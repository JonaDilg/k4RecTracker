#include "Clustering_Pixels.h"

DECLARE_COMPONENT(Clustering_Pixels)

Clustering_Pixels::Clustering_Pixels(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("TrackerHitCollectionName", {"UNDEFINED_TrackerHitCollectionName"}),
                        KeyValues("TrackerHitSimTrackerHitLinkCollectionName", {"UNDEFINED_TrackerHitSimTrackerHitLinkCollectionName"}),
                        KeyValues("HeaderName", {"UNDEFINED_HeaderName"}),},
                       {KeyValues("TrackerClusterHitCollectionName", {"UNDEFINED_TrackerClusterHitCollectionName"}),
                        KeyValues("SimTrkHitRelationsCollection", {"UNDEFINED_SimTrkHitRelationsCollection"})}) {
  info() << "Constructed successfully" << endmsg;
}

StatusCode Clustering_Pixels::initialize() {
  info() << "INITIALIZING..." << endmsg;

  InitializeServicesAndGeometry();

  CheckGaudiProperties();

  if (m_debugHistograms)
    InitHistograms();

  verbose() << " - Initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode Clustering_Pixels::finalize() {
  info() << "FINALIZING Clustering_Pixels ..." << endmsg;

  PrintCountersSummary();

  verbose() << " - finalized successfully" << endmsg;
  return StatusCode::SUCCESS;
} 

/* -- event loop -- */

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> Clustering_Pixels::operator()
  (const edm4hep::TrackerHitPlaneCollection& hits, const edm4hep::TrackerHitSimTrackerHitLinkCollection& hitLinks, const edm4hep::EventHeaderCollection& headers) const {

  /* Initial check: returns false if event has no simHits, throws if there is a problem with the setup */
  if (!CheckInitialSetup(hits, headers))
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());

  /* output collections */
  auto clusters = edm4hep::TrackerHitPlaneCollection();
  auto clusterHitLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();

  /* link navigator (to find which hits are linked to which sim hits) */
  const auto hitLinkNavigator = podio::LinkNavigator(hitLinks);

  /* sort hits to their respective sensor
   * Note: do NOT copy the objects behind the pointers, this messes up their associations under the hood of PODIO!
   * instead, simply access via member access -> */
  // std::map<int, std::vector<const edm4hep::TrackerHitPlane*>> hitsSensorMap;
  std::map<int, std::vector<std::shared_ptr<const edm4hep::TrackerHitPlane>>> hitsSensorMap;

  for (const auto& hit : hits) {
    ++m_counter_hitsRead;

    if (!CheckHitCuts(hit))
      continue; // counters updated in CheckHitCuts()
    ++m_counter_hitsAccepted;
    
    /* We need a unique identifier for each sensor, so that we can group hits by sensor
    * this is because we only want to cluster hits on the same sensor
    * TODO what I am doing here is only a temporary solution, needs to be replaced with a proper way to get a unique identifier for each sensor
    * (this works for IDEA with ARCADIA and ATLASpix layers, because ATLASpix has 4 and ARCADIA has 1 sensor per module)*/
    const int module = m_cellIDdecoder->get(hit.getCellID(), "module");
    const int sensor = m_cellIDdecoder->get(hit.getCellID(), "sensor");
    const int sensorUid = module*4 + sensor; // unique identifier for each sensor
    
    /* create managed copy, will lead to segfaults as soon as this loop goes out of scope */
    auto hitPtr = std::make_shared<edm4hep::TrackerHitPlane>(hit);
    hitsSensorMap[sensorUid].push_back(hitPtr);
  }

  /* Check if we have any hits that are on layers we want to cluster */
  if (hitsSensorMap.empty()) {
    debug() << " - No hits passed the cuts, returning empty output collections." << endmsg;
    ++m_counter_eventsRejected_noHitsOnSpecifiedLayers;
    return std::make_tuple(std::move(clusters), std::move(clusterHitLinks));
  }
  verbose() << " - Sorted " << hits.size() << " hits to " << hitsSensorMap.size() << " sensors. Looping over sensors:" << endmsg;


  
  /* loop over sensors, create clusters */
  for (auto& [sensorUid, SensorHits] : hitsSensorMap) {

    if (SensorHits.empty())
      throw GaudiException("Internal error: SensorHits is empty in Clustering_Pixels::operator()", "Clustering_Pixels::operator()", StatusCode::FAILURE);
    if (SensorHits.front() == nullptr)
      throw GaudiException("Internal error: first hit pointer is null in Clustering_Pixels::operator()", "Clustering_Pixels::operator()", StatusCode::FAILURE);

    int nHits = SensorHits.size();
    const int layer = m_cellIDdecoder->get(SensorHits.front()->getCellID(), "layer");
    const dd4hep::DDSegmentation::CellID cellID = SensorHits.front()->getCellID();

    verbose() << "   - Processing hits on (SensorUID " << sensorUid << ", layer " << layer << ", cellID " << cellID << "): found " << nHits << " hits." << endmsg;

    /* sort pixels by energy deposition (last element has highest charge) */
    std::sort(
      SensorHits.begin(), SensorHits.end(),
      [](const std::shared_ptr<const edm4hep::TrackerHitPlane>& a, const std::shared_ptr<const edm4hep::TrackerHitPlane>& b) {
        return a->getEDep() < b->getEDep();
      });

    /* TODO: set up conversion from global coordinates to pixel indices (i_u, i_v) for this sensor 
    * This would require knowing the pixel pitch and sensor size. the pitch is not directly accessible as of now.
    * (for now, i simply use the local position )*/

    /* place all accepted hits into struct, compute their local coordinates */
    std::vector<HitData> SensorHitsData;
    SensorHitsData.reserve(nHits);

    for (const auto& hit : SensorHits) {
      HitData data;
      dd4hep::rec::Vector3D globalPos = ConvertVector(hit->getPosition());
      dd4hep::rec::Vector3D localPos = TransformGlobalToLocal(globalPos, hit->getCellID());
      data.u = localPos.x();
      data.v = localPos.y();
      data.eDep = hit->getEDep();
      data.time = hit->getTime();
      data.hitPtr = hit;
      SensorHitsData.push_back(data);
    }

    std::vector<HitData> clusterHitsData;
    clusterHitsData.reserve(9); // reserve space for 3x3 cluster
    int nClusters = 0;

    /* clustering algorithm
     * this follows the structure of the Corryvreckan module ClusteringSpatial 
     * (see https://gitlab.cern.ch/corryvreckan/corryvreckan/-/tree/master/src/modules/ClusteringSpatial?ref_type=heads) */
    while (!SensorHitsData.empty()) {
      
      /* TODO: implement charge-threshold and time-cut */

      verbose() << "      - Starting new cluster with " << SensorHitsData.size() << " unclustered hits remaining." << endmsg;
      ++nClusters;
      clusterHitsData.clear();
      clusterHitsData.push_back(SensorHitsData.back()); // modifying the last object in vector is much faster than modifying the first
      SensorHitsData.pop_back();

      /* look for neighbors, add them, look for neighbors of newly added neighbors, ...
       * until no new neighbors are found */
      bool addedNeighbor = true;
      while (addedNeighbor) {
        addedNeighbor = false;
        /* loop over hits and see if any are closer than r_u^2+r_v^2
         * this is a brute-force algorithm. My naive guess is that our low occupancy makes this a non-issue. ~ Jona, 2025-11*/
        for (auto iterator = SensorHitsData.begin(); iterator != SensorHitsData.end(); ) {
          const float du = iterator->u - clusterHitsData.back().u; // local u coordinate of last added hit
          const float dv = iterator->v - clusterHitsData.back().v; // local v coordinate of last added hit

          /* use ellipse equation (du/r_u)^2 + (dv/r_v)^2 <= 1 to check if hit is within neighbor radius in u and v*/
          if ((du*du*m_neighborRadiusHelper[1] + dv*dv*m_neighborRadiusHelper[0]) <= m_neighborRadiusHelper[2]) {
            verbose() << "         - Found neighbor hit at distance (du, dv) = (" << du << " mm, " << dv << " mm), adding to cluster." << endmsg;
            clusterHitsData.push_back(*iterator);
            iterator = SensorHitsData.erase(iterator); // erase returns the next iterator
            addedNeighbor = true;
          }
          else {
            ++iterator;
            verbose() << "         - Hit at distance (du, dv) = (" << du << " mm, " << dv << " mm) is too far away." << endmsg;
          }
        }
      }

      /* compute cluster position and timestamp from clusterHitsData */
      float clusterPos[2] = {0.,0.}, clusterEDep = 0, clusterTime = 0;
      for (const auto& hitData : clusterHitsData) {
        clusterPos[0] += hitData.u * hitData.eDep;
        clusterPos[1] += hitData.v * hitData.eDep;
        clusterEDep += hitData.eDep;
        /* TODO: weighed average is a crude estimate for cluster timestamp. Think of something better. */
        clusterTime += hitData.time * hitData.eDep; 
      }
      clusterPos[0] /= clusterEDep;
      clusterPos[1] /= clusterEDep;
      clusterTime /= clusterEDep;

      /* collect simHits that are linked to any hit in this cluster. (see
      * https://github.com/AIDASoft/podio/blob/master/doc/links.md 
      * https://github.com/AIDASoft/podio/blob/1a678f5f46273e3e5a2ea3ff16eb8d41990e7c70/include/podio/LinkNavigator.h#L87) 
      * Problem: I am now copying the simHits, and creating hitLinks from the clusterHits to these copied simHits.
      * I am not sure if this works correctly, but I can imagine it does. 
      * There seems to be no better way to do this (ie. I cannot find one) ~ Jona 2025-12 */
      std::vector<edm4hep::SimTrackerHit> clusterSimHits;

      for (const auto& hitData : clusterHitsData) {
        std::vector<
          podio::detail::links::WeightedObject<edm4hep::SimTrackerHit>,std::allocator<podio::detail::links::WeightedObject<edm4hep::SimTrackerHit>>
        > linkedWs = hitLinkNavigator.getLinked(*hitData.hitPtr); // full type just so I know what is going on ~ Jona 2025-12

        for (const auto& linkedW : linkedWs) {
          const edm4hep::SimTrackerHit& linkedSimHit = linkedW.o; // might be .object depending on podio version (i think?) ~ Jona 2025-12

          /* check if we already have this simHit in clusterSimHits, if not, add it. */
          bool isDuplicate = false;
          for (const auto& simHit : clusterSimHits) {
            if (simHit == linkedSimHit) {
              isDuplicate = true;
              break;
            }
          }
          if (!isDuplicate) {
            clusterSimHits.push_back(linkedSimHit);
          }
        }
      }
      const int nLinkedSimHits = clusterSimHits.size();
      verbose() << "      - Cluster has " << nLinkedSimHits << " unique linked simHits." << endmsg;


      
      /* create cluster hit */

      auto cluster = clusters.create();
      cluster.setCellID(cellID);
      dd4hep::rec::Vector3D clusterLocalPos(clusterPos[0], clusterPos[1], 0.f);
      dd4hep::rec::Vector3D clusterGlobalPos = TransformLocalToGlobal(clusterLocalPos, cellID);
      cluster.setPosition(ConvertVector(clusterGlobalPos));
      cluster.setEDep(clusterEDep);
      cluster.setTime(clusterTime);

      /* get pointers to all simHits that are linked to by digiHits in the cluster */
      



      /* TODO: loop over hits in this cluster, find all linked simHits */


      // auto clusterLink = clusterLinks.create();
      // clusterLink.setFrom(cluster);


      
      // for (const auto& hitData : clusterHitsData) {
      //   clusterLink.setTo(hitData.hitPtr);
      // }

      if (m_debugHistograms.value()) {
        FillHistograms_perCluster(cluster, clusterHitsData, clusterSimHits, cellID, layer);
      }
    } // end clustering loop

    verbose() << "   - Created " << nClusters << " clusters." << endmsg;
  }

  debug() << "FINISHED event." << endmsg;
  return std::make_tuple(std::move(clusters), std::move(clusterHitLinks));
} // operator()


/* -- Initialization / finalization functions -- */

void Clustering_Pixels::InitializeServicesAndGeometry() {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc)
    throw GaudiException("Unable to get UniqueIDGenSvc", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  m_geometryService = serviceLocator()->service(m_geometryServiceName);
  if (!m_geometryService)
    throw GaudiException("Unable to retrieve the GeoSvc", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  std::string cellIDstr = m_geometryService->constantAsString(m_encodingStringVariable.value());
  m_cellIDdecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDstr);
  if (!m_cellIDdecoder)
    throw GaudiException("Unable to retrieve the cellID decoder", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  m_detector = m_geometryService->getDetector();
  if (!m_detector)
    throw GaudiException("Unable to retrieve the DD4hep detector from GeoSvc", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  const dd4hep::rec::SurfaceManager* simSurfaceManager = m_detector->extension<dd4hep::rec::SurfaceManager>();
  if (!simSurfaceManager)
    throw GaudiException("Unable to retrieve the SurfaceManager from the DD4hep detector", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  m_simSurfaceMap = simSurfaceManager->map(m_subDetName.value());
  if (!m_simSurfaceMap)
    throw GaudiException("Unable to retrieve the simSurface map for subdetector " + m_subDetName.value(), "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  m_volumeManager = m_detector->volumeManager();
  if (!m_volumeManager.isValid())
    throw GaudiException("Unable to retrieve the VolumeManager from the DD4hep detector", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  /* IDEA / Allegro: 
   *   - subDet is child of detector, eg. Vertex 
   *   - subDetChild is child of subDet, eg. VertexBarrel
   *   - subDetChildChild are layers 
   * If this is not the case, layers might be direct children of subDet: */

  if (m_subDetChildName.value() != m_undefinedString) {
    /* IDEA/Allegro setup */
    const dd4hep::DetElement subDet = m_detector->detector(m_subDetName.value());
    m_subDetector = subDet.child(m_subDetChildName.value());
  }
  else {
    /* layers are direct children of subDet */
    m_subDetector = m_detector->detector(m_subDetName.value());
  }

  if (!m_subDetector)
    throw GaudiException("Unable to retrieve the subdetector DetElement " + m_subDetName.value(), "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  debug() << " - Successfully retrieved all necessary services and detector elements, starting to check Gaudi properties." << endmsg;
}

void Clustering_Pixels::CheckGaudiProperties() {
  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "Clustering_Pixels::checkGaudiProperties()", StatusCode::FAILURE);

  if (m_layersToRun.value().empty()) {
    /* TODO: Find all layers in geometry & set m_layerCount accordingly. right now this doesnt work.*/
    throw GaudiException("Finding all layers in geometry is not implemented yet. Please specify the layers to be digitized via the LayersToDigitize property.", "Clustering_Pixels::setupSensorParameters()", StatusCode::FAILURE);

    /* Collect layers in this subdetector from geometry */


    m_layerN = 0; // placeholder, needs to be set properly when finding all layers in geometry is implemented
  }
  else { // if layers are specified as Gaudi property
    m_layerN = m_layersToRun.size();
    verbose() << " - Digitizing " << m_layerN << " layers, as specified by LayersToRun property." << endmsg;
  }

  /* pre-compute helper for checking if a hit is within the neighbor radius, using the ellipse equation (du/r_u)^2 + (dv/r_v)^2 <= 1:
   * is neighbor if: du^2*[1] + dv^2*[0]) <= [2] */
  m_neighborRadiusHelper[0] = m_neighborRadius_u.value()*m_neighborRadius_u.value() / 1e6; // convert from um to mm and square
  m_neighborRadiusHelper[1] = m_neighborRadius_v.value()*m_neighborRadius_v.value() / 1e6;
  m_neighborRadiusHelper[2] = m_neighborRadiusHelper[0] * m_neighborRadiusHelper[1] + 0.001f; // add small numeric limit to avoid floating point issues. We can simply add a small epsilon, because we know the order of magnitude of the distances we are working with.
}

void Clustering_Pixels::InitHistograms() {
  
  /* -- per-layer histograms -- */
  for (int layer : m_layersToRun.value()) {
    std::array< std::unique_ptr< Gaudi::Accumulators::StaticHistogram<1,Gaudi::Accumulators::atomicity::full,float>>, hist1dArrayLen>  hist1d;

    hist1d.at(hist1d_clusterSize).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize",
        "Cluster size - Layer " + std::to_string(layer) + ";Pixels per cluster;Entries",
        {50, -0.5f, 49.5f}
      }
    );

    hist1d.at(hist1d_clusterCharge).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterCharge",
        "Cluster charge - Layer " + std::to_string(layer) + ";Total cluster charge [e];Entries",
        {1000, 0.f, 50.f*500.f} // for 50 mu sensor thickness
      }
    );

    hist1d.at(hist1d_simHitsPerCluster).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHitsPerCluster",
        "Number of individual particles that contributed to this cluster - Layer " + std::to_string(layer) + ";SimHits per cluster;Entries",
        {30, -0.5f, 29.5f}
      }
    );

    hist1d.at(hist1d_residual).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual",
        "Total residual - Layer " + std::to_string(layer) + ";Residual [um];Entries",
        {2000, -1000.f, 1000.f}
      }
    );
    hist1d.at(hist1d_residual_z).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_Global_z",
        "Residual in z (global) - Layer " + std::to_string(layer) + ";Residual z [um];Entries",
        {2000, -1000.f, 1000.f}
      }
    );
    hist1d.at(hist1d_residual_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_Local_u",
        "Residual in u (local) - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        {2000, -1000.f, 1000.f}
      }
    );
    hist1d.at(hist1d_residual_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_Local_v",
        "Residual in v (local) - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        {2000, -1000.f, 1000.f}
      }
    );
    hist1d.at(hist1d_residual_w).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_Local_w",
        "Residual in w (local) - Layer " + std::to_string(layer) + ";Residual w [um];Entries",
        {2000, -1000.f, 1000.f}
      }
    );



    m_hist1d.emplace(layer, std::move(hist1d));

    /* -- 2D histograms -- */

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticHistogram<2,Gaudi::Accumulators::atomicity::full,float>>, hist2dArrayLen>  hist2d;

    hist2d.at(hist2d_clusterSize_vs_clusterCharge).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_clusterCharge",
        "Cluster size vs. cluster charge - Layer " + std::to_string(layer) + ";Total cluster charge [e];Pixels per cluster",
        {1000, 0.f, 50.f*500.f},
        {50, -0.5f, 49.5f} // for 50 mu sensor thickness
      }
    );

    hist2d.at(hist2d_residual_local).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_2D_local",
        "Residual in local (sensor) coordinates - Layer " + std::to_string(layer) + ";u [um];v [um]",
        {2000, -1000.f, 1000.f},
        {2000, -1000.f, 1000.f}
      }
    );



    
    m_hist2d.emplace(layer, std::move(hist2d));

    /* -- 1D Profile histograms -- */

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticProfileHistogram<1,Gaudi::Accumulators::atomicity::full,float>>, histProfile1dArrayLen>  histProfile1d;

    histProfile1d.at(histProfile1d_clusterSize_vs_hit_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_hit_z",
        "Cluster size - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Pixels per cluster",
        {2000, -500.f, 500.f}
      }
    );
    histProfile1d.at(histProfile1d_clusterSize_vs_hit_cosTheta).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_hit_cosTheta",
        "Cluster size - Layer " + std::to_string(layer) + ";cosTheta of simHit position;Pixels per cluster",
        {200, 0.f, 1.f}
      }
    );

    histProfile1d.at(histProfile1d_clusterSize_vs_module_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_module_z",
        "Cluster size - Layer " + std::to_string(layer) + ";module global z position [mm];Pixels per cluster",
        {2000, -500.f, 500.f}
      }
    );
    histProfile1d.at(histProfile1d_clusterSize_vs_module_ID).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_module_ID",
        "Cluster size - Layer " + std::to_string(layer) + ";Module ID;Pixels per cluster",
        {2000, -0.5f, 1999.5f}
      }
    );

    histProfile1d.at(histProfile1d_residual_vs_hit_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_vs_hit_z",
        "Total residual - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Residual [um]",
        {2000, -500.f, 500.f}
      }
    );
    histProfile1d.at(histProfile1d_residual_u_vs_hit_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_u_vs_hit_z",
        "Local residual in u - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Residual in u (sensor frame) [um]",
        {2000, -500.f, 500.f}
      }
    );
    histProfile1d.at(histProfile1d_residual_v_vs_hit_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_v_vs_hit_z",
        "Local residual in v - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Residual in v (sensor frame) [um]",
        {2000, -500.f, 500.f}
      }
    );
    histProfile1d.at(histProfile1d_residual_vs_hit_cosTheta).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_vs_hit_cosTheta",
        "Total residual - Layer " + std::to_string(layer) + ";SimHit global cosTheta;Residual [um]",
        {200, 0.f, 1.f}
      }
    );
    histProfile1d.at(histProfile1d_residual_vs_clusterSize).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/residual_vs_clusterSize",
        "Total residual - Layer " + std::to_string(layer) + ";Cluster size [pixels];Residual [um]",
        {50, -0.5f, 49.5f}
      }
    );



    m_histProfile1d.emplace(layer, std::move(histProfile1d));
  }


}

void Clustering_Pixels::PrintCountersSummary() const {
  const int colWidths[] = {65, 10};  
  info() << " Counters summary: " << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events read"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsRead.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events rejected (no hits)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsRejected_noHits.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events rejected (no hits above seed threshold)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsRejected_noHitsAboveSeedThreshold.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Tracker hits read"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_hitsRead.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Tracker hits rejected (on layer that is ignored)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_hitsRejected_layerIgnored.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Tracker hits accepted"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_hitsAccepted.value() << " |" << endmsg;
}

/* -- Core algorithm functions -- */
bool Clustering_Pixels::CheckInitialSetup(const edm4hep::TrackerHitPlaneCollection& hits, const edm4hep::EventHeaderCollection& headers) const {
  info() << "PROCESSING event (run " << headers.at(0).getRunNumber() << ", event " << headers.at(0).getEventNumber() << ", found " << hits.size() << " hits" << endmsg;
  ++m_counter_eventsRead;

  // early sanity checks to avoid segfaults from null pointers
  if (!m_simSurfaceMap)
    throw GaudiException("SimSurfaceMap is null in operator(). Did initialize() succeed?", "Clustering_Pixels::CheckInitialSetup()", StatusCode::FAILURE);

  if (hits.size()==0) {
    debug() << " - No SimTrackerHits in collection, returning empty output collections" << endmsg;
    ++m_counter_eventsRejected_noHits;
    return false;
  }

  return true;
}


/* -- Transformation between global frame and local sensor frame -- */

TGeoHMatrix Clustering_Pixels::ComputeTransformationMatrix(const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = m_volumeManager.lookupDetElement(cellID).nominal().worldTransformation(); // given in cm

  // rotation is unitless, but need to convert translation from cm to mm
  double* translationComponent = transformationMatrix.GetTranslation();
  translationComponent[0] = translationComponent[0] * 10; // convert to mm
  translationComponent[1] = translationComponent[1] * 10;
  translationComponent[2] = translationComponent[2] * 10;
  transformationMatrix.SetTranslation(translationComponent);

  SetProperDirectFrame(transformationMatrix); // Change coordinates to have z orthogonal to sensor, with direct (right-handed) frame

  return transformationMatrix;
}

void Clustering_Pixels::SetProperDirectFrame(TGeoHMatrix& transformationMatrix) const {
  /** Change the transformationMatrix to have a direct frame with z orthogonal to sensor simSurface
   * 
   * copied from https://github.com/jessy-daniel/k4RecTracker/blob/New_Detailed_VTXdigitizer/VTXdigiDetailed/src/VTXdigitizerDetailed.cpp
   * I don't really understand what it does. Defaults to do nothing though, that's ok with me. ~ Jona, 2025-09
   */
  
  std::string localNormalVectorDir = m_localNormalVectorDir;
  bool IsDirect = true; // Is the origin frame direct ?
  if (localNormalVectorDir[0]=='-') {
    IsDirect = false;
    localNormalVectorDir = localNormalVectorDir[1];
  }  

  // If the orthogonal direction is along X or Y in local frame, rotate the frame to have Z orthogonal to sensors instead
  if (localNormalVectorDir=="x") {
    TGeoRotation rot("rot",90.,90.,0.); // X->Z / Y->X / Z->Y
    transformationMatrix.Multiply(rot);
  }
  if (localNormalVectorDir=="y") {
    TGeoRotation rot("rot",0.,-90.,0.); // X->X / Y->Z / Z->-Y
    transformationMatrix.Multiply(rot);
  }

  // If the frame isn't direct, make it direct by reflecting the x axis. This is necessary to correctly calculte the drift in X-Y due to B-field
  if (!IsDirect) {
    transformationMatrix.ReflectX(false);
  }
} // setProperDirectFrame()

dd4hep::rec::Vector3D Clustering_Pixels::TransformGlobalToLocal(const dd4hep::rec::Vector3D& globalPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = ComputeTransformationMatrix(cellID);

  double localPos[3] = {0, 0, 0};

  transformationMatrix.MasterToLocal(globalPos, localPos);

  return dd4hep::rec::Vector3D(localPos[0], localPos[1], localPos[2]);
}

dd4hep::rec::Vector3D Clustering_Pixels::TransformLocalToGlobal(const dd4hep::rec::Vector3D& localPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = ComputeTransformationMatrix(cellID);

  double globalPos[3] = {0, 0, 0};

  transformationMatrix.LocalToMaster(localPos, globalPos);

  return dd4hep::rec::Vector3D(globalPos[0], globalPos[1], globalPos[2]);
}

dd4hep::rec::Vector3D Clustering_Pixels::ConvertVector(edm4hep::Vector3d vec) const {
  // return dd4hep::rec::Vector3D(vec.x*dd4hep::mm, vec.y*dd4hep::mm, vec.z*dd4hep::mm);
  return dd4hep::rec::Vector3D(vec.x, vec.y, vec.z);
}
edm4hep::Vector3d Clustering_Pixels::ConvertVector(dd4hep::rec::Vector3D vec) const {
  // return edm4hep::Vector3d(vec[0]/dd4hep::mm, vec[1]/dd4hep::mm, vec[2]/dd4hep::mm);
  return edm4hep::Vector3d(vec.x(), vec.y(), vec.z());
}

void Clustering_Pixels::FillHistograms_perCluster(const edm4hep::MutableTrackerHitPlane& cluster, const std::vector<Clustering_Pixels::HitData>& clusterHitsData,  std::vector<edm4hep::SimTrackerHit> clusterSimHits, const dd4hep::DDSegmentation::CellID& cellID, const int layer) const {
  verbose() << "   - Filling debug histograms for cluster in layer " << layer << " with " << clusterHitsData.size() << " hits." << endmsg;

  const dd4hep::rec::Vector3D clusterGlobalPos = ConvertVector(cluster.getPosition()); // re-computing them here. This is slightly slower than passing them, but cleaner.
  const dd4hep::rec::Vector3D clusterLocalPos = TransformGlobalToLocal(clusterGlobalPos, cellID);

  ++(*m_hist1d.at(layer).at(hist1d_simHitsPerCluster))[clusterSimHits.size()]; 

  /* ClusterSize */
  const int clusterSize = clusterHitsData.size();
  ++(*m_hist1d.at(layer).at(hist1d_clusterSize))[static_cast<float>(clusterSize)]; 

  /* ClusterCharge */
  float clusterE = 0.0f; // in keV
  for (const auto& hit : clusterHitsData) {
    clusterE += hit.eDep;
  }
  const float clusterCharge = clusterE * m_chargePerkeV; // in electrons
  verbose() << "      - Cluster eDep " << clusterE << " keV, charge " << clusterCharge << " e." << endmsg;
  ++(*m_hist1d.at(layer).at(hist1d_clusterCharge))[clusterCharge];
  ++(*m_hist2d.at(layer).at(hist2d_clusterSize_vs_clusterCharge))[{clusterCharge, clusterSize}];

  /* plot sim-hit dependent variables once for every simHit that contributed to the cluster */
  for (const auto& simHit : clusterSimHits) {
    const dd4hep::rec::Vector3D simHitGlobalPos = ConvertVector(simHit.getPosition());
    const dd4hep::rec::Vector3D simHitLocalPos = TransformGlobalToLocal(simHitGlobalPos, cellID);

    const float simHit_z = simHitGlobalPos.z(); // in mm
    (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_vs_hit_z))[simHit_z] += clusterSize;

    const float simHit_cosTheta = std::abs(simHitGlobalPos.z() / simHitGlobalPos.r());
    (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_vs_hit_cosTheta))[simHit_cosTheta] += clusterSize;

    /* Residual = cluster pos. - truth pos. */
    const dd4hep::rec::Vector3D residualLocal = clusterLocalPos - simHitLocalPos;
    const dd4hep::rec::Vector3D residualGlobal = clusterGlobalPos - simHitGlobalPos;

    ++(*m_hist1d.at(layer).at(hist1d_residual))[residualGlobal.r()*1000.f]; // convert to um
    ++(*m_hist1d.at(layer).at(hist1d_residual_z))[residualGlobal.z()*1000.f];
    ++(*m_hist1d.at(layer).at(hist1d_residual_u))[residualLocal.x()*1000.f];
    ++(*m_hist1d.at(layer).at(hist1d_residual_v))[residualLocal.y()*1000.f];
    ++(*m_hist1d.at(layer).at(hist1d_residual_w))[residualLocal.z()*1000.f];

    ++(*m_hist2d.at(layer).at(hist2d_residual_local))[{residualLocal.x()*1000.f, residualLocal.y()*1000.f}];

    (*m_histProfile1d.at(layer).at(histProfile1d_residual_vs_hit_z))[simHit_z] += residualGlobal.r()*1000.f; // convert mm to um
    (*m_histProfile1d.at(layer).at(histProfile1d_residual_u_vs_hit_z))[simHit_z] += residualLocal.x()*1000.f;
    (*m_histProfile1d.at(layer).at(histProfile1d_residual_v_vs_hit_z))[simHit_z] += residualLocal.y()*1000.f;
    (*m_histProfile1d.at(layer).at(histProfile1d_residual_vs_hit_cosTheta))[simHit_cosTheta] += residualGlobal.r()*1000.f;
    (*m_histProfile1d.at(layer).at(histProfile1d_residual_vs_clusterSize))[clusterSize] += residualGlobal.r()*1000.f;
  }

  const float module_z = TransformLocalToGlobal(dd4hep::rec::Vector3D(0.f, 0.f, 0.f), cellID).z();
  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_vs_module_z))[module_z] += clusterSize;

  const int moduleID = m_cellIDdecoder->get(cellID, "module");
  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_vs_module_ID))[static_cast<float>(moduleID)] += clusterSize;

}

/* -- Other Helper functions -- */

int Clustering_Pixels::ComputeBinIndex(float x, float binX0, float binWidth, int binN) const {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  if (binN <= 0) throw GaudiException("computeBinIndex: binN must be positive", "Clustering_Pixels::computeBinIndex()", StatusCode::FAILURE);
  if (binWidth <= 0.0) throw GaudiException("computeBinIndex: binWidth must be positive", "Clustering_Pixels::computeBinIndex()", StatusCode::FAILURE);

  float relativePos = (x - binX0) / binWidth;
  if (relativePos < 0.0f || relativePos > static_cast<float>(binN))
    return -1;
  return static_cast<int>(relativePos);
} // computeBinIndex()

// std::tuple<int, int> Clustering_Pixels::computePixelIndices(const edm4hep::TrackerHitPlane& hit) const {
  
//   int i_u = computeBinIndex(
//     segmentPos.x(),
//     -0.5*length_u,
//     m_pixelPitch_u.at(layerIndex),
//     m_pixelCount_u.value().at(layerIndex));

//   int i_v = computeBinIndex(
//     segmentPos.y(),
//     -0.5*length_v,
//     m_pixelPitch_v.at(layerIndex),
//     m_pixelCount_v.value().at(layerIndex));
//   return std::make_tuple(i_u, i_v);
// } // computePixelIndices()

bool Clustering_Pixels::CheckHitCuts(const edm4hep::TrackerHitPlane& hit) const {
  // verbose() << "   - Checking hit cuts for hit with EDep = " << hit.getEDep() << " keV" << endmsg; // gives segFault ???

  const int layer = m_cellIDdecoder->get(hit.getCellID(), "layer");
  if (m_layerN != -1) { // if layers are specified
    if (std::find(m_layersToRun.value().begin(), m_layersToRun.value().end(), layer) == m_layersToRun.value().end()) {
      debug() << "   - Hit in layer " << layer << " skipped, as this layer is not in LayersToRun property." << endmsg;
      ++m_counter_hitsRejected_layerIgnored;
      return false;
    }
  }

  return true;
}

void Clustering_Pixels::CreateCluster() const {

}



