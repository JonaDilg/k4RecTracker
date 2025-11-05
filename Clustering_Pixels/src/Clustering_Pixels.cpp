#include "Clustering_Pixels.h"

DECLARE_COMPONENT(Clustering_Pixels)

Clustering_Pixels::Clustering_Pixels(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("TrackerHitCollectionName", {"UNDEFINED_TrackerHitCollectionName"}),
                        KeyValues("HeaderName", {"UNDEFINED_HeaderName"}),},
                       {KeyValues("TrackerClusterHitCollectionName", {"UNDEFINED_TrackerClusterHitCollectionName"}),
                        KeyValues("SimTrkHitRelationsCollection", {"UNDEFINED_SimTrkHitRelationsCollection"})}) {
  info() << "Constructed successfully" << endmsg;
}

StatusCode Clustering_Pixels::initialize() {
  info() << "INITIALIZING Clustering_Pixels ..." << endmsg;

  InitializeServicesAndGeometry();

  CheckGaudiProperties();

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
  (const edm4hep::TrackerHitPlaneCollection& hits, const edm4hep::EventHeaderCollection& headers) const {

  /* Initial check: returns false if event has no simHits, throws if there is a problem with the setup */
  if (!CheckInitialSetup(hits, headers))
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());

  /* output collections */
  auto clusters = edm4hep::TrackerHitPlaneCollection();
  auto clusterHitLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();

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
    const int layer = m_cellIDdecoder->get(hit.getCellID(), "layer");
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
  for (auto& [sensorUid, hitsOnSensor] : hitsSensorMap) {

    if (hitsOnSensor.empty())
      throw GaudiException("Internal error: hitsOnSensor is empty in Clustering_Pixels::operator()", "Clustering_Pixels::operator()", StatusCode::FAILURE);
    if (hitsOnSensor.front() == nullptr)
      throw GaudiException("Internal error: first hit pointer is null in Clustering_Pixels::operator()", "Clustering_Pixels::operator()", StatusCode::FAILURE);

    int nHits = hitsOnSensor.size();
    const int layer = m_cellIDdecoder->get(hitsOnSensor.front()->getCellID(), "layer");
    const int cellID = hitsOnSensor.front()->getCellID();

    verbose() << "   - Processing hits on (SensorUID " << sensorUid << ", layer " << layer << ", cellID " << cellID << "): found " << nHits << " hits." << endmsg;

    /* sort pixels by energy deposition (last element has highest charge) */
    std::sort(
      hitsOnSensor.begin(), hitsOnSensor.end(),
      [](const std::shared_ptr<const edm4hep::TrackerHitPlane>& a, const std::shared_ptr<const edm4hep::TrackerHitPlane>& b) {
        return a->getEDep() < b->getEDep();
      });

    /* TODO: set up conversion from global coordinates to pixel indices (i_u, i_v) for this sensor */

    /* TODO: create data structure that has (i_u, i_v, energy, timestamp, ptr to hit) for each hit */

    std::vector<std::shared_ptr<const edm4hep::TrackerHitPlane>> hitsInCluster; 
    hitsInCluster.reserve(9); // reserve space for 3x3 cluster
    int nClusters = 0;

    /* clustering algorithm
     * this follows the structure of the Corryvreckan module ClusteringSpatial 
     * (see https://gitlab.cern.ch/corryvreckan/corryvreckan/-/tree/master/src/modules/ClusteringSpatial?ref_type=heads) */
    while (!hitsOnSensor.empty()) {
      
      /* TODO: implement charge-threshold and time-cut */
      
      verbose() << "      - Starting new cluster with " << hitsOnSensor.size() << " unclustered hits remaining." << endmsg;
      ++nClusters;
      hitsInCluster.clear();
      hitsInCluster.push_back(hitsOnSensor.back()); // modifying the last object in vector is much faster than modifying the first
      hitsOnSensor.pop_back();

      /* look for neighbors, add them, look for neighbors of newly added neighbors, ...
       * until no new neighbors are found */
      bool addedNeighbor = true;
      while (addedNeighbor) {
        addedNeighbor = false;

        /* TODO: get position in terms of pixels (or in terms of u,v position?)*/

        /* TODO: find neighboring hits and add them to hitsInCluster */
      }

      /* TODO: get cluster position and timestamp from hitsInCluster */

      /* TODO: create hit & links for the cluster */
    }

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
  
  const dd4hep::Detector* detector = m_geometryService->getDetector();
  if (!detector)
    throw GaudiException("Unable to retrieve the DD4hep detector from GeoSvc", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  const dd4hep::rec::SurfaceManager* simSurfaceManager = detector->extension<dd4hep::rec::SurfaceManager>();
  if (!simSurfaceManager)
    throw GaudiException("Unable to retrieve the SurfaceManager from the DD4hep detector", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  m_simSurfaceMap = simSurfaceManager->map(m_subDetName.value());
  if (!m_simSurfaceMap)
    throw GaudiException("Unable to retrieve the simSurface map for subdetector " + m_subDetName.value(), "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  m_volumeManager = detector->volumeManager();
  if (!m_volumeManager.isValid())
    throw GaudiException("Unable to retrieve the VolumeManager from the DD4hep detector", "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  /* subDetector not needed as of now. Keep for future reference ~ Jona 2025-10 */
  // const dd4hep::DetElement subDetector = detector->detector(m_subDetName.value());
  // if (!subDetector.isValid())
  //   throw GaudiException("Unable to retrieve the DetElement for subdetector " + m_subDetName.value(), "Clustering_Pixels::initializeServicesAndGeometry()", StatusCode::FAILURE);

  debug() << " - Successfully retrieved all necessary services and detector elements, starting to check Gaudi properties." << endmsg;
}

void Clustering_Pixels::CheckGaudiProperties() {
  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "Clustering_Pixels::checkGaudiProperties()", StatusCode::FAILURE);

  if (m_layersToRun.value().empty()) {
    /* TODO: Find all layers in geometry & set m_layerCount accordingly. right now this doesnt work.*/
    throw GaudiException("Finding all layers in geometry is not implemented yet. Please specify the layers to be digitized via the LayersToDigitize property.", "VTXdigi_Allpix2::setupSensorParameters()", StatusCode::FAILURE);

    m_layerN = 0; // placeholder, needs to be set properly when finding all layers in geometry is implemented
  }
  else { // if layers are specified as Gaudi property
    m_layerN = m_layersToRun.size();
    verbose() << " - Digitizing " << m_layerN << " layers, as specified by LayersToRun property." << endmsg;
  }
}

void Clustering_Pixels::SetupDebugHistograms() {
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
  info() << "PROCESSING event. run " << headers.at(0).getRunNumber() << ", event " << headers.at(0).getEventNumber() << ". Found " << hits.size() << " hits." << endmsg;
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



