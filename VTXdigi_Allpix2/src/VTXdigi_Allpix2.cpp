#include "VTXdigi_Allpix2.h"

DECLARE_COMPONENT(VTXdigi_Allpix2)

/* Notes ~ Jona 2025-09
 * - all lengths in mm (edm4hep already does this)
 *   BUT: dd4hep uses cm internally, so convert when passing values to/from dd4hep via dd4hep::mm = 0.1
 *        mm -> cm: a [cm] = dd4hep::mm * a [mm]
 *        cm -> mm: a [mm] = 1/dd4hep::mm * a [cm]

 * - Vectors can be given in 
 *        a) dd4hep::rec::Vector3D <- fully featured vector, overloads operators *+- etc
 *        b) edm4hep::Vector3d <- natively used by edm4hep (where simHit, digiHit are from)
 *      -> generally use dd4hep::rec::Vector3D, convert via ConvertVector() where edm4hep::Vector3d is needed
 * - Indices named i_ ... refer to pixels. Indices named j_ ... refer to in-pixel bins (for charge deposition)
 * - Reference frames: 
 *        - global detector frame, use (x,y,z)
 *             - z along beamline
 *        - local sensor frame: (u,v,w)
 *             - u,v span sensor plane, (for ARCADIA in barrel: v along z)
 *             - w normal to sensor plane
 * - Energies in keV, but deposited energy is always converted to the electron charge equivalent [e-], 3.65 eV per eh-pair
 * - Charges are given as either 
 *        - "raw charge" (as given by Geant4/Allpix2, before thresholding and noise) or 
 *        - "measured charge" (after thresholding and noise)
 */

VTXdigi_Allpix2::VTXdigi_Allpix2(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("SimTrackHitCollectionName", {"UNDEFINED_SimTrackHitCollectionName"}),
                        KeyValues("HeaderName", {"UNDEFINED_HeaderName"}),},
                       {KeyValues("TrackerHitCollectionName", {"UNDEFINED_TrackerHitCollectionName"}),
                        KeyValues("SimTrkHitRelationsCollection", {"UNDEFINED_SimTrkHitRelationsCollection"})}) {
  info() << "Constructed successfully" << endmsg;
}

StatusCode VTXdigi_Allpix2::initialize() {
  info() << "INITIALIZING ..." << endmsg;

  InitServicesAndGeometry();

  InitGaudiProperties();

  InitDetectorGeometry();

  InitLookupTable();
  
  if (m_debugHistograms.value())
    InitHistograms();

  if (!m_debugCsvFileName.value().empty())
    InitCsvOutput();

  /* TODO: check that pixel pitch and count match sensor size in geometry (I am not sure how to get the sensor size from the geometry though, does seem to be a bit more general, subDetectors have children that might or might not be layers). this is currently done for every event, but could be done in the initialization phase, which is of course far more efficient ~ Jona 2025-09
   *  TODO load lookup table for charge transport and diffusion
   *  TODO check that k4geo has same pitch as lookup table */
  info() << " - Initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode VTXdigi_Allpix2::finalize() {
  info() << "FINALIZING VTXdigi_Allpix2..." << endmsg;

  if (m_debugCsvFile.is_open()) {
    m_debugCsvFile.close();
    debug() << " - Closed debug CSV file" << endmsg;
  }

  PrintCounterSummary();

  debug() << " - finalized successfully" << endmsg;
  return StatusCode::SUCCESS;
} 


/* -- Event loop -- */

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> VTXdigi_Allpix2::operator()
  (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {
  /* Initial check: returns false if event has no simHits, throws if there is a problem with the setup */
  if (!CheckEventSetup(simHits, headers))
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());

  /* Output collections */
  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitsLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();

  const u_int32_t rngSeed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  TRandom2 rngEngine = TRandom2(rngSeed);
  verbose() << " - RNG engine initialized, using seed " << rngSeed << endmsg;

  /* Loop over sim hits, digitize them, create output collections */
  for (const auto& simHit : simHits) {
    ++m_counter_simHitsRead;

    /* check if the simHit is on a relevant layer. Needs to be done here ( before GatherHitInfo() ), to avoid out-of-bounds access in Gaudi properties */
    if (!CheckLayerCut(simHit)) continue;

    HitInfo hitInfo; // contains all relevant quantities of the simHit
    HitPosition hitPos; // contains all relevant positions related to the simHit
    std::tie(hitInfo, hitPos) = GatherHitInfoAndPositions(simHit, headers);

    { /* debug statements */
      debug() << "   - Processing simHit (event " << hitInfo.eventNumber() << ", layer " << hitInfo.layerIndex() << ", cellID " << hitInfo.cellID() << " (" << std::bitset<24>(hitInfo.cellID()) << "))" << endmsg;
      debug() << "   - Momentum = " << hitInfo.simMomentum() << " GeV/c, dep. charge = " << hitInfo.charge() << " e-, path length = " << hitPos.path.r() << " mm, Geant4 path length = " << hitInfo.simPathLength() << " mm" << endmsg;
      verbose() << "   - Position (global) " << hitPos.global.x() << " mm, " << hitPos.global.y() << " mm, " << hitPos.global.z() << " mm" << endmsg;
      verbose() << "   - Position (local) " << hitPos.local.x() << " mm, " << hitPos.local.y() << " mm, " << hitPos.local.z() << " mm (in sensor frame)" << endmsg;
      verbose() << "   - Distance to surface " << hitInfo.simSurface()->distance(dd4hep::mm * hitPos.global) << " mm" << endmsg; // dd4hep expects cm, simHitGlobalPosition is in mm. dd4hep::cm=0.1
      verbose() << "   - Entry point (sensor local) : " << hitPos.entry.x() << " mm, " << hitPos.entry.y() << " mm, " << hitPos.entry.z() << " mm" << endmsg;
    }

    if (!CheckSimHitCuts(hitInfo, hitPos))
      continue;
    ++m_counter_simHitsAccepted;
    
    /* Loop over path segments, share each segments charge to pixels */
    PixelChargeMatrix pixelChargeMatrix = DepositAndCollectCharge(hitInfo, hitPos);
    
    /* Generate noise for each pixel (only now, after we know which pixels are fired) */
    pixelChargeMatrix.GenerateNoise(rngEngine, m_electronicNoise); // generate noise only after we know how large the matrix is in the end

    /* find pixels with charge above threshold, create digiHits */
    AnalyseSharedCharge(hitInfo, hitPos, pixelChargeMatrix, simHit, digiHits, digiHitsLinks);

    if (m_debugHistograms)
      FillHistograms_PerSimHit(hitInfo, hitPos, pixelChargeMatrix);
      
    if (m_debugCsv) {
      int i_u_debug, i_v_debug;
      std::tie(i_u_debug, i_v_debug) = ComputePixelIndices(hitPos.local, m_sensorLength.at(0), m_sensorLength.at(1));
      AppendSimHitToCsv(hitInfo, hitPos, i_u_debug, i_v_debug);
    }
  } // loop over sim hits
  
  debug() << "FINISHED event." << endmsg;
  verbose() << " - Returning collections: digiHits.size()=" << digiHits.size() << ", digiHitsLinks.size()=" << digiHitsLinks.size() << endmsg;
  
  return std::make_tuple(std::move(digiHits), std::move(digiHitsLinks));
} // event loop


/* -- Initialization / finalization functions -- */

void VTXdigi_Allpix2::InitServicesAndGeometry() {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc)
    throw GaudiException("Unable to get UniqueIDGenSvc", "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_geometryService = serviceLocator()->service(m_geometryServiceName);
  if (!m_geometryService)
    throw GaudiException("Unable to retrieve the GeoSvc", "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  std::string cellIDstr = m_geometryService->constantAsString(m_encodingStringVariable.value());
  m_cellIDdecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDstr);
  if (!m_cellIDdecoder)
    throw GaudiException("Unable to retrieve the cellID decoder", "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  m_detector = m_geometryService->getDetector();
  if (!m_detector)
    throw GaudiException("Unable to retrieve the DD4hep detector from GeoSvc", "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  const dd4hep::rec::SurfaceManager* simSurfaceManager = m_detector->extension<dd4hep::rec::SurfaceManager>();
  if (!simSurfaceManager)
    throw GaudiException("Unable to retrieve the SurfaceManager from the DD4hep detector", "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  m_simSurfaceMap = simSurfaceManager->map(m_subDetName.value());
  if (!m_simSurfaceMap)
    throw GaudiException("Unable to retrieve the simSurface map for subdetector " + m_subDetName.value(), "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_volumeManager = m_detector->volumeManager();
  if (!m_volumeManager.isValid())
    throw GaudiException("Unable to retrieve the VolumeManager from the DD4hep detector", "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);

  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);

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
    throw GaudiException("Unable to retrieve the subdetector DetElement " + m_subDetName.value(), "VTXdigi_Allpix2::InitServicesAndGeometry()", StatusCode::FAILURE);

  debug() << " - Successfully retrieved all necessary services and detector elements, starting to check Gaudi properties." << endmsg;
}

void VTXdigi_Allpix2::InitGaudiProperties() {
  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "VTXdigi_Allpix2::InitGaudiProperties()", StatusCode::FAILURE);

  if (m_pixelCount.size() != 2)
    throw GaudiException("Property PixelCount expects 2 entries for pix count in u and v (eg. [1024,64])", "VTXdigi_Allpix2::InitGaudiProperties()", StatusCode::FAILURE);

  if (m_maxClusterSize.value().size() != 2) 
    throw GaudiException("Property MaximumClusterSize must have exactly 2 entries (u and v).", "VTXdigi_Allpix2::InitGaudiProperties()", StatusCode::FAILURE);
  for (int i = 0; i < 2; ++i) {
    if (m_maxClusterSize.value().at(i) < 1) 
      throw GaudiException("Property MaximumClusterSize entries must be positive odd numbers.", "VTXdigi_Allpix2::InitGaudiProperties()", StatusCode::FAILURE);
    if ((m_maxClusterSize.value().at(i)-1) % 2 != 0) 
      throw GaudiException("Property MaximumClusterSize entries must be positive odd numbers.", "VTXdigi_Allpix2::InitGaudiProperties()", StatusCode::FAILURE);
  }
}

void VTXdigi_Allpix2::InitDetectorGeometry() {
  /* require members to be set:
   *  - m_sensorThickness
   *  - Gaudi-property: m_pixelCount
   * Sets members:
   *  - m_pixelPitch
   *  - m_sensorThickness
   *  - m_pixelCount
   *  - m_layerCount
   *  - m_layerToIndex
   *  - m_sensorLength (TODO)
   *  - */

  /* Find the relevant layers */
  if (!m_layersToDigitize.value().empty()) {
    /* relevant layers are specified as Gaudi property */

    m_layerCount = m_layersToDigitize.size();
    debug() << " - Digitizing " << m_layerCount << " layers as specified in LayersToDigitize property." << endmsg;

    m_layerToIndex.clear();
    for (size_t layerIndex = 0; layerIndex < m_layersToDigitize.size(); layerIndex++) {
      int layer = m_layersToDigitize.value().at(layerIndex);
      m_layerToIndex.insert({layer, static_cast<int>(layerIndex)});
    }

    std::string layerListStr;
    for (const auto& [layer, index] : m_layerToIndex)
      layerListStr += std::to_string(layer) + " (" + std::to_string(index) + "), ";
    debug() << "   - Layers to be digitized [layer (index)]: " << layerListStr << endmsg;

  }
  else { 
    debug() << " - No layers specified to be digitized, will digitize all layers found in geometry for this subDetector. This will fail if the segmentation varies across layers." << endmsg;

    m_layerToIndex.clear();
    m_layerCount = 0;

    for (const auto& [layerName, layer] : m_subDetector.children()) {
      ++m_layerCount;

      dd4hep::VolumeID layerVolumeID = layer.volumeID();
      int layerNumber = m_cellIDdecoder->get(layerVolumeID, "layer");

      m_layerToIndex.insert({layerNumber, m_layerCount-1});
    }

    std::string layerListStr;
    for (const auto& [layer, index] : m_layerToIndex)
      layerListStr += std::to_string(layer) + " (" + std::to_string(index) + "), ";

    debug() << "   - Found " << m_layerCount << " layers in subDetector " << m_subDetName.value() << "; [layer (index)]: " << layerListStr <<endmsg;
  }

  /* Get pixel pitch from the segmentation
   * The readout is also given by the collection name
   * -> loop over readout, find the one matching our collection
   * -> retrieve pixel pitch from segmentation */
  debug() << " - Retrieving sensor pixel count from segmentation for subDetector \"" << m_subDetName.value() << "\"..." << endmsg;
  
  std::string simHitCollectionName;
  if (this->getProperty("SimTrackHitCollectionName", simHitCollectionName).isFailure())
    throw GaudiException("Could not retrieve SimTrackHitCollectionName property while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
  verbose() << "   - retrieved property SimTrackHitCollectionName \"" << simHitCollectionName << "\". Looking for matching readouts..." << endmsg;

  dd4hep::Detector::HandleMap readoutHandleMap = m_detector->readouts();
  int readoutCount = 0;
  std::string matchedReadoutKey;
  for (const auto& [readoutKey, readoutHandle] : readoutHandleMap) {
    if (simHitCollectionName.find(readoutKey) != std::string::npos && readoutCount == 0) {
      verbose() << "     - Readout \"" << readoutKey << "\" MATCHES the SimTrackHitCollectionName \"" << simHitCollectionName << "\"." << endmsg;
      ++readoutCount;
      if (readoutCount != 1) {
        warning() << "Found multiple (" << readoutCount << ") readouts matching SimTrackHitCollectionName \"" << simHitCollectionName << "\" in detector while checking geometry consistency. Used the first one found. Enable verbose messages for more info." << endmsg;
        continue;
      }
      matchedReadoutKey = readoutKey;

      const dd4hep::Segmentation& segmentation = m_detector->readout(readoutKey).segmentation();
      if (!segmentation.isValid())
        throw GaudiException("Segmentation for readout " + readoutKey + " is not valid while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
      const auto cellDimensions = segmentation.cellDimensions(0); // this assumes all cells have the same dimensions (ie. only one sensor type in this readout)
      m_pixelPitch.at(0) = cellDimensions.at(0) * 10; // convert cm to mm
      m_pixelPitch.at(1) = cellDimensions.at(1) * 10;
      verbose() << "       - Cell dimensions for cellID 0: (" << m_pixelPitch.at(0) << " x " << m_pixelPitch.at(1) << ") mm2" << endmsg;
      /* TODO: get pixel pitch (and count) from segmentation*/
    }
    else {
      verbose() << "     - Readout \"" << readoutKey << "\" does NOT MATCH SimTrackHitCollectionName \"" << simHitCollectionName << "\". skipping." << endmsg;
    }
  } // end loop over readouts
  if (readoutCount == 0)
    throw GaudiException("Could not find any readout matching SimTrackHitCollectionName " + simHitCollectionName + " in detector while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
    
  /* Get sensor length in u,v and active thickness from a random (ie. the first) sensor in the geometry */
  debug() << " - Retrieving sensor dimensions from geometry for subDetector \"" << m_subDetName.value() << "\"..." << endmsg;
  const auto& firstLayer = m_subDetector.children().begin()->second;
  if (!firstLayer)
    throw GaudiException("No layers found in subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
  if (firstLayer.children().empty())
    throw GaudiException("No modules found in first layer of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
  const auto& firstModule = firstLayer.children().begin()->second;

  if (!firstModule)
    throw GaudiException("No modules found in first layer of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
  if (firstModule.children().empty())
    throw GaudiException("No sensors found in first module of first layer of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
  const auto& firstSensor = firstModule.children().begin()->second;
  if (!firstSensor)
    throw GaudiException("No sensors found in first module of first layer of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);

  const dd4hep::VolumeID firstSensorVolumeID = firstSensor.volumeID();
  
  const auto firstSensorItSimSurface = m_simSurfaceMap->find(firstSensorVolumeID);
  if (firstSensorItSimSurface == m_simSurfaceMap->end())
    throw GaudiException("Could not find SimSurface for first sensor (volumeID " + std::to_string(firstSensorVolumeID) + ") in subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
  dd4hep::rec::ISurface* firstSimSurface = firstSensorItSimSurface->second;
  if (!firstSimSurface)
    throw GaudiException("SimSurface pointer for first sensor (volumeID " + std::to_string(firstSensorVolumeID) + ") in subDetector " + m_subDetName.value() + " is null while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);

  float sensorInnerThickness = firstSimSurface->innerThickness() * 10; // in mm
  float sensorOuterThickness = firstSimSurface->outerThickness() * 10; // in mm
  m_sensorThickness = sensorInnerThickness + sensorOuterThickness;

  m_sensorLength.at(0) = firstSimSurface->length_along_u() * 10; // convert to mm
  m_sensorLength.at(1) = firstSimSurface->length_along_v() * 10;

  float pixelCountU = m_sensorLength.at(0) / m_pixelPitch.at(0);
  float pixelCountV = m_sensorLength.at(1) / m_pixelPitch.at(1);
  if (abs(pixelCountU - std::round(pixelCountU)) > 0.001 ||
      abs(pixelCountV - std::round(pixelCountV)) > 0.001) {
    throw GaudiException("Sensor side length (" + std::to_string(m_sensorLength.at(0)) + " x " + std::to_string(m_sensorLength.at(1)) + ") mm and pixel pitch (" + std::to_string(m_pixelPitch.at(0)) + " x " + std::to_string(m_pixelPitch.at(1)) + ") mm result in a non-integer pixel count (" + std::to_string(pixelCountU) + " x " + std::to_string(pixelCountV) + ") in subDetector " + m_subDetName.value() + ".", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
  }
  m_pixelCount.at(0) = std::round(pixelCountU);
  m_pixelCount.at(1) = std::round(pixelCountV);

  debug() << " - Found sensor parameters: area (" << m_sensorLength.at(0) << " x " << m_sensorLength.at(1) << ") mm, thickness (" << sensorInnerThickness << " + " << sensorOuterThickness << ") = " << m_sensorThickness << " mm, pixel pitch (" << m_pixelPitch.at(0) << " x " << m_pixelPitch.at(1) << ") mm, pixel count (" << m_pixelCount.at(0) << " x " << m_pixelCount.at(1) << "). Start looping through layers and checking for consistency in sensor geometry..." << endmsg;

  /* Check how the local sensor coordinates (u,v,w) are defined wrt. the global detector coordinates (x,y,z)
  * This is needed st. the B-field that is applied in the Allpix2 simulation for the LUT goes in the correct direction. 
  * Note: The local coordinate system in DD4hep is different to the one we use.
  *       The transformation matrix that dd4hep provides needs to be adjusted accordingly. */
 const float parallelTolerance = 0.2; // tolerance to consider two vectors parallel (in gradient). Necessary for tilted sensors.
 
 const dd4hep::rec::Vector3D globalSensorPos = TransformLocalToGlobal(dd4hep::rec::Vector3D(0.,0.,0.), firstSensorVolumeID);
 
 const dd4hep::rec::Vector3D globalSensorPosU = TransformLocalToGlobal(dd4hep::rec::Vector3D(1.,0.,0.), firstSensorVolumeID);
 const dd4hep::rec::Vector3D globalSensorUDir = globalSensorPosU - globalSensorPos;
 
 const dd4hep::rec::Vector3D globalSensorPosV = TransformLocalToGlobal(dd4hep::rec::Vector3D(0.,1.,0.), firstSensorVolumeID);
 const dd4hep::rec::Vector3D globalSensorVDir = globalSensorPosV - globalSensorPos;

 const dd4hep::rec::Vector3D globalSensorPosW = TransformLocalToGlobal(dd4hep::rec::Vector3D(0.,0.,1.), firstSensorVolumeID);
 const dd4hep::rec::Vector3D globalSensorWDir = globalSensorPosW - globalSensorPos; // sensor normal vector in global coordinates

  if (simHitCollectionName.find("barrel") != std::string::npos || 
  simHitCollectionName.find("Barrel") != std::string::npos ||
  simHitCollectionName.find("BARREL") != std::string::npos) {
    /* For a barrel sensor, we expect:
     *  - u lies in x-y plane, pointing in positive r-phi direction
     *  - v goes parallel to z (also in the same direct)
     *  - w is sensor normal, points outward (from IP) */
    debug() << " - Detected barrel sensors. Checking transformation from local to global coordinates." << endmsg;
    verbose() << "   - This was inferred by matching \"barrel\" to SimTrackHitCollectionName \"" << simHitCollectionName << "\".  First sensor has global position (" << globalSensorPos.x() << ", " << globalSensorPos.y() << ", " << globalSensorPos.z() << ") mm." << endmsg;
    
    /* Check u plane (expected to lie in x-y plane) */ 
    if ( (globalSensorUDir.z()*globalSensorUDir.z()) > parallelTolerance*parallelTolerance*
         (globalSensorUDir.x()*globalSensorUDir.x() + globalSensorUDir.y()*globalSensorUDir.y()) ) {
      warning() << "Local sensor u-direction (1,0,0) has a z component in global coordinates (" << globalSensorUDir.x() << ", " << globalSensorUDir.y() << ", " << globalSensorUDir.z() << ") mm. Thus it does not lie in x-y plane as is expected, likely not matching the definition in Allpix2. This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* Check u sign (expected to point in positive r-phi direction) */
    const float phiSensor = std::atan2(globalSensorPos.y(), globalSensorPos.x());
    const float phiUDir = std::atan2(globalSensorUDir.y(), globalSensorUDir.x());
    const float dphi = phiUDir - phiSensor; // angle between position vector and u-direction vector in global coordinates. if this is in [0, pi], a positive change in u corresponds to a positive rotation in phi
    if (dphi < 0 || dphi > 3.14159265) {
      warning() << "Sensor u-direction points opposite to the global r-phi direction for barrel sensor at position (" << globalSensorPos.x() << ", " << globalSensorPos.y() << ", " << globalSensorPos.z() << ") mm. The transformation from the global detector coordinates to local sensor does (quite likely) not match the definition in Allpix2. This means the matrix might be rotated or mirrored, and thus the B-field in the simulation might not be applied in the correct direction." << endmsg;
    }

    /* Check v axis (expected to be parallel to z) */
    if ( (globalSensorVDir.x()*globalSensorVDir.x() + globalSensorVDir.y()*globalSensorVDir.y()) > parallelTolerance*parallelTolerance*
         globalSensorVDir.z()*globalSensorVDir.z() ){
      warning() << "Local sensor v-direction (0,1,0) has x- or y-components in global coordinates (" << globalSensorVDir.x() << ", " << globalSensorVDir.y() << ", " << globalSensorVDir.z() << "), and is not parallel to z as is expected.  This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* Check v sign (expected to point in positive z-direction) */
    if (globalSensorVDir.z() < 0) {      
      warning() << "Local sensor v-direction (0,1,0) points opposite to the global z-axis for barrel sensor at position (" << globalSensorPos.x() << ", " << globalSensorPos.y() << ", " << globalSensorPos.z() << ") mm. The transformation from the global detector coordinates to local sensor does (quite likely) not match the definition in Allpix2. This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* Check w plane (expected lie in xy-plane) */
    if (abs(globalSensorWDir.z())/sqrt(globalSensorWDir.x()*globalSensorWDir.x() + globalSensorWDir.y()*globalSensorWDir.y()) > parallelTolerance) {
      warning() << "Local sensor normal (0,0,1) has a z component in global coordinates (" << globalSensorWDir.x() << ", " << globalSensorWDir.y() << ", " << globalSensorWDir.z() << ") mm. Thus it does not lie in x-y plane as is expected, likely not matching the definition in Allpix2. This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* Check w sign (expected to point away from z-axis) */
    if (globalSensorPos.x()*globalSensorWDir.x() + globalSensorPos.y()*globalSensorWDir.y() < 0.) {
      warning() << "Sensor normal vector points inward towards the z-axis for barrel sensor at position (" << globalSensorPos.x() << ", " << globalSensorPos.y() << ", " << globalSensorPos.z() << ") mm. The transformation from the global detector coordinates to local sensor does (quite likely) not match the definition in Allpix2. This means the matrix might be rotated or mirrored, and thus the B-field in the simulation might not be applied in the correct direction." << endmsg;
    }

    /* Check that w has a bigger radial component than u (ie. u points more outwards)
     * Idea: The one with the greater dot product with SensorPos must have a greater radial component */
    const float dotProductU = globalSensorPos.dot(globalSensorUDir);
    const float dotProductW = globalSensorPos.dot(globalSensorWDir);
    if (dotProductW < dotProductU) {
      warning() << "Global sensor normal vector (" << globalSensorWDir.x() << ", " << globalSensorWDir.y() << ", " << globalSensorWDir.z() << ") has a smaller radial component than the u-direction (" << globalSensorUDir.x() << ", " << globalSensorUDir.y() << ", " << globalSensorUDir.z() << "). The transformation from the global detector coordinates to local sensor does (quite likely) not match the definition in Allpix2. This means the matrix might be rotated or mirrored, and thus the B-field in the simulation might not be applied in the correct direction." << endmsg;
    }
    /* TODO: I am sure there is a better way to do these checks, please implement it if you have ideas / find problems. ~ Jona, 2025-11*/
    
    debug() << "   - Barrel sensor coordinate system check completed, local coordinates comply with Allpix2 assumptions." << endmsg;
  }
  if (simHitCollectionName.find("endcap") != std::string::npos || 
           simHitCollectionName.find("Endcap") != std::string::npos ||
           simHitCollectionName.find("ENDCAP") != std::string::npos ||
           simHitCollectionName.find("disk") != std::string::npos ||
           simHitCollectionName.find("Disk") != std::string::npos ||
           simHitCollectionName.find("DISK") != std::string::npos) {
    /* For an endcap sensor, we expect:
     *  - u lies in x-y plane, pointing in positive r direction
     *  - v also lies in x-y planegoes along positive r-phi direction 
     *  - w is sensor normal, points outward (along z) */
    debug() << " - Detected endcap sensors. Checking transformation from local to global coordinates." << endmsg;
    verbose() << "   - This was inferred by matching \"endcap\" or \"disk\" to SimTrackHitCollectionName \"" << simHitCollectionName << "\".  First sensor has global position (" << globalSensorPos.x() << ", " << globalSensorPos.y() << ", " << globalSensorPos.z() << ") mm." << endmsg;
    
    /* Check u plane (expected to lie in x-y plane) */ 
    if (globalSensorUDir.z()*globalSensorUDir.z() > parallelTolerance*parallelTolerance*
         (globalSensorUDir.x()*globalSensorUDir.x() + globalSensorUDir.y()*globalSensorUDir.y()) ) {
      warning() << "Local sensor u-direction (1,0,0) has a z component in global coordinates (" << globalSensorUDir.x() << ", " << globalSensorUDir.y() << ", " << globalSensorUDir.z() << ") mm. Thus it does not lie in x-y plane as is expected, likely not matching the definition in Allpix2. This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* Check v plane (expected to lie in x-y plane) */
    if (globalSensorVDir.z()*globalSensorVDir.z() > parallelTolerance*parallelTolerance*
         (globalSensorVDir.x()*globalSensorVDir.x() + globalSensorVDir.y()*globalSensorVDir.y()) ) {
      warning() << "Local sensor v-direction (0,1,0) has a z component in global coordinates (" << globalSensorVDir.x() << ", " << globalSensorVDir.y() << ", " << globalSensorVDir.z() << ") mm. Thus it does not lie in x-y plane as is expected, likely not matching the definition in Allpix2. This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* Check u sign (expected to point away from z-axis) */
    if (globalSensorPos.x()*globalSensorUDir.x() + globalSensorPos.y()*globalSensorUDir.y() < 0.) {
      warning() << "Sensor u-direction points inward towards the z-axis for endcap sensor at position (" << globalSensorPos.x() << ", " << globalSensorPos.y() << ", " << globalSensorPos.z() << ") mm. The transformation from the global detector coordinates to local sensor does (quite likely) not match the definition in Allpix2. This means the matrix might be rotated or mirrored, and thus the B-field in the simulation might not be applied in the correct direction." << endmsg;
    }

    /* Check that u has a bigger radial component than v (ie. u points more outwards)
     * Idea: The one with the greater dot product with SensorPos must have a greater radial component */
    const float dotProductU = globalSensorPos.dot(globalSensorUDir);
    const float dotProductV = globalSensorPos.dot(globalSensorVDir);
    if (dotProductV > dotProductU) {
      warning() << "Global sensor normal vector (" << globalSensorVDir.x() << ", " << globalSensorVDir.y() << ", " << globalSensorVDir.z() << ") has a smaller radial component than the u-direction (" << globalSensorUDir.x() << ", " << globalSensorUDir.y() << ", " << globalSensorUDir.z() << "). The transformation from the global detector coordinates to local sensor does (quite likely) not match the definition in Allpix2. This means the matrix might be rotated or mirrored, and thus the B-field in the simulation might not be applied in the correct direction." << endmsg;
    }

    /* Check w axis (expected to be parallel to z) */
    if ( sqrt(globalSensorWDir.x()*globalSensorWDir.x() + globalSensorWDir.y()*globalSensorWDir.y()) / abs(globalSensorWDir.z()) > parallelTolerance ) {
      warning() << "Local sensor normal (0,0,1) has x- or y-components in global coordinates (" << globalSensorWDir.x() << ", " << globalSensorWDir.y() << ", " << globalSensorWDir.z() << "), and is not parallel to z as is expected.  This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* Check w sign (expected to point in positive z-direction) */
    if (globalSensorWDir.z() < 0) {      
      warning() << "Local sensor normal (0,0,1) points opposite to the global z-axis for endcap sensor at position (" << globalSensorPos.x() << ", " << globalSensorPos.y() << ", " << globalSensorPos.z() << ") mm. The transformation from the global detector coordinates to local sensor does (quite likely) not match the definition in Allpix2. This means the matrix might be rotated or mirrored and the B-field simulated in AP2 might be applied in a wrong direction." << endmsg;
    }

    /* TODO: I am sure there is a better way to do these checks, please implement it if you have ideas / find problems. ~ Jona, 2025-11*/
    
    debug() << "   - Barrel sensor coordinate system check completed, local coordinates comply with Allpix2 assumptions." << endmsg;

  }
  else {
    warning() << " - Could not determine if sensors are in barrel or endcap from SimTrackHitCollectionName \"" << simHitCollectionName << "\" by matching it to \"barrel\", \"endcap\", or \"disk\". The transformation from the global detector coordinates to local sensor coordinates is not guaranteed to match the definition in Allpix2. This means the matrix might be rotated or mirrored, and thus the B-field in the simulation might not be applied in the correct direction." << endmsg;
  }

  /* check that every sensor in the relevant layers matches the dimensions of the first one we found.*/
  debug() << " - Looping over layers, modules and sensors in subDetector \"" << m_subDetName.value() << "\" to check for consistent sensor dimensions..." << endmsg;
  for (const auto& [layerName, layer] : m_subDetector.children()) {
    dd4hep::VolumeID layerVolumeID = layer.volumeID();
    int layerNumber = m_cellIDdecoder->get(layerVolumeID, "layer");
    
    const int nModules = layer.children().size();
    if (m_layerToIndex.find(layerNumber) == m_layerToIndex.end()) {
      debug() << "   - Skipping layer " << layerName << " (layerNumber " << layerNumber << ", volumeID " << layerVolumeID << " with " << nModules << " modules) as it is not in the LayersToDigitize list." << endmsg;
      continue;
    }
    else {
      debug() << "   - Found layer \"" << layerName << "\" (layerNumber " << layerNumber << ", volumeID " << layerVolumeID << ", " << nModules << " modules) for subDetector \"" << m_subDetName.value() << "\"." << endmsg;
    }

    /* loop over modules and sensors, check dimensions for each. 
     * this takes far less than a second for the complete IDEA vertex detector ~ Jona, 2025-11 */
    int moduleN = 0, sensorN = 0;
    for (const auto& [moduleName, module] : layer.children()) {
      ++moduleN;
      
      for (const auto& [sensorName, sensor] : module.children()) {
        dd4hep::VolumeID sensorVolumeID = sensor.volumeID();
        ++sensorN;

        for (const auto& [pixelName, pixel]: sensor.children()) {
          // nothing to do here yet
          verbose() << "     - Found pixel \"" << pixelName << "\" in sensor \"" << sensorName << "\"." << endmsg;
        }

        const auto itSimSurface = m_simSurfaceMap->find(sensorVolumeID);
        if (itSimSurface == m_simSurfaceMap->end())
          throw GaudiException("Could not find SimSurface for sensor " + sensorName + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layerNumber) + " of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
        dd4hep::rec::ISurface* simSurface = itSimSurface->second;
        if (!simSurface)
          throw GaudiException("SimSurface pointer for sensor " + sensorName + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layerNumber) + " of subDetector " + m_subDetName.value() + " is null while checking geometry consistency.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);

        const float sensorLength_u = simSurface->length_along_u() * 10; // convert to mm
        const float sensorLength_v = simSurface->length_along_v() * 10; 
        const float sensorThickness = (simSurface->innerThickness() + simSurface->outerThickness()) * 10; // in mm

        /* sanity check: does this sensor have the same dimensions as the one we looked at in the beginning? */
        if (std::abs(sensorLength_u - m_sensorLength.at(0)) > 0.001 ||
            std::abs(sensorLength_v - m_sensorLength.at(1)) > 0.001 ||
            std::abs(sensorThickness - m_sensorThickness) > 0.001) {
          throw GaudiException("Sensor dimension mismatch found in sensor " + sensorName + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layerNumber) + " of subDetector " + m_subDetName.value() + ": expected dimensions of (" + std::to_string(m_sensorLength.at(0)) + " x " + std::to_string(m_sensorLength.at(1)) + " x " + std::to_string(m_sensorThickness) + ") mm3, but found (" + std::to_string(sensorLength_u) + " x " + std::to_string(sensorLength_v) + " x " + std::to_string(sensorThickness) + ") mm3. This algorithm expects exactly one type of sensor per subDetector. Use different instances of the algorithm if different layers consist of different sensors.", "VTXdigi_Allpix2::InitDetectorGeometry()", StatusCode::FAILURE);
        }
        else {
          verbose() << "     - Found sensor: " << sensorName << ", volumeID: " << sensorVolumeID << " (sensor " << sensorN << " in layer " << layerNumber << "). Dimensions are consistent." << endmsg;
        }
      } // loop over sensors per module
    } // loop over modules per layer
  } // loop over layers

  info() << " - Retrieved sensor geometry parameters for subDetector \"" << m_subDetName.value() << "\":" << endmsg;
  info() << "    - Pixel pitch (" << m_pixelPitch.at(0) << " x " << m_pixelPitch.at(1) << ") mm and count (" << m_pixelCount.at(0) << " x " << m_pixelCount.at(1) << ")" << endmsg;
  info() << "    - Sensor active thickness: " << m_sensorThickness << " mm and active area (" << (m_pixelPitch.at(0) * m_pixelCount.at(0)) << " x " << (m_pixelPitch.at(1) * m_pixelCount.at(1)) << ") mm2" << endmsg;

  std::string layerListStr = "";
  for (const auto& [layer, index] : m_layerToIndex)
    layerListStr += std::to_string(layer) + ", ";
  info() << " - Relevant layers: " + layerListStr + " (All sensors in these layers have been checked to be consistent with the retrieved geometry)" << endmsg;
}

void VTXdigi_Allpix2::InitLookupTable() {
  /* requires members to be set:
   *  - m_pixelPitch
   *  - m_sensorThickness
   * Checks the consistency of the members:
   *  - Gaudi-property: m_GlobalSharingMatrix *or* m_MatrixFileName
   * Sets members:
   *  - m_matrixSize
   *  - m_inPixelBinCount
   *  - m_LUT */
  debug() << " - Importing lookup table..." << endmsg;
  
  /* TODO: check if member m_matrixSize is even needed. It is also saved in m_LUT*/

  /* sanity check */
  if (m_globalSharingMatrix.value().empty() && m_LUTFileName.value().empty())
    throw GaudiException("No lookup table specified! Please provide either a global matrix or a LUT file.", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);
  if (m_globalSharingMatrix.value().size() > 0 && m_LUTFileName.value().length() > 0)
    throw GaudiException("Please provide either a global matrix or a matrix file, not both!", "VTXdigi_Allpix2::InitGaudiProperties()", StatusCode::FAILURE);

  /* load matrices */
  if (!m_globalSharingMatrix.value().empty()) { // use global matrix (mostly for debugging)
    debug() << "   - Global matrix specified in <config>, applying this for all layers" << endmsg;

    for (int i=0; i<3; i++)
      m_inPixelBinCount[i] = 10;

    m_matrixSize = static_cast<int>(std::sqrt(m_globalSharingMatrix.value().size()));
    verbose() << "   - Using matrix size of " << m_matrixSize << " (from global matrix with " << m_globalSharingMatrix.value().size() << " entries)" << endmsg;

    m_LUT = std::make_unique<LookupTable>(m_inPixelBinCount[0], m_inPixelBinCount[1], m_inPixelBinCount[2], m_matrixSize);

    /* matrix size is and validity of weights checked in LookupTable::SetAllMatrices() */
    m_LUT->SetAllMatrices(m_globalSharingMatrix.value());
  }
  else { // load from file
    info() << " - Loading matrices from file: \"" << m_LUTFileName.value() << "\". Expecting Allpix2 format, defined in ChargePropagationWriter module." << endmsg;
    /* This implements parsing the default Allpix2 LUT file format. A general version was implemented in a previous commit (2025-11), but removed because the amount of options made it unneccessarily hard to validate and use. 
     * See https://indico.cern.ch/event/1489052/contributions/6475539/attachments/3063712/5418424/Allpix_workshop_Lemoine.pdf (slide 10) for more info on fields in the LUT file */

    const int headerLines = 5; // allpix2 LUT files have 5 header lines, and then one matrix per line

    debug() << "   - Opening LUT file: \"" << m_LUTFileName.value() << "\"." << endmsg;
    std::ifstream LUTFile(m_LUTFileName.value());
    if (!LUTFile.is_open())
      throw GaudiException("Could not open LUT file: " + m_LUTFileName.value(), "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);

    std::string line;
    int lineNumber = 0; // next line to be read (0-indexed)
    float matricesWeightSum = 0.f;

    /* loading pixel-pitch, thickness and in-pixel bin counts from header (all in 5th line) */
    for (; lineNumber < 5; ++lineNumber)
      std::getline(LUTFile, line);
    std::istringstream headerStringStream(line);
    std::string headerEntry;
    std::vector<std::string> headerLineEntries;

    while (std::getline(headerStringStream, headerEntry, ' ')) {
      if (!headerEntry.empty())
        headerLineEntries.push_back(headerEntry);
    }

    if (headerLineEntries.size() != 11)
      throw GaudiException("Invalid number of entries in LUT file at header line 5: found " + std::to_string(headerLineEntries.size()) + " entries, expected 11.", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);

    for (int i=0; i<3; i++) {
      m_inPixelBinCount.at(i) = std::stoi(headerLineEntries.at(7+i));
    }
    debug() << "   - found in-pixel bin count of (" << m_inPixelBinCount.at(0) << ", " << m_inPixelBinCount.at(1) << ", " << m_inPixelBinCount.at(2) << ") from LUT file header." << endmsg;

    /* Sensor pitch & thickness are given in the header. Compare them to values from the detector geometry we retrieved in InitDetectorGeometry(). */
    const float fileSensorThickness = std::stof(headerLineEntries.at(0)) / 1000.f; // convert from um to mm
    const float filePixelPitchU = std::stof(headerLineEntries.at(1)) / 1000.f;
    const float filePixelPitchV = std::stof(headerLineEntries.at(2)) / 1000.f;

    if (std::abs(fileSensorThickness - m_sensorThickness) > 0.001) {
      throw GaudiException("Sensor thickness mismatch between LUT file and detector geometry: LUT file specifies " + std::to_string(fileSensorThickness) + " mm, but geometry has " + std::to_string(m_sensorThickness) + " mm.", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);
    }
    if (std::abs(filePixelPitchU - m_pixelPitch.at(0)) > 0.001 ||
        std::abs(filePixelPitchV - m_pixelPitch.at(1)) > 0.001) {
      throw GaudiException("Pixel pitch mismatch between LUT file and detector geometry: LUT file specifies (" + std::to_string(filePixelPitchU) + " x " + std::to_string(filePixelPitchV) + ") mm2, but geometry has (" + std::to_string(m_pixelPitch.at(0)) + " x " + std::to_string(m_pixelPitch.at(1)) + ") mm2.", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);
    }
    debug() << "   - found pixel-pitch of (" << m_pixelPitch.at(0) << " x " << m_pixelPitch.at(1) << ") mm2 and thickness of " << m_sensorThickness << " mm in the LUT file header. These match the pitch and thickness in the detector geometry." << endmsg;

    /* get the matrix size (5x5, 7x7, ...) from the length of the first line after the header */

    if (std::getline(LUTFile, line)) {
      m_matrixSize = static_cast<int>(std::sqrt(std::count(line.begin(), line.end(), ' ') - 2)); // not very robust, but works for valid Allpix2 files. first 3 entries are bin indices
    }
    else {
      throw GaudiException("Could not read first line after header in LUT file: " + m_LUTFileName.value(), "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);
    }
    debug() << "   - Detected matrix size of " << m_matrixSize << " from first line." << endmsg;

    /* set up mapping from Allpix2 LUT format
    *   (row-major, starts on bottom left)
    * to the format expected by the LookupTable class 
    *   (row-major, starts on top-left) */
    std::vector<int> valueIndices(m_matrixSize*m_matrixSize, 0); // i: index in local format; valueIndices[i]: index in Allpix2 format
    for (int i_u = 0; i_u < m_matrixSize; i_u++) {
      for (int i_v = 0; i_v < m_matrixSize; i_v++) {
        int i_allpix2 = i_u + (m_matrixSize - 1 - i_v) * m_matrixSize;
        int i_local = i_u + i_v * m_matrixSize;
        valueIndices[i_local] = i_allpix2;
      }
    }

    /* reset the file and go to the beginning of the matrix lines (after header) */
    LUTFile.clear();
    LUTFile.seekg(0, std::ios::beg);
    for (int i=0; i<headerLines; ++i)
      std::getline(LUTFile, line);
    
    /* loop over lines that contain a matrix each, set the matrices */
    int matrixSize = 0;
    m_LUT = std::make_unique<LookupTable>(m_inPixelBinCount[0], m_inPixelBinCount[1], m_inPixelBinCount[2], m_matrixSize);

    debug() << "   - Loading matrices from LUT file as lines ..." << endmsg;
    while (std::getline(LUTFile, line)) {
      if (line.empty() || line[0] == '#')
        throw GaudiException("Empty or comment line found in LUT file at line " + std::to_string(lineNumber+1) + ". All lines (after 5 header lines) must contain valid matrix data.", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);
      
      std::istringstream stringStream(line);
      std::vector<std::string> lineEntries;
      std::string entryString;

      /* read the line & do sanity checks */
      while (std::getline(stringStream, entryString, ' ')) {
        if (!entryString.empty())
          lineEntries.push_back(entryString);
      }
      if (matrixSize == 0)
        matrixSize = static_cast<int>(std::sqrt(lineEntries.size() - 3)); // first 3 entries are bin indices
      else if (static_cast<int>(std::sqrt(lineEntries.size() - 3)) != matrixSize)
        throw GaudiException("Invalid number of entries in LUT file at line " + std::to_string(lineNumber+1) + ": found " + std::to_string(lineEntries.size()) + " entries, expected (3 indices + matrixSize^2) = " + std::to_string(3 + m_matrixSize * m_matrixSize) + ".", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);

      /* in-pixel bin that this line corresponds to, given in first 3 columns */
      int j_u = std::stoi(lineEntries[0]) - 1; // Allpix2 input is 1-indexed, sane people use 0-indexing
      int j_v = std::stoi(lineEntries[1]) - 1;
      int j_w = std::stoi(lineEntries[2]) - 1;

      if (j_u < 0 || j_u >= m_inPixelBinCount[0] ||
          j_v < 0 || j_v >= m_inPixelBinCount[1] ||
          j_w < 0 || j_w >= m_inPixelBinCount[2]) {
        throw GaudiException("Invalid in-pixel bin indices in LUT file at line " + std::to_string(lineNumber+1) + ": got (" + std::to_string(j_u) + ", " + std::to_string(j_v) + ", " + std::to_string(j_w) + "), but expected ranges are [0, " + std::to_string(m_inPixelBinCount[0]-1) + "], [0, " + std::to_string(m_inPixelBinCount[1]-1) + "], [0, " + std::to_string(m_inPixelBinCount[2]-1) + "].", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);
      }

      /* Parse matrix values & set it */
      std::vector<float> matrixWeights(m_matrixSize*m_matrixSize, 0.);
      float matrixWeightSum = 0.f;
      for (int i = 0; i < m_matrixSize*m_matrixSize; i++) {
        float entry = std::stof(lineEntries[3 + valueIndices[i]]); // NaN check done on sum
        matrixWeights[i] = entry;
        matrixWeightSum += entry;
      } 
      verbose() << "     - Parsed matrix for in-pixel bin (" << j_u << ", " << j_v << ", " << j_w << ") with entry sum " << std::to_string(matrixWeightSum) << ", setting it now..." << endmsg;
      matricesWeightSum += matrixWeightSum;
      m_LUT->SetMatrix(j_u, j_v, j_w, matrixWeights);

      lineNumber++;
    } // loop over lines containing a matrix each

    if (lineNumber - headerLines != m_inPixelBinCount[0] * m_inPixelBinCount[1] * m_inPixelBinCount[2])
      throw GaudiException("Invalid number of matrices loaded from file: expected " + std::to_string(m_inPixelBinCount[0] * m_inPixelBinCount[1] * m_inPixelBinCount[2]) + " matrices (inferred from InPixelBinCount = [" + std::to_string(m_inPixelBinCount[0]) + ", " + std::to_string(m_inPixelBinCount[1]) + ", " + std::to_string(m_inPixelBinCount[2]) + "]), but got " + std::to_string(lineNumber - headerLines) + ".", "VTXdigi_Allpix2::InitLookupTable()", StatusCode::FAILURE);

    matricesWeightSum /= static_cast<float>(lineNumber - headerLines);
    info() << " -   Loaded " << (lineNumber - headerLines) << " matrices from file. " << matricesWeightSum*100 << " percent of charge deposited in the sensor volume is collected by the pixels (the rest is lost, eg. due to being outside of depletion or due to trapping)." << endmsg;
  } // if load LUT from file
} // InitLookupTable()

void VTXdigi_Allpix2::InitHistograms() {
  warning() << " - You enabled creating debug histograms by setting `DebugHistograms = True`. This is NOT MULTITHREADING SAFE and will cause crashes if multithreading is used." << endmsg;
  verbose () << " - Creating debug histograms ..." << endmsg;

  /* -- Global Histograms (collect from all layers) -- */

  m_histGlobal.at(histGlobal_simHitE).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
      "SimHit/DepositedE",
      "SimHit deposited energy in the sensor;Deposited energy [keV]", 
      {1000, 0.f, m_sensorThickness*2000.f}}); 
  m_histGlobal.at(histGlobal_simHitCharge).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "SimHit/DepositedCharge",
      "SimHit deposited charge in the sensor;Number of dep. charges [e-]",
      {1000, 0.f, m_sensorThickness*500000.f}});

  m_histGlobal.at(histGlobal_clusterSize_raw).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
      "DigiHit_Global/PixelsPerSimHit_Raw", 
      "Number of pixels that receive any charge per simhit);Number of pixels [pix]",
      {50, -0.5f, 49.5f}
    }
  );
  m_histGlobal.at(histGlobal_clusterSize_measured).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "DigiHit_Global/PixelsPerSimHit_Measured",
      "Number of digiHits per simHit (ie. number of pixels above threshold, similar to cluster size);Number of digiHits [pix]",
      {50, -0.5f, 49.5f}
    }
  );

  m_histGlobal.at(histGlobal_pathLength).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
      "SimHit/PathLength",
      "Path length in sensor active volume, as computed by VTXdigi reconstruction;Path length [um]",
      {500, 0.f, m_sensorThickness*10*1000.f}
    }
  ); // in um, max 10x sensor thickness (for very shallow angles)
  m_histGlobal.at(histGlobal_pathLengthGeant4).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "SimHit/PathLength-Geant4",
      "Path length in sensor active volume, as given by Geant4;Path length [um]",
      {500, 0.f, m_sensorThickness*10*1000.f}
    }
  ); // in um, max 10x sensor thickness (for very shallow angles)

  m_histGlobal.at(histGlobal_EntryPointX).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "SimHit/PathBeginning_x",
      "SimHit computed path starting point U (in local sensor frame);u [mm]",
      {400, -4.f, 4.f}
    }
  );
  m_histGlobal.at(histGlobal_EntryPointY).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "SimHit/PathBeginning_y",
      "SimHit Entry Point V (in local sensor frame);v [mm]",
      {2000, -20.f, 20.f}
    }
  );
  m_histGlobal.at(histGlobal_EntryPointZ).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
      "SimHit/PathBeginning_z",
      "SimHit Entry Point W (in local sensor frame);w [mm]",
      {1000, -300.f, 300.f}
    }
  );

  m_histGlobal.at(histGlobal_DisplacementU).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "DigiHit_Global/LocalDisplacementU",
      "Displacement in u (local sensor frame): digiHit_u - simHit_u;Displacement in u [um]",
      {800, -200.f, 200.f}
    }
  );
  m_histGlobal.at(histGlobal_DisplacementV).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "DigiHit_Global/LocalDisplacementV",
      "Displacement in v (local sensor frame): digiHit_v - simHit_v;Displacement in v [um]",
      {800, -200.f, 200.f}
    }
  );
  m_histGlobal.at(histGlobal_DisplacementR).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "DigiHit_Global/GlobalDisplacementR",
      "Displacement distance (global frame): | digiHit - simHit |;Displacement distance [um]",
      {300, 0.f, 300.f}
    }
  );

  m_histGlobal.at(histGlobal_chargeCollectionEfficiency_raw).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "DigiHit_Global/ChargeCollectionEfficiency_rawCharge",
      "Charge collection efficiency (raw charge, eg. before noise & threshold);# e- (digitised) / # e- (simHit)",
      {500, 0.f, 2.f}
    }
  );
  
  m_histGlobal.at(histGlobal_chargeCollectionEfficiency).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "DigiHit_Global/ChargeCollectionEfficiency", "Charge collection efficiency;# e- (digitised) / # e- (simHit)",
      {500, 0.f, 2.f}
    }
  );

  m_histGlobal.at(histGlobal_pixelChargeMatrix_size_u).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "Internal/PixelChargeMatrix_Size_u",
      "Final Size of the Pixel Charge Matrix in u (local);pixel charge matrix size u [pix]",
      {100, -0.5f, 99.5f}
    }
  );
  m_histGlobal.at(histGlobal_pixelChargeMatrix_size_v).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "Internal/PixelChargeMatrix_Size_v",
      "Final Size of the Pixel Charge Matrix in v (local);pixel charge matrix size v [pix]",
      {100, -0.5f, 99.5f}
    }
  );
  m_histGlobal.at(histGlobal_simHitPDG).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
      "SimHit/PDG",
      "PDG number of the SimHit particle;PDG number",
      {1401, -700.5f, 700.5f}
    }
  );



  m_histGlobal2d.at(histGlobal2d_pathLength_vs_G4PathLength).reset(
    new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this,
      "SimHit/PathLength_vs_G4PathLength",
      "Path length in sensor active volume: VTXdigi reconstruction vs Geant4;Path length Geant4 [um];Path length VTXdigi [um]",
      {500, 0.f, m_sensorThickness*10*1000.f},
      {500, 0.f, m_sensorThickness*10*1000.f}
    }
  );

  /* -- Per-layer Histograms -- */
  
  m_hist1d.resize(m_layersToDigitize.size());
  m_hist2d.resize(m_layersToDigitize.size()); // resize vector to hold all histograms
  m_histProfile1d.resize(m_layersToDigitize.size());
  m_histWeighted2d.resize(m_layersToDigitize.size());

  for (int layer : m_layersToDigitize) {
    int layerIndex = m_layerToIndex.at(layer);

    /* -- 1d Histograms -- */

    m_hist1d.at(layerIndex).at(hist1d_DigiHitCharge_raw).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/TotalCharge_Raw",
        "Sum of digiHit charge per simHit (raw, ie. before noise & threshold) - Layer "+std::to_string(layer)+";Sum of digiHit charge [e-]",
        {1000, 0.f, m_sensorThickness*500000.f}
      }
    );
    m_hist1d.at(layerIndex).at(hist1d_DigiHitCharge_measured).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/TotalCharge_Measured",
        "Sum of digiHit charge per simHit (measured, ie. after noise & threshold) - Layer "+std::to_string(layer)+";Sum of digiHit charge [e-]",
        {1000, 0.f, m_sensorThickness*500000.f}
      }
    );

    m_hist1d.at(layerIndex).at(hist1d_ClusterSize_raw).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/PixelsPerSimHit_Raw",
        "Number of pixels that receive any charge per simhit - Layer "+std::to_string(layer)+";Number of pixels [pix]",
        {50, -0.5f, 49.5f}
      }
    );
    m_hist1d.at(layerIndex).at(hist1d_ClusterSize_measured).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/PixelsPerSimHit_Measured",
        "Number of digiHits per simHit (ie. number of pixels above threshold, similar to cluster size) - Layer "+std::to_string(layer)+";Number of digiHits [pix]",
        {50, -0.5f, 49.5f}
      }
    );

    m_hist1d.at(layerIndex).at(hist1d_IncidentAngle_ThetaLocal).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
        "SimHit/IncidentAngle_LocalTheta_Layer"+std::to_string(layer),
        "Local theta: particle momentum incident angle to sensor normal - Layer "+std::to_string(layer)+";Polar angle [deg]",
        {360, 0.f, 180.f}
      }
    );
    m_hist1d.at(layerIndex).at(hist1d_IncidentAngle_PhiLocal).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
        "SimHit/IncidentAngle_LocalPhi_Layer"+std::to_string(layer),
        "Local phi: particle momentum incident angle to sensor normal - Layer "+std::to_string(layer)+";Azimuthal angle [deg]",
        {360, -180.f, 180.f}
      }
    );

    m_hist1d.at(layerIndex).at(hist1d_SimHitMomentum).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this, 
        "SimHit/SimHitMomentum_Layer"+std::to_string(layer),
        "Momentum of the simHit particle at sensor position - Layer "+std::to_string(layer)+";Momentum [MeV/c]",
        {10000, 0.f, 5000.f}
      }
    );

    /* -- Profile Histograms -- */

    m_histProfile1d.at(layerIndex).at(histProfile1d_clusterSize_vs_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
        "DigiHit_Layer"+std::to_string(layer)+"/PixelsPerSimHit_vs_z",
        "Number of digiHits per simhit vs z position of the simHit - Layer "+std::to_string(layer)+";SimHit z [mm];Number of digiHits [pix]",
        {2000, -500.f, 500.f}
      }
    );
    m_histProfile1d.at(layerIndex).at(histProfile1d_clusterSize_vs_module_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
        "DigiHit_Layer"+std::to_string(layer)+"/PixelsPerSimHit_vs_Module_z",
        "Number of digiHits per simhit vs z position of the module - Layer "+std::to_string(layer)+";SimHit z [mm];Number of digiHits [pix]",
        {2000, -500.f, 500.f}
      }
    );
    m_histProfile1d.at(layerIndex).at(histProfile1d_clusterSize_vs_moduleID).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float>{this,
        "DigiHit_Layer"+std::to_string(layer)+"/PixelsPerSimHit_vs_ModuleID",
        "Number of digiHits per simhit vs Module ID - Layer "+std::to_string(layer)+";Module ID;Number of digiHits [pix]",
        {2000, -0.5f, 1999.5f}
      }
    );

    /* -- 2D Histograms -- */

    m_hist2d.at(layerIndex).at(hist2d_hitMap_simHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/HitMap_simHits",
        "SimHit Hitmap - Layer " + std::to_string(layer) + ";u [pix]; v [pix]",
        {static_cast<unsigned int>(m_pixelCount.at(0)), -0.5f, static_cast<float>(m_pixelCount.at(0)+0.5f)},
        {static_cast<unsigned int>(m_pixelCount.at(1)), -0.5f, static_cast<float>(m_pixelCount.at(1)+0.5f)}
      }
    );
    m_hist2d.at(layerIndex).at(hist2d_hitMap_digiHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/HitMap_digiHits",
        "DigiHit Hitmap - Layer " + std::to_string(layer) + ";u [pix]; v [pix]",
        {static_cast<unsigned int>(m_pixelCount.at(0)), -0.5f, static_cast<float>(m_pixelCount.at(0)+0.5f)},
        {static_cast<unsigned int>(m_pixelCount.at(1)), -0.5f, static_cast<float>(m_pixelCount.at(1)+0.5f)}
      }
    );
    m_hist2d.at(layerIndex).at(hist2d_pathLength_vs_simHit_v).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "SimHit/PathLength_vs_z_2D_Layer"+std::to_string(layer),
        "Path length in sensor active volume (as computed here and used for charge sharing) vs. simHit z position (global) - Layer " + std::to_string(layer) + ";SimHit z [mm];Path length [um]",
        {2000, -500.f, 500.f},
        {200, 0., 20*m_sensorThickness*1000} // in um
      }
    );
    m_hist2d.at(layerIndex).at(hist2d_pixelChargeMatrixSize).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/PixelChargeMatrixSize",
        "Pixel Charge Matrix Size - Layer " + std::to_string(layer) + ";u [pix];v [pix]",
        {100, -0.5f, 99.5f},
        {100, -0.5f, 99.5f}
      }
    );
    m_hist2d.at(layerIndex).at(hist2d_IncidentAngle).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "SimHit/IncidentAngle_2D_Layer"+std::to_string(layer),
        "Incident particle angle to sensor normal - Layer " + std::to_string(layer) + ";Local phi [deg];Local theta [deg]",
        {360, -180.f, 180.f},
        {360, 0.f, 180.f}
      }
    );
    
    m_hist2d.at(layerIndex).at(hist2d_clusterSize_vs_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/PixelsPerSimHit_vs_z_2D",
        "Number of digiHits per simhit vs z position of the simHit (global) - Layer " + std::to_string(layer) + ";SimHit z [mm];Number of digiHits [pix]",
        {2000, -500.f, 500.f},
        {20, -0.5f, 19.5f}
      }
    );
    m_hist2d.at(layerIndex).at(hist2d_clusterSize_vs_module_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/PixelsPerSimHit_vs_z_2D",
        "Number of digiHits per simhit vs z position of the module (global) - Layer " + std::to_string(layer) + ";Module z [mm];Number of digiHits [pix]",
        {2000, -500.f, 500.f},
        {20, -0.5f, 19.5f}
      }
    );
    m_hist2d.at(layerIndex).at(hist2d_averageCluster_binary).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this, 
        "DigiHit_Layer"+std::to_string(layer)+"/AverageCluster_binary",
        "Layout of the pixels above threshold per simHit (central pixel is defined by the simHit) - Layer " + std::to_string(layer) +";u [pix]; v [pix]",
        {static_cast<unsigned int>(m_maxClusterSize.value().at(0)),
          -static_cast<float>(m_maxClusterSize.value().at(0))/2,
          static_cast<float>(m_maxClusterSize.value().at(0))/2},
        {static_cast<unsigned int>(m_maxClusterSize.value().at(1)),
          -static_cast<float>(m_maxClusterSize.value().at(1))/2,
          static_cast<float>(m_maxClusterSize.value().at(1))/2}
      }
    );
    m_hist2d.at(layerIndex).at(hist2d_totalCharge_vs_simHitCharge).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this,
        "DigiHit_Layer"+std::to_string(layer)+"/TotalDigiHitCharge_vs_SimHitCharge",
        "Sum of digiHit charge (per simHit) vs simHit deposited charge - Layer " + std::to_string(layer) + ";SimHit deposited charge [e-];Sum of digiHit charge [e-]",
        {1000, 0.f, m_sensorThickness*500000.f},
        {1000, 0.f, m_sensorThickness*500000.f}
      }
    ); 

    /* -- Weighted 2D Histograms -- */

    m_histWeighted2d.at(layerIndex).at(histWeighted2d_averageCluster_analog).reset(
      new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this,
        "DigiHit_Layer"+std::to_string(layer)+"/AverageCluster_analog",
        "Layout of the charge collection per simHit (central pixel is defined by the simHit) - Layer " + std::to_string(layer) +";u [pix]; v [pix]",
        {static_cast<unsigned int>(m_maxClusterSize.value().at(0)),
          -static_cast<float>(m_maxClusterSize.value().at(0))/2,
          static_cast<float>(m_maxClusterSize.value().at(0))/2},
        {static_cast<unsigned int>(m_maxClusterSize.value().at(1)),
          -static_cast<float>(m_maxClusterSize.value().at(1))/2,
          static_cast<float>(m_maxClusterSize.value().at(1))/2}
      }
    );
    m_histWeighted2d.at(layerIndex).at(histWeighted2d_chargeOriginU).reset(
      new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this,
        "DigiHit_Layer"+std::to_string(layer)+"/ChargeOrigin_u",
        "Charge collected from in-pix bins in u - Layer " + std::to_string(layer) + ";u pos. relative to collecting pixel center [pix]; w pos. [um]",
        {static_cast<unsigned int>(m_matrixSize*m_inPixelBinCount[0]),
          -static_cast<float>(m_matrixSize)/2,
          static_cast<float>(m_matrixSize)/2},
        {static_cast<unsigned int>(m_inPixelBinCount[2]),
          -m_sensorThickness*1000.f/2.f, 
          m_sensorThickness*1000.f/2.f}
      }
    );
    m_histWeighted2d.at(layerIndex).at(histWeighted2d_chargeOriginV).reset(
      new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, float>{this,
        "DigiHit_Layer"+std::to_string(layer)+"/ChargeOrigin_v",
        "Charge collected from in-pix bins in v - Layer " + std::to_string(layer) + ";v pos. relative to collecting pixel center [pix]; w pos. [um]",
        {static_cast<unsigned int>(m_matrixSize*m_inPixelBinCount[1]),
          -static_cast<float>(m_matrixSize)/2,
          static_cast<float>(m_matrixSize)/2},
        {static_cast<unsigned int>(m_inPixelBinCount[2]),
          -m_sensorThickness*1000.f/2.f, 
          m_sensorThickness*1000.f/2.f}
      }
    );
  }
}

void VTXdigi_Allpix2::InitCsvOutput() {
  m_debugCsv = true;

  warning() << " - You enabled the CSV output by setting `DebugCsvFileName` to a path. This is NOT MULTITHREADING SAFE and will cause crashes if multithreading is used." << endmsg;
  m_debugCsvFile.open(m_debugCsvFileName.value());

  if (!m_debugCsvFile.is_open()) {
    throw GaudiException("Failed to open debug CSV file", "VTXdigi_Allpix2::initialize", StatusCode::FAILURE);
  } else {
    m_debugCsvFile << "eventNumber,layerIndex,segmentCount,sensorThickness,pix_u,pix_v,simHitPos_u,simHitPos_v,simHitPos_w,simHitEntryPos_u,simHitEntryPos_v,simHitEntryPos_w,simHitPath_u,simHitPath_v,simHitPath_w,pathLengthGeant4,pathLength,rawChargeDeposition,debugFlag\n";
    m_debugCsvFile.flush();
    debug() << "   - writing to file: " << m_debugCsvFileName.value() << endmsg;
  }
}

void VTXdigi_Allpix2::PrintCounterSummary() const {
  const int colWidths[] = {65, 10};  
  info() << " Counters summary: " << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events read"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsRead.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events rejected (no simHits)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsRejected_noSimHits.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events accepted"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsAccepted.value() << " |" << endmsg;
  if (m_counter_eventsRead.value() != m_counter_eventsRejected_noSimHits.value() + m_counter_eventsAccepted.value())
    warning() << " | Number of accepted and rejected events does not add up to total number of processed events!" << endmsg;
  else info() << " | - event numbers add up." << endmsg;
  
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits read"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRead.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits rejected (in layer not to be digitized)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRejected_LayerNotToBeDigitized.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits rejected (below charge cut)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRejected_ChargeCut.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits rejected (path outside sensor)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRejected_OutsideSensor.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits accepted"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsAccepted.value() << " |" << endmsg;
  const long int simHitsRejected = (
    m_counter_simHitsRejected_LayerNotToBeDigitized.value() + 
    m_counter_simHitsRejected_ChargeCut.value() + 
    m_counter_simHitsRejected_OutsideSensor.value() );
  if (m_counter_simHitsRead.value() != simHitsRejected + m_counter_simHitsAccepted.value())
    warning() << " | Number of accepted and rejected simHits does not add up to total number of processed simHits!" << endmsg;
  else info() << " | - simHit numbers add up." << endmsg;

  info() << " | " << std::setw(colWidths[0]) << std::left << "Calculated path did not pass through the sensor active volume"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_acceptedButNoSegmentsInSensor.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "DigiHits created"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_digiHitsCreated.value() << " |" << endmsg;
}

/* -- Core algorithm functions -- */

bool VTXdigi_Allpix2::CheckEventSetup(const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {
  info() << "PROCESSING event (run " << headers.at(0).getRunNumber() << ", event " << headers.at(0).getEventNumber() << ", found " << simHits.size() << " simHits)" << endmsg;
  ++m_counter_eventsRead;

  /* early sanity checks to avoid segfaults from null pointers */
  if (!m_LUT)
    throw GaudiException("LookupTable is null in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::CheckEventSetup()", StatusCode::FAILURE);

  if (!m_simSurfaceMap)
    throw GaudiException("SimSurfaceMap is null in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::CheckEventSetup()", StatusCode::FAILURE);

  if (simHits.size()==0) {
    debug() << " - No SimTrackerHits in collection, returning empty output collections" << endmsg;
    ++m_counter_eventsRejected_noSimHits;
    return false;
  }

  ++m_counter_eventsAccepted;
  return true;
}

bool VTXdigi_Allpix2::CheckLayerCut(const edm4hep::SimTrackerHit& simHit) const {
  const int layer = m_cellIDdecoder->get(simHit.getCellID(), "layer");
  if (m_layersToDigitize.value().size()>0) { 
    if (std::find(m_layersToDigitize.value().begin(), m_layersToDigitize.value().end(),  layer) == m_layersToDigitize.value().end()) {
      verbose() << "   - DISMISSED SimHit in layer " << layer << ". (not in the list of layers to digitize)" << endmsg;
      ++m_counter_simHitsRejected_LayerNotToBeDigitized;
      return false;
    }
  }
  else {
    verbose() << "   - All layers are digitized, as property \"LayersToDigitize\" is not set ." << endmsg;
  }
  return true;
}

std::tuple<VTXdigi_Allpix2::HitInfo, VTXdigi_Allpix2::HitPosition> VTXdigi_Allpix2::GatherHitInfoAndPositions(const edm4hep::SimTrackerHit& simHit, const edm4hep::EventHeaderCollection& headers) const {
  HitInfo hitInfo(*this, simHit, headers);

  HitPosition hitPos;
  hitPos.global = ConvertVector(simHit.getPosition()); // global simHit position (from Geant4)
  hitPos.local = TransformGlobalToLocal(hitPos.global, hitInfo.cellID()); // the same position, but in the local sensor frame (u,v,w)
  std::tie(hitPos.entry, hitPos.path) = ConstructSimHitPath(hitInfo, hitPos, simHit); // simHit path through the sensor
  /* -> hitPos is now fully defined */

  hitInfo.setNSegments(std::max(1, int( hitPos.path.r() / m_targetPathSegmentLength ))); // both in mm
  /* -> with this, hitInfo is now fully defined as well */

  return std::make_tuple(hitInfo, hitPos);
}

bool VTXdigi_Allpix2::CheckSimHitCuts (const HitInfo& hitInfo, const HitPosition& hitPos) const {

  // DISMISS if entry point is outside sensor
  if (m_cutPathOutsideSensor) {

    if (abs(hitPos.entry.z()) > m_sensorThickness / 2 + m_numericLimit_float) { // entry point is outside sensor thickness
      verbose() << "   - DISMISSED simHit (entry point is outside sensor thickness (local w = " << hitPos.entry.z()*1000 << " um, sensor thickness = " << m_sensorThickness*1000 << " um)." << endmsg;
      ++m_counter_simHitsRejected_OutsideSensor;
      return false;
    }

    if (abs(hitPos.entry.z()+hitPos.path.z()) > m_sensorThickness / 2 + m_numericLimit_float) { // exit point is outside sensor thickness
      verbose() << "   - DISMISSED simHit (exit point is outside sensor thickness (local w = " << (hitPos.entry.z()+hitPos.path.z())*1000 << " um, sensor thickness = " << m_sensorThickness*1000 << " um)." << endmsg;
      ++m_counter_simHitsRejected_OutsideSensor;
      return false;
    }
  }

  // DISMISS if outside minimum charge cut
  if (hitInfo.charge() < m_cutDepositedCharge.value()) {
    ++m_counter_simHitsRejected_ChargeCut;
    verbose() << "   - DISMISSED simHit (charge below cut, " << hitInfo.charge() << " e-)" << endmsg;
    return false;
  }

  debug() << "   - ACCEPTED SimHit in layer " << m_cellIDdecoder->get(hitInfo.cellID(), "layer") << ". Charge =" << hitInfo.charge() << " e-. Starting loop over segments" << endmsg;
  return true;
}

std::tuple<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D> VTXdigi_Allpix2::ConstructSimHitPath(HitInfo& hitInfo, HitPosition& hitPos, const edm4hep::SimTrackerHit& simHit) const {
  /* Get the simHitPath (the path of the simulated particle through a sensor) of a given simHit. Describet by the path direction (unit vector) and entry point (local frame).
   *  
   * The way I do this differs from what Jessy does in their VTXdigi_Detailed.
   *  My approach ensures that the simHitPath always extends exactly from one sensor surface to the other, and that the simHitPos always lies on this path. 
   * Then, the paths are shortened to the PathLength that Geant4 gives us, if they are longer (the Geant4 pathLength is always longer than the linear approx., because it accounts for B-field and multiple scattering inside the sensor). We need to do this to avoid unphysically long paths, where the charge-per-pathLength would be unphysically tiny
   *  I am not sure if my approach is better or worse, it definitely produces less problems. The problem is that assuming paths are linear is generally not compatible with setting the path-vectors length to the path-length (as it is given by Geant4, where paths include multiple scattering and B-field effects).
   */
  
  /* get simHitMomentum (this vector gives the exact angles of the particle's path through the sensor. We assume that path is linear.). */
  double simHitGlobalMomentum_double[3] = {simHit.getMomentum().x * dd4hep::GeV, simHit.getMomentum().y * dd4hep::GeV, simHit.getMomentum().z * dd4hep::GeV}; // need floats for the TransfomationMatrix functions

  TGeoHMatrix transformationMatrix = ComputeTransformationMatrix(hitInfo.cellID());

  double simHitLocalMomentum_double[3] = {0.0, 0.0, 0.0};
  transformationMatrix.MasterToLocalVect(simHitGlobalMomentum_double, simHitLocalMomentum_double);

  /* Scale the simHitPath such that it extends from one sensor surface to the other */
  dd4hep::rec::Vector3D simHitPath(
    simHitLocalMomentum_double[0],
    simHitLocalMomentum_double[1],
    simHitLocalMomentum_double[2]);

  float scaleFactor = m_sensorThickness / std::abs(simHitLocalMomentum_double[2]);
  simHitPath = scaleFactor * simHitPath ; // now, simHitPath extends for one sensor surface to the other surface
    
  /* calculate the path's entry position into the sensor, by placing it such that the path passes through the simHit position */
  if (abs(hitPos.local.z()) > (0.5*m_sensorThickness + m_numericLimit_float)) {
    warning() << "SimHit position is outside the sensor volume (local w = " << hitPos.local.z() << " mm, sensor thickness = " << m_sensorThickness << " mm). This should never happen. Forcing it to w=0." << endmsg;
    hitPos.local.z() = 0.;
  }
  scaleFactor = 0.; // this will now hold the fraction of the simHitPath between simHitPos and entry point, in terms of [0,1] on simHitPath
  if (simHitPath.z() >= 0.) {
    const float shiftDist_w = 0.5 * m_sensorThickness + hitPos.local.z();
    scaleFactor = shiftDist_w / simHitPath.z();
  }
  else {
    const float shiftDist_w = -0.5 * m_sensorThickness + hitPos.local.z();
    scaleFactor = shiftDist_w / simHitPath.z();
  }
  dd4hep::rec::Vector3D simHitEntryPos = hitPos.local - scaleFactor * simHitPath; // entry pos is now on one of the two sensor surfaces

  /* If the path passes through the sensor edges, clip it to the edges in u and v direction */
  float t_min = 0.f, t_max = 1.f;
  std::tie(t_min, t_max) = ComputePathClippingFactors(t_min, t_max, simHitEntryPos.x(), simHitPath.x(), m_sensorLength.at(0));
  std::tie(t_min, t_max) = ComputePathClippingFactors(t_min, t_max, simHitEntryPos.y(), simHitPath.y(), m_sensorLength.at(1));
  if (t_min != 0.f || t_max != 1.f) { // check if clipping is even necessary, for performance
    hitInfo.setDebugFlag();
    if (0. <= t_min && t_min <= t_max && t_max <= 1.) {
      verbose() << "   - Clipping simHitPath to sensor edges with factors t_min = " << t_min << ", t_max = " << t_max << ". PathLength changed to " << (t_max - t_min) * simHitPath.r() << " mm from " << simHitPath.r() << " mm" << endmsg;
      
      simHitEntryPos = simHitEntryPos + t_min * simHitPath;
      simHitPath = (t_max - t_min) * simHitPath;
    } 
    else {
      warning() << "ConstructSimHitPath(): Cannot clip simHitPath to sensor edges. The path lies completely outside the sensor volume. Clipping t_min = " << t_min << ", t_max = " << t_max << "." << endmsg;
      verbose() << "   - before clipping: EntryPos: (" << simHitEntryPos.x() << " mm, " << simHitEntryPos.y() << " mm, " << simHitEntryPos.z() << " mm), exitPos: (" << simHitPath.x() << " mm, " << simHitPath.y() << " mm, " << simHitPath.z() << " mm)" << endmsg;
    }
  }

  /* if pathLength given by Geant4 is more than X% shorter than the length we calculate, shorten our calculated path accordingly. */
  if (simHitPath.r() > m_pathLengthShorteningFactorGeant4 * hitInfo.simPathLength()) {
    verbose() << "   - Shortening simHitPath from " << simHitPath.r() << " mm to Geant4 pathLength of " << hitInfo.simPathLength() << " mm, because it's length is more than " << m_pathLengthShorteningFactorGeant4 << " of the Geant4 path length." << endmsg;
    
    hitInfo.setDebugFlag();
    /* make sure the path stays as centered around the simHitPos as possible
    * find out where the simHitPos lies on the current path:
    * project (simHitPos - simHitEntryPos) onto simHitPath -> gives distance from entryPos to simHitPos along the path (or to the point on path closest to simHitPos)
    * dotProduct (simHitPos-simHitEntryPos, simHitPath) = |simHitPos - simHitEntryPos| * |simHitPath| * cos(angle between them) */
    const float t_simHit = ( (hitPos.local - simHitEntryPos).dot(simHitPath) ) / ( simHitPath.r() * simHitPath.r() ); 
    
    /* clamp new path (centered at t_center) to [0,1] */
    const float t_length_half = hitInfo.simPathLength() / simHitPath.r() / 2.f; // half-length of new path, in terms of t [0,1] on old path
    const float t_center = std::max(t_length_half, std::min(t_simHit, 1.f - t_length_half)); // center of new path, clamped to [t_length_half, 1 - t_length_half] while not exceeding [0,1]

    t_min = t_center - t_length_half;
    t_max = t_center + t_length_half;

    simHitEntryPos = simHitEntryPos + t_min * simHitPath;
    simHitPath = (t_max - t_min) * simHitPath;
  }

  verbose() << "   - Calculated SimHitPath with length " << simHitPath.r() << " mm" << " (the Geant4 path length is " << hitInfo.simPathLength() << " mm)" << endmsg;
  return std::make_tuple(simHitEntryPos, simHitPath);
}

VTXdigi_Allpix2::PixelChargeMatrix VTXdigi_Allpix2::DepositAndCollectCharge(HitInfo& hitInfo, const HitPosition& hitPos) const {
  const auto& [i_u_simHit, i_v_simHit] = ComputePixelIndices(hitPos.local, m_sensorLength.at(0), m_sensorLength.at(1));
  
  PixelChargeMatrix pixelChargeMatrix(
    i_u_simHit,
    i_v_simHit,
    m_maxClusterSize.value().at(0),
    m_maxClusterSize.value().at(1)
  );

  const float segmentCharge = hitInfo.charge() / hitInfo.nSegments();
  
  /* get first segment */
  SegmentIndices segment = ComputeSegmentIndices(hitInfo, hitPos.entry, hitPos.path, 0);
  int segmentsInBin = 1;
  SegmentIndices nextSegment;

  /* loop over segments */
  for (int nextSegmentIndex = 1; nextSegmentIndex < hitInfo.nSegments(); ++nextSegmentIndex) {

    nextSegment = ComputeSegmentIndices(hitInfo, hitPos.entry, hitPos.path, nextSegmentIndex);

    if (segment == nextSegment) {
      verbose() << "       - Segment lies in the same pixel and in-pixel bin as previous segment, continuing." << endmsg;
      ++segmentsInBin;
      continue;
    } 
    else {
      verbose() << "       - Crossed bin-boundary wrt. last segment. Sharing " << segmentCharge*segmentsInBin << " e- from " << segmentsInBin << " segments. The last segment of these has nextSegmentIndex " << nextSegmentIndex-1 << "." << endmsg;
      DistributeSegmentCharge(hitInfo, pixelChargeMatrix, segment, segmentCharge); // write charge for this set of segments into pixelChargeMatrix (avoids copying the matrix in memory every time we write into it)

      segment = nextSegment;
      segmentsInBin = 1;
    }
  } // loop over segments

  /* write out last set of segments */
  verbose() << "       - Reached last segment. Sharing " << segmentCharge*segmentsInBin << " e- from last " << segmentsInBin << " segments." << endmsg;
  DistributeSegmentCharge(hitInfo, pixelChargeMatrix, segment, segmentCharge); // write charge for this segment into pixelChargeMatrix (done like this to avoid copying the matrix in memory every time we write into it)

  return pixelChargeMatrix;
}


void VTXdigi_Allpix2::DistributeSegmentCharge(HitInfo& hitInfo, PixelChargeMatrix& pixelChargeMatrix, const SegmentIndices& segment, const float segmentCharge) const {
  if (segment.i_u == -1) { // ComputeSegmentIndices() returns -1 if any dimension is outside sensor volume
      warning() << "Applying Kernel: Bin lies outside sensor volume. Dismissing." << endmsg;
      return; 
    }

  int i_u_target, i_v_target;

  /* each matrix entry shares charge from the source-pixel i_u, i_v to one target-pixel.
   * The matrix is centered on the source pixel, so loop over all target pixels it covers.
   * i_x_previous defines the source-pixel */
  for (int i_m = -1*(m_matrixSize-1)/2; i_m<=(m_matrixSize-1)/2; i_m++) {
    i_u_target = segment.i_u + i_m;
    if (i_u_target<0 || i_u_target>=m_pixelCount.at(0))
      continue; // target pixel outside pixel matrix in u

    for (int i_n = -1*(m_matrixSize-1)/2; i_n<=(m_matrixSize-1)/2; i_n++) {
      i_v_target = segment.i_v + i_n;

      if (i_v_target<0 || i_v_target>=m_pixelCount.at(1))
        continue; // target pixel outside pixel matrix in v

      const float weight = m_LUT->GetWeight(segment, i_m, i_n);
      if (weight < m_numericLimit_float) 
        continue; // skip zero entries
        
      const float sharedCharge = weight * segmentCharge;
      pixelChargeMatrix.FillRawCharge(i_u_target, i_v_target, sharedCharge);

      if (m_debugHistograms)
        FillHistograms_PerSegment(hitInfo, segment, i_m, i_n, sharedCharge);
    }
  }
}

void VTXdigi_Allpix2::AnalyseSharedCharge(const HitInfo& hitInfo, const HitPosition& hitPos, const PixelChargeMatrix& pixelChargeMatrix, const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks) const {
  /** Process the shared charge in pixelChargeMatrix, create digiHits and fill digiHits and digiHitsLinks collections */

  debug() << "       - Processed all segments. Looping over pixelChargeMatrix." << endmsg;
  
  int nPixelsReceivedCharge = 0;
  int nPixelsFired = 0;


  for (int i_u = pixelChargeMatrix.GetRangeMin_u(); i_u <= pixelChargeMatrix.GetRangeMax_u(); ++i_u) {
    for (int i_v = pixelChargeMatrix.GetRangeMin_v(); i_v <= pixelChargeMatrix.GetRangeMax_v(); ++i_v) {

      float pixelChargeRaw = pixelChargeMatrix.GetRawCharge(i_u, i_v);
      if (pixelChargeRaw < m_numericLimit_float) 
        continue; 
      nPixelsReceivedCharge++;

      /* skip pixels with zero or negative charge BEFORE noise addition, st. each pixel without charge is treated the same.
        * Otherwise noise would appear only around simHits (because pixelChargeMatrix only holds pixels around simHit).
        * This would bias the results. */

      float pixelChargeMeasured = pixelChargeRaw + pixelChargeMatrix.GetNoise(i_u, i_v);

      if (pixelChargeMeasured < m_pixelThreshold) {
        verbose() << "       - Pixel (" << i_u << ", " << i_v << ") received a measured/raw charge of " << pixelChargeMeasured << " / " << pixelChargeRaw << " e-. This is below threshold (" << m_pixelThreshold << " e-). Discarding pixel." << endmsg;
        continue;
      }

      nPixelsFired++;
      dd4hep::rec::Vector3D pixelCenterLocal = ComputePixelCenter_Local(i_u, i_v, *hitInfo.simSurface());
      dd4hep::rec::Vector3D pixelCenterGlobal = TransformLocalToGlobal(pixelCenterLocal, hitInfo.cellID());

      debug() << "       - Pixel (" << i_u << ", " << i_v << ") at (" << pixelCenterLocal.x() << ", " << pixelCenterLocal.y() << ", " << pixelCenterLocal[2] << ") mm received a measured/raw charge of " << pixelChargeMeasured << "/" << pixelChargeRaw << " e-, center at global position " << pixelCenterGlobal[0] << " mm, " << pixelCenterGlobal[1] << " mm, " << pixelCenterGlobal[2] << " mm" << endmsg;
      CreateDigiHit(simHit, digiHits, digiHitsLinks, pixelCenterGlobal, pixelChargeMeasured);
      ++m_counter_digiHitsCreated;

      if (m_debugHistograms)
        FillHistograms_PerPixelHit(hitInfo, hitPos, i_u, i_v, pixelChargeMeasured);
    }
  } // end loop over pixels, create digiHits

  if (m_debugHistograms) {
    ++(*m_histGlobal.at(histGlobal_clusterSize_raw))[static_cast<float>(nPixelsReceivedCharge)]; // cluster size = number of pixels fired per simHit
    ++(*m_histGlobal.at(histGlobal_clusterSize_measured))[static_cast<float>(nPixelsFired)];
    
    ++(*m_hist1d.at(hitInfo.layerIndex()).at(hist1d_ClusterSize_raw))[static_cast<float>(nPixelsReceivedCharge)];
    ++(*m_hist1d.at(hitInfo.layerIndex()).at(hist1d_ClusterSize_measured))[static_cast<float>(nPixelsFired)];

    (*m_histProfile1d.at(hitInfo.layerIndex()).at(histProfile1d_clusterSize_vs_z))[hitPos.global.z()] += nPixelsFired;
    ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_clusterSize_vs_z))[{hitPos.global.z(), nPixelsFired}];
    const float moduleZ = TransformLocalToGlobal(dd4hep::rec::Vector3D(0.f, 0.f, 0.f), hitInfo.cellID()).z();
    (*m_histProfile1d.at(hitInfo.layerIndex()).at(histProfile1d_clusterSize_vs_module_z))[moduleZ] += nPixelsFired;
    ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_clusterSize_vs_module_z))[{moduleZ, nPixelsFired}];
    const int moduleID = m_cellIDdecoder->get(hitInfo.cellID(), "module");
    (*m_histProfile1d.at(hitInfo.layerIndex()).at(histProfile1d_clusterSize_vs_moduleID))[moduleID] += nPixelsFired;
  }
    
  if (nPixelsFired <= 0) {
    ++m_counter_acceptedButNoSegmentsInSensor;
    debug() << "   - No pixels fired for this simHit" << endmsg;
  }
  else {
    verbose() << "   - From this simHit, " << nPixelsFired << " pixels received charge." << endmsg;
  }
} // AnalyseSharedCharge()

void VTXdigi_Allpix2::FillHistograms_PerSimHit(HitInfo& hitInfo, const HitPosition& hitPos, const PixelChargeMatrix& pixelChargeMatrix) const {
  /* Is executed once per simHit (that passes all cuts etc.), after it has been processed  */
  verbose() << "   - Filling 1D histograms" << endmsg;
  ++(*m_histGlobal.at(histGlobal_simHitCharge))[hitInfo.charge()]; // in e
  ++(*m_histGlobal.at(histGlobal_simHitE))[hitInfo.charge() / m_chargePerkeV]; // in keV
  ++(*m_histGlobal.at(histGlobal_EntryPointX))[hitPos.entry.x()]; // in mm
  ++(*m_histGlobal.at(histGlobal_EntryPointY))[hitPos.entry.y()]; // in mm
  ++(*m_histGlobal.at(histGlobal_EntryPointZ))[hitPos.entry.z()*1000]; // convert mm to um

  ++(*m_histGlobal.at(histGlobal_pathLength))[hitPos.path.r()*1000]; // convert mm to um
  ++(*m_histGlobal.at(histGlobal_pathLengthGeant4))[hitInfo.simPathLength()*1000]; // convert mm to um
  ++(*m_histGlobal2d.at(histGlobal2d_pathLength_vs_G4PathLength))[{hitInfo.simPathLength()*1000 , hitPos.path.r()*1000}]; // both in um

  ++(*m_histGlobal.at(histGlobal_simHitPDG))[static_cast<float>(hitInfo.simPdg())];

  ++(*m_hist1d.at(hitInfo.layerIndex()).at(hist1d_SimHitMomentum))[hitInfo.simMomentum()*1000]; // Convert GeV to MeV

  const float clusterChargeRaw = pixelChargeMatrix.GetTotalRawCharge();
  const float clusterChargeMeasured = pixelChargeMatrix.GetTotalMeasuredCharge(m_pixelThreshold);
  ++(*m_hist1d.at(hitInfo.layerIndex()).at(hist1d_DigiHitCharge_raw))[ clusterChargeRaw ];
  ++(*m_hist1d.at(hitInfo.layerIndex()).at(hist1d_DigiHitCharge_measured))[ clusterChargeMeasured ];
  ++(*m_histGlobal.at(histGlobal_chargeCollectionEfficiency_raw))[ clusterChargeRaw / hitInfo.charge() ]; // in e-, only if at least one pixel fired
  ++(*m_histGlobal.at(histGlobal_chargeCollectionEfficiency))[ clusterChargeMeasured / hitInfo.charge()];
  ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_totalCharge_vs_simHitCharge))[{hitInfo.charge(), clusterChargeMeasured}];

  int pix_u, pix_v;
  std::tie(pix_u, pix_v) = ComputePixelIndices(hitPos.local, m_sensorLength.at(0), m_sensorLength.at(1));
  verbose() << "   - Filling 2D histograms for layer " << m_cellIDdecoder->get(hitInfo.cellID(), "layer") << " (index " << hitInfo.layerIndex() << ")" << endmsg;
  ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_hitMap_simHits))[{pix_u, pix_v}];
  ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_pathLength_vs_simHit_v))[{hitPos.global.z(), hitPos.path.r()*1000}]; // in um and mm

  ++(*m_histGlobal.at(histGlobal_pixelChargeMatrix_size_u))[pixelChargeMatrix.GetSize_u()];
  ++(*m_histGlobal.at(histGlobal_pixelChargeMatrix_size_v))[pixelChargeMatrix.GetSize_v()];
  
  ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_pixelChargeMatrixSize))[{pixelChargeMatrix.GetSize_u(), pixelChargeMatrix.GetSize_v()}];

  const float pathAngleTheta = acos(hitPos.path.z() / hitPos.path.r()) / 3.14159265f * 180.f; // in degrees
  const float pathAnglePhi = atan2(hitPos.path.y(), hitPos.path.x()) / 3.14159265f * 180.f; // in degrees. v is .y(), and parallel to global 

  ++(*m_hist1d.at(hitInfo.layerIndex()).at(hist1d_IncidentAngle_ThetaLocal))[pathAngleTheta];
  ++(*m_hist1d.at(hitInfo.layerIndex()).at(hist1d_IncidentAngle_PhiLocal))[pathAnglePhi];
  ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_IncidentAngle))[{pathAnglePhi, pathAngleTheta}];
}

void VTXdigi_Allpix2::FillHistograms_PerSegment(const HitInfo& hitInfo, const SegmentIndices& segment, int i_m, int i_n, const float sharedCharge) const {

  /* work out the distance between target pixel center and the origin bin, in terms of pixels */
  const float dist_u = -i_m - 0.5f + (segment.j_u + 0.5f) / static_cast<float>(m_inPixelBinCount[0]); // in pixels
  const float dist_v = -i_n - 0.5f + (segment.j_v + 0.5f) / static_cast<float>(m_inPixelBinCount[1]); // in pixels
  const float pos_w = ( (segment.j_w + 0.5f) * (m_sensorThickness / m_inPixelBinCount[2]) - m_sensorThickness/2.f) * 1000.f; // conversion from mm to um
  
  (*m_histWeighted2d.at(hitInfo.layerIndex()).at(histWeighted2d_chargeOriginU))[{dist_u, pos_w}] += sharedCharge; // in e-
  (*m_histWeighted2d.at(hitInfo.layerIndex()).at(histWeighted2d_chargeOriginV))[{dist_v, pos_w}] += sharedCharge; // in e-
}

void VTXdigi_Allpix2::FillHistograms_PerPixelHit(const HitInfo& hitInfo, const HitPosition& hitPos, int i_u, int i_v, float pixelChargeMeasured) const {
  /* executed for every digiHit (ie. for every pixel that is above threshold after charge collection and noise) */

  const auto& [i_u_simHit, i_v_simHit] = ComputePixelIndices(hitPos.local, m_sensorLength.at(0), m_sensorLength.at(1));

  ++(*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_hitMap_digiHits))[{i_u, i_v}];
  (*m_histWeighted2d.at(hitInfo.layerIndex()).at(histWeighted2d_averageCluster_analog))[{i_u - i_u_simHit, i_v - i_v_simHit}] += pixelChargeMeasured; // in e-
  ++ (*m_hist2d.at(hitInfo.layerIndex()).at(hist2d_averageCluster_binary))[{i_u - i_u_simHit, i_v - i_v_simHit}];

  dd4hep::rec::Vector3D pixelCenterLocal = ComputePixelCenter_Local(i_u, i_v, *hitInfo.simSurface());
  dd4hep::rec::Vector3D pixelCenterGlobal = TransformLocalToGlobal(pixelCenterLocal, hitInfo.cellID());

  
  ++ (*m_histGlobal.at(histGlobal_DisplacementU))[ (pixelCenterLocal.x() - hitPos.local.x()) * 1000]; // convert mm to um
  ++ (*m_histGlobal.at(histGlobal_DisplacementV))[ (pixelCenterLocal.y() - hitPos.local.y()) * 1000]; // convert mm to um
  ++ (*m_histGlobal.at(histGlobal_DisplacementR))[sqrt( (pixelCenterGlobal.x() - hitPos.global.x())*(pixelCenterGlobal.x() - hitPos.global.x()) + (pixelCenterGlobal.y() - hitPos.global.y())*(pixelCenterGlobal.y() - hitPos.global.y()) + (pixelCenterGlobal.z() - hitPos.global.z())*(pixelCenterGlobal.z() - hitPos.global.z()) ) * 1000]; // convert mm to um
}

/* -- Pixel and in-pixel binning magic (logic) -- */

int VTXdigi_Allpix2::ComputeBinIndex(float x, float binX0, float binWidth, int binN) const {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  if (binN <= 0) throw GaudiException("ComputeBinIndex: binN must be positive", "VTXdigi_Allpix2::ComputeBinIndex()", StatusCode::FAILURE);
  if (binWidth <= 0.0) throw GaudiException("ComputeBinIndex: binWidth must be positive", "VTXdigi_Allpix2::ComputeBinIndex()", StatusCode::FAILURE);

  float relativePos = (x - binX0) / binWidth; // shift to [0, binN]
  if (relativePos < 0.0f || relativePos >= static_cast<float>(binN))
    return -1;
  return static_cast<int>(relativePos);
} // ComputeBinIndex()

std::tuple<int, int> VTXdigi_Allpix2::ComputePixelIndices(const dd4hep::rec::Vector3D& segmentPos, const float length_u, const float length_v) const {
  
  int i_u = ComputeBinIndex(
    segmentPos.x(),
    -0.5*length_u,
    m_pixelPitch.at(0),
    m_pixelCount.at(0));

  int i_v = ComputeBinIndex(
    segmentPos.y(),
    -0.5*length_v,
    m_pixelPitch.at(1),
    m_pixelCount.at(1));
  return std::make_tuple(i_u, i_v);
} // ComputePixelIndices()

std::tuple<int, int, int> VTXdigi_Allpix2::ComputeInPixelIndices(const dd4hep::rec::Vector3D& segmentPos, const float length_u, const float length_v) const {
  int j_u, j_v, j_w;

  verbose() << "         - Computing in-pixel indices. Number of in-pix bins: (" << m_inPixelBinCount[0] << ", " << m_inPixelBinCount[1] << ", " << m_inPixelBinCount[2] << ")" << endmsg;

  float shiftedPos_u = segmentPos.x() + 0.5 * length_u; // shift to [0, length_u]
  float pitch_u = m_pixelPitch.at(0);
  float inPixelPos_u = std::fmod(shiftedPos_u, pitch_u);
  if (inPixelPos_u < 0.0) inPixelPos_u += pitch_u; // ensure positive remainder

  j_u = ComputeBinIndex(inPixelPos_u, 0.0, pitch_u / m_inPixelBinCount[0], m_inPixelBinCount[0]);

  float shiftedPos_v = segmentPos.y() + 0.5 * length_v;
  float pitch_v = m_pixelPitch.at(1);
  float inPixelPos_v = std::fmod(shiftedPos_v, pitch_v);
  if (inPixelPos_v < 0.0) inPixelPos_v += pitch_v;

  j_v = ComputeBinIndex(inPixelPos_v, 0.0, pitch_v / m_inPixelBinCount[1], m_inPixelBinCount[1]);

  // vertical (w) binning: shift to [0, thickness]
  float shiftedPos_w = segmentPos.z() + 0.5 * m_sensorThickness;
  j_w = ComputeBinIndex(shiftedPos_w, 0.0, m_sensorThickness / m_inPixelBinCount[2], m_inPixelBinCount[2]);

  return std::make_tuple(j_u, j_v, j_w);
} // ComputeInPixelIndices()

VTXdigi_Allpix2::SegmentIndices VTXdigi_Allpix2::ComputeSegmentIndices(HitInfo& hitInfo, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const int segmentIndex) const {

  verbose() << "     - Processing segment (" << segmentIndex << " out of " << hitInfo.nSegments() << ", length " << simHitPath.r() << " mm)" << endmsg;
  if (segmentIndex < 0 || segmentIndex >= hitInfo.nSegments()) {
    error() << "ComputeSegmentIndices(): Invalid segment number " << hitInfo.nSegments() << " or segment index " << segmentIndex << endmsg;
    throw std::runtime_error("VTXdigi_Allpix2::ComputeSegmentIndices(): Invalid segment number or segment index");
  }
  
  SegmentIndices segment;
  float pathFraction = (static_cast<float>(segmentIndex)+0.5) / hitInfo.nSegments();
  const dd4hep::rec::Vector3D segmentRelativePos = pathFraction * simHitPath;
  const dd4hep::rec::Vector3D segmentPos = simHitEntryPos + segmentRelativePos;
  verbose() << "       - Local position (u,v,w): (" << segmentPos.x() << " mm, " << segmentPos.y() << " mm, " << segmentPos.z() << " mm)" << endmsg;

  /* Compute pixel indices */
  std::tie(segment.i_u, segment.i_v) = ComputePixelIndices(segmentPos, m_sensorLength.at(0), m_sensorLength.at(1));
  if (segment.i_u==-1 || segment.i_v==-1) {
    warning() << "ComputeSegmentIndices(): Segment lies outside sensor area (in u or v). Dismissing." << endmsg;
    // hitInfo.setDebugFlag();
    SegmentIndices emptySegment;
    return emptySegment;
  }
  
  /* Compute in-pixel indices */
  std::tie(segment.j_u, segment.j_v, segment.j_w) = ComputeInPixelIndices(segmentPos, m_sensorLength.at(0), m_sensorLength.at(1));
  if (segment.j_u==-1 || segment.j_v==-1 || segment.j_w==-1) {
    warning() << "ComputeSegmentIndices(): Segment lies inside sensor area (in u and v), but vertically outside sensor volume. Dismissing." << endmsg;
    SegmentIndices emptySegment;
    return emptySegment;
  }
  verbose() << "       - Pixel indices (" << segment.i_u << ", " << segment.i_v << "), In-pixel indices (" << segment.j_u << ", " << segment.j_v << ", " << segment.j_w << ")" << endmsg;

  return segment;
}

dd4hep::rec::Vector3D VTXdigi_Allpix2::ComputePixelCenter_Local(const int i_u, const int i_v, const dd4hep::rec::ISurface& simSurface) const {
  /* returns the position of the center of pixel i_u, i_v) in the global
   * u - short ARCADIA axis, corresponds to x in sensor local frame
   * v - long ARCADIA axis, corresponds to y in sensor local frame */

  float length_u = simSurface.length_along_u() * 10; // convert to mm (works, checked 2025-10-17)
  float length_v = simSurface.length_along_v() * 10; // convert to mm 

  float posU = -0.5 * length_u + (i_u + 0.5) * m_pixelPitch.at(0); // in mm
  float posV = -0.5 * length_v + (i_v + 0.5) * m_pixelPitch.at(1); // in mm

  return dd4hep::rec::Vector3D(posU, posV, 0.); 
}


/* -- Transformation between global frame and local sensor frame -- */

TGeoHMatrix VTXdigi_Allpix2::ComputeTransformationMatrix(const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = m_volumeManager.lookupDetElement(cellID).nominal().worldTransformation(); // given in cm

  /* rotation is unitless, but need to convert translation from cm to mm */
  double* translationComponent = transformationMatrix.GetTranslation();
  translationComponent[0] = translationComponent[0] * 10; // convert to mm
  translationComponent[1] = translationComponent[1] * 10;
  translationComponent[2] = translationComponent[2] * 10;
  transformationMatrix.SetTranslation(translationComponent);

  /* Rotate and mirror the transformationMatrix to have a direct frame with z orthogonal to sensor simSurface
   * (a) right-handed and (b) has z perpendicular to sensor plane. see the checks in InitDetectorGeometry() for details.
   * This is necessary to correctly calculate the drift in X-Y due to B-field, as our coordinate system needs to match the one in AP2.
   * I basically copied this from https://github.com/jessy-daniel/k4RecTracker/blob/New_Detailed_VTXdigitizer/VTXdigiDetailed/src/VTXdigitizerDetailed.cpp */
  std::string localNormalVectorDir = m_localNormalVectorDir;
  bool IsDirect = true; // Is the origin frame direct ?
  if (localNormalVectorDir[0]=='-') {
    IsDirect = false;
    localNormalVectorDir = localNormalVectorDir[1];
  }  

  /* If the orthogonal direction is along X or Y in local frame, rotate the frame to have Z orthogonal to sensors instead */
  if (localNormalVectorDir=="x") {
    TGeoRotation rot("rot",90.,90.,0.); // X->Z / Y->X / Z->Y
    transformationMatrix.Multiply(rot);
  }
  else if (localNormalVectorDir=="y") {
    TGeoRotation rot("rot",0.,-90.,0.); // X->X / Y->Z / Z->-Y
    transformationMatrix.Multiply(rot);
  }
  else if (localNormalVectorDir!="z") {
    throw std::runtime_error("VTXdigi_Allpix2::ComputeTransformationMatrix(): Invalid localNormalVectorDir. Must be 'x', 'y', 'z', '-x', '-y' or '-z' (use x for IDEA vertex barrel).");
  }

  /* If the frame isn't direct, make it direct by reflecting the x axis. This is necessary to correctly calculte the drift in X-Y due to B-field
   * TODO: this is does not cover all cases, ie. needing to relect y or z. Maybe this is not necessary due to constraints in dd4hep. ~ Jona, 2025-11 */
  if (!IsDirect) {
    transformationMatrix.ReflectX(false);
  }

  return transformationMatrix;
}
TGeoHMatrix VTXdigi_Allpix2::ComputeTransformationMatrix(const edm4hep::SimTrackerHit& simHit) const {
  return ComputeTransformationMatrix(simHit.getCellID());
}

dd4hep::rec::Vector3D VTXdigi_Allpix2::TransformGlobalToLocal(const dd4hep::rec::Vector3D& globalPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = ComputeTransformationMatrix(cellID);

  double localPos[3] = {0, 0, 0};

  transformationMatrix.MasterToLocal(globalPos, localPos);

  return dd4hep::rec::Vector3D(localPos[0], localPos[1], localPos[2]);
}

dd4hep::rec::Vector3D VTXdigi_Allpix2::TransformLocalToGlobal(const dd4hep::rec::Vector3D& localPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = ComputeTransformationMatrix(cellID);

  double globalPos[3] = {0, 0, 0};

  transformationMatrix.LocalToMaster(localPos, globalPos);

  return dd4hep::rec::Vector3D(globalPos[0], globalPos[1], globalPos[2]);
}


/* -- Other Helper functions -- */

dd4hep::rec::Vector3D VTXdigi_Allpix2::ConvertVector(edm4hep::Vector3d vec) const {
  // return dd4hep::rec::Vector3D(vec.x*dd4hep::mm, vec.y*dd4hep::mm, vec.z*dd4hep::mm);
  return dd4hep::rec::Vector3D(vec.x, vec.y, vec.z);
}
edm4hep::Vector3d VTXdigi_Allpix2::ConvertVector(dd4hep::rec::Vector3D vec) const {
  // return edm4hep::Vector3d(vec[0]/dd4hep::mm, vec[1]/dd4hep::mm, vec[2]/dd4hep::mm);
  return edm4hep::Vector3d(vec.x(), vec.y(), vec.z());
}

std::tuple<float, float> VTXdigi_Allpix2::ComputePathClippingFactors(float t_min, float t_max, const float entryPos_ax, const float pathLength_ax, const float sensorLength_ax) const {
  /** Calculate the clipping factor for a single dimension (u or v)
   * 
   * Used to calculate the total clipping factor for the simHitPath
   */

  const bool posDirection = (entryPos_ax + pathLength_ax > entryPos_ax) ? true : false; // whether the path goes in positive direction along this axis (false = negative direction)
  
  const float minPos = std::min(entryPos_ax, entryPos_ax + pathLength_ax);
  if (2*minPos < -sensorLength_ax) {
    // path extends out of sensor in - direction

    float t = (-minPos - 0.5*sensorLength_ax) / abs(pathLength_ax); // fraction of pathLength_ax that is outside sensor (in this axis)
    
    if (posDirection)
    t_min = std::max(t_min, t);
    else
    t_max = std::min(t_max, 1 - t);
    
    debug() << "   - SimHitPath extends outside sensor on NEG. side, in " << (posDirection ? "POS" : "NEG") << ". direction. (min at " << minPos << " mm, edge at " << -0.5*sensorLength_ax << " mm) => " << -0.5*sensorLength_ax - minPos << "mm outside of the sensor, " << t*100 << " percent of the path length. Clipping to [" << t_min << ", " << t_max << "]" << endmsg;
  }
  
  const float maxPos = std::max(entryPos_ax, entryPos_ax + pathLength_ax);
  if (2*maxPos > sensorLength_ax) {

    float t = (maxPos - 0.5*sensorLength_ax) / abs(pathLength_ax);
    
    if (posDirection)
    t_max = std::min(t_max, 1 - t);
    else
    t_min = std::max(t_min, t);

    debug() << "   - SimHitPath extends outside sensor on POS. side, in " << (posDirection ? "POS" : "NEG") << ". direction. (min at " << maxPos << " mm, edge at " << -0.5*sensorLength_ax << " mm) => " << maxPos - 0.5*sensorLength_ax << "mm outside of the sensor, " << t*100 << " percent of the path length. Clipping to [" << t_min << ", " << t_max << "]" << endmsg;

  }

  return std::make_tuple(t_min, t_max);
}

void VTXdigi_Allpix2::CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const dd4hep::rec::Vector3D& position, const float charge) const {
  // overload to allow passing dd4hep::rec::Vector3D as position. ~ Jona 2025-09

  auto digiHit = digiHits.create();
  digiHit.setCellID(simHit.getCellID());
  verbose() << "         - Creating digiHit in cellID " << simHit.getCellID() << " (" << std::bitset<24>(simHit.getCellID()) << "), setting eDep to " << charge/m_chargePerkeV << " keV" << endmsg;
  digiHit.setEDep(charge / m_chargePerkeV); // convert e- to keV
  digiHit.setPosition(ConvertVector(position));
  // TODO: check if position is within sensor bounds & force it onto sensor simSurface ~ Jona 2025-09
  digiHit.setTime(simHit.getTime());
  
  auto digiHitLink = digiHitsLinks.create();
  digiHitLink.setFrom(digiHit);
  digiHitLink.setTo(simHit);
}

void VTXdigi_Allpix2::AppendSimHitToCsv(const HitInfo& hitInfo, const HitPosition& hitPos, const int i_u, const int i_v) const {
  if (!m_debugCsvFile.is_open()) {
    error() << "AppendSimHitToCsv(): DebugCsv file is not open." << endmsg;
    return;
  }

  // eventNumber,layerIndex,segmentCount,sensorThickness,pix_u,pix_v,
  // simHitPos_u,simHitPos_v,simHitPos_w,
  // simHitEntryPos_u,simHitEntryPos_v,simHitEntryPos_w,
  // simHitPath_u,simHitPath_v,simHitPath_w,
  // pathLengthGeant4,pathLength,chargeDeposition,debugFlag

  m_debugCsvFile << std::to_string(hitInfo.eventNumber()) << "," << hitInfo.layerIndex() << "," << hitInfo.nSegments() << ","  << m_sensorThickness << "," << i_u << "," << i_v << ",";

  m_debugCsvFile << hitPos.local.x() << "," << hitPos.local.y() << "," << hitPos.local.z() << ",";
  m_debugCsvFile << hitPos.entry.x() << "," << hitPos.entry.y() << "," << hitPos.entry.z() << ",";
  m_debugCsvFile << hitPos.path.x() << "," << hitPos.path.y() << "," << hitPos.path.z() << ",";
  m_debugCsvFile << hitInfo.simPathLength() << "," << hitPos.path.r() << "," << hitInfo.charge() <<  "," << hitInfo.debugFlag() << "\n";
  m_debugCsvFile.flush();
  verbose() << "Wrote simHit with event " << hitInfo.eventNumber() << ", layerIndex " << hitInfo.layerIndex() << ", segmentN " << hitInfo.nSegments() << ", i_u " << i_u << ", i_v " << i_v << " to debug CSV file." << endmsg;
  return;
}
