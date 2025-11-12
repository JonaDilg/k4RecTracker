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
 *      -> generally use dd4hep::rec::Vector3D, convert via convertVector() where edm4hep::Vector3d is needed
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

  initializeServicesAndGeometry();

  checkGaudiProperties();
  
  /* Sensor parameters are given as either
  *  - a single value, which is then applied to all layers
  *    -> rewrite the vector to contain the same value for all layers (simplifies code later)
  *  - a vector of values, one per layer (for each layer that will be digitized)
  *    -> check that the number of entries matches the number of layers in the geometry
  * Need to create a layer index mapping as well, to map from layer number (from cellID) to index in the sensor parameter vectors (to access the correct sensor parameters for each layer). */ 
  setupSensorParameters();

  loadKernels(); // also sets m_pixelPitch_u, m_pixelPitch_v, m_sensorThickness based on kernel file header. this needs to be updated if different sensor types are used across layers
  
  verbose() << " - Using the following sensor parameters for each layer:" << endmsg;
  for (const auto& [layer, index] : m_layerToIndex) {
    verbose() << "   - Layer " << layer << " (index " << index << "): pitch_u = " << m_pixelPitch_u.at(index) << " mm, pitch_v = " << m_pixelPitch_v.at(index) << " mm, count_u = " << m_pixelCount_u.value().at(index) << ", count_v = " << m_pixelCount_v.value().at(index) << ", thickness = " << m_sensorThickness.at(index) << " mm" << endmsg;
  }

  if (m_debugHistograms.value())
    setupDebugHistograms();

  if (!m_debugCsvFileName.value().empty())
    setupDebugCsvOutput();

  /* TODO: check that pixel pitch and count match sensor size in geometry (I am not sure how to get the sensor size from the geometry though, does seem to be a bit more general, subDetectors have children that might or might not be layers). this is currently done for every event, but could be done in the initialization phase, which is of course far more efficient ~ Jona 2025-09
   *  TODO load lookup table for charge transport and diffusion
   *  TODO check that k4geo has same pitch as lookup table */
  info() << " - Initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode VTXdigi_Allpix2::finalize() {
  info() << "FINALIZING ..." << endmsg;

  if (m_debugCsvFile.is_open()) {
    m_debugCsvFile.close();
    verbose() << " - Closed debug CSV file" << endmsg;
  }

  printCountersSummary();

  verbose() << " - finalized successfully" << endmsg;
  return StatusCode::SUCCESS;
} 


/* -- event loop -- */

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> VTXdigi_Allpix2::operator()
  (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {

  // Initial check: returns false if event has no simHits, throws if there is a problem with the setup
  if (!CheckInitialSetup(simHits, headers))
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());

  // output collections
  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitsLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();

  // initialize random number generator for this event
  const u_int32_t rngSeed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  TRandom2 rngEngine = TRandom2(rngSeed);
  verbose() << " - RNG engine initialized, using seed " << rngSeed << endmsg;

  // loop over sim hits, digitize them, create output collections
  for (const auto& simHit : simHits) {
    ++m_counter_simHitsRead;

    // check if the simHit is on a relevant layer. Needs to be done here ( before GatherHitInfo() ), to avoid out-of-bounds access in Gaudi properties
    if (!CheckLayerCut(simHit))
      continue;

    HitInfo hitInfo; // contains all relevant quantities of the simHit
    HitPosition hitPos; // contains all relevant positions related to the simHit
    std::tie(hitInfo, hitPos) = GatherHitInfoAndPosition(simHit, headers);

    { // debug statements
      debug() << "   - Processing simHit in event " << hitInfo.eventNumber() << ", layer " << hitInfo.layerIndex() << ", cellID " << hitInfo.cellID() << " (" << std::bitset<24>(hitInfo.cellID()) << ")" << endmsg;
      debug() << "   - Momentum = " << hitInfo.simMomentum() << " GeV/c, dep. charge = " << hitInfo.charge() << " e-, path length = " << hitPos.path.r() << " mm, Geant4 path length = " << hitInfo.simPathLength() << " mm" << endmsg;
      verbose() << "   - Position (global) " << hitPos.global.x() << " mm, " << hitPos.global.y() << " mm, " << hitPos.global.z() << " mm" << endmsg;
      verbose() << "   - Position (local) " << hitPos.local.x() << " mm, " << hitPos.local.y() << " mm, " << hitPos.local.z() << " mm (in sensor frame)" << endmsg;
      verbose() << "   - Distance to surface " << hitInfo.simSurface()->distance(dd4hep::mm * hitPos.global) << " mm" << endmsg; // dd4hep expects cm, simHitGlobalPosition is in mm. dd4hep::cm=0.1
      verbose() << "   - Entry point (sensor local) : " << hitPos.entry.x() << " mm, " << hitPos.entry.y() << " mm, " << hitPos.entry.z() << " mm" << endmsg;
    }

    if (!CheckSimHitCuts(hitInfo, hitPos))
      continue;
    ++m_counter_simHitsAccepted;
    
    // Loop over path segments, share each segments charge to pixels
    PixelChargeMatrix pixelChargeMatrix = DepositAndCollectCharge(hitInfo, hitPos);
    
    // Generate noise for each pixel (only now, after we know which pixels are fired)
    pixelChargeMatrix.GenerateNoise(rngEngine, m_electronicNoise); // only generate noise for pixels, after we know how large the matrix is in the end

    // find pixels with charge above threshold, create digiHits
    AnalyseSharedCharge(hitInfo, hitPos, pixelChargeMatrix, simHit, digiHits, digiHitsLinks);

    if (m_debugHistograms)
      FillGeneralDebugHistograms(hitInfo, hitPos, pixelChargeMatrix);
      
    if (m_debugCsv) {
      int i_u_debug, i_v_debug;
      std::tie(i_u_debug, i_v_debug) = computePixelIndices(hitPos.local, hitInfo.layerIndex(), hitInfo.length(0), hitInfo.length(1));
      appendSimHitToCsv(hitInfo, hitPos, i_u_debug, i_v_debug);
    }
  } // end loop over sim hits
  
  debug() << "FINISHED event." << endmsg;
  verbose() << " - Returning collections: digiHits.size()=" << digiHits.size() << ", digiHitsLinks.size()=" << digiHitsLinks.size() << endmsg;
  
  return std::make_tuple(std::move(digiHits), std::move(digiHitsLinks));
} // operator()


/* -- Initialization / finalization functions -- */

void VTXdigi_Allpix2::initializeServicesAndGeometry() {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc)
    throw GaudiException("Unable to get UniqueIDGenSvc", "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);

  m_geometryService = serviceLocator()->service(m_geometryServiceName);
  if (!m_geometryService)
    throw GaudiException("Unable to retrieve the GeoSvc", "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  std::string cellIDstr = m_geometryService->constantAsString(m_encodingStringVariable.value());
  m_cellIDdecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDstr);
  if (!m_cellIDdecoder)
    throw GaudiException("Unable to retrieve the cellID decoder", "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  const dd4hep::Detector* detector = m_geometryService->getDetector();
  if (!detector)
    throw GaudiException("Unable to retrieve the DD4hep detector from GeoSvc", "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  const dd4hep::rec::SurfaceManager* simSurfaceManager = detector->extension<dd4hep::rec::SurfaceManager>();
  if (!simSurfaceManager)
    throw GaudiException("Unable to retrieve the SurfaceManager from the DD4hep detector", "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);
  
  m_simSurfaceMap = simSurfaceManager->map(m_subDetName.value());
  if (!m_simSurfaceMap)
    throw GaudiException("Unable to retrieve the simSurface map for subdetector " + m_subDetName.value(), "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);

  m_volumeManager = detector->volumeManager();
  if (!m_volumeManager.isValid())
    throw GaudiException("Unable to retrieve the VolumeManager from the DD4hep detector", "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);

  /* subDetector not needed as of now. Keep for future reference ~ Jona 2025-10 */
  // const dd4hep::DetElement subDetector = detector->detector(m_subDetName.value());
  // if (!subDetector.isValid())
  //   throw GaudiException("Unable to retrieve the DetElement for subdetector " + m_subDetName.value(), "VTXdigi_Allpix2::initializeServicesAndGeometry()", StatusCode::FAILURE);

  debug() << " - Successfully retrieved all necessary services and detector elements, starting to check Gaudi properties." << endmsg;
}

void VTXdigi_Allpix2::checkGaudiProperties() {
  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);

  if (m_pixelCount_u.value().empty())
    throw GaudiException("Property PixelCount_u is not set! Give either a single value \"[<value>]\" or one per layer.", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);
  if (m_pixelCount_v.value().empty())
    throw GaudiException("Property PixelCount_v is not set! Give either a single value \"[<value>]\" or one per layer.", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);

  if (m_maxClusterSize.value().size() != 2) 
    throw GaudiException("Property MaximumClusterSize must have exactly 2 entries (u and v).", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);
  for (int i = 0; i < 2; ++i) {
    if (m_maxClusterSize.value().at(i) < 1) 
      throw GaudiException("Property MaximumClusterSize entries must be positive odd numbers.", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);
    if ((m_maxClusterSize.value().at(i)-1) % 2 != 0) 
      throw GaudiException("Property MaximumClusterSize entries must be positive odd numbers.", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);
  }

  if (m_globalKernel.value().size() > 0 && m_kernelFileName.value().length() > 0)
    throw GaudiException("Please provide either a global kernel or a kernel file, not both!", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);
  else if (m_globalKernel.value().size() == 0 && m_kernelFileName.value().length() == 0)
    throw GaudiException("No charge sharing kernel specified! Please provide either a global kernel or a file containing charge sharing kernels.", "VTXdigi_Allpix2::checkGaudiProperties()", StatusCode::FAILURE);
}

void VTXdigi_Allpix2::setupSensorParameters() {

  if (m_layersToDigitize.value().empty()) {
    verbose() << " - No layers specified to be digitized, will digitize all layers found in geometry." << endmsg;

    /* TODO: Find all layers in geometry & set m_layerCount accordingly. right now this doesnt work.*/
    throw GaudiException("Finding all layers in geometry is not implemented yet. Please specify the layers to be digitized via the LayersToDigitize property.", "VTXdigi_Allpix2::setupSensorParameters()", StatusCode::FAILURE);

    m_layerCount = 5; // placeholder, needs to be set properly when finding all layers in geometry is implemented
  }
  else { // if layers are specified as Gaudi property
    m_layerCount = m_layersToDigitize.size();
    verbose() << " - Digitizing " << m_layerCount << " layers as specified in LayersToDigitize property." << endmsg;

    /* create layer number to internal index map
     * needed in case the user specifies non-consecutive layer numbers */
    m_layerToIndex.clear();
    for (size_t layerIndex = 0; layerIndex < m_layersToDigitize.size(); layerIndex++) {
      int layer = m_layersToDigitize.value().at(layerIndex);
      m_layerToIndex.insert({layer, static_cast<int>(layerIndex)});
    }
    debug() << " - Digitizing the following layers: " << endmsg;
    for (const auto& [layer, index] : m_layerToIndex) {
      debug() << "   - Layer " << layer << " (index " << index << ")" << endmsg;
    }

    /* check sensor & pixel parameters */
    if (m_pixelCount_u.value().size() == 1 && m_pixelCount_v.value().size()) {
      // single value for all layers, rewrite to vector of correct size
      m_pixelCount_u.value().resize(m_layerCount, m_pixelCount_u.value().at(0));
      m_pixelCount_v.value().resize(m_layerCount, m_pixelCount_v.value().at(0));
      verbose() << " - found single values for pixel counts, applying them to all layers" << endmsg;
    }
    else if (
      m_pixelCount_u.value().size() == static_cast<long unsigned int>(m_layerCount) && 
      m_pixelCount_v.value().size() == static_cast<long unsigned int>(m_layerCount)) {
        // one entry per layer
        verbose() << " - found pixel count values for all " << m_layerCount << " layers" << endmsg;
      }
    else {
      throw GaudiException("Pixel counts must either have a single value (that is then applied to all layers) or one value per layer to be digitized.", "VTXdigi_Allpix2::setupSensorParameters()", StatusCode::FAILURE);
    }
  }
}

void VTXdigi_Allpix2::loadKernels() {
  verbose() << " - Importing charge sharing kernels..." << endmsg;
  
  /* sanity check */
  if (m_globalKernel.value().empty() && m_kernelFileName.value().empty())
    throw GaudiException("No charge sharing kernel specified! Please provide either a global kernel or a kernel file.", "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);

  /* load kernels */
  if (!m_globalKernel.value().empty()) { // use global kernel (mostly for debugging)
    debug() << " - Global kernel specified in <config>, applying this for all layers" << endmsg;

    for (int i=0; i<3; i++)
      m_inPixelBinCount[i] = 10;

    m_kernelSize = static_cast<int>(std::sqrt(m_globalKernel.value().size()));
    verbose() << " -   Using kernel size of " << m_kernelSize << " (from global kernel with " << m_globalKernel.value().size() << " entries)" << endmsg;

    m_chargeSharingKernels = std::make_unique<ChargeSharingKernels>(m_inPixelBinCount[0], m_inPixelBinCount[1], m_inPixelBinCount[2], m_kernelSize);

    /* kernel size is checked in ChargeSharingKernel::SetAllKernels() */
    m_chargeSharingKernels->SetAllKernels(m_globalKernel.value());
  }
  else { // load from file
    debug() << " - Loading kernels from file: " << m_kernelFileName.value() << endmsg;
    /* This implements parsing the default Allpix2 kernel file format. A general version was implemented in a previous commit (2025-11), but removed because the amount of options made it unneccessarily hard to validate and use. 
     * See https://indico.cern.ch/event/1489052/contributions/6475539/attachments/3063712/5418424/Allpix_workshop_Lemoine.pdf (slide 10) for more info on fields in the kernel file */

    const int headerLines = 5; // allpix2 kernel files have 5 header lines, and then a kernel per line

    info() << " -   Opening kernel file: " << m_kernelFileName.value() << ". Expecting Allpix2 format, defined in ChargePropagationWriter module." << endmsg;
    std::ifstream kernelFile(m_kernelFileName.value());
    if (!kernelFile.is_open())
      throw GaudiException("Could not open kernel file: " + m_kernelFileName.value(), "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);

    std::string line;
    int lineNumber = 0; // next line to be read (0-indexed)
    float kernelsSum = 0.f;

    /* loading pixel-pitch, thickness and in-pixel bin counts from header (all in 5th line) */
    for (; lineNumber < 5; ++lineNumber)
      std::getline(kernelFile, line);

    std::istringstream headerStringStream(line);
    std::string headerEntry;
    std::vector<std::string> headerLineEntries;

    while (std::getline(headerStringStream, headerEntry, ' ')) {
      if (!headerEntry.empty())
        headerLineEntries.push_back(headerEntry);
    }

    if (headerLineEntries.size() != 11)
      throw GaudiException("Invalid number of entries in kernel file at header line 5: found " + std::to_string(headerLineEntries.size()) + " entries, expected 11.", "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);

    
    
    float thickness = std::stof(headerLineEntries.at(0)) / 1000.f; // convert from um to mm
    float pitch_u = std::stof(headerLineEntries.at(1)) / 1000.f;
    float pitch_v = std::stof(headerLineEntries.at(2)) / 1000.f;
    for (int i=0; i<3; i++)
      m_inPixelBinCount[i] = std::stoi(headerLineEntries.at(7+i));

    /* Initialize sensor properties. This implementation needs to be changed in case of different sensor types across layers */
    m_pixelPitch_u.resize(m_layerCount, pitch_u);
    m_pixelPitch_v.resize(m_layerCount, pitch_v);
    m_sensorThickness.resize(m_layerCount, thickness);

    verbose() << " -   found pixel-pitch of (" << pitch_u << " x " << pitch_v << ") mm2, thickness of " << thickness << " mm, and in-pixel bin count of (" << m_inPixelBinCount[0] << ", " << m_inPixelBinCount[1] << ", " << m_inPixelBinCount[2] << ")." << endmsg;

    /* get the kernel size (5x5, 7x7, ...) from the length of the first line after the header */

    if (std::getline(kernelFile, line)) {
      m_kernelSize = static_cast<int>(std::sqrt(std::count(line.begin(), line.end(), ' ') - 2)); // not very robust, but works for valid Allpix2 files. first 3 entries are bin indices
    }
    else {
      throw GaudiException("Could not read first kernel line after header in kernel file: " + m_kernelFileName.value(), "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);
    }
    verbose() << " - Detected kernel size of " << m_kernelSize << " from first kernel line." << endmsg;

    /* set up mapping from Allpix2 kernel format
    *   (row-major, starts on bottom left)
    * to the format expected by the ChargeSharingKernels class 
    *   (row-major, starts on top-left) */
    std::vector<int> valueIndices(m_kernelSize*m_kernelSize, 0); // i: index in local format; valueIndices[i]: index in Allpix2 format
    for (int i_u = 0; i_u < m_kernelSize; i_u++) {
      for (int i_v = 0; i_v < m_kernelSize; i_v++) {
        int i_allpix2 = i_u + (m_kernelSize - 1 - i_v) * m_kernelSize;
        int i_local = i_u + i_v * m_kernelSize;
        valueIndices[i_local] = i_allpix2;
      }
    }

    /* reset the file and go to the beginning of the kernel data (after header) */
    kernelFile.clear();
    kernelFile.seekg(0, std::ios::beg);
    for (int i=0; i<headerLines; ++i)
      std::getline(kernelFile, line);
    
    /* loop over lines that contain a kernel each, set the kernels */
    int kernelSize = 0;
    m_chargeSharingKernels = std::make_unique<ChargeSharingKernels>(m_inPixelBinCount[0], m_inPixelBinCount[1], m_inPixelBinCount[2], m_kernelSize);

    debug() << " -   Loading kernels from file lines ..." << endmsg;
    while (std::getline(kernelFile, line)) {
      if (line.empty() || line[0] == '#')
        throw GaudiException("Empty or comment line found in kernel file at line " + std::to_string(lineNumber+1) + ". All lines (after 5 header lines) must contain valid kernel data.", "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);
      
      std::istringstream stringStream(line);
      std::vector<std::string> lineEntries;
      std::string entryString;

      /* read the line & do sanity checks */
      while (std::getline(stringStream, entryString, ' ')) {
        if (!entryString.empty())
          lineEntries.push_back(entryString);
      }

      if (kernelSize == 0)
        kernelSize = static_cast<int>(std::sqrt(lineEntries.size() - 3)); // first 3 entries are bin indices
      else if (kernelSize != static_cast<int>(std::sqrt(lineEntries.size() - 3)))
        throw GaudiException("Invalid number of entries in kernel file at line " + std::to_string(lineNumber+1) + ": found " + std::to_string(lineEntries.size()) + " entries, expected " + std::to_string(3 + m_kernelSize * m_kernelSize) + " = 3 indices + kernelSize^2.", "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);

      /* in-pixel bin of this kernel, given in first 3 columns */
      int j_u = std::stoi(lineEntries[0]) - 1; // Allpix2 input is 1-indexed, sane people use 0-indexing
      int j_v = std::stoi(lineEntries[1]) - 1;
      int j_w = std::stoi(lineEntries[2]) - 1;

      if (j_u < 0 || j_u >= m_inPixelBinCount[0] ||
          j_v < 0 || j_v >= m_inPixelBinCount[1] ||
          j_w < 0 || j_w >= m_inPixelBinCount[2]) {
        throw GaudiException("Invalid in-pixel bin indices in kernel file at line " + std::to_string(lineNumber+1) + ": got (" + std::to_string(j_u) + ", " + std::to_string(j_v) + ", " + std::to_string(j_w) + "), but expected ranges are [0, " + std::to_string(m_inPixelBinCount[0]-1) + "], [0, " + std::to_string(m_inPixelBinCount[1]-1) + "], [0, " + std::to_string(m_inPixelBinCount[2]-1) + "].", "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);
      }

      /* Parse kernel values & do sanity checks */
      std::vector<float> kernelValues(m_kernelSize*m_kernelSize, 0.);
      float kernelValueSum = 0.f;
      for (int i = 0; i < m_kernelSize*m_kernelSize; i++) {
        float entry = std::stof(lineEntries[3 + valueIndices[i]]); // NaN check done on sum
        kernelValues[i] = entry;
        kernelValueSum += entry;
      } 

      if (std::isnan(kernelValueSum))
        throw GaudiException("NaN encountered in kernel for in-pixel bin (" + std::to_string(j_u) + ", " + std::to_string(j_v) + ", " + std::to_string(j_w) + ") at line " + std::to_string(lineNumber) + " in kernel file.", "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);

      if (kernelValueSum < 0.f || kernelValueSum > 1.f + m_numericLimit_float) {
        warning() << "Kernel for in-pixel bin (" << j_u << ", " << j_v << ", " << j_w << ") has invalid sum (" << kernelValueSum << ", must be in [0, 1]). Normalising to 1." << endmsg;
        for (auto& value : kernelValues)
          value /= kernelValueSum;
      }

      /* finalize this line, actually set the kernel */
      debug() << " -     Parsed kernel for in-pixel bin (" << j_u << ", " << j_v << ", " << j_w << ") with entry sum " << std::to_string(kernelValueSum) << ", setting it now..." << endmsg;
      kernelsSum += kernelValueSum;
      m_chargeSharingKernels->SetKernel(j_u, j_v, j_w, kernelValues);
      lineNumber++;

    } // loop over lines with kernels

    if (lineNumber - headerLines != m_inPixelBinCount[0] * m_inPixelBinCount[1] * m_inPixelBinCount[2])
      throw GaudiException("Invalid number of kernels loaded from file: expected " + std::to_string(m_inPixelBinCount[0] * m_inPixelBinCount[1] * m_inPixelBinCount[2]) + " kernels (inferred from InPixelBinCount = [" + std::to_string(m_inPixelBinCount[0]) + ", " + std::to_string(m_inPixelBinCount[1]) + ", " + std::to_string(m_inPixelBinCount[2]) + "]), but got " + std::to_string(lineNumber - headerLines) + ".", "VTXdigi_Allpix2::loadKernels()", StatusCode::FAILURE);

    kernelsSum /= static_cast<float>(lineNumber - headerLines);
    info() << " - Loaded " << (lineNumber - headerLines) << " kernels from file. " << kernelsSum*100 << " percent of charge deposited in the sensor volume is collected by the pixels (the rest is lost, eg. due to being outside of depletion or due to trapping)." << endmsg;
  } // if load from file
}

void VTXdigi_Allpix2::setupDebugHistograms() {
  error() << " - You enabled creating debug histograms by setting `DebugHistograms = True`. This is NOT MULTITHREADING SAFE and will cause crashes if multithreading is used." << endmsg;
  verbose () << " - Creating debug histograms ..." << endmsg;

  /* -- Global Histograms (collect from all layers) -- */

  m_histogramsGlobal.at(histGlobal_simHitE).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, 
      "simHitE", "SimHit Energy;keV", 
      {1000, 0, m_sensorThickness.at(0)*2000}}); 
  m_histogramsGlobal.at(histGlobal_simHitCharge).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this,
      "simHitCharge", "SimHit Raw Charge;dep. charges [e-]",
      {1000, 0, m_sensorThickness.at(0)*500000}});

  m_histogramsGlobal.at(histGlobal_clusterSize_raw).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, 
      "ClusterSize_Raw", 
      "Global Cluster Size of Raw Charges (ie. number of pixels that receive any charge from this simHit);Cluster size [pix]",
      {50, -0.5, 49.5}
    }
  );

  m_histogramsGlobal.at(histGlobal_clusterSize_measured).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this,
      "ClusterSize_Measured",
      "Global Cluster Size of Measured Charges (ie. number of pixels that are above threshold due to the simHit's deposited charge);Cluster size [pix]",
      {50, -0.5, 49.5}
    }
  );

  m_histogramsGlobal.at(histGlobal_pathLength).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "PathLength", "Path Length in Sensor Active Volume;Path length [um]", {500, 0., m_sensorThickness.at(0)*10*1000}}); // in um, max 10x sensor thickness (for very shallow angles)
  m_histogramsGlobal.at(histGlobal_pathLengthGeant4).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "PathLength-Geant4", "Path Length in Sensor Active Volume\nAs Given by Geant4;Path length [um]", {500, 0., m_sensorThickness.at(0)*10*1000}}); // in um, max 10x sensor thickness (for very shallow angles)

  m_histogramsGlobal.at(histGlobal_EntryPointX).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "EntryPointX", "SimHit Entry Point U (in local sensor frame);u [mm]", {400, -4, 4}});
  m_histogramsGlobal.at(histGlobal_EntryPointY).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "EntryPointY", "SimHit Entry Point V (in local sensor frame);v [mm]", {2000, -20, 20}});
  m_histogramsGlobal.at(histGlobal_EntryPointZ).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "EntryPointZ", "SimHit Entry Point W (in local sensor frame);w [mm]", {1000, -300, 300}});

  m_histogramsGlobal.at(histGlobal_DisplacementU).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "LocalDisplacementU", "Displacement U (local sensor frame): digiHit_u - simHit_u;dU [um]", {800, -200, 200}});
  m_histogramsGlobal.at(histGlobal_DisplacementV).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "LocalDisplacementV", "Displacement V (local sensor frame): digiHit_v - simHit_v;dV [um]", {800, -200, 200}});
  m_histogramsGlobal.at(histGlobal_DisplacementR).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "GlobalDisplacementR", "Displacement R (global frame): | digiHit - simHit |;dR [um]", {300, 0., 300}});

  m_histogramsGlobal.at(histGlobal_chargeCollectionEfficiency_raw).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "ChargeCollectionEfficiency_rawCharge", "Charge collection efficiency (raw charge, eg. before noise & threshold);# e- (digitised) / # e- (simHit)", {500, 0., 2.}});
  
  m_histogramsGlobal.at(histGlobal_chargeCollectionEfficiency).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "ChargeCollectionEfficiency", "Charge collection efficiency;# e- (digitised) / # e- (simHit)", {500, 0., 2.}});

  

  m_histogramsGlobal.at(histGlobal_pixelChargeMatrix_size_u).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "PixelChargeMatrix_Size_u", "Final Size of the Pixel Charge Matrix in u (local);pixel charge matrix size u [pix]", {100, -0.5, 99.5}});
  m_histogramsGlobal.at(histGlobal_pixelChargeMatrix_size_v).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, "PixelChargeMatrix_Size_v", "Final Size of the Pixel Charge Matrix in v (local);pixel charge matrix size v [pix]", {100, -0.5, 99.5}});


  
  /* -- Per-layer Histograms -- */
  
  m_histograms1d.resize(m_layersToDigitize.size());
  m_histograms2d.resize(m_layersToDigitize.size()); // resize vector to hold all histograms
  m_histWeighted2d.resize(m_layersToDigitize.size());

  for (int layer : m_layersToDigitize) {
    int layerIndex = m_layerToIndex.at(layer);

    m_histograms1d.at(layerIndex).at(hist1d_ClusterSize_raw).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_ClusterSize_Raw",
        "Cluster Size of Raw Charges (ie. number of pixels that receive any charge from this simHit);Cluster size [pix]",
        {50, -0.5, 49.5}
      }
    );
    m_histograms1d.at(layerIndex).at(hist1d_ClusterSize_measured).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_ClusterSize_Measured",
        "Cluster Size of Measured Charges (ie. number of pixels that are above threshold due to the simHit's deposited charge);Cluster size [pix]",
        {50, -0.5, 49.5}
      }
    );

    m_histograms1d.at(layerIndex).at(hist1d_IncidentAngle_ThetaLocal).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_IncidentAngle_LocalTheta",
        "Local theta: particle momentum incident angle to sensor normal;Polar angle [deg]",
        {360, -180., 180.}
      }
    );
    m_histograms1d.at(layerIndex).at(hist1d_IncidentAngle_PhiLocal).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_IncidentAngle_LocalPhi",
        "Local phi: particle momentum incident angle to sensor normal;Azimuthal angle [deg]",
        {360, -180., 180.}
      }
    );

    m_histograms1d.at(layerIndex).at(hist1d_SimHitMomentum).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_simHitMomentum",
        "SimHit Momentum;Momentum [MeV/c]",
        {10000, 0., 5000.}
      }
    );

    verbose () << " - Creating 2D histograms for layer " << layer << " (index " << layerIndex << ")" << endmsg;

    m_histograms2d.at(layerIndex).at(hist2d_hitMap_simHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_HitMap_simHits",
        "SimHit Hitmap, Layer " + std::to_string(layer) + ";u [pix]; v [pix]",
        {static_cast<unsigned int>(m_pixelCount_u.value().at(0)), -0.5, static_cast<double>(m_pixelCount_u.value().at(0)+0.5)},
        {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)}
      }
    );
    m_histograms2d.at(layerIndex).at(hist2d_hitMap_digiHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_HitMap_digiHits",
        "DigiHit Hitmap, Layer " + std::to_string(layer) + ";u [pix]; v [pix]",
        {static_cast<unsigned int>(m_pixelCount_u.value().at(0)), -0.5, static_cast<double>(m_pixelCount_u.value().at(0)+0.5)},
        {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)}
      }
    );
    m_histograms2d.at(layerIndex).at(hist2d_pathLength_vs_simHit_v).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_TrackLength_vs_simHit_v",
        "Track Length in Sensor Active Volume vs. SimHit v position (local), Layer " + std::to_string(layer) + ";SimHit v [pix];Track length [um]",
          {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)},
        {200, 0., 10*m_sensorThickness.at(layerIndex)*1000}, // in um
      
      }
    );
    m_histograms2d.at(layerIndex).at(hist2d_pixelChargeMatrixSize).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_PixelChargeMatrixSize",
        "Pixel Charge Matrix Size, Layer " + std::to_string(layer) + ";u [pix];v [pix]",
        {100, -0.5, 99.5},
        {100, -0.5, 99.5}
      }
    );
    m_histograms2d.at(layerIndex).at(hist2d_IncidentAngle).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full>{this, 
        "Layer"+std::to_string(layer)+"_IncidentAngle_2D",
        "Incident particle angle to sensor normal, Layer " + std::to_string(layer) + ";Local theta [deg];Local phi [deg]",
        {360, -180., 180.},
        {360, -180., 180.}
      }
    );
    
    m_histWeighted2d.at(layerIndex).at(histWeighted2d_averageCluster).reset(
      new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, double>{this,
        "Layer"+std::to_string(layer)+"_AverageCluster",
        "Average Cluster Shape, after applying Threshold & Noise, Layer " + std::to_string(layer) + "(charge per hit normalised to 1);u [pix]; v [pix]",
        {static_cast<unsigned int>(m_maxClusterSize.value().at(0)),
          -static_cast<double>(m_maxClusterSize.value().at(0))/2,
          static_cast<double>(m_maxClusterSize.value().at(0))/2},
        {static_cast<unsigned int>(m_maxClusterSize.value().at(1)),
          -static_cast<double>(m_maxClusterSize.value().at(1))/2,
          static_cast<double>(m_maxClusterSize.value().at(1))/2}
      }
    );
    m_histWeighted2d.at(layerIndex).at(histWeighted2d_chargeOriginU).reset(
      new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, double>{this,
        "Layer"+std::to_string(layer)+"_ChargeOriginU",
        "Charge collected from in-pix bins in u, Layer " + std::to_string(layer) + ";u pos. relative to collecting pixel center [pix]; w pos. [um]",
        {static_cast<unsigned int>(m_kernelSize*m_inPixelBinCount[0]),
          -static_cast<double>(m_kernelSize)/2,
          static_cast<double>(m_kernelSize)/2},
        {static_cast<unsigned int>(m_inPixelBinCount[2]),
          -m_sensorThickness.at(layerIndex)*1000/2, 
          m_sensorThickness.at(layerIndex)*1000/2}
      }
    );
    m_histWeighted2d.at(layerIndex).at(histWeighted2d_chargeOriginV).reset(
      new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, double>{this,
        "Layer"+std::to_string(layer)+"_ChargeOriginV",
        "Charge collected from in-pix bins in v, Layer " + std::to_string(layer) + ";v pos. relative to collecting pixel center [pix]; w pos. [um]",
        {static_cast<unsigned int>(m_kernelSize*m_inPixelBinCount[1]),
          -static_cast<double>(m_kernelSize)/2,
          static_cast<double>(m_kernelSize)/2},
        {static_cast<unsigned int>(m_inPixelBinCount[2]),
          -m_sensorThickness.at(layerIndex)*1000/2, 
          m_sensorThickness.at(layerIndex)*1000/2}
      }
    );
  }
}

void VTXdigi_Allpix2::setupDebugCsvOutput() {
  m_debugCsv = true;

  error() << " - You enabled the CSV output by setting `DebugCsvFileName` to a path. This is NOT MULTITHREADING SAFE and will cause crashes if multithreading is used." << endmsg;
  m_debugCsvFile.open(m_debugCsvFileName.value());

  if (!m_debugCsvFile.is_open()) {
    throw GaudiException("Failed to open debug CSV file", "VTXdigi_Allpix2::initialize", StatusCode::FAILURE);
  } else {
    m_debugCsvFile << "eventNumber,layerIndex,segmentCount,sensorThickness,pix_u,pix_v,simHitPos_u,simHitPos_v,simHitPos_w,simHitEntryPos_u,simHitEntryPos_v,simHitEntryPos_w,simHitPath_u,simHitPath_v,simHitPath_w,pathLengthGeant4,pathLength,rawChargeDeposition,debugFlag\n";
    m_debugCsvFile.flush();
    debug() << "   - writing to file: " << m_debugCsvFileName.value() << endmsg;
  }
}

void VTXdigi_Allpix2::printCountersSummary() const {
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

bool VTXdigi_Allpix2::CheckInitialSetup(const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {
  info() << "PROCESSING event. run " << headers.at(0).getRunNumber() << "; event " << headers.at(0).getEventNumber() << "; with " << simHits.size() << " simHits" << endmsg;
  ++m_counter_eventsRead;

  // early sanity checks to avoid segfaults from null pointers
  if (!m_chargeSharingKernels)
    throw GaudiException("ChargeSharingKernels is null in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::CheckInitialSetup()", StatusCode::FAILURE);

  if (!m_simSurfaceMap)
    throw GaudiException("SimSurfaceMap is null in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::CheckInitialSetup()", StatusCode::FAILURE);

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

std::tuple<VTXdigi_Allpix2::HitInfo, VTXdigi_Allpix2::HitPosition> VTXdigi_Allpix2::GatherHitInfoAndPosition(const edm4hep::SimTrackerHit& simHit, const edm4hep::EventHeaderCollection& headers) const {
  HitInfo hitInfo(*this, simHit, headers);

  HitPosition hitPos;
  hitPos.global = convertVector(simHit.getPosition()); // global simHit position (from Geant4)
  hitPos.local = transformGlobalToLocal(hitPos.global, hitInfo.cellID()); // the same position, but in the local sensor frame (u,v,w)
  std::tie(hitPos.entry, hitPos.path) = constructSimHitPath(hitInfo, hitPos, simHit); // simHit path through the sensor
  // -> hitPos is now fully defined

  hitInfo.setNSegments(std::max(1, int( hitPos.path.r() / m_targetPathSegmentLength ))); // both in mm
  // -> with this, hitInfo is now fully defined as well

  return std::make_tuple(hitInfo, hitPos);
}

bool VTXdigi_Allpix2::CheckSimHitCuts (const HitInfo& hitInfo, const HitPosition& hitPos) const {

  // DISMISS if entry point is outside sensor
  if (m_cutPathOutsideSensor) {

    if (abs(hitPos.entry.z()) > hitInfo.thickness() / 2 + m_numericLimit_float) { // entry point is outside sensor thickness
      verbose() << "   - DISMISSED simHit (entry point is outside sensor thickness (local w = " << hitPos.entry.z()*1000 << " um, sensor thickness = " << hitInfo.thickness()*1000 << " um)." << endmsg;
      ++m_counter_simHitsRejected_OutsideSensor;
      return false;
    }

    if (abs(hitPos.entry.z()+hitPos.path.z()) > hitInfo.thickness() / 2 + m_numericLimit_float) { // exit point is outside sensor thickness
      verbose() << "   - DISMISSED simHit (exit point is outside sensor thickness (local w = " << (hitPos.entry.z()+hitPos.path.z())*1000 << " um, sensor thickness = " << hitInfo.thickness()*1000 << " um)." << endmsg;
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

std::tuple<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D> VTXdigi_Allpix2::constructSimHitPath (HitInfo& hitInfo, HitPosition& hitPos, const edm4hep::SimTrackerHit& simHit) const {
  /* Get the simHitPath (the path of the simulated particle through a sensor) of a given simHit. Describet by the path direction (unit vector) and entry point (local frame).
   *  
   * The way I do this differs from what Jessy does in their vtxDigi Detailed.
   *  My approach ensures that the simHitPath always extends exactly from one sensor surface to the other, and that the simHitPos always lies on this path. 
   * Then, the paths are shortened to the PathLength that Geant4 gives us, if they are longer (the Geant4 pathLength is always longer than the linear approx., because it accounts for B-field and multiple scattering inside the sensor). We need to do this to avoid unphysically long paths, where the charge-per-pathLength would be unphysically tiny
   *  I am not sure if this is better or worse, it definitely produces less problems. 
   *  I think the assumption of a linear path through the sensor is not compatible with using the path-length as normalisation for the path-vector
   */
  
  // get simHitMomentum (this vector gives the exact agles of the particle's path through the sensor. We assume that path is linear.). Transform to sensor's local coordinates.
  double simHitGlobalMomentum_double[3] = {simHit.getMomentum().x * dd4hep::GeV, simHit.getMomentum().y * dd4hep::GeV, simHit.getMomentum().z * dd4hep::GeV}; // need floats for the TransfomationMatrix functions

  TGeoHMatrix transformationMatrix = computeTransformationMatrix(hitInfo.cellID());

  double simHitLocalMomentum_double[3] = {0.0, 0.0, 0.0};
  transformationMatrix.MasterToLocalVect(simHitGlobalMomentum_double, simHitLocalMomentum_double);

  dd4hep::rec::Vector3D simHitPath(
    simHitLocalMomentum_double[0],
    simHitLocalMomentum_double[1],
    simHitLocalMomentum_double[2]);

  float scaleFactor = hitInfo.thickness() / std::abs(simHitLocalMomentum_double[2]);
  simHitPath = scaleFactor * simHitPath ; // now, simHitPath extends for one sensor surface to the other surface
  // edge case with curlers having horizontal momentum is handled later, by cutting the path to the length supplied by Geant4.
  
  // calculate the path's entry position into the sensor, by placing it such that the path passes through the simHit position
  if (abs(hitPos.local.z()) > (0.5*hitInfo.thickness() + m_numericLimit_float)) {
    warning() << "SimHit position is outside the sensor volume (local w = " << hitPos.local.z() << " mm, sensor thickness = " << hitInfo.thickness() << " mm). This should never happen. Forcing it to w=0." << endmsg;
    hitPos.local.z() = 0.;
  }
  
  // calculate how far the entry point is from the simHitPos, in terms of the simHitPath
  scaleFactor = 0.;
  if (simHitPath.z() >= 0.) {
    const float shiftDist_w = 0.5 * hitInfo.thickness() + hitPos.local.z();
    scaleFactor = shiftDist_w / simHitPath.z();
  }
  else {
    const float shiftDist_w = -0.5 * hitInfo.thickness() + hitPos.local.z();
    scaleFactor = shiftDist_w / simHitPath.z();
  }
  
  // Entry position is on either surface of the sensor, upstream from the simHit position (exactly by the fraction we just calculated)
  dd4hep::rec::Vector3D simHitEntryPos = hitPos.local - scaleFactor * simHitPath;

  // Clip the path to the sensor edges in u and v direction (if it passes through the side of the sensor)
  float t_min = 0.f, t_max = 1.f;
  std::tie(t_min, t_max) = computePathClippingFactors(t_min, t_max, simHitEntryPos.x(), simHitPath.x(), hitInfo.length(0));
  std::tie(t_min, t_max) = computePathClippingFactors(t_min, t_max, simHitEntryPos.y(), simHitPath.y(), hitInfo.length(1));

  // Apply clipping
  if (t_min != 0.f || t_max != 1.f) { // check if clipping is even necessary, for performance
    hitInfo.setDebugFlag();
    if (0. <= t_min && t_min <= t_max && t_max <= 1.) {
      verbose() << "   - Clipping simHitPath to sensor edges with factors t_min = " << t_min << ", t_max = " << t_max << ". PathLength changed to " << (t_max - t_min) * simHitPath.r() << " mm from " << simHitPath.r() << " mm" << endmsg;
      
      simHitEntryPos = simHitEntryPos + t_min * simHitPath;
      simHitPath = (t_max - t_min) * simHitPath;

      // if pathLength given by Geant4 is more than X% shorter than the length we calculate, shorten our calculated path accordingly.
      if (simHitPath.r() > m_pathLengthShorteningFactorGeant4 * hitInfo.simPathLength()) {
        verbose() << "   - Shortening simHitPath from " << simHitPath.r() << " mm to Geant4 pathLength of " << hitInfo.simPathLength() << " mm, because it's length is more than " << m_pathLengthShorteningFactorGeant4 << " of the Geant4 path length." << endmsg;
        
        /* make sure the path stays as centered around the simHitPos as possible (regarding edges)  
        * find out where the simHitPos lies on the current path:
        * project (simHitPos - simHitEntryPos) onto simHitPath -> gives distance from entryPos to simHitPos along the path (or to the point on path closest to simHitPos)
        * dotProduct (simHitPos-simHitEntryPos, simHitPath) = |simHitPos - simHitEntryPos| * |simHitPath| * cos(angle between them) */
        const float t_simHit = ( (hitPos.local - simHitEntryPos).dot(simHitPath) ) / ( simHitPath.r() * simHitPath.r() ); 
        
        // clamp new path (centered at t_center) to [0,1]
        const float t_length_half = hitInfo.simPathLength() / simHitPath.r() / 2.f; // half-length of new path, in terms of t [0,1] on old path
        const float t_center = std::max(t_length_half, std::min(t_simHit, 1.f - t_length_half)); // center of new path, clamped to [t_length_half, 1 - t_length_half] while not exceeding [0,1]

        t_min = t_center - t_length_half;
        t_max = t_center + t_length_half;

        simHitEntryPos = simHitEntryPos + t_min * simHitPath;
        simHitPath = (t_max - t_min) * simHitPath;
      }
    } 
    else {
      warning() << "constructSimHitPath(): Cannot clip simHitPath to sensor edges. The path lies completely outside the sensor volume. Clipping t_min = " << t_min << ", t_max = " << t_max << "." << endmsg;
      verbose() << "   - before clipping: EntryPos: (" << simHitEntryPos.x() << " mm, " << simHitEntryPos.y() << " mm, " << simHitEntryPos.z() << " mm), exitPos: (" << simHitPath.x() << " mm, " << simHitPath.y() << " mm, " << simHitPath.z() << " mm)" << endmsg;
      // hitInfo.setDebugFlag();
      // return std::make_tuple(dd4hep::rec::Vector3D(0.0, 0.0, 0.0), dd4hep::rec::Vector3D(0.0, 0.0, 0.0));
    }
  }

  verbose() << "   - Calculated SimHitPath with length " << simHitPath.r() << " mm" << " (the Geant4 path length is " << hitInfo.simPathLength() << " mm)" << endmsg;
  return std::make_tuple(simHitEntryPos, simHitPath);
}

VTXdigi_Allpix2::PixelChargeMatrix VTXdigi_Allpix2::DepositAndCollectCharge(HitInfo& hitInfo, const HitPosition& hitPos) const {
  const auto& [i_u_simHit, i_v_simHit] = computePixelIndices(hitPos.local, hitInfo.layerIndex(), hitInfo.length(0), hitInfo.length(1));
  
  PixelChargeMatrix pixelChargeMatrix(
    i_u_simHit,
    i_v_simHit,
    m_maxClusterSize.value().at(0),
    m_maxClusterSize.value().at(1)
  );

  const float segmentCharge = hitInfo.charge() / hitInfo.nSegments();
  
  /* get first segment */
  SegmentIndices segment = computeSegmentIndices(hitInfo, hitPos.entry, hitPos.path, 0);
  int segmentsInBin = 1;
  SegmentIndices nextSegment;

  /* loop over segments */
  for (int nextSegmentIndex = 1; nextSegmentIndex < hitInfo.nSegments(); ++nextSegmentIndex) {

    nextSegment = computeSegmentIndices(hitInfo, hitPos.entry, hitPos.path, nextSegmentIndex);

    if (segment == nextSegment) {
      verbose() << "       - Segment lies in the same pixel and in-pixel bin as previous segment, continuing." << endmsg;
      ++segmentsInBin;
      continue;
    } 
    else {
      verbose() << "       - Crossed bin-boundary wrt. last segment. Sharing " << segmentCharge*segmentsInBin << " e- from " << segmentsInBin << " segments. The last segment of these has nextSegmentIndex " << nextSegmentIndex-1 << "." << endmsg;
      collectSegmentCharge(hitInfo, pixelChargeMatrix, segment, segmentCharge); // write charge for this set of segments into pixelChargeMatrix (avoids copying the matrix in memory every time we write into it)

      segment = nextSegment;
      segmentsInBin = 1;
    }
  } // loop over segments

  /* write out last set of segments */
  verbose() << "       - Reached last segment. Sharing " << segmentCharge*segmentsInBin << " e- from last " << segmentsInBin << " segments." << endmsg;
  collectSegmentCharge(hitInfo, pixelChargeMatrix, segment, segmentCharge); // write charge for this segment into pixelChargeMatrix (done like this to avoid copying the matrix in memory every time we write into it)

  return pixelChargeMatrix;
}


void VTXdigi_Allpix2::collectSegmentCharge(HitInfo& hitInfo, PixelChargeMatrix& pixelChargeMatrix, const SegmentIndices& segment, const float segmentCharge) const {
  if (segment.i_u == -1) { // computeSegmentIndices() returns -1 if any dimension is outside sensor volume
      warning() << "Applying Kernel: Bin lies outside sensor volume. Dismissing." << endmsg;
      return; 
    }

  int i_u_target, i_v_target;

  /* each kernel entry shares charge from the source-pixel i_u, i_v to one target-pixel.
      * The kernel is centered on the source pixel, so loop over all target pixels covered by the kernel.
      * i_x_previous defines the source-pixel */
  for (int i_m = -1*(m_kernelSize-1)/2; i_m<=(m_kernelSize-1)/2; i_m++) {
    i_u_target = segment.i_u + i_m;
    if (i_u_target<0 || i_u_target>=hitInfo.pixCount(0))
      continue; // target pixel outside pixel matrix in u

    for (int i_n = -1*(m_kernelSize-1)/2; i_n<=(m_kernelSize-1)/2; i_n++) {
      i_v_target = segment.i_v + i_n;

      if (i_v_target<0 || i_v_target>=hitInfo.pixCount(1))
        continue; // target pixel outside pixel matrix in v

      const float kernelEntry = m_chargeSharingKernels->GetWeight(segment, i_m, i_n);
      if (kernelEntry < m_numericLimit_float) 
        continue; // skip zero entries
        
      const float sharedCharge = kernelEntry * segmentCharge;
      pixelChargeMatrix.FillRawCharge(i_u_target, i_v_target, sharedCharge);

      if (m_debugHistograms)
        fillDebugHistograms_segmentLoop(hitInfo, segment, i_m, i_n, sharedCharge);
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
      dd4hep::rec::Vector3D pixelCenterLocal = computePixelCenter_Local(i_u, i_v, hitInfo.layerIndex(), *hitInfo.simSurface());
      dd4hep::rec::Vector3D pixelCenterGlobal = transformLocalToGlobal(pixelCenterLocal, hitInfo.cellID());

      debug() << "       - Pixel (" << i_u << ", " << i_v << ") at (" << pixelCenterLocal.x() << ", " << pixelCenterLocal.y() << ", " << pixelCenterLocal[2] << ") mm received a measured/raw charge of " << pixelChargeMeasured << "/" << pixelChargeRaw << " e-, center at global position " << pixelCenterGlobal[0] << " mm, " << pixelCenterGlobal[1] << " mm, " << pixelCenterGlobal[2] << " mm" << endmsg;
      createDigiHit(simHit, digiHits, digiHitsLinks, pixelCenterGlobal, pixelChargeMeasured);
      ++m_counter_digiHitsCreated;

      if (m_debugHistograms)
        fillDebugHistograms_targetPixelLoop(hitInfo, hitPos, i_u, i_v, pixelChargeMeasured);
    }
  } // end loop over pixels, create digiHits

  if (m_debugHistograms) {
    ++(*m_histogramsGlobal.at(histGlobal_clusterSize_raw))[static_cast<double>(nPixelsReceivedCharge)]; // cluster size = number of pixels fired per simHit
    ++(*m_histogramsGlobal.at(histGlobal_clusterSize_measured))[static_cast<double>(nPixelsFired)];
    
    ++(*m_histograms1d.at(hitInfo.layerIndex()).at(hist1d_ClusterSize_raw))[static_cast<double>(nPixelsReceivedCharge)];
    ++(*m_histograms1d.at(hitInfo.layerIndex()).at(hist1d_ClusterSize_measured))[static_cast<double>(nPixelsFired)];
  }
    

  if (nPixelsFired <= 0) {
    ++m_counter_acceptedButNoSegmentsInSensor;
    debug() << "   - No pixels fired for this simHit" << endmsg;
  }
  else {
    verbose() << "   - From this simHit, " << nPixelsFired << " pixels received charge." << endmsg;
  }
} // AnalyseSharedCharge()

void VTXdigi_Allpix2::FillGeneralDebugHistograms(HitInfo& hitInfo, const HitPosition& hitPos, const PixelChargeMatrix& pixelChargeMatrix) const {
  /* mostly implements per-simhit histograms */
  verbose() << "   - Filling 1D histograms" << endmsg;
  ++(*m_histogramsGlobal.at(histGlobal_simHitCharge))[hitInfo.charge()]; // in e
  ++(*m_histogramsGlobal.at(histGlobal_simHitE))[hitInfo.charge() / m_chargePerkeV]; // in keV
  ++(*m_histogramsGlobal.at(histGlobal_EntryPointX))[hitPos.entry.x()]; // in mm
  ++(*m_histogramsGlobal.at(histGlobal_EntryPointY))[hitPos.entry.y()]; // in mm
  ++(*m_histogramsGlobal.at(histGlobal_EntryPointZ))[hitPos.entry.z()*1000]; // convert mm to um
  ++(*m_histogramsGlobal.at(histGlobal_pathLength))[hitPos.path.r()*1000]; // convert mm to um
  ++(*m_histogramsGlobal.at(histGlobal_pathLengthGeant4))[hitInfo.simPathLength()*1000]; // convert mm to um
  ++(*m_histogramsGlobal.at(histGlobal_chargeCollectionEfficiency_raw))[ pixelChargeMatrix.GetTotalRawCharge() / hitInfo.charge() ]; // in e-, only if at least one pixel fired
  ++(*m_histogramsGlobal.at(histGlobal_chargeCollectionEfficiency))[ pixelChargeMatrix.GetTotalMeasuredCharge(m_pixelThreshold) / hitInfo.charge()];
  
  ++(*m_histograms1d.at(hitInfo.layerIndex()).at(hist1d_SimHitMomentum))[hitInfo.simMomentum()*1000]; // Convert GeV to MeV

  int pix_u, pix_v;
  std::tie(pix_u, pix_v) = computePixelIndices(hitPos.local, hitInfo.layerIndex(), hitInfo.length(0), hitInfo.length(1));
  verbose() << "   - Filling 2D histograms for layer " << m_cellIDdecoder->get(hitInfo.cellID(), "layer") << " (index " << hitInfo.layerIndex() << ")" << endmsg;
  ++(*m_histograms2d.at(hitInfo.layerIndex()).at(hist2d_hitMap_simHits))[{pix_u, pix_v}]; // in mm
  ++(*m_histograms2d.at(hitInfo.layerIndex()).at(hist2d_pathLength_vs_simHit_v))[{pix_v, hitPos.path.r()*1000}]; // in um and mm

  ++(*m_histogramsGlobal.at(histGlobal_pixelChargeMatrix_size_u))[pixelChargeMatrix.GetSize_u()];
  ++(*m_histogramsGlobal.at(histGlobal_pixelChargeMatrix_size_v))[pixelChargeMatrix.GetSize_v()];
  
  ++(*m_histograms2d.at(hitInfo.layerIndex()).at(hist2d_pixelChargeMatrixSize))[{pixelChargeMatrix.GetSize_u(), pixelChargeMatrix.GetSize_v()}];

  float pathAnglePhi = atan2(hitPos.path.x(), hitPos.path.z()) / 3.14159265 * 180.f; // in degrees
  float pathAngleTheta = atan2(hitPos.path.y(), hitPos.path.z()) / 3.14159265 * 180.f; // in degrees
  ++(*m_histograms1d.at(hitInfo.layerIndex()).at(hist1d_IncidentAngle_ThetaLocal))[pathAngleTheta];
  ++(*m_histograms1d.at(hitInfo.layerIndex()).at(hist1d_IncidentAngle_PhiLocal))[pathAnglePhi];
  ++(*m_histograms2d.at(hitInfo.layerIndex()).at(hist2d_IncidentAngle))[{pathAngleTheta, pathAnglePhi}];
}

void VTXdigi_Allpix2::fillDebugHistograms_segmentLoop(const HitInfo& hitInfo, const SegmentIndices& segment, int i_m, int i_n, const float sharedCharge) const {

  /* work out the distance between target pixel center and the origin bin, in terms of pixels */
  const float dist_u = -i_m - 0.5f + (segment.j_u + 0.5f) / static_cast<float>(m_inPixelBinCount[0]); // in pixels
  const float dist_v = -i_n - 0.5f + (segment.j_v + 0.5f) / static_cast<float>(m_inPixelBinCount[1]); // in pixels
  const float pos_w = ( (segment.j_w + 0.5f) * (hitInfo.thickness() / m_inPixelBinCount[2]) - hitInfo.thickness()/2.f) * 1000.f; // conversion from mm to um
  
  (*m_histWeighted2d.at(hitInfo.layerIndex()).at(histWeighted2d_chargeOriginU))[{dist_u, pos_w}] += sharedCharge; // in e-
  (*m_histWeighted2d.at(hitInfo.layerIndex()).at(histWeighted2d_chargeOriginV))[{dist_v, pos_w}] += sharedCharge; // in e-
}

void VTXdigi_Allpix2::fillDebugHistograms_targetPixelLoop(const HitInfo& hitInfo, const HitPosition& hitPos, int i_u, int i_v, float pixelChargeMeasured) const {
  const auto& [i_u_simHit, i_v_simHit] = computePixelIndices(hitPos.local, hitInfo.layerIndex(), hitInfo.length(0), hitInfo.length(1));

  (*m_histWeighted2d.at(hitInfo.layerIndex()).at(histWeighted2d_averageCluster))[{i_u - i_u_simHit, i_v - i_v_simHit}] += pixelChargeMeasured; // in e-

  dd4hep::rec::Vector3D pixelCenterLocal = computePixelCenter_Local(i_u, i_v, hitInfo.layerIndex(), *hitInfo.simSurface());
  dd4hep::rec::Vector3D pixelCenterGlobal = transformLocalToGlobal(pixelCenterLocal, hitInfo.cellID());
  
  ++ (*m_histogramsGlobal.at(histGlobal_DisplacementU))[ (pixelCenterLocal.x() - hitPos.local.x()) * 1000]; // convert mm to um
  ++ (*m_histogramsGlobal.at(histGlobal_DisplacementV))[ (pixelCenterLocal.y() - hitPos.local.y()) * 1000]; // convert mm to um
  ++ (*m_histogramsGlobal.at(histGlobal_DisplacementR))[sqrt( (pixelCenterGlobal.x() - hitPos.global.x())*(pixelCenterGlobal.x() - hitPos.global.x()) + (pixelCenterGlobal.y() - hitPos.global.y())*(pixelCenterGlobal.y() - hitPos.global.y()) + (pixelCenterGlobal.z() - hitPos.global.z())*(pixelCenterGlobal.z() - hitPos.global.z()) ) * 1000]; // convert mm to um

}


/* -- Pixel and in-pixel binning magic (logic) -- */

int VTXdigi_Allpix2::computeBinIndex(float x, float binX0, float binWidth, int binN) const {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  if (binN <= 0) throw GaudiException("computeBinIndex: binN must be positive", "VTXdigi_Allpix2::computeBinIndex()", StatusCode::FAILURE);
  if (binWidth <= 0.0) throw GaudiException("computeBinIndex: binWidth must be positive", "VTXdigi_Allpix2::computeBinIndex()", StatusCode::FAILURE);

  float relativePos = (x - binX0) / binWidth; // shift to [0, binN]
  if (relativePos < 0.0f || relativePos >= static_cast<float>(binN))
    return -1;
  return static_cast<int>(relativePos);
} // computeBinIndex()

std::tuple<int, int> VTXdigi_Allpix2::computePixelIndices(const dd4hep::rec::Vector3D& segmentPos, const int layerIndex, const float length_u, const float length_v) const {
  
  int i_u = computeBinIndex(
    segmentPos.x(),
    -0.5*length_u,
    m_pixelPitch_u.at(layerIndex),
    m_pixelCount_u.value().at(layerIndex));

  int i_v = computeBinIndex(
    segmentPos.y(),
    -0.5*length_v,
    m_pixelPitch_v.at(layerIndex),
    m_pixelCount_v.value().at(layerIndex));
  return std::make_tuple(i_u, i_v);
} // computePixelIndices()

std::tuple<int, int, int> VTXdigi_Allpix2::computeInPixelIndices(const dd4hep::rec::Vector3D& segmentPos, const int layerIndex, const float length_u, const float length_v) const {
  int j_u, j_v, j_w;

  verbose() << "         - Computing in-pixel indices. Number of in-pix bins: (" << m_inPixelBinCount[0] << ", " << m_inPixelBinCount[1] << ", " << m_inPixelBinCount[2] << ")" << endmsg;

  float shiftedPos_u = segmentPos.x() + 0.5 * length_u; // shift to [0, length_u]
  float pitch_u = m_pixelPitch_u.at(layerIndex);
  float inPixelPos_u = std::fmod(shiftedPos_u, pitch_u);
  if (inPixelPos_u < 0.0) inPixelPos_u += pitch_u; // ensure positive remainder

  j_u = computeBinIndex(inPixelPos_u, 0.0, pitch_u / m_inPixelBinCount[0], m_inPixelBinCount[0]);

  float shiftedPos_v = segmentPos.y() + 0.5 * length_v;
  float pitch_v = m_pixelPitch_v.at(layerIndex);
  float inPixelPos_v = std::fmod(shiftedPos_v, pitch_v);
  if (inPixelPos_v < 0.0) inPixelPos_v += pitch_v;

  j_v = computeBinIndex(inPixelPos_v, 0.0, pitch_v / m_inPixelBinCount[1], m_inPixelBinCount[1]);

  // vertical (w) binning: shift to [0, thickness]
  float shiftedPos_w = segmentPos.z() + 0.5 * m_sensorThickness.at(layerIndex);
  j_w = computeBinIndex(shiftedPos_w, 0.0, m_sensorThickness.at(layerIndex) / m_inPixelBinCount[2], m_inPixelBinCount[2]);

  return std::make_tuple(j_u, j_v, j_w);
} // computeInPixelIndices()

VTXdigi_Allpix2::SegmentIndices VTXdigi_Allpix2::computeSegmentIndices(HitInfo& hitInfo, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const int segmentIndex) const {

  verbose() << "     - Processing segment " << segmentIndex << " out of " << hitInfo.nSegments() << " with length " << simHitPath.r() << " mm" << endmsg;
  if (segmentIndex < 0 || segmentIndex >= hitInfo.nSegments()) {
    error() << "computeSegmentIndices(): Invalid segment number " << hitInfo.nSegments() << " or segment index " << segmentIndex << endmsg;
    throw std::runtime_error("VTXdigi_Allpix2::computeSegmentIndices(): Invalid segment number or segment index");
  }
  
  SegmentIndices segment;
  float pathFraction = (static_cast<float>(segmentIndex)+0.5) / hitInfo.nSegments();
  const dd4hep::rec::Vector3D segmentRelativePos = pathFraction * simHitPath;
  const dd4hep::rec::Vector3D segmentPos = simHitEntryPos + segmentRelativePos;
  verbose() << "       - Local position (u,v,w): (" << segmentPos.x() << " mm, " << segmentPos.y() << " mm, " << segmentPos.z() << " mm)" << endmsg;

  /* Compute pixel indices */
  std::tie(segment.i_u, segment.i_v) = computePixelIndices(segmentPos, hitInfo.layerIndex(), hitInfo.length(0), hitInfo.length(1));
  if (segment.i_u==-1 || segment.i_v==-1) {
    warning() << "computeSegmentIndices(): Segment lies outside sensor area (in u or v). Dismissing." << endmsg;
    // hitInfo.setDebugFlag();
    SegmentIndices emptySegment;
    return emptySegment;
  }
  
  /* Compute in-pixel indices */
  std::tie(segment.j_u, segment.j_v, segment.j_w) = computeInPixelIndices(segmentPos, hitInfo.layerIndex(), hitInfo.length(0), hitInfo.length(1));
  if (segment.j_u==-1 || segment.j_v==-1 || segment.j_w==-1) {
    warning() << "computeSegmentIndices(): Segment lies inside sensor area (in u and v), but vertically outside sensor volume. Dismissing." << endmsg;
    SegmentIndices emptySegment;
    return emptySegment;
  }
  verbose() << "       - Pixel indices (" << segment.i_u << ", " << segment.i_v << "), In-pixel indices (" << segment.j_u << ", " << segment.j_v << ", " << segment.j_w << ")" << endmsg;

  return segment;
}

dd4hep::rec::Vector3D VTXdigi_Allpix2::computePixelCenter_Local(const int i_u, const int i_v, const int layerIndex, const dd4hep::rec::ISurface& simSurface) const {
  // returns the position of the center of pixel i_u, i_v) in the global
  // u - short ARCADIA axis, corresponds to x in sensor local frame
  // v - long ARCADIA axis, corresponds to y in sensor local frame

  float length_u = simSurface.length_along_u() * 10; // convert to mm (works, checked 2025-10-17)
  float length_v = simSurface.length_along_v() * 10; // convert to mm 

  float posU = -0.5 * length_u + (i_u + 0.5) * m_pixelPitch_u.at(layerIndex); // in mm
  float posV = -0.5 * length_v + (i_v + 0.5) * m_pixelPitch_v.at(layerIndex); // in mm

  return dd4hep::rec::Vector3D(posU, posV, 0.); 
}


/* -- Transformation between global frame and local sensor frame -- */

TGeoHMatrix VTXdigi_Allpix2::computeTransformationMatrix(const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = m_volumeManager.lookupDetElement(cellID).nominal().worldTransformation(); // given in cm

  // rotation is unitless, but need to convert translation from cm to mm
  double* translationComponent = transformationMatrix.GetTranslation();
  translationComponent[0] = translationComponent[0] * 10; // convert to mm
  translationComponent[1] = translationComponent[1] * 10;
  translationComponent[2] = translationComponent[2] * 10;
  transformationMatrix.SetTranslation(translationComponent);

  setProperDirectFrame(transformationMatrix); // Change coordinates to have z orthogonal to sensor, with direct (right-handed) frame

  return transformationMatrix;
}
TGeoHMatrix VTXdigi_Allpix2::computeTransformationMatrix(const edm4hep::SimTrackerHit& simHit) const {
  return computeTransformationMatrix(simHit.getCellID());
}

void VTXdigi_Allpix2::setProperDirectFrame(TGeoHMatrix& transformationMatrix) const {
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

dd4hep::rec::Vector3D VTXdigi_Allpix2::transformGlobalToLocal(const dd4hep::rec::Vector3D& globalPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = computeTransformationMatrix(cellID);

  double localPos[3] = {0, 0, 0};

  transformationMatrix.MasterToLocal(globalPos, localPos);

  return dd4hep::rec::Vector3D(localPos[0], localPos[1], localPos[2]);
}

dd4hep::rec::Vector3D VTXdigi_Allpix2::transformLocalToGlobal(const dd4hep::rec::Vector3D& localPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = computeTransformationMatrix(cellID);

  double globalPos[3] = {0, 0, 0};

  transformationMatrix.LocalToMaster(localPos, globalPos);

  return dd4hep::rec::Vector3D(globalPos[0], globalPos[1], globalPos[2]);
}


/* -- Other Helper functions -- */

dd4hep::rec::Vector3D VTXdigi_Allpix2::convertVector(edm4hep::Vector3d vec) const {
  // return dd4hep::rec::Vector3D(vec.x*dd4hep::mm, vec.y*dd4hep::mm, vec.z*dd4hep::mm);
  return dd4hep::rec::Vector3D(vec.x, vec.y, vec.z);
}
edm4hep::Vector3d VTXdigi_Allpix2::convertVector(dd4hep::rec::Vector3D vec) const {
  // return edm4hep::Vector3d(vec[0]/dd4hep::mm, vec[1]/dd4hep::mm, vec[2]/dd4hep::mm);
  return edm4hep::Vector3d(vec.x(), vec.y(), vec.z());
}

std::tuple<float, float> VTXdigi_Allpix2::computePathClippingFactors(float t_min, float t_max, const float entryPos_ax, const float pathLength_ax, const float sensorLength_ax) const {
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

void VTXdigi_Allpix2::createDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const dd4hep::rec::Vector3D& position, const float charge) const {
  // overload to allow passing dd4hep::rec::Vector3D as position. ~ Jona 2025-09

  auto digiHit = digiHits.create();
  digiHit.setCellID(simHit.getCellID());
  verbose() << "         - Creating digiHit in cellID " << simHit.getCellID() << " (" << std::bitset<24>(simHit.getCellID()) << "), setting eDep to " << charge/m_chargePerkeV << " keV" << endmsg;
  digiHit.setEDep(charge / m_chargePerkeV); // convert e- to keV
  digiHit.setPosition(convertVector(position));
  // TODO: check if position is within sensor bounds & force it onto sensor simSurface ~ Jona 2025-09
  digiHit.setTime(simHit.getTime());
  
  auto digiHitLink = digiHitsLinks.create();
  digiHitLink.setFrom(digiHit);
  digiHitLink.setTo(simHit);
}

void VTXdigi_Allpix2::appendSimHitToCsv(const HitInfo& hitInfo, const HitPosition& hitPos, const int i_u, const int i_v) const {
  if (!m_debugCsvFile.is_open()) {
    error() << "appendSimHitToCsv(): DebugCsv file is not open." << endmsg;
    return;
  }

  // eventNumber,layerIndex,segmentCount,sensorThickness,pix_u,pix_v,
  // simHitPos_u,simHitPos_v,simHitPos_w,
  // simHitEntryPos_u,simHitEntryPos_v,simHitEntryPos_w,
  // simHitPath_u,simHitPath_v,simHitPath_w,
  // pathLengthGeant4,pathLength,chargeDeposition,debugFlag

  m_debugCsvFile << std::to_string(hitInfo.eventNumber()) << "," << hitInfo.layerIndex() << "," << hitInfo.nSegments() << ","  << hitInfo.thickness() << "," << i_u << "," << i_v << ",";

  m_debugCsvFile << hitPos.local.x() << "," << hitPos.local.y() << "," << hitPos.local.z() << ",";
  m_debugCsvFile << hitPos.entry.x() << "," << hitPos.entry.y() << "," << hitPos.entry.z() << ",";
  m_debugCsvFile << hitPos.path.x() << "," << hitPos.path.y() << "," << hitPos.path.z() << ",";
  m_debugCsvFile << hitInfo.simPathLength() << "," << hitPos.path.r() << "," << hitInfo.charge() <<  "," << hitInfo.debugFlag() << "\n";
  m_debugCsvFile.flush();
  verbose() << "Wrote simHit with event " << hitInfo.eventNumber() << ", layerIndex " << hitInfo.layerIndex() << ", segmentN " << hitInfo.nSegments() << ", i_u " << i_u << ", i_v " << i_v << " to debug CSV file." << endmsg;
  return;
}



