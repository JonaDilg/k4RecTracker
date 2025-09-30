#include "VTXdigi_Allpix2.h"

DECLARE_COMPONENT(VTXdigi_Allpix2)


/* Notes ~ Jona 2025-09
 * - all lengths in mm (edm4hep already does this)
 *   BUT: dd4hep uses cm internally, so convert when interfacing with dd4hep via dd4hep::mm = 0.1
 *        mm -> cm: a [cm] = dd4hep::mm * a [mm]
 *        cm -> mm: a [mm] = 1/dd4hep::mm * a [cm]
 * - energies in keV
 * - Vectors can be given in 
 *        a) dd4hep::rec::Vector3D <- fully featured vector, overloads operators *+- etc
 *        b) edm4hep::Vector3d <- natively used by edm4hep (where simHit, digiHit are from)
 *      -> generally use dd4hep::rec::Vector3D, overload functions where edm4hep::Vector3d is needed
 */

VTXdigi_Allpix2::VTXdigi_Allpix2(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("SimTrackHitCollectionName", {"SimTrackerHits"}),
                        KeyValues("HeaderName", {"EventHeader"}),},
                       {KeyValues("TrackerHitCollectionName", {"VTXTrackerHits"}),
                        KeyValues("SimTrkHitRelationsCollection", {"VTXTrackerHitRelations"})}) {

  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }  

  // -- Debugging Histograms --

  m_histograms[hist_hitE].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "hitE", "SimHit Energy;keV", {1000, 0, 200}});
  m_histograms[hist_clusterSize].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "ClusterSize", "Cluster Size (ie. number of digiHits per accepted simHit)", {20, -0.5, 19.5}});

  m_histograms[hist_EntryPointX].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "EntryPointX", "SimHit Entry Point U (in local sensor frame);mm", {400, -4, 4}});
  m_histograms[hist_EntryPointY].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "EntryPointY", "SimHit Entry Point V (in local sensor frame);mm", {2000, -20, 20}});
  m_histograms[hist_EntryPointZ].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "EntryPointZ", "SimHit Entry Point W (in local sensor frame);um", {1000, -300, 300}});

  m_histograms[hist_DisplacementU].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "LocalDisplacementU", "Displacement U (local sensor frame): digiHit_u - simHit_u;um", {200, -200, 200}});
  m_histograms[hist_DisplacementV].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "LocalDisplacementV", "Displacement V (local sensor frame): digiHit_v - simHit_v;um", {200, -200, 200}});
  m_histograms[hist_DisplacementR].reset(
    new Gaudi::Accumulators::StaticRootHistogram<1>{this, "GlobalDisplacementR", "Displacement R (global frame): | digiHit - simHit |;um", {300, 0., 300}});

  info() << "Constructed successfully" << endmsg;
}

StatusCode VTXdigi_Allpix2::initialize() {
  debug() << "Initializing..." << endmsg;
  m_geometryService = serviceLocator()->service(m_geometryServiceName);
  if (!m_geometryService) {
    error() << " - Unable to retrieve the GeoSvc. Abort." << endmsg;
    return StatusCode::FAILURE;
  }

  // retrieve the volume manager
  m_volumeManager = m_geometryService->getDetector()->volumeManager();

  // retrieve the cellID encoding string from GeoSvc
  std::string cellIDstr = m_geometryService->constantAsString(m_encodingStringVariable.value());
  m_cellIDdecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDstr);

  // get surface map
  const auto detector = m_geometryService->getDetector();
  const auto surfaceManager = detector->extension<dd4hep::rec::SurfaceManager>();
  // dd4hep::DetElement subDetector = detector->detector(m_subDetName.value()); // not used? ~ Jona 2025-09
  m_surfaceMap = surfaceManager->map(m_subDetName.value());

  /* pixel pitch and count are given as either
  *  - a single value, which is then applied to all layers
  *    -> rewrite the vector to contain the same value for all layers (simplifies code later)
  *  - a vector of values, one per layer
  *    -> check that the number of entries matches the number of layers in the geometry
  */ 
  if (m_pixelPitch_u.value().size() == 1 && m_pixelPitch_v.value().size() == 1 && m_pixelCount_u.value().size() == 1 && m_pixelCount_v.value().size() == 1) {
    // single value for all layers, rewrite to vector of correct size
    m_pixelPitch_u.value().resize(m_layerCount, m_pixelPitch_u.value()[0]);
    m_pixelPitch_v.value().resize(m_layerCount, m_pixelPitch_v.value()[0]);
    m_pixelCount_u.value().resize(m_layerCount, m_pixelCount_u.value()[0]);
    m_pixelCount_v.value().resize(m_layerCount, m_pixelCount_v.value()[0]);
    debug() << " - found single pixel pitch and pixel count values, applying them to all layers" << endmsg;
  }
  else if (m_pixelPitch_u.value().size() == m_layerCount && m_pixelPitch_v.value().size() == m_layerCount && m_pixelCount_u.value().size() == m_layerCount && m_pixelCount_v.value().size() == m_layerCount) {
      // one entry per layer
      debug() << " - found pixel pitch and pixel count values for all " << m_layerCount << " layers" << endmsg;
    }
  else {
    error() << " - Pixel pitch and count must be given either as a single value for all layers (in brackets though: []) or as a vector with one entry per layer. Abort." << endmsg;
    return StatusCode::FAILURE;
  }


  // -- import charge sharing kernels --

  debug() << " - Importing charge sharing kernels..." << endmsg;
  m_chargeSharingKernels = std::make_unique<ChargeSharingKernels>(m_InPixelBinCount_u.value(), m_InPixelBinCount_v.value(), m_InPixelBinCount_w.value(), m_KernelSize.value());

  // TODO: read the kernels from file. Instead, set gaussian kernels for now

  const std::vector<double> inputKernel_Gaussian = { // 3x3 kernel, gaussian
    1./16., 2./16., 1./16.,
    2./16., 4./16., 2./16.,
    1./16., 2./16., 1./16.
  };

  const std::vector<double> inputKernel_Basic = { // 3x3 kernel, double gaussian
    0,0,0,
    0,1,0,
    0,0,0
  };

  for (int i_u=0; i_u<m_InPixelBinCount_u.value(); i_u++) {
    for (int i_v=0; i_v<m_InPixelBinCount_v.value(); i_v++) {
      for (int i_w=0; i_w<m_InPixelBinCount_w.value(); i_w++) {
        m_chargeSharingKernels->SetKernel(i_u, i_v, i_w, inputKernel_Basic);
      }
    }
  }

  // TODO: check that pixel pitch and count match sensor size in geometry (I am not sure how to get the sensor size from the geometry though, does seem to be a bit more general, subDetectors have children that might or might not be layers) ~ Jona 2025-09
  // TODO load lookup table for charge transport and diffusion
  // TODO check that k4geo has same pitch as lookup table 
  info() << " - Initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}


// -- event loop -- 

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> VTXdigi_Allpix2::operator()
  (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {

  debug() << "PROCESSING event. run " << headers[0].getRunNumber() << "; event " << headers[0].getEventNumber() << "; with " << simHits.size() << " simHits" << endmsg;

  if (simHits.size()==0) {
    info() << " - No SimTrackerHits in collection, returning empty output collections" << endmsg;
    // TODO: Bookkeeping: count empty events, fill histogram etc.
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());
  }
  
  // create output collections
  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitsLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();
  unsigned int nHitsRead=0, nHitsCreated=0, nHitsDismissed=0;

  // loop over sim hits, digitize them, create output collections
  for (const auto& simHit : simHits) {

    // -- gather information about the simHit --    
    dd4hep::DDSegmentation::CellID cellID = simHit.getCellID(); // TODO: apply 64-bit length mask as done in DDPlanarDigi.cpp? ~ Jona 2025-09
    dd4hep::rec::Vector3D simHitGlobalPos = Vector3dConvert(simHit.getPosition());
    dd4hep::rec::Vector3D simHitLocalPos = GlobalToLocal(simHitGlobalPos, cellID);
    int layer = m_cellIDdecoder->get(cellID, "layer");
    double eDeposition = simHit.getEDep() * (dd4hep::GeV / dd4hep::keV); // convert to keV
    
    const dd4hep::rec::ISurface* surface = m_surfaceMap->find(cellID)->second;

    debug() << " - FOUND simHit in cellID " << simHit.getCellID() << ", Edep = " << eDeposition << " keV, path length = " << simHit.getPathLength()*1000 << " um" << endmsg;
    debug() << "   at position " << simHit.getPosition().x << " mm, " << simHit.getPosition().y << " mm, " << simHit.getPosition().z << " mm" << endmsg;
    debug() << "   with distance " << surface->distance(dd4hep::mm * simHitGlobalPos) << " mm to surface" << endmsg; // dd4hep expects cm, simHitGlobalPosition is in mm. dd4hep::cm=0.1
    
    // -- apply cuts --
    if (!ApplySimHitCuts(layer, eDeposition)) {
      
      nHitsDismissed++;
      continue;
    }

    // initial sensor size check
    // do NOT do this before cuts, as otherwise we might check layers that are not used
    if (m_initialSensorSizeCheckPassed.find(layer) == m_initialSensorSizeCheckPassed.end())
      InitialSensorSizeCheck(simHit);

    nHitsRead++;
    debug() << "   - ACCEPTED SimHit in layer " << layer << ". Edep =" << eDeposition << " keV" << endmsg;
    

    // -- get simHit path, segment into ionization points. (approximates paths as straight line) --

    // 1) split the path into equally long segments
    float simHitPathLength = simHit.getPathLength(); // in mm
    int segmentNumber = int( simHitPathLength / m_targetPathSegmentLength ); // both in mm
    segmentNumber = std::max(1, segmentNumber); // want at least one segment
    float segmentLength = simHitPathLength / segmentNumber;
    float segmentEdep = eDeposition / segmentNumber;
    debug() << "   - Segmenting path of length " << simHitPathLength*1000 << " um into " << segmentNumber << " segments of length " << segmentLength*1000 << " um" << endmsg;

    // 2) get path direction (unit vector), entry point (local frame)
    auto [simHitPath, simHitEntryPoint] = GetSimHitPath(simHit); // both are of type dd4hep::rec::Vector3D, given in mm
    debug() << "   - Entry point (sensor local) : " << simHitEntryPoint[0] << " mm, " << simHitEntryPoint[1] << " mm, " << simHitEntryPoint[2] << " mm" << endmsg;

    // get pixel edges in local frame
    /* u - short ARCADIA axis, corresponds to x in sensor local frame
     * v - long ARCADIA axis, corresponds to y in sensor local frame
     * w - orthogonal to sensor plane, corresponds to z in sensor local frame
    */
    const double length_u = surface->length_along_u() * 10; // convert to mm
    const double length_v = surface->length_along_v() * 10; // convert to mm
    const int size_u = m_pixelCount_u[layer];
    const int size_v = m_pixelCount_v[layer];

    // -- loop over the segments, assign charges to pixels according to kernel--

    debug() << "   - Looping over segments, assigning charges to pixels" << endmsg;
    std::vector<double> pixelCharge(size_u * size_v, 0.0); // array to store the charge assigned to each pixel, reset for each simHit. Index = i_u + i_v*size_u (Use 1-d vector to avoid overhead)
    int i_u, i_v; // pixel indices
    dd4hep::rec::Vector3D segmentPos; // in mm
    for (int i=0; i<segmentNumber; i++) {
      segmentPos = ((static_cast<double>(i)/segmentNumber)*simHitPath) + simHitEntryPoint; // in mm
      debug() << "     - Processing segment " << i << endmsg;
      
      // get pixel indices (i_u and i_v) for this segment
      i_u = GetBinIndex(segmentPos[0], -0.5*length_u, m_pixelPitch_u[layer], m_pixelCount_u[layer]);
      i_v = GetBinIndex(segmentPos[1], -0.5*length_v, m_pixelPitch_v[layer], m_pixelCount_v[layer]);

      // get the in-pixel position (for charge sharing kernel)
      

      // debug() << "       - Pixel indices (" << i_u << ", " << i_v << ")" << endmsg;
      if (i_u==-1 || i_v==-1) {
        debug() << "       - Dismissed segment (is outside sensor area)" << endmsg;
        continue; // segment is outside sensor area
      }
      
      // assign charge to pixel
      pixelCharge[i_u + i_v*size_u] += segmentEdep; // convert energy to number of electrons (3.61 eV per electron-hole pair in Si)
    }
    
    
    //  -- find fired pixels, create digiHit --
    int nPixelsFired = 0;
    
    debug() << "   - Processed all segments. Creating digiHits." << endmsg;
    for (i_u = 0; i_u<m_pixelCount_u[layer]; i_u++) {
      for (i_v = 0; i_v<m_pixelCount_v[layer]; i_v++) {
        if (pixelCharge[i_u + i_v*size_u]>0) {
          nPixelsFired++;
          
          dd4hep::rec::Vector3D pixelCenter = GetPixelCenter_Local(i_u, i_v, layer, *surface);
          pixelCenter = LocalToGlobal(pixelCenter, cellID); // convert to global frame
          
          debug() << "     - Pixel (" << i_u << ", " << i_v << ") received " << pixelCharge[i_u + i_v*size_u] << " keV, center at global position " << pixelCenter[0] << " mm, " << pixelCenter[1] << " mm, " << pixelCenter[2] << " mm" << endmsg;
          CreateDigiHit(simHit, digiHits, digiHitsLinks, pixelCharge[i_u + i_v*size_u], pixelCenter);
          nHitsCreated++;

          ++ (*m_histograms[hist_DisplacementU])[1000 * (pixelCenter[0] - simHitGlobalPos[0])]; // in um
          ++ (*m_histograms[hist_DisplacementV])[1000 * (pixelCenter[1] - simHitGlobalPos[1])]; // in um
          ++ (*m_histograms[hist_DisplacementR])[1000 * sqrt( (pixelCenter[0] - simHitGlobalPos[0])*(pixelCenter[0] - simHitGlobalPos[0]) + (pixelCenter[1] - simHitGlobalPos[1])*(pixelCenter[1] - simHitGlobalPos[1]) + (pixelCenter[2] - simHitGlobalPos[2])*(pixelCenter[2] - simHitGlobalPos[2]) )]; // in um
        }
      }
    }
    debug() << "   - From this simHit, " << nPixelsFired << " pixels received charge." << endmsg;


    // -- fill debugging histograms --

    ++(*m_histograms[hist_hitE])[eDeposition]; // in keV
    ++(*m_histograms[hist_EntryPointX])[simHitEntryPoint[0]]; // in mm
    ++(*m_histograms[hist_EntryPointY])[simHitEntryPoint[1]]; // in mm
    ++(*m_histograms[hist_EntryPointZ])[simHitEntryPoint[2]*1000]; // in um
    ++(*m_histograms[hist_clusterSize])[static_cast<double>(nPixelsFired)]; // cluster size = number of pixels fired per simHit
    
    nHitsRead++;
  } // end loop over sim hits

  debug() << "FINISHED event. Processed " << nHitsRead << " simHits. From this, created " << nHitsCreated << " digiHits and dismissed " << nHitsDismissed << " simHits." << endmsg;

  return std::make_tuple(std::move(digiHits), std::move(digiHitsLinks));
} // operator()





// -- Transformation between global frame and local sensor frame --

TGeoHMatrix VTXdigi_Allpix2::GetTransformationMatrix(const dd4hep::DDSegmentation::CellID& cellID) const {

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
TGeoHMatrix VTXdigi_Allpix2::GetTransformationMatrix(const edm4hep::SimTrackerHit& simHit) const {
  return GetTransformationMatrix(simHit.getCellID());
}

void VTXdigi_Allpix2::SetProperDirectFrame(TGeoHMatrix& transformationMatrix) const {
  /** Change the transformationMatrix to have a direct frame with z orthogonal to sensor surface
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
} // SetProperDirectFrame()

dd4hep::rec::Vector3D VTXdigi_Allpix2::GlobalToLocal(const dd4hep::rec::Vector3D& globalPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = GetTransformationMatrix(cellID);

  double localPos[3] = {0, 0, 0};

  // transformationMatrix.MasterToLocal(dd4hep::mm*globalPos, localPos); // globalPos in mm, but matrix in cm

  // return dd4hep::rec::Vector3D(localPos[0]/dd4hep::mm, localPos[1]/dd4hep::mm, localPos[2]/dd4hep::mm);

  transformationMatrix.MasterToLocal(globalPos, localPos);

  return dd4hep::rec::Vector3D(localPos[0], localPos[1], localPos[2]);
}

dd4hep::rec::Vector3D VTXdigi_Allpix2::LocalToGlobal(const dd4hep::rec::Vector3D& localPos, const dd4hep::DDSegmentation::CellID& cellID) const {

  TGeoHMatrix transformationMatrix = GetTransformationMatrix(cellID);

  double globalPos[3] = {0, 0, 0};

  // transformationMatrix.LocalToMaster(dd4hep::mm*localPos, globalPos); // both in mm

  // return dd4hep::rec::Vector3D(globalPos[0]/dd4hep::mm, globalPos[1]/dd4hep::mm, globalPos[2]/dd4hep::mm);

  transformationMatrix.LocalToMaster(localPos, globalPos);

  return dd4hep::rec::Vector3D(globalPos[0], globalPos[1], globalPos[2]);
}



// -- Helper functions --

dd4hep::rec::Vector3D VTXdigi_Allpix2::Vector3dConvert(edm4hep::Vector3d vec) const {
  // return dd4hep::rec::Vector3D(vec.x*dd4hep::mm, vec.y*dd4hep::mm, vec.z*dd4hep::mm);
  return dd4hep::rec::Vector3D(vec.x, vec.y, vec.z);
}
edm4hep::Vector3d VTXdigi_Allpix2::Vector3dConvert(dd4hep::rec::Vector3D vec) const {
  // return edm4hep::Vector3d(vec[0]/dd4hep::mm, vec[1]/dd4hep::mm, vec[2]/dd4hep::mm);
  return edm4hep::Vector3d(vec[0], vec[1], vec[2]);
}

void VTXdigi_Allpix2::InitialSensorSizeCheck(const edm4hep::SimTrackerHit& simHit) const {
  /** Check that the pixel pitch and count match the sensor size in the geometry
   * 
   * This is only done once for each layer, for the first simHit in the first event.
   * If it fails, an error is raised and no hits are digitized.
   */
  
  const std::uint64_t cellID = simHit.getCellID(); 
  int layer = m_cellIDdecoder->get(cellID, "layer");
  const dd4hep::rec::ISurface* surface = m_surfaceMap->find(cellID)->second;
  debug() << "STARTING initial sensor size check for layer " << layer << endmsg;

  const double length_u = surface->length_along_u() * 10; // convert to mm
  const double length_v = surface->length_along_v() * 10; // convert to mm

  // check sensor size in u, v
  if ( abs(length_u - m_pixelPitch_u[layer]*m_pixelCount_u[layer]) > 0.00001 ) { // numverical tolerance of 0.01 um
    error() << "Sensor size in u (" << length_u << " mm) does not match pixel pitch * count (" << m_pixelPitch_u[layer] << " mm * " << m_pixelCount_u[layer] << ") = " << (m_pixelPitch_u[layer] * m_pixelCount_u[layer]) << " mm for layer " << layer << " (defined in options file). Abort." << endmsg;
    throw std::runtime_error("VTXdigi_Allpix2::InitialSensorSizeCheck: Sensor size in u does not match pixel pitch * count");
  }
  if ( abs(length_v - m_pixelPitch_v[layer]*m_pixelCount_v[layer]) > 0.00001 ) { // numverical tolerance of 0.01 um
    error() << "Sensor size in v (" << length_v << " mm) does not match pixel pitch * count (" << m_pixelPitch_v[layer] << " mm * " << m_pixelCount_v[layer] << ") = " << (m_pixelPitch_v[layer] * m_pixelCount_v[layer]) << " mm for layer " << layer << " (defined in options file). Abort." << endmsg;
    throw std::runtime_error("VTXdigi_Allpix2::InitialSensorSizeCheck: Sensor size in v does not match pixel pitch * count");
  }
  
  m_initialSensorSizeCheckPassed.insert(layer);
  debug() << "PASSED initial sensor size check for layer " << layer << ": " << endmsg;
  return;
}
dd4hep::rec::Vector3D VTXdigi_Allpix2::GetPixelCenter_Local(const int& i_u, const int& i_v, const int& layer, const dd4hep::rec::ISurface& surface) const {
  // returns the position of the center of pixel i_u, i_v) in the global
  // u - short ARCADIA axis, corresponds to x in sensor local frame
  // v - long ARCADIA axis, corresponds to y in sensor local frame

  const double length_u = surface.length_along_u() * 10; // convert to mm
  const double length_v = surface.length_along_v() * 10; // convert to mm

  double posU = (-0.5*length_u) + (i_u + 0.5)*m_pixelPitch_u[layer]; // in mm
  double posV = (-0.5*length_v) + (i_v + 0.5)*m_pixelPitch_v[layer]; // in mm

  return dd4hep::rec::Vector3D(posU, posV, 0); 
}

bool VTXdigi_Allpix2::ApplySimHitCuts (int layer, double eDeposition) const {
  // DISMISS if layer is not in the list of layers to digitize
  if (m_layersToDigitize.value().size()>0) {
    // only digitize hits in the layers specified by the user
    if (std::find(m_layersToDigitize.value().begin(), m_layersToDigitize.value().end(), layer) == m_layersToDigitize.value().end()) {
      debug() << "   - DISMISSED SimHit in layer " << layer << ". (not in the list of layers to digitize)" << endmsg;
      return false;
    }
  }
  // DISMISS if outside minimum energy cut
  if (eDeposition < m_minDepositedEnergy) {
    debug() << "   - DISMISSED simHit (energy below cut, " << eDeposition << " keV)" << endmsg;
    return false;
  }
  return true;
}

std::tuple<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D> VTXdigi_Allpix2::GetSimHitPath (const edm4hep::SimTrackerHit& simHit) const {
  /** Get the path direction (unit vector), entry point (local frame), segmentEdep, segmentLength and segmentNumber for a given simHit
   * 
   */

  // Get the global position of the hit (defined by default in Geant4 as the mean between the entry and exit point in the active material) and apply unit transformation (translation matrix is stored in cm)

  double simHitGlobalCentralPos[3] = {simHit.getPosition().x, simHit.getPosition().y, simHit.getPosition().z}; // given in mm

  double simHitGlobalMomentum[3] = {simHit.getMomentum().x * dd4hep::GeV, simHit.getMomentum().y * dd4hep::GeV, simHit.getMomentum().z * dd4hep::GeV};

  // convert to local detector frame
  const dd4hep::DDSegmentation::CellID& cellID = simHit.getCellID();

  TGeoHMatrix transformationMatrix = GetTransformationMatrix(cellID);
  double simHitLocalCentralPos[3] = {0, 0, 0};
  double simHitLocalMomentum[3] = {0, 0, 0};

  // get the simHit coordinate in cm in the sensor reference frame

  transformationMatrix.MasterToLocal(simHitGlobalCentralPos, simHitLocalCentralPos);
  transformationMatrix.MasterToLocalVect(simHitGlobalMomentum, simHitLocalMomentum);

  // create vector of direction, normalized to path length
  dd4hep::rec::Vector3D simHitPath(simHitLocalMomentum[0], simHitLocalMomentum[1], simHitLocalMomentum[2]);
  simHitPath = (simHit.getPathLength() / simHitPath.r())*simHitPath; // normalize to path length

  // get vector of entry point position
  dd4hep::rec::Vector3D simHitEntryPoint(simHitLocalCentralPos[0], simHitLocalCentralPos[1], simHitLocalCentralPos[2]); // This is the central point of the path as of now
  simHitEntryPoint = simHitEntryPoint - 0.5*simHitPath; // entry point is half a path-length upstream of central position

  // TODO: Maybe add a check that the point is inside the material or go to the closest material and do the same for exit point and then redefine the length (Jessy had this comment already) ~ Jona, 2025-09

  return std::make_tuple(simHitPath, simHitEntryPoint);
}

int VTXdigi_Allpix2::GetBinIndex(float x, float binX0, float binWidth, int binN) const {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  int binIndex = static_cast<int>((x - binX0) / binWidth);
  // debug() << "         - GetBinIndex: x = " << x << ", binX0 = " << binX0 << ", binWidth = " << binWidth << ", binN = " << binN << " => binIndex = " << binIndex << endmsg;
  if (binIndex < 0 || binIndex >= binN) 
    return -1;
  return binIndex;
} // GetBinIndex()

void VTXdigi_Allpix2::CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const double eDeposition, const edm4hep::Vector3d& position) const {

  auto digiHit = digiHits.create();
  digiHit.setCellID(simHit.getCellID());
  digiHit.setEDep(eDeposition);
  digiHit.setPosition(position);
  // TODO: check if position is within sensor bounds & force it onto sensor surface ~ Jona 2025-09
  digiHit.setTime(simHit.getTime());

  auto digiHitLink = digiHitsLinks.create();
  digiHitLink.setFrom(digiHit);
  digiHitLink.setTo(simHit);

  // debug() << "CREATED digiHit with ID " << digiHit.getCellID() << " and Edep " << eDeposition << " keV" << endmsg;
  return;
}
void VTXdigi_Allpix2::CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const double eDeposition, const dd4hep::rec::Vector3D& position) const {
  // overload to allow passing dd4hep::rec::Vector3D as position. ~ Jona 2025-09

  CreateDigiHit(simHit, digiHits, digiHitsLinks, eDeposition, edm4hep::Vector3d(position[0], position[1], position[2]));
  return;
}

