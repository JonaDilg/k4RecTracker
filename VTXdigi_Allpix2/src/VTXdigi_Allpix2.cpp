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
                       {KeyValues("SimTrackHitCollectionName", {"SimTrackerHits"}),
                        KeyValues("HeaderName", {"EventHeader"}),},
                       {KeyValues("TrackerHitCollectionName", {"VTXTrackerHits"}),
                        KeyValues("SimTrkHitRelationsCollection", {"VTXTrackerHitRelations"})}) {
  info() << "Constructed successfully" << endmsg;
}

StatusCode VTXdigi_Allpix2::initialize() {
  // i think this is executed after the Gaudi properties are loaded (different to the constructor, where they are not necessarily loaded yet). ~ Jona 2025-09
  info() << "INITIALIZING ..." << endmsg;

  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }  

  m_geometryService = serviceLocator()->service(m_geometryServiceName);
  if (!m_geometryService) {
    error() << "Unable to retrieve the GeoSvc. Abort." << endmsg;
    return StatusCode::FAILURE;
  }
  
  std::string cellIDstr = m_geometryService->constantAsString(m_encodingStringVariable.value());
  m_cellIDdecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDstr);
  if (!m_cellIDdecoder) {
    error() << "Unable to retrieve the cellID decoder. Abort." << endmsg;
    return StatusCode::FAILURE;
  }
  
  const dd4hep::Detector* detector = m_geometryService->getDetector();
  if (!detector) {
    error() << " - Unable to retrieve the DD4hep detector from GeoSvc. Abort." << endmsg;
    return StatusCode::FAILURE;
  }
  
  const dd4hep::rec::SurfaceManager* simSurfaceManager = detector->extension<dd4hep::rec::SurfaceManager>();
  if (!simSurfaceManager) {
    error() << " - Unable to retrieve the SurfaceManager from the DD4hep detector. Abort." << endmsg;
    return StatusCode::FAILURE;
  }
  
  m_simSurfaceMap = simSurfaceManager->map(m_subDetName.value());
  if (!m_simSurfaceMap) {
    error() << " - Unable to retrieve the simSurface map for subdetector " << m_subDetName.value() << ". Abort." << endmsg;
    return StatusCode::FAILURE;
  }
  
  m_volumeManager = detector->volumeManager();
  if (!m_volumeManager.isValid()) {
    error() << " - Unable to retrieve the VolumeManager from the DD4hep detector. Abort." << endmsg;
    return StatusCode::FAILURE;
  }

  /* subDetector not needed as of now. Keep for future reference ~ Jona 2025-10 */
  // const dd4hep::DetElement subDetector = detector->detector(m_subDetName.value());
  // if (!subDetector.isValid()) {
  //   error() << " - Unable to retrieve the DetElement for subdetector " << m_subDetName.value() << ". Abort." << endmsg;
  //   return StatusCode::FAILURE;
  // }


  debug() << " - Successfully retrieved all necessary services and detector elements, starting to check Gaudi properties." << endmsg;

  // check some config parameters
  if (m_maxClusterSize.value().size() != 2) {
    error() << "Property MaximumClusterSize must contain exactly two entries: [max in u, max in v]. Abort." << endmsg;
    return StatusCode::FAILURE;
  }
  for (int i = 0; i < 2; ++i) {
    if (m_maxClusterSize.value().at(i) < 1) {
      error() << "Property MaximumClusterSize entries must be >= 1. Abort." << endmsg;
      return StatusCode::FAILURE;
    }
    if ((m_maxClusterSize.value().at(i)-1) % 2 != 0) {
      error() << "Property MaximumClusterSize entries must be odd numbers. Abort." << endmsg;
      return StatusCode::FAILURE;
    }
  }

  /* pixel pitch and count are given as either
  *  - a single value, which is then applied to all layers
  *    -> rewrite the vector to contain the same value for all layers (simplifies code later)
  *  - a vector of values, one per layer
  *    -> check that the number of entries matches the number of layers in the geometry
  */ 

  // check that m_layersToDigitize is consistent with m_layerCount
  if (!m_layersToDigitize.empty()) {
    if (m_layersToDigitize.size() != static_cast<size_t>(m_layerCount.value())) {
      error() << "LayersToDigitize size does not match LayerCount. Abort. (LayersToDigitize is optional, set [] to digitize all layers)" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // create layer number to internal index map
  m_layerToIndex.clear();
  for (size_t layerIndex = 0; layerIndex < m_layersToDigitize.size(); layerIndex++) {
    int layer = m_layersToDigitize.value().at(layerIndex);
    m_layerToIndex.insert({layer, static_cast<int>(layerIndex)});
  }
  verbose() << " - Digitizing the following layers: " << endmsg;
  for (const auto& [layer, index] : m_layerToIndex) {
    verbose() << "   - Layer " << layer << " (index " << index << ")" << endmsg;
  }

  if (m_pixelPitch_u.value().size() == 1 && m_pixelPitch_v.value().size() == 1 && m_pixelCount_u.value().size() == 1 && m_pixelCount_v.value().size() == 1 && m_sensorThickness.value().size() == 1) {
    // single value for all layers, rewrite to vector of correct size
    m_pixelPitch_u.value().resize(m_layerCount, m_pixelPitch_u.value().at(0));
    m_pixelPitch_v.value().resize(m_layerCount, m_pixelPitch_v.value().at(0));
    m_pixelCount_u.value().resize(m_layerCount, m_pixelCount_u.value().at(0));
    m_pixelCount_v.value().resize(m_layerCount, m_pixelCount_v.value().at(0));
    m_sensorThickness.value().resize(m_layerCount, m_sensorThickness.value().at(0));
    verbose() << " - found single values for pixel pitch and pixel count, applying them to all layers" << endmsg;
  }
  else if (m_pixelPitch_u.value().size() == m_layerCount && m_pixelPitch_v.value().size() == m_layerCount && m_pixelCount_u.value().size() == m_layerCount && m_pixelCount_v.value().size() == m_layerCount && m_sensorThickness.value().size() == m_layerCount) {
      // one entry per layer
      verbose() << " - found pixel pitch and pixel count values for all " << m_layerCount << " layers" << endmsg;
    }
  else {
    error() << "Pixel pitch and count, and sensor thickness must be given either as a single value for all layers (in brackets though: []) or as a vector with one entry per layer. Abort." << endmsg;
    return StatusCode::FAILURE;
  }

  verbose() << " - Using the following sensor parameters for each layer:" << endmsg;
  for (const auto& [layer, index] : m_layerToIndex) {
    verbose() << "   - Layer " << layer << " (index " << index << "): pitch_u = " << m_pixelPitch_u.value().at(index) << " mm, pitch_v = " << m_pixelPitch_v.value().at(index) << " mm, count_u = " << m_pixelCount_u.value().at(index) << ", count_v = " << m_pixelCount_v.value().at(index) << ", thickness = " << m_sensorThickness.value().at(index) << " mm" << endmsg;
  }

  // check that given pixel pitch and count match the sensor size in the geometry
  

  // -- import charge sharing kernels --

  verbose() << " - Importing charge sharing kernels..." << endmsg;
  m_chargeSharingKernels = std::make_unique<ChargeSharingKernels>(m_inPixelBinCount_u.value(), m_inPixelBinCount_v.value(), m_inPixelBinCount_w.value(), m_kernelSize.value());

  // TODO: read the kernels from file. Instead, set gaussian kernels for now

  if (!m_globalKernel.value().empty()) {
    debug() << " - Global kernel specified in <config>, applying this for all layers" << endmsg;

    if (m_globalKernel.size() != static_cast<size_t>(m_kernelSize * m_kernelSize)) {
      error() << "Global kernel size does not match KernelSize property. Abort." << endmsg;
      return StatusCode::FAILURE;
    }

    // reorder global kernel to match 
    
    for (int i_u=0; i_u<m_inPixelBinCount_u.value(); i_u++) {
      for (int i_v=0; i_v<m_inPixelBinCount_v.value(); i_v++) {
        for (int i_w=0; i_w<m_inPixelBinCount_w.value(); i_w++) {
          m_chargeSharingKernels->SetKernel(i_u, i_v, i_w, m_globalKernel.value());
        }
      }
    }
  }


  // -- Debugging Histograms --

  if (m_debugHistograms.value()) {

    error() << " - You enabled creating debug histograms by setting `DebugHistograms = True`. This is NOT MULTITHREADING SAFE and will cause crashes if multithreading is used." << endmsg;
    verbose () << " - Creating debug histograms ..." << endmsg;

    // -- 1D Histograms --

    m_histograms.at(hist_simHitE).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, 
        "simHitE", "SimHit Energy;keV", 
        {1000, 0, m_sensorThickness.value().at(0)*2000}}); 
    m_histograms.at(hist_simHitCharge).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this,
        "simHitCharge", "SimHit Raw Charge;dep. charges [e-]",
        {1000, 0, m_sensorThickness.value().at(0)*500000}});

    m_histograms.at(hist_clusterSize_raw).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "ClusterSize_Raw", "Cluster Size of Raw Charges (ie. number of pixels that receive any charge from this simHit);Cluster size [pix]", {50, -0.5, 49.5}});
    m_histograms.at(hist_clusterSize_measured).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "ClusterSize_Measured", "Cluster Size of Measured Charges (ie. number of pixels that are above threshold due to the simHit's deposited charge);Cluster size [pix]", {50, -0.5, 49.5}});

    m_histograms.at(hist_pathLength).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "PathLength", "Path Length in Sensor Active Volume;Path length [um]", {500, 0., m_sensorThickness.value().at(0)*10*1000}}); // in um, max 10x sensor thickness (for very shallow angles)
    m_histograms.at(hist_pathLengthGeant4).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "PathLength-Geant4", "Path Length in Sensor Active Volume\nAs Given by Geant4;Path length [um]", {500, 0., m_sensorThickness.value().at(0)*10*1000}}); // in um, max 10x sensor thickness (for very shallow angles)

    m_histograms.at(hist_EntryPointX).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "EntryPointX", "SimHit Entry Point U (in local sensor frame);u [mm]", {400, -4, 4}});
    m_histograms.at(hist_EntryPointY).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "EntryPointY", "SimHit Entry Point V (in local sensor frame);v [mm]", {2000, -20, 20}});
    m_histograms.at(hist_EntryPointZ).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "EntryPointZ", "SimHit Entry Point W (in local sensor frame);w [mm]", {1000, -300, 300}});

    m_histograms.at(hist_DisplacementU).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "LocalDisplacementU", "Displacement U (local sensor frame): digiHit_u - simHit_u;dU [um]", {800, -200, 200}});
    m_histograms.at(hist_DisplacementV).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "LocalDisplacementV", "Displacement V (local sensor frame): digiHit_v - simHit_v;dV [um]", {800, -200, 200}});
    m_histograms.at(hist_DisplacementR).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "GlobalDisplacementR", "Displacement R (global frame): | digiHit - simHit |;dR [um]", {300, 0., 300}});

    m_histograms.at(hist_HitRawChargeDifference).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "HitRawChargeDifference", "Hit Raw Charge Difference: #e- (digiHits, raw) - # e- (simHit);dep. charges [e-]", {5000, -5000, 5000}});

    m_histograms.at(hist_pixelChargeMatrix_size_u).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "PixelChargeMatrix_Size_u", "Final Size of the Pixel Charge Matrix in u (local);pixel charge matrix size u [pix]", {100, -0.5, 99.5}});
    m_histograms.at(hist_pixelChargeMatrix_size_v).reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "PixelChargeMatrix_Size_v", "Final Size of the Pixel Charge Matrix in v (local);pixel charge matrix size v [pix]", {100, -0.5, 99.5}});

    // -- 2D Histograms --

    m_histograms2d.resize(m_layersToDigitize.size()); // resize vector to hold all histograms
    m_histWeighted2d.resize(m_layersToDigitize.size());

    for (int layer : m_layersToDigitize) {
      int layerIndex = m_layerToIndex.at(layer);
      verbose () << " - Creating 2D histograms for layer " << layer << " (index " << layerIndex << ")" << endmsg;

      m_histograms2d.at(layerIndex).at(hist2d_hitMap_simHits).reset(
        new Gaudi::Accumulators::StaticRootHistogram<2>{this, 
          "Layer"+std::to_string(layer)+"_HitMap_simHits",
          "SimHit Hitmap, Layer " + std::to_string(layer) + ";u [pix]; v [pix]",
          {static_cast<unsigned int>(m_pixelCount_u.value().at(0)), -0.5, static_cast<double>(m_pixelCount_u.value().at(0)+0.5)},
          {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)}
        }
      );
      m_histograms2d.at(layerIndex).at(hist2d_hitMap_digiHits).reset(
        new Gaudi::Accumulators::StaticRootHistogram<2>{this, 
          "Layer"+std::to_string(layer)+"_HitMap_digiHits",
          "DigiHit Hitmap, Layer " + std::to_string(layer) + ";u [pix]; v [pix]",
          {static_cast<unsigned int>(m_pixelCount_u.value().at(0)), -0.5, static_cast<double>(m_pixelCount_u.value().at(0)+0.5)},
          {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)}
        }
      );
      m_histograms2d.at(layerIndex).at(hist2d_hitMap_simHitDebug).reset(
        new Gaudi::Accumulators::StaticRootHistogram<2>{this, 
          "Layer"+std::to_string(layer)+"_HitMap_simHitsDebug",
          "SimHit Debugging Hitmap, Layer " + std::to_string(layer) + " showing hits where local w != -25 um;u [pix]; v [pix]",
          {static_cast<unsigned int>(m_pixelCount_u.value().at(0)), -0.5, static_cast<double>(m_pixelCount_u.value().at(0)+0.5)},
          {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)}
        }
      );
      m_histograms2d.at(layerIndex).at(hist2d_hitMap_digiHitDebug).reset(
        new Gaudi::Accumulators::StaticRootHistogram<2>{this, 
          "Layer"+std::to_string(layer)+"_HitMap_digiHitsDebug",
          "DigiHit Debugging Hitmap, Layer " + std::to_string(layer) + " showing hits where local w != -25 um;u [pix]; v [pix]",
          {static_cast<unsigned int>(m_pixelCount_u.value().at(0)), -0.5, static_cast<double>(m_pixelCount_u.value().at(0)+0.5)},
          {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)}
        }
      );
      m_histograms2d.at(layerIndex).at(hist2d_pathLength_vs_simHit_v).reset(
        new Gaudi::Accumulators::StaticRootHistogram<2>{this, 
          "Layer"+std::to_string(layer)+"_TrackLength_vs_simHit_v",
          "Track Length in Sensor Active Volume vs. SimHit v position (local), Layer " + std::to_string(layer) + ";SimHit v [pix];Track length [um]",
            {static_cast<unsigned int>(m_pixelCount_v.value().at(0)), -0.5, static_cast<double>(m_pixelCount_v.value().at(0)+0.5)},
          {200, 0., 10*m_sensorThickness.value().at(layerIndex)*1000}, // in um
        
        }
      );
      m_histograms2d.at(layerIndex).at(hist_pixelChargeMatrixSize).reset(
        new Gaudi::Accumulators::StaticRootHistogram<2>{this, 
          "Layer"+std::to_string(layer)+"_PixelChargeMatrixSize",
          "Pixel Charge Matrix Size, Layer " + std::to_string(layer) + ";u [pix];v [pix]",
          {100, -0.5, 99.5},
          {100, -0.5, 99.5}
        }
      );


      m_histWeighted2d.at(layerIndex).at(histWeighted2d_averageCluster_rawCharge).reset(
        new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, double>{this,
          "Layer"+std::to_string(layer)+"_AverageCluster_rawCharge",
          "Average Cluster Shape, without Threshold & Noise, Layer " + std::to_string(layer) + "(charge per hit normalised to 1);u [pix]; v [pix]",
          {static_cast<unsigned int>(m_maxClusterSize.value().at(0)),
            -static_cast<double>(m_maxClusterSize.value().at(0))/2,
            static_cast<double>(m_maxClusterSize.value().at(0))/2},
          {static_cast<unsigned int>(m_maxClusterSize.value().at(1)),
            -static_cast<double>(m_maxClusterSize.value().at(1))/2,
            static_cast<double>(m_maxClusterSize.value().at(1))/2}
        }
      );

      m_histWeighted2d.at(layerIndex).at(histWeighted2d_averageCluster_measuredCharge).reset(
        new Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, double>{this,
          "Layer"+std::to_string(layer)+"_AverageCluster_measuredCharge",
          "Average Cluster Shape, after applying Threshold & Noise, Layer " + std::to_string(layer) + "(charge per hit normalised to 1);u [pix]; v [pix]",
          {static_cast<unsigned int>(m_maxClusterSize.value().at(0)),
            -static_cast<double>(m_maxClusterSize.value().at(0))/2,
            static_cast<double>(m_maxClusterSize.value().at(0))/2},
          {static_cast<unsigned int>(m_maxClusterSize.value().at(1)),
            -static_cast<double>(m_maxClusterSize.value().at(1))/2,
            static_cast<double>(m_maxClusterSize.value().at(1))/2}
        }
      );
    }
  } // if debug histograms

  // -- Debugging CSV output --
  if (!m_debugCsvFileName.value().empty()) {
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
  } else { 
    m_debugCsv = false;
    verbose() << " - Debug CSV output disabled" << endmsg;
  }

  // TODO: check that pixel pitch and count match sensor size in geometry (I am not sure how to get the sensor size from the geometry though, does seem to be a bit more general, subDetectors have children that might or might not be layers) ~ Jona 2025-09
  // TODO load lookup table for charge transport and diffusion
  // TODO check that k4geo has same pitch as lookup table 
  info() << " - Initialized successfully" << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode VTXdigi_Allpix2::finalize() {
  info() << "FINALIZING ..." << endmsg;

  if (m_debugCsvFile.is_open()) {
    m_debugCsvFile.close();
    verbose() << " - Closed debug CSV file" << endmsg;
  }

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
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits rejected (surface distance too large)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRejected_SurfaceDistToLarge.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits rejected (outside sensor)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRejected_OutsideSensor.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Simhits accepted"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsAccepted.value() << " |" << endmsg;
  const long int simHitsRejected = (
    m_counter_simHitsRejected_LayerNotToBeDigitized.value() + 
    m_counter_simHitsRejected_ChargeCut.value() + 
    m_counter_simHitsRejected_SurfaceDistToLarge.value() + 
    m_counter_simHitsRejected_OutsideSensor.value() );
  if (m_counter_simHitsRead.value() != simHitsRejected + m_counter_simHitsAccepted.value())
    warning() << " | Number of accepted and rejected simHits does not add up to total number of processed simHits!" << endmsg;
  else info() << " | - simHit numbers add up." << endmsg;

  info() << " | " << std::setw(colWidths[0]) << std::left << "Calculated path did not pass through the sensor active volume"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_acceptedButNoSegmentsInSensor.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "DigiHits created"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_digiHitsCreated.value() << " |" << endmsg;
  


  verbose() << " - finalized successfully" << endmsg;
  return StatusCode::SUCCESS;
} 


// -- event loop -- 

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> VTXdigi_Allpix2::operator()
  (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {
  const long int eventNumber = headers.at(0).getEventNumber();

  info() << "PROCESSING event. run " << headers.at(0).getRunNumber() << "; event " << eventNumber << "; with " << simHits.size() << " simHits" << endmsg;
  ++m_counter_eventsRead;

  // early sanity checks to avoid segfaults from null pointers
  if (!m_chargeSharingKernels) 
    throw GaudiException("ChargeSharingKernels is null in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::operator()", StatusCode::FAILURE);
  if (simHits.size()==0) {
    debug() << " - No SimTrackerHits in collection, returning empty output collections" << endmsg;
    ++m_counter_eventsRejected_noSimHits;
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());
  }
  ++m_counter_eventsAccepted;

  // create output collections
  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitsLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();
  unsigned int nHitsRead=0, nHitsCreated=0, nHitsRejected=0; // per-event counters

  // initialize random number generator for this event
  const u_int32_t rngSeed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  TRandom2 rngEngine = TRandom2(rngSeed);
  verbose() << " - RNG engine initialized, using seed " << rngSeed << endmsg;

  // loop over sim hits, digitize them, create output collections
  for (const auto& simHit : simHits) {
    nHitsRead++;
    ++m_counter_simHitsRead;

    HitInfo hitInfo; // struct to hold information about the simHit (and be passed to functions)
    
    // gather information about the simHit
    hitInfo.cellID = simHit.getCellID(); // TODO: apply 64-bit length mask as done in DDPlanarDigi.cpp? ~ Jona 2025-09
    const dd4hep::rec::Vector3D simHitGlobalPos = convertVector(simHit.getPosition());
    
    hitInfo.simLocalPos = transformGlobalToLocal(simHitGlobalPos, hitInfo.cellID);
    
    // DISMISS if layer is not in the list of layers to digitize. Cant be done with other cuts later, would lead to out-of-bounds access in m_pixelCount_u etc
    if (m_layersToDigitize.value().size()>0) {
      const int layer = m_cellIDdecoder->get(hitInfo.cellID, "layer");
      if (std::find(m_layersToDigitize.value().begin(), m_layersToDigitize.value().end(), layer) == m_layersToDigitize.value().end()) {
        verbose() << "   - DISMISSED SimHit in layer " << layer << ". (not in the list of layers to digitize)" << endmsg;
        ++m_counter_simHitsRejected_LayerNotToBeDigitized;
        nHitsRejected++;
        continue;
      }
    }
    
    debug() << "  - FOUND SimHit in layer " << m_cellIDdecoder->get(hitInfo.cellID, "layer") << endmsg;
    hitInfo.layerIndex = m_layerToIndex.at(m_cellIDdecoder->get(hitInfo.cellID, "layer")); // will throw if layer not in map, but this is already checked above. Use index instead of layer to avoid segfaults if layers are not numbered consecutively

    hitInfo.simCharge = simHit.getEDep() * (dd4hep::GeV / dd4hep::keV) * m_chargePerkeV; // in electrons
    hitInfo.simPathLength = simHit.getPathLength(); // in mm
    
    if (!m_simSurfaceMap) throw GaudiException("SimSurface map is null in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::operator()", StatusCode::FAILURE);
    const auto itSimSurface = m_simSurfaceMap->find(hitInfo.cellID);
    if (itSimSurface == m_simSurfaceMap->end()) throw GaudiException("No simSurface found for cellID " + std::to_string(hitInfo.cellID) + " in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::operator()", StatusCode::FAILURE);
    hitInfo.simSurface = itSimSurface->second;
    if (!hitInfo.simSurface) throw GaudiException("SimSurface pointer for cellID " + std::to_string(hitInfo.cellID) + " is null in operator(). Did initialize() succeed?", "VTXdigi_Allpix2::operator()", StatusCode::FAILURE);
    
    
    // get sensor size in local frame
    hitInfo.length[0] = hitInfo.simSurface->length_along_u() * 10; // convert to mm
    hitInfo.length[1] = hitInfo.simSurface->length_along_v() * 10; // convert to mm

    hitInfo.thickness = m_sensorThickness.value().at(hitInfo.layerIndex); // default value, will be overwritten below

    hitInfo.pixPitch[0] = m_pixelPitch_u.value().at(hitInfo.layerIndex);
    hitInfo.pixPitch[1] = m_pixelPitch_v.value().at(hitInfo.layerIndex);

    hitInfo.pixCount[0] = m_pixelCount_u.value().at(hitInfo.layerIndex);
    hitInfo.pixCount[1] = m_pixelCount_v.value().at(hitInfo.layerIndex);

    
    { // check if pitch & count match sensor size (doing this repeatedly for every simHit is not optimal, but anything significantly better I can come up with is not multithreading-safe ~ Jona 2025-10)
      const float length_u = hitInfo.pixPitch[0] * hitInfo.pixCount[0];
      if ( abs(length_u - hitInfo.length[0]) > m_numericLimit_float )
      throw std::runtime_error("VTXdigi_Allpix2::operator(): Sensor size from geometry u does not match pixel pitch * count in layer " + std::to_string(m_cellIDdecoder->get(hitInfo.cellID, "layer")));
      const float length_v = hitInfo.pixPitch[1] * hitInfo.pixCount[1];
      if ( abs(length_v - hitInfo.length[1]) > m_numericLimit_float )
      throw std::runtime_error("VTXdigi_Allpix2::operator(): Sensor size from geometry v does not match pixel pitch * count in layer " + std::to_string(m_cellIDdecoder->get(hitInfo.cellID, "layer")));
    }
    
    const auto& [simHitEntryPos, simHitPath] = ConstructSimHitPath(hitInfo, simHit); // both are of type dd4hep::rec::Vector3D, given in mm
    
    { // debug statements
      debug() << "   - dep. charge = " << hitInfo.simCharge << " e-, path length = " << simHitPath.r() << " mm, Geant4 path length = " << hitInfo.simPathLength << " mm" << endmsg;
      verbose() << "   - Position (global) " << simHit.getPosition().x << " mm, " << simHit.getPosition().y << " mm, " << simHit.getPosition().z << " mm" << endmsg;
      debug() << "   - Position (local) " << hitInfo.simLocalPos.x() << " mm, " << hitInfo.simLocalPos.y() << " mm, " << hitInfo.simLocalPos.z() << " mm (in sensor frame)" << endmsg;
      verbose() << "   - Distance to surface " << hitInfo.simSurface->distance(dd4hep::mm * simHitGlobalPos) << " mm" << endmsg; // dd4hep expects cm, simHitGlobalPosition is in mm. dd4hep::cm=0.1
      verbose() << "   - Entry point (sensor local) : " << simHitEntryPos.x() << " mm, " << simHitEntryPos.y() << " mm, " << simHitEntryPos.z() << " mm" << endmsg;
    }
    
    // apply cuts
    if (!ApplySimHitCuts(hitInfo, simHitEntryPos, simHitPath, simHitGlobalPos)) {
      nHitsRejected++; // the specific m_counter_... gets incremented inside ApplySimHitCuts()
      continue;
    }

    /* sanity check: pixel pitch * pixel count must match sensor size from geometry. 
     * Do this after applying cuts, as non-digitized layers would lead to out-of-bounds access in m_pixelCount_u etc.
     * Doing this for every simHit is not optimal (because pixPitch and pixCount per layer are constant across all events)
     * but accessing the surface-based sensor size is not easy to implement without using a cellID that is known to lie on that layer. */
    if ( abs(hitInfo.pixPitch[0] * hitInfo.pixCount[0] - hitInfo.length[0]) > m_numericLimit_float
      || abs(hitInfo.pixPitch[1] * hitInfo.pixCount[1] - hitInfo.length[1]) > m_numericLimit_float) {
      throw GaudiException("Gaudi properties (pixelPitch * pixelCount) does not match sensor size from detector geometry in layer " + std::to_string(m_cellIDdecoder->get(hitInfo.cellID, "layer")) + " in operator().", "VTXdigi_Allpix2::operator()", StatusCode::FAILURE);
    }
    
    // get more info about the simHit path segmentation
    hitInfo.nSegments = int( hitInfo.simPathLength / m_targetPathSegmentLength ); // both in mm
    hitInfo.nSegments = std::max(1, hitInfo.nSegments); // want at least one segment
    float segmentLength = hitInfo.simPathLength / hitInfo.nSegments;
    hitInfo.segmentCharge = hitInfo.simCharge / hitInfo.nSegments;
    debug() << "   - ACCEPTED SimHit in layer " << m_cellIDdecoder->get(hitInfo.cellID, "layer") << ". Charge =" << hitInfo.simCharge << " e-. Starting loop over segments" << endmsg;
    ++m_counter_simHitsAccepted;
    // cannot increase counter for accepted simHits just yet, only after segments loop succeeded
    
    
    // -- loop over the segments, assign charges to pixels according to kernel--
    verbose() << "   - Segmenting path of length " << hitInfo.simPathLength << " mm into " << hitInfo.nSegments << " segments of length " << segmentLength << " mm" << endmsg;
    
    hitInfo.pixCount[0] = m_pixelCount_u.value().at(hitInfo.layerIndex);
    hitInfo.pixCount[1] = m_pixelCount_v.value().at(hitInfo.layerIndex);

    // int i_u_simHit, i_v_simHit;
    const auto& [i_u_simHit, i_v_simHit] = computePixelIndices(hitInfo.simLocalPos, hitInfo.layerIndex, hitInfo.length[0], hitInfo.length[1]);
    PixelChargeMatrix pixelChargeMatrix(
      i_u_simHit, i_v_simHit,
      m_maxClusterSize.value().at(0), m_maxClusterSize.value().at(1)
      );

    int i_u, i_v; // segment being processed: pixel indices
    int j_u, j_v, j_w; // in-pixel indices
    int i_u_next, i_v_next; // next segment
    int j_u_next, j_v_next, j_w_next;

    // first segment (segment=0) is treated outside the loop, to initialize j_u_previous
    verbose() << " Pre-loading first segment" << endmsg;
    std::tie(i_u, i_v, j_u, j_v, j_w) = computeSegmentIndices(hitInfo, simHitEntryPos, simHitPath, 0);
    int nSegmentsInBin = 1; // number of segments in the same in-pixel bin (used to apply kernel only once per bin, not per segment)

    for (int segment=1; segment<hitInfo.nSegments; segment++) {

      std::tie(i_u_next, i_v_next, j_u_next, j_v_next, j_w_next) = computeSegmentIndices(hitInfo, simHitEntryPos, simHitPath, segment);

      if ((i_u_next == i_u && i_v_next == i_v && j_u_next == j_u && j_v_next == j_v && j_w_next == j_w) && segment < hitInfo.nSegments-1) {
        nSegmentsInBin++;
        verbose() << "       - Segment lies in the same pixel and in-pixel bin as previous segment, increasing counter to " << nSegmentsInBin << endmsg;
        continue; // still in the same pixel and in-pixel bin
      } // if sement in same bin
      else {
        verbose() << "       - Crossed bin-boundary or reached last segment. Sharing " << hitInfo.segmentCharge*nSegmentsInBin << " from previous bin." << endmsg;

        if (i_u == -1 && i_v == -1 && j_u == -1 && j_v == -1 && j_w == -1) {
          warning() << "Applying Kernel: Segments lie outside sensor area. Dismissing." << endmsg;
          continue; // segment is outside sensor area
        }
        // new pixel or in-pixel bin, apply kernel for previous bin

        /* each kernel entry shares charge from the source-pixel i_u, i_v to one target-pixel.
         * The kernel is centered on the source pixel, so loop over all target pixels covered by the kernel.
         * i_x_previous defines the source-pixel */

        int i_u_target, i_v_target;

        for (int i_m = -1*(m_kernelSize.value()-1)/2; i_m<=(m_kernelSize.value()-1)/2; i_m++) {
          i_u_target = i_u + i_m;
          if (i_u_target<0 || i_u_target>=hitInfo.pixCount[0]) continue;

          for (int i_n = -1*(m_kernelSize.value()-1)/2; i_n<=(m_kernelSize.value()-1)/2; i_n++) {
            i_v_target = i_v + i_n;
            if (i_v_target<0 || i_v_target>=hitInfo.pixCount[1]) continue;

            if (j_u == -1 || j_v == -1 || j_w == -1) {
              warning() << "Applying Kernel: Segment lies outside sensor area. Dismissing." << endmsg;
              continue; 
            }

            try {
              const float kernelEntry = m_chargeSharingKernels->GetWeight(j_u, j_v, j_w, i_m, i_n);
              if (kernelEntry < m_numericLimit_float) continue; // skip zero entries
              pixelChargeMatrix.FillRawCharge(
                i_u_target, i_v_target,
                kernelEntry * nSegmentsInBin * hitInfo.segmentCharge);
            } catch (const std::exception& e) {
              warning() << "Exception while accessing kernel entry and saving charge for pixel (" << i_u_target << "," << i_v_target << "), in pix bin(" << j_u << "," << j_v << "," << j_w << "): " << e.what() << ". Skipping this bin." << endmsg;
              continue;
            }
          }
        }

        nSegmentsInBin = 1;
        i_u = i_u_next;
        i_v = i_v_next;
        j_u = j_u_next;
        j_v = j_v_next;
        j_w = j_w_next;
      } // if segment in different bin
    } // loop over segments
    
    //  -- find fired pixels, create digiHit --
    int nPixelsFired = 0;
    
    debug() << "   - Processed all segments. Looping over pixelChargeMatrix." << endmsg;

    // Generate noise for each pixel (only now, after we know which )
    pixelChargeMatrix.GenerateNoise(rngEngine, m_electronicNoise);

    for (i_u = pixelChargeMatrix.GetRangeMin_u(); i_u <= pixelChargeMatrix.GetRangeMax_u(); ++i_u) {

      for (i_v = pixelChargeMatrix.GetRangeMin_v(); i_v <= pixelChargeMatrix.GetRangeMax_v(); ++i_v) {

        float pixelChargeRaw = pixelChargeMatrix.GetRawCharge(i_u, i_v);
        if (pixelChargeRaw < m_numericLimit_float) 
          continue; 
        /* skip pixels with zero or negative charge BEFORE noise addition, st. each pixel without charge is treated the same
         * (otherwise noise would happen only around simHits, because pixelChargeMatrix only holds pixels around simHit. This would bias the results.) */

        float pixelChargeMeasured = pixelChargeRaw + pixelChargeMatrix.GetNoise(i_u, i_v);

        if (pixelChargeMeasured < m_pixelThreshold) {
          verbose() << "     - Pixel (" << i_u << ", " << i_v << ") received a measured/raw charge of " << pixelChargeMeasured << " / " << pixelChargeRaw << " e-. This is below threshold (" << m_pixelThreshold << " e-). Discarding pixel." << endmsg;
          continue; // pixel did not pass threshold
        }

        nPixelsFired++;
        dd4hep::rec::Vector3D pixelCenterLocal = computePixelCenter_Local(i_u, i_v, hitInfo.layerIndex, *hitInfo.simSurface);
        dd4hep::rec::Vector3D pixelCenterGlobal = transformLocalToGlobal(pixelCenterLocal, hitInfo.cellID);

        debug() << "     - Pixel (" << i_u << ", " << i_v << ") at (" << pixelCenterLocal.x() << ", " << pixelCenterLocal.y() << ", " << pixelCenterLocal[2] << ") mm received a measured/raw charge of " << pixelChargeMeasured << "/" << pixelChargeRaw << " e-, center at global position " << pixelCenterGlobal[0] << " mm, " << pixelCenterGlobal[1] << " mm, " << pixelCenterGlobal[2] << " mm" << endmsg;
        CreateDigiHit(simHit, digiHits, digiHitsLinks, pixelCenterGlobal, pixelChargeMeasured);
        nHitsCreated++;
        ++m_counter_digiHitsCreated;
      }
    } // loop over pixels, create digiHits

    if (nPixelsFired <= 0) {
      ++m_counter_acceptedButNoSegmentsInSensor;
      debug() << "   - No pixels fired for this simHit (no segments lie inside the sensor volume)" << endmsg;
      continue;
    }
    verbose() << "   - From this simHit, " << nPixelsFired << " pixels received charge." << endmsg;

    // -- fill debugging histograms --

    if (m_debugHistograms) {
      verbose() << "   - Filling 1D histograms" << endmsg;
      ++(*m_histograms.at(hist_simHitCharge))[hitInfo.simCharge]; // in e
      ++(*m_histograms.at(hist_simHitE))[hitInfo.simCharge / m_chargePerkeV]; // in keV
      ++(*m_histograms.at(hist_EntryPointX))[simHitEntryPos.x()]; // in mm
      ++(*m_histograms.at(hist_EntryPointY))[simHitEntryPos.y()]; // in mm
      ++(*m_histograms.at(hist_EntryPointZ))[simHitEntryPos.z()*1000]; // convert mm to um
      ++(*m_histograms.at(hist_pathLength))[simHitPath.r()*1000]; // convert mm to um
      ++(*m_histograms.at(hist_pathLengthGeant4))[hitInfo.simPathLength*1000]; // convert mm to um
      ++(*m_histograms.at(hist_HitRawChargeDifference))[ pixelChargeMatrix.GetTotalRawCharge() - hitInfo.simCharge]; // in e-, only if at least one pixel fired
      
      int pix_u, pix_v;
      std::tie(pix_u, pix_v) = computePixelIndices(hitInfo.simLocalPos, hitInfo.layerIndex, hitInfo.length[0], hitInfo.length[1]);
      verbose() << "   - Filling 2D histograms for layer " << m_cellIDdecoder->get(hitInfo.cellID, "layer") << " (index " << hitInfo.layerIndex << ")" << endmsg;
      ++(*m_histograms2d.at(hitInfo.layerIndex).at(hist2d_hitMap_simHits))[{pix_u, pix_v}]; // in mm
      ++(*m_histograms2d.at(hitInfo.layerIndex).at(hist2d_pathLength_vs_simHit_v))[{pix_v, simHitPath.r()*1000}]; // in um and mm
      
      ++(*m_histograms.at(hist_pixelChargeMatrix_size_u))[pixelChargeMatrix.GetSize_u()];
      ++(*m_histograms.at(hist_pixelChargeMatrix_size_v))[pixelChargeMatrix.GetSize_v()];
      
      ++(*m_histograms2d.at(hitInfo.layerIndex).at(hist_pixelChargeMatrixSize))[{pixelChargeMatrix.GetSize_u(), pixelChargeMatrix.GetSize_v()}];
      
      
      if (hitInfo.debugFlag)
      ++(*m_histograms2d.at(hitInfo.layerIndex).at(hist2d_hitMap_simHitDebug))[{pix_u, pix_v}];
      
      int nPixelsReceivedCharge = 0;

      // fill in-sensor-dependent histograms (not efficient with histograms enabled, but makes code slightly faster when histograms are disabled)
      for (i_u = pixelChargeMatrix.GetRangeMin_u(); i_u <= pixelChargeMatrix.GetRangeMax_u(); ++i_u) {
        for (i_v = pixelChargeMatrix.GetRangeMin_v(); i_v <= pixelChargeMatrix.GetRangeMax_v(); ++i_v) {
          


          float pixelChargeRaw = pixelChargeMatrix.GetRawCharge(i_u, i_v);
          if (pixelChargeRaw < m_numericLimit_float) continue; 

          nPixelsReceivedCharge++;

          (*m_histWeighted2d.at(hitInfo.layerIndex).at(histWeighted2d_averageCluster_rawCharge))[{i_u - i_u_simHit, i_v - i_v_simHit}] += pixelChargeRaw / pixelChargeMatrix.GetTotalRawCharge(); // in e-
          
          float pixelChargeMeasured = pixelChargeRaw + pixelChargeMatrix.GetNoise(i_u, i_v);
          if (pixelChargeMeasured < m_pixelThreshold) continue;
          
          dd4hep::rec::Vector3D pixelCenterLocal = computePixelCenter_Local(i_u, i_v, hitInfo.layerIndex, *hitInfo.simSurface);
          dd4hep::rec::Vector3D pixelCenterGlobal = transformLocalToGlobal(pixelCenterLocal, hitInfo.cellID);
          
          ++ (*m_histograms.at(hist_DisplacementU))[ (pixelCenterLocal.x() - hitInfo.simLocalPos.x()) * 1000]; // convert mm to um
          ++ (*m_histograms.at(hist_DisplacementV))[ (pixelCenterLocal.y() - hitInfo.simLocalPos.y()) * 1000]; // convert mm to um
          ++ (*m_histograms.at(hist_DisplacementR))[sqrt( (pixelCenterGlobal.x() - simHitGlobalPos.x())*(pixelCenterGlobal.x() - simHitGlobalPos.x()) + (pixelCenterGlobal.y() - simHitGlobalPos.y())*(pixelCenterGlobal.y() - simHitGlobalPos.y()) + (pixelCenterGlobal.z() - simHitGlobalPos.z())*(pixelCenterGlobal.z() - simHitGlobalPos.z()) ) * 1000]; // convert mm to um
          
          (*m_histWeighted2d.at(hitInfo.layerIndex).at(histWeighted2d_averageCluster_measuredCharge))[{i_u - i_u_simHit, i_v - i_v_simHit}] += pixelChargeMeasured / pixelChargeMatrix.GetTotalMeasuredCharge() ; // in e-
        }
      }
      ++(*m_histograms.at(hist_clusterSize_raw))[static_cast<double>(nPixelsReceivedCharge)]; // cluster size = number of pixels fired per simHit
      ++(*m_histograms.at(hist_clusterSize_measured))[static_cast<double>(nPixelsFired)]; // cluster size = number of pixels fired per simHit
    }
    
    if (m_debugCsv) {
      int i_u_debug, i_v_debug;
      std::tie(i_u_debug, i_v_debug) = computePixelIndices(hitInfo.simLocalPos, hitInfo.layerIndex, hitInfo.length[0], hitInfo.length[1]);
      if (m_debugCsv)
      appendSimHitToCsv(hitInfo, eventNumber, simHitEntryPos, simHitPath, i_u_debug, i_v_debug);
    }
  } // end loop over sim hits
  
  debug() << "FINISHED event. Processed " << nHitsRead << " simHits. From this, created " << nHitsCreated << " digiHits and dismissed " << nHitsRejected << " simHits." << endmsg;
  
  verbose() << " - Returning collections: digiHits.size()=" << digiHits.size() << ", digiHitsLinks.size()=" << digiHitsLinks.size() << endmsg;
  
  return std::make_tuple(std::move(digiHits), std::move(digiHitsLinks));
} // operator()





// -- Core algorithm functions -- 

bool VTXdigi_Allpix2::ApplySimHitCuts (HitInfo& hitInfo, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const dd4hep::rec::Vector3D& simHitGlobalPos) const {

  // DISMISS if simHitPosition is not on the DD4hep sensor simSurface
  if (m_cutDistanceToSurface) {
    if (abs(hitInfo.simSurface->distance(dd4hep::mm * simHitGlobalPos)) > m_numericLimit_float) {
      verbose() << "   - DISMISSED simHit (is not on the DD4hep sensor simSurface (distance = " << hitInfo.simSurface->distance(dd4hep::mm * simHitGlobalPos) * 1000 << " um)." << endmsg; // convert to um
      ++m_counter_simHitsRejected_SurfaceDistToLarge;
      return false;
    }
  }
  
  // DISMISS if entry point is outside sensor
  if (m_cutPathOutsideSensor) {
    if (abs(simHitEntryPos.z()) > hitInfo.thickness / 2 + m_numericLimit_float) { // entry point is outside sensor thickness
      verbose() << "   - DISMISSED simHit (entry point is outside sensor thickness (local w = " << simHitEntryPos.z()*1000 << " um, sensor thickness = " << hitInfo.thickness*1000 << " um)." << endmsg;
      ++m_counter_simHitsRejected_OutsideSensor;
      return false;
    }
    if (abs(simHitEntryPos.z()+simHitPath.z()) > hitInfo.thickness / 2 + m_numericLimit_float) { // exit point is outside sensor thickness
      verbose() << "   - DISMISSED simHit (exit point is outside sensor thickness (local w = " << (simHitEntryPos.z()+simHitPath.z())*1000 << " um, sensor thickness = " << hitInfo.thickness*1000 << " um)." << endmsg;
      ++m_counter_simHitsRejected_OutsideSensor;
      return false;
    }
  }

  // DISMISS if outside minimum charge cut
  if (hitInfo.simCharge < m_cutDepositedCharge.value()) {
    ++m_counter_simHitsRejected_ChargeCut;
    verbose() << "   - DISMISSED simHit (charge below cut, " << hitInfo.simCharge << " e-)" << endmsg;
    return false;
  }
  
  return true;
}

std::tuple<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D> VTXdigi_Allpix2::ConstructSimHitPath (HitInfo& hitInfo, const edm4hep::SimTrackerHit& simHit) const {
  /** Get the simHitPath (the path of the simulated particle through a sensor) of a given simHit. Describet by the path direction (unit vector) and entry point (local frame).
   *  
   * The way I do this differs from what Jessy does in their vtxDigi Detailed.
   *  My approach ensures that the simHitPath always extends exactly from one sensor surface to the other, and that the simHitPos always lies on this path. 
   * Then, the paths are shortened to the PathLength that Geant4 gives us, if they are longer (the Geant4 pathLength is always longer than the linear approx., because it accounts for B-field and multiple scattering inside the sensor). We need to do this to avoid unphysically long paths, where the charge-per-pathLength would be unphysically tiny
   *  I am not sure if this is better or worse, it definitely produces less problems. 
   *  I think the assumption of a linear path through the sensor is not compatible with using the path-length as normalisation for the path-vector
   */

  // const dd4hep::DDSegmentation::CellID& cellID = simHit.getCellID();
  
  // get simHitMomentum (this vector gives the exact agles of the particle's path through the sensor. We assume that path is linear.). Transform to sensor's local coordinates.
  double simHitGlobalMomentum_double[3] = {simHit.getMomentum().x * dd4hep::GeV, simHit.getMomentum().y * dd4hep::GeV, simHit.getMomentum().z * dd4hep::GeV}; // need floats for the TransfomationMatrix functions

  TGeoHMatrix transformationMatrix = computeTransformationMatrix(hitInfo.cellID);

  double simHitLocalMomentum_double[3] = {0.0, 0.0, 0.0};
  transformationMatrix.MasterToLocalVect(simHitGlobalMomentum_double, simHitLocalMomentum_double);

  dd4hep::rec::Vector3D simHitPath(
    simHitLocalMomentum_double[0],
    simHitLocalMomentum_double[1],
    simHitLocalMomentum_double[2]);

  float scaleFactor = hitInfo.thickness / std::abs(simHitLocalMomentum_double[2]);
  simHitPath = scaleFactor * simHitPath ; // now, simHitPath extends for one sensor surface to the other surface
  // edge case with curlers having horizontal momentum is handled later, by cutting the path to the length supplied by Geant4.
  
  // calculate the path's entry position into the sensor, by placing it such that the path passes through the simHit position
  if (abs(hitInfo.simLocalPos.z()) > (0.5*hitInfo.thickness + m_numericLimit_float)) {
    warning() << "SimHit position is outside the sensor volume (local w = " << hitInfo.simLocalPos.z() << " mm, sensor thickness = " << hitInfo.thickness << " mm). This should never happen. Forcing it to w=0." << endmsg;
    hitInfo.simLocalPos.z() = 0.;
  }
  
  // calculate how far the entry point is from the simHitPos, in terms of the simHitPath
  scaleFactor = 0.;
  if (simHitPath.z() >= 0.) {
    const float shiftDist_w = 0.5 * hitInfo.thickness + hitInfo.simLocalPos.z();
    scaleFactor = shiftDist_w / simHitPath.z();
  }
  else {
    const float shiftDist_w = -0.5 * hitInfo.thickness + hitInfo.simLocalPos.z();
    scaleFactor = shiftDist_w / simHitPath.z();
  }
  
  // Entry position is on either surface of the sensor, upstream from the simHit position (exactly by the fraction we just calculated)
  dd4hep::rec::Vector3D simHitEntryPos = hitInfo.simLocalPos - scaleFactor * simHitPath;

  // Clip the path to the sensor edges in u and v direction (if it passes through the side of the sensor)
  float t_min = 0.f, t_max = 1.f;
  std::tie(t_min, t_max) = computePathClippingFactors(t_min, t_max, simHitEntryPos.x(), simHitPath.x(), hitInfo.length[0]);
  std::tie(t_min, t_max) = computePathClippingFactors(t_min, t_max, simHitEntryPos.y(), simHitPath.y(), hitInfo.length[1]);

  // Apply clipping
  if (t_min != 0.f || t_max != 1.f) { // check if clipping is even necessary, for performance
    hitInfo.debugFlag = true;
    if (0. <= t_min && t_min <= t_max && t_max <= 1.) {
      verbose() << "   - Clipping simHitPath to sensor edges with factors t_min = " << t_min << ", t_max = " << t_max << ". PathLength changed to " << (t_max - t_min) * simHitPath.r() << " mm from " << simHitPath.r() << " mm" << endmsg;
      
      simHitEntryPos = simHitEntryPos + t_min * simHitPath;
      simHitPath = (t_max - t_min) * simHitPath;

      // if pathLength given by Geant4 is more than X% shorter than the length we calculate, shorten our calculated path accordingly.
      if (simHitPath.r() > m_pathLengthShorteningFactorGeant4 * hitInfo.simPathLength) {
        verbose() << "   - Shortening simHitPath from " << simHitPath.r() << " mm to Geant4 pathLength of " << hitInfo.simPathLength << " mm, because it's length is more than " << m_pathLengthShorteningFactorGeant4 << " of the Geant4 path length." << endmsg;
        
        /* make sure the path stays as centered around the simHitPos as possible (regarding edges)  
        * find out where the simHitPos lies on the current path:
        * project (simHitPos - simHitEntryPos) onto simHitPath -> gives distance from entryPos to simHitPos along the path (or to the point on path closest to simHitPos)
        * dotProduct (simHitPos-simHitEntryPos, simHitPath) = |simHitPos - simHitEntryPos| * |simHitPath| * cos(angle between them) */
        const float t_simHit = ( (hitInfo.simLocalPos - simHitEntryPos).dot(simHitPath) ) / ( simHitPath.r() * simHitPath.r() ); 
        
        // clamp new path (centered at t_center) to [0,1]
        const float t_length_half = hitInfo.simPathLength / simHitPath.r() / 2.f; // half-length of new path, in terms of t [0,1] on old path
        const float t_center = std::max(t_length_half, std::min(t_simHit, 1.f - t_length_half)); // center of new path, clamped to [t_length_half, 1 - t_length_half] while not exceeding [0,1]

        t_min = t_center - t_length_half;
        t_max = t_center + t_length_half;

        simHitEntryPos = simHitEntryPos + t_min * simHitPath;
        simHitPath = (t_max - t_min) * simHitPath;
      }
    } 
    else {
      warning() << "ConstructSimHitPath(): Cannot clip simHitPath to sensor edges. The path lies completely outside the sensor volume. Clipping t_min = " << t_min << ", t_max = " << t_max << "." << endmsg;
      verbose() << "   - before clipping: EntryPos: (" << simHitEntryPos.x() << " mm, " << simHitEntryPos.y() << " mm, " << simHitEntryPos.z() << " mm), exitPos: (" << simHitPath.x() << " mm, " << simHitPath.y() << " mm, " << simHitPath.z() << " mm)" << endmsg;
      // hitInfo.debugFlag = true;
      // return std::make_tuple(dd4hep::rec::Vector3D(0.0, 0.0, 0.0), dd4hep::rec::Vector3D(0.0, 0.0, 0.0));
    }
  }

  verbose() << "   - Calculated SimHitPath with length " << simHitPath.r() << " mm" << " (the Geant4 path length is " << hitInfo.simPathLength << " mm)" << endmsg;
  return std::make_tuple(simHitEntryPos, simHitPath);
}

void VTXdigi_Allpix2::CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const dd4hep::rec::Vector3D& position, const float charge) const {
  // overload to allow passing dd4hep::rec::Vector3D as position. ~ Jona 2025-09

  auto digiHit = digiHits.create();
  digiHit.setCellID(simHit.getCellID());
  digiHit.setEDep(charge / m_chargePerkeV); // convert e- to keV
  digiHit.setPosition(convertVector(position));
  // TODO: check if position is within sensor bounds & force it onto sensor simSurface ~ Jona 2025-09
  digiHit.setTime(simHit.getTime());
  
  auto digiHitLink = digiHitsLinks.create();
  digiHitLink.setFrom(digiHit);
  digiHitLink.setTo(simHit);
}


// -- Pixel and in-pixel binning magic (logic) --

int VTXdigi_Allpix2::computeBinIndex(float x, float binX0, float binWidth, int binN) const {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  if (binN <= 0) throw GaudiException("computeBinIndex: binN must be positive", "VTXdigi_Allpix2::computeBinIndex()", StatusCode::FAILURE);
  if (binWidth <= 0.0) throw GaudiException("computeBinIndex: binWidth must be positive", "VTXdigi_Allpix2::computeBinIndex()", StatusCode::FAILURE);

  float relativePos = (x - binX0) / binWidth;
  if (relativePos < 0.0f || relativePos > static_cast<float>(binN))
    return -1;
  return static_cast<int>(relativePos);
} // computeBinIndex()

std::tuple<int, int> VTXdigi_Allpix2::computePixelIndices(const dd4hep::rec::Vector3D& segmentPos, const int layerIndex, const float length_u, const float length_v) const {
  
  int i_u = computeBinIndex(
    segmentPos.x(),
    -0.5*length_u,
    m_pixelPitch_u.value().at(layerIndex),
    m_pixelCount_u.value().at(layerIndex));

  int i_v = computeBinIndex(
    segmentPos.y(),
    -0.5*length_v,
    m_pixelPitch_v.value().at(layerIndex),
    m_pixelCount_v.value().at(layerIndex));
  return std::make_tuple(i_u, i_v);
} // computePixelIndices()

std::tuple<int, int, int> VTXdigi_Allpix2::computeInPixelIndices(const dd4hep::rec::Vector3D& segmentPos, const int layerIndex, const float length_u, const float length_v) const {
  int j_u, j_v, j_w;

  // compute in-pixel position
  float shiftedPos_u = segmentPos.x() + 0.5 * length_u; // expected in [0, length_u]
  float pitch_u = m_pixelPitch_u.value().at(layerIndex);
  float inPixelPos_u = std::fmod(shiftedPos_u, pitch_u);
  if (inPixelPos_u < 0.0) inPixelPos_u += pitch_u; // ensure positive remainder

  j_u = computeBinIndex(inPixelPos_u, 0.0, pitch_u / m_inPixelBinCount_u.value(), m_inPixelBinCount_u.value());

  float shiftedPos_v = segmentPos.y() + 0.5 * length_v;
  float pitch_v = m_pixelPitch_v.value().at(layerIndex);
  float inPixelPos_v = std::fmod(shiftedPos_v, pitch_v);
  if (inPixelPos_v < 0.0) inPixelPos_v += pitch_v;

  j_v = computeBinIndex(inPixelPos_v, 0.0, pitch_v / m_inPixelBinCount_v.value(), m_inPixelBinCount_v.value());

  // vertical (w) binning: shift to [0, thickness]
  float shiftedPos_w = segmentPos.z() + 0.5 * m_sensorThickness.value().at(layerIndex);
  j_w = computeBinIndex(shiftedPos_w, 0.0, m_sensorThickness.value().at(layerIndex) / m_inPixelBinCount_w.value(), m_inPixelBinCount_w.value());

  return std::make_tuple(j_u, j_v, j_w);
} // computeInPixelIndices()

std::tuple<int, int, int, int, int> VTXdigi_Allpix2::computeSegmentIndices(HitInfo& hitInfo, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const int segment) const {
  /* Get the segment position, calculate pixel indices and in-pixel indices.
   * 
   * If segment is outside sensor area, return (-1,-1,-1,-1,-1)
   */
 verbose() << "     - The pathLength is " << simHitPath.r()*1000 << " um" << endmsg;
  if (segment < 0 || segment >= hitInfo.nSegments) {
    error() << "computeSegmentIndices(): Invalid segment number " << hitInfo.nSegments << " or segment index " << segment << endmsg;
    throw std::runtime_error("VTXdigi_Allpix2::computeSegmentIndices(): Invalid segment number or segment index");
  }

  float pathFraction = (static_cast<float>(segment)+0.5) / hitInfo.nSegments;
  const dd4hep::rec::Vector3D segmentRelativePos = pathFraction * simHitPath; // in mm
  verbose() << "     - Processing segment " << segment << " out of " << hitInfo.nSegments << " at fraction " << pathFraction << " with length " << segmentRelativePos.r()*1000 << " um" << endmsg;
  verbose() << "     - The pathLength is " << simHitPath.r()*1000 << " um" << endmsg;


  const dd4hep::rec::Vector3D segmentPos = simHitEntryPos + segmentRelativePos; // in mm
  verbose() << "         - Local position (u,v,w): (" << segmentPos.x() << " mm, " << segmentPos.y() << " mm, " << segmentPos.z() << " mm)" << endmsg;
  int i_u, i_v; // pixel indices
  std::tie(i_u, i_v) = computePixelIndices(segmentPos, hitInfo.layerIndex, hitInfo.length[0], hitInfo.length[1]);
  verbose() << "         - Pixel indices (" << i_u << ", " << i_v << ")" << endmsg;
  if (i_u==-1 || i_v==-1) {
    warning() << "computeSegmentIndices(): Segment lies outside sensor area (in u or v). Dismissing." << endmsg;
    // hitInfo.debugFlag = true;
    return std::make_tuple(-1, -1, -1, -1, -1);
  }
    // segment is outside sensor area

  int j_u, j_v, j_w; // in-pixel-binning indices
  std::tie(j_u, j_v, j_w) = computeInPixelIndices(segmentPos, hitInfo.layerIndex, hitInfo.length[0], hitInfo.length[1]);
  verbose() << "         - In-pixel indices (" << j_u << ", " << j_v << ", " << j_w << ")" << endmsg;
  if (j_u==-1 || j_v==-1 || j_w==-1) {
    warning() << "computeSegmentIndices(): Segment lies inside sensor area (in u and v), but vertically outside sensor volume. Dismissing." << endmsg;
    // hitInfo.debugFlag = true;
    return std::make_tuple(-1, -1, -1, -1, -1);
  }

  return std::make_tuple(i_u, i_v, j_u, j_v, j_w);
}

dd4hep::rec::Vector3D VTXdigi_Allpix2::computePixelCenter_Local(const int i_u, const int i_v, const int layerIndex, const dd4hep::rec::ISurface& simSurface) const {
  // returns the position of the center of pixel i_u, i_v) in the global
  // u - short ARCADIA axis, corresponds to x in sensor local frame
  // v - long ARCADIA axis, corresponds to y in sensor local frame

  float length_u = simSurface.length_along_u() * 10; // convert to mm (works, checked 2025-10-17)
  float length_v = simSurface.length_along_v() * 10; // convert to mm 

  float posU = -0.5 * length_u + (i_u + 0.5) * m_pixelPitch_u.value().at(layerIndex); // in mm
  float posV = -0.5 * length_v + (i_v + 0.5) * m_pixelPitch_v.value().at(layerIndex); // in mm

  return dd4hep::rec::Vector3D(posU, posV, 0.); 
}

// -- Transformation between global frame and local sensor frame --

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

// -- Other Helper functions --

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

void VTXdigi_Allpix2::appendSimHitToCsv(const HitInfo& hitInfo, const long int eventNumber, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const int i_u, const int i_v) const {
  if (!m_debugCsvFile.is_open()) {
    error() << "appendSimHitToCsv(): DebugCsv file is not open." << endmsg;
    return;
  }

  // eventNumber,layerIndex,segmentCount,sensorThickness,pix_u,pix_v,
  // simHitPos_u,simHitPos_v,simHitPos_w,
  // simHitEntryPos_u,simHitEntryPos_v,simHitEntryPos_w,
  // simHitPath_u,simHitPath_v,simHitPath_w,
  // pathLengthGeant4,pathLength,chargeDeposition,debugFlag

  m_debugCsvFile << std::to_string(eventNumber) << "," << hitInfo.layerIndex << "," << hitInfo.nSegments << ","  << hitInfo.thickness << "," << i_u << "," << i_v << ",";

  m_debugCsvFile << hitInfo.simLocalPos.x() << "," << hitInfo.simLocalPos.y() << "," << hitInfo.simLocalPos.z() << ",";
  m_debugCsvFile << simHitEntryPos.x() << "," << simHitEntryPos.y() << "," << simHitEntryPos.z() << ",";
  m_debugCsvFile << simHitPath.x() << "," << simHitPath.y() << "," << simHitPath.z() << ",";
  m_debugCsvFile << hitInfo.simPathLength << "," << simHitPath.r() << "," << hitInfo.simCharge <<  "," << hitInfo.debugFlag << "\n";
  m_debugCsvFile.flush();
  verbose() << "Wrote simHit with event " << eventNumber << ", layerIndex " << hitInfo.layerIndex << ", segmentN " << hitInfo.nSegments << ", i_u " << i_u << ", i_v " << i_v << " to debug CSV file." << endmsg;
  return;
}