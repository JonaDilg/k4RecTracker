#pragma once

// #include "Gaudi/Algorithm.h"
// #include "GaudiKernel/IRndmGenSvc.h"
// #include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "Gaudi/Accumulators/RootHistogram.h" // added by Jona
#include "Gaudi/Accumulators/Histogram.h" // added by Jona

#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h" // added by Jona
#include "edm4hep/TrackerHitPlaneCollection.h" // added by Jona
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h" // added by Jona

// K4FWCORE
// #include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

#include "k4FWCore/Transformer.h" // added by Jona
#include "k4Interface/IUniqueIDGenSvc.h" // added by Jona

// DD4HEP
#include "DDRec/SurfaceManager.h"

#include "TRandom2.h"

#include <string> // added by Jona
#include <vector>
#include <cmath> // for std::fmod

// for debugging-csv
#include <iostream>
#include <fstream>  // for std::ofstream
#include <sstream>

#include "GAUDI_VERSION.h"

#if GAUDI_MAJOR_VERSION < 39 // don't really have a clue why we need this, assume it's important ~ Jona 2025-09
namespace Gaudi::Accumulators {
template <unsigned int ND, atomicity Atomicity = atomicity::full, typename Arithmetic = double>
using StaticRootHistogram =
    Gaudi::Accumulators::RootHistogramingCounterBase<ND, Atomicity, Arithmetic, naming::histogramString>;
}
#endif

/** @class VTXdigi_Allpix2
 *
 * Creates trackerHits from simHits. Produces clusters from simHits, outputs either the cluster centre or all hits in the cluster as digitized hits.
 *
 *  @author Jona Dilg, Armin Ilg
 *  @date   2025-09
 *
 */

struct VTXdigi_Allpix2 final 
  : k4FWCore::MultiTransformer 
    <std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> (const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

  VTXdigi_Allpix2(const std::string& name, ISvcLocator* svcLoc);
  
  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator() (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const override;

private: 

  /* -- Classes used internally -- */

  /** @brief Class to store relevant quantities about a simHit 
   * @note Content is defined in constructor, cannot be modified afterwards. The exceptions are nSegments and the debugFlag. */
  class HitInfo; 
  friend class HitInfo;

  /** @brief Struct to store the position information of a simHit, including its path through the sensor */
  struct HitPosition {
    dd4hep::rec::Vector3D entry;
    dd4hep::rec::Vector3D path;
    dd4hep::rec::Vector3D global;
    dd4hep::rec::Vector3D local;
  };

  /** @brief Struct to store the indices of a path segment */
  struct SegmentIndices {
    int i_u = -1, i_v = -1;
    int j_u = -1, j_v = -1, j_w = -1;

    inline bool operator==(const SegmentIndices& other) const {
      return (i_u == other.i_u) && (i_v == other.i_v) && (j_u == other.j_u) && (j_v == other.j_v) && (j_w == other.j_w);
    }
  };

  /** @brief Class to store the charge deposited in each pixel around the simHit. 
   * @note Matrix size starts from property "MaximumClusterSize" and expands dynamically if charge is shared to pixels outside of that range */
  class PixelChargeMatrix;
  
  /** @brief Class to store and access charge sharing kernels */
  class ChargeSharingKernels;

  /* ---- Initialization & finalization functions ---- */

  void initializeServicesAndGeometry();
  void checkGaudiProperties();
  void setupSensorParameters();
  void loadKernels();
  void setupDebugHistograms();
  void setupDebugCsvOutput();

  void printCountersSummary() const;

  /* ---- Core algorithm functions ---- */

  /** @brief Initial checks before processing an event.
   * @return true if setup is OK and event can be processed, false if event has no simHits 
   * @note throws if chargeSharingKernels or surfaceMap are invalid */
  bool CheckInitialSetup(const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const;

  /** @brief Apply layer cuts to a simHit: reject hits not on layers to be digitized, if Gaudi property "LayersToDigitize" is set.
   * @return true if the hit is accepted, false if dismissed. */
  bool CheckLayerCut(const edm4hep::SimTrackerHit& simHit) const;

  /** @brief Gather hit information and position from a simHit.
   * @return A tuple containing the hit information and position. */
  std::tuple<HitInfo, HitPosition> GatherHitInfoAndPosition(const edm4hep::SimTrackerHit& simHit, const edm4hep::EventHeaderCollection& headers) const;

  /** @brief Apply cuts to a simHit. 
   * @return true if the hit is accepted, false if dismissed. */
  bool CheckSimHitCuts(const HitInfo& hitInfo, const HitPosition& hitPos) const;

  /** @brief Loop over the path segments and share each segments deposited charge among its neighbors according to the charge sharing kernels.
   * @return A PixelChargeMatrix containing the charge deposited in each pixel around the simHit.
   * @note Noise has not yet been generated for the pixelChargeMatrix. */
  PixelChargeMatrix DepositAndCollectCharge(HitInfo& hitInfo, const HitPosition& hitPos) const;

  /** @brief Find pixels with charge above threshold in pixelChargeMatrix, create digiHits and fill digiHits and digiHitsLinks collections */
  void AnalyseSharedCharge(const HitInfo& hitInfo, const HitPosition& hitPos, const PixelChargeMatrix& pixelChargeMatrix, const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks) const;
  
  /** @brief Fill debug histograms.
   * @note Contains a very similar loop to AnalyseSharedCharge, but does not create digiHits. This optimises performance with debug histograms disabled.*/
  void FillGeneralDebugHistograms(HitInfo& hitInfo, const HitPosition& hitPos, const PixelChargeMatrix& pixelChargeMatrix) const;


  /* ---- Helper functions (called by core algorithm functions) ---- */

  /* -- Pixel and in-pixel binning magic (logic) -- */

  /** @brief Given a histogram definition (x0, binWidth, nBins) and a value x, compute the bin index i in which x falls.
   * @return Int, -1 if x is out of range.
   * @note Bins are 0-indexed (vs ROOT's 1-indexing) */
  int computeBinIndex(float x, float binX0, float binWidth, int binN) const;

  /** @brief Compute the pixel indices (i_u, i_v) for a given (local) position inside the sensor */
  std::tuple<int, int> computePixelIndices(const dd4hep::rec::Vector3D& pos, const int layerIndex, const float length_u, const float length_v) const;
  
  /** @brief Compute the in-pixel indices (j_u, j_v, j_w) for a given (local) position inside the pixel and layer index
   *  @note Assumption: each layer has only 1 type if sensor */
  std::tuple<int, int, int> computeInPixelIndices(const dd4hep::rec::Vector3D& pos, const int layerIndex, const float length_u, const float length_v) const;

  /** @brief Compute the local position of a pixel center (u,v,0) */
  dd4hep::rec::Vector3D computePixelCenter_Local(const int i_u, const int i_v, const int layerIndex, const dd4hep::rec::ISurface& simSurface) const;
  
  /** @brief Find the in-pixel indices (j_u, j_v, j_w) for a given (local) position inside the pixel and layer index
   *  @note Assumption: each layer has only 1 type if sensor */
  SegmentIndices computeSegmentIndices(HitInfo& hitInfo, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const int segment) const;


  /* -- Transformation between global frame and local sensor frame -- */

  /** @brief Computes the transformation matrix to the local sensor frame */
  TGeoHMatrix computeTransformationMatrix(const edm4hep::SimTrackerHit& simHit) const;
  /** @copydoc computeTransformationMatrix(const edm4hep::SimTrackerHit&) */
  TGeoHMatrix computeTransformationMatrix(const dd4hep::DDSegmentation::CellID& cellID) const;

  /** @brief Adjust the transformation matrix to the local sensor frame.
   * Such that the coordinate system is x-u, y-v, z-n and right-handed. Uses user-input "LocalNormalVectorDir"
   */
  void setProperDirectFrame(TGeoHMatrix& sensorTransformMatrix) const;
  
  /** @brief Transform a global position to the local sensor frame of a hit, given its cellID */
  dd4hep::rec::Vector3D transformGlobalToLocal(const dd4hep::rec::Vector3D& globalPos, const dd4hep::DDSegmentation::CellID& cellID) const;

  /** @brief Transform a local sensor position to the global frame, given the sensors cellID */
  dd4hep::rec::Vector3D transformLocalToGlobal(const dd4hep::rec::Vector3D& localPos, const dd4hep::DDSegmentation::CellID& cellID) const;


  /* -- Other helper functions -- */

  /** @brief Simply convert EDM4HEP vector to DD4HEP vector, and vice versa */
  dd4hep::rec::Vector3D convertVector(edm4hep::Vector3d vec) const;
  /** @copydoc convertVector(edm4hep::Vector3d) */
  edm4hep::Vector3d convertVector(dd4hep::rec::Vector3D vec) const;

  /** @brief Calculate the entry point and path vector of a simHit (in local sensor frame).
   * @return A tuple containing (0) the entryPoint into the sensor and (1) the path vector through the sensor */
  std::tuple<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D> constructSimHitPath(HitInfo& hitInfo, HitPosition& hitPos, const edm4hep::SimTrackerHit& simHit) const;

  /** @brief Calculate the clipping factors for the path of a simHit in local sensor frame. Called in constructSimHitPath() */
  std::tuple<float, float> computePathClippingFactors(float t_min, float t_max, const float entryPos_ax, const float pathLength_ax, const float sensorLength_ax) const;

  void collectSegmentCharge(HitInfo& hitInfo, PixelChargeMatrix& pixelChargeMatrix, const SegmentIndices& segment, const float segmentCharge) const;

  /** @brief Create a digitized hit*/
  void createDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const dd4hep::rec::Vector3D& position, const float charge) const;

  /** @brief Write simHit & path information to the debugging CSV file. Definitely not thread-safe. */
  void appendSimHitToCsv(const HitInfo& hitInfo, const HitPosition& hitPos, const int i_u, const int i_v) const;

  void fillDebugHistograms_segmentLoop(const HitInfo& hitInfo, const SegmentIndices& segment, int i_m, int i_n, const float sharedCharge) const;
  void fillDebugHistograms_targetPixelLoop(const HitInfo& hitInfo, const HitPosition& hitPos, int i_u, int i_v, float pixelChargeMeasured) const;
  

  /* -- Properties -- */

  const std::string m_undefinedString = "UNDEFINED";

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", m_undefinedString, "Name of the subdetector"};
  Gaudi::Property<std::string> m_geometryServiceName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"}; // what is this for?
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  
  Gaudi::Property<bool> m_cutPathOutsideSensor{this, "CutPathOutsideSensor", false, "Whether to cut simHits where the simHitPath is not on or inside the sensor."};
  Gaudi::Property<float> m_cutDepositedCharge{this, "CutDepositedCharge", 0.0, "Minimum charge (e-) of SimTrackerHit to be digitized"};
  
  Gaudi::Property<std::vector<int>> m_layersToDigitize{this, "LayersToDigitize", {}, "Which layers to digitize (0-indexed). If empty, all layers are digitized."};

  Gaudi::Property<std::vector<int>> m_pixelCount_u{this, "PixelCount_u", {}, "Number of pixels in direction of u; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<int>> m_pixelCount_v{this, "PixelCount_v", {}, "Number of pixels in direction of v; either one per layer or one for all layers"};
  
  Gaudi::Property<float> m_targetPathSegmentLength{this, "TargetPathSegmentLength", 0.002, "Length of the path segments, that the simHits path through a sensor is divided into. In mm. Defines the precision of the charge deposition along the path."};
  Gaudi::Property<float> m_pathLengthShorteningFactorGeant4{this, "PathLengthShorteningFactorGeant4", 1.05, "Relative path length (to Geant4 path length), above which the path is shortened to the Geant4 length. (Geant4 length includes multiple scattering and curling in B-field, which are lost in our linear approximation of the path)."};
  
  Gaudi::Property<std::vector<int>> m_maxClusterSize{this, "MaximumClusterSize", {15,15}, "Default maximum cluster size [max in u, max in v] in terms of pixels. The area is expanded dynamically if charge is shared to pixels outside of this range, but this is computationally expensive."};
  Gaudi::Property<float> m_pixelThreshold{this, "PixelThreshold", 100, "Threshold in electrons for a pixel to fire (1 eh-pair = 3.65 eV)"};
  Gaudi::Property<float> m_electronicNoise{this, "PixelElectronicNoise", 20, "Electronic noise in electrons (1 eh-pair = 3.65 eV). Defines the width of the Gaussian noise added to each pixel."};
  Gaudi::Property<float> m_timeSmearFactor{this, "PixelTimeSmear", 0., "Gaussian width for the time smearing applied on the pixel time. Applied for each digiHit individually."};
  
  Gaudi::Property<std::string> m_localNormalVectorDir{this, "LocalNormalVectorDir", "", "Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo). If defined correctly, the local frame is transformed such that z is orthogonal to the sensor plane."};

  /* Kernel import properties */
  Gaudi::Property<std::vector<float>> m_globalKernel{this, "GlobalKernel", {}, "Flat vector containing the global charge sharing kernel in row-major order (ie. row-by-row), starting on top left. Length must be KernelSize*KernelSize"};
  Gaudi::Property<std::string> m_kernelFileName{this, "KernelFileName", "", "Name of the file supplying the charge sharing kernel."};

  /* Debugging */
  Gaudi::Property<bool> m_debugHistograms{this, "DebugHistograms", false, "Whether to create and fill debug histograms. Not recommended for multithreading, might lead to crashes."};
  Gaudi::Property<std::string> m_debugCsvFileName{this, "DebugCsvFileName", "", "Name of a CSV file to output detailed per-hit debug information. If empty, no debug output is created."};

  /* -- Services -- */
  
  SmartIF<IGeoSvc> m_geometryService;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_cellIDdecoder; // Decoder for the cellID

  dd4hep::VolumeManager m_volumeManager; // volume manager to get the physical cell sensitive volume
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  const dd4hep::rec::SurfaceMap* m_simSurfaceMap;

  /* -- Constants -- */

  const float m_numericLimit_float = 0.00001f; // limit for floating point comparisons
  const double m_numericLimit_double = 0.000000001; // limit for floating point
  const float m_chargePerkeV = 273.97f; // number of electron-hole pairs created per keV of deposited energy in silicon. eh-pair ~ 3.65 eV
  
  /* -- Detector & sensor parameters -- */

  int m_layerCount;
  std::unordered_map<int, int> m_layerToIndex; // layer number (from cellID & m_layersToDigitize) to internal index (0...N-1, where N is the number of layers to digitize)

  std::vector<double> m_sensorThickness;
  std::vector<double> m_pixelPitch_u;
  std::vector<double> m_pixelPitch_v;

  int m_inPixelBinCount[3];
  int m_kernelSize;
  std::unique_ptr<ChargeSharingKernels> m_chargeSharingKernels; // the charge sharing kernel

  /* -- Counters -- */

  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRead{this, "Events read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRejected_noSimHits{this, "Events rejected (no simHits)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsAccepted{this, "Events accepted"};
  
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRead{this, "SimTrackerHits read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRejected_LayerNotToBeDigitized{this, "SimTrackerHits rejected (layer not to be digitized)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRejected_ChargeCut{this, "SimTrackerHits rejected (charge cut)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRejected_OutsideSensor{this, "SimTrackerHits rejected (outside sensor)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsAccepted{this, "SimTrackerHits accepted"};
  mutable Gaudi::Accumulators::Counter<> m_counter_acceptedButNoSegmentsInSensor{this, "SimTrackerHits rejected (no segments in sensor)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_digiHitsCreated{this, "DigiHits created"};

  /* -- Debugging -- */

  bool m_debugCsv = false; // create & write debug CSV file?
  mutable std::ofstream m_debugCsvFile; // debug output file. Definitely not thread-safe.

  enum { 
    hist_simHitE, 
    hist_simHitCharge,
    hist_clusterSize_raw, 
    hist_clusterSize_measured,
    hist_EntryPointX, 
    hist_EntryPointY, 
    hist_EntryPointZ, 
    hist_DisplacementU, 
    hist_DisplacementV, 
    hist_DisplacementR, 
    hist_pathLength, 
    hist_pathLengthGeant4,
    hist_chargeCollectionEfficiency_raw, 
    hist_chargeCollectionEfficiency,
    hist_pixelChargeMatrix_size_u,
    hist_pixelChargeMatrix_size_v,
    histArrayLen}; // histogram indices. histArrayLen must be last
  std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<1>>, histArrayLen> m_histograms; 

  enum { 
    hist2d_hitMap_simHits,
    hist2d_hitMap_digiHits, 
    hist2d_pathLength_vs_simHit_v,
    hist2d_pixelChargeMatrixSize,
    hist2d_pathAngleToSensorNormal,
    hist2dArrayLen }; // 2D histogram indices. hist2dArrayLen must be last
  std::vector<std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<2>>, hist2dArrayLen>> m_histograms2d;

  enum {
    histWeighted2d_averageCluster,
    histWeighted2d_chargeOriginU,
    histWeighted2d_chargeOriginV,
    histWeighted2dArrayLen };
  std::vector<std::array<std::unique_ptr<Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, double>>, histWeighted2dArrayLen>> m_histWeighted2d;

  /* TODO: implement having a set of kernels per layer (array of unique_ptr ?) */
};

class VTXdigi_Allpix2::HitInfo {
  bool m_debugFlag = false; // set to true for hits that should be logged in the debug CSV output

  long int m_eventNumber;
  dd4hep::DDSegmentation::CellID m_cellID;
  dd4hep::rec::ISurface* m_simSurface;

  int m_layerIndex;
  float m_charge;
  float m_simPathLength;

  float m_length[2]; // sensor length in u and v direction [mm]
  float m_thickness; // sensor thickness [mm]
  float m_pixPitch[2]; // pixel pitch in u and v direction [mm]
  int m_pixCount[2]; // number of pixels in u and v direction

  int m_nSegments;

  public:
    HitInfo() = default;

    HitInfo(const VTXdigi_Allpix2& vtxdigi_AP2, const edm4hep::SimTrackerHit& simHit, const edm4hep::EventHeaderCollection& headers) {
      m_eventNumber = headers.at(0).getEventNumber();
      m_cellID = simHit.getCellID(); // TODO: apply 64-bit length mask as done in DDPlanarDigi.cpp? ~ Jona 2025-09

      const auto itSimSurface = vtxdigi_AP2.m_simSurfaceMap->find(m_cellID);
      if (itSimSurface == vtxdigi_AP2.m_simSurfaceMap->end())
        throw std::runtime_error("VTXdigi_Allpix2::HitInfo constructor (called from VTXdigi_Allpix2::operator()): No simSurface found for cellID " + std::to_string(m_cellID) + ". Did initialize() succeed?");

      m_simSurface = itSimSurface->second;
      if (!m_simSurface)
        throw std::runtime_error("VTXdigi_Allpix2::HitInfo constructor (called from VTXdigi_Allpix2::operator()): SimSurface pointer for cellID " + std::to_string(m_cellID) + " is null. Did initialize() succeed?");

      /* Use layer index instead of layer to avoid segfaults, if layers are not numbered consecutively
       * This will throw if layer not in map (ie. layers not-to-be digitized have not beed dismissed yet).*/
      m_layerIndex = vtxdigi_AP2.m_layerToIndex.at(vtxdigi_AP2.m_cellIDdecoder->get(m_cellID, "layer"));

      m_charge = simHit.getEDep() * (dd4hep::GeV / dd4hep::keV) * vtxdigi_AP2.m_chargePerkeV; // in electrons
      m_simPathLength = simHit.getPathLength(); // in mm

      // sensor dimensions in local frame, in mm
      m_length[0] = m_simSurface->length_along_u() * 10; // convert to mm
      m_length[1] = m_simSurface->length_along_v() * 10; // convert to mm

      m_thickness = vtxdigi_AP2.m_sensorThickness.at(m_layerIndex);

      m_pixPitch[0] = vtxdigi_AP2.m_pixelPitch_u.at(m_layerIndex);
      m_pixPitch[1] = vtxdigi_AP2.m_pixelPitch_v.at(m_layerIndex);

      m_pixCount[0] = vtxdigi_AP2.m_pixelCount_u.value().at(m_layerIndex);
      m_pixCount[1] = vtxdigi_AP2.m_pixelCount_v.value().at(m_layerIndex);

      /* sanity check: pixel pitch * pixel count must match sensor size from geometry. 
      * Do this after applying cuts, as non-digitized layers would lead to out-of-bounds access in m_pixelCount_u etc.
      * Doing this for every simHit is not optimal (because pixPitch and pixCount per layer are constant across all events)
      * but accessing the surface-based sensor size from initialize() is not easy to implement (without access to a cellID that is known to lie on that layer). */
      if ( abs(m_pixPitch[0] * m_pixCount[0] - m_length[0]) > vtxdigi_AP2.m_numericLimit_float
        || abs(m_pixPitch[1] * m_pixCount[1] - m_length[1]) > vtxdigi_AP2.m_numericLimit_float) {
        throw GaudiException("Gaudi properties (pixelPitch * pixelCount) do not match sensor size from detector geometry in layer " + std::to_string(vtxdigi_AP2.m_cellIDdecoder->get(m_cellID, "layer")) + " in operator().", "VTXdigi_Allpix2::GatherHitInfo()", StatusCode::FAILURE);
      }
    }


    inline void setDebugFlag() { m_debugFlag = true; }
    inline bool debugFlag() const { return m_debugFlag; }

    inline long int eventNumber() const { return m_eventNumber; }
    inline dd4hep::DDSegmentation::CellID cellID() const { return m_cellID; }
    inline dd4hep::rec::ISurface* simSurface() const { return m_simSurface; }

    inline int layerIndex() const { return m_layerIndex; }
    inline float charge() const { return m_charge; }
    inline float simPathLength() const { return m_simPathLength; }

    inline float length(int axis) const { return m_length[axis]; } // axis: 0 = u, 1 = v
    inline float thickness() const { return m_thickness; }
    inline float pixPitch(int axis) const { return m_pixPitch[axis]; } // axis: 0 = u, 1 = v
    inline int pixCount(int axis) const { return m_pixCount[axis]; } // axis: 0 = u, 1 = v  

    inline void setNSegments(int n) { m_nSegments = n; }
    inline int nSegments() const { return m_nSegments; }
  }; // class HitInfo

class VTXdigi_Allpix2::PixelChargeMatrix {
  /* Stores the charge deposited in a (size_u x size_v) pixel matrix around a given origin pixel.
    * In case charge is added outside the matrix bounds, the matrix range is expanded in that direction.
    * The size of the matrix is defined via the Gaudi property MaximumClusterSize.
  */

  std::vector<float> m_pixelCharge;
  std::vector<float> m_pixelChargeNoise;

  const int m_overExpansionStep = 0; // number of extra pixels to expand the matrix by, when a charge is added outside the current bounds
  
  int m_range_u[2], m_range_v[2]; // Inclusive matrix range. -> size = range[1] - range[0] + 1
  /* range CAN theoretically extend into negative values, this ensures graceful handling of hits outside inditial bounds. These might be discarded later, if outside of sensor. */
  int m_origin[2]; // origin pixel indices

  public:

    PixelChargeMatrix(int i_origin_u, int i_origin_v, int size_u, int size_v) 
      : m_origin{ i_origin_u, i_origin_v } {
        // At first, the range is centered around the origin
      m_range_u[0] = m_origin[0] - (size_u - 1) / 2;
      m_range_u[1] = m_origin[0] + (size_u - 1) / 2;
      m_range_v[0] = m_origin[1] - (size_v - 1) / 2;
      m_range_v[1] = m_origin[1] + (size_v - 1) / 2;

      m_pixelCharge.resize(size_u * size_v, 0.f);
    }

    void Reset() {
      std::fill(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
    }
    
    inline int GetOriginU() const { return m_origin[0]; }
    inline int GetOriginV() const { return m_origin[1]; }
    inline int GetRangeMin_u() const { return m_range_u[0]; }
    inline int GetRangeMax_u() const { return m_range_u[1]; }
    inline int GetRangeMin_v() const { return m_range_v[0]; }
    inline int GetRangeMax_v() const { return m_range_v[1]; }
    inline int GetSize_u() const { return m_range_u[1] - m_range_u[0] + 1; }
    inline int GetSize_v() const { return m_range_v[1] - m_range_v[0] + 1; }

    inline std::tuple<int, int> GetSize() const {
      return std::make_tuple(GetSize_u(), GetSize_v());
    }
    inline std::tuple<int, int> GetOrigin() const {
      return std::make_tuple(m_origin[0], m_origin[1]);
    }

    float GetTotalRawCharge() const {
      return std::accumulate(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
    }

    float GetRawCharge(int i_u, int i_v) const {
      if (_OutOfBounds(i_u, i_v)) {
        throw std::runtime_error("PixelChargeMatrix::GetCharge: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") out of range");
      }
      return m_pixelCharge[_FindIndex(i_u, i_v)];
    }

    void FillRawCharge(int i_u, int i_v, float charge) {
      if (_OutOfBounds(i_u, i_v)) {
        // this might occur due to the limited size of the matrix or delta-electrons. 
        _ExpandMatrix(i_u, i_v);
      }
      if (_OutOfBounds(i_u, i_v)) {
        throw std::runtime_error("PixelChargeMatrix::FillCharge: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") still out of range ( " + std::to_string(m_range_u[0]) + ", " + std::to_string(m_range_u[1]) + ", " + std::to_string(m_range_v[0]) + ", " + std::to_string(m_range_v[1]) + ") after expansion");
      }
      m_pixelCharge[_FindIndex(i_u, i_v)] += charge;
    }

    // functins related to noise

    void GenerateNoise(TRandom2& rndEngine, float sigma) {
      /* Generate noise values for each pixel in the matrix.
       * @param sigma Standard deviation of the Gaussian distribution.
       * @param randGen Random number generator to use.
      */
      m_pixelChargeNoise.resize(GetSize_u() * GetSize_v(), 0.f);
      for (int i = 0; i < GetSize_u() * GetSize_v(); ++i) {
        m_pixelChargeNoise[i] = rndEngine.Gaus(0.f, sigma);
      }
    }

    float GetNoise(int i_u, int i_v) const {
      if (_OutOfBounds(i_u, i_v)) {
        throw std::runtime_error("PixelChargeMatrix::GetNoise: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") out of range");
      }
      if (m_pixelChargeNoise.size() != m_pixelCharge.size()) {
        throw std::runtime_error("PixelChargeMatrix::GetNoise: noise has not been generated yet.");
      }
      return m_pixelChargeNoise[_FindIndex(i_u, i_v)];
    }

    float GetMeasuredCharge(int i_u, int i_v) const {
      if (_OutOfBounds(i_u, i_v)) {
        throw std::runtime_error("PixelChargeMatrix::GetCharge: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") out of range");
      }
      if (m_pixelChargeNoise.size() != m_pixelCharge.size()) {
        throw std::runtime_error("PixelChargeMatrix::GetMeasuredCharge: noise has not been generated yet.");
      }
      return ( 
        m_pixelCharge[_FindIndex(i_u, i_v)] +
        m_pixelChargeNoise[_FindIndex(i_u, i_v)]);
    }

    // float GetTotalMeasuredCharge() const {
    //   if (m_pixelChargeNoise.size() != m_pixelCharge.size()) {
    //     throw std::runtime_error("PixelChargeMatrix::TotalMeasuredCharge: noise has not been generated yet.");
    //   }

    //   return (
    //     std::accumulate(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f) + 
    //     std::accumulate(m_pixelChargeNoise.begin(), m_pixelChargeNoise.end(), 0.f));
    // }

    float GetTotalMeasuredCharge(float threshold) const {
      if (m_pixelChargeNoise.size() != m_pixelCharge.size())
        throw std::runtime_error("PixelChargeMatrix::TotalMeasuredCharge: noise has not been generated yet.");
      
      float totalCharge = 0.f;
      for (int i_u = m_range_u[0]; i_u <= m_range_u[1]; ++i_u) {
        for (int i_v = m_range_v[0]; i_v <= m_range_v[1]; ++i_v) {
          float measuredCharge = GetMeasuredCharge(i_u, i_v);
          if (measuredCharge >= threshold)
            totalCharge += measuredCharge;
        }
      }
      return totalCharge;
    }

  private:
    
    inline int _FindIndex(int i_u, int i_v) const {
      int i_u_rel = i_u - m_range_u[0];
      int i_v_rel = i_v - m_range_v[0];
      return i_u_rel + i_v_rel * GetSize_u();
    }

    inline bool _OutOfBounds(int i_u, int i_v) const {
      return (
        i_u < m_range_u[0]
        || i_u > m_range_u[1]
        || i_v < m_range_v[0]
        || i_v > m_range_v[1]);
    }

    void _ExpandMatrix(int i_u, int i_v) {
      /* Expand the matrix to include (i_u, i_v) and an excess of m_expansionStep pixels.
       */

      int rangeNew_u[2] = { m_range_u[0], m_range_u[1] }; 
      int rangeNew_v[2] = { m_range_v[0], m_range_v[1] };

      // TODO: optimise the expansion logic below, to expand to two directions at once if hit is close to corner ~ Jona 2025-10

      // expand in u direction?
      if (i_u < m_range_u[0]) {
        rangeNew_u[0] = i_u - m_overExpansionStep; // expand to include i_u plus some excess
      } else if (i_u > m_range_u[1]) {
        rangeNew_u[1] = i_u + m_overExpansionStep;
      }
      // expand in v direction?
      if (i_v < m_range_v[0]) {
        rangeNew_v[0] = i_v - m_overExpansionStep;
      } else if (i_v > m_range_v[1]) {
        rangeNew_v[1] = i_v + m_overExpansionStep;
      }

      const int newSizeU = rangeNew_u[1] - rangeNew_u[0] + 1;
      const int newSizeV = rangeNew_v[1] - rangeNew_v[0] + 1;

      std::vector<float> newPixelCharge(newSizeU * newSizeV, 0.f);

      // copy old charges into new array, row by row
      
      for (int row = 0; row < GetSize_v(); ++row) {
        std::copy(
          m_pixelCharge.begin() + row * GetSize_u(),
          m_pixelCharge.begin() + (row + 1) * GetSize_u(),
          newPixelCharge.begin() + (row + (m_range_v[0] - rangeNew_v[0])) * newSizeU + (m_range_u[0] - rangeNew_u[0])
        );
      }

      m_pixelCharge.swap(newPixelCharge);
      std::copy(rangeNew_u, rangeNew_u + 2, m_range_u);
      std::copy(rangeNew_v, rangeNew_v + 2, m_range_v);
    }
  }; // class PixelChargeMatrix

class VTXdigi_Allpix2::ChargeSharingKernels {
    int m_binCountU, m_binCountV, m_binCountW;
    int m_kernelSize;
    /** Vector of kernels, one per in-pixel bin
     *  Kernel-indexing is col-major (ie. i = col * size + row) */
    std::vector<std::vector<float>> m_kernels; 

  public:

    ChargeSharingKernels(int& binCountU, int& binCountV, int& binCountW, int& kernelSize) : m_binCountU(binCountU), m_binCountV(binCountV), m_binCountW(binCountW), m_kernelSize(kernelSize) {
      if (kernelSize < 3 || kernelSize % 2 == 0)
        throw std::runtime_error("ChargeSharingKernel: Kernel size must be an odd integer >= 3, but is " + std::to_string(kernelSize) + ".");

      m_kernels.resize(binCountU * binCountV * binCountW);
      for (auto& kernel : m_kernels) {
        kernel.resize(kernelSize*kernelSize, 0.f);
      }
    }

    /** Set the charge sharing kernel for a specific in-pixel bin 
     * @param weights A flat vector containing the kernel values in row-major order (ie. row-by-row) (length must be kernelSize*kernelSize)
     */
    void SetKernel(const int j_u, const int j_v, const int j_w, const std::vector<float>& weights) {
      if (static_cast<int>(weights.size()) != m_kernelSize*m_kernelSize)
        throw std::runtime_error("ChargeSharingKernel::SetKernel: weights size (" + std::to_string(weights.size()) + ") does not match kernel size (" + std::to_string(m_kernelSize*m_kernelSize) + ")");

      // check if Kernel is normalized to 1
      float sum = 0.f;
      for (int row = 0; row < m_kernelSize; ++row) {
        for (int col = 0; col < m_kernelSize; ++col) {
          sum += weights.at(row*m_kernelSize + col);
        }
      }
      if (std::abs(sum) > 1.000001f)
        throw GaudiException("Supplied ChargeSharingKernel weight sum needs to be <= 1, but is " + std::to_string(sum) + ", .", "VTXdigi_Allpix2::ChargeSharingKernels::SetKernel()", StatusCode::FAILURE);

      const int index = _FindIndex(j_u, j_v, j_w);

      for (int row = 0; row < m_kernelSize; ++row) {
        for (int col = 0; col < m_kernelSize; ++col) {
          // weights are given in row-major order, starting at top left. 
          // We store kernels in col-major order, starting at bottom left (lowest bin index)
          // m_kernels.at(index).at(col * m_kernelSize + row) = weights.at(row*m_kernelSize + col);
          m_kernels.at(index).at(col * m_kernelSize + row) = weights.at((m_kernelSize-1-row)*m_kernelSize + col);
        }
      }
    }

    void SetAllKernels(const std::vector<float>& weights) {
      for (int j_u = 0; j_u < m_binCountU; ++j_u) {
        for (int j_v = 0; j_v < m_binCountV; ++j_v) {
          for (int j_w = 0; j_w < m_binCountW; ++j_w) {
            SetKernel(j_u, j_v, j_w, weights);
          }
        }
      }
    }

    /** Access kernel as const reference */
    const std::vector<float>& GetKernel(int j_u, int j_v, int j_w) const {
      if (j_u < 0 || j_u >= m_binCountU)
        throw std::runtime_error("ChargeSharingKernel::GetKernel: j_u (= " + std::to_string(j_u) + ") out of range");
      if (j_v < 0 || j_v >= m_binCountV)
        throw std::runtime_error("ChargeSharingKernel::GetKernel: j_v (= " + std::to_string(j_v) + ") out of range");
      if (j_w < 0 || j_w >= m_binCountW)
        throw std::runtime_error("ChargeSharingKernel::GetKernel: j_w (= " + std::to_string(j_w) + ") out of range");
      return m_kernels[_FindIndex(j_u, j_v, j_w)];
    }

    /** Access a specific entry of a kernel
     * @param j_u, j_v, j_w In-pixel indices
     * @param j_col, j_row Column and row of the kernel entry to access (go from -(kernelsize-1)/2 to +(kernelsize-1)/2)
    */
    float GetWeight(int j_u, int j_v, int j_w, int j_col, int j_row) const {
      if (abs(j_col) > (m_kernelSize-1)/2)
        throw std::runtime_error("GetKernelEntry(): j_col (=" + std::to_string(j_col) + ") out of range of kernel size.");
      if (abs(j_row) > (m_kernelSize-1)/2)
        throw std::runtime_error("GetKernelEntry(): j_row (=" + std::to_string(j_row) + ") out of range of kernel size.");
      const auto& kernel = GetKernel(j_u, j_v, j_w);
      return kernel[(j_col + (m_kernelSize-1)/2) * m_kernelSize + (j_row + (m_kernelSize-1)/2)];
    }
    float GetWeight(const VTXdigi_Allpix2::SegmentIndices& segment, int i_col, int i_row) const {
      return GetWeight(segment.j_u, segment.j_v, segment.j_w, i_col, i_row);
    }

    inline int GetSize() const {
      return m_kernelSize;
    }
    inline int GetBinCountU() const {
      return m_binCountU;
    }
    inline int GetBinCountV() const {
      return m_binCountV;
    }
    inline int GetBinCountW() const {
      return m_binCountW;
    }

  private:

    int _FindIndex (int j_u, int j_v, int j_w) const {
      if (j_u < 0 || j_u >= m_binCountU)
        throw std::runtime_error("ChargeSharingKernel::_FindIndex: j_u (= " + std::to_string(j_u) + ") out of range");
      if (j_v < 0 || j_v >= m_binCountV)
        throw std::runtime_error("ChargeSharingKernel::_FindIndex: j_v (= " + std::to_string(j_v) + ") out of range");
      if (j_w < 0 || j_w >= m_binCountW)
        throw std::runtime_error("ChargeSharingKernel::_FindIndex: j_w (= " + std::to_string(j_w) + ") out of range");

      return j_u + m_binCountU * (j_v + m_binCountV * j_w); 
    }
  }; // class ChargeSharingKernel
