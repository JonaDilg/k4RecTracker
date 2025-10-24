#pragma once

// GAUDI
// #include "Gaudi/Algorithm.h"
// #include "GaudiKernel/IRndmGenSvc.h"
// #include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "Gaudi/Accumulators/RootHistogram.h" // added by Jona
#include "Gaudi/Accumulators/Histogram.h" // added by Jona

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"
// #include "edm4hep/TrackerHit3DCollection.h"
// #include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include "edm4hep/EventHeaderCollection.h" // added by Jona
#include "edm4hep/TrackerHitPlaneCollection.h" // added by Jona
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h" // added by Jona

// K4FWCORE
// #include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

#include "k4FWCore/Transformer.h" // added by Jona
#include "k4Interface/IUniqueIDGenSvc.h" // added by Jona

// DD4HEP
// #include "DD4hep/Detector.h" // for dd4hep::VolumeManager
#include "DDRec/SurfaceManager.h"
// #include "DDRec/Vector3D.h"

//#include "DDSegmentation/BitFieldCoder.h"


#include "TRandom2.h" // added by Jona
#include "TMatrixD.h" // added by Jona - for storing the kernel

#include <string> // added by Jona
#include <vector>
#include <cmath> // for std::fmod

// for debugging-csv
#include <iostream>
#include <fstream>  // for std::ofstream

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
  struct HitInfo {
    bool debugFlag = false;
    dd4hep::DDSegmentation::CellID cellID;

    dd4hep::rec::ISurface* simSurface;

    int layerIndex;
    float length[2]; // sensor length in u and v direction [mm]
    float thickness; // sensor thickness [mm]
    int pixCount[2]; // number of pixels in u and v direction
    float pixPitch[2]; // pixel pitch in u and v direction [mm]

    float simCharge;
    float simPathLength;
    dd4hep::rec::Vector3D simLocalPos;

    float charge;
    float pathLength;
    
    int nSegments;
    float segmentCharge;
  };
  
  // -- Core algorithm functions -- 
  
  /** Check that the pixel pitch and count match the sensor size in the geometry
   */
  void InitialSensorSizeCheck(const edm4hep::SimTrackerHit& simHit) const;
  // TODO: remove this. 

  /** @brief Apply cuts to a simHit. 
   * @return true if the hit is accepted, false if dismissed. 
   */
  bool ApplySimHitCuts(HitInfo& hitInfo, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const dd4hep::rec::Vector3D& simHitGlobalPos) const;
  
  /** @brief Calculate the entry point and path vector of a simHit in local sensor frame. */
  std::tuple<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D> ConstructSimHitPath(HitInfo& hitInfo, const edm4hep::SimTrackerHit& simHit) const;
  
  /** @brief Create a digitized hit*/
  void CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const dd4hep::rec::Vector3D& position, const float charge) const;


  // -- Pixel and in-pixel binning magic (logic) --

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
  std::tuple<int, int, int, int, int> computeSegmentIndices(HitInfo& hitInfo, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const int segment) const;


  // -- Transformation between global frame and local sensor frame -- 

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


  // -- Other helper functions --

  /** @brief Convert EDM4HEP vector to DD4HEP vector, and vice versa */
  dd4hep::rec::Vector3D convertVector(edm4hep::Vector3d vec) const;
  /** @copydoc convertVector(edm4hep::Vector3d) */
  edm4hep::Vector3d convertVector(dd4hep::rec::Vector3D vec) const;

  /** @brief Calculate the clipping factors for the path of a simHit in local sensor frame. Used in FindSimHitPath() */
  std::tuple<float, float> computePathClippingFactors(float t_min, float t_max, const float entryPos_ax, const float pathLength_ax, const float sensorLength_ax) const;

  /** @brief Write simHit & path information to a CSV file */
  void writeSimHitToCsv(const edm4hep::SimTrackerHit& simHit, const dd4hep::rec::Vector3D& simHitPos, const dd4hep::rec::Vector3D& simHitEntryPos, const dd4hep::rec::Vector3D& simHitPath, const int segmentNumber, const int layerIndex, const int i_u, const int i_v) const;
  

  // -- Properties --

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "VXD", "Name of the subdetector"};
  Gaudi::Property<bool> m_isStrip{this, "IsStrip", false, "Whether the hits are 1D strip hits"};
  Gaudi::Property<std::string> m_geometryServiceName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"}; // what is this for?
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  
  
  Gaudi::Property<bool> m_cutDistanceToSurface{this, "CutDistanceToSurface", false, "Whether to cut simHits where the simHitPosition is not exactly on the DD4hep sensor simSurface."};
  Gaudi::Property<bool> m_cutPathOutsideSensor{this, "CutPathOutsideSensor", false, "Whether to cut simHits where the simHitPath is not on or inside the sensor."};
  Gaudi::Property<float> m_cutDepositedCharge{this, "CutDepositedCharge", 0.0, "Minimum charge (e-) of SimTrackerHit to be digitized"};
  
  // Sensor pitch and size (enter either a single value for all layers or a vector, containing one value per layer)
  Gaudi::Property<std::vector<double>> m_pixelPitch_u{this, "PixelPitch_u", {0.025}, "Pixel pitch in direction of u in mm; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<double>> m_pixelPitch_v{this, "PixelPitch_v", {0.025}, "Pixel pitch in direction of v in mm; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<int>> m_pixelCount_u{this, "PixelCount_u", {1024}, "Number of pixels in direction of u; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<int>> m_pixelCount_v{this, "PixelCount_v", {1024}, "Number of pixels in direction of v; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<double>> m_sensorThickness{this, "SensorThickness", {0.05}, "Sensor thickness in mm; either one per layer or one for all layers"};
  
  Gaudi::Property<int> m_layerCount{this, "LayerCount", 5, "Number of layers in the subdetector. Used to validate the size of the pixel pitch and count vectors."};
  Gaudi::Property<std::vector<int>> m_layersToDigitize{this, "LayersToDigitize", {}, "Which layers to digitize (0-indexed). If empty, all layers are digitized."};

  Gaudi::Property<std::vector<int>> m_maxClusterSize{this, "MaximumClusterSize", {15,15}, "Maximum cluster size [max in u, max in v] in terms of pixels. Charges in pixels further that 1/2 max from the simHit pixel are discarded. This greatly improves performance. "};

  Gaudi::Property<float> m_targetPathSegmentLength{this, "TargetPathSegmentLength", 0.002, "Length of the path segments, that the simHits path through a sensor is divided into. In mm. Defines the precision of the charge deposition along the path."};
  Gaudi::Property<float> m_pathLengthShorteningFactorGeant4{this, "PathLengthShorteningFactorGeant4", 1.05, "Relative path length (to Geant4 path length), above which the path is shortened to the Geant4 length. (Geant4 length includes multiple scattering and curling in B-field, which are lost in our linear approximation of the path)."};
  
  Gaudi::Property<int> m_kernelSize{this, "KernelSize", 3, "Size of the charge spreading kernel imported from Allpix2 (ie. 3 for 3x3 kernel) Must be an odd integer >= 3."};
  Gaudi::Property<float> m_pixelThreshold{this, "PixelThreshold", 100, "Threshold in electrons for a pixel to fire (1 eh-pair = 3.65 eV)"};

  Gaudi::Property<int> m_inPixelBinCount_u{this, "InPixelBinCount_u", 3, "Number of bins per pixel in u direction for charge deposition. Must agree with the imported Kernel map."};
  Gaudi::Property<int> m_inPixelBinCount_v{this, "InPixelBinCount_v", 3, "Number of bins per pixel in v direction for charge deposition. Must agree with the imported Kernel map."};
  Gaudi::Property<int> m_inPixelBinCount_w{this, "InPixelBinCount_w", 3, "Number of bins per pixel in w (vertical) direction for charge deposition. Must agree with the imported Kernel map."};

  // Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo). Defaults to no transformation.
  Gaudi::Property<std::string> m_localNormalVectorDir{this, "LocalNormalVectorDir", "", "Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo). If defined correctly, the local frame is transformed such that z is orthogonal to the sensor plane."};

  Gaudi::Property<std::vector<float>> m_globalKernel{this, "GlobalKernel", {}, "Flat vector containing the global charge sharing kernel in row-major order (ie. row-by-row), starting on top left. Length must be KernelSize*KernelSize"};

  Gaudi::Property<std::string> m_debugCsvName{this, "DebugCsvName", "", "Name of a CSV file to output detailed per-hit debug information. If empty, no debug output is created."};


  // -- Services --
  
  SmartIF<IGeoSvc> m_geometryService;
  // Decoder for the cellID
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_cellIDdecoder;

  dd4hep::VolumeManager m_volumeManager; // volume manager to get the physical cell sensitive volume
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  // -- Global variables --

  mutable long int m_eventNumber; // current event number

  const float m_numericLimit_float = 0.00001f; // limit for floating point comparisons
  const double m_numericLimit_double = 0.000000001; // limit for floating point

  const float m_chargePerkeV = 273.97f; // number of electron-hole pairs created per keV of deposited energy in silicon. eh-pair ~ 3.65 eV
  
  std::unordered_map<int, int> m_layerToIndex; // layer number (from cellID & m_layersToDigitize) to internal index (0...N-1, where N is the number of layers to digitize)

  mutable bool m_debugFlag; // per-hit debug flag, set to true to add hit to debug histograms
  mutable std::ofstream m_debugCsvFile; // debug output file

  const dd4hep::rec::SurfaceMap* m_simSurfaceMap;
  
  enum { 
    counter_eventsRead, 
    counter_eventsRejected_noSimHits, 
    counter_eventsAccepted, 
    counter_simHitsRead, 
    counter_simHitsRejected_LayerNotToBeDigitized, 
    counter_simHitsRejected_ChargeCut, 
    counter_simHitsRejected_SurfaceDistToLarge, 
    counter_simHitsRejected_OutsideSensor,
    counter_simHitsRejected_NoSegmentsInSensor,
    counter_simHitsAccepted, 
    counter_digiHitsCreated, 
    counterArrayLen }; // counter indices.

  mutable std::array<long unsigned int, counterArrayLen> m_counters = {0}; // array of counters, size counterArrayLen

  enum { 
    hist_hitE, 
    hist_hitCharge,
    hist_clusterSize, 
    hist_EntryPointX, 
    hist_EntryPointY, 
    hist_EntryPointZ, 
    hist_DisplacementU, 
    hist_DisplacementV, 
    hist_DisplacementR, 
    hist_pathLength, 
    hist_pathLengthGeant4,
    hist_HitChargeDifference, 
    histArrayLen}; // histogram indices. histArrayLen must be last
  std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<1>>, histArrayLen> m_histograms; 

enum { 
  hist2d_hitMap_simHits,
  hist2d_hitMap_digiHits, 
  hist2d_hitMap_simHitDebug, 
  hist2d_hitMap_digiHitDebug, 
  hist2d_pathLength_vs_simHit_v,
  hist2dArrayLen }; // 2D histogram indices. hist2dArrayLen must be last
  std::vector<std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<2>>, hist2dArrayLen>> m_histograms2d;


  mutable std::unordered_set<int> m_initialSensorSizeCheckPassed; // whether the check that the pixel pitch and count match the sensor size in the geometry has been done and passed

  enum {
    histWeighted2d_averageCluster,
    histWeighted2dArrayLen };
  std::vector<std::array<std::unique_ptr<Gaudi::Accumulators::StaticWeightedHistogram<2, Gaudi::Accumulators::atomicity::full, double>>, histWeighted2dArrayLen>> m_histWeighted2d;



  struct PixelChargeMatrix {
    /* Stores the charge deposited in a (size_u x size_v) pixel matrix around a given origin pixel.
     * In case charge is added outside the matrix bounds, the matrix range is expanded in that direction.
     * The size of the matrix is defined via the Gaudi property MaximumClusterSize.
    */
    PixelChargeMatrix(int i_origin_u, int i_origin_v, int size_u, int size_v) 
      : m_origin{ i_origin_u, i_origin_v } {
        // At first, the range is centered around the origin
      m_range_u[0] = m_origin[0] - (size_u - 1) / 2;
      m_range_u[1] = m_origin[0] + (size_u - 1) / 2;
      m_range_v[0] = m_origin[1] - (size_v - 1) / 2;
      m_range_v[1] = m_origin[1] + (size_v - 1) / 2;

      m_pixelCharge.resize(size_u * size_v, 0.f);
    }
    
    inline int GetOriginU() const { return m_origin[0]; }
    inline int GetOriginV() const { return m_origin[1]; }
    inline int GetRangeMin_u() const { return m_range_u[0]; }
    inline int GetRangeMax_u() const { return m_range_u[1]; }
    inline int GetRangeMin_v() const { return m_range_v[0]; }
    inline int GetRangeMax_v() const { return m_range_v[1]; }
    inline int GetSizeU() const { return m_range_u[1] - m_range_u[0] + 1; }
    inline int GetSizeV() const { return m_range_v[1] - m_range_v[0] + 1; }

    inline std::tuple<int, int> GetSize() const {
      return std::make_tuple(GetSizeU(), GetSizeV());
    }
    inline std::tuple<int, int> GetOrigin() const {
      return std::make_tuple(m_origin[0], m_origin[1]);
    }
    
    float TotalCharge() const {
      return std::accumulate(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
    }

    float GetCharge(int i_u, int i_v) const {
      if (_OutOfBounds(i_u, i_v)) {
        throw std::runtime_error("PixelChargeMatrix::GetCharge: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") out of range");
      }
      return m_pixelCharge[_FindIndex(i_u, i_v)];
    }

    void FillCharge(int i_u, int i_v, float charge) {
      if (_OutOfBounds(i_u, i_v)) {
        // this might occur due to the limited size of the matrix or delta-electrons. 
        _ExpandMatrix(i_u, i_v);
      }
      if (_OutOfBounds(i_u, i_v)) {
        throw std::runtime_error("PixelChargeMatrix::FillCharge: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") still out of range ( " + std::to_string(m_range_u[0]) + ", " + std::to_string(m_range_u[1]) + ", " + std::to_string(m_range_v[0]) + ", " + std::to_string(m_range_v[1]) + ") after expansion");
      }
      m_pixelCharge[_FindIndex(i_u, i_v)] += charge;
    }
    // TODO: resize array if out of bounds (instead of simply discarding the charge)? This would be slow, but more robust. Is probably acceptable if it only occurs rarely, ie. the size_u and size_v are large enough for >99% of cases. ~ Jona 2025-10

    void Reset() {
      std::fill(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
    }

  private:
    const int m_minExpansionStep = 5; // number of pixels to expand the matrix by, if out of bounds occurs
  
    int m_range_u[2], m_range_v[2]; // Inclusive matrix range. -> size = range[1] - range[0] + 1
    // CAN theoretically extend into negative values, this ensures graceful handling of hits outside inditial bounds. These might be discarded later, if outside of sensor.
    int m_origin[2]; // origin pixel indices
    std::vector<float> m_pixelCharge;

    inline int _FindIndex(int i_u, int i_v) const {
      int i_u_rel = i_u - m_range_u[0];
      int i_v_rel = i_v - m_range_v[0];
      return i_u_rel + i_v_rel * GetSizeU();
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
        rangeNew_u[0] = i_u - m_minExpansionStep; // expand to include i_u plus some excess
      } else if (i_u > m_range_u[1]) {
        rangeNew_u[1] = i_u + m_minExpansionStep;
      }
      // expand in v direction?
      if (i_v < m_range_v[0]) {
        rangeNew_v[0] = i_v - m_minExpansionStep;
      } else if (i_v > m_range_v[1]) {
        rangeNew_v[1] = i_v + m_minExpansionStep;
      }

      const int newSizeU = rangeNew_u[1] - rangeNew_u[0] + 1;
      const int newSizeV = rangeNew_v[1] - rangeNew_v[0] + 1;

      std::vector<float> newPixelCharge(newSizeU * newSizeV, 0.f);

      // copy old charges into new array, row by row
      
      for (int row = 0; row < GetSizeV(); ++row) {
        std::copy(
          m_pixelCharge.begin() + row * GetSizeU(),
          m_pixelCharge.begin() + (row + 1) * GetSizeU(),
          newPixelCharge.begin() + (row + (m_range_v[0] - rangeNew_v[0])) * newSizeU + (m_range_u[0] - rangeNew_u[0])
        );
      }

      m_pixelCharge.swap(newPixelCharge);
      std::copy(rangeNew_u, rangeNew_u + 2, m_range_u);
      std::copy(rangeNew_v, rangeNew_v + 2, m_range_v);
    }
  };

  // class to store the charge sharing kernel for each in-pixel bin
  struct ChargeSharingKernels {
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
        throw std::runtime_error("GetKernelEntry: j_col (=" + std::to_string(j_col) + ") out of range");
      if (abs(j_row) > (m_kernelSize-1)/2)
        throw std::runtime_error("GetKernelEntry: j_row (=" + std::to_string(j_row) + ") out of range");
      const auto& kernel = GetKernel(j_u, j_v, j_w);
      return kernel[(j_col + (m_kernelSize-1)/2) * m_kernelSize + (j_row + (m_kernelSize-1)/2)];
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
    int m_binCountU, m_binCountV, m_binCountW;
    int m_kernelSize;
    /** Vector of kernels, one per in-pixel bin
     *  Kernel-indexing is col-major (ie. i = col * size + row) */
    std::vector<std::vector<float>> m_kernels; 

    int _FindIndex (int j_u, int j_v, int j_w) const {
      if (j_u < 0 || j_u >= m_binCountU)
        throw std::runtime_error("ChargeSharingKernel::_FindIndex: j_u (= " + std::to_string(j_u) + ") out of range");
      if (j_v < 0 || j_v >= m_binCountV)
        throw std::runtime_error("ChargeSharingKernel::_FindIndex: j_v (= " + std::to_string(j_v) + ") out of range");
      if (j_w < 0 || j_w >= m_binCountW)
        throw std::runtime_error("ChargeSharingKernel::_FindIndex: j_w (= " + std::to_string(j_w) + ") out of range");

      return j_u + m_binCountU * (j_v + m_binCountV * j_w); 
    }

  }; // struct ChargeSharingKernel

  std::unique_ptr<ChargeSharingKernels> m_chargeSharingKernels; // the charge sharing kernel
  // TODO: implement having a kernel per layer (array of unique_ptr ?)
};