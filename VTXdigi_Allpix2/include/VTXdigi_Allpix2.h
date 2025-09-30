#pragma once

// GAUDI
// #include "Gaudi/Algorithm.h"
// #include "GaudiKernel/IRndmGenSvc.h"
// #include "GaudiKernel/RndmGenerators.h"
#include "Gaudi/Property.h"

#include "Gaudi/Accumulators/RootHistogram.h" // added by Jona

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

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator() (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const override;

private:


  // -- Transformation between global frame and local sensor frame

  /** Get the transformation matrix to the local sensor frame */
  TGeoHMatrix GetTransformationMatrix(const edm4hep::SimTrackerHit& simHit) const;
  TGeoHMatrix GetTransformationMatrix(const dd4hep::DDSegmentation::CellID& cellID) const;

  /** adjust the transformation matrix to the local sensor frame, so that it is x-u, y-v, z-n and right-handed. Uses user-input "LocalNormalVectorDir" */
  void SetProperDirectFrame(TGeoHMatrix& sensorTransformMatrix) const;

  /** Transform a global position to the local sensor frame of a hit, given its cellID */
  dd4hep::rec::Vector3D GlobalToLocal(const dd4hep::rec::Vector3D& globalPos, const dd4hep::DDSegmentation::CellID& cellID) const;

  /** Transform a local sensor position to the global frame, given the sensors cellID */
  dd4hep::rec::Vector3D LocalToGlobal(const dd4hep::rec::Vector3D& localPos, const dd4hep::DDSegmentation::CellID& cellID) const;


  // -- Helper functions --

  /** Convert EDM4HEP vector to DD4HEP vector, and vice versa */
  dd4hep::rec::Vector3D Vector3dConvert(edm4hep::Vector3d vec) const;
  edm4hep::Vector3d Vector3dConvert(dd4hep::rec::Vector3D vec) const;

  /** Check that the pixel pitch and count match the sensor size in the geometry
   */
  void InitialSensorSizeCheck(const edm4hep::SimTrackerHit& simHit) const;
    
  /** Apply cuts to a simHit. Returns true if the hit is accepted, false if dismissed. */
  bool ApplySimHitCuts (int layer, double eDeposited) const;

  /** Calculate the entry point and path vector of a simHit in local sensor frame. */
  std::tuple<dd4hep::rec::Vector3D, dd4hep::rec::Vector3D> GetSimHitPath(const edm4hep::SimTrackerHit& simHit) const;

  /**Given a histogram definition (x0, binWidth, nBins) and a value x, return the bin index in which x falls.
   * Returns -1 if x is out of range.
   * Bins are 0-indexed (vs ROOT's 1-indexing) */
  int GetBinIndex(float x, float binX0, float binWidth, int binN) const;

  /** Create a digitized hit */
  void CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const double eDeposition, const dd4hep::rec::Vector3D& position) const;
  void CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitsLinks, const double eDeposition, const edm4hep::Vector3d& position) const;
  
  /** Get the local position of a pixel center (u,v,0) */
  dd4hep::rec::Vector3D GetPixelCenter_Local(const int& i_u, const int& i_v, const int& layer, const dd4hep::rec::ISurface& surface) const;



  // -- Properties --

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "VXD", "Name of the subdetector"};
  Gaudi::Property<bool> m_isStrip{this, "IsStrip", false, "Whether the hits are 1D strip hits"};
  Gaudi::Property<double> m_minDepositedEnergy{this, "MinDepositedEnergy", 0.0, "Minimum energy (GeV) of SimTrackerHit to be digitized"};
  Gaudi::Property<double> m_targetPathSegmentLength{this, "TargetPathSegmentLength", 0.002, "Length of the path segments, that the simHits path through a sensor is divided into. In mm. Defines the precision of the charge deposition along the path."};

  Gaudi::Property<std::string> m_geometryServiceName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"}; // what is this for?

  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  Gaudi::Property<std::vector<int>> m_layersToDigitize{this, "LayersToDigitize", {}, "Which layers to digitize (0-indexed). If empty, all layers are digitized."};

  // Sensor pitch and size (enter either a single value for all layers or a vector, containing one value per layer)
  Gaudi::Property<std::vector<float>> m_pixelPitch_u{this, "PixelPitch_u", {0.025}, "Pixel pitch in direction of u in mm; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<float>> m_pixelPitch_v{this, "PixelPitch_v", {0.025}, "Pixel pitch in direction of v in mm; either one per layer or one for all layers"};

  Gaudi::Property<std::vector<int>> m_pixelCount_u{this, "PixelCount_u", {1024}, "Number of pixels in direction of u; either one per layer or one for all layers"};
  Gaudi::Property<std::vector<int>> m_pixelCount_v{this, "PixelCount_v", {1024}, "Number of pixels in direction of v; either one per layer or one for all layers"};

  Gaudi::Property<int> m_layerCount{this, "LayerCount", 5, "Number of layers in the subdetector. Used to validate the size of the pixel pitch and count vectors."};

  Gaudi::Property<int> m_KernelSize{this, "KernelSize", 3, "Size of the charge spreading kernel imported from Allpix2 (ie. 3 for 3x3 kernel) Must be an odd integer >= 3."};
  Gaudi::Property<float> m_Threshold{this, "Threshold", 0.6, "Threshold in keV for a pixel to fire (1 keV = 274 eh-pairs)"};

  Gaudi::Property<int> m_InPixelBinCount_u{this, "InPixelBinCount_u", 3, "Number of bins per pixel in u direction for charge deposition. Must agree with the imported Kernel map."};
  Gaudi::Property<int> m_InPixelBinCount_v{this, "InPixelBinCount_v", 3, "Number of bins per pixel in v direction for charge deposition. Must agree with the imported Kernel map."};
  Gaudi::Property<int> m_InPixelBinCount_w{this, "InPixelBinCount_w", 3, "Number of bins per pixel in w (vertical) direction for charge deposition. Must agree with the imported Kernel map."};

  // Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo). Defaults to no transformation.
  Gaudi::Property<std::string> m_localNormalVectorDir{this, "LocalNormalVectorDir", "", "Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo). If defined correctly, the local frame is transformed such that z is orthogonal to the sensor plane."};

  // -- Services --
  
  SmartIF<IGeoSvc> m_geometryService;
  // Decoder for the cellID
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_cellIDdecoder;

  dd4hep::VolumeManager m_volumeManager; // volume manager to get the physical cell sensitive volume
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  // -- Global variables --

  const dd4hep::rec::SurfaceMap* m_surfaceMap;

  enum { hist_hitE, hist_clusterSize, hist_EntryPointX, hist_EntryPointY, hist_EntryPointZ, hist_DisplacementU, hist_DisplacementV, hist_DisplacementR, histArrayLen }; // histogram indices. histArrayLen must be last

  std::array<std::unique_ptr<Gaudi::Accumulators::StaticRootHistogram<1>>, histArrayLen> m_histograms; 

  mutable std::unordered_set<int> m_initialSensorSizeCheckPassed; // whether the check that the pixel pitch and count match the sensor size in the geometry has been done and passed

  // class to store the charge sharing kernel for each in-pixel bin
  struct ChargeSharingKernels {
    ChargeSharingKernels(int& binCountU, int& binCountV, int& binCountW, int& kernelSize) : m_binCountU(binCountU), m_binCountV(binCountV), m_binCountW(binCountW), m_kernelSize(kernelSize) {
      if (kernelSize < 3 || kernelSize % 2 == 0)
        throw std::runtime_error("ChargeSharingKernel: Kernel size must be an odd integer >= 3, but is " + std::to_string(kernelSize) + ".");

      m_kernels.resize(binCountU * binCountV * binCountW);
      for (auto& kernel : m_kernels) {
        kernel.resize(kernelSize*kernelSize, 0.);
      }
    }

    /** Set the charge sharing kernel for a specific in-pixel bin 
     * @param values A flat vector containing the kernel values in row-major order (ie. row-by-row) (length must be kernelSize*kernelSize)
     */
    void SetKernel(const int j_u, const int j_v, const int j_w, const std::vector<double>& values) {
      if (static_cast<int>(values.size()) != m_kernelSize*m_kernelSize)
        throw std::runtime_error("ChargeSharingKernel::SetKernel: values size (" + std::to_string(values.size()) + ") does not match kernel size (" + std::to_string(m_kernelSize*m_kernelSize) + ")");

      const int index = _index(j_u, j_v, j_w);
      std::vector<double>& kernel = m_kernels[index];

      for (int row = 0; row < m_kernelSize; ++row) {
        for (int col = 0; col < m_kernelSize; ++col) {
          kernel[col * m_kernelSize + row] = values[row*m_kernelSize + col];
        }
      }
    }

    /** Access kernel as const reference */
    const std::vector<double>& GetKernel(int j_u, int j_v, int j_w) const {
      return m_kernels[_index(j_u, j_v, j_w)];
    }

    /** Access a specific entry of a kernel */
    double GetKernelEntry(int j_u, int j_v, int j_w, int j_col, int j_row) const {
      if (j_col < 0 || j_col >= m_kernelSize || j_row < 0 || j_row >= m_kernelSize)
        throw std::runtime_error("GetKernelEntry: col/row out of range");
      const auto& kernel = GetKernel(j_u, j_v, j_w);
      return kernel[j_col * m_kernelSize + j_row];
    }

  private:
    int m_binCountU, m_binCountV, m_binCountW;
    int m_kernelSize;
    /** Vector of kernels, one per in-pixel bin
     *  Kernel-indexing is col-major (ie. i = col * size + row) */
    std::vector<std::vector<double>> m_kernels; 

    int _index (int j_u, int j_v, int j_w) const {
      if (j_u < 0 || j_u >= m_binCountU)
        throw std::runtime_error("ChargeSharingKernel::SetKernel: j_u (" + std::to_string(j_u) + ") out of range");
      if (j_v < 0 || j_v >= m_binCountV)
        throw std::runtime_error("ChargeSharingKernel::SetKernel: j_v (" + std::to_string(j_v) + ") out of range");
      if (j_w < 0 || j_w >= m_binCountW)
        throw std::runtime_error("ChargeSharingKernel::SetKernel: j_w (" + std::to_string(j_w) + ") out of range");

      return j_u + m_binCountU * (j_v + m_binCountV * j_w); 
    }

  }; // struct ChargeSharingKernel

  std::unique_ptr<ChargeSharingKernels> m_chargeSharingKernels; // the charge sharing kernel
  // TODO: implement having a kernel per layer (array of unique_ptr ?)
};