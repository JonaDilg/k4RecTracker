#pragma once

#include "Gaudi/Property.h"

#include "Gaudi/Accumulators/RootHistogram.h" // added by Jona
#include "Gaudi/Accumulators/Histogram.h" // added by Jona

#include "podio/LinkCollection.h"
#include "podio/LinkNavigator.h"

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

/** @class Clustering_Pixels
 *  @author Jona Dilg, Armin Ilg
 *  @date   2025-09 */

struct Clustering_Pixels final 
  : k4FWCore::MultiTransformer 
    <std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> (const edm4hep::TrackerHitPlaneCollection&,  const edm4hep::TrackerHitSimTrackerHitLinkCollection&, const edm4hep::EventHeaderCollection&)> {

  Clustering_Pixels(const std::string& name, ISvcLocator* svcLoc);
  
  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator() (const edm4hep::TrackerHitPlaneCollection& hits, const edm4hep::TrackerHitSimTrackerHitLinkCollection& hitLinks, const edm4hep::EventHeaderCollection& headers) const override;

private:

  /* -- Classes -- */

  struct HitData {
    float u = 0, v = 0;
    float eDep = 0;
    float time = 0;
    std::shared_ptr<const edm4hep::TrackerHitPlane> hitPtr = nullptr;
  };

  /* -- Initialization / finalization functions -- */

  void InitializeServicesAndGeometry();
  void CheckGaudiProperties();
  void InitHistograms();
  void PrintCountersSummary() const;

  /* -- Core algorithm functions -- */

  bool CheckInitialSetup(const edm4hep::TrackerHitPlaneCollection& hits, const edm4hep::EventHeaderCollection& headers) const;


  void FillHistograms_perCluster(const edm4hep::MutableTrackerHitPlane& cluster, const std::vector<Clustering_Pixels::HitData>& clusterHitsData,  std::vector<edm4hep::SimTrackerHit> clusterSimHits, const dd4hep::DDSegmentation::CellID& cellID, const int layer) const;


  /* -- Helper functions -- */

  int ComputeBinIndex(float x, float binX0, float binWidth, int binN) const;

  TGeoHMatrix ComputeTransformationMatrix(const dd4hep::DDSegmentation::CellID& cellID) const;
  inline TGeoHMatrix ComputeTransformationMatrix(const edm4hep::TrackerHitPlane& hit) const {
    return ComputeTransformationMatrix(hit.getCellID());
  };

  void SetProperDirectFrame(TGeoHMatrix& transformationMatrix) const;

  dd4hep::rec::Vector3D TransformGlobalToLocal(const dd4hep::rec::Vector3D& globalPos, const dd4hep::DDSegmentation::CellID& cellID) const;

  dd4hep::rec::Vector3D TransformLocalToGlobal(const dd4hep::rec::Vector3D& localPos, const dd4hep::DDSegmentation::CellID& cellID) const;

  dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec) const;
  edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec) const;

  bool CheckHitCuts(const edm4hep::TrackerHitPlane& hit) const;

  void CreateCluster() const;

  /* -- Services --*/

  SmartIF<IGeoSvc> m_geometryService;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_cellIDdecoder; // Decoder for the cellID

  dd4hep::VolumeManager m_volumeManager; // volume manager to get the physical cell sensitive volume
  const dd4hep::Detector* m_detector = nullptr; // pointer to the DD4hep detector
  dd4hep::DetElement m_subDetector; // subdetector DetElement. contains layers as children
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  const dd4hep::rec::SurfaceMap* m_simSurfaceMap;

  /* -- Constants -- */

  const float m_chargePerkeV = 273.97f; // number of electron-hole pairs created per keV of deposited energy in silicon. eh-pair ~ 3.65 eV

  /* -- Detector & sensor parameters -- */

  int m_layerN = -1; // number of layers to digitize. -1 means all layers

  float m_neighborRadiusHelper[3]; // filled in initialize(). Used for calculating distannces between hits


  /* -- Counters -- */

  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRead{this, "Events read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRejected_noHits{this, "Events rejected (no hits)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRejected_noHitsOnSpecifiedLayers{this, "Events rejected (no hits on specified layers)"};

  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRejected_noHitsAboveSeedThreshold{this, "Events rejected (no hits above seed threshold)"};

  mutable Gaudi::Accumulators::Counter<> m_counter_hitsRead{this, "Tracker hits read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_hitsRejected_layerIgnored{this, "Tracker hits rejected (on layer that is ignored)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_hitsAccepted{this, "Tracker hits accepted"};

  /* -- Gaudi properties -- */

  const std::string m_undefinedString = "UNDEFINED";

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", m_undefinedString, "Name of the subdetector"};
  Gaudi::Property<std::string> m_subDetChildName{this, "SubDetectorChildName", m_undefinedString, "Name of the subdetector child (eg. \"VertexBarrel\"), if applicable. If undefined, the subdetector itself is assumed to contain layers as children."};

  Gaudi::Property<std::string> m_geometryServiceName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"}; // what is this for?
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  Gaudi::Property<std::string> m_localNormalVectorDir{this, "LocalNormalVectorDir", "", "Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo). If defined correctly, the local frame is transformed such that z is orthogonal to the sensor plane."};
  Gaudi::Property<std::vector<int>> m_layersToRun{this, "LayersToRun", {}, "Which layers to digitize (0-indexed). If empty, all layers are digitized."};

  /* Notes on the radius:
   *   - I use radius instead of bins to avoid the complexity of converting position to bin. this is not necessarily the best decision, but it simplifies the code.
   *   - use pixel-pitch to only use direct neighbors, 1.5*pitch to also include diagonal neighbors. Higher numbers give a larger area, always elliptic. */
  Gaudi::Property<float> m_neighborRadius_u{this, "NeighborRadius_u", 37.5, "Radius (in um) in which to look for neighbors of a pixel, in u direction"};
  Gaudi::Property<float> m_neighborRadius_v{this, "NeighborRadius_v", 37.5, "Radius (in um) in which to look for neighbors of a pixel, in v direction"};

  /* Debugging */
  Gaudi::Property<bool> m_debugHistograms{this, "DebugHistograms", false, "Whether to create and fill debug histograms. Not recommended for multithreading, might lead to crashes. Requires RootHistSvc."};

  enum {
    hist1d_clusterSize,
    hist1d_clusterCharge,
    hist1d_simHitsPerCluster,
    hist1d_residual,
    hist1d_residual_z,
    hist1d_residual_u,
    hist1d_residual_v,
    hist1d_residual_w,
    hist1dArrayLen }; // 1D histogram indices
  mutable std::unordered_map<
    int,
    std::array<
      std::unique_ptr<
        Gaudi::Accumulators::StaticHistogram<
          1, 
          Gaudi::Accumulators::atomicity::full,
          float
        >
      >,
      hist1dArrayLen
    >
  > m_hist1d;

  enum {
    hist2d_clusterSize_vs_clusterCharge,
    hist2d_residual_local,
    hist2dArrayLen };
  mutable std::unordered_map<
    int,
    std::array<
      std::unique_ptr<
        Gaudi::Accumulators::StaticHistogram<
          2, 
          Gaudi::Accumulators::atomicity::full,
          float
        >
      >,
      hist2dArrayLen
    >
  > m_hist2d;

  enum{
    histProfile1d_clusterSize_vs_hit_z,
    histProfile1d_clusterSize_vs_hit_cosTheta,
    histProfile1d_clusterSize_vs_module_z,
    histProfile1d_clusterSize_vs_module_ID,
    histProfile1d_residual_vs_hit_z,
    histProfile1d_residual_u_vs_hit_z,
    histProfile1d_residual_v_vs_hit_z,
    histProfile1d_residual_vs_hit_cosTheta,
    histProfile1d_residual_vs_clusterSize,
    histProfile1dArrayLen };
  mutable std::unordered_map<
    int,
    std::array<
      std::unique_ptr<
        Gaudi::Accumulators::StaticProfileHistogram<
          1, 
          Gaudi::Accumulators::atomicity::full,
          float
        >
      >,
      histProfile1dArrayLen
    >
  > m_histProfile1d;
};

