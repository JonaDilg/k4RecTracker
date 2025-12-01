#pragma once

#include "IChargeCollector.h" 

#include "Gaudi/Property.h"

#include "Gaudi/Accumulators/RootHistogram.h" // added by Jona
#include "Gaudi/Accumulators/Histogram.h" // added by Jona

#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h" // added by Jona
#include "edm4hep/TrackerHitPlaneCollection.h" // added by Jona
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h" // added by Jona

// K4FWCORE
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

/** @class VTXdigi_Modular
 *
 * Creates trackerHits from simHits. Produces clusters from simHits, outputs either the cluster centre or all hits in the cluster as digitized hits.
 *
 *  @author Jona Dilg, Armin Ilg
 *  @date   2025-09
 */

/* -- Forward declarations -- */
namespace VTXdigi_details {
  class IChargeCollector;
}

struct VTXdigi_Modular final 
  : k4FWCore::MultiTransformer 
    <std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> (const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

  VTXdigi_Modular(const std::string& name, ISvcLocator* svcLoc);
  
  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator() (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const override;

private: 

  friend class IChargeCollector;
  
  /* -- Properties -- */

  const std::string m_undefinedString = "UNDEFINED";

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", m_undefinedString, "Name of the subdetector (eg. \"Vertex\")"};
  Gaudi::Property<std::string> m_subDetChildName{this, "SubDetectorChildName", m_undefinedString, "Name of the subdetector child (eg. \"VertexBarrel\"), if applicable. If undefined, the subdetector itself is assumed to contain layers as children."};

  Gaudi::Property<std::string> m_geometryServiceName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"}; // what is this for?
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  /* -- Properties and members related to the various charge collection algorithms-- */

  Gaudi::Property<std::string> m_chargeCollectionMethod{this, "ChargeCollectionMethod", "Drift", "Method used for charge collection: \"Fast\", \"Drift\", \"LookupTable\", etc."};
  
  /* -- Services, geometry variables -- */
  
  SmartIF<IGeoSvc> m_geoService;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_cellIDdecoder;

  dd4hep::VolumeManager m_volumeManager; // volume manager to get the physical cell sensitive volume
  const dd4hep::Detector* m_detector = nullptr;
  dd4hep::DetElement m_subDetector; // subdetector DetElement. contains layers as children
  SmartIF<IUniqueIDGenSvc> m_uidSvc;
  
  const dd4hep::rec::SurfaceMap* m_SurfaceMap;
  
  /* -- Constants -- */
  
  const float m_chargePerkeV = 273.97f; // number of electron-hole pairs created per keV of deposited energy in silicon. eh-pair ~ 3.65 eV
  
  /* -- Member variables -- */

  std::unique_ptr<VTXdigi_details::IChargeCollector> m_chargeCollector;
}; // class VTXdigi_Modular





