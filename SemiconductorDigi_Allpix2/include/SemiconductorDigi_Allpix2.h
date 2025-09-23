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

#include <string> // added by Jona
#include <vector>

/** @class SemiconductorDigi_Allpix2
 *
 * Creates TrackerHits from SimTrackerHits. Produces clusters from simHits, outputs either the cluster centre or all hits in the cluster as digitized hits.
 *
 *  @author Jona Dilg, Armin Ilg
 *  @date   2025-09
 *
 */



struct SemiconductorDigi_Allpix2 final : k4FWCore::MultiTransformer <std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> (const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

  SemiconductorDigi_Allpix2(const std::string& name, ISvcLocator* svcLoc);
  
  StatusCode initialize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator()(const edm4hep::SimTrackerHitCollection& simTrackerHits, const edm4hep::EventHeaderCollection& headers) const override;

private:
  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "VXD", "Name of the subdetector"};
  Gaudi::Property<bool> m_isStrip{this, "IsStrip", false, "Whether the hits are 1D strip hits"};

  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;
};

StatusCode SemiconductorDigi_Allpix2::initialize() {
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> SemiconductorDigi_Allpix2::operator()
  (const edm4hep::SimTrackerHitCollection& simTrackerHits, const edm4hep::EventHeaderCollection& headers) const {
    
  auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  debug() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run " << headers[0].getRunNumber() << endmsg;

  return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());
}