#include "VTXdigi_Allpix2.h"

DECLARE_COMPONENT(VTXdigi_Allpix2)

VTXdigi_Allpix2::VTXdigi_Allpix2(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                          KeyValues("SimTrackHitCollectionName", {"SimTrackerHits"}),
                          KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("TrackerHitCollectionName", {"VTXTrackerHits"}),
                        KeyValues("SimTrkHitRelationsCollection", {"VTXTrackerHitRelations"})}) {

  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }  
}


StatusCode VTXdigi_Allpix2::initialize() {
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> VTXdigi_Allpix2::operator()
  (const edm4hep::SimTrackerHitCollection& simTrackerHits, const edm4hep::EventHeaderCollection& headers) const {
    
  auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  debug() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run " << headers[0].getRunNumber() << endmsg;

  return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());
}

