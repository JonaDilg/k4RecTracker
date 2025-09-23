#include "SemiconductorDigi_Allpix2.h"

DECLARE_COMPONENT(SemiconductorDigi_Allpix2)

SemiconductorDigi_Allpix2::SemiconductorDigi_Allpix2(const std::string& name, ISvcLocator* svcLoc)
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

