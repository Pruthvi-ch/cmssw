// \class JetTracksAssociatorAtVertex JetTracksAssociatorAtVertex.cc
//
// Original Author:  Andrea Rizzi
//         Created:  Wed Apr 12 11:12:49 CEST 2006
// Accommodated for Jet Package by: Fedor Ratnikov Jul. 30, 2007
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

#include "JetTracksAssociatorAtVertex.h"

JetTracksAssociatorAtVertex::JetTracksAssociatorAtVertex(const edm::ParameterSet& fConfig)
    : mJets{consumes<edm::View<reco::Jet> >(fConfig.getParameter<edm::InputTag>("jets"))},
      mTracks{consumes<reco::TrackCollection>(fConfig.getParameter<edm::InputTag>("tracks"))},
      mAssociator(fConfig.getParameter<double>("coneSize")),
      mAssociatorAssigned(fConfig.getParameter<double>("coneSize")),
      useAssigned(fConfig.getParameter<bool>("useAssigned")),
      pvSrc() {
  if (useAssigned) {
    pvSrc = consumes<reco::VertexCollection>(fConfig.getParameter<edm::InputTag>("pvSrc"));
  }
  produces<reco::JetTracksAssociation::Container>();
}

void JetTracksAssociatorAtVertex::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag(""));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<double>("coneSize", 0.4);
  desc.add<bool>("useAssigned", false);
  desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlinePrimaryVertices"));
  descriptions.addWithDefaultLabel(desc);
}

void JetTracksAssociatorAtVertex::produce(edm::Event& fEvent, const edm::EventSetup& fSetup) {
  edm::Handle<edm::View<reco::Jet> > jets_h;
  fEvent.getByToken(mJets, jets_h);
  edm::Handle<reco::TrackCollection> tracks_h;
  fEvent.getByToken(mTracks, tracks_h);

  auto jetTracks = std::make_unique<reco::JetTracksAssociation::Container>(reco::JetRefBaseProd(jets_h));

  // format inputs
  std::vector<edm::RefToBase<reco::Jet> > allJets;
  allJets.reserve(jets_h->size());
  for (unsigned i = 0; i < jets_h->size(); ++i)
    allJets.push_back(jets_h->refAt(i));
  std::vector<reco::TrackRef> allTracks;
  allTracks.reserve(tracks_h->size());
  // run algo
  for (unsigned i = 0; i < tracks_h->size(); ++i) {
    allTracks.push_back(reco::TrackRef(tracks_h, i));
  }
  if (!useAssigned) {
    mAssociator.produce(&*jetTracks, allJets, allTracks);
  } else {
    edm::Handle<reco::VertexCollection> pvHandle;
    fEvent.getByToken(pvSrc, pvHandle);
    const reco::VertexCollection& vertices = *pvHandle.product();
    mAssociatorAssigned.produce(&*jetTracks, allJets, allTracks, vertices);
  }

  // store output
  fEvent.put(std::move(jetTracks));
}
