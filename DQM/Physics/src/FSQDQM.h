#ifndef FSQDQM_H
#define FSQDQM_H

#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <map>
#include <TMath.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TComplex.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TComplex.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

//Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
//Candidate
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/CaloTauFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TProfile.h>
#include "DQMServices/Core/interface/MonitorElement.h" 
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

class DQMStore;

class FSQDQM : public edm ::EDAnalyzer{

 public:
 FSQDQM(const edm::ParameterSet& ps);
  virtual ~FSQDQM();

 protected:

  void beginJob();
  void beginRun(edm::Run const& run, 
                edm::EventSetup const& eSetup);
  void analyze(edm::Event const& e, 
               edm::EventSetup const& eSetup);
  void beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                            edm::EventSetup const& context) ;
  void endRun(edm::Run const& run, 
              edm::EventSetup const& eSetup);
 private:
  //  float deltaPhi(float phi1, float phi2);
  void bookHistos(DQMStore* bei);

  DQMStore* bei_; 
  /*
  HLTConfigProvider hltConfigProvider_;
  bool isValidHltConfig_;

  edm::InputTag theTriggerResultsCollection;
  edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  edm::InputTag HLTriggerResults_;
  std::vector<std::string> triggers_;
  std::string labelVtx_, labelBS_, labelTrack_,labelPFJet_,labelCastorJet_;

  edm::EDGetTokenT<edm::TriggerResults>                 tok_trigRes_;
  edm::EDGetTokenT<reco::BeamSpot>                        tok_bs_;
  edm::EDGetTokenT<reco::VertexCollection>                tok_Vtx_;
  edm::EDGetTokenT<reco::TrackCollection>                tok_track_;
  std::vector<std::string> all_triggers;
  edm::EDGetTokenT<reco::PFJetCollection> tok_pfjet_;
  edm::EDGetTokenT<reco::BasicJetCollection> tok_castorjet_;
  
  std::vector<int>          hltresults;
  unsigned int runNumber_, eventNumber_ , lumiNumber_, bxNumber_;

  //Histograms
  
  MonitorElement *PFJetpt;
  MonitorElement *PFJeteta;
  MonitorElement *PFJetphi;


  MonitorElement *CastorJetphi;
  MonitorElement *CastorJetMulti;
  MonitorElement *PFJetMulti;
  MonitorElement *PFJetRapidity;
  MonitorElement *Track_HP_Phi;
  MonitorElement *Track_HP_Eta;
  MonitorElement *Track_HP_Pt;
  MonitorElement *Track_HP_ptErr_over_pt;
  MonitorElement *Track_HP_dzvtx_over_dzerr;
  MonitorElement *Track_HP_dxyvtx_over_dxyerror;
  MonitorElement *NPV;
  MonitorElement *PV_chi2;
  MonitorElement *PV_d0;
  MonitorElement *PV_numTrks;
  MonitorElement *PV_sumTrks;
  MonitorElement *h_ptsum_towards;
  MonitorElement *h_ptsum_transverse;
  MonitorElement *h_ptsum_away;
  MonitorElement *h_ntracks_towards;
  MonitorElement *h_ntracks_transverse;
  MonitorElement *h_ntracks_away;
  
  TProfile *h_leadingtrkpt_ntrk_away;
  TProfile *h_leadingtrkpt_ntrk_towards;
  TProfile *h_leadingtrkpt_ntrk_transverse;
  TProfile *h_leadingtrkpt_ptsum_away;
  TProfile *h_leadingtrkpt_ptsum_towards;
  TProfile *h_leadingtrkpt_ptsum_transverse;
  */
};
#endif
