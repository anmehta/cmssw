#include "DQM/Physics/src/FSQDQM.h"
#include <memory>
using namespace std;
using namespace reco;
using namespace edm;

struct SortByPt

{

  bool operator () (const TLorentzVector& a, const TLorentzVector& b) const {
    
    return a.Pt() > b.Pt();
    
  }
  
};
/*
float FSQDQM::deltaPhi(float phi1, float phi2){
  float result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}
*/

FSQDQM::FSQDQM(const edm::ParameterSet& iConfig)
{
  edm::LogInfo("FSQDQM") << " Creating FSQDQM "<< "\n";
  /*  
  HLTriggerResults_ = iConfig.getParameter<edm::InputTag>("HLTriggerResults");
  labelBS_          = iConfig.getParameter<std::string>("LabelBeamSpot");
  labelVtx_         = iConfig.getParameter<std::string>("LabelVertex");
  labelPFJet_         = iConfig.getParameter<std::string>("LabelPFJet");
  labelCastorJet_         = iConfig.getParameter<std::string>("LabelCastorJet");
  theTriggerResultsCollection = iConfig.getParameter<InputTag>("triggerResultsCollection");
  labelTrack_       = iConfig.getParameter<std::string>("LabelTrack");
  tok_trigRes_      = consumes<edm::TriggerResults>(HLTriggerResults_);
  tok_bs_           = consumes<reco::BeamSpot>(labelBS_);
  tok_pfjet_          = consumes<reco::PFJetCollection>(labelPFJet_);
  tok_track_        = consumes<reco::TrackCollection>(labelTrack_);
  tok_Vtx_          = consumes<reco::VertexCollection>(labelVtx_);
  tok_castorjet_          = consumes<reco::BasicJetCollection>(labelCastorJet_);
  isValidHltConfig_ = false;
  */ 
}


FSQDQM::~FSQDQM()
{
 
  edm::LogInfo("FSQDQM") << " Deleting FSQDQM "
			   << "\n";
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
void FSQDQM::beginJob(){
  cout<<"Enterin FSQDQM::beginJob: "<<endl;
  bei_->setCurrentFolder("Physics/FSQ/TESTFSQ");
  cout<<"created dir "<<endl;
  bookHistos(bei_);
  cout<<"histos booked"<<endl;

  

  //cout<<"...leaving FSQDQM::beginJob. "<<endl;
}
//
// -- Begin Run
//
void FSQDQM::beginRun(Run const& run, edm::EventSetup const& eSetup) {
  edm::LogInfo ("FSQDQM") <<"[FSQDQM]: Begining of Run";
  // passed as parameter to HLTConfigProvider::init(), not yet used
  //bool isConfigChanged = false;
  
  // isValidHltConfig_ used to short-circuit analyze() in case of problems
  //const std::string hltProcessName( "HLT" );
  // isValidHltConfig_ = hltConfigProvider_.init( run, eSetup, hltProcessName, isConfigChanged );

}
//
// -- Begin  Luminosity Block
//

void FSQDQM::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
				    edm::EventSetup const& context) {
  cout<<"Entering FSQDQM::beginLuminosityBlock: "<<endl;
  
  edm::LogInfo ("FSQDQM") <<"[FSQDQM]: Begin of LS transition";

  cout<<"...leaving FSQDQM::beginLuminosityBlock. "<<endl;
}


//
//  -- Book histograms
//
void FSQDQM::bookHistos(DQMStore* bei){
  //  bei->cd();
  bei->setCurrentFolder("Physics/FSQ/TESTFSQ");

  //void FSQDQM::bookHistograms(DQMStore::IBooker & ibooker,edm::Run const &, edm::EventSetup const & ){
  //  ibooker.setCurrentFolder("Physics/FSQ");
  /*
  PFJetpt = bei->book1D("PFJetpt",";p_{T}(PFJet)", 100,0.0 , 100);
  PFJeteta = bei->book1D("PFJeteta", ";#eta(PFJet)", 50, -2.5, 2.5);
  PFJetphi = bei->book1D("PFJetphi", ";#phi(PFJet)", 50, -3.14,3.14);
  PFJetMulti         = bei->book1D("PFJetMulti", ";No. of PFJets", 10, -0.5, 9.5);
  PFJetRapidity      = bei->book1D("PFJetRapidity",";PFJetRapidity", 50, -6.0,6.0);
  CastorJetphi = bei->book1D("CastorJetphi", ";#phi(CastorJet)", 50, -3.14,3.14);
  CastorJetMulti         = bei->book1D("CastorJetMulti", ";No. of CastorJets", 10, -0.5, 9.5);
  Track_HP_Phi       =bei->book1D("Track_HP_Phi",";#phi(HPtrack)", 50, -3.14,3.14);
  Track_HP_Eta       =bei->book1D("Track_HP_Eta", ";#eta(HPtrack)", 50, -2.5, 2.5);
  Track_HP_Pt        =bei->book1D("Track_HP_Pt",  ";p_{T}(HPtrack)",500, 0.0 , 50);
  Track_HP_ptErr_over_pt=bei->book1D("Track_HP_ptErr_over_pt",";{p_{T}Err}/p_{T}",100,0,0.1);
  Track_HP_dzvtx_over_dzerr=bei->book1D("Track_HP_dzvtx_over_dzerr",";dZerr/dZ",100,-10,10);
  Track_HP_dxyvtx_over_dxyerror =bei->book1D("Track_HP_dxyvtx_over_dxyerror",";dxyErr/dxy",100,-10,10);
  NPV      = bei->book1D("NPV",";NPV",10, -0.5, 9.5);
  PV_chi2  = bei->book1D("PV_chi2",";PV_chi2",100, 0.0, 2.0);
  PV_d0    = bei->book1D("PV_d0",";PV_d0",100, -10.0, 10.0);
  PV_numTrks         = bei->book1D("PV_numTrks",";PV_numTrks",100, -0.5, 99.5);
  PV_sumTrks=bei->book1D("PV_sumTrks",";PV_sumTrks",100,0,100);
  h_ptsum_towards      = bei->book1D("h_ptsum_towards",";h_ptsum_towards",100,0,100);
  h_ptsum_transverse      = bei->book1D("h_ptsum_transverse",";h_ptsum_transverse",100,0,100);
  h_ptsum_away      = bei->book1D("h_ptsum_away",";h_ptsum_away",100,0,100);
  h_ntracks_towards      = bei->book1D("h_ntracks_towards",";h_ntracks_towards",50,-0.5,49.5);
  h_ntracks_transverse      = bei->book1D("h_ntracks_transverse",";h_ntracks_transverse",50,-0.5,49.5);
  h_ntracks_away      = bei->book1D("h_ntracks_away",";h_ntracks_away",50,-0.5,49.5);

  h_leadingtrkpt_ntrk_away = fs->make<TProfile>("h_leadingtrkpt_ntrk_away","h_leadingtrkpt_ntrk_away",50,0,50,0,30);
  h_leadingtrkpt_ntrk_towards = fs->make<TProfile>("h_leadingtrkpt_ntrk_towards","h_leadingtrkpt_ntrk_towards",50,0,50,0,30);
  h_leadingtrkpt_ntrk_transverse = fs->make<TProfile>("h_leadingtrkpt_ntrk_transverse","h_leadingtrkpt_ntrk_transverse",50,0,50,0,30);
  h_leadingtrkpt_ptsum_away = fs->make<TProfile>("h_leadingtrkpt_ptsum_away","h_leadingtrkpt_ptsum_away",50,0,50,0,30);
  h_leadingtrkpt_ptsum_towards = fs->make<TProfile>("h_leadingtrkpt_ptsum_towards","h_leadingtrkpt_ptsum_towards",50,0,50,0,30);
  h_leadingtrkpt_ptsum_transverse = fs->make<TProfile>("h_leadingtrkpt_ptsum_transverse","h_leadingtrkpt_ptsum_transverse",50,0,50,0,30);
  */
  bei->cd();
}//booking histos
//
// member functions
//

// ------------ method called for each event  ------------
void
FSQDQM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /*
  using namespace edm;
  using namespace std;
  using namespace reco;

  runNumber_   = iEvent.id().run();
  eventNumber_ = iEvent.id().event();
  lumiNumber_  = iEvent.id().luminosityBlock();
  bxNumber_ = iEvent.bunchCrossing();
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_Vtx_, vertices);
  //	 reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
  NPV->Fill(vertices->size());
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  
  for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
    if (!(vtx->chi2()==0 && vtx->ndof()==0)  &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0  && fabs(vtx->position().Z())<=24.0){
      bestvz = vtx->z(); 
      bestvx = vtx->x(); 
      bestvy = vtx->y();
      bestvzError = vtx->zError(); 
      bestvxError = vtx->xError(); 
      bestvyError = vtx->yError();
      PV_chi2->Fill(vtx->normalizedChi2());
      PV_d0->Fill(sqrt(vtx->x() * vtx->x() + vtx->y() * vtx->y()));
      PV_numTrks->Fill(vtx->tracksSize());
      double vertex_sumTrks = 0.0;
      for(reco::Vertex::trackRef_iterator iTrack= vtx->tracks_begin(); iTrack != vtx->tracks_end(); iTrack++)
	{
	  vertex_sumTrks += (*iTrack)->pt();
	}
      
      
      PV_sumTrks->Fill(vertex_sumTrks);
    }
  }
	 
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(tok_bs_, beamSpotH);
  std::vector<Jet> recoPFJets;
  recoPFJets.clear();
  
  int  nPFCHSJet=0;
  edm::Handle<PFJetCollection> pfjetchscoll;
  iEvent.getByToken(tok_pfjet_, pfjetchscoll);
  const reco::PFJetCollection *pfchsjets = pfjetchscoll.product();
  reco::PFJetCollection::const_iterator pfjetchsclus = pfchsjets->begin();
  for(pfjetchsclus = pfchsjets->begin(); pfjetchsclus!= pfchsjets->end() ; ++pfjetchsclus){
    PFJetpt->Fill( pfjetchsclus->pt());
    PFJeteta->Fill( pfjetchsclus->eta());
    PFJetphi->Fill( pfjetchsclus->phi());
    PFJetRapidity->Fill( pfjetchsclus->rapidity());
    nPFCHSJet++;
  }
  PFJetMulti->Fill( nPFCHSJet);

  
  std::vector<Jet> recoCastorJets;
  recoCastorJets.clear();
  
  
  edm::Handle<BasicJetCollection> castorJets;
  iEvent.getByToken(tok_castorjet_, castorJets);
  for (unsigned ijet=0; ijet<castorJets->size();ijet++) {
    recoCastorJets.push_back((*castorJets)[ijet]);
  }
  for (unsigned ijet=0; ijet<recoCastorJets.size(); ijet++) {
    cout<<recoCastorJets[ijet].pt()<<endl;
    CastorJetphi->Fill( recoCastorJets[ijet].phi());

    CastorJetMulti->Fill(recoCastorJets.size());
  }
  
	 edm::Handle<reco::TrackCollection> itracks;
	 iEvent.getByToken(tok_track_, itracks);
	 reco::TrackBase::TrackQuality hiPurity = reco::TrackBase::qualityByName("highPurity");
	 std::vector<TLorentzVector> T_trackRec_P4;
	 
	 int ntracks = 0;
	 int ntracks_towards = 0;
	 int ntracks_transverse = 0;
	 int ntracks_away = 0;
	 float ptsum = 0;
	 float dphi=0;
	 float ptsum_towards = 0;
	 float ptsum_transverse = 0;
	 float ptsum_away = 0;
	 
	 T_trackRec_P4.clear();
	 for(reco::TrackCollection::const_iterator iT = itracks->begin(); iT != itracks->end(); ++iT){
	     if(iT->quality(hiPurity)){
	       math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
	       double dzvtx = iT->dz(bestvtx);
	       double dxyvtx = iT->dxy(bestvtx);
	       double dzerror = sqrt(iT->dzError()*iT->dzError()+bestvzError*bestvzError);
	       double dxyerror = sqrt(iT->d0Error()*iT->d0Error()+bestvxError*bestvyError);
	       if((iT->ptError())/iT->pt() < 0.05 && dzvtx < 3.0 && dxyvtx < 3.0){
	       TLorentzVector trk; 
	       trk.SetPtEtaPhiE(iT->pt(),iT->eta(),iT->phi(),iT->p());
	       T_trackRec_P4.push_back(trk);
               Track_HP_Eta->Fill(iT->eta());
	       Track_HP_Phi->Fill(iT->phi());
	       Track_HP_Pt->Fill(iT->pt());
	       Track_HP_ptErr_over_pt->Fill((iT->ptError())/iT->pt());
	       Track_HP_dzvtx_over_dzerr->Fill(dzvtx/dzerror);
	       Track_HP_dxyvtx_over_dxyerror->Fill(dxyvtx/dxyerror); 
	       }
	     }
	 }
	 std::sort(T_trackRec_P4.begin(), T_trackRec_P4.end(), SortByPt());
	 for(unsigned int itrk=0;itrk<T_trackRec_P4.size();itrk++){
	     ++ntracks;
	     ptsum= ptsum + T_trackRec_P4.at(itrk).Pt();
	     dphi = deltaPhi(T_trackRec_P4.at(itrk).Phi(),T_trackRec_P4.at(0).Phi());
	     if(fabs(dphi) < 1.05){
	       ++ntracks_towards;
	       ptsum_towards = ptsum_towards + T_trackRec_P4.at(itrk).Pt();}
	     if(fabs(dphi) > 1.05 && fabs(dphi) < 2.09){
	       ++ntracks_transverse;
	       ptsum_transverse = ptsum_transverse + T_trackRec_P4.at(itrk).Pt();}
	     if(fabs(dphi) > 2.09){
	       ++ntracks_away; 
		   ptsum_away = ptsum_away + T_trackRec_P4.at(itrk).Pt();}
	   }       
	 
	 
	 
	 
	 h_ptsum_towards->Fill(ptsum_towards);
	 h_ptsum_transverse->Fill(ptsum_transverse);
	 h_ptsum_away->Fill(ptsum_away);
	 h_ntracks_towards->Fill(ntracks_towards);
	 h_ntracks_transverse->Fill(ntracks_transverse);
	 h_ntracks_away->Fill(ntracks_away);
	 if(T_trackRec_P4.size()>0){
	 h_leadingtrkpt_ntrk_towards->Fill(T_trackRec_P4.at(0).Pt(),ntracks_towards/8.37);
	 h_leadingtrkpt_ntrk_transverse->Fill(T_trackRec_P4.at(0).Pt(),ntracks_transverse/8.37);
	 h_leadingtrkpt_ntrk_away->Fill(T_trackRec_P4.at(0).Pt(),ntracks_away/8.37);
	 h_leadingtrkpt_ptsum_towards->Fill(T_trackRec_P4.at(0).Pt(),ptsum_towards/8.37);
	 h_leadingtrkpt_ptsum_transverse->Fill(T_trackRec_P4.at(0).Pt(),ptsum_transverse/8.37);
	 h_leadingtrkpt_ptsum_away->Fill(T_trackRec_P4.at(0).Pt(),ptsum_away/8.37);
	 }*/
}//analyze

void FSQDQM::endRun(edm::Run const& run, edm::EventSetup const& eSetup) {
  cout<<"Entering FSQDQM::endRun: "<<endl;

  // edm::LogVerbatim ("FSQDQM") <<"[FSQDQM]: End of Run, saving  DQM output
  // ";
  // int iRun = run.run();

  cout<<"...leaving FSQDQM::endRun. "<<endl;
}


//define this as a plug-in
//DEFINE_FWK_MODULE(FSQDQM);

//  LocalWords:  TH1F ptsum fs ntracks
