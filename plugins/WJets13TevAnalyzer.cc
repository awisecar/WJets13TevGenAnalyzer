// -*- C++ -*-
//
// Package:    WJets13TeV/WJets13TevAnalyzer
// Class:      WJets13TevAnalyzer
// 
/**\class WJets13TevAnalyzer WJets13TevAnalyzer.cc WJets13TeV/WJets13TevAnalyzer/plugins/WJets13TevAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Apichart Hortiangtham
//         Created:  Tue, 29 Nov 2016 13:48:53 GMT
//
//

#define DEBUG 0

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <TLorentzVector.h>

using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class WJets13TevAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit WJets13TevAnalyzer(const edm::ParameterSet&);
      ~WJets13TevAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
	edm::EDGetTokenT<reco::GenParticleCollection> theGenParLabel;
	//edm::EDGetTokenT<pat::PackedGenParticleCollection> theGenParLabel;
	edm::EDGetTokenT<std::vector<reco::GenJet>> theGenJetLabel;
	edm::EDGetTokenT<GenEventInfoProduct> theGEvntInfLabel;
	edm::EDGetTokenT<LHEEventProduct> theLHEEvtProdLabel;
	
	edm::Service<TFileService> fs;
	
	int nevent_run;
	//double jetPtCut;
	
	std::vector<TH1D *> hWeights;
	TH1D* _countEventSample;
	TH1D* _countEventBonzai;
	TH1D* _sumEventWeightSample;
	TH1D* _sumEventWeightBonzai;
	
	TH1D* _hist_excl_WJetMult;
	TH1D* _hist_excl_WJetMult_jetcut;
	TH1D* _hist_excl_WJetMult_vetomu;
	TH1D* _hist_excl_WJetMult_vetomujetcut;
	TH1D* _hist_NMuon;
	TH1D* _hist_NNeutrino;
	TH1D* _hist_NPhoton;
	TH1D* _hist_MuonPt_0j;
	TH1D* _hist_MuonEta_0j;
	TH1D* _hist_MT_0j;
	
	
	TTree *outputTree;
	
	vector<double> mcEventWeight_;
	vector<double> EvtWeights_;
	vector<double> pdfInfo_;
	
	vector<int> genLepId_;
	vector<int> genLepSt_;
	vector<int> genLepQ_;
	vector<float> genLepPt_;
	vector<float> genLepEta_;
	vector<float> genLepPhi_;
	vector<float> genLepE_;
	vector<int> genLepMomId_;
	vector<bool> genLepTauProd_;
	vector<bool> genLepPrompt_;
	
	vector<float> genPhoPt_;
	vector<float> genPhoEta_;
	vector<float> genPhoPhi_;
	vector<float> genPhoE_;
	vector<int> genPhoSt_;
	vector<bool> genPhoPrompt_;
	vector<int> genPhoMotherId_;
	
	vector<float> genJetPt_;
	vector<float> genJetEta_;
	vector<float> genJetPhi_;
	vector<float> genJetE_;
	//vector<double> genJetChF_;

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
WJets13TevAnalyzer::WJets13TevAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");
	
	theGenParLabel = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
	//theGenParLabel = consumes<pat::PackedGenParticleCollection>(edm::InputTag("packedGenParticles"));
	theGenJetLabel = consumes<std::vector<reco::GenJet>>(edm::InputTag("slimmedGenJets"));
	theGEvntInfLabel = consumes<GenEventInfoProduct>(edm::InputTag ("generator"));
	theLHEEvtProdLabel = consumes<LHEEventProduct>(edm::InputTag ("externalLHEProducer"));
	
	//jetPtCut = 30.0;
}


WJets13TevAnalyzer::~WJets13TevAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
WJets13TevAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	
	if (DEBUG) cout << __LINE__ << " NEW EVENT " << endl;
	using namespace edm;
	using namespace std;
	
	pdfInfo_.clear();
	mcEventWeight_.clear();
	EvtWeights_.clear();
	
	genLepId_.clear();
	genLepSt_.clear();
	genLepQ_.clear();
	genLepPt_.clear();
	genLepEta_.clear();
	genLepPhi_.clear();
	genLepE_.clear();
	genLepTauProd_.clear();
	genLepPrompt_.clear();
	genLepMomId_.clear();
	
	// -- gen photons ---
	genPhoPt_.clear();
	genPhoEta_.clear();
	genPhoPhi_.clear();
	genPhoE_.clear();
	genPhoSt_.clear();
	genPhoPrompt_.clear();
	genPhoMotherId_.clear();
	
	
	// -- gen jets ---
	genJetPt_.clear();
	genJetEta_.clear();
	genJetPhi_.clear();
	genJetE_.clear();
	//genJetChF_.clear();
	
	// ------------------------
	
	nevent_run++;
	
	edm::Handle<reco::GenParticleCollection> genParticles;
	//edm::Handle<pat::PackedGenParticleCollection> genParticles;
	iEvent.getByToken(theGenParLabel, genParticles);
	
	edm::Handle<reco::GenParticleCollection> genParticles_2;
	iEvent.getByToken(theGenParLabel, genParticles_2);
	
	edm::Handle<std::vector<reco::GenJet>> genjets;
	iEvent.getByToken(theGenJetLabel, genjets);
	
	edm::Handle<GenEventInfoProduct> genEventInfo;
	iEvent.getByToken(theGEvntInfLabel, genEventInfo);
	
	edm::Handle<LHEEventProduct> LHEEvtProd;
	iEvent.getByToken(theLHEEvtProdLabel, LHEEvtProd);
	
	// this has only one element, the central weight
	EvtWeights_ = genEventInfo->weights();
	// this also get the central weight
	double theWeight = genEventInfo->weight();
	
	int GMuCnt(0);
	int GMNuCnt(0);
	int GPhoCnt(0);
	
	double id1(0);
	double id2(0);
	double x1(0);
	double x2(0);
	double qscale(0);
	
	qscale = genEventInfo->qScale();
	if (genEventInfo->pdf()){
		x1  = genEventInfo->pdf()->x.first;
		x2  = genEventInfo->pdf()->x.second;
		id1 = genEventInfo->pdf()->id.first;
		id2 = genEventInfo->pdf()->id.second;
	}
	pdfInfo_.push_back((double) id1);
	pdfInfo_.push_back((double) id2);
	pdfInfo_.push_back((double) x1);
	pdfInfo_.push_back((double) x2);
	pdfInfo_.push_back((double) qscale);
	
	//---look at weight
	if (DEBUG) cout << " weight " << theWeight << " orig "  << LHEEvtProd->originalXWGTUP() << " size: " << LHEEvtProd->weights().size() /*<< " size: " << genEventInfo->weights().size() <<  " " << genEventInfo->weights()[0] */ << endl;
	
	//-------- check weights()[i].id
	//for (unsigned int i=0; i<LHEEvtProd->weights().size(); i++) {
	//    std::cout << " line: "<<  __LINE__ << " i " << i << " id " << LHEEvtProd->weights()[i].id << " wgt " << LHEEvtProd->weights()[i].wgt << std::endl;
	//}
	
	if (DEBUG) cout << __LINE__ << " mcEventWeight_.size() " <<  mcEventWeight_.size() << endl;
	//-------- Need the following to get sum of weight of all events
	mcEventWeight_.push_back((double) theWeight);
	for (unsigned int ih = 0 ; ih < 111; ih++){
		int pdfset = 1001;
		if (ih<9) {
			pdfset=1001+ih; //pdf sets for direct scale variation
		}
		else if (ih>=9  && ih< 111){
			pdfset=2000+ih-8;
		}
		else{
			cout << "ERROR: non existing pdfset" << endl;
		}
		std::string whichWeightId = std::to_string(pdfset);
		
		double weight = theWeight;
		for (unsigned int i=0; i<LHEEvtProd->weights().size(); i++) {
			if (LHEEvtProd->weights()[i].id == whichWeightId){
				weight *= LHEEvtProd->weights()[i].wgt/LHEEvtProd->originalXWGTUP();
				//if (DEBUG) std::cout << " i " << i << " id " << whichWeightId << " wgt " << LHEEvtProd->weights()[i].wgt << " new weight " << weight << std::endl;
				break;
			}
			//if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl;
		}
		hWeights[ih] ->Fill (0., weight);
		mcEventWeight_.push_back((double) weight);
	}
	if (DEBUG) cout << __LINE__ << " mcEventWeight_.size() " <<  mcEventWeight_.size() << endl;
	
	
	//--- GenParticles ---
	for (reco::GenParticleCollection::const_iterator p = genParticles->begin();p != genParticles->end(); ++p ) {
		
		int st = p->status();
		int id = p->pdgId();
		int motherId = p->numberOfMothers() ? p->mother(0)->pdgId() : -1 ;
		//int motherid = 0;
		//if (p->numberOfMothers()!=0) motherid = p->mother(0)->pdgId();
		//cout << " motherid " << motherid << endl;
		//if (id == -22) cout << " gen id " << id << endl;
		//if (st == -1 ) cout << " gen st " << st << endl;
		
		bool isPrompt = p->isPromptFinalState();
		bool isDirectTauDecay = p->isDirectPromptTauDecayProductFinalState();
		
		//----- lepton
		if (st==1 && (fabs(id)==11 || fabs(id)==13 || fabs(id)==15
					  || fabs(id)==12|| fabs(id)==14|| fabs(id)==16)){
			
			bool isFromTau(0);
			const reco::Candidate *tau_mother = 0;
			if(p->numberOfMothers()!=0) {
				tau_mother= p->mother(0);
				//cout << " mother is " << tau_mother->pdgId() << endl;
				if(fabs(tau_mother->pdgId()) == 15) isFromTau = true;
			}
			
			TLorentzVector genLep1(0,0,0,0);
			genLep1.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
			double charge(0);
			
			if ( fabs(id)==12|| fabs(id)==14|| fabs(id)==16) charge = 0;
			else if (id > 0) charge = -1;
			else if (id < 0) charge = 1;
			
			if (DEBUG) cout << " Lep id " << id << " " << fabs(id) << " pt " << genLep1.Pt() << " eta " << genLep1.Eta() << " isPromptFinalState " << isPrompt << " isFromTau " << isFromTau << " isDirectPromptTauDecayProductFinalState " << isDirectTauDecay << endl;
			
			//---- dressed lepton
			TLorentzVector genR1DressLep1(0,0,0,0);
			TLorentzVector genR1Pho1(0,0,0,0);
			for (reco::GenParticleCollection::const_iterator p_2 = genParticles_2->begin(); p_2 != genParticles_2->end(); ++p_2 ){
				if ( !(p_2->numberOfMothers() && p_2->status() == 1 && p_2->pdgId() == 22 && p_2->energy() >= 0.000001) ) continue;
				TLorentzVector thisPho1(0,0,0,0);
				thisPho1.SetPtEtaPhiE(p_2->pt(),p_2->eta(),p_2->phi(),p_2->energy());
				double dR = genLep1.DeltaR(thisPho1);
				if(dR<0.1){
					genR1Pho1+=thisPho1;
				}
			}
			genR1DressLep1 = genLep1+genR1Pho1;
			
			if (DEBUG){
				if (genLep1.Pt() < 15 && genR1DressLep1.Pt() > 25 && fabs(genR1DressLep1.Eta()) < 2.4 && fabs(id)==13){
					cout << " This event has been missing: " << " genR1Pho1.Pt() " << genR1Pho1.Pt() << " genLep1.Pt() " << genLep1.Pt() << " genR1DressLep1.Pt() " << genR1DressLep1.Pt() << " genR1DressLep1.Eta() " << genR1DressLep1.Eta() << endl;
				}
			}
			
			//if (genLep1.Pt() > 15 && fabs(id)==13) GMuCnt++;
			if (genR1DressLep1.Pt() > 15 && fabs(id)==13) GMuCnt++;
			if (fabs(id)==14) GMNuCnt++;
			
			genLepId_.push_back(id);
			genLepSt_.push_back(st);
			genLepQ_.push_back(charge);
			genLepPt_.push_back(genLep1.Pt());
			genLepEta_.push_back(genLep1.Eta());
			genLepPhi_.push_back(genLep1.Phi());
			genLepE_.push_back(genLep1.Energy());
			genLepTauProd_.push_back(isFromTau);
			genLepPrompt_.push_back(isPrompt);
			genLepMomId_.push_back(motherId);
		}
		
		//----- photon
		if (st==1 && fabs(id) == 22){
			
			TLorentzVector genPho(0, 0, 0, 0);
			genPho.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
			
			if (DEBUG) cout << " Pho id " << id << " " << fabs(id) << " pt " << genPho.Pt() << " eta " << genPho.Eta() << " isPromptFinalState " << isPrompt << " isFromTau " << "NAN" << " isDirectPromptTauDecayProductFinalState " << isDirectTauDecay << endl;
			
			// count PromptFinalState photon ET > 10GeV
			if (genPho.Et() > 10 && isPrompt) GPhoCnt++;
			
			genPhoPt_.push_back(genPho.Pt());
			genPhoEta_.push_back(genPho.Eta());
			genPhoPhi_.push_back(genPho.Phi());
			genPhoE_.push_back(genPho.Energy());
			genPhoSt_.push_back(st);
			genPhoPrompt_.push_back(isPrompt);
			genPhoMotherId_.push_back(motherId);
			
		}
		
	} // end loop over GenParticles
	
	if (DEBUG) cout << "  GMuCnt : " << GMuCnt << " GMNuCnt " << GMNuCnt << endl;
	if (DEBUG) cout << "  GPhoCnt : " << GPhoCnt << endl;
	
	
	//--- GenJets ---
	int ngenjet=0;
	int ngenjetNoCut=0;
	for (reco::GenJetCollection::const_iterator jet = genjets->begin(); jet!=genjets->end(); jet++){
		TLorentzVector thisJet(0,0,0,0);
		thisJet.SetPtEtaPhiE(jet->pt(),jet->eta(),jet->phi(),jet->energy());
		if (thisJet.Pt() < 10.) continue;
		
		genJetPt_.push_back(thisJet.Pt());
		genJetEta_.push_back(thisJet.Eta());
		genJetPhi_.push_back(thisJet.Phi());
		genJetE_.push_back(thisJet.Energy());
		
		ngenjetNoCut++;
		if (jet->pt()< 30) continue;
		if (fabs(jet->rapidity())>2.4) continue;
		ngenjet++;
	}
	//if (genjets.failedToGet()) cout << " genjets.failedToGet() " << genjets.failedToGet() << " ngenjetNoCut " << ngenjetNoCut << endl;
	if (DEBUG) cout << " ngenjetNoCut " << ngenjetNoCut << " size " << genJetPt_.size() << endl;
	
	
	
	//------- Filling Histograms --------------
	_hist_excl_WJetMult->Fill(ngenjetNoCut, theWeight);
	_hist_excl_WJetMult_jetcut->Fill(ngenjet, theWeight);
	_hist_NMuon->Fill(GMuCnt, theWeight);
	_hist_NNeutrino->Fill(GMNuCnt, theWeight);
	_hist_NPhoton->Fill(GPhoCnt, theWeight);
	//    _hist_MuonPt_0j->Fill(muDressed.Pt(), theWeight);
	//    _hist_MuonEta_0j->Fill(muDressed.Pt(), theWeight);
	//    _hist_MT_0j->Fill(WmT, theWeight);
	_countEventSample->Fill(0.);
	_sumEventWeightSample->Fill(0., theWeight);

	
	//---- if you want to veto events regarding GMuCnt, do it here
	//if ( ! ( GMuCnt > 0 && GMNuCnt > 0) ) return;
	if (GMuCnt < 1) return;
	//-----
	
	_hist_excl_WJetMult_vetomu->Fill(ngenjetNoCut, theWeight);
	_hist_excl_WJetMult_vetomujetcut->Fill(ngenjet, theWeight);
	
	
	
	//----- sum of weight after veto events
	for (unsigned int ih = 0 ; ih < 111; ih++){
		int pdfset = 1001;
		if (ih<9) {
			pdfset=1001+ih; //pdf sets for direct scale variation
		}
		else if (ih>=9  && ih< 111){
			pdfset=2000+ih-8;
		}
		else{
			cout << "ERROR: non existing pdfset" << endl;
		}
		std::string whichWeightId = std::to_string(pdfset);
		
		double weight = theWeight;
		for (unsigned int i=0; i<LHEEvtProd->weights().size(); i++) {
			if (LHEEvtProd->weights()[i].id == whichWeightId){
				weight *= LHEEvtProd->weights()[i].wgt/LHEEvtProd->originalXWGTUP();
				//if (DEBUG) std::cout << " i " << i << " id " << whichWeightId << " wgt " << LHEEvtProd->weights()[i].wgt << " new weight " << weight << std::endl;
				break;
			}
			//if (DEBUG) std::cout << __LINE__ << " " << theWeight << " " << evtWeights.back() << std::endl;
		}
		hWeights[ih] ->Fill (1., weight);
	}
	_countEventBonzai->Fill(0.);
	_sumEventWeightBonzai->Fill(0., theWeight);
	
	//----- save tree
	outputTree->Fill();
 
	
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
	
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
WJets13TevAnalyzer::beginJob()
{
	_countEventSample = fs->make<TH1D>("countEventSample", "countEventSample", 2, -0.5, 1.5);
	_countEventBonzai = fs->make<TH1D>("countEventBonzai", "countEventBonzai", 2, -0.5, 1.5);
	_sumEventWeightSample = fs->make<TH1D>("sumEventWeightSample", "sumEventWeightSample", 2, -0.5, 1.5);
	_sumEventWeightBonzai = fs->make<TH1D>("sumEventWeightBonzai", "sumEventWeightBonzai", 2, -0.5, 1.5);

	_hist_excl_WJetMult = fs->make<TH1D>("njetWJet_excl", "njetWJet_excl", 8, -0.5, 7.5);
	_hist_excl_WJetMult_jetcut = fs->make<TH1D>("njetWJet_excl_jetcut", "njetWJet_excl_jetcut", 8, -0.5, 7.5);
	_hist_excl_WJetMult_vetomu = fs->make<TH1D>("njetWJet_excl_vetomu", "njetWJet_excl_vetomu", 8, -0.5, 7.5);
	_hist_excl_WJetMult_vetomujetcut = fs->make<TH1D>("njetWJet_excl_vetomujetcut", "njetWJet_excl_vetomujetcut", 8, -0.5, 7.5);
	_hist_NMuon     = fs->make<TH1D>("NMuon", "NMuon", 8, -0.5, 7.5);
	_hist_NNeutrino = fs->make<TH1D>("NNeutrino", "NNeutrino", 8, -0.5, 7.5);
	_hist_NPhoton   = fs->make<TH1D>("NPhoton", "NPhoton", 8, -0.5, 7.5);
	_hist_MuonPt_0j     = fs->make<TH1D>("muPt_inc0jet", "muPt_inc0jet", 40, 0, 200);
	_hist_MuonEta_0j    = fs->make<TH1D>("muEta_inc0jet", "muEta_inc0jet", 24,-2.4, 2.4);
	_hist_MT_0j         = fs->make<TH1D>("MT_inc0jet", "MT_inc0jet", 200,0.,400);
	
	outputTree = new TTree("EventTree"," EventTree");
	
	outputTree->Branch("pdfInfo_", &pdfInfo_);
	outputTree->Branch("mcEventWeight_", &mcEventWeight_);
	outputTree->Branch("EvtWeights", &EvtWeights_);
	
	
	outputTree->Branch("GLepBareId", &genLepId_);
	outputTree->Branch("GLepBareSt", &genLepSt_);
	outputTree->Branch("GLepBareQ", &genLepQ_);
	outputTree->Branch("GLepBarePt", &genLepPt_);
	outputTree->Branch("GLepBareEta", &genLepEta_);
	outputTree->Branch("GLepBarePhi", &genLepPhi_);
	outputTree->Branch("GLepBareE", &genLepE_);
	outputTree->Branch("GLepBareTauProd", &genLepTauProd_);
	outputTree->Branch("GLepBarePrompt", &genLepPrompt_);
	outputTree->Branch("GLepBareMomId", &genLepMomId_);
	
	outputTree->Branch("GPhotPt", &genPhoPt_);
	outputTree->Branch("GPhotEta", &genPhoEta_);
	outputTree->Branch("GPhotPhi", &genPhoPhi_);
	outputTree->Branch("GPhotE", &genPhoE_);
	outputTree->Branch("GPhotSt", &genPhoSt_);
	outputTree->Branch("GPhotPrompt", &genPhoPrompt_);
	outputTree->Branch("GPhotMotherId", &genPhoMotherId_);
	
	outputTree->Branch("GJetAk04Pt", &genJetPt_ );
	outputTree->Branch("GJetAk04Eta", &genJetEta_);
	outputTree->Branch("GJetAk04Phi", &genJetPhi_);
	outputTree->Branch("GJetAk04E", &genJetE_ );
	//outputTree->Branch("genJetChF_", &genJetChF_ );
	
	std::stringstream name;
	TH1D* histo = NULL;
	for (unsigned int ih = 0; ih < 9 ; ih++) {
		name.str("");
		name<< "hWeights_"<<1001+ih;
		histo=fs->make<TH1D> (name.str().c_str(),  name.str().c_str(), 3, -0.5, 2.5);
		hWeights.push_back(histo);
		name.str("");
	}
	for (unsigned int ih = 1; ih < 103 ; ih++) {
		name.str("");
		name<< "hWeights_"<<2000+ih;
		histo=fs->make<TH1D> (name.str().c_str(),  name.str().c_str(), 3, -0.5, 2.5);
		hWeights.push_back(histo);
		name.str("");
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void
WJets13TevAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
WJets13TevAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WJets13TevAnalyzer);
