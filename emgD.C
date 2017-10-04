#include "TH1.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>

R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"

using namespace std;

int idbg=0;
float ConeSize=0.4;
float D0SigCut=3;
float D0Cut=0.2;
float HTCUT = 1000.;
float PT1CUT = 400;
float PT2CUT = 200;
float PT3CUT = 200;
float PT4CUT = 100;
float JETETACUT = 2;
float ALPHAMAXCUT = 0.1;

std::ofstream myfile;



//------------------------------------------------------------------------------

float DeltaR(float eta1, float phi1, float eta2, float phi2) {

	float dR=0.;
	float deta = std::fabs(eta1-eta2);
	float dphi = std::fabs(phi1-phi2);
	if(dphi>3.14159) dphi = 2.*3.14159-dphi;
	dR=std::sqrt(deta*deta+dphi*dphi);

	return dR;
}

bool isQuark(int id){
	return abs(id) >= 1 && abs(id) <= 8;
}

struct Variables
{
	vector<string> names;
	int addVariable(string name){
		names.push_back(name);
		return names.size() - 1;
	}
	int fJetPT = addVariable("fJetPT");
	int fJetAM = addVariable("fJetAM");
	int fJetAMp = addVariable("fJetAMp");
	int fJetD0max = addVariable("fJetD0max");
	int fJetD0med = addVariable("fJetD0med");
	int fJetTHmed = addVariable("fJetTHmed");
	int fnJet = addVariable("fnJet");
	int fnTRK = addVariable("fnTRK");
	int ftrkPT = addVariable("ftrkPT");
	int ftrkTH = addVariable("ftrkTH");
	int ftrkD0 = addVariable("ftrkD0");
	int ftrkD0Error = addVariable("ftrkD0Error");
	int ftrkD0sig = addVariable("ftrkD0sig");
	int fMissingET = addVariable("fMissingET");
	int felectronPT = addVariable("felectronPT");
	int fmuonPT = addVariable("fmuonPT");
	int fHT = addVariable("fHT");
	int fdqd0 = addVariable("fdqd0");
	int fdd0 = addVariable("fdd0");

	int fhtnm1 = addVariable("fhtnm1");
	int fjpt1nm1 = addVariable("fjpt1nm1");
	int fjpt2nm1 = addVariable("fjpt2nm1");
	int fjpt3nm1 = addVariable("fjpt3nm1");
	int fjpt4nm1 = addVariable("fjpt4nm1");
	int famnm1 = addVariable("famnm1");

	int fnBJet = addVariable("fnBJet");
	int fBJetPT = addVariable("fBJetPT");

	int fptTop = addVariable("fptTop");
	int fptTopW = addVariable("fptTopW");

	//My additions
	int fangleWdecay = addVariable("fangleWdecay");
	int ptOfb = addVariable("ptOfb");
	int ptOfe = addVariable("ptOfe");
	int ptOfmu = addVariable("ptOfmu");
	int ptOftau = addVariable("ptOftau");

	int angleBvsWQuark = addVariable("angleBvsWQuark");
	int angleBvsWLepton = addVariable("angleBvsWLepton");

	int sumNeutrinoPt = addVariable("sumNeutrinoPt");

	int wAngleVsPt = addVariable("wAngleVsPt");

	int D0avgDist = addVariable("D0avgDist");

	int nTracksTopJets = addVariable("nTracksTopJets");
	int nTracksDarkJets = addVariable("nTracksDarkJets");

} variables;

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

TNtuple *AnalyseEvents(ExRootTreeReader *treeReader)
{
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchTRK = treeReader->UseBranch("Track");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
	TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");



	Long64_t allEntries = treeReader->GetEntries();

	cout << "** Chain contains " << allEntries << " events" << endl;

	GenParticle *prt;
	GenParticle *prt2;
	GenParticle *prtT;
	Track *trk;
	Jet *jet;
	MissingET *met;
	ScalarHT *ht;
	Electron *electron;
	Muon *muon;

	Long64_t entry;

	Int_t i;
	float dR;

	//Create NTuple
	std::stringstream ss;//Step 1: build up string that has variable names for the ntuple
	for(size_t i = 0; i < variables.names.size(); ++i) {
		if(i != 0) ss << ":";
		ss << variables.names[i];
	}
	TNtuple *ntuple = new TNtuple("variables","variables",ss.str().c_str());

	// Loop over all events
	
	int ijloop = allEntries;
	if(idbg>0) ijloop = 10;
	for(entry = 0; entry < ijloop; ++entry)
	{
		//Set up array to hold this row of the ntuple
		float values[variables.names.size()];


		if(idbg>0) myfile<<std::endl;
		if(idbg>0) myfile<<"event "<<entry<<std::endl;
		// Load selected branches with data from specified event
		treeReader->ReadEntry(entry);

		// Analyse gen particles
		int ngn = branchParticle->GetEntriesFast();
		int firstdq = -1;//These variables keep track of the first dark/down quarks found in the event. Not sure which variable is which yet.
		int firstadq = -1;
		int firstq = -1;
		int firstaq = -1;
		vector<int> pointtops;//gather particles relevant to study into here
		for(int i=0;i<ngn;i++ ) {//For the start, loop over particles and gather all relevant particles into pointtops

			prt = (GenParticle*) branchParticle->At(i);
			int id=(prt->PID);

			//look for W
			//find the initial daughters of the mediator
			if((id==4900101)&&(firstdq<0)) {
				firstdq=i;
				if(idbg>0) myfile<<" first dark quark"<<std::endl;
				firstq=i-1;
				prt2 = (GenParticle*) branchParticle->At(firstq);
			}
			if((id==-4900101)&&(firstadq<0)) {
				firstadq=i;
				if(idbg>0) myfile<<" first dark antiquark"<<std::endl;
				firstaq=i-1;
				prt2 = (GenParticle*) branchParticle->At(firstaq);
			}
			if(idbg>20) {
				myfile<<"genparticle "<<i<<" has pid "<<prt->PID<<" and pt of "<<prt->PT<<" status "<<prt->Status<<" mothers "<<prt->M1<<" "<<prt->M2<<std::endl;
			}

			//to help with background studies, find top and anti top
			if(abs(id)==6) {
				if(idbg>0) {
					std::cout<<" top at particle "<<i<<std::endl;
					std::cout<<" daugters are particles "<<prt->D1<<" "<<prt->D2<<std::endl;
				}
				prtT = (GenParticle*) branchParticle->At(prt->D1);
				int idpid1=abs(prtT->PID);
				if(idbg>0) std::cout<<"daughter 1 has pid "<<idpid1<<std::endl;
				prtT = (GenParticle*) branchParticle->At(prt->D2);
				int idpid2=abs(prtT->PID);
				if(idbg>0) std::cout<<"daughter 2 has pid "<<idpid2<<std::endl;

				//find the one that decays to W
				if((idpid1==24)||(idpid2==24) ) {
					pointtops.push_back(i);
					if(idbg>0) std::cout<<"choosing this top"<<std::endl;
				}
			}

		}//end looping over particles


		//make some plots about the tops in the events
		for(int i=0;i<pointtops.size();i++) {
			//get W
			GenParticle* w = NULL;
			GenParticle* b = NULL;


			prt = (GenParticle*) branchParticle->At(pointtops[i]);
			values[variables.fptTop] = prt->PT;
			prtT = (GenParticle*) branchParticle->At(prt->D1);

			//check first daughter particle
			if(abs(prtT->PID) == 24){//check for W
				values[variables.fptTopW] = prtT->PT;
				w = prtT;
			} else if (abs(prtT->PID) == 5){//check for b
				b = prtT;
				values[variables.ptOfb] = b->PT;
			}

			//check second daughter particle
			prtT = (GenParticle*) branchParticle->At(prt->D2);
			if(abs(prtT->PID)==24){//check for W
				values[variables.fptTopW] = prtT->PT;
				w = prtT;
			} else if (abs(prtT->PID) == 5){//check for b
				b = prtT;
				values[variables.ptOfb] = b->PT;
			}

			if (w){//if a W was found
				GenParticle *daughters[] = {(GenParticle*) branchParticle->At(w->D1), (GenParticle*)  branchParticle->At(w->D2)};
				bool d0IsQ = isQuark(daughters[0]->PID);
				bool d1IsQ = isQuark(daughters[1]->PID);
				bool d0IsLep = daughters[0]->PID == 11 || daughters[0]->PID == 13 || daughters[0]->PID == 15;
				bool d1IsLep = daughters[1]->PID == 11 || daughters[1]->PID == 13 || daughters[1]->PID == 15;
				bool d0IsNeutrino = daughters[0]->PID == 12 || daughters[0]->PID == 14 || daughters[0]->PID == 16;
				bool d1IsNeutrino = daughters[1]->PID == 12 || daughters[1]->PID == 14 || daughters[1]->PID == 16;

				if (isQuark(daughters[0]->PID) && isQuark(daughters[1]->PID)){
					float angle = DeltaR(daughters[0]->Eta, daughters[0]->Phi, daughters[1]->Eta, daughters[1]->Phi);
					values[variables.fangleWdecay] = angle;
					values[variables.wAngleVsPt] = w->PT, angle;
				}

				for (GenParticle *daughter: daughters){
					int pid = abs(daughter->PID);
					if (pid == 11){//if electron
						values[variables.ptOfe] = daughter->PT;
					} else if (pid == 13){
						values[variables.ptOfmu] = daughter->PT;
					} else if (pid == 15){
						values[variables.ptOftau] = daughter->PT;
					}
				}

				if (b){
					GenParticle *quark = d0IsQ? daughters[0]: d1IsQ ? daughters[1]: NULL;
					GenParticle *lepton = d0IsLep? daughters[0]: d1IsLep? daughters[1]: NULL;

					if (quark){
						values[variables.angleBvsWQuark] = DeltaR(quark->Eta, quark->Phi, b->Eta, b->Phi);
					}
					if (lepton){
						values[variables.angleBvsWLepton] = DeltaR(lepton->Eta, lepton->Phi, b->Eta, b->Phi);
					}
				}

				if (d0IsNeutrino || d1IsNeutrino){
					float x = 0, y = 0, z = 0;
					if (d0IsNeutrino){
						x += daughters[0]->Px;
						y += daughters[0]->Py;
						z += daughters[0]->Pz;
					}
					if (d1IsNeutrino){
						x += daughters[1]->Px;
						y += daughters[1]->Py;
						z += daughters[1]->Pz;
					}
					float magnitude = sqrt(x*x + y*y + z*z);
					values[variables.sumNeutrinoPt] = magnitude;
				}

			}
		}


		// Analyse tracks
		int ntrk = branchTRK->GetEntriesFast();
		vector<float> trkTheta(ntrk);
		values[variables.fnTRK] = ntrk;
		for(int i=0;i<ntrk;i++ ) {//Loop over all tracks
			trk = (Track*) branchTRK->At(i);
			// doing this at the generator level because I am too lazy to figure out the formulas
			// for the reconstructed
			// this would not be the right formula for pileup or if there
			// were a realistic vertex z distribution
			prt = (GenParticle*) trk->Particle.GetObject();
			float x1=prt->X;
			float y1=prt->Y;
			float z1=prt->Z;
			float px1=prt->Px;
			float py1=prt->Py;
			float pz1=prt->Pz;
			trkTheta[i]=0.;
			if((fabs(prt->X)>0.001)||(fabs(prt->Y)>0.001)) {//Calculate the theta? angle of each track
				float costt = (x1*px1+y1*py1+z1*pz1)/sqrt(x1*x1+y1*y1+z1*z1)/sqrt(px1*px1+py1*py1+pz1*pz1);
				trkTheta[i]=acos(costt);
			}
			values[variables.ftrkTH] = trkTheta[i];
			values[variables.ftrkPT] = trk->PT;
			values[variables.ftrkD0] = trk->D0;
			values[variables.ftrkD0Error] = fabs(trk->ErrorD0);  // for some reason, delphse pulls this from a caussian with a mean of zero, so half the time it is neg, which makes no sense to me
			//      std::cout<<"track d0 d0error "<<trk->D0<<" "<<trk->ErrorD0<<std::endl;
			if((trk->ErrorD0)>0) values[variables.ftrkD0sig] = fabs(trk->D0/(trk->ErrorD0));
		}


		// plots for jets and calculate displaced jet variables
		int njet = branchJet->GetEntriesFast();
		values[variables.fnJet] = njet;

		vector<float> alphaMax(njet);  // not really alpha max but best we can do here
		vector<float> alphaMaxp(njet);
		vector<float> D0Max(njet);
		vector<float> D0Med(njet);
		vector<float> THMed(njet);
		vector<int> ntrk1(njet);
		vector<bool> goodjet(njet);
		float allpT,cutpT,cutpTp;
		int ntrkj;
		vector<bool> adkq(njet);
		vector<bool> adq(njet);
		if(idbg>0) myfile<<" number of jets is "<<njet<<std::endl;
		int nbjets = 0;

		for(int i=0;i<njet;i++) {//Loop over jets
			jet = (Jet*) branchJet->At(i);

			bool isBJet = (jet->BTag>>0) & 0x1;
			// btag working points are accessed by bit-shifting
			// use 0 for loose, 1 for medium, and 2 for tight
			if (isBJet) {
				nbjets++;
				values[variables.fBJetPT] = jet->PT;
			}

			if(idbg>0) myfile<<"jet "<<i<<"  with pt, eta, phi of "<<jet->PT<<" "<<jet->Eta<<" "<<jet->Phi<<std::endl;
			values[variables.fJetPT] = jet->PT;
			adkq[i]=false;
			adq[i]=false;
			//see if it matches a dark or down quark
			if(firstdq>0) {
				prt2 = (GenParticle*) branchParticle->At(firstdq);
				float dr1=DeltaR(jet->Eta,jet->Phi,prt2->Eta,prt2->Phi);
				if(dr1<0.04) adkq[i]=true;
			}
			if(firstadq>0) {
				prt2 = (GenParticle*) branchParticle->At(firstadq);
				float dr1=DeltaR(jet->Eta,jet->Phi,prt2->Eta,prt2->Phi);
				if(dr1<0.04) adkq[i]=true;
			}
			if(firstq>0) {
				prt2 = (GenParticle*) branchParticle->At(firstq);
				float dr1=DeltaR(jet->Eta,jet->Phi,prt2->Eta,prt2->Phi);
				if(dr1<0.04) adq[i]=true;
			}
			if(firstaq>0) {
				prt2 = (GenParticle*) branchParticle->At(firstaq);
				float dr1=DeltaR(jet->Eta,jet->Phi,prt2->Eta,prt2->Phi);
				if(dr1<0.04) adq[i]=true;
			}

			// calculate track based variables
			alphaMax[i]=1.;
			alphaMaxp[i]=1.;
			goodjet[i]=false;
			D0Max[i]=0.;
			D0Med[i]=0.;
			THMed[i]=0.;
			allpT=0.;
			cutpT=0;
			cutpTp=0;
			ntrkj=0;

			for(int j=0;j<ntrk;j++) {
				trk = (Track*) branchTRK->At(j);
				dR=DeltaR(jet->Eta,jet->Phi,trk->Eta,trk->Phi);  //Find all tracks that are in the jet cone. Seems like this information should be directly acessable from jet...
				if(dR<ConeSize) {
					if((trk->PT)>1) {
						if(adkq[i]) {
							values[variables.fdqd0] = trk->D0;
						}
						if(adq[i]) {
							values[variables.fdd0] = trk->D0;
						}
						values[variables.D0avgDist] = trk->D0;
						ntrkj+=1;
						if((trk->D0)>D0Max[i]) D0Max[i]=(trk->D0);
						D0Med[i]=D0Med[i]+(trk->D0);
						THMed[i]=THMed[i]+trkTheta[j];
						allpT+=(trk->PT);
						if((fabs(trk->ErrorD0))>0) {  // this is not implemented by default.  Hope I did it right!
							if(fabs((trk->D0)/(trk->ErrorD0))<D0SigCut) {
								cutpT+=trk->PT;
							}}
						if(fabs((trk->D0))<D0Cut) {
							cutpTp+=trk->PT;
						}
						if(i<4) {
							if(idbg>3) myfile<<"   contains track "<<j<<" with pt, eta, phi of "<<trk->PT<<" "<<trk->Eta<<" "<<trk->Phi<<" d0 of "<<trk->D0<<
									//" and D0error of "<<trk->ErrorD0<<
									std::endl;
							prt = (GenParticle*) trk->Particle.GetObject();
							if(idbg>3) myfile<<"     which matches to get particle with XY of "<<prt->X<<" "<<prt->Y<<std::endl;

						}  // end first 4 jets
					}  //end pT cut of 1 GeV
				} //end in cone

			}  //end loop over tracks

			if(allpT>0) {
				alphaMax[i]=cutpT/allpT;
				alphaMaxp[i]=cutpTp/allpT;
			}
			if(alphaMax[i]>0.99999) alphaMax[i]=0.99999;
			if(alphaMaxp[i]>0.99999) alphaMaxp[i]=0.99999;

			ntrk1[i]=ntrkj;
			if(ntrkj>0) {
				D0Med[i]=D0Med[i]/ntrkj;
				THMed[i]=THMed[i]/ntrkj;
			}
			if((fabs(jet->Eta)<JETETACUT)&&(ntrk1[i]>0)) goodjet[i]=true;

			if(idbg>0) myfile<<"alpha max is "<<alphaMax[i]<<std::endl;
		} // end loop over all jets
		values[variables.fnBJet] = nbjets;


		// Analyse missing ET
		if(branchMissingET->GetEntriesFast() > 0)
		{
			met = (MissingET*) branchMissingET->At(0);
			values[variables.fMissingET] = met->MET;
		}


		// Analyse  HT
		if(branchScalarHT->GetEntriesFast() > 0)
		{
			ht = (ScalarHT*) branchScalarHT->At(0);
			values[variables.fHT] = ht->HT;
		}


		// Loop over all electrons in event
		for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
		{
			electron = (Electron*) branchElectron->At(i);
			values[variables.felectronPT] = electron->PT;
		}


		// Loop over all muons in event
		for(i = 0; i < branchMuon->GetEntriesFast(); ++i)
		{
			muon = (Muon*) branchMuon->At(i);
			values[variables.fmuonPT] = muon->PT;
		}



		//count number of the 4 leading jets with alpha max < a cut
		int nalpha=0;
		int iloop=min(4,njet);
		for(int i=0;i<iloop;i++) {
			values[variables.fJetAM] = alphaMax[i];
			values[variables.fJetAMp] = alphaMaxp[i];
			values[variables.fJetD0max] = D0Max[i];
			values[variables.fJetD0med] = D0Med[i];
			values[variables.fJetTHmed] = THMed[i];
			if(alphaMax[i]<ALPHAMAXCUT) {
				nalpha+=1;
				if(idbg>0) myfile<<" jet "<<i<<" passes alphamax cut with alphamax of "<<alphaMax[i]<<std::endl;
			}
		}

		// do pseudo emerging jets analysis

		// see if passes cuts
		bool Pnjet=false;
		bool Pht=false;
		bool Ppt1=false;
		bool Ppt2=false;
		bool Ppt3=false;
		bool Ppt4=false;
		bool Pam=false;
		if(njet>3) Pnjet=true;
		if(njet>3) {
			if((ht->HT)>HTCUT) Pht=true;
			jet = (Jet*) branchJet->At(0);
			if(((jet->PT)>PT1CUT)&&goodjet[0]) Ppt1=true;
			jet = (Jet*) branchJet->At(1);
			if(((jet->PT)>PT2CUT)&&goodjet[1]) Ppt2=true;
			jet = (Jet*) branchJet->At(2);
			if(((jet->PT)>PT3CUT)&&goodjet[2]) Ppt3=true;
			jet = (Jet*) branchJet->At(3);
			if(((jet->PT)>PT4CUT)&&goodjet[3]) Ppt4=true;
			if(nalpha>1) Pam=true;
		}


		//n-1 plots

		if(Pnjet&&Ppt1&&Ppt2&&Ppt3&&Ppt4&&Pam) values[variables.fhtnm1] = ht->HT;
		jet = (Jet*) branchJet->At(0);
		if(Pnjet&&Pht&&Ppt2&&Ppt3&&Ppt4&&Pam) values[variables.fjpt1nm1] = jet->PT;
		jet = (Jet*) branchJet->At(1);
		if(Pnjet&&Pht&&Ppt1&&Ppt3&&Ppt4&&Pam) values[variables.fjpt2nm1] = jet->PT;
		jet = (Jet*) branchJet->At(2);
		if(Pnjet&&Pht&&Ppt1&&Ppt2&&Ppt4&&Pam) values[variables.fjpt3nm1] = jet->PT;
		jet = (Jet*) branchJet->At(3);
		if(Pnjet&&Pht&&Ppt1&&Ppt2&&Ppt3&&Pam) values[variables.fjpt4nm1] = jet->PT;

		if(Pnjet&&Pht&&Ppt1&&Ppt2&&Ppt3&&Ppt4) {
			values[variables.famnm1] = alphaMax[0];
			values[variables.famnm1] = alphaMax[1];
			values[variables.famnm1] = alphaMax[2];
			values[variables.famnm1] = alphaMax[3];
		}
		

		//Now that we have collected all of our values for this event, we need to put the event into the ntuple
		ntuple->Fill(values);
	}
	return ntuple;
}

void emgD(const char *inputFile)
{
	gSystem->Load("libDelphes");

	TChain *chain = new TChain("Delphes");
	chain->Add(inputFile);

	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	//ExRootResult *result = new ExRootResult();

	//Data *plots = new Data;

	myfile.open("debug.txt");

	//BookHistograms(result, plots);
	
	TFile outFile("results.root","RECREATE");

	TNtuple *ntuple = AnalyseEvents(treeReader);
	ntuple->Write();
	outFile.Close();

	//plots->Count->LabelsDeflate();
	//plots->Count->LabelsOption("v");

//	//PrintHistograms(result, plots);

	//result->Write("results.root");

	myfile.close();

	cout << "** Exiting..." << endl;

	//delete plots;
	//delete result;
	delete treeReader;
	delete chain;
}

//------------------------------------------------------------------------------
