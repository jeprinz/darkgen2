
#include "tdrstyle.C"
#include "TH1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFrame.h"
#include <iostream>
#include <vector>
#include <map>

using namespace std;

class HistoOverlay{
public:
	vector<TH1*> histos;
	TCanvas *canv;
	TLegend *lgd;

	HistoOverlay(string name):
		histos(){
		setupCanvasAndLegend(name, 1);

	}
	~HistoOverlay(){
		delete lgd;
		delete canv;
	}
	const TH1F* addHisto(string name, string title, string key, TFile *file, float crossSection){
			TH1F *histo = static_cast<TH1F*>(file->Get(name.c_str())->Clone());
			histo->SetDirectory(0);
			double aaA = histo->Integral();
			std::cout<<" first entries is "<<aaA<<std::endl;
			histo->Scale(crossSection/aaA);

			float max = histo->GetMaximum();
			histo->SetMaximum(max*1.3);

			if (histos.size() == 0){//First histogram
				histo->GetYaxis()->SetTitle("percent ");
				histo->GetYaxis()->SetTitleSize(0.05);
				histo->GetXaxis()->SetTitle(title.c_str());
				histo->GetXaxis()->SetTitleSize(0.05);
			}

			histo->SetLineWidth(3);
			histo->SetStats(1);

			lgd->AddEntry(histo, key.c_str(), "l");

			if (histos.size() == 0){//For first histogram
				histo->Draw();
			} else { //For later ones
				histo->Draw("same");
			}
			histo->SetLineColor(histos.size() + 2);
			histos.push_back(histo);
			return histo;
	}

	void print(string fileName){

			float max = 0.0;
			for (TH1* histo: histos){
				if (histo->GetMaximum() > max) {
					max = histo->GetMaximum();
				}
			}

			if (histos.size()){
				histos[0]->SetMaximum(max*1.3);
			}

			gStyle->SetOptStat(0);
			lgd->Draw();



			canv->Update();
			canv->RedrawAxis();
			canv->GetFrame()->Draw();
			lgd->Draw();


			canv->Print((fileName+".png").c_str(),".png");
	}
private:
	void setupCanvasAndLegend(string name, bool dolog){
			int W = 800;
			int H = 700;
			// references for T, B, L, R
			float T = 0.08*H;
			float B = 0.22*H;
			float L = 0.12*W;
			float R = 0.04*W;

			canv = new TCanvas(name.c_str(),name.c_str(),50,50,W,H);

			canv->SetFillColor(0);
			canv->SetBorderMode(0);
			canv->SetFrameFillStyle(0);
			canv->SetFrameBorderMode(0);
			canv->SetLeftMargin( L/W );
			canv->SetRightMargin( R/W );
			canv->SetTopMargin( T/H );
			canv->SetBottomMargin( B/H );
			canv->SetTickx(0);
			canv->SetTicky(0);

			float x1_l = 1.2;
			float y1_l = 0.80;

			float dx_l = 0.60;
			float dy_l = 0.1;
			float x0_l = x1_l-dx_l;
			float y0_l = y1_l-dy_l;

			lgd = new TLegend(x0_l,y0_l,x1_l, y1_l);
			lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);


			if (dolog) canv->SetLogy();
	}
};

#define SIGNAL_CROSS 0.00615134 * 3
#define BACKGROUND_CROSS 888.0

#if 0
#define SIGNAL_CROSS 1
#define BACKGROUND_CROSS 1
#endif

int dolog=1;
void Overlay(string signalFileName, string backgroundFileName)
{ 

	vector<vector<string>> names = vector<vector<string>>({
//		vector<string>({"fangleWdecay", "angle between quarks from W decay"}),
//		vector<string>({"ptOfb", "Pt of b from top decay"}),
//		vector<string>({"ptOfe", "pt of electrons from W decay"}),
//		vector<string>({"ptOfmu", "pt of mu from W decay"}),
//		vector<string>({"ptOftau", "pt of tau from W decay"}),
//		vector<string>({"angleBvsWQuark", "angle between b and quark from W decay"}),
//		vector<string>({"angleBvsWLepton", "angle between b and lepton from W decay"}),
//		vector<string>({"sumNeutrinoPt", "magnitude of vector sum of neutrino pt from W decay"}),
//		vector<string>({"ptTop", "pt of top quarks"}),
//		vector<string>({"jet_D0max", "max track impact parameter in jets"}),
//		vector<string>({"jet_D0med", "med track impact parameter in jets"}),
		vector<string>({"D0avgDist", "average D0 distribution"}),
		vector<string>({"fdqd0","impact parameter for dark quark jets"}),
		vector<string>({"fdd0","impact parameter for down quark jets"})
	});

	map<string, TH1*> histos;

	TFile *signalFile = new TFile(signalFileName.c_str());
	TFile *backgroundFile = new TFile(backgroundFileName.c_str());

	for (vector<string> nameAndTitle: names){
		HistoOverlay histos(nameAndTitle[0]);
		histos.addHisto(nameAndTitle[0], nameAndTitle[1], "signal", signalFile, SIGNAL_CROSS);
		histos.addHisto(nameAndTitle[0], nameAndTitle[1], "background", backgroundFile, BACKGROUND_CROSS);
		histos.print(nameAndTitle[0]);

	}

	if(0){
		//Now make the extra graph
		HistoOverlay extra("neutrinos and leptons");
		extra.addHisto("ptOfe", "pt of electrons", "electrons signal", signalFile, SIGNAL_CROSS);
		extra.addHisto("ptOfe", "pt of electrons", "electrons background", backgroundFile, BACKGROUND_CROSS);
		extra.addHisto("sumNeutrinoPt", "pt of neutrinos", "neutrinos signal", signalFile, SIGNAL_CROSS);
		extra.addHisto("sumNeutrinoPt", "pt of neutrinos", "neutrinos background", backgroundFile, BACKGROUND_CROSS);
		extra.print("neutrinos-and-leptons-overlay");
	}
	if(0){
		//Overlay generator with detector level electrons
		HistoOverlay extra("generator vs detector level electrons");
		extra.addHisto("ptOfe", "pt of electrons", "e- gen signal", signalFile, SIGNAL_CROSS);
		extra.addHisto("ptOfe", "pt of electrons", "e- gen background", backgroundFile, BACKGROUND_CROSS);
		extra.addHisto("electronPT", "pt of electrons", "e- detector signal", signalFile, SIGNAL_CROSS);
		extra.addHisto("electronPT", "pt of electrons", "e- detector background", backgroundFile, BACKGROUND_CROSS);
		extra.print("generator-vs-detector-electrons");
	}
	if(1){
		double total, dark, down;
		{
			HistoOverlay dummyCounter("Impact parameter for various things");
			total = dummyCounter.addHisto("D0avgDist", "sum track impact parameter in all jets", "signal all med", signalFile, SIGNAL_CROSS)->GetEntries();
			dark = dummyCounter.addHisto("fdqd0", "impact parameter for dark quark tracks", "dark signal", signalFile, SIGNAL_CROSS)->GetEntries();
			down = dummyCounter.addHisto("fdd0", "impact parameter for down quark tracks", "down signal", signalFile, SIGNAL_CROSS)->GetEntries();
		}

		HistoOverlay extra("Impact parameter for various things");
		extra.addHisto("D0avgDist", "sum track impact parameter in all jets", "signal all tracks", signalFile, SIGNAL_CROSS);
		extra.addHisto("D0avgDist", "sum track impact parameter in all jets", "background all tracks", backgroundFile, BACKGROUND_CROSS);
		extra.addHisto("fdqd0", "impact parameter for dark quark tracks", "dark signal", signalFile, SIGNAL_CROSS * dark / total);
		extra.addHisto("fdd0", "impact parameter for down quark tracks", "down signal", signalFile, SIGNAL_CROSS * down / total);
		extra.print("impact-parameters");
	}
}

void Overlay(){
	string signalFileName = "signal/results.root";
	string backgroundFileName = "background/results.root";
	cout << "Using default arguments of signal file: " << signalFileName << " and background file: " << backgroundFileName << endl;
	Overlay(signalFileName, backgroundFileName);
}
