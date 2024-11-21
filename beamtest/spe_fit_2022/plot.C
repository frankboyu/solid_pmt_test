#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>
#include <TLegend.h>

using namespace std;

void plot(){
	gStyle->SetPalette(1);
	TFile *f8 = new TFile("beamtest_level1_4577_merged.root");
	TTree *tree_beam = (TTree*) f8->Get("T");
	const int rebinfac=4;
    double SPE[16]={296.6,250.3,243.2,198.6,248.4,271.2,278.9,269.7,87.8,88.8,111.5,114.8,290.6,311.4,262.22,307.2};//single p.e. value
    TString name_array[16]={"CerA0","CerA1","CerA2","CerA3","CerB0","CerB1","CerB2","CerB3","CerC0","CerC1","CerC2","CerC3","CerD0","CerD1","CerD2","CerD3"};
    TString Cer_array[16] = {"Cer[0]","Cer[1]","Cer[2]","Cer[3]","Cer[4]","Cer[5]","Cer[6]","Cer[7]","Cer[8]","Cer[9]","Cer[10]","Cer[11]","Cer[12]","Cer[13]","Cer[14]","Cer[15]"};  
	TH1F *hist_Cer[16]; 
    for(int t=0;t<16;t++)
	{
        tree_beam->Draw(Form("%s>>hist_Cer_%d(1500,0,1500)",Cer_array[t].Data(),t),"TrigType==8 && (Shower_t*0.88+Shower_l+Shower_r*0.91)>12000 && (PreSh_t*0.78+PreSh_l+PreSh_r*1.14)>2000","goff");
	 	hist_Cer[t] = (TH1F*)gROOT->FindObject(Form("hist_Cer_%d",t));

		ofstream outf(Form("20220001_%d000/Histo_Ch0%d.txt", t/4+1, t%4+1));
		for (int s=0; s<hist_Cer[t]->GetNbinsX(); s++)
		{
			outf << hist_Cer[t]->GetBinLowEdge(s+1) << "    " << hist_Cer[t]->GetBinContent(s+1) << endl;
		}
    }




}
