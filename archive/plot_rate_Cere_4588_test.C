
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

void plot_rate_Cere_4588_test(){
	gStyle->SetPalette(1);
	TFile *f8 = new TFile("beamtest_level1_merged.root");
	TTree *tree_beam = (TTree*) f8->Get("T");
	const int rebinfac=4;
        double SPE[16]={296.6,250.3,243.2,198.6,248.4,271.2,278.9,269.7,87.8,88.8,111.5,114.8,290.6,311.4,262.22,307.2};//single p.e. value
        TString name_array[16]={"CerA0","CerA1","CerA2","CerA3","CerB0","CerB1","CerB2","CerB3","CerC0","CerC1","CerC2","CerC3","CerD0","CerD1","CerD2","CerD3"};
        TString Cer_array[16] = {"Cer[0]","Cer[1]","Cer[2]","Cer[3]","Cer[4]","Cer[5]","Cer[6]","Cer[7]","Cer[8]","Cer[9]","Cer[10]","Cer[11]","Cer[12]","Cer[13]","Cer[14]","Cer[15]"};  
	TH1F *hist_Cer[16]; 
        for(int t=0;t<16;t++){
         tree_beam->Draw(Form("%s>>hist_Cer_%d(100,0,1500)",Cer_array[t].Data(),t),"TrigType==8 && (Shower_t*0.88+Shower_l+Shower_r*0.91)>12000 && (PreSh_t*0.78+PreSh_l+PreSh_r*1.14)>2000","goff");
	 hist_Cer[t] = (TH1F*)gROOT->FindObject(Form("hist_Cer_%d",t));
        }
	tree_beam->Draw("(Cer[0]/296.6+Cer[1]/250.3+Cer[2]/243.2+Cer[3]/198.6+Cer[4]/248.4+Cer[5]/271.2+Cer[6]/278.9+Cer[7]/269.7+Cer[8]/87.8+Cer[9]/88.8+Cer[10]/111.5+Cer[11]/114.8+Cer[12]/290.6+Cer[13]/311.4+Cer[14]/262.22+Cer[15]/307.2)>>hist_Cersum(100,0,50)","TrigType==8 && Shower_l>12000 && PreSh_l>3000","goff");
	TH1F *hist_Cer_sum = (TH1F*)gROOT->FindObject("hist_Cersum");
	TCanvas *c[2];
	c[0] = new TCanvas("c[0]","c[0]",1000,1000);
	c[0]->Divide(4,4);
	for(int j=0;j<16;j++){
		c[0]->cd(j+1);
		gPad->SetLogy();
		hist_Cer[j]->SetTitle(name_array[j]);
		hist_Cer[j]->Draw();
	}
	c[1] = new TCanvas("c[1]","c[1]",1000,1000);
	hist_Cer_sum->SetLineColor(4);
	hist_Cer_sum->GetXaxis()->SetTitle("Npe");
	hist_Cer_sum->Draw();
	double rate_hist_Cer_sum = hist_Cer_sum->Integral();
	//
	for(int a=0;a<2;a++){
		c[a]->SaveAs(Form("c%d.pdf",a));
	}
	gSystem->Exec("pdfunite ./c*.pdf ./CerABCD_run4588_highE_electron_left_cuts.pdf");
	gSystem->Exec(Form("rm -rf ./c*.pdf"));

}
