using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TFitResult.h>

int Xmax    = 1000;
int Xmin    = 0;
int Nbin    = 1000;
double Size =  1.0;

double parlo[10]      = {   50.0,      0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,   -20.0};
double parhi[10]      = {  100.0,     20.0,   7.0, 60.0,  1.0,  1.0,  1.0,  1.0, 12.0,    20.0};
double p0save[10]     = {   75.0,     10.0,   1.0, 20.0,  0.2,  0.2,  0.2,  0.5,  5.0,     0.0};
char par_names[10][7] = {"scale", "#sigma", "#mu", "v1", "a2", "c1", "c2", "c3", "#xi", "ped0"};

int Mmax     =  4;

const bool DEBUG = true;

double fitf(double *, double *);
double T(double, double, double, double, double, double, double);
double G(double, double, double);

double fit0(double *, double *);
double fit1(double *, double *);
double fit2(double *, double *);
double fit3(double *, double *);

int pmt_fit_adc(int date, int time, int readout)
{
	
    gStyle->SetOptStat(0);
	gStyle->SetOptFit(0112);    
    
    ifstream file(Form("/group/solid/www/solid/html/files/temp/%d_%d/Histo_Ch0%d.txt",date, time, readout));
	
	TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
	canvas->SetGrid(); canvas->SetLogy();

    TH1 *pmt_hist = new TH1D("pmt_hist", "pmt_hist", Nbin, Xmin, Xmax);
    int channel, count;
    while (file >> channel >> count) 
    {
        if (channel >= 1000)
            break;
        pmt_hist->SetBinContent(channel, count);

    }

    int xpedpk = pmt_hist->GetMaximumBin();
//         bool xpedpk_flag;
//         for(int bin_peak = pmt_hist->FindFirstBinAbove(); bin_peak <= 1000; bin_peak++)
//         {
//             xpedpk_flag = 1;
//             for(int ii = -10; ii <= 10; ii++)
//             {
//                 if (ii==0)
//                     continue;
//                 if (pmt_hist->GetBinContent(bin_peak+ii)<=0 || pmt_hist->GetBinContent(bin_peak+ii)>pmt_hist->GetBinContent(bin_peak))
//                 {	
//                     xpedpk_flag = 0;
//                     break;
//                 }
//             }
//             if(xpedpk_flag == 1)
//             {
//                 xpedpk = bin_peak;
//                 break;
//             }
//         }
    if(DEBUG) std::cout << "Pedestal peak position: " << xpedpk << std::endl;

    //CREATE THE NORMALIZED AND PEDESTAL SUBTRACTED HISTOGRAM
    Xmin = -100.0;
    Xmax = 500.0;
    Nbin = Xmax - Xmin;
    TH1D* hist_norm = new TH1D(Form("hist_norm_%d_%d_Ch0%d", date, time, readout), Form("hist_norm_%d_%d_Ch0%d", date, time, readout), Nbin, Xmin, Xmax );

    double norm_factor = 1. / ((double)pmt_hist->Integral());
    for(int bin = 1; bin <= hist_norm->GetXaxis()->GetNbins(); bin++)
    {
        double content      = (double)pmt_hist->GetBinContent(Xmin+bin+xpedpk);
        double norm_content = norm_factor * content;
        hist_norm->SetBinContent(bin, norm_content);

        if(content != 0)
        {
            hist_norm->SetBinError(bin, norm_factor*TMath::Sqrt(content));
        }
        else
        {
            hist_norm->SetBinError(bin, norm_factor);
        }
    }

    hist_norm->SetLineColor(kBlack);
    hist_norm->GetXaxis()->SetTitle( "s = adc - ped [adc channels]" );
    hist_norm->GetYaxis()->SetTitle( "dN/ds p.d.f. [a.u.]" );


    //INITIALZE THE FITTING FUNCTION
    TF1 *fit_pmt = new TF1("fit_pmt", fitf, -20., Xmax, 10);

    if( DEBUG ) std::cout << "Initial parameters:" << std::endl;
    for(int i = 0; i < 10; i++) 
    {
        fit_pmt->SetParameter(i, p0save[i]);
        fit_pmt->SetParLimits(i,  parlo[i],  parhi[i]);
        fit_pmt->SetParName(i, par_names[i]);
        if( DEBUG ) 
            cout << "\t" << par_names[i] << ": " << p0save[i] << endl;
    }
    fit_pmt->FixParameter(4, 0.0);
    fit_pmt->FixParameter(5, 0.0);
    fit_pmt->FixParameter(6, 0.0);
    fit_pmt->FixParameter(7, 0.0);

    double minx = -30;
    double maxx = 300;
    //double minx = Xmin;
    //double maxx = Xmax;    

    fit_pmt->SetRange(minx, maxx);

    //PERFORM THE FIT
    if(DEBUG) std::cout << "Fitting results: " << std::endl;
    TFitResultPtr result = hist_norm->Fit( fit_pmt, "SR0" );
    if(DEBUG) std::cout << " chi^2/NDF: " << result->Chi2() / result->Ndf() << std::endl;

    //DRAW THE FITTING RESULTS
    hist_norm->Draw();
    fit_pmt->SetNpx(1000);
    fit_pmt->Draw("same");

    TF1 *fd0 = new TF1( "fd0", fit0, minx, maxx, 10 );
    fd0->SetParameters( fit_pmt->GetParameters() );
    fd0->SetLineColor( kCyan );
    fd0->SetFillColor( kCyan );
    fd0->SetFillStyle( 3004 );
    fd0->Draw( "same" );

    TF1 *fd1 = new TF1( "fd1", fit1, minx, maxx, 10 );
    fd1->SetParameters( fit_pmt->GetParameters() );
    fd1->SetFillColor( kBlue );
    fd1->SetLineColor( kBlue );
    fd1->SetFillStyle( 3004 );
    fd1->Draw( "same" );

    TF1 *fd2 = new TF1( "fd2", fit2, minx, maxx, 10 );
    fd2->SetParameters( fit_pmt->GetParameters() );
    fd2->SetFillColor( kMagenta );
    fd2->SetLineColor( kMagenta );
    fd2->SetFillStyle( 3005 );
    fd2->Draw( "same" );

    TF1 *fd3 = new TF1( "fd3", fit3, minx, maxx, 10 );
    fd3->SetParameters( fit_pmt->GetParameters() );
    fd3->SetFillColor( kGreen );
    fd3->SetLineColor( kGreen );
    fd3->SetFillStyle( 3005 );
    fd3->Draw( "same" );

    canvas->SetLogy(true);
    canvas->Print(Form("/group/solid/www/solid/html/files/temp/%d_%d/Fit_Ch0%d.png", date, time, readout));

    //PRINT THE FITTING PARAMETERS
    ofstream outf(Form("/group/solid/www/solid/html/files/temp/%d_%d/Fit_Ch0%d.dat", date, time, readout));
    outf << "xpedpk" << "     " << xpedpk << "     " << "0.0" << endl;
    for( int i=0; i<10; i++ ) {
        double loc_par      = fit_pmt->GetParameter(i);
        double loc_parE     = fit_pmt->GetParError(i);
        string loc_par_name = par_names[i];
        outf << loc_par_name << "     " << loc_par << "     " << loc_parE << endl;
    }
    outf << "chi^2/NDF" << "     " << result->Chi2() / result->Ndf() << "     " << "0.0" << endl;
    outf.close();
    
	return 0;
	
}




double fitf(double *xx, double *par)
{
	double scale = par[0], sigma = par[1], mu = par[2], v1 = par[3], a2 = par[4], c1 = par[5], 
	       c2 = par[6], c3 = par[7], xi = par[8], ped = par[9];
	
	
	double a    = xx[0] / scale;
	double peda = ped   / scale;
	a          -= peda;
	
	
	double v2 = v1 * c1;
	double v3 = v1 * c3;
	double a3 = (1 - a2) * c2;
	double v = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;
	
	double A = 1;
 	double B = -1 * (2*v*a + 64/xi);
	double C = v*v * (a*a - 64*sigmaa*sigmaa);
	double D = B*B - 4*A*C;
	int nmax = 0;
	
	if(D>0)
	{
		nmax = TMath::Ceil((-1*B + TMath::Sqrt(D)) / (2*A));
		if(nmax<0)
		{
			nmax = 0;
			
		}
	}
	
	
	double result = 0;
	
	for(int n = 0; n < nmax; n++)
	{
		double seff = TMath::Sqrt(sigmaa*sigmaa + n/(v*v * xi));

		if(TMath::Abs(a - n/v) > 8*seff)
		{
			continue;
		}

		double q0 = 1;
		if(n!=0)
		{
			q0 = 0;
		}
		
		double pt = 0;

		for(int m = 1; m <= Mmax; m++)
		{
			pt += TMath::Poisson(m, mu) * T(n, m, v1, a2, v2, a3, v3);
		}
		result += G(a, n/v, seff) * (TMath::Exp(-1*mu) * q0 + pt);
		
	}
	return result / scale * Size;
}

double fit0(double *xx, double *par)
{
	double scale = par[0], sigma = par[1], mu = par[2], v1 = par[3], a2 = par[4], c1 = par[5], 
	       c2 = par[6], c3 = par[7], xi = par[8], ped = par[9];
	
	
	double a    = xx[0] / scale;
	double peda = ped   / scale;
	a          -= peda;
	
	
	double v2 = v1 * c1;
	double v3 = v1 * c3;
	double a3 = (1 - a2) * c2;
	double v = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;

	double A = 1;
 	double B = -1 * (2*v*a + 64/xi);
	double C = v*v * (a*a - 64*sigmaa*sigmaa);
	double D = B*B - 4*A*C;
	int nmax = 0;
	if(D>0)
	{
		nmax = TMath::Ceil((-1*B + TMath::Sqrt(D)) / (2*A));
		if(nmax<0)
		{
			nmax = 0;
			
		}
	}
	
	double result = 0;

	for(int n = 0; n < nmax; n++)
	{
		double seff = TMath::Sqrt(sigmaa*sigmaa + n/(v*v * xi));

		if(TMath::Abs(a - n/v) > 8*seff)
		{
			continue;
		}

		double q0 = 1;
		if(n!=0)
		{
			q0 = 0;
		}
		
		double pt = 0;
		
		int m = 0;
		pt += TMath::Poisson(m, mu) * T(n, m, v1, a2, v2, a3, v3);
		
		result += G(a, n/v, seff) * (TMath::Exp(-1*mu) * q0);
		
	}
	return result / scale * Size;
}

double fit1(double *xx, double *par)
{
	double scale = par[0], sigma = par[1], mu = par[2], v1 = par[3], a2 = par[4], c1 = par[5], 
	       c2 = par[6], c3 = par[7], xi = par[8], ped = par[9];
	
	
	double a    = xx[0] / scale;
	double peda = ped   / scale;
	a          -= peda;
	
	
	double v2 = v1 * c1;
	double v3 = v1 * c3;
	double a3 = (1 - a2) * c2;
	double v = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;

	double A = 1;
 	double B = -1 * (2*v*a + 64/xi);
	double C = v*v * (a*a - 64*sigmaa*sigmaa);
	double D = B*B - 4*A*C;
	int nmax = 0;
	if(D>0)
	{
		nmax = TMath::Ceil((-1*B + TMath::Sqrt(D)) / (2*A));
		if(nmax<0)
		{
			nmax = 0;
			
		}
	}
	
	double result = 0;

	for(int n = 0; n < nmax; n++)
	{
		double seff = TMath::Sqrt(sigmaa*sigmaa + n/(v*v * xi));
		
		if(TMath::Abs(a - n/v) > 8*seff)
		{
			continue;
		}

		double q0 = 1;
		if(n!=0)
		{
			q0 = 0;
		}
		
		double pt = 0;
		
		int m = 1;
		pt += TMath::Poisson(m, mu) * T(n, m, v1, a2, v2, a3, v3);
		
		result += G(a, n/v, seff) * (pt);
		
	}
	return result / scale * Size;
}

double fit2(double *xx, double *par)
{
	double scale = par[0], sigma = par[1], mu = par[2], v1 = par[3], a2 = par[4], c1 = par[5], 
	       c2 = par[6], c3 = par[7], xi = par[8], ped = par[9];
	
	
	double a    = xx[0] / scale;
	double peda = ped   / scale;
	a          -= peda;
	
	
	double v2 = v1 * c1;
	double v3 = v1 * c3;
	double a3 = (1 - a2) * c2;
	double v = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;

	double A = 1;
 	double B = -1 * (2*v*a + 64/xi);
	double C = v*v * (a*a - 64*sigmaa*sigmaa);
	double D = B*B - 4*A*C;
	int nmax = 0;
	if(D>0)
	{
		nmax = TMath::Ceil((-1*B + TMath::Sqrt(D)) / (2*A));
		if(nmax<0)
		{
			nmax = 0;
			
		}
	}
	
	double result = 0;

	for(int n = 0; n < nmax; n++)
	{
		double seff = TMath::Sqrt(sigmaa*sigmaa + n/(v*v * xi));

		if(TMath::Abs(a - n/v) > 8*seff)
		{
			continue;
		}

		double q0 = 1;
		if(n!=0)
		{
			q0 = 0;
		}
		
		double pt = 0;
		
		int m = 2;
		pt += TMath::Poisson(m, mu) * T(n, m, v1, a2, v2, a3, v3);
		
		result += G(a, n/v, seff) * (pt);
		
	}
	return result / scale * Size;
}

double fit3(double *xx, double *par)
{
	double scale = par[0], sigma = par[1], mu = par[2], v1 = par[3], a2 = par[4], c1 = par[5], 
	       c2 = par[6], c3 = par[7], xi = par[8], ped = par[9];
	
	
	double a    = xx[0] / scale;
	double peda = ped   / scale;
	a          -= peda;
	
	
	double v2 = v1 * c1;
	double v3 = v1 * c3;
	double a3 = (1 - a2) * c2;
	double v = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;

	double A = 1;
 	double B = -1 * (2*v*a + 64/xi);
	double C = v*v * (a*a - 64*sigmaa*sigmaa);
	double D = B*B - 4*A*C;
	int nmax = 0;
	if(D>0)
	{
		nmax = TMath::Ceil((-1*B + TMath::Sqrt(D)) / (2*A));
		if(nmax<0)
		{
			nmax = 0;
			
		}
	}
	
	double result = 0;

	for(int n = 0; n < nmax; n++)
	{
		double seff = TMath::Sqrt(sigmaa*sigmaa + n/(v*v * xi));

		if(TMath::Abs(a - n/v) > 8*seff)
		{
			continue;
		}

		double q0 = 1;
		if(n!=0)
		{
			q0 = 0;
		}
		
		double pt = 0;
		
		int m = 3;
		pt += TMath::Poisson(m, mu) * T(n, m, v1, a2, v2, a3, v3);
		
		result += G(a, n/v, seff) * (pt);
		
	}
	return result / scale * Size;
}

double T(double n, double m, double v1, double a2, double v2, double a3, double v3)
{
	double a1 = 1-a2-a3;

	double result = 0;
	for(double i1 = 0; i1 <= m; i1++)
	{
		for(double i2 = 0; i2 <= m-i1; i2++)
		{
			double i3 = m - i1 - i2;
			double vc = v1*i1 + v2*i2 + v3*i3;
			result += TMath::Factorial(m) / 
				(TMath::Factorial(i1)*TMath::Factorial(i2)*TMath::Factorial(i3))
				 *TMath::Power(a1, i1) * TMath::Power(a2, i2) * TMath::Power(a3, i3) 
				 *TMath::Poisson(n, vc);

		}
	}
	return result;
}

double G(double a, double nv, double s)//TMath::Gaus is giving a different answer than I am expecting
{
	return 1/(s *TMath::Sqrt(2*TMath::Pi())) * TMath::Exp(-1 * (a - nv)*(a-nv) /(2*s*s));
}
