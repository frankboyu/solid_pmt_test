using namespace std;

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <vector>
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

int    num_paras     =   10;
int    m_max         =    5;
int    hist_pmt_min  =    0;
int    hist_pmt_max  = 1000;
int    hist_norm_min = -100;
int    hist_norm_max =  500;
int    fit_min       =  -50;
int    fit_max       =  500;
int    peak_low      =  200;
int    peak_high     =  400;

double parlo[10]      = {  -10.0,   10.0,      0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,    1.0};
double parhi[10]      = {   10.0,  150.0,     20.0,   5.0, 50.0,  1.0,  1.0,  1.0,  1.0,   30.0};
double p0save[10]     = {    0.0,   70.0,     10.0,   0.1, 10.0,  0.2,  0.2,  0.2,  0.2,    5.0};
char par_names[10][7] = { "ped0", "scale", "#sigma", "#mu", "v1", "a2", "c1", "c2", "c3", "#xi"};

bool DEBUG = true;

double fitf(double *, double *);
double fit0(double *, double *);
double fit1(double *, double *);
double fit2(double *, double *);
double fit3(double *, double *);
double T(double, double, double, double, double, double, double);
double G(double, double, double);

int pmt_fit_benchtest(int date, int time, int readout, int iteration)
{
    //SET UP THE CANVAS
    gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);    
	TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 1000);
	canvas->SetGrid(); 
	// canvas->SetLogy();

    //CREATE THE ORIGINAL HISTOGRAM
    TH1 *hist_pmt = new TH1D("hist_pmt", "hist_pmt", hist_pmt_max-hist_pmt_min, hist_pmt_min, hist_pmt_max);
    int channel, count;
    ifstream file(Form("/group/solid/www/solid/html/files/temp/%d_%d/Histo_Ch0%d.txt",date, time, readout));
    while (file >> channel >> count) 
    {
        if (channel >= hist_pmt_max)
            break;
        hist_pmt->SetBinContent(channel, count);
    }

	//FIND THE PEDESTAL PEAK POSITION
    int xpedpk = hist_pmt->GetMaximumBin();
    bool PEAK_FOUND;
    for(int peak_position = peak_low; peak_position <= peak_high; peak_position++)
    {
        PEAK_FOUND = 1;
        for(int ii = -20; ii <= 20; ii++)
        {
            if (ii==0)
                continue;
            if (hist_pmt->GetBinContent(peak_position+ii)<=0 || hist_pmt->GetBinContent(peak_position+ii)>hist_pmt->GetBinContent(peak_position))
            {	
                PEAK_FOUND = 0;
                break;
            }
        }
        if(PEAK_FOUND)
        {
            xpedpk = peak_position;
            break;
        }
    }
    if(DEBUG) std::cout << "Pedestal peak position: " << xpedpk << std::endl;

    //CREATE THE NORMALIZED AND PEDESTAL SUBTRACTED HISTOGRAM
    TH1D* hist_norm = new TH1D(Form("hist_norm_%d_%d_Ch0%d_Iter0%d", date, time, readout, iteration), Form("hist_norm_%d_%d_Ch0%d_Iter0%d", date, time, readout, iteration), hist_norm_max-hist_norm_min, hist_norm_min, hist_norm_max);

    double norm_factor = 1. / ((double)hist_pmt->Integral());
    for(int bin = 1; bin <= hist_norm->GetXaxis()->GetNbins(); bin++)
    {
        double content      = (double)hist_pmt->GetBinContent(hist_norm_min+bin+xpedpk);
        double norm_content = norm_factor * content;
        hist_norm->SetBinContent(bin, norm_content);

        if(content != 0) 
            hist_norm->SetBinError(bin, norm_factor*TMath::Sqrt(content));
        else  
            hist_norm->SetBinError(bin, norm_factor);
    }

    hist_norm->SetLineColor(kBlack);
    hist_norm->GetXaxis()->SetTitle( "s = adc - ped [adc channels]" );
    hist_norm->GetYaxis()->SetTitle( "dN/ds p.d.f. [a.u.]" );

	//GET THE INITIAL PARAMETERS
    if(iteration == 1)
	{
		TF1 *ped_fit = new TF1("ped_fit", "gaus", -3.0, 3.0);
		ped_fit->SetParameters(hist_norm->GetMaximum(), 0.0, 10.0);
		ped_fit->SetParLimits(2, 0.0, 20.0);
		
		hist_norm->Fit("ped_fit", "R0L");
		
		double ped0   = ped_fit->GetParameter(1);
		double sigma0 = ped_fit->GetParameter(2);
		double ypcal  = ped_fit->Integral(hist_norm_min, hist_norm_max);	
		double mu0    = -1. * TMath::Log(ypcal / hist_norm->Integral());
		double xave   = hist_norm->GetMean();
		double scale0 = xave/mu0;

		p0save[1] = scale0;
		p0save[2] = sigma0;
		p0save[3] = mu0;
	}
	else
	{
		ifstream file_paras(Form("/group/solid/www/solid/html/files/temp/%d_%d/Fit_Ch0%d.dat",date, time, readout));
		string line;
		while(getline(file_paras, line))
		{
			istringstream iss(line);
			double number, numbers[24];
			int index = 0;
			while (iss >> number) 
			{
				numbers[index] = number;
				index++;
			}
			for(int i = 0; i < 10; i++) 
			{
				p0save[i] = numbers[i+4];
			}
			file_paras.close();
		}
	}
	
	//INITIALZE THE FITTING FUNCTION
    TF1 *fit_pmt = new TF1("fit_pmt", fitf, fit_min, fit_max, 10);
	fit_pmt->SetRange(fit_min, fit_max);

    if( DEBUG ) std::cout << "Initial parameters:" << std::endl;
    for(int i = 0; i < 10; i++) 
    {
        fit_pmt->SetParameter(i, p0save[i]);
        fit_pmt->SetParLimits(i,  parlo[i],  parhi[i]);
        fit_pmt->SetParName(i, par_names[i]);
    }
	if(num_paras == 5)
	{
		fit_pmt->FixParameter(5, 0.0);
    	fit_pmt->FixParameter(6, 0.0);
    	fit_pmt->FixParameter(7, 0.0);
    	fit_pmt->FixParameter(8, 0.0);
    	fit_pmt->FixParameter(9, 1.0);
	}
	if(num_paras == 9)
	{
    	fit_pmt->FixParameter(9, 1.0);
	}
	if(iteration%2 == 0)
	{
		fit_pmt->FixParameter(4, fit_pmt->GetParameter(4));
		fit_pmt->FixParameter(5, fit_pmt->GetParameter(5));
		fit_pmt->FixParameter(6, fit_pmt->GetParameter(6));
		fit_pmt->FixParameter(7, fit_pmt->GetParameter(7));
		fit_pmt->FixParameter(8, fit_pmt->GetParameter(8));
		fit_pmt->FixParameter(9, fit_pmt->GetParameter(9));
	}
	if(DEBUG)
	{
		for(int i = 0; i < 10; i++) 
		{
			cout << "\t" << par_names[i] << ": " << fit_pmt->GetParameter(i) << endl;
		}
	}
    
    //PERFORM THE FIT
    if(DEBUG) std::cout << "Fitting results: " << std::endl;
    TFitResultPtr result = hist_norm->Fit( fit_pmt, "SR0" );
    if(DEBUG) std::cout << " chi^2/NDF: " << result->Chi2() / result->Ndf() << std::endl;

    //DRAW THE FITTING RESULTS
    hist_norm->Draw();
    fit_pmt->SetNpx(1000);
    fit_pmt->Draw("same");

    TF1 *fd0 = new TF1( "fd0", fit0, fit_min, fit_max, 10 );
    fd0->SetParameters( fit_pmt->GetParameters() );
    fd0->SetLineColor( kCyan );
    fd0->SetFillColor( kCyan );
    fd0->SetFillStyle( 3004 );
    fd0->Draw( "same" );

    TF1 *fd1 = new TF1( "fd1", fit1, fit_min, fit_max, 10 );
    fd1->SetParameters( fit_pmt->GetParameters() );
    fd1->SetFillColor( kBlue );
    fd1->SetLineColor( kBlue );
    fd1->SetFillStyle( 3004 );
    fd1->Draw( "same" );

    TF1 *fd2 = new TF1( "fd2", fit2, fit_min, fit_max, 10 );
    fd2->SetParameters( fit_pmt->GetParameters() );
    fd2->SetFillColor( kMagenta );
    fd2->SetLineColor( kMagenta );
    fd2->SetFillStyle( 3005 );
    fd2->Draw( "same" );

    TF1 *fd3 = new TF1( "fd3", fit3, fit_min, fit_max, 10 );
    fd3->SetParameters( fit_pmt->GetParameters() );
    fd3->SetFillColor( kGreen );
    fd3->SetLineColor( kGreen );
    fd3->SetFillStyle( 3005 );
    fd3->Draw( "same" );

    // canvas->SetLogy(true);
    canvas->Print(Form("/group/solid/www/solid/html/files/temp/%d_%d/Fit_Ch0%d.png", date, time, readout));

    //PRINT THE FITTING PARAMETERS
	ofstream outf(Form("/group/solid/www/solid/html/files/temp/%d_%d/Fit_Ch0%d.dat", date, time, readout));
	outf << date << "    " << time << "    " << readout << "    " << result->Chi2() / result->Ndf() << "    ";
	outf << xpedpk + fit_pmt->GetParameter(0) << "    ";
	for(int i=1; i<10; i++)    
	{
		outf << fit_pmt->GetParameter(i) << "    ";
	}
	for(int i=0; i<10; i++)    
	{
		outf << fit_pmt->GetParError(i) << "    ";
	}	
	outf << endl;
	outf.close();
    
	//END OF MAIN FUNCTION
    return 0;
}

double fitf(double *xx, double *par)
{
	double ped = par[0], scale = par[1], sigma = par[2], mu = par[3], v1 = par[4], a2 = par[5], c1 = par[6], 
	       c2 = par[7], c3 = par[8], xi = par[9];
	
	double peda   = ped   / scale;
	double a      = xx[0] / scale - peda;
	
	double v2     = v1 * c1;
	double v3     = v1 * c3;
	double a3     = (1 - a2) * c2;
	double v      = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;
	
	double A      = 1;
 	double B      = -1 * (2*v*a + 64/xi);
	double C      = v*v * (a*a - 64*sigmaa*sigmaa);
	double D      = B*B - 4*A*C;
	int    nmax   = 0;
	
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

		for(int m = 1; m <= m_max; m++)
		{
			pt += TMath::Poisson(m, mu) * T(n, m, v1, a2, v2, a3, v3);
		}
		result += G(a, n/v, seff) * (TMath::Exp(-1*mu) * q0 + pt);
		
	}
	return result / scale;
}

double fit0(double *xx, double *par)
{
	double ped = par[0], scale = par[1], sigma = par[2], mu = par[3], v1 = par[4], a2 = par[5], c1 = par[6], 
	       c2 = par[7], c3 = par[8], xi = par[9];

	double peda   = ped   / scale;
	double a      = xx[0] / scale - peda;
	
	double v2     = v1 * c1;
	double v3     = v1 * c3;
	double a3     = (1 - a2) * c2;
	double v      = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;
	
	double A      = 1;
 	double B      = -1 * (2*v*a + 64/xi);
	double C      = v*v * (a*a - 64*sigmaa*sigmaa);
	double D      = B*B - 4*A*C;
	int    nmax   = 0;

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
	return result / scale;
}

double fit1(double *xx, double *par)
{
	double ped = par[0], scale = par[1], sigma = par[2], mu = par[3], v1 = par[4], a2 = par[5], c1 = par[6], 
	       c2 = par[7], c3 = par[8], xi = par[9];
	
	double peda   = ped   / scale;
	double a      = xx[0] / scale - peda;
	
	double v2     = v1 * c1;
	double v3     = v1 * c3;
	double a3     = (1 - a2) * c2;
	double v      = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;
	
	double A      = 1;
 	double B      = -1 * (2*v*a + 64/xi);
	double C      = v*v * (a*a - 64*sigmaa*sigmaa);
	double D      = B*B - 4*A*C;
	int    nmax   = 0;

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
	return result / scale;
}

double fit2(double *xx, double *par)
{
	double ped = par[0], scale = par[1], sigma = par[2], mu = par[3], v1 = par[4], a2 = par[5], c1 = par[6], 
	       c2 = par[7], c3 = par[8], xi = par[9];
	
	double peda   = ped   / scale;
	double a      = xx[0] / scale - peda;
	
	double v2     = v1 * c1;
	double v3     = v1 * c3;
	double a3     = (1 - a2) * c2;
	double v      = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;
	
	double A      = 1;
 	double B      = -1 * (2*v*a + 64/xi);
	double C      = v*v * (a*a - 64*sigmaa*sigmaa);
	double D      = B*B - 4*A*C;
	int    nmax   = 0;	

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
	return result / scale;
}

double fit3(double *xx, double *par)
{
	double ped = par[0], scale = par[1], sigma = par[2], mu = par[3], v1 = par[4], a2 = par[5], c1 = par[6], 
	       c2 = par[7], c3 = par[8], xi = par[9];
	
	double peda   = ped   / scale;
	double a      = xx[0] / scale - peda;
	
	double v2     = v1 * c1;
	double v3     = v1 * c3;
	double a3     = (1 - a2) * c2;
	double v      = v1*(1-a2-a3) + v2*a2 + v3*a3;	
	double sigmaa = sigma/scale;
	
	double A      = 1;
 	double B      = -1 * (2*v*a + 64/xi);
	double C      = v*v * (a*a - 64*sigmaa*sigmaa);
	double D      = B*B - 4*A*C;
	int    nmax   = 0;	

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
	return result / scale;
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
