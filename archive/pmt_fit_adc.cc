
#include "fit.h"


int Xmax    = 3500;
int Xmin    = -500;
int Nbin    = 4000;
double Size =  1.0;

double parlo[9]  = {  0.0,  0.0,   0.0,   7.5,  0.1,  0.1,  0.1,  0.4,  2.5};
double parhi[9]  = {750.0,  5.0,   0.8,  40.0,  0.5,  0.4,  0.9,  0.8, 12.0};
double p0save[9] = {  0.0,  1.5, 0.007,  20.0,  0.2,  0.2,  0.2,  0.5,  5.0};

double nfpmt =  0;
int Mmax     =  4;

char par_names[9][7] = {"scale", "#sigma", "#mu", "v1", "a2", "c1", "c2", "c3", "#xi"};

const bool DEBUG = true;


int main(int argc, char **arg) 
{
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	
	if(argc != 3) return 0;
	
	string pmt  = "HA0053";
	int anode   = stoi(argv[1);
	int hv      = stoi(argv[2]);
	
	if(anode < 1 || anode > 64) return 0;
	
	//---------   Read in the file   ----------//
	
	//TFile *file = new TFile("inFile.root", "READ");
	TFile *file = new TFile(Form("hv_data/%dV/run_000001.spectra.root",hv), "READ");
	
	TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
	canvas->SetGrid(); canvas->SetLogy();
	
	{
		if(DEBUG) cout << "Channel: " << anode << endl;
		
		TH1 *pmt_hist = new TH1D(Form("pmt_hist_%d", anode), Form("pmt_hist_%d", anode), 
			Nbin, Xmin, Xmax);
		
		file->GetObject(Form("pix%02d_adc", anode), pmt_hist);
		
		int xpedpk = pmt_hist->GetMaximumBin();
		
		Xmin = pmt_hist->FindFirstBinAbove() - xpedpk;
		Xmax = pmt_hist->FindLastBinAbove()  - xpedpk;
		if( (Xmax-Xmin)%2 ) Xmax += 1;
		Nbin = Xmax - Xmin;
		
		double norm_factor = 1. / ((double)pmt_hist->Integral());
		
		if(DEBUG) {
			std::cout << "min: " << Xmin << ", max: " << Xmax << std::endl;
			std::cout << "integral: " << (pmt_hist->Integral()) << endl;
		}
		
		TH1D* hist_norm = new TH1D( Form("hist_norm_%d", anode), 
			Form("hist_norm_%d", anode), Nbin, Xmin, Xmax );
		
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
		
		//hist_norm->Rebin(2);
		//hist_norm->Scale(0.5);
		hist_norm->SetLineColor(kBlack);
		hist_norm->SetTitle( Form("MaPMT HA00%02d Anode #%d", pmt.c_str(), anode) );
		hist_norm->GetXaxis()->SetTitle( "s = adc - ped [adc channels]" );
		hist_norm->GetYaxis()->SetTitle( "dN/ds p.d.f. [a.u.]" );
		
		hist_norm->Draw();
		
		if(DEBUG) std::cout << "Anode: " << anode << " xpedpk: " << xpedpk << std::endl;
		
		//-----  Gaussian fit to pedestal  -----//
		
		double ped_min = hist_norm->GetBinCenter(hist_norm->GetMaximumBin()-3);
		double ped_max = hist_norm->GetBinCenter(hist_norm->GetMaximumBin()+3);
		
		TF1 *ped_fit = new TF1( "ped_fit", "gaus", ped_min, ped_max );
		ped_fit->SetParameters( hist_norm->GetMaximum(), 0., 0.002 );
		ped_fit->SetParLimits( 2, 0., 0.1 );
		
		hist_norm->Fit("ped_fit", "R0L");
		
		double ped0   =  ped_fit->GetParameter(1);
		double sigma0 =  ped_fit->GetParameter(2);
		if(DEBUG) printf("Pedestal = %.4f; Width = %.4f", ped0, sigma0);
		
		ped_fit->Draw("same");
		
		// Main fit inital Parameters
		
		double scale =  p0save[0];
		double sigma =  p0save[1];
		double mu    =  p0save[2];
		double v1    =  p0save[3];
		double a2    =  p0save[4];
		double c1    =  p0save[5];
		double c2    =  p0save[6];
		double c3    =  p0save[7];
		double xi    =  p0save[8];
		
		// Guess Some Parameters:
		
		double ypcal  = ped_fit->Integral( Xmin, Xmax ) / Size;	
		double mu0    = -1. * TMath::Log(ypcal / hist_norm->Integral());
		double xave   = hist_norm->GetMean();
		double scale0 = xave/mu0;
		
		
		if(DEBUG) {
			cout << "Integral, ypcal, mu0, xave, scale0: " << hist_norm->Integral() 
			     << ", "<< ypcal << ", " << mu0 << ", " << xave 
			     << ", " << scale0 << std::endl << endl;
		}
		
		double nguess[6] = {9, 4, 4, 4, 4, 1};
		double step[6] = {0, 0, 0, 0, 0, 0};
		double ngt = 1;
		double ngtt[6] = {1, 1, 1, 1, 1, 1};
		
		/*
		for(int i = 0; i < 6; i++) {
			if(nguess[i] != 0) {
				step[i] = (parhi[i+3] - parlo[i+3]) / (nguess[i] + 1);
			}
			else {
				step[i] = 0;
			}
			ngt = ngt*nguess[i];
			ngtt[i+1]   = ngtt[i]*nguess[i];
			p0save[i+3] = parlo[i+3];
		}
		*/
		
		// initial guess for the 'scale' based on the hv setting:
		
		if(hv == 950) {
			if( scale0 < 0. )        scale0 =  50.;
			else if( scale0 > 200. ) scale0 = 200.;
			
			if(mu0 < 0.005)      mu0 = 0.005;
			else if(mu0 > 0.025) mu0 = 0.025;
		}
		if(hv == 1000) {
			if( scale0 < 0. )        scale0 = 100.;
			else if( scale0 > 500. ) scale0 = 500.;
			
			if( mu0 < 0.005 )      mu0 = 0.005;
			else if( mu0 > 0.025 ) mu0 = 0.025;
		}
		
		
		p0save[0] = scale0;
		p0save[1] = sigma0;
		p0save[2] = mu0;
		
		double chisqrmin = 1.e100;
		
		double p0[9] = { scale0, sigma0, mu0, 0, 0, 0, 0, 0, 2.5 };
		double ycal = 0;
		
		
		TF1 *fit_pmt = new TF1("fit_pmt", fitf, -20., Xmax, 10);
		
		fit_pmt->SetParameter(0,scale0);
		fit_pmt->SetParameter(1,sigma0);
		fit_pmt->SetParameter(2,mu0);
		
		if( DEBUG ) std::cout << "Parameters:" << std::endl;
		for(int i = 0; i<3; i++) {
			fit_pmt->SetParLimits(i, parlo[i], parhi[i]);
			fit_pmt->SetParName(i, par_names[i]);
			if( DEBUG ) 
				cout << "\t" << par_names[i] << ": " << p0save[i] << endl;
		}
		for(int i = 3; i<9; i++) {
			fit_pmt->SetParameter(i, p0save[i]);
			fit_pmt->SetParLimits(i,  parlo[i],  parhi[i]);
			fit_pmt->SetParName(i, par_names[i]);
			if( DEBUG ) 
				cout << "\t" << par_names[i] << ": " << p0save[i] << endl;
		}
		
		fit_pmt->FixParameter(9, ped0);
		
		
		canvas->Clear();
		hist_norm->Draw();
		
		double minx = ped0 + 5.0*sigma0;
		double maxx = (double)hist_norm->GetBinCenter(hist_norm->GetXaxis()->GetNbins()+1);
		
		fit_pmt->SetRange(minx,maxx);
		fit_pmt->FixParameter(1,fit_pmt->GetParameter(1));
		
		if( hv==1000 && anode==16 ) {
			fit_pmt->SetParameter( 0, 65. );
			fit_pmt->SetParameter( 2, 0.0086 );
			fit_pmt->SetParameter( 3, 17.2 );
			fit_pmt->SetParameter( 4, 0.22 );
			fit_pmt->SetParameter( 5, 0.16 );
			fit_pmt->SetParameter( 6, 0.28 );
			fit_pmt->SetParameter( 7, 0.56 );
		}
		if( hv==1000 && anode==36 ) {
			fit_pmt->SetParameter( 0, 65. );
			fit_pmt->SetParameter( 2, 0.0085 );
			fit_pmt->SetParameter( 3, 15.3 );
			fit_pmt->SetParameter( 4, 0.29 );
			fit_pmt->SetParameter( 5, 0.15 );
			fit_pmt->SetParameter( 6, 0.31 );
			fit_pmt->SetParameter( 7, 0.57 );
		}
		if( hv==1000 && anode==39 ) {
			fit_pmt->SetParameter( 0, 65. );
			fit_pmt->SetParameter( 2, 0.0094 );
			fit_pmt->SetParameter( 3, 15.8 );
			fit_pmt->SetParameter( 4, 0.23 );
			fit_pmt->SetParameter( 5, 0.16 );
			fit_pmt->SetParameter( 6, 0.24 );
			fit_pmt->SetParameter( 7, 0.53 );
		}
		if( hv==1000 && anode==47 ) {
			fit_pmt->SetParameter( 0, 65. );
			fit_pmt->SetParameter( 2, 0.0102 );
			fit_pmt->SetParameter( 3, 17.5 );
			fit_pmt->SetParameter( 4, 0.22 );
			fit_pmt->SetParameter( 5, 0.11 );
			fit_pmt->SetParameter( 6, 0.26 );
			fit_pmt->SetParameter( 7, 0.50 );
		}
		if( hv==950 && anode==33 ) {
			fit_pmt->SetParameter( 0, 50. );
			fit_pmt->SetParameter( 2, 0.0083 );
			fit_pmt->SetParameter( 3, 17.5 );
			fit_pmt->SetParameter( 4, 0.22 );
			fit_pmt->SetParameter( 5, 0.11 );
			fit_pmt->SetParameter( 6, 0.26 );
			fit_pmt->SetParameter( 7, 0.50 );
		}
		
		if(hv < 925) {
			fit_pmt->SetParLimits(2,0.9*fit_pmt->GetParameter(2),1.1*fit_pmt->GetParameter(2));
		}
		
		//fit_pmt->FixParameter(8,2.5);
		
		TFitResultPtr result = hist_norm->Fit( fit_pmt, "SR0" );
		double chi2 = result->Chi2() / result->Ndf();
		
		fit_pmt->SetRange(ped0 - 10.*sigma0, maxx);
		
		hist_norm->Draw();
		fit_pmt->SetNpx(1000);
		fit_pmt->Draw("same");
		
		TF1 *fd0 = new TF1( "fd0", fit0, -10., maxx, 10 );
		fd0->SetParameters( fit_pmt->GetParameters() );
		fd0->SetLineColor( kCyan );
		fd0->SetFillColor( kCyan );
		fd0->SetFillStyle( 3004 );
		fd0->Draw( "same" );
		
		TF1 *fd1 = new TF1( "fd1", fit1, -10., maxx, 10 );
		fd1->SetParameters( fit_pmt->GetParameters() );
		fd1->SetFillColor( kBlue );
		fd1->SetLineColor( kBlue );
		fd1->SetFillStyle( 3004 );
		fd1->Draw( "same" );
		
		TF1 *fd2 = new TF1( "fd2", fit2, -10., maxx, 10 );
		fd2->SetParameters( fit_pmt->GetParameters() );
		fd2->SetFillColor( kMagenta );
		fd2->SetLineColor( kMagenta );
		fd2->SetFillStyle( 3005 );
		fd2->Draw( "same" );
		
		TF1 *fd3 = new TF1( "fd3", fit3, -10., maxx, 10 );
		fd3->SetParameters( fit_pmt->GetParameters() );
		fd3->SetFillColor( kGreen );
		fd3->SetLineColor( kGreen );
		fd3->SetFillStyle( 3005 );
		fd3->Draw( "same" );
		
		canvas->SetLogy(true);
    //canvas->Print(Form("~/RICH/pavel/run19_fits/run19_hspe%d_LOGGED.pdf", anode));
		//canvas->SetLogy(false);
		//canvas->Print(pdf_name);
		
		canvas->Print(Form("canvas_hv%d_%d.png",hv,anode));
		/*
		//ofstream outf( Form("HA%04d.%02d.w%d.dat",pmt.c_str(),anode,run_num) );
		ofstream outf(Form("output_hv%d_%d.dat",hv,anode));
		
		//outf << chi2 << endl;
		for( int i=0; i<10; i++ ) {
			double loc_par  = fit_pmt->GetParameter(i);
			double loc_parE = fit_pmt->GetParError(i);
			outf << loc_par << "     " << loc_parE << endl;
		}
		outf.close();
		*/
	}
	
	
	
	return 0;
	
}

double fitf(double *xx, double *par)
{
	//c1=v2/v1
	//c3=v3/v1
	//c2=a4/(1-a2)
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
	
	//Decide upper limit of n
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

	//Decide upper limit of n
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

	//Decide upper limit of n
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

	//Decide upper limit of n
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

	//Decide upper limit of n
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














/*void fitfn(double *xx, double *par)
{
// I remember I thought I could make fitf faster 
// 	but I don't know what I was trying to do

	double scale = par[0], sigma = par[1], mu = par[2], v1 = par[3], a2 = par[4], 
	       c1 = par[5], c2 = par [6], c3 = par[7], xi = par[8], n = par[9];

	double a = xx[0]/scale, v2 = c1*v1, a3 = c2*(1-a2), v3 = c3*v1, sigmaa = sigma/scale;
	double v = v1*(1-a2-a3) + v2*a2 + v3*a3;

	double seff = TMath::Sqrt(sigmaa*sigmaa + n/(v*v * xi));

	//if(TMath::Abs(a-n/v)>8*seff)
	//{
	//	continue;
	//}

	double q0=1;
	
	if(n!=0)
	{
		q0=0;
	}
	double pt=0;
	for(int m = 1; m <= Mmax; m++)
	{	
		if(TMath::Abs(a-n/v)>8*seff)
        	{
                	continue;
        	}
		pt += TMath::Poisson(m, mu) * T(n, m, v1, a2, v2, a3, v3);
	}
	return (G(a, n/v, seff)*(TMath::Exp(-1*mu)*q0+pt)) / scale * Size;

}*/

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

