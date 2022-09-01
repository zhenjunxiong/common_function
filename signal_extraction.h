#include<iostream>
#include<iomanip>
#include<fstream>
#include"TLatex.h"
#include"TFile.h"
#include"TStyle.h"
#include"AliHistogramManager.h"
#include"TMatrixDSym.h"
#include"TFitResultPtr.h"
#include"/home/zjx/include/common_files/Common_Seting.h"
#include "/home/zjx/include/useful.h" 

using namespace std;

double mass_low      = 2.0;
double mass_high     = 4.0;
const double massMax_jpsi  = 3.16;
const double massMin_jpsi  = 2.92;
char buf[1024];
TH1F *significance;
TH1F *hSignal4Fit;
TH1F *hMCJpsi;
TH1F *hFitJpsi ;
TF1  *fBkgFitFunction;
TF1  *fGlobalFitFunction;

void fit_directly(TH1F*,TH1F*,string,string,int,double);
void fit_likesign(TH1F*,TH1F*,TH1F*,string,string,int,double);
Double_t GlobalFitFunction(Double_t *x,Double_t* par);
void SetYRangeHist(TH1 *his);

void fit_directly(TH1F* hMass,TH1F* hunlike,string outputpath,string outputname,int pol_n,double result_imformation[5])
{
    gStyle->SetOptStat(0);
    TH1F* hMC = (TH1F*)hMass->Clone();
    TH1F* hDataUnlike = (TH1F*)hunlike->Clone();
    hDataUnlike->Sumw2();
	// hDataLike->Sumw2();
	hMC->Sumw2();
    if(hMC->Integral() !=0)
	{
		hMC->Scale(1./hMC->Integral());
	}
	else  return -1;
    //set style
	hMC->SetLineColor(kBlue);
	hMC->SetLineWidth(3);
	hMC->SetLineStyle(3);

    hMCJpsi = (TH1F*) hMC->Clone("hMCJpsi");
	TH1F *hMCTemp = (TH1F*) hMC->Clone("hMCTemp");
    //set fit function parameters
    sprintf(buf,"pol%d",pol_n);
	fBkgFitFunction = new TF1("fBkgFitFunction",buf,0,10);
	fGlobalFitFunction = new TF1("GlobalFitFunction",GlobalFitFunction,0.0,10.0,1+fBkgFitFunction->GetNpar());
	fGlobalFitFunction->SetLineColor(2);
	fGlobalFitFunction->SetNpx(10000);
	fGlobalFitFunction->SetLineWidth(3);
	fBkgFitFunction->SetLineColor(4);
	fBkgFitFunction->SetNpx(10000);
	fBkgFitFunction->SetLineWidth(3);
    //initialize fit result
	int nTotalPars = fGlobalFitFunction->GetNpar();
	double nJpsi            = -999.;
	double nJpsiErr         = -999.;
	double nBkgShape        = -999.;
	double nBkgShapeErr     = -999.;
	double nTotalFit        = -999.;
	double nTotalFitErr     = -999.;
    double nTotalBkg        = -999.;
    double nTotalBkgErr     = -999.;
	double nBkgFit          = -999.;
	double nBkgFitErr       = -999.;
	double StoB             = -999.;
	double StoBErr          = -999.;
	double significance     = -999.;
	double significance_err = -999.;
	double scale_factor     = -999.;
	double scale_factor_err = -999.;

	TH1F *hRawJpsi   = new TH1F(*hDataUnlike);
	hRawJpsi->SetName("hRawJpsi");
	TH1F *hRawJpsiTemp = (TH1F*)hRawJpsi->Clone("hRawJpsiTemp");    
    hSignal4Fit = new TH1F(*hDataUnlike);
    hSignal4Fit->Fit(fGlobalFitFunction,"R","",mass_low,mass_high);

    for(int ibin=0;ibin<hMCTemp->GetNbinsX();ibin++)
    {
    hMCTemp->SetBinError(ibin+1,hMCTemp->GetBinContent(ibin+1)*fGlobalFitFunction->GetParError(0));
    hMCTemp->SetBinContent(ibin+1,hMCTemp->GetBinContent(ibin+1)*fGlobalFitFunction->GetParameter(0));
    }

    TFitResultPtr r=hSignal4Fit->Fit(fGlobalFitFunction,"RS","same",mass_low,mass_high);
    TMatrixDSym covTotal = r->GetCovarianceMatrix();
    TMatrixDSym covBG;
    covTotal.GetSub(1,nTotalPars-1,1,nTotalPars-1,covBG);

    double p[nTotalPars]; double perr[nTotalPars];
    for(int ip=0; ip<nTotalPars; ip++)
    {
        p[ip]=fGlobalFitFunction->GetParameter(ip);
        perr[ip]=fGlobalFitFunction->GetParError(ip);
    }

    for(int ip=0; ip<nTotalPars-1; ip++)
    {
        fBkgFitFunction->SetParameter(ip,p[ip+1]);//   p[1], p[2],p[3],p[4]);
        fBkgFitFunction->SetParError(ip,perr[ip+1]);//   p[1], p[2],p[3],p[4]);
    }

    // count the number of jpsi and bkg 
    nTotalFit     = fGlobalFitFunction->Integral(massMin_jpsi,massMax_jpsi);
    nTotalFitErr  = fGlobalFitFunction->IntegralError(massMin_jpsi,massMax_jpsi);
    nJpsi         = hSignal4Fit  ->IntegralAndError(hSignal4Fit  ->FindBin(massMin_jpsi+1.e-9),hSignal4Fit  ->FindBin(massMax_jpsi-1.e-9),nJpsiErr);
    //nBkgShape     = hBkgShape->IntegralAndError(hBkgShape->FindBin(massMin_jpsi+1.e-9),hBkgShape->FindBin(massMax_jpsi-1.e-9),nBkgShapeErr);
    nBkgFit       = fBkgFitFunction->Integral(massMin_jpsi,massMax_jpsi);
    nBkgFitErr    = fBkgFitFunction->IntegralError(massMin_jpsi,massMax_jpsi,fBkgFitFunction->GetParameters(), covBG.GetMatrixArray());

    double binWidth = hDataUnlike->GetBinWidth(1);
    // cout<<"binWidth == "<<binWidth<<endl;
    nTotalFit/=binWidth;
    nTotalFitErr/=binWidth;
    nBkgFit=nBkgFit/binWidth;
    nBkgFitErr/=binWidth;

    nJpsi = nJpsi-nBkgFit; 
    nJpsiErr = sqrt(pow(nJpsiErr,2)+pow(nBkgFitErr,2));
    // get the total background
    nTotalBkg = nBkgFit;
    nTotalBkgErr = sqrt(pow(nBkgFitErr,2));
    // compute S/B
    StoB    = nJpsi/nTotalBkg;
    StoBErr = StoB * sqrt(pow(nJpsiErr/nJpsi,2)+pow(nTotalBkgErr/nTotalBkg,2));
    // compute significance
    significance = nJpsi/sqrt(nJpsi+nTotalBkg);
	significance_err =0.;

    //draw the fit result
    TCanvas *c1 = new TCanvas("c1","",2500,2700);
    c1->Divide(1,2,1e-6,1e-6);
    c1->cd(1);
    gPad->SetBottomMargin(0);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.14);

    hDataUnlike ->GetYaxis()->SetRangeUser(0.2*hDataUnlike ->GetMinimum(),1.8*hDataUnlike ->GetMaximum());

    SetHistogram(hDataUnlike   ,20,2,3,1,2,2,"{m}_{ee} (GeV/c^{2})","Counts per 40 MeV/#font[12]{c}^{2}","");
    SetAxisSize(hDataUnlike    ,0.06,0.06,1.1,1.2,0.05,0.055,0.01,0.01);

    hDataUnlike->Draw("PE");
    fBkgFitFunction->Draw("same");
    fGlobalFitFunction->Draw("same");

    drawLine(massMin_jpsi,massMin_jpsi,0,0.6*hDataUnlike->GetMaximum(),6,7,2);
    drawLine(massMax_jpsi,massMax_jpsi,0,0.6*hDataUnlike->GetMaximum(),6,7,2);

    //<----------------------------------------------------->//
    sprintf(buf,"pp, #sqrt{s}=13 TeV");
    drawLatex(0.17,0.84,buf,22,0.049,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Events:1727.0 M");
    drawLatex(0.17,0.76,buf,22,0.047,1);
    //<----------------------------------------------------->//

    TLegend *legend = new TLegend(0.65,0.72,0.85,0.9);
    legend->AddEntry(hDataUnlike,"Unlike-sign same event","p");
    legend->AddEntry(fGlobalFitFunction,"Fit unlike-sign","L");
    legend->AddEntry(fBkgFitFunction,"Fit Bkg","L");
    drawLegend(legend,22,0.042,0,0);

    c1->cd(2);
    gPad->SetTopMargin(0);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);

    hFitJpsi = new TH1F(*hDataUnlike);
    // sprintf(buf,"hFitJpsi_%d",ipt);
    // hFitJpsi->SetName(buf);
    double temp_max = 0.;
    double temp_min = 0.;
    for(int ibin=1;ibin<=hFitJpsi->GetNbinsX();ibin++)
    {
        double temp = hDataUnlike->GetBinContent(ibin) - fBkgFitFunction->Eval(hFitJpsi->GetBinCenter(ibin));
        hFitJpsi->SetBinContent(ibin,temp);
        if (temp>temp_max)
            temp_max = temp;
        if (temp<temp_min)
            temp_min = temp;
    }
    if(temp_max>0)
        hFitJpsi->SetMaximum(1.6*temp_max);
    else 
        hFitJpsi->SetMinimum(0.6*temp_max);
    if(temp_min>0)
        hFitJpsi->SetMinimum(0.6*temp_min);
    else 
        hFitJpsi->SetMinimum(1.6*temp_min);
    //SetYRangeHist(hFitJpsi);
    SetAxisSize(hFitJpsi,0.06,0.07,1.1,1.1,0.07,0.06,0.01,0.01);
    SetHistogram(hFitJpsi,20,1,3,1,1,2,"#font[12]{m}_{ee} (GeV/#font[12]{c}^{2})","","");
    hFitJpsi->Draw("PE");
    hMCTemp->Draw("histsame");

    //<----------------------------------------------------->//
    sprintf(buf,"N_{J/#psi}: %i #pm %i",int(nJpsi),int(nJpsiErr));
    drawLatex(0.17,0.9,buf,22,0.047,1);
    //<----------------------------------------------------->//
    sprintf(buf,"S/B: %5.3f #pm %5.3f",StoB,StoBErr);
    drawLatex(0.17,0.82,buf,22,0.047,1);
    //<----------------------------------------------------->//
    sprintf(buf,"S/#sqrt{S+B}: %5.1f",significance);
    drawLatex(0.17,0.74,buf,22,0.047,1);
    //<----------------------------------------------------->//
    sprintf(buf,"#chi^{2}/ndf = %2.2f",fGlobalFitFunction->GetChisquare()/fGlobalFitFunction->GetNDF());
    drawLatex(0.72,0.9,buf,22,0.04,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Fit all: %i #pm %i",int(nTotalFit),int(nTotalFitErr));
    drawLatex(0.72,0.82,buf,22,0.04,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Fit bkg: %i #pm %i",int(nBkgFit),int(nBkgFitErr));
    drawLatex(0.72,0.75,buf,22,0.04,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Fit all- Fit bkg: %i #pm %i",int(nTotalFit-nBkgFit),int(sqrt(pow(nTotalFitErr,2)+pow(nBkgFitErr,2))));
    drawLatex(0.72,0.68,buf,22,0.04,1);

    drawLine(massMin_jpsi,massMin_jpsi,hFitJpsi->GetMinimum(),1e5*hFitJpsi->GetMaximum(),6,7,2);
    drawLine(massMax_jpsi,massMax_jpsi,hFitJpsi->GetMinimum(),1e5*hFitJpsi->GetMaximum(),6,7,2);
    drawLine(mass_low,mass_high,0,0,6,7,2);

    sprintf(buf,"%s/fit_result_%s.pdf",outputpath.c_str(),outputname.c_str());
    c1->SaveAs(buf);
    delete c1;
    //get fit result array
    result_imformation[0] = nJpsi;
    result_imformation[1] = nJpsiErr;
    result_imformation[2] = StoB;
    result_imformation[3] = StoBErr;
    result_imformation[4] = significance;
}

void fit_likesign(TH1F* hMass,TH1F* hunlike,TH1F* hlike,string outputpath,string outputname,int pol_n,double result_imformation[5])
{
    gStyle->SetOptStat(0);
    TH1F* hMC = (TH1F*)hMass->Clone();
    TH1F* hDataUnlike = (TH1F*)hunlike->Clone();
    TH1F* hDataLike = (TH1F*)hlike->Clone();
    hDataUnlike->Sumw2();
	hDataLike->Sumw2();
	hMC->Sumw2();
    if(hMC->Integral() !=0)
	{
		hMC->Scale(1./hMC->Integral());
	}
	else  return -1;
    //set style
	hMC->SetLineColor(kBlue);
	hMC->SetLineWidth(3);
	hMC->SetLineStyle(3);

    hMCJpsi = (TH1F*) hMC->Clone("hMCJpsi");
	TH1F *hMCTemp = (TH1F*) hMC->Clone("hMCTemp");
    //set fit function parameters
    sprintf(buf,"pol%d",pol_n);
	fBkgFitFunction = new TF1("fBkgFitFunction",buf,0,10);
	fGlobalFitFunction = new TF1("GlobalFitFunction",GlobalFitFunction,0.0,10.0,1+fBkgFitFunction->GetNpar());
	fGlobalFitFunction->SetLineColor(2);
	fGlobalFitFunction->SetNpx(10000);
	fGlobalFitFunction->SetLineWidth(3);
	fBkgFitFunction->SetLineColor(4);
	fBkgFitFunction->SetNpx(10000);
	fBkgFitFunction->SetLineWidth(3);
    //initialize fit result
	int nTotalPars = fGlobalFitFunction->GetNpar();
	double nJpsi            = -999.;
	double nJpsiErr         = -999.;
	double nBkgShape        = -999.;
	double nBkgShapeErr     = -999.;
	double nTotalFit        = -999.;
	double nTotalFitErr     = -999.;
    double nTotalBkg        = -999.;
    double nTotalBkgErr     = -999.;
	double nBkgFit          = -999.;
	double nBkgFitErr       = -999.;
	double StoB             = -999.;
	double StoBErr          = -999.;
	double significance     = -999.;
	double significance_err = -999.;
	double scale_factor     = -999.;
	double scale_factor_err = -999.;

    scale_factor = 1.;
    TH1F *hBkgShape = new TH1F(*hDataLike);
    hBkgShape->Scale(scale_factor);
    TH1F *hRawJpsi   = new TH1F(*hDataUnlike);
    hRawJpsi->SetName("hRawJpsi");
    hRawJpsi->Add(hBkgShape,-1);
    TH1F *hRawJpsiTemp = (TH1F*)hRawJpsi->Clone("hRawJpsiTemp");

    hSignal4Fit = new TH1F(*hRawJpsi);
    hSignal4Fit->Fit(fGlobalFitFunction,"R","",mass_low,mass_high);

    for(int ibin=0;ibin<hMCTemp->GetNbinsX();ibin++)
    {
    hMCTemp->SetBinError(ibin+1,hMCTemp->GetBinContent(ibin+1)*fGlobalFitFunction->GetParError(0));
    hMCTemp->SetBinContent(ibin+1,hMCTemp->GetBinContent(ibin+1)*fGlobalFitFunction->GetParameter(0));
    }

    TFitResultPtr r=hSignal4Fit->Fit(fGlobalFitFunction,"RS","same",mass_low,mass_high);
    TMatrixDSym covTotal = r->GetCovarianceMatrix();
    TMatrixDSym covBG;
    covTotal.GetSub(1,nTotalPars-1,1,nTotalPars-1,covBG);

    double p[nTotalPars]; double perr[nTotalPars];
    for(int ip=0; ip<nTotalPars; ip++)
    {
        p[ip]=fGlobalFitFunction->GetParameter(ip);
        perr[ip]=fGlobalFitFunction->GetParError(ip);
    }

    for(int ip=0; ip<nTotalPars-1; ip++)
    {
        fBkgFitFunction->SetParameter(ip,p[ip+1]);//   p[1], p[2],p[3],p[4]);
        fBkgFitFunction->SetParError(ip,perr[ip+1]);//   p[1], p[2],p[3],p[4]);
    }

    // count the number of jpsi and bkg 
    nTotalFit     = fGlobalFitFunction->Integral(massMin_jpsi,massMax_jpsi);
    nTotalFitErr  = fGlobalFitFunction->IntegralError(massMin_jpsi,massMax_jpsi);
    nJpsi         = hSignal4Fit  ->IntegralAndError(hSignal4Fit  ->FindBin(massMin_jpsi+1.e-9),hSignal4Fit  ->FindBin(massMax_jpsi-1.e-9),nJpsiErr);
    nBkgShape     = hBkgShape->IntegralAndError(hBkgShape->FindBin(massMin_jpsi+1.e-9),hBkgShape->FindBin(massMax_jpsi-1.e-9),nBkgShapeErr);
    nBkgFit       = fBkgFitFunction->Integral(massMin_jpsi,massMax_jpsi);
    nBkgFitErr    = fBkgFitFunction->IntegralError(massMin_jpsi,massMax_jpsi,fBkgFitFunction->GetParameters(), covBG.GetMatrixArray());

    double binWidth = hDataUnlike->GetBinWidth(1);
    // cout<<"binWidth == "<<binWidth<<endl;
    nTotalFit/=binWidth;
    nTotalFitErr/=binWidth;
    nBkgFit=nBkgFit/binWidth;
    nBkgFitErr/=binWidth;

    nJpsi = nJpsi-nBkgFit; 
    nJpsiErr = sqrt(pow(nJpsiErr,2)+pow(nBkgFitErr,2));
    // get the total background
    nTotalBkg = nBkgShape+nBkgFit;
    nTotalBkgErr = sqrt(pow(nBkgFitErr,2)+pow(nBkgShapeErr,2));
    // compute S/B
    StoB    = nJpsi/nTotalBkg;
    StoBErr = StoB * sqrt(pow(nJpsiErr/nJpsi,2)+pow(nTotalBkgErr/nTotalBkg,2));
    // compute significance
    significance = nJpsi/sqrt(nJpsi+2*nTotalBkg);
	significance_err =0.;

    //draw the fit result
    TCanvas *c1 = new TCanvas("c1","",2500,2700);
    c1->Divide(1,2,1e-6,1e-6);
    c1->cd(1);
    gPad->SetBottomMargin(0);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.14);

    hDataUnlike ->GetYaxis()->SetRangeUser(0.2*hDataUnlike ->GetMinimum(),1.8*hDataUnlike ->GetMaximum());

    SetHistogram(hDataUnlike   ,20,2,3,1,2,2,"{m}_{ee} (GeV/c^{2})","Counts per 40 MeV/#font[12]{c}^{2}","");
    SetHistogram(hBkgShape     ,24,4,3,1,4,2,"{m}_{ee} (GeV/c^{2})","Counts per 40 MeV/#font[12]{c}^{2}","");
    SetAxisSize(hDataUnlike    ,0.06,0.06,1.1,1.2,0.05,0.055,0.01,0.01);

    hDataUnlike->Draw("PE");
    hBkgShape->Draw("samePE");

    drawLine(massMin_jpsi,massMin_jpsi,0,0.6*hDataUnlike->GetMaximum(),6,7,2);
    drawLine(massMax_jpsi,massMax_jpsi,0,0.6*hDataUnlike->GetMaximum(),6,7,2);

    //<----------------------------------------------------->//
    sprintf(buf,"pp, #sqrt{s}=13 TeV");
    drawLatex(0.17,0.84,buf,22,0.049,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Events:1727.0 M");
    drawLatex(0.17,0.76,buf,22,0.047,1);
    //<----------------------------------------------------->//

    TLegend *legend = new TLegend(0.65,0.72,0.85,0.9);
    legend->AddEntry(hDataUnlike,"Unlike-sign same event","p");
    legend->AddEntry(hBkgShape,"Like-sign same event","p");
    legend->AddEntry(hRawJpsi,"Raw J/#psi signal","p");
    drawLegend(legend,22,0.042,0,0);

    TLegend *legend_2 = new TLegend(0.65,0.54,0.85,0.72);
    legend_2->AddEntry(hMCTemp,"MC","L");
    sprintf(buf,"MC+Pol%d Fit",pol_n);
    legend_2->AddEntry(fGlobalFitFunction,buf,"L");
    sprintf(buf,"Pol%d Residual",pol_n);
    legend_2->AddEntry(fBkgFitFunction,buf,"L");
    drawLegend(legend_2,22,0.041,0,0);

    c1->cd(2);
    gPad->SetTopMargin(0);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);

    SetYRangeHist(hRawJpsi);
    SetAxisSize(hRawJpsi,0.06,0.07,1.1,1.1,0.07,0.06,0.01,0.01);
    SetHistogram(hRawJpsi,20,1,3,1,1,2,"#font[12]{m}_{ee} (GeV/#font[12]{c}^{2})","","");

    hRawJpsi->Draw("PE");
    hMCTemp->Draw("histsame");
    fBkgFitFunction->Draw("same");
    fGlobalFitFunction->Draw("same");

    //<----------------------------------------------------->//
    sprintf(buf,"N_{J/#psi}: %i #pm %i",int(nJpsi),int(nJpsiErr));
    drawLatex(0.17,0.9,buf,22,0.047,1);
    //<----------------------------------------------------->//
    sprintf(buf,"S/B: %5.3f #pm %5.3f",StoB,StoBErr);
    drawLatex(0.17,0.82,buf,22,0.047,1);
    //<----------------------------------------------------->//
//				sprintf(buf,"S/#sqrt{S+2B}: %5.1f #pm %5.2f",significance,significance_err);
    sprintf(buf,"S/#sqrt{S+2B}: %5.1f ",significance);
    drawLatex(0.17,0.74,buf,22,0.047,1);
    //<----------------------------------------------------->//
    sprintf(buf,"#chi^{2}/ndf = %2.2f",fGlobalFitFunction->GetChisquare()/fGlobalFitFunction->GetNDF());
    drawLatex(0.72,0.9,buf,22,0.04,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Fit all: %i #pm %i",int(nTotalFit),int(nTotalFitErr));
    drawLatex(0.72,0.82,buf,22,0.04,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Fit bkg: %i #pm %i",int(nBkgFit),int(nBkgFitErr));
    drawLatex(0.72,0.75,buf,22,0.04,1);
    //<----------------------------------------------------->//
    sprintf(buf,"Fit all- Fit bkg: %i #pm %i",int(nTotalFit-nBkgFit),int(sqrt(pow(nTotalFitErr,2)+pow(nBkgFitErr,2))));
    drawLatex(0.72,0.68,buf,22,0.04,1);
    //<----------------------------------------------------->//
    //sprintf(buf,"Fit J/#psi %i #pm %i",int(mc_jpsi),int(mc_jpsi_err));
    //drawLatex(0.72,0.61,buf,22,0.04,1);
    //<----------------------------------------------------->//

    drawLine(massMin_jpsi,massMin_jpsi,hRawJpsi->GetMinimum(),1e5*hRawJpsi->GetMaximum(),6,7,2);
    drawLine(massMax_jpsi,massMax_jpsi,hRawJpsi->GetMinimum(),1e5*hRawJpsi->GetMaximum(),6,7,2);
    drawLine(mass_low,mass_high,0,0,6,7,2);

    sprintf(buf,"%s/fit_result_%s.pdf",outputpath.c_str(),outputname.c_str());
    c1->SaveAs(buf);
    delete c1;
    //get fit result array
    result_imformation[0] = nJpsi;
    result_imformation[1] = nJpsiErr;
    result_imformation[2] = StoB;
    result_imformation[3] = StoBErr;
    result_imformation[4] = significance;
}

Double_t GlobalFitFunction(Double_t *x, Double_t* par) {

	Double_t val = hMCJpsi->GetBinContent(hMCJpsi->FindBin(x[0]));
	val *= par[0];
	for(Int_t i=0;i<fBkgFitFunction->GetNpar();++i) {
		fBkgFitFunction->SetParameter(i, par[i+1]);
	}
	val += fBkgFitFunction->Eval(x[0]);
	return val;
}

void SetYRangeHist(TH1 *his)
{
	if(his->GetMinimum() > 0)
	{
		his->GetYaxis()->SetRangeUser(0,1.6*his->GetMaximum());
	}
	else if(his->GetMinimum() <= 0)
	{
		his->GetYaxis()->SetRangeUser(1.6*his->GetMinimum(),1.6*his->GetMaximum());
	}

}