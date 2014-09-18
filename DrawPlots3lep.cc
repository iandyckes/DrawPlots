//  object#:   no number -> ee, 2 -> MuMu
//
//  e.g. samplehisto -> histogram of ee dilepton mass
//       samplehisto2 -> histogram of MuMu dilepton mass

#include <vector>
#include <iostream>
#include <stdexcept>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "THStack.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

//----------Global Variables---------------//
int mll_minbin=0;  //  Why not 50?
int mll_maxbin=200;
int mll_binsize=5;
int met_minbin=0;  //
int met_maxbin=200;
int met_binsize=2;
//-----------------------------------------//


void h_format(TH1F *histo)
{

  TString name=histo->GetName();
  Color_t color=kMagenta+4;

  if(name.Contains("tar") && name.Contains("mll"))
	{histo->Rebin(5);}
  else if (name.Contains("inc") && name.Contains("mll"))
	{histo->Rebin(5);}
  else if (name.Contains("met"))
	{histo->Rebin(2);
	  histo->SetMarkerSize(0.5);}

  if(name.Contains("zjets"))     //long because first, needs more formatting
	{
	  color=kCyan;
	  histo->SetFillColor(color);

	  if(name.Contains("mll"))
		{ 
		  histo->GetXaxis()->SetRangeUser(50,200); 
		  histo->GetYaxis()->SetRangeUser(1e0,5e6);//1e1, 5e8

		  if(name.Contains("ee"))
			{
			  histo->GetXaxis()->SetTitle("m_{ee} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 5 GeV");
			}
		  else if(name.Contains("mumu"))
			{
			  histo->GetXaxis()->SetTitle("m_{mumu} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 5 GeV");
			}
		}
	  else if(name.Contains("met"))
		{  
		  histo->GetXaxis()->SetRangeUser(0,200);  //350
		  histo->GetYaxis()->SetRangeUser(1e-1,5e6);//1e0,5e7
		 
		  if(name.Contains("ee"))
			{
			  histo->GetXaxis()->SetTitle("met_{ee} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 2 GeV");
			}
		  else if(name.Contains("mumu"))
			{
			  histo->GetXaxis()->SetTitle("met_{mumu} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 2 GeV");
			}
		}

	}
  else if(name.Contains("singlet"))
	{color=kGreen-3;
	  histo->SetFillColor(color);}
  else if(name.Contains("wz"))
	{color=kBlue;
	  histo->SetFillColor(color);}
  else if(name.Contains("zz"))
	{color=kRed;
	  histo->SetFillColor(color);}
  else if(name.Contains("ttbar"))
	{color=kYellow;
	  histo->SetFillColor(color);}
  else if(name.Contains("ww"))
	{color=kBlue-3;
	  histo->SetFillColor(color);}
  else if(name.Contains("ttv"))
	{color=kGray+1;
	  histo->SetFillColor(color);}
  else if(name.Contains("vvv"))
	{color=kGray+3;
	  histo->SetFillColor(color);}
  else if(name.Contains("data"))
	{
	  color=kBlack;
	  histo->SetMarkerStyle(8);
	}

  histo->SetMarkerColor(color);
  histo->SetLineColor(kBlack);
 
}

/*
void overflow(TH1F *histo)
{
  double uncertainty=0.;

  histo->SetBinContent(histo->GetNbinsX(),histo->GetBinContent(histo->GetNbinsX())+histo->IntegralAndError(histo->GetNbinsX(),-1,uncertainty));

  histo->SetBinError((histo->GetNbinsX()),sqrt(pow(histo->GetBinError((histo->GetNbinsX())),2) + pow(uncertainty,2)));
}
*/
void overflow(TH1F *histo)
{
  double uncertainty=0.;
  int overflowbin=-1;
  TString name=histo->GetName();

  if(name.Contains("mll"))
	{overflowbin=((mll_maxbin - mll_minbin)/mll_binsize);}
  else if(name.Contains("met"))
	{overflowbin=((met_maxbin - met_minbin)/met_binsize)+1;}

  histo->SetBinContent(overflowbin,histo->GetBinContent(overflowbin)+histo->IntegralAndError(overflowbin,-1,uncertainty));

  histo->SetBinError(overflowbin,sqrt(pow(histo->GetBinError((overflowbin)),2) + pow(uncertainty,2)));
}

int DrawPlots3lep()
{

  //load histos
  TFile *InputFile = new TFile("zjets.root","read");
  TH1F *h_mll_ee_inc_zjets = (TH1F*) InputFile->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_zjets");
  TH1F *h_mll_mumu_inc_zjets = (TH1F*) InputFile->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_zjets");
  TH1F *h_mll_ee_tar_zjets = (TH1F*) InputFile->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_zjets");
  TH1F *h_mll_mumu_tar_zjets = (TH1F*) InputFile->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_zjets");
  TH1F *h_met_ee_inc_zjets = (TH1F*) InputFile->Get("h_met_ee_inc")->Clone("h_met_ee_inc_zjets");
  TH1F *h_met_mumu_inc_zjets = (TH1F*) InputFile->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_zjets");
  TH1F *h_met_ee_tar0_zjets = (TH1F*) InputFile->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_zjets");
  TH1F *h_met_mumu_tar0_zjets = (TH1F*) InputFile->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_zjets");
  TH1F *h_met_ee_tar2_zjets = (TH1F*) InputFile->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_zjets");
  TH1F *h_met_mumu_tar2_zjets = (TH1F*) InputFile->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_zjets");
  TH1F *h_met_emu_inc_zjets = (TH1F*) InputFile->Get("h_met_emu_inc")->Clone("h_met_emu_inc_zjets");
  TH1F *h_met_emu_tar0_zjets = (TH1F*) InputFile->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_zjets");
  TH1F *h_met_emu_tar2_zjets = (TH1F*) InputFile->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_zjets");

  TFile *InputFile_t = new TFile("singlet.root","read");
  TH1F *h_mll_ee_inc_singlet = (TH1F*) InputFile_t->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_singlet"); 
  TH1F *h_mll_mumu_inc_singlet = (TH1F*) InputFile_t->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_singlet"); 
  TH1F *h_mll_ee_tar_singlet = (TH1F*) InputFile_t->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_singlet");
  TH1F *h_mll_mumu_tar_singlet = (TH1F*) InputFile_t->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_singlet");
  TH1F *h_met_ee_inc_singlet = (TH1F*) InputFile_t->Get("h_met_ee_inc")->Clone("h_met_ee_inc_singlet");
  TH1F *h_met_mumu_inc_singlet = (TH1F*) InputFile_t->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_singlet");
  TH1F *h_met_ee_tar0_singlet = (TH1F*) InputFile_t->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_singlet");
  TH1F *h_met_mumu_tar0_singlet = (TH1F*) InputFile_t->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_singlet");
  TH1F *h_met_ee_tar2_singlet = (TH1F*) InputFile_t->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_singlet");
  TH1F *h_met_mumu_tar2_singlet = (TH1F*) InputFile_t->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_singlet");
  TH1F *h_met_emu_inc_singlet = (TH1F*) InputFile_t->Get("h_met_emu_inc")->Clone("h_met_emu_inc_singlet");
  TH1F *h_met_emu_tar0_singlet = (TH1F*) InputFile_t->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_singlet");
  TH1F *h_met_emu_tar2_singlet = (TH1F*) InputFile_t->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_singlet");

  TFile *InputFile_wz = new TFile("wz.root","read");
  TH1F *h_mll_ee_inc_wz = (TH1F*) InputFile_wz->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_wz"); 
  TH1F *h_mll_mumu_inc_wz = (TH1F*) InputFile_wz->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_wz");
  TH1F *h_mll_ee_tar_wz = (TH1F*) InputFile_wz->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_wz");
  TH1F *h_mll_mumu_tar_wz = (TH1F*) InputFile_wz->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_wz"); 
  TH1F *h_met_ee_inc_wz = (TH1F*) InputFile_wz->Get("h_met_ee_inc")->Clone("h_met_ee_inc_wz");
  TH1F *h_met_mumu_inc_wz = (TH1F*) InputFile_wz->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_wz");
  TH1F *h_met_ee_tar0_wz = (TH1F*) InputFile_wz->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_wz");
  TH1F *h_met_mumu_tar0_wz = (TH1F*) InputFile_wz->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_wz");
  TH1F *h_met_ee_tar2_wz = (TH1F*) InputFile_wz->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_wz");
  TH1F *h_met_mumu_tar2_wz = (TH1F*) InputFile_wz->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_wz");
  TH1F *h_met_emu_inc_wz = (TH1F*) InputFile_wz->Get("h_met_emu_inc")->Clone("h_met_emu_inc_wz");
  TH1F *h_met_emu_tar0_wz = (TH1F*) InputFile_wz->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_wz");
  TH1F *h_met_emu_tar2_wz = (TH1F*) InputFile_wz->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_wz");

  TFile *InputFile_zz = new TFile("zz.root","read");
  TH1F *h_mll_ee_inc_zz = (TH1F*) InputFile_zz->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_zz"); 
  TH1F *h_mll_mumu_inc_zz = (TH1F*) InputFile_zz->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_zz");
  TH1F *h_mll_ee_tar_zz = (TH1F*) InputFile_zz->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_zz");
  TH1F *h_mll_mumu_tar_zz = (TH1F*) InputFile_zz->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_zz");
  TH1F *h_met_ee_inc_zz = (TH1F*) InputFile_zz->Get("h_met_ee_inc")->Clone("h_met_ee_inc_zz");
  TH1F *h_met_mumu_inc_zz = (TH1F*) InputFile_zz->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_zz");
  TH1F *h_met_ee_tar0_zz = (TH1F*) InputFile_zz->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_zz");
  TH1F *h_met_mumu_tar0_zz = (TH1F*) InputFile_zz->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_zz");
  TH1F *h_met_ee_tar2_zz = (TH1F*) InputFile_zz->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_zz");
  TH1F *h_met_mumu_tar2_zz = (TH1F*) InputFile_zz->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_zz");
  TH1F *h_met_emu_inc_zz = (TH1F*) InputFile_zz->Get("h_met_emu_inc")->Clone("h_met_emu_inc_zz");
  TH1F *h_met_emu_tar0_zz = (TH1F*) InputFile_zz->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_zz");
  TH1F *h_met_emu_tar2_zz = (TH1F*) InputFile_zz->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_zz");

  TFile *InputFile_ttbar = new TFile("ttbar.root","read");
  TH1F *h_mll_ee_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_ttbar");
  TH1F *h_mll_mumu_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_ttbar");
  TH1F *h_mll_ee_tar_ttbar = (TH1F*) InputFile_ttbar->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_ttbar");
  TH1F *h_mll_mumu_tar_ttbar = (TH1F*) InputFile_ttbar->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_ttbar");
  TH1F *h_met_ee_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_ee_inc")->Clone("h_met_ee_inc_ttbar");
  TH1F *h_met_mumu_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_ttbar");
  TH1F *h_met_ee_tar0_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_ttbar");
  TH1F *h_met_mumu_tar0_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_ttbar");
  TH1F *h_met_ee_tar2_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_ttbar");
  TH1F *h_met_mumu_tar2_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_ttbar");
  TH1F *h_met_emu_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_emu_inc")->Clone("h_met_emu_inc_ttbar");
  TH1F *h_met_emu_tar0_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_ttbar");
  TH1F *h_met_emu_tar2_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_ttbar");

  TFile *InputFile_ww = new TFile("ww.root","read");
  TH1F *h_mll_ee_inc_ww = (TH1F*) InputFile_ww->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_ww");
  TH1F *h_mll_mumu_inc_ww = (TH1F*) InputFile_ww->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_ww");
  TH1F *h_mll_ee_tar_ww = (TH1F*) InputFile_ww->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_ww");
  TH1F *h_mll_mumu_tar_ww = (TH1F*) InputFile_ww->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_ww");
  TH1F *h_met_ee_inc_ww = (TH1F*) InputFile_ww->Get("h_met_ee_inc")->Clone("h_met_ee_inc_ww");
  TH1F *h_met_mumu_inc_ww = (TH1F*) InputFile_ww->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_ww");
  TH1F *h_met_ee_tar0_ww = (TH1F*) InputFile_ww->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_ww");
  TH1F *h_met_mumu_tar0_ww = (TH1F*) InputFile_ww->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_ww");
  TH1F *h_met_ee_tar2_ww = (TH1F*) InputFile_ww->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_ww");
  TH1F *h_met_mumu_tar2_ww = (TH1F*) InputFile_ww->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_ww");
  TH1F *h_met_emu_inc_ww = (TH1F*) InputFile_ww->Get("h_met_emu_inc")->Clone("h_met_emu_inc_ww");
  TH1F *h_met_emu_tar0_ww = (TH1F*) InputFile_ww->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_ww");
  TH1F *h_met_emu_tar2_ww = (TH1F*) InputFile_ww->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_ww");

  TFile *InputFile_ttv = new TFile("ttv.root","read");
  TH1F *h_mll_ee_inc_ttv = (TH1F*) InputFile_ttv->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_ttv");
  TH1F *h_mll_mumu_inc_ttv = (TH1F*) InputFile_ttv->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_ttv");
  TH1F *h_mll_ee_tar_ttv = (TH1F*) InputFile_ttv->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_ttv");
  TH1F *h_mll_mumu_tar_ttv = (TH1F*) InputFile_ttv->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_ttv");
  TH1F *h_met_ee_inc_ttv = (TH1F*) InputFile_ttv->Get("h_met_ee_inc")->Clone("h_met_ee_inc_ttv");
  TH1F *h_met_mumu_inc_ttv = (TH1F*) InputFile_ttv->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_ttv");
  TH1F *h_met_ee_tar0_ttv = (TH1F*) InputFile_ttv->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_ttv");
  TH1F *h_met_mumu_tar0_ttv = (TH1F*) InputFile_ttv->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_ttv");
  TH1F *h_met_ee_tar2_ttv = (TH1F*) InputFile_ttv->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_ttv");
  TH1F *h_met_mumu_tar2_ttv = (TH1F*) InputFile_ttv->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_ttv");
  TH1F *h_met_emu_inc_ttv = (TH1F*) InputFile_ttv->Get("h_met_emu_inc")->Clone("h_met_emu_inc_ttv");
  TH1F *h_met_emu_tar0_ttv = (TH1F*) InputFile_ttv->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_ttv");
  TH1F *h_met_emu_tar2_ttv = (TH1F*) InputFile_ttv->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_ttv");

  TFile *InputFile_vvv = new TFile("vvv.root","read");
  TH1F *h_mll_ee_inc_vvv = (TH1F*) InputFile_vvv->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_vvv");
  TH1F *h_mll_mumu_inc_vvv = (TH1F*) InputFile_vvv->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_vvv");
  TH1F *h_mll_ee_tar_vvv = (TH1F*) InputFile_vvv->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_vvv");
  TH1F *h_mll_mumu_tar_vvv = (TH1F*) InputFile_vvv->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_vvv");
  TH1F *h_met_ee_inc_vvv = (TH1F*) InputFile_vvv->Get("h_met_ee_inc")->Clone("h_met_ee_inc_vvv");
  TH1F *h_met_mumu_inc_vvv = (TH1F*) InputFile_vvv->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_vvv");
  TH1F *h_met_ee_tar0_vvv = (TH1F*) InputFile_vvv->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_vvv");
  TH1F *h_met_mumu_tar0_vvv = (TH1F*) InputFile_vvv->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_vvv");
  TH1F *h_met_ee_tar2_vvv = (TH1F*) InputFile_vvv->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_vvv");
  TH1F *h_met_mumu_tar2_vvv = (TH1F*) InputFile_vvv->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_vvv");
  TH1F *h_met_emu_inc_vvv = (TH1F*) InputFile_vvv->Get("h_met_emu_inc")->Clone("h_met_emu_inc_vvv");
  TH1F *h_met_emu_tar0_vvv = (TH1F*) InputFile_vvv->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_vvv");
  TH1F *h_met_emu_tar2_vvv = (TH1F*) InputFile_vvv->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_vvv");

  TFile *InputFile_data = new TFile("data.root","read");
  TH1F *h_mll_ee_inc_data = (TH1F*) InputFile_data->Get("h_mll_ee_inc")->Clone("h_mll_ee_inc_data");
  TH1F *h_mll_mumu_inc_data = (TH1F*) InputFile_data->Get("h_mll_mumu_inc")->Clone("h_mll_mumu_inc_data");
  TH1F *h_mll_ee_tar_data = (TH1F*) InputFile_data->Get("h_mll_ee_tar")->Clone("h_mll_ee_tar_data");
  TH1F *h_mll_mumu_tar_data = (TH1F*) InputFile_data->Get("h_mll_mumu_tar")->Clone("h_mll_mumu_tar_data");
  TH1F *h_met_ee_inc_data = (TH1F*) InputFile_data->Get("h_met_ee_inc")->Clone("h_met_ee_inc_data");
  TH1F *h_met_mumu_inc_data = (TH1F*) InputFile_data->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_data");
  TH1F *h_met_ee_tar0_data = (TH1F*) InputFile_data->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_data");
  TH1F *h_met_mumu_tar0_data = (TH1F*) InputFile_data->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_data");
  TH1F *h_met_ee_tar2_data = (TH1F*) InputFile_data->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_data");
  TH1F *h_met_mumu_tar2_data = (TH1F*) InputFile_data->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_data");
  TH1F *h_met_emu_inc_data = (TH1F*) InputFile_data->Get("h_met_emu_inc")->Clone("h_met_emu_inc_data");
  TH1F *h_met_emu_tar0_data = (TH1F*) InputFile_data->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_data");
  TH1F *h_met_emu_tar2_data = (TH1F*) InputFile_data->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_data");

  //-------------------------------------------------------------------------------------------------//

  //--------------------------------------Put in Vectors---------------------------------------------//


  vector<TH1F*> v_mll_ee_inc;
  v_mll_ee_inc.push_back(h_mll_ee_inc_zjets);
  v_mll_ee_inc.push_back(h_mll_ee_inc_wz);
  v_mll_ee_inc.push_back(h_mll_ee_inc_zz);
  v_mll_ee_inc.push_back(h_mll_ee_inc_ttbar);
  v_mll_ee_inc.push_back(h_mll_ee_inc_singlet);
  v_mll_ee_inc.push_back(h_mll_ee_inc_ww);
  v_mll_ee_inc.push_back(h_mll_ee_inc_ttv);
  v_mll_ee_inc.push_back(h_mll_ee_inc_vvv);
  v_mll_ee_inc.push_back(h_mll_ee_inc_data);
  
  vector<TH1F*> v_mll_mumu_inc;
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_zjets);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_wz);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_zz);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_ttbar);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_singlet);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_ww);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_ttv);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_vvv);
  v_mll_mumu_inc.push_back(h_mll_mumu_inc_data);

  vector<TH1F*> v_mll_ee_tar;
  v_mll_ee_tar.push_back(h_mll_ee_tar_zjets);
  v_mll_ee_tar.push_back(h_mll_ee_tar_wz);
  v_mll_ee_tar.push_back(h_mll_ee_tar_zz);
  v_mll_ee_tar.push_back(h_mll_ee_tar_ttbar);
  v_mll_ee_tar.push_back(h_mll_ee_tar_singlet);
  v_mll_ee_tar.push_back(h_mll_ee_tar_ww);
  v_mll_ee_tar.push_back(h_mll_ee_tar_ttv);
  v_mll_ee_tar.push_back(h_mll_ee_tar_vvv);
  v_mll_ee_tar.push_back(h_mll_ee_tar_data);

  vector<TH1F*> v_mll_mumu_tar;
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_zjets);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_wz);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_zz);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_ttbar);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_singlet);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_ww);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_ttv);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_vvv);
  v_mll_mumu_tar.push_back(h_mll_mumu_tar_data);

  vector<TH1F*> v_met_ee_inc;
  v_met_ee_inc.push_back(h_met_ee_inc_zjets);
  v_met_ee_inc.push_back(h_met_ee_inc_wz);
  v_met_ee_inc.push_back(h_met_ee_inc_zz);
  v_met_ee_inc.push_back(h_met_ee_inc_ttbar);
  v_met_ee_inc.push_back(h_met_ee_inc_singlet);
  v_met_ee_inc.push_back(h_met_ee_inc_ww);
  v_met_ee_inc.push_back(h_met_ee_inc_ttv);
  v_met_ee_inc.push_back(h_met_ee_inc_vvv);
  v_met_ee_inc.push_back(h_met_ee_inc_data);

  vector<TH1F*> v_met_mumu_inc;
  v_met_mumu_inc.push_back(h_met_mumu_inc_zjets);
  v_met_mumu_inc.push_back(h_met_mumu_inc_wz);
  v_met_mumu_inc.push_back(h_met_mumu_inc_zz);
  v_met_mumu_inc.push_back(h_met_mumu_inc_ttbar);
  v_met_mumu_inc.push_back(h_met_mumu_inc_singlet);
  v_met_mumu_inc.push_back(h_met_mumu_inc_ww);
  v_met_mumu_inc.push_back(h_met_mumu_inc_ttv);
  v_met_mumu_inc.push_back(h_met_mumu_inc_vvv);
  v_met_mumu_inc.push_back(h_met_mumu_inc_data);

  vector<TH1F*> v_met_emu_inc;
  v_met_emu_inc.push_back(h_met_emu_inc_zjets);
  v_met_emu_inc.push_back(h_met_emu_inc_wz);
  v_met_emu_inc.push_back(h_met_emu_inc_zz);
  v_met_emu_inc.push_back(h_met_emu_inc_ttbar);
  v_met_emu_inc.push_back(h_met_emu_inc_singlet);
  v_met_emu_inc.push_back(h_met_emu_inc_ww);
  v_met_emu_inc.push_back(h_met_emu_inc_ttv);
  v_met_emu_inc.push_back(h_met_emu_inc_vvv);
  v_met_emu_inc.push_back(h_met_emu_inc_data);

  vector<TH1F*> v_met_ee_tar0;
  v_met_ee_tar0.push_back(h_met_ee_tar0_zjets);
  v_met_ee_tar0.push_back(h_met_ee_tar0_wz);
  v_met_ee_tar0.push_back(h_met_ee_tar0_zz);
  v_met_ee_tar0.push_back(h_met_ee_tar0_ttbar);
  v_met_ee_tar0.push_back(h_met_ee_tar0_singlet);
  v_met_ee_tar0.push_back(h_met_ee_tar0_ww);
  v_met_ee_tar0.push_back(h_met_ee_tar0_ttv);
  v_met_ee_tar0.push_back(h_met_ee_tar0_vvv);
  v_met_ee_tar0.push_back(h_met_ee_tar0_data);

  vector<TH1F*> v_met_mumu_tar0;
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_zjets);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_wz);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_zz);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_ttbar);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_singlet);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_ww);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_ttv);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_vvv);
  v_met_mumu_tar0.push_back(h_met_mumu_tar0_data);

  vector<TH1F*> v_met_emu_tar0;
  v_met_emu_tar0.push_back(h_met_emu_tar0_zjets);
  v_met_emu_tar0.push_back(h_met_emu_tar0_wz);
  v_met_emu_tar0.push_back(h_met_emu_tar0_zz);
  v_met_emu_tar0.push_back(h_met_emu_tar0_ttbar);
  v_met_emu_tar0.push_back(h_met_emu_tar0_singlet);
  v_met_emu_tar0.push_back(h_met_emu_tar0_ww);
  v_met_emu_tar0.push_back(h_met_emu_tar0_ttv);
  v_met_emu_tar0.push_back(h_met_emu_tar0_vvv);
  v_met_emu_tar0.push_back(h_met_emu_tar0_data);

  vector<TH1F*> v_met_ee_tar2;
  v_met_ee_tar2.push_back(h_met_ee_tar2_zjets);
  v_met_ee_tar2.push_back(h_met_ee_tar2_wz);
  v_met_ee_tar2.push_back(h_met_ee_tar2_zz);
  v_met_ee_tar2.push_back(h_met_ee_tar2_ttbar);
  v_met_ee_tar2.push_back(h_met_ee_tar2_singlet);
  v_met_ee_tar2.push_back(h_met_ee_tar2_ww);
  v_met_ee_tar2.push_back(h_met_ee_tar2_ttv);
  v_met_ee_tar2.push_back(h_met_ee_tar2_vvv);
  v_met_ee_tar2.push_back(h_met_ee_tar2_data);

  vector<TH1F*> v_met_mumu_tar2;
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_zjets);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_wz);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_zz);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_ttbar);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_singlet);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_ww);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_ttv);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_vvv);
  v_met_mumu_tar2.push_back(h_met_mumu_tar2_data);

  vector<TH1F*> v_met_emu_tar2;
  v_met_emu_tar2.push_back(h_met_emu_tar2_zjets);
  v_met_emu_tar2.push_back(h_met_emu_tar2_wz);
  v_met_emu_tar2.push_back(h_met_emu_tar2_zz);
  v_met_emu_tar2.push_back(h_met_emu_tar2_ttbar);
  v_met_emu_tar2.push_back(h_met_emu_tar2_singlet);
  v_met_emu_tar2.push_back(h_met_emu_tar2_ww);
  v_met_emu_tar2.push_back(h_met_emu_tar2_ttv);
  v_met_emu_tar2.push_back(h_met_emu_tar2_vvv);
  v_met_emu_tar2.push_back(h_met_emu_tar2_data);
 

  //-------------------------------------------------------------------------------------------------//

  //loop through ee and mumu and subract eu*R.  Before formatting, overflow, summing, and stacking
  int size=v_mll_ee_inc.size();
  
  double Ree=0.440978;
  double Rmumu=0.524251;

  for(int i=0; i<size; i++)
	{
	  v_met_ee_inc[i]->Add(v_met_emu_inc[i],-Ree);      //don't forget minus
	  v_met_mumu_inc[i]->Add(v_met_emu_inc[i],-Rmumu);  //can use clones to not permanently change
	  v_met_ee_tar0[i]->Add(v_met_emu_tar0[i],-Ree);      
	  v_met_mumu_tar0[i]->Add(v_met_emu_tar0[i],-Rmumu);  
	  v_met_ee_tar2[i]->Add(v_met_emu_tar2[i],-Ree);      
	  v_met_mumu_tar2[i]->Add(v_met_emu_tar2[i],-Rmumu); 
	}
  

  //---------------------------------------Format and Overflow---------------------------------------//

	
  for(int i=0; i < size; i++)
	{
	  // cout<<"i = "<<i<<endl;
	  // cout<<"length = "<<v_mll_mumu_inc.size()<<endl;
	  
	  h_format(v_mll_ee_inc[i]);
	  h_format(v_mll_mumu_inc[i]);
	  h_format(v_mll_ee_tar[i]);
	  h_format(v_mll_mumu_tar[i]);
	  h_format(v_met_ee_inc[i]);
	  h_format(v_met_mumu_inc[i]);
	  h_format(v_met_emu_inc[i]);
	  h_format(v_met_ee_tar0[i]);
	  h_format(v_met_mumu_tar0[i]);
	  h_format(v_met_emu_tar0[i]);
	  h_format(v_met_ee_tar2[i]);
	  h_format(v_met_mumu_tar2[i]);
	  h_format(v_met_emu_tar2[i]);

	  overflow(v_mll_ee_inc[i]);
	  overflow(v_mll_mumu_inc[i]);
	  overflow(v_mll_ee_tar[i]);
	  overflow(v_mll_mumu_tar[i]);
	  overflow(v_met_ee_inc[i]);
	  overflow(v_met_mumu_inc[i]);
	  overflow(v_met_emu_inc[i]);
	  overflow(v_met_ee_tar0[i]);
	  overflow(v_met_mumu_tar0[i]);
	  overflow(v_met_emu_tar0[i]);
	  overflow(v_met_ee_tar2[i]);
	  overflow(v_met_mumu_tar2[i]);
	  overflow(v_met_emu_tar2[i]);
	}

  cout<<"\n End Formatting \n"<<endl;
  //-------------------------------------------------------------------------------------------------//

  //-------------------------------------  Sum  -----------------------------------------------------//
  
  TH1F *sum = (TH1F*) v_mll_ee_inc[0]->Clone("sum");
  TH1F *sum2 = (TH1F*) v_mll_mumu_inc[0]->Clone("sum2");
  TH1F *tarsum = (TH1F*) v_mll_ee_tar[0]->Clone("tarsum");
  TH1F *tarsum2 = (TH1F*) v_mll_mumu_tar[0]->Clone("tarsum2");
  TH1F *h_metsum_ee_inc = (TH1F*) v_met_ee_inc[0]->Clone("h_metsum_ee_inc");
  TH1F *h_metsum_mumu_inc = (TH1F*) v_met_mumu_inc[0]->Clone("h_metsum_mumu_inc");
  TH1F *h_metsum_emu_inc = (TH1F*) v_met_emu_inc[0]->Clone("h_metsum_emu_inc");
  TH1F *h_metsum_ee_tar0 = (TH1F*) v_met_ee_tar0[0]->Clone("h_metsum_ee_tar0");
  TH1F *h_metsum_mumu_tar0 = (TH1F*) v_met_mumu_tar0[0]->Clone("h_metsum_mumu_tar0");
  TH1F *h_metsum_emu_tar0 = (TH1F*) v_met_emu_tar0[0]->Clone("h_metsum_emu_tar0");
  TH1F *h_metsum_ee_tar2 = (TH1F*) v_met_ee_tar2[0]->Clone("h_metsum_ee_tar2");
  TH1F *h_metsum_mumu_tar2 = (TH1F*) v_met_mumu_tar2[0]->Clone("h_metsum_mumu_tar2");
  TH1F *h_metsum_emu_tar2 = (TH1F*) v_met_emu_tar2[0]->Clone("h_metsum_emu_tar2");

  
  for(int i=1; i < size-1; i++)    //start at 1 so skip cloned, -1 so you don't add the data.
	{
	  // cout<<"i = "<<i<<endl;
	  // cout<<"length = "<<v_mll_mumu_inc.size()<<endl;
	 
	  sum->Add(v_mll_ee_inc[i]);
	  sum2->Add(v_mll_mumu_inc[i]);
	  tarsum->Add(v_mll_ee_tar[i]);
	  tarsum2->Add(v_mll_mumu_tar[i]);
	  h_metsum_ee_inc->Add(v_met_ee_inc[i]);
	  h_metsum_mumu_inc->Add(v_met_mumu_inc[i]);
	  h_metsum_emu_inc->Add(v_met_emu_inc[i]);
	  h_metsum_ee_tar0->Add(v_met_ee_tar0[i]);
	  h_metsum_mumu_tar0->Add(v_met_mumu_tar0[i]);
	  h_metsum_emu_tar0->Add(v_met_emu_tar0[i]);
	  h_metsum_ee_tar2->Add(v_met_ee_tar2[i]);
	  h_metsum_mumu_tar2->Add(v_met_mumu_tar2[i]);
	  h_metsum_emu_tar2->Add(v_met_emu_tar2[i]);
	}

  overflow(sum);
  overflow(sum2);
  overflow(tarsum);
  overflow(tarsum2);
  overflow(h_metsum_ee_inc);
  overflow(h_metsum_mumu_inc);
  overflow(h_metsum_emu_inc);
  overflow(h_metsum_ee_tar0);
  overflow(h_metsum_mumu_tar0);
  overflow(h_metsum_emu_tar0);
  overflow(h_metsum_ee_tar2);
  overflow(h_metsum_mumu_tar2);
  overflow(h_metsum_emu_tar2);

  cout<<"\n End Summing \n"<<endl;

  //--------------------------------------------------------------------------------//

  //Get rid of stat box
  gStyle->SetOptStat(0);

  //-------------------------Stacked Histograms--------------------------------------//

  //Declare Stacked Histogram
  THStack *hs = new THStack("hs","ee");
  THStack *hs2= new THStack("hs","mumu");
  THStack *hstar = new THStack("hstar","tar ee");
  THStack *hstar2 = new THStack("hstar2","tar mumu");
  THStack *hs_met_ee_inc = new THStack("hs_met_ee_inc","met ee inc");
  THStack *hs_met_mumu_inc = new THStack("hs_met_mumu_inc","met mumu inc");
  THStack *hs_met_emu_inc = new THStack("hs_met_emu_inc","met emu inc");
  THStack *hs_met_ee_tar0 = new THStack("hs_met_ee_tar0","met ee tar0");
  THStack *hs_met_mumu_tar0 = new THStack("hs_met_mumu_tar0","met mumu tar0");
  THStack *hs_met_emu_tar0 = new THStack("hs_met_emu_tar0","met emu tar0");
  THStack *hs_met_ee_tar2 = new THStack("hs_met_ee_tar2","met ee tar2");
  THStack *hs_met_mumu_tar2 = new THStack("hs_met_mumu_tar2","met mumu tar2");
  THStack *hs_met_emu_tar2 = new THStack("hs_met_emu_tar2","met emu tar2");

  //Add to Stacked Histograms
  for(int i=size-2; i >= 0; i--)  //Reverse,don't add data  >=0????????????????????????????????????????????????????
	{
	  hs->Add(v_mll_ee_inc[i]);
	  hs2->Add(v_mll_mumu_inc[i]);
	  hstar->Add(v_mll_ee_tar[i]);
	  hstar2->Add(v_mll_mumu_tar[i]);
	  hs_met_ee_inc->Add(v_met_ee_inc[i]);
	  hs_met_mumu_inc->Add(v_met_mumu_inc[i]);
	  hs_met_emu_inc->Add(v_met_emu_inc[i]);
	  hs_met_ee_tar0->Add(v_met_ee_tar0[i]);
	  hs_met_mumu_tar0->Add(v_met_mumu_tar0[i]);
	  hs_met_emu_tar0->Add(v_met_emu_tar0[i]);
	  hs_met_ee_tar2->Add(v_met_ee_tar2[i]);
	  hs_met_mumu_tar2->Add(v_met_mumu_tar2[i]);
	  hs_met_emu_tar2->Add(v_met_emu_tar2[i]);
	}
  
  //-----------------------------------------------------------------------------------//
  //----------------------------------Drawing------------------------------------------//
  //-----------------------------------------------------------------------------------//  
  
  //-----------------------mll ee inc------------------------- -//
  //Canvas
  TCanvas *c1=new TCanvas("c1","mll ee inc ",800,800);
  TPad *pad_h = new TPad("pad_h","Histo Pad",0., 0, 1., 0.8);
  TPad *pad_r = new TPad("pad_r","Ratio Pad",0., 0.8, 1., 1.);

  pad_h->Draw();
  pad_r->Draw();

  pad_h->SetLogy();

  pad_h->cd();  //change pad

  v_mll_ee_inc[0]->Draw();

  hs->Draw("histsame");
 
  v_mll_ee_inc[size-1]->Draw("same e1 x0");

  pad_h->RedrawAxis();  //c1->pad_h

  //--------------------------Legend-----------------------------//
  TLegend *leg = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg->SetLineColor(kWhite);
  leg->SetTextFont(42); 
  leg->SetTextSize(0.026);
  leg->SetShadowColor(kWhite); 
  leg->SetFillColor(kWhite);
  leg->AddEntry("h_mll_ee_inc_data","data","ep"); 
  leg->AddEntry("h_mll_ee_inc_zjets","Z+jets","F");
  leg->AddEntry("h_mll_ee_inc_wz","WZ","F");
  leg->AddEntry("h_mll_ee_inc_zz","ZZ","F");
  leg->AddEntry("h_mll_ee_inc_ttbar","ttbar","F");
  leg->AddEntry("h_mll_ee_inc_singlet","single top","F"); 
  leg->AddEntry("h_mll_ee_inc_ww","WW","F");
  leg->AddEntry("h_mll_ee_inc_ttv","ttV","F");
  leg->AddEntry("h_mll_ee_inc_vvv","VVV","F");
  

  leg->Draw();
  //-------------------------------------------------------------//

  //---------------------Ratio Pad-------------------------------//

  pad_r->cd();

  TH1F *h_mll_ee_inc_data_clone = (TH1F*) v_mll_ee_inc[size-1]->Clone("h_mll_ee_inc_data_clone");
  
  //  overflow(sum);

  //Divide
  h_mll_ee_inc_data_clone->Divide(h_mll_ee_inc_data_clone,sum);

  h_mll_ee_inc_data_clone->GetXaxis()->SetRangeUser(50,200); 
  h_mll_ee_inc_data_clone->GetYaxis()->SetRangeUser(0,2);
  h_mll_ee_inc_data_clone->GetYaxis()->SetNdivisions(4);
  h_mll_ee_inc_data_clone->GetYaxis()->SetLabelSize(.12);
  h_mll_ee_inc_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r->SetGridy();

  h_mll_ee_inc_data_clone->Draw();  //"e1"?

  //-------------------------------------------------------------//


  //------------------------MuMu---------------------------------//
  //Canvas
  TCanvas *c2=new TCanvas("c2","mll mumu inc",800,800);
  TPad *pad_h2 = new TPad("pad_h2","Histo Pad2",0., 0, 1., 0.8);
  TPad *pad_r2 = new TPad("pad_r2","Ratio Pad2",0., 0.8, 1., 1.);

  pad_h2->Draw();
  pad_r2->Draw();

  pad_h2->SetLogy();

  pad_h2->cd();  //change pad

  v_mll_mumu_inc[0]->Draw();
  hs2->Draw("histsame");
 
  v_mll_mumu_inc[size-1]->Draw("same e1 x0");

  pad_h2->RedrawAxis(); 


  //-------------------------------Legend-----------------------------------------//
  TLegend *leg2 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg2->SetLineColor(kWhite);
  leg2->SetTextFont(42); 
  leg2->SetTextSize(0.026);
  leg2->SetShadowColor(kWhite); 
  leg2->SetFillColor(kWhite); 
  leg2->AddEntry("h_mll_mumu_inc_data","data","ep");
  leg2->AddEntry("h_mll_mumu_inc_zjets","Z+jets","F");
  leg2->AddEntry("h_mll_mumu_inc_wz","WZ","F");
  leg2->AddEntry("h_mll_mumu_inc_zz","ZZ","F");
  leg2->AddEntry("h_mll_mumu_inc_ttbar","ttbar","F");
  leg2->AddEntry("h_mll_mumu_inc_singlet","single top","F");
  leg2->AddEntry("h_mll_mumu_inc_ww","WW","F");
  leg2->AddEntry("h_mll_mumu_inc_ttv","ttV","F");
  leg2->AddEntry("h_mll_mumu_inc_vvv","VVV","F");

  leg2->Draw();
  //------------------------------------------------------------------------------//

  //---------------------Ratio Pad-------------------------------//

  pad_r2->cd();

  TH1F *h_mll_mumu_inc_data_clone = (TH1F*) v_mll_mumu_inc[size-1]->Clone("h_mll_mumu_inc_data_clone");

  //  overflow(sum2);
  
  //Divide
  h_mll_mumu_inc_data_clone->Divide(h_mll_mumu_inc_data_clone,sum2);

  h_mll_mumu_inc_data_clone->GetXaxis()->SetRangeUser(50,200); 
  h_mll_mumu_inc_data_clone->GetYaxis()->SetRangeUser(0,2);
  h_mll_mumu_inc_data_clone->GetYaxis()->SetNdivisions(4);
  h_mll_mumu_inc_data_clone->GetYaxis()->SetLabelSize(.12);
  h_mll_mumu_inc_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r2->SetGridy();

  h_mll_mumu_inc_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //-------------------Targeted ee-------------------------------//
  
  //Canvas
  TCanvas *c3=new TCanvas("c3","mll ee targeted",800,800);
  TPad *pad_tarh = new TPad("pad_tarh","Targeted Histo Pad",0., 0, 1., 0.8);
  TPad *pad_tarr = new TPad("pad_tarr","Targeted Ratio Pad",0., 0.8, 1., 1.);

  pad_tarh->Draw();
  pad_tarr->Draw();

  pad_tarh->SetLogy();

  pad_tarh->cd();  //change pad

  v_mll_ee_tar[0]->Draw();
  hstar->Draw("histsame");
 
  v_mll_ee_tar[size-1]->Draw("same e1 x0");

  pad_tarh->RedrawAxis();

  //--------------------------Legend-----------------------------//
  TLegend *tarleg = new TLegend(0.78, 0.63, 0.87, 0.89);
  tarleg->SetLineColor(kWhite);
  tarleg->SetTextFont(42); 
  tarleg->SetTextSize(0.026);
  tarleg->SetShadowColor(kWhite); 
  tarleg->SetFillColor(kWhite);
  tarleg->AddEntry("h_mll_ee_tar_data","data","ep"); 
  tarleg->AddEntry("h_mll_ee_tar_zjets","Z+jets","F");
  tarleg->AddEntry("h_mll_ee_tar_wz","WZ","F");
  tarleg->AddEntry("h_mll_ee_tar_zz","ZZ","F");
  tarleg->AddEntry("h_mll_ee_tar_ttbar","ttbar","F");
  tarleg->AddEntry("h_mll_ee_tar_singlet","single top","F"); 
  tarleg->AddEntry("h_mll_ee_tar_ww","WW","F");
  tarleg->AddEntry("h_mll_ee_tar_ttv","ttV","F");
  tarleg->AddEntry("h_mll_ee_tar_vvv","VVV","F");
  
  tarleg->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_tarr->cd();

  TH1F *h_mll_ee_tar_data_clone = (TH1F*) v_mll_ee_tar[size-1]->Clone("h_mll_ee_tar_data_clone");
  
  //  overflow(tarsum);

  //Divide
  h_mll_ee_tar_data_clone->Divide(h_mll_ee_tar_data_clone,tarsum);

  h_mll_ee_tar_data_clone->GetXaxis()->SetRangeUser(50,200); 
  h_mll_ee_tar_data_clone->GetYaxis()->SetRangeUser(0,2);
  h_mll_ee_tar_data_clone->GetYaxis()->SetNdivisions(4);
  h_mll_ee_tar_data_clone->GetYaxis()->SetLabelSize(.12);
  h_mll_ee_tar_data_clone->GetYaxis()->SetTitle("data/MC");
  
  pad_tarr->SetGridy();

  h_mll_ee_tar_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //-------------------Targeted MuMu-----------------------------//
  //Canvas
  TCanvas *c4=new TCanvas("c4","mll mumu targeted",800,800);
  TPad *pad_tarh2 = new TPad("pad_tarh2","Targeted Histo Pad",0., 0, 1., 0.8);
  TPad *pad_tarr2 = new TPad("pad_tarr2","Targeted Ratio Pad",0., 0.8, 1., 1.);

  pad_tarh2->Draw();
  pad_tarr2->Draw();

  pad_tarh2->SetLogy();

  pad_tarh2->cd();  //change pad

  v_mll_mumu_tar[0]->Draw();
  hstar2->Draw("histsame");
 
  v_mll_mumu_tar[size-1]->Draw("same e1 x0");

  pad_tarh2->RedrawAxis();

  //--------------------------Legend-----------------------------//
  TLegend *tarleg2 = new TLegend(0.78, 0.63, 0.87, 0.89);
  tarleg2->SetLineColor(kWhite);
  tarleg2->SetTextFont(42); 
  tarleg2->SetTextSize(0.026);
  tarleg2->SetShadowColor(kWhite); 
  tarleg2->SetFillColor(kWhite);
  tarleg2->AddEntry("h_mll_mumu_tar_data","data","ep"); 
  tarleg2->AddEntry("h_mll_mumu_tar_zjets","Z+jets","F");
  tarleg2->AddEntry("h_mll_mumu_tar_wz","WZ","F");
  tarleg2->AddEntry("h_mll_mumu_tar_zz","ZZ","F");
  tarleg2->AddEntry("h_mll_mumu_tar_ttbar","ttbar","F");
  tarleg2->AddEntry("h_mll_mumu_tar_singlet","single top","F"); 
  tarleg2->AddEntry("h_mll_mumu_tar_ww","WW","F");
  tarleg2->AddEntry("h_mll_mumu_tar_ttv","ttV","F");
  tarleg2->AddEntry("h_mll_mumu_tar_vvv","VVV","F");
  
  tarleg2->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_tarr2->cd();

  TH1F *h_mll_ee_tar_data_clone2 = (TH1F*) v_mll_mumu_tar[size-1]->Clone("h_mll_ee_tar_data_clone2");
  

  //  overflow(tarsum2);

  //Divide
  h_mll_ee_tar_data_clone2->Divide(h_mll_ee_tar_data_clone2,tarsum2);

  h_mll_ee_tar_data_clone2->GetXaxis()->SetRangeUser(50,200); 
  h_mll_ee_tar_data_clone2->GetYaxis()->SetRangeUser(0,2);
  h_mll_ee_tar_data_clone2->GetYaxis()->SetNdivisions(4);
  h_mll_ee_tar_data_clone2->GetYaxis()->SetLabelSize(.12);
  h_mll_ee_tar_data_clone2->GetYaxis()->SetTitle("data/MC");
  
  pad_tarr2->SetGridy();

  h_mll_ee_tar_data_clone2->Draw();

  //-------------------------------------------------------------//
  
  //---------------------MET ee INC------------------------------//
  //Canvas
  TCanvas *c5=new TCanvas("c5","Met ee inc",800,800);
  TPad *pad_h5 = new TPad("pad_h5","Histo Pad5",0., 0, 1., 0.8);
  TPad *pad_r5 = new TPad("pad_r5","Ratio Pad5",0., 0.8, 1., 1.);

  pad_h5->Draw();
  pad_r5->Draw();

  pad_h5->SetLogy();  //c1->pad_h

  pad_h5->cd();  //change pad

  v_met_ee_inc[0]->Draw();
  hs_met_ee_inc->Draw("histsame");
 
  v_met_ee_inc[size-1]->Draw("same e1 x0");

  pad_h5->RedrawAxis();  //c1->pad_h


  //-------------------------------Legend---------------------------------------//
  TLegend *leg5 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg5->SetLineColor(kWhite);
  leg5->SetTextFont(42); 
  leg5->SetTextSize(0.026);
  leg5->SetShadowColor(kWhite); 
  leg5->SetFillColor(kWhite); 
  leg5->AddEntry("h_met_ee_inc_data","data","ep");
  leg5->AddEntry("h_met_ee_inc_zjets","Z+jets","F");
  leg5->AddEntry("h_met_ee_inc_wz","WZ","F");
  leg5->AddEntry("h_met_ee_inc_zz","ZZ","F");
  leg5->AddEntry("h_met_ee_inc_ttbar","ttbar","F");
  leg5->AddEntry("h_met_ee_inc_singlet","single top","F");
  leg5->AddEntry("h_met_ee_inc_ww","WW","F");
  leg5->AddEntry("h_met_ee_inc_ttv","ttV","F");
  leg5->AddEntry("h_met_ee_inc_vvv","VVV","F");

  leg5->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r5->cd();

  TH1F *h_met_ee_inc_data_clone = (TH1F*) v_met_ee_inc[size-1]->Clone("h_met_ee_inc_data_clone");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_met_ee_inc_data_clone->Divide(h_met_ee_inc_data_clone,h_metsum_ee_inc);

  h_met_ee_inc_data_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_ee_inc_data_clone->GetYaxis()->SetRangeUser(0,3);
  h_met_ee_inc_data_clone->GetYaxis()->SetNdivisions(4);
  h_met_ee_inc_data_clone->GetYaxis()->SetLabelSize(.12);
  h_met_ee_inc_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r5->SetGridy();

  h_met_ee_inc_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //---------------------MET MuMu INC----------------------------//
  //Canvas
  TCanvas *c6=new TCanvas("c6","Met mumu inc",800,800);
  TPad *pad_h6 = new TPad("pad_h6","Histo Pad6",0., 0, 1., 0.8);
  TPad *pad_r6 = new TPad("pad_r6","Ratio Pad6",0., 0.8, 1., 1.);

  pad_h6->Draw();
  pad_r6->Draw();

  pad_h6->SetLogy();  //c1->pad_h

  pad_h6->cd();  //change pad

  v_met_mumu_inc[0]->Draw();
  hs_met_mumu_inc->Draw("histsame");
 
  v_met_mumu_inc[size-1]->Draw("same e1 x0");

  pad_h6->RedrawAxis();  //c1->pad_h

  //-------------------------------Legend-------------------------------------//
  TLegend *leg6 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg6->SetLineColor(kWhite);
  leg6->SetTextFont(42); 
  leg6->SetTextSize(0.026);
  leg6->SetShadowColor(kWhite); 
  leg6->SetFillColor(kWhite); 
  leg6->AddEntry("h_met_mumu_inc_data","data","ep");
  leg6->AddEntry("h_met_mumu_inc_zjets","Z+jets","F");
  leg6->AddEntry("h_met_mumu_inc_wz","WZ","F");
  leg6->AddEntry("h_met_mumu_inc_zz","ZZ","F");
  leg6->AddEntry("h_met_mumu_inc_ttbar","ttbar","F");
  leg6->AddEntry("h_met_mumu_inc_singlet","single top","F");
  leg6->AddEntry("h_met_mumu_inc_ww","WW","F");
  leg6->AddEntry("h_met_mumu_inc_ttv","ttV","F");
  leg6->AddEntry("h_met_mumu_inc_vvv","VVV","F");

  leg6->Draw();
  //----------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r6->cd();

  TH1F *h_met_mumu_inc_data_clone = (TH1F*) v_met_mumu_inc[size-1]->Clone("h_met_mumu_inc_data_clone");

  //  overflow(h_metsum_mumu_inc);
  
  //Divide
  h_met_mumu_inc_data_clone->Divide(h_met_mumu_inc_data_clone,h_metsum_mumu_inc);

  h_met_mumu_inc_data_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_mumu_inc_data_clone->GetYaxis()->SetRangeUser(0,3);
  h_met_mumu_inc_data_clone->GetYaxis()->SetNdivisions(4);
  h_met_mumu_inc_data_clone->GetYaxis()->SetLabelSize(.12);
  h_met_mumu_inc_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r6->SetGridy();

  h_met_mumu_inc_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //---------------------MET ee Tar0-----------------------------//
  //Canvas
  TCanvas *c7=new TCanvas("c7","Met ee tar0",800,800);
  TPad *pad_h7 = new TPad("pad_h7","Histo Pad7",0., 0, 1., 0.8);
  TPad *pad_r7 = new TPad("pad_r7","Ratio Pad7",0., 0.8, 1., 1.);

  pad_h7->Draw();
  pad_r7->Draw();

  pad_h7->SetLogy();  //c1->pad_h

  pad_h7->cd();  //change pad

  v_met_ee_tar0[0]->Draw();
  hs_met_ee_tar0->Draw("histsame");
 
  v_met_ee_tar0[size-1]->Draw("same e1 x0");

  pad_h7->RedrawAxis();  //c1->pad_h

  //-------------------------------Legend------------------------------------//
  TLegend *leg7 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg7->SetLineColor(kWhite);
  leg7->SetTextFont(42); 
  leg7->SetTextSize(0.026);
  leg7->SetShadowColor(kWhite); 
  leg7->SetFillColor(kWhite); 
  leg7->AddEntry("h_met_ee_tar0_data","data","ep");
  leg7->AddEntry("h_met_ee_tar0_zjets","Z+jets","F");
  leg7->AddEntry("h_met_ee_tar0_wz","WZ","F");
  leg7->AddEntry("h_met_ee_tar0_zz","ZZ","F");
  leg7->AddEntry("h_met_ee_tar0_ttbar","ttbar","F");
  leg7->AddEntry("h_met_ee_tar0_singlet","single top","F");
  leg7->AddEntry("h_met_ee_tar0_ww","WW","F");
  leg7->AddEntry("h_met_ee_tar0_ttv","ttV","F");
  leg7->AddEntry("h_met_ee_tar0_vvv","VVV","F");

  leg7->Draw();
  //--------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r7->cd();

  TH1F *h_met_ee_tar0_data_clone = (TH1F*) v_met_ee_tar0[size-1]->Clone("h_met_ee_tar0_data_clone");

  //  overflow(h_metsum_ee_tar0);
  
  //Divide
  h_met_ee_tar0_data_clone->Divide(h_met_ee_tar0_data_clone,h_metsum_ee_tar0);

  h_met_ee_tar0_data_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_ee_tar0_data_clone->GetYaxis()->SetRangeUser(0,3);
  h_met_ee_tar0_data_clone->GetYaxis()->SetNdivisions(4);
  h_met_ee_tar0_data_clone->GetYaxis()->SetLabelSize(.12);
  h_met_ee_tar0_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r7->SetGridy();

  h_met_ee_tar0_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //---------------------MET MuMu Tar0---------------------------//
  //Canvas
  TCanvas *c8 =new TCanvas("c8","Met mumu tar0",800,800);
  TPad *pad_h8 = new TPad("pad_h8","Histo Pad8",0., 0, 1., 0.8);
  TPad *pad_r8 = new TPad("pad_r8","Ratio Pad8",0., 0.8, 1., 1.);

  pad_h8->Draw();
  pad_r8->Draw();

  pad_h8->SetLogy();  //c1->pad_h

  pad_h8->cd();  //change pad

  v_met_mumu_tar0[0]->Draw();
  hs_met_mumu_tar0->Draw("histsame");
 
  v_met_mumu_tar0[size-1]->Draw("same e1 x0");

  pad_h8->RedrawAxis();  //c1->pad_h

  //-------------------------------Legend-------------------------------------//
  TLegend *leg8 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg8->SetLineColor(kWhite);
  leg8->SetTextFont(42); 
  leg8->SetTextSize(0.026);
  leg8->SetShadowColor(kWhite); 
  leg8->SetFillColor(kWhite); 
  leg8->AddEntry("h_met_mumu_tar0_data","data","ep");
  leg8->AddEntry("h_met_mumu_tar0_zjets","Z+jets","F");
  leg8->AddEntry("h_met_mumu_tar0_wz","WZ","F");
  leg8->AddEntry("h_met_mumu_tar0_zz","ZZ","F");
  leg8->AddEntry("h_met_mumu_tar0_ttbar","ttbar","F");
  leg8->AddEntry("h_met_mumu_tar0_singlet","single top","F");
  leg8->AddEntry("h_met_mumu_tar0_ww","WW","F");
  leg8->AddEntry("h_met_mumu_tar0_ttv","ttV","F");
  leg8->AddEntry("h_met_mumu_tar0_vvv","VVV","F");

  leg8->Draw();
  //--------------------------------------------------------------------------//
  //--------------------Ratio Pad-------------------------------//

  pad_r8->cd();

  TH1F *h_met_mumu_tar0_data_clone = (TH1F*) v_met_mumu_tar0[size-1]->Clone("h_met_mumu_tar0_data_clone");

  //  overflow(h_metsum_mumu_tar0);
  
  //Divide
  h_met_mumu_tar0_data_clone->Divide(h_met_mumu_tar0_data_clone,h_metsum_mumu_tar0);

  h_met_mumu_tar0_data_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_mumu_tar0_data_clone->GetYaxis()->SetRangeUser(0,3);
  h_met_mumu_tar0_data_clone->GetYaxis()->SetNdivisions(4);
  h_met_mumu_tar0_data_clone->GetYaxis()->SetLabelSize(.12);
  h_met_mumu_tar0_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r8->SetGridy();

  h_met_mumu_tar0_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //---------------------MET ee Tar2------------------------------//
  //Canvas
  TCanvas *c9=new TCanvas("c9","Met ee tar2",800,800);
  TPad *pad_h9 = new TPad("pad_h9","Histo Pad9",0., 0, 1., 0.8);
  TPad *pad_r9 = new TPad("pad_r9","Ratio Pad9",0., 0.8, 1., 1.);

  pad_h9->Draw();
  pad_r9->Draw();

  pad_h9->SetLogy();  //c1->pad_h

  pad_h9->cd();  //change pad

  v_met_ee_tar2[0]->Draw();
  hs_met_ee_tar2->Draw("histsame");
 
  v_met_ee_tar2[size-1]->Draw("same e1 x0");

  pad_h9->RedrawAxis();  //c1->pad_h

  //-------------------------------Legend-------------------------------------//
  TLegend *leg9 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg9->SetLineColor(kWhite);
  leg9->SetTextFont(42); 
  leg9->SetTextSize(0.026);
  leg9->SetShadowColor(kWhite); 
  leg9->SetFillColor(kWhite); 
  leg9->AddEntry("h_met_ee_tar2_data","data","ep");
  leg9->AddEntry("h_met_ee_tar2_zjets","Z+jets","F");
  leg9->AddEntry("h_met_ee_tar2_wz","WZ","F");
  leg9->AddEntry("h_met_ee_tar2_zz","ZZ","F");
  leg9->AddEntry("h_met_ee_tar2_ttbar","ttbar","F");
  leg9->AddEntry("h_met_ee_tar2_singlet","single top","F");
  leg9->AddEntry("h_met_ee_tar2_ww","WW","F");
  leg9->AddEntry("h_met_ee_tar2_ttv","ttV","F");
  leg9->AddEntry("h_met_ee_tar2_vvv","VVV","F");

  leg9->Draw();
  //----------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r9->cd();

  TH1F *h_met_ee_tar2_data_clone = (TH1F*) v_met_ee_tar2[size-1]->Clone("h_met_ee_tar2_data_clone");

  //  overflow(h_metsum_ee_tar2);
  
  //Divide
  h_met_ee_tar2_data_clone->Divide(h_met_ee_tar2_data_clone,h_metsum_ee_tar2);

  h_met_ee_tar2_data_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_ee_tar2_data_clone->GetYaxis()->SetRangeUser(0,3);
  h_met_ee_tar2_data_clone->GetYaxis()->SetNdivisions(4);
  h_met_ee_tar2_data_clone->GetYaxis()->SetLabelSize(.12);
  h_met_ee_tar2_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r9->SetGridy();

  h_met_ee_tar2_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //---------------------MET MuMu Tar2------------------------------//
  //Canvas
  TCanvas *c10 =new TCanvas("c10","Met mumu tar2",800,800);
  TPad *pad_h10 = new TPad("pad_h10","Histo Pad10",0., 0, 1., 0.8);
  TPad *pad_r10 = new TPad("pad_r10","Ratio Pad10",0., 0.8, 1., 1.);

  pad_h10->Draw();
  pad_r10->Draw();

  pad_h10->SetLogy();  //c1->pad_h

  pad_h10->cd();  //change pad

  v_met_mumu_tar2[0]->Draw();
  hs_met_mumu_tar2->Draw("histsame");
 
  v_met_mumu_tar2[size-1]->Draw("same e1 x0");

  pad_h10->RedrawAxis();  //c1->pad_h

  //-------------------------------Legend-----------------------------------------//
  TLegend *leg10 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg10->SetLineColor(kWhite);
  leg10->SetTextFont(42); 
  leg10->SetTextSize(0.026);
  leg10->SetShadowColor(kWhite); 
  leg10->SetFillColor(kWhite); 
  leg10->AddEntry("h_met_mumu_tar2_data","data","ep");
  leg10->AddEntry("h_met_mumu_tar2_zjets","Z+jets","F");
  leg10->AddEntry("h_met_mumu_tar2_wz","WZ","F");
  leg10->AddEntry("h_met_mumu_tar2_zz","ZZ","F");
  leg10->AddEntry("h_met_mumu_tar2_ttbar","ttbar","F");
  leg10->AddEntry("h_met_mumu_tar2_singlet","single top","F");
  leg10->AddEntry("h_met_mumu_tar2_ww","WW","F");
  leg10->AddEntry("h_met_mumu_tar2_ttv","ttV","F");
  leg10->AddEntry("h_met_mumu_tar2_vvv","VVV","F");

  leg10->Draw();
  //----------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r10->cd();

  TH1F *h_met_mumu_tar2_data_clone = (TH1F*) v_met_mumu_tar2[size-1]->Clone("h_met_mumu_tar2_data_clone");

  //  overflow(h_metsum_mumu_tar2);;
  
  //Divide
  h_met_mumu_tar2_data_clone->Divide(h_met_mumu_tar2_data_clone,h_metsum_mumu_tar2);

  h_met_mumu_tar2_data_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_mumu_tar2_data_clone->GetYaxis()->SetRangeUser(0,3);
  h_met_mumu_tar2_data_clone->GetYaxis()->SetNdivisions(4);
  h_met_mumu_tar2_data_clone->GetYaxis()->SetLabelSize(.12);
  h_met_mumu_tar2_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r10->SetGridy();

  h_met_mumu_tar2_data_clone->Draw();

  //-------------------------------------------------------------//
  
  //----------------------------------------------------------------------------//
  //-----------------------------Overlay Plots----------------------------------//
  //----------------------------------------------------------------------------//
  
  //---------------------MET Tar0 background---------------------------//
  //Canvas
  TCanvas *c11=new TCanvas("c11","Met ee vs mumu tar0 bkgd",800,800);
  TPad *pad_h11 = new TPad("pad_h11","Histo Pad11",0., 0, 1., 0.8);
  TPad *pad_r11 = new TPad("pad_r11","Ratio Pad11",0., 0.8, 1., 1.);

  pad_h11->Draw();
  pad_r11->Draw();

  pad_h11->SetLogy();

  pad_h11->cd();  

  //THStack *hs_tar0_MC_comparison = new THStack("hs_tar0_MC_comparison","hs_tar0_MC_comparison");  
  
  h_metsum_ee_tar0->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_mumu_tar0->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_mumu_tar0->SetLineColor(kRed);
  h_metsum_ee_tar0->SetLineColor(kBlue);

  //hs_tar0_MC_comparison->Add(h_metsum_ee_tar0);
  //hs_tar0_MC_comparison->Add(h_metsum_mumu_tar0);
  h_metsum_mumu_tar0->Draw("e1");
  //hs_tar0_MC_comparison->Draw("histsame");
  h_metsum_ee_tar0->Draw("same e1");

  pad_h11->RedrawAxis();
  
  //-------------------------------Legend------------------------------------//
  TLegend *leg11 = new TLegend(0.78, 0.75, 0.87, 0.83);
  leg11->SetLineColor(kWhite);
  leg11->SetTextFont(42); 
  leg11->SetTextSize(0.026);
  leg11->SetShadowColor(kWhite); 
  leg11->SetFillColor(kWhite); 
  leg11->AddEntry("h_metsum_ee_tar0","ee MC","l e");
  leg11->AddEntry("h_metsum_mumu_tar0","mumu MC","e l");


  leg11->Draw();
  
  //--------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//
  
  pad_r11->cd();

  TH1F *h_metsum_ee_tar0_clone = (TH1F*) h_metsum_ee_tar0->Clone("h_metsum_ee_tar0_clone");

  //  overflow(h_metsum_ee_tar0);
  
  //Divide
  h_metsum_ee_tar0_clone->Divide(h_metsum_ee_tar0_clone,h_metsum_mumu_tar0);

  h_metsum_ee_tar0_clone->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_ee_tar0_clone->GetYaxis()->SetRangeUser(0,2);
  h_metsum_ee_tar0_clone->GetYaxis()->SetNdivisions(4);
  h_metsum_ee_tar0_clone->GetYaxis()->SetLabelSize(.12);
  h_metsum_ee_tar0_clone->GetYaxis()->SetTitle("ee/mumu MC");

  pad_r11->SetGridy();

  h_metsum_ee_tar0_clone->Draw();
  
  //-------------------------------------------------------------//

  //---------------------MET Tar0 Data---------------------------//
  //Canvas
  TCanvas *c12=new TCanvas("c12","Met ee vs mumu tar0 data",800,800);
  TPad *pad_h12 = new TPad("pad_h12","Histo Pad12",0., 0, 1., 0.8);
  TPad *pad_r12 = new TPad("pad_r12","Ratio Pad12",0., 0.8, 1., 1.);

  pad_h12->Draw();
  pad_r12->Draw();

  pad_h12->SetLogy();  //c1->pad_h

  pad_h12->cd();  //change pad

  //  THStack *hs_tar0_data_comparison = new THStack("hs_tar0_data_comparison","hs_tar0_data_comparison");  
  
  v_met_ee_tar0[size-1]->SetLineColor(kBlue); 
  v_met_mumu_tar0[size-1]->SetLineColor(kRed); 
  v_met_ee_tar0[size-1]->GetXaxis()->SetRangeUser(0,200); 
  v_met_mumu_tar0[size-1]->GetXaxis()->SetRangeUser(0,200); 

  //  hs_tar0_data_comparison->Add(v_met_ee_tar0[size-1]);
  //  hs_tar0_data_comparison->Add(v_met_mumu_tar0[size-1]);
  v_met_mumu_tar0[size-1]->Draw();
  v_met_ee_tar0[size-1]->Draw("same e1");
  //  hs_tar0_data_comparison->Draw("histsame");

  pad_h12->RedrawAxis();
  
  //-------------------------------Legend------------------------------------//
  TLegend *leg12 = new TLegend(0.78, 0.75, 0.87, 0.83);
  leg12->SetLineColor(kWhite);
  leg12->SetTextFont(42); 
  leg12->SetTextSize(0.026);
  leg12->SetShadowColor(kWhite); 
  leg12->SetFillColor(kWhite); 
  leg12->AddEntry("h_met_ee_tar0_data","ee data","l e");
  leg12->AddEntry("h_met_mumu_tar0_data","mumu data","l e");

  leg12->Draw();
  
  //--------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//
  
  pad_r12->cd();

  TH1F *h_met_ee_tar0_datacomp_clone = (TH1F*) v_met_ee_tar0[size-1]->Clone("h_met_ee_tar0_datacomp_clone");

  //  overflow(h_metsum_ee_tar0);
  
  //Divide
  h_met_ee_tar0_datacomp_clone->Divide(h_met_ee_tar0_datacomp_clone,v_met_mumu_tar0[size-1]);

  h_met_ee_tar0_datacomp_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_ee_tar0_datacomp_clone->GetYaxis()->SetRangeUser(0,2);
  h_met_ee_tar0_datacomp_clone->GetYaxis()->SetNdivisions(4);
  h_met_ee_tar0_datacomp_clone->GetYaxis()->SetLabelSize(.12);
  h_met_ee_tar0_datacomp_clone->GetYaxis()->SetTitle("ee/mumu MC");

  pad_r12->SetGridy();

  h_met_ee_tar0_datacomp_clone->Draw();
  
  //-------------------------------------------------------------//

  //------------------------------------------------------------------------------//  
   
  //---------------------MET Tar2 background---------------------------//
  //Canvas
  TCanvas *c13=new TCanvas("c13","Met ee vs mumu tar2 bkgd",800,800);
  TPad *pad_h13 = new TPad("pad_h13","Histo Pad13",0., 0, 1., 0.8);
  TPad *pad_r13 = new TPad("pad_r13","Ratio Pad13",0., 0.8, 1., 1.);

  pad_h13->Draw();
  pad_r13->Draw();

  pad_h13->SetLogy();  //c1->pad_h

  pad_h13->cd();  //change pad

  //  THStack *hs_tar2_MC_comparison = new THStack("hs_tar2_MC_comparison","hs_tar2_MC_comparison");  
  
  h_metsum_ee_tar2->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_mumu_tar2->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_mumu_tar2->SetLineColor(kRed);
  h_metsum_ee_tar2->SetLineColor(kBlue);

  //  hs_tar2_MC_comparison->Add(h_metsum_ee_tar2);
  // hs_tar2_MC_comparison->Add(h_metsum_mumu_tar2);
  h_metsum_mumu_tar2->Draw();
  h_metsum_ee_tar2->Draw("same e1");
  //  hs_tar2_MC_comparison->Draw("histsame");

  pad_h13->RedrawAxis();
  
  //-------------------------------Legend------------------------------------//
  TLegend *leg13 = new TLegend(0.78, 0.75, 0.87, 0.83);
  leg13->SetLineColor(kWhite);
  leg13->SetTextFont(42); 
  leg13->SetTextSize(0.026);
  leg13->SetShadowColor(kWhite); 
  leg13->SetFillColor(kWhite); 
  leg13->AddEntry("h_metsum_ee_tar2","ee MC","l e");
  leg13->AddEntry("h_metsum_mumu_tar2","mumu MC","l e");

  leg13->Draw();
  
  //--------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//
  
  pad_r13->cd();

  TH1F *h_metsum_ee_tar2_clone = (TH1F*) h_metsum_ee_tar2->Clone("h_metsum_ee_tar2_clone");

  //  overflow(h_metsum_ee_tar2);
  
  //Divide
  h_metsum_ee_tar2_clone->Divide(h_metsum_ee_tar2_clone,h_metsum_mumu_tar2);

  h_metsum_ee_tar2_clone->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_ee_tar2_clone->GetYaxis()->SetRangeUser(0,2);
  h_metsum_ee_tar2_clone->GetYaxis()->SetNdivisions(4);
  h_metsum_ee_tar2_clone->GetYaxis()->SetLabelSize(.12);
  h_metsum_ee_tar2_clone->GetYaxis()->SetTitle("ee/mumu MC");

  pad_r13->SetGridy();

  h_metsum_ee_tar2_clone->Draw();
  
  //-------------------------------------------------------------//
  
  //---------------------MET Tar2 Data---------------------------//
  //Canvas
  TCanvas *c14=new TCanvas("c14","Met ee vs mumu tar2 data",800,800);
  TPad *pad_h14 = new TPad("pad_h14","Histo Pad14",0., 0, 1., 0.8);
  TPad *pad_r14 = new TPad("pad_r14","Ratio Pad14",0., 0.8, 1., 1.);

  pad_h14->Draw();
  pad_r14->Draw();

  pad_h14->SetLogy();  //c1->pad_h

  pad_h14->cd();  //change pad

  //  THStack *hs_tar2_data_comparison = new THStack("hs_tar2_data_comparison","hs_tar2_data_comparison");  
  
  v_met_ee_tar2[size-1]->SetLineColor(kBlue); 
  v_met_mumu_tar2[size-1]->SetLineColor(kRed); 
  v_met_ee_tar2[size-1]->GetXaxis()->SetRangeUser(0,200); 
  v_met_mumu_tar2[size-1]->GetXaxis()->SetRangeUser(0,200); 

  //  hs_tar2_data_comparison->Add(v_met_ee_tar2[size-1]);
  //  hs_tar2_data_comparison->Add(v_met_mumu_tar2[size-1]);
  v_met_mumu_tar2[size-1]->Draw();
  v_met_ee_tar2[size-1]->Draw("same e1");
  //  hs_tar2_data_comparison->Draw("histsame");

  pad_h14->RedrawAxis();
  
  //-------------------------------Legend------------------------------------//
  TLegend *leg14 = new TLegend(0.78, 0.75, 0.87, 0.83);
  leg14->SetLineColor(kWhite);
  leg14->SetTextFont(42); 
  leg14->SetTextSize(0.026);
  leg14->SetShadowColor(kWhite); 
  leg14->SetFillColor(kWhite); 
  leg14->AddEntry("h_met_ee_tar2_data","ee data","l e");
  leg14->AddEntry("h_met_mumu_tar2_data","mumu data","l e");

  leg14->Draw();
  
  //--------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//
  
  pad_r14->cd();

  TH1F *h_met_ee_tar2_datacomp_clone = (TH1F*) v_met_ee_tar2[size-1]->Clone("h_met_ee_tar2_datacomp_clone");

  //  overflow(h_metsum_ee_tar2);
  
  //Divide
  h_met_ee_tar2_datacomp_clone->Divide(h_met_ee_tar2_datacomp_clone,v_met_mumu_tar2[size-1]);

  h_met_ee_tar2_datacomp_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_ee_tar2_datacomp_clone->GetYaxis()->SetRangeUser(0,2);
  h_met_ee_tar2_datacomp_clone->GetYaxis()->SetNdivisions(4);
  h_met_ee_tar2_datacomp_clone->GetYaxis()->SetLabelSize(.12);
  h_met_ee_tar2_datacomp_clone->GetYaxis()->SetTitle("ee/mumu MC");

  pad_r14->SetGridy();

  h_met_ee_tar2_datacomp_clone->Draw();
  
  //-------------------------------------------------------------//

  //---------------------------------------------------------------------------------------------------------------//

  //---------------------MET Inc background---------------------------//
  //Canvas
  TCanvas *c15=new TCanvas("c15","Met ee vs mumu inc bkgd",800,800);
  TPad *pad_h15 = new TPad("pad_h15","Histo Pad15",0., 0, 1., 0.8);
  TPad *pad_r15 = new TPad("pad_r15","Ratio Pad15",0., 0.8, 1., 1.);

  pad_h15->Draw();
  pad_r15->Draw();

  pad_h15->SetLogy();  //c1->pad_h

  pad_h15->cd();  //change pad

  //  THStack *hs_inc_MC_comparison = new THStack("hs_inc_MC_comparison","hs_inc_MC_comparison");  
  
  h_metsum_ee_inc->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_mumu_inc->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_mumu_inc->SetLineColor(kRed);
  h_metsum_ee_inc->SetLineColor(kBlue);

  //  hs_inc_MC_comparison->Add(h_metsum_ee_inc);
  //  hs_inc_MC_comparison->Add(h_metsum_mumu_inc);
  h_metsum_mumu_inc->Draw();
  h_metsum_ee_inc->Draw("same e1");
  //  hs_inc_MC_comparison->Draw("histsame");

  pad_h15->RedrawAxis();
  
  //-------------------------------Legend------------------------------------//
  TLegend *leg15 = new TLegend(0.78, 0.75, 0.87, 0.83);
  leg15->SetLineColor(kWhite);
  leg15->SetTextFont(42); 
  leg15->SetTextSize(0.026);
  leg15->SetShadowColor(kWhite); 
  leg15->SetFillColor(kWhite); 
  leg15->AddEntry("h_metsum_ee_inc","ee MC","l e");
  leg15->AddEntry("h_metsum_mumu_inc","mumu MC","l e");

  leg15->Draw();
  
  //--------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//
  
  pad_r15->cd();

  TH1F *h_metsum_ee_inc_clone = (TH1F*) h_metsum_ee_inc->Clone("h_metsum_ee_inc_clone");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_metsum_ee_inc_clone->Divide(h_metsum_ee_inc_clone,h_metsum_mumu_inc);

  h_metsum_ee_inc_clone->GetXaxis()->SetRangeUser(0,200); 
  h_metsum_ee_inc_clone->GetYaxis()->SetRangeUser(0,2);
  h_metsum_ee_inc_clone->GetYaxis()->SetNdivisions(4);
  h_metsum_ee_inc_clone->GetYaxis()->SetLabelSize(.12);
  h_metsum_ee_inc_clone->GetYaxis()->SetTitle("ee/mumu MC");

  pad_r15->SetGridy();

  h_metsum_ee_inc_clone->Draw();
  
  //-------------------------------------------------------------//

  //---------------------MET Inc Data---------------------------//
  //Canvas
  TCanvas *c16=new TCanvas("c16","Met ee vs mumu inc data",800,800);
  TPad *pad_h16 = new TPad("pad_h16","Histo Pad16",0., 0, 1., 0.8);
  TPad *pad_r16 = new TPad("pad_r16","Ratio Pad16",0., 0.8, 1., 1.);

  pad_h16->Draw();
  pad_r16->Draw();

  pad_h16->SetLogy();  //c1->pad_h

  pad_h16->cd();  //change pad

  //  THStack *hs_inc_data_comparison = new THStack("hs_inc_data_comparison","hs_inc_data_comparison");  
  
  v_met_ee_inc[size-1]->SetLineColor(kBlue); 
  v_met_mumu_inc[size-1]->SetLineColor(kRed); 
  v_met_ee_inc[size-1]->GetXaxis()->SetRangeUser(0,200); 
  v_met_mumu_inc[size-1]->GetXaxis()->SetRangeUser(0,200); 

  //  hs_inc_data_comparison->Add(v_met_ee_inc[size-1]);
  //  hs_inc_data_comparison->Add(v_met_mumu_inc[size-1]);
  v_met_mumu_inc[size-1]->Draw();
  v_met_ee_inc[size-1]->Draw("same e1");
  //  hs_inc_data_comparison->Draw("histsame");

  pad_h16->RedrawAxis();
  
  //-------------------------------Legend------------------------------------//
  TLegend *leg16 = new TLegend(0.78, 0.75, 0.87, 0.83);
  leg16->SetLineColor(kWhite);
  leg16->SetTextFont(42); 
  leg16->SetTextSize(0.026);
  leg16->SetShadowColor(kWhite); 
  leg16->SetFillColor(kWhite); 
  leg16->AddEntry("h_met_ee_inc_data","ee data","e l");
  leg16->AddEntry("h_met_mumu_inc_data","mumu data","e ;");

  leg16->Draw();
  
  //--------------------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//
  
  pad_r16->cd();

  TH1F *h_met_ee_inc_datacomp_clone = (TH1F*) v_met_ee_inc[size-1]->Clone("h_met_ee_inc_datacomp_clone");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_met_ee_inc_datacomp_clone->Divide(h_met_ee_inc_datacomp_clone,v_met_mumu_inc[size-1]);

  h_met_ee_inc_datacomp_clone->GetXaxis()->SetRangeUser(0,200); 
  h_met_ee_inc_datacomp_clone->GetYaxis()->SetRangeUser(0,2);
  h_met_ee_inc_datacomp_clone->GetYaxis()->SetNdivisions(4);
  h_met_ee_inc_datacomp_clone->GetYaxis()->SetLabelSize(.12);
  h_met_ee_inc_datacomp_clone->GetYaxis()->SetTitle("ee/mumu MC");

  pad_r16->SetGridy();

  h_met_ee_inc_datacomp_clone->Draw();
  
  //-------------------------------------------------------------//



  //---------------------------------------//
  /*
  c5->SaveAs("./3lep_pics/ee_inc.jpg");
  c6->SaveAs("./3lep_pics/mumu_inc.jpg");
  c7->SaveAs("./3lep_pics/ee_tar0.jpg");
  c8->SaveAs("./3lep_pics/mumu_tar0.jpg");
  c9->SaveAs("./3lep_pics/ee_tar2.jpg");
  c10->SaveAs("./3lep_pics/mumu_tar2.jpg");
  c11->SaveAs("./3lep_pics/ee_vs_mumu_tar0_MC.jpg");
  c12->SaveAs("./3lep_pics/ee_vs_mumu_tar0_data.jpg");
  c13->SaveAs("./3lep_pics/ee_vs_mumu_tar2_MC.jpg");
  c14->SaveAs("./3lep_pics/ee_vs_mumu_tar2_data.jpg");
  c15->SaveAs("./3lep_pics/ee_vs_mumu_inc_MC.jpg");
  c16->SaveAs("./3lep_pics/ee_vs_mumu_inc_data.jpg");
  */
  
  return 0;
}
 
