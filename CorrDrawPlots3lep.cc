//  object#:   no number -> ee, 2 -> MuMu
//
//  e.g. samplehisto -> histogram of ee dilepton mass
//       samplehisto2 -> histogram of MuMu dilepton mass

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "THStack.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TFitResultPtr.h"
#include "TMatrixDSym.h"

using namespace std;

//----------Global Variables---------------//
int mll_minbin=0;  //  Why not 50?
int mll_maxbin=200;
int mll_binsize=5;
int met_minbin=0;  //
int met_maxbin=200;
int met_binsize=2;
int corecutoff=60;
int sumet_minbin=0;
int sumet_maxbin=2500;
int sumet_binsize=5;
//-----------------------------------------//


void h_format(TH1F *histo)
{

  TString name=histo->GetName();
  Color_t color=kMagenta+4;

  if(name.Contains("tar") && name.Contains("mll"))
	{histo->Rebin(5);}
  else if (name.Contains("inc") && name.Contains("mll"))
	{histo->Rebin(5);}
  else if (name.Contains("met") && !name.Contains("sum"))
	{histo->Rebin(2);
	  histo->SetMarkerSize(0.5);}
  else if (name.Contains("met") && name.Contains("sum"))
	{histo->Rebin(5);
	  histo->SetMarkerSize(0.5);}

  if(name.Contains("zjets"))     //long because first, needs more formatting
	{
	  color=kCyan;
	  histo->SetFillColor(color);

	  if(name.Contains("mll"))
		{ 
		  histo->GetXaxis()->SetRangeUser(50,200); 
		  histo->GetYaxis()->SetRangeUser(1e0,5e7);//1e1, 5e8

		  if(name.Contains("ee"))
			{
			  histo->GetXaxis()->SetTitle("m_{ll} (GeV)"); 
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
			  histo->GetXaxis()->SetTitle("met_{ll} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 2 GeV");
			}
		  else if(name.Contains("mumu"))
			{
			  histo->GetXaxis()->SetTitle("met_{mumu} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 2 GeV");
			}
		}
	  if(name.Contains("sumet"))
		{  
		  histo->GetXaxis()->SetRangeUser(0,2500);  //350
		  histo->GetYaxis()->SetRangeUser(1e-1,1e5);//1e0,5e7
		 
		  if(name.Contains("ee"))
			{
			  histo->GetXaxis()->SetTitle("summet_{ll} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 2 GeV");
			}
		  else if(name.Contains("mumu"))
			{
			  histo->GetXaxis()->SetTitle("summet_{mumu} (GeV)"); 
			  histo->GetYaxis()->SetTitle("Events / 2 GeV");
			}
		}
	  else if(name.Contains("nvtx"))
		{
		  histo->GetXaxis()->SetRangeUser(0,35);  //350
		  histo->GetYaxis()->SetRangeUser(0,10e6);//1e0,5e7
		  histo->GetXaxis()->SetTitle("number of vertices"); 
		  histo->GetYaxis()->SetTitle("Events");
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


void overflow(TH1F *histo)
{
  double uncertainty=0.;
  int overflowbin=-1;
  TString name=histo->GetName();

  if(name.Contains("mll"))
	{overflowbin=((mll_maxbin - mll_minbin)/mll_binsize);}
  else if(name.Contains("met") && !name.Contains("sum"))
	{overflowbin=((met_maxbin - met_minbin)/met_binsize)+1;}
  // else if(name.Contains("met") && name.Contains("sum"))
  //	{overflowbin=((sumet_maxbin - sumet_minbin)/sumet_binsize)+1;}

  histo->SetBinContent(overflowbin,histo->GetBinContent(overflowbin)+histo->IntegralAndError(overflowbin,-1,uncertainty));

  histo->SetBinError(overflowbin,sqrt(pow(histo->GetBinError((overflowbin)),2) + pow(uncertainty,2)));
}


void tails(TH1D *histo)
{
  int min=histo->FindBin(-corecutoff);
  int max=histo->FindBin(corecutoff);
  for(int i = min; i<= max; i++)
	{
	  histo->SetBinContent(i,0.);
	  histo->SetBinError(i,0.);
	}
}

void core(TH1D *histo)
{
  int min=histo->FindBin(-corecutoff);
  int max=histo->FindBin(corecutoff);
  int last=histo->GetNbinsX();
  for(int i = 1; i< min; i++)
	{
	  histo->SetBinContent(i,0.);
	  histo->SetBinError(i,0.);
	}
  for(int i = max + 1; i<= last ; i++)
	{
	  histo->SetBinContent(i,0.);
	  histo->SetBinError(i,0.);
	}
}


/*
double fitf(double *v, double *par)
{
  double arg = 0;

  if (par[2] != 0)
	{
	  arg = (v[0]-par[1])/par[2];;
	}

  double fitval = par[0]*exp(-0.5*arg*arg);

  return fitval;
}

void fit(TH1D *histo)
{
  
  TF1 *func = new TF1("func", fitf,-350,350,3);  //func-fit
  func->SetParameters(histo->GetMaximum(), histo->GetMean(), histo->GetRMS());
  func->SetParNames("Constant","Mean_value", "RMS_width");
  func->SetLineColor(kBlack);
  histo->Fit("func","0RS");//"0"  func->fit

  cout<<"\nChi-square = "<<func->GetChisquare()<<",    Prob = "<<func->GetProb()<<endl;

  func->DrawClone("Same");
}
*/

void betterfit(TH1D *histo)
{
  TFitResultPtr fit_data = (TFitResultPtr) histo->Fit("gaus","0S");
  TF1 *fit_func = (TF1*) histo->GetFunction("gaus");
  int npar=3;
  cout<<"\nChisquare = "<<fit_func->GetChisquare()<<endl;
}

using std::pair;
pair <double, double> getweightedavg(double x, double dx, double y, double dy)
{
  double val = (x/pow(dx,2)+y/pow(dy,2))/(1/pow(dx,2)+1/pow(dy,2));
  double err = 1/sqrt(1/pow(dx,2)+1/pow(dy,2));
  return make_pair(val, err);
}

void sumparameters(TH1D *core, TH1D *tails)
{
  TF1 *coref = (TF1*) core->GetFunction("gaus");
  TF1 *tailsf = (TF1*) tails->GetFunction("gaus");
  pair<double, double> constant = getweightedavg( coref->GetParameter(0),coref->GetParError(0),tailsf->GetParameter(0), tailsf->GetParError(0));
  pair<double, double> mean  = getweightedavg( coref->GetParameter(1),coref->GetParError(1),tailsf->GetParameter(1), tailsf->GetParError(1));
  pair<double, double> sigma = getweightedavg( coref->GetParameter(2),coref->GetParError(2),tailsf->GetParameter(2), tailsf->GetParError(2));

  int w=11;
  cout<<setw(30)<<"\nWeighted Sum Parameters"<<endl;
  cout<<scientific<<setw(w)<<"\nConstant = "<<setw(w)<<constant.first<<" +- "<<constant.second<<endl;
  cout<<scientific<<setw(w)<<"Mean = "<<setw(w)<<mean.first<<" +- "<<mean.second<<endl;
  cout<<scientific<<setw(w)<<"Sigma = "<<setw(w)<<sigma.first<<" +- "<<sigma.second<<endl;

  TF1 *sum = new TF1("sum","[0]*exp(-0.5*pow( (x-[1])/[2] ,2))",-350,350);
  sum->SetParameter(0, constant.first);
  sum->SetParameter(1, mean.first);
  sum->SetParameter(2, sigma.first);
  sum->SetLineColor(kBlack);
  sum->Draw("same");
}

void all_fitting(TH1D *histo)
{
  TH1D *histo_core = (TH1D*) histo->Clone("histo_core");
  core(histo_core);
  TH1D *histo_tails = (TH1D*) histo->Clone("histo_tails");
  tails(histo_tails);

  TString name=histo->GetName();

  cout<<"\n"<<setw(50)<<name<<endl;
  betterfit(histo);

  cout<<"\nCore"<<endl;
  betterfit(histo_core);

  cout<<"\nTails"<<endl;
  betterfit(histo_tails);

  sumparameters(histo_core, histo_tails);

  //Save to root file.  Change name depending on name of histogram passed to function
  TFile *OutputFile = new TFile(name + ".root","recreate");
  OutputFile->cd();
  histo->Write();
  histo_core->Write();
  histo_tails->Write();
  OutputFile->Close();
}

double find_width(TH1F *histo)
{
  int maxbin = 0;
  maxbin = histo->GetMaximumBin();
  double max = 0;
  max = histo->GetBinContent(maxbin);
  double half_max = 0;
  half_max = pow( 10, 0.5*TMath::Log10(max) );

  int half_maxbin = 0;
  int half_maxbin_low = 0;  
  int half_maxbin_high = 0;  
  int percent_low = 1.00;  //determine if bin content is close enought to half_max
  int percent_high = 0.95;  //determine if bin content is close enought to half_max


  int i = maxbin;
  bool again = true;
  while(again && i<=histo->GetNbinsX())
	{
	  if(TMath::Abs( half_max - histo->GetBinContent(i) ) < percent_low*half_max)
		{
		  half_maxbin_low = i;
		  half_maxbin_high = i;
		  again=false;
		  i++;
		  while(TMath::Abs( half_max - histo->GetBinContent(i) ) < percent_high*half_max && i<=histo->GetNbinsX())
			{half_maxbin_high=i;
			  i++;}
		}
	  else
		{
		  i++;
		}
	}
  
  half_maxbin = (half_maxbin_high + half_maxbin_low)/2;
  double width = 0;
  width = (half_maxbin - maxbin)*histo->GetBinWidth(maxbin);
  cout<<"maxbin = "<<maxbin<<endl;
  cout<<"max = "<<max<<endl;
  cout<<"half_max = "<<half_max<<endl;
  cout<<"half_maxbin_low = "<<half_maxbin_low<<endl;
  cout<<"half_maxbin_high = "<<half_maxbin_high<<endl;
  cout<<"half_maxbin = "<<half_maxbin<<endl;
  return width;
}

//----------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------//

int CorrDrawPlots3lep()
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
  TH1F *h_nvtx_scaled_zjets = (TH1F*) InputFile->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_zjets");
  TH1F *h_nvtx_unscaled_zjets = (TH1F*) InputFile->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_zjets");
  TH2F *h_met_ll_inc_zjets = (TH2F*) InputFile->Get("h_met_ll_inc")->Clone("h_met_ll_inc_zjets");
  TH2F *h_met_ll_tar0_zjets = (TH2F*) InputFile->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_zjets");
  TH2F *h_met_ll_tar2_zjets = (TH2F*) InputFile->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_zjets");
  TH2F *h_met_emu2_inc_zjets = (TH2F*) InputFile->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_zjets");
  TH2F *h_met_emu2_tar0_zjets = (TH2F*) InputFile->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_zjets");
  TH2F *h_met_emu2_tar2_zjets = (TH2F*) InputFile->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_zjets");
  TH1F *h_phi_zjets = (TH1F*) InputFile->Get("h_phi")->Clone("h_phi_zjets");
  TH1F *h_sumet_ee_inc_zjets = (TH1F*) InputFile->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_zjets");
  TH1F *h_sumet_mumu_inc_zjets = (TH1F*) InputFile->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_zjets");
  TH1F *h_sumet_emu_inc_zjets = (TH1F*) InputFile->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_zjets");

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
  TH1F *h_nvtx_scaled_singlet = (TH1F*) InputFile_t->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_singlet");
  TH1F *h_nvtx_unscaled_singlet = (TH1F*) InputFile_t->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_singlet");
  TH2F *h_met_ll_inc_singlet = (TH2F*) InputFile_t->Get("h_met_ll_inc")->Clone("h_met_ll_inc_singlet");
  TH2F *h_met_ll_tar0_singlet = (TH2F*) InputFile_t->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_singlet");
  TH2F *h_met_ll_tar2_singlet = (TH2F*) InputFile_t->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_singlet");
  TH2F *h_met_emu2_inc_singlet = (TH2F*) InputFile_t->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_singlet");
  TH2F *h_met_emu2_tar0_singlet = (TH2F*) InputFile_t->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_singlet");
  TH2F *h_met_emu2_tar2_singlet = (TH2F*) InputFile_t->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_singlet");
  TH1F *h_phi_singlet = (TH1F*) InputFile_t->Get("h_phi")->Clone("h_phi_singlet");
  TH1F *h_sumet_ee_inc_singlet = (TH1F*) InputFile_t->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_singlet");
  TH1F *h_sumet_mumu_inc_singlet = (TH1F*) InputFile_t->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_singlet");
  TH1F *h_sumet_emu_inc_singlet = (TH1F*) InputFile_t->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_singlet");

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
  TH1F *h_nvtx_scaled_wz = (TH1F*) InputFile_wz->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_wz");
  TH1F *h_nvtx_unscaled_wz = (TH1F*) InputFile_wz->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_wz");
  TH2F *h_met_ll_inc_wz = (TH2F*) InputFile_wz->Get("h_met_ll_inc")->Clone("h_met_ll_inc_wz");
  TH2F *h_met_ll_tar0_wz = (TH2F*) InputFile_wz->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_wz");
  TH2F *h_met_ll_tar2_wz = (TH2F*) InputFile_wz->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_wz");
  TH2F *h_met_emu2_inc_wz = (TH2F*) InputFile_wz->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_wz");
  TH2F *h_met_emu2_tar0_wz = (TH2F*) InputFile_wz->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_wz");
  TH2F *h_met_emu2_tar2_wz = (TH2F*) InputFile_wz->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_wz");
  TH1F *h_phi_wz = (TH1F*) InputFile_wz->Get("h_phi")->Clone("h_phi_wz");
  TH1F *h_sumet_ee_inc_wz = (TH1F*) InputFile_wz->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_wz");
  TH1F *h_sumet_mumu_inc_wz = (TH1F*) InputFile_wz->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_wz");
  TH1F *h_sumet_emu_inc_wz = (TH1F*) InputFile_wz->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_wz");

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
  TH1F *h_nvtx_scaled_zz = (TH1F*) InputFile_zz->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_zz");
  TH1F *h_nvtx_unscaled_zz = (TH1F*) InputFile_zz->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_zz");
  TH2F *h_met_ll_inc_zz = (TH2F*) InputFile_zz->Get("h_met_ll_inc")->Clone("h_met_ll_inc_zz");
  TH2F *h_met_ll_tar0_zz = (TH2F*) InputFile_zz->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_zz");
  TH2F *h_met_ll_tar2_zz = (TH2F*) InputFile_zz->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_zz");
  TH2F *h_met_emu2_inc_zz = (TH2F*) InputFile_zz->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_zz");
  TH2F *h_met_emu2_tar0_zz = (TH2F*) InputFile_zz->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_zz");
  TH2F *h_met_emu2_tar2_zz = (TH2F*) InputFile_zz->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_zz");
  TH1F *h_phi_zz = (TH1F*) InputFile_zz->Get("h_phi")->Clone("h_phi_zz");
  TH1F *h_sumet_ee_inc_zz = (TH1F*) InputFile_zz->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_zz");
  TH1F *h_sumet_mumu_inc_zz = (TH1F*) InputFile_zz->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_zz");
  TH1F *h_sumet_emu_inc_zz = (TH1F*) InputFile_zz->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_zz");

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
  TH1F *h_nvtx_scaled_ttbar = (TH1F*) InputFile_ttbar->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_ttbar");
  TH1F *h_nvtx_unscaled_ttbar = (TH1F*) InputFile_ttbar->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_ttbar");
  TH2F *h_met_ll_inc_ttbar = (TH2F*) InputFile_ttbar->Get("h_met_ll_inc")->Clone("h_met_ll_inc_ttbar");
  TH2F *h_met_ll_tar0_ttbar = (TH2F*) InputFile_ttbar->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_ttbar");
  TH2F *h_met_ll_tar2_ttbar = (TH2F*) InputFile_ttbar->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_ttbar");
  TH2F *h_met_emu2_inc_ttbar = (TH2F*) InputFile_ttbar->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_ttbar");
  TH2F *h_met_emu2_tar0_ttbar = (TH2F*) InputFile_ttbar->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_ttbar");
  TH2F *h_met_emu2_tar2_ttbar = (TH2F*) InputFile_ttbar->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_ttbar");
  TH1F *h_phi_ttbar = (TH1F*) InputFile_ttbar->Get("h_phi")->Clone("h_phi_ttbar");
  TH1F *h_sumet_ee_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_ttbar");
  TH1F *h_sumet_mumu_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_ttbar");
  TH1F *h_sumet_emu_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_ttbar");

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
  TH1F *h_nvtx_scaled_ww = (TH1F*) InputFile_ww->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_ww");
  TH1F *h_nvtx_unscaled_ww = (TH1F*) InputFile_ww->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_ww");
  TH2F *h_met_ll_inc_ww = (TH2F*) InputFile_ww->Get("h_met_ll_inc")->Clone("h_met_ll_inc_ww");
  TH2F *h_met_ll_tar0_ww = (TH2F*) InputFile_ww->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_ww");
  TH2F *h_met_ll_tar2_ww = (TH2F*) InputFile_ww->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_ww");
  TH2F *h_met_emu2_inc_ww = (TH2F*) InputFile_ww->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_ww");
  TH2F *h_met_emu2_tar0_ww = (TH2F*) InputFile_ww->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_ww");
  TH2F *h_met_emu2_tar2_ww = (TH2F*) InputFile_ww->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_ww");
  TH1F *h_phi_ww = (TH1F*) InputFile_ww->Get("h_phi")->Clone("h_phi_ww");
  TH1F *h_sumet_ee_inc_ww = (TH1F*) InputFile_ww->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_ww");
  TH1F *h_sumet_mumu_inc_ww = (TH1F*) InputFile_ww->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_ww");
  TH1F *h_sumet_emu_inc_ww = (TH1F*) InputFile_ww->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_ww");

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
  TH1F *h_nvtx_scaled_ttv = (TH1F*) InputFile_ttv->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_ttv");
  TH1F *h_nvtx_unscaled_ttv = (TH1F*) InputFile_ttv->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_ttv");
  TH2F *h_met_ll_inc_ttv = (TH2F*) InputFile_ttv->Get("h_met_ll_inc")->Clone("h_met_ll_inc_ttv");
  TH2F *h_met_ll_tar0_ttv = (TH2F*) InputFile_ttv->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_ttv");
  TH2F *h_met_ll_tar2_ttv = (TH2F*) InputFile_ttv->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_ttv");
  TH2F *h_met_emu2_inc_ttv = (TH2F*) InputFile_ttv->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_ttv");
  TH2F *h_met_emu2_tar0_ttv = (TH2F*) InputFile_ttv->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_ttv");
  TH2F *h_met_emu2_tar2_ttv = (TH2F*) InputFile_ttv->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_ttv");
  TH1F *h_phi_ttv = (TH1F*) InputFile_ttv->Get("h_phi")->Clone("h_phi_ttv");
  TH1F *h_sumet_ee_inc_ttv = (TH1F*) InputFile_ttv->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_ttv");
  TH1F *h_sumet_mumu_inc_ttv = (TH1F*) InputFile_ttv->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_ttv");
  TH1F *h_sumet_emu_inc_ttv = (TH1F*) InputFile_ttv->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_ttv");

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
  TH1F *h_nvtx_scaled_vvv = (TH1F*) InputFile_vvv->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_vvv");
  TH1F *h_nvtx_unscaled_vvv = (TH1F*) InputFile_vvv->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_vvv");
  TH2F *h_met_ll_inc_vvv = (TH2F*) InputFile_vvv->Get("h_met_ll_inc")->Clone("h_met_ll_inc_vvv");
  TH2F *h_met_ll_tar0_vvv = (TH2F*) InputFile_vvv->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_vvv");
  TH2F *h_met_ll_tar2_vvv = (TH2F*) InputFile_vvv->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_vvv");
  TH2F *h_met_emu2_inc_vvv = (TH2F*) InputFile_vvv->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_vvv");
  TH2F *h_met_emu2_tar0_vvv = (TH2F*) InputFile_vvv->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_vvv");
  TH2F *h_met_emu2_tar2_vvv = (TH2F*) InputFile_vvv->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_vvv");
  TH1F *h_phi_vvv = (TH1F*) InputFile_vvv->Get("h_phi")->Clone("h_phi_vvv");
  TH1F *h_sumet_ee_inc_vvv = (TH1F*) InputFile_vvv->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_vvv");
  TH1F *h_sumet_mumu_inc_vvv = (TH1F*) InputFile_vvv->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_vvv");
  TH1F *h_sumet_emu_inc_vvv = (TH1F*) InputFile_vvv->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_vvv");

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
  TH1F *h_nvtx_scaled_data = (TH1F*) InputFile_data->Get("h_nvtx_scaled")->Clone("h_nvtx_scaled_data");
  TH1F *h_nvtx_unscaled_data = (TH1F*) InputFile_data->Get("h_nvtx_unscaled")->Clone("h_nvtx_unscaled_data");
  TH2F *h_met_ll_inc_data = (TH2F*) InputFile_data->Get("h_met_ll_inc")->Clone("h_met_ll_inc_data");
  TH2F *h_met_ll_tar0_data = (TH2F*) InputFile_data->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar0_data");
  TH2F *h_met_ll_tar2_data = (TH2F*) InputFile_data->Get("h_met_ll_tar_njets0")->Clone("h_met_ll_tar2_data");
  TH2F *h_met_emu2_inc_data = (TH2F*) InputFile_data->Get("h_met_emu2_inc")->Clone("h_met_emu2_inc_data");
  TH2F *h_met_emu2_tar0_data = (TH2F*) InputFile_data->Get("h_met_emu2_tar_njets0")->Clone("h_met_emu2_tar0_data");
  TH2F *h_met_emu2_tar2_data = (TH2F*) InputFile_data->Get("h_met_emu2_tar_njets2")->Clone("h_met_emu2_tar2_data");
  TH1F *h_phi_data = (TH1F*) InputFile_data->Get("h_phi")->Clone("h_phi_data");
  TH1F *h_sumet_ee_inc_data = (TH1F*) InputFile_data->Get("h_sumet_ee_inc")->Clone("h_sumet_ee_inc_data");
  TH1F *h_sumet_mumu_inc_data = (TH1F*) InputFile_data->Get("h_sumet_mumu_inc")->Clone("h_sumet_mumu_inc_data");
  TH1F *h_sumet_emu_inc_data = (TH1F*) InputFile_data->Get("h_sumet_emu_inc")->Clone("h_sumet_emu_inc_data");

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

  vector<TH2F*> v_met_ll_inc;
  v_met_ll_inc.push_back(h_met_ll_inc_zjets);
  v_met_ll_inc.push_back(h_met_ll_inc_wz);
  v_met_ll_inc.push_back(h_met_ll_inc_zz);
  v_met_ll_inc.push_back(h_met_ll_inc_ttbar);
  v_met_ll_inc.push_back(h_met_ll_inc_singlet);
  v_met_ll_inc.push_back(h_met_ll_inc_ww);
  v_met_ll_inc.push_back(h_met_ll_inc_ttv);
  v_met_ll_inc.push_back(h_met_ll_inc_vvv);
  v_met_ll_inc.push_back(h_met_ll_inc_data);

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

  vector<TH2F*> v_met_emu2_inc;
  v_met_emu2_inc.push_back(h_met_emu2_inc_zjets);
  v_met_emu2_inc.push_back(h_met_emu2_inc_wz);
  v_met_emu2_inc.push_back(h_met_emu2_inc_zz);
  v_met_emu2_inc.push_back(h_met_emu2_inc_ttbar);
  v_met_emu2_inc.push_back(h_met_emu2_inc_singlet);
  v_met_emu2_inc.push_back(h_met_emu2_inc_ww);
  v_met_emu2_inc.push_back(h_met_emu2_inc_ttv);
  v_met_emu2_inc.push_back(h_met_emu2_inc_vvv);
  v_met_emu2_inc.push_back(h_met_emu2_inc_data);

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

  vector<TH2F*> v_met_ll_tar0;
  v_met_ll_tar0.push_back(h_met_ll_tar0_zjets);
  v_met_ll_tar0.push_back(h_met_ll_tar0_wz);
  v_met_ll_tar0.push_back(h_met_ll_tar0_zz);
  v_met_ll_tar0.push_back(h_met_ll_tar0_ttbar);
  v_met_ll_tar0.push_back(h_met_ll_tar0_singlet);
  v_met_ll_tar0.push_back(h_met_ll_tar0_ww);
  v_met_ll_tar0.push_back(h_met_ll_tar0_ttv);
  v_met_ll_tar0.push_back(h_met_ll_tar0_vvv);
  v_met_ll_tar0.push_back(h_met_ll_tar0_data);

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

  vector<TH2F*> v_met_emu2_tar0;
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_zjets);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_wz);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_zz);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_ttbar);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_singlet);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_ww);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_ttv);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_vvv);
  v_met_emu2_tar0.push_back(h_met_emu2_tar0_data);

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

  vector<TH2F*> v_met_ll_tar2;
  v_met_ll_tar2.push_back(h_met_ll_tar2_zjets);
  v_met_ll_tar2.push_back(h_met_ll_tar2_wz);
  v_met_ll_tar2.push_back(h_met_ll_tar2_zz);
  v_met_ll_tar2.push_back(h_met_ll_tar2_ttbar);
  v_met_ll_tar2.push_back(h_met_ll_tar2_singlet);
  v_met_ll_tar2.push_back(h_met_ll_tar2_ww);
  v_met_ll_tar2.push_back(h_met_ll_tar2_ttv);
  v_met_ll_tar2.push_back(h_met_ll_tar2_vvv);
  v_met_ll_tar2.push_back(h_met_ll_tar2_data);

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

  vector<TH2F*> v_met_emu2_tar2;
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_zjets);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_wz);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_zz);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_ttbar);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_singlet);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_ww);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_ttv);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_vvv);
  v_met_emu2_tar2.push_back(h_met_emu2_tar2_data);

  vector<TH1F*> v_nvtx_scaled;
  v_nvtx_scaled.push_back(h_nvtx_scaled_zjets);
  v_nvtx_scaled.push_back(h_nvtx_scaled_wz);
  v_nvtx_scaled.push_back(h_nvtx_scaled_zz);
  v_nvtx_scaled.push_back(h_nvtx_scaled_ttbar);
  v_nvtx_scaled.push_back(h_nvtx_scaled_singlet);
  v_nvtx_scaled.push_back(h_nvtx_scaled_ww);
  v_nvtx_scaled.push_back(h_nvtx_scaled_ttv);
  v_nvtx_scaled.push_back(h_nvtx_scaled_vvv);
  v_nvtx_scaled.push_back(h_nvtx_scaled_data);

  vector<TH1F*> v_nvtx_unscaled;
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_zjets);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_wz);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_zz);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_ttbar);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_singlet);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_ww);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_ttv);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_vvv);
  v_nvtx_unscaled.push_back(h_nvtx_unscaled_data);
 
  vector<TH1F*> v_phi;
  v_phi.push_back(h_phi_zjets);
  v_phi.push_back(h_phi_wz);
  v_phi.push_back(h_phi_zz);
  v_phi.push_back(h_phi_ttbar);
  v_phi.push_back(h_phi_singlet);
  v_phi.push_back(h_phi_ww);
  v_phi.push_back(h_phi_ttv);
  v_phi.push_back(h_phi_vvv);
  v_phi.push_back(h_phi_data);

  vector<TH1F*> v_sumet_ee_inc;
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_zjets);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_wz);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_zz);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_ttbar);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_singlet);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_ww);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_ttv);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_vvv);
  v_sumet_ee_inc.push_back(h_sumet_ee_inc_data);

  vector<TH1F*> v_sumet_mumu_inc;
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_zjets);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_wz);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_zz);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_ttbar);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_singlet);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_ww);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_ttv);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_vvv);
  v_sumet_mumu_inc.push_back(h_sumet_mumu_inc_data);

  vector<TH1F*> v_sumet_emu_inc;
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_zjets);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_wz);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_zz);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_ttbar);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_singlet);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_ww);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_ttv);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_vvv);
  v_sumet_emu_inc.push_back(h_sumet_emu_inc_data);
  //------------------------------------FSBG-------------------------------------------------------------//

  //loop through ee and mumu and subract eu*R.  Before formatting, overflow, summing, and stacking
  int size=v_mll_ee_inc.size();
  
  double Ree = 0.440978;
  double Rmumu = 0.524251;
  double R = Ree+Rmumu;

  for(int i=0; i<size; i++)
	{
	  v_met_ee_inc[i]->Add(v_met_emu_inc[i],-Ree);      //don't forget minus
	  v_met_mumu_inc[i]->Add(v_met_emu_inc[i],-Rmumu);  //can use clones to not permanently change
	  v_met_ee_tar0[i]->Add(v_met_emu_tar0[i],-Ree);      
	  v_met_mumu_tar0[i]->Add(v_met_emu_tar0[i],-Rmumu);  
	  v_met_ee_tar2[i]->Add(v_met_emu_tar2[i],-Ree);      
	  v_met_mumu_tar2[i]->Add(v_met_emu_tar2[i],-Rmumu);
	  v_met_ll_inc[i]->Add(v_met_emu2_inc[i], -R);
	  v_met_ll_tar0[i]->Add(v_met_emu2_tar0[i], -R);
	  v_met_ll_tar2[i]->Add(v_met_emu2_tar2[i], -R);
	  v_sumet_ee_inc[i]->Add(v_sumet_emu_inc[i],-Ree);      //don't forget minus
	  v_sumet_mumu_inc[i]->Add(v_sumet_emu_inc[i],-Rmumu);  //can use clones to not permanently change
	}
  
  //----------------------add ee and mumu----------------------------//
  for(int i=0; i<size; i++)
	{
	  v_mll_ee_inc[i]->Add(v_mll_mumu_inc[i]);
	  v_mll_ee_tar[i]->Add(v_mll_mumu_tar[i]);
	  v_met_ee_inc[i]->Add(v_met_mumu_inc[i]);     
	  v_met_ee_tar0[i]->Add(v_met_mumu_tar0[i]);      
	  v_met_ee_tar2[i]->Add(v_met_mumu_tar2[i]);      
	  v_sumet_ee_inc[i]->Add(v_sumet_mumu_inc[i]);
	}
  /*
  vector<TH1F*> v_met_inc = v_met_ee_inc;    //maybe use these? better name, but makes it harder to separate ee and mumu
  vector<TH1F*> v_met_tar0 = v_met_ee_tar0;
  vector<TH1F*> v_met_tar2 = v_met_ee_tar2;
  */
  //---------------------------------------Format and Overflow---------------------------------------//

	
  for(int i=0; i < size; i++)
	{
	  // cout<<"i = "<<i<<endl;
	  // cout<<"length = "<<v_mll_mumu_inc.size()<<endl;
	  
	  h_format(v_mll_ee_inc[i]);       //no point in formatting emu's, or mumu's after I get rid of them...
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

	  h_format(v_nvtx_scaled[i]);
	  h_format(v_nvtx_unscaled[i]);	  
	  //overflow(v_nvtx_scaled[i]);   //figure this out later
	  // overflow(v_nvtx_unscaled[i])

	   h_format(v_sumet_ee_inc[i]);
	   // overflow(v_sumet_ee_inc[i]);  //do in scanchain
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
  TH2F *h_metsum_ll_inc = (TH2F*) v_met_ll_inc[0]->Clone("h_metsum_ll_inc");
  TH1F *h_metsum_emu_inc = (TH1F*) v_met_emu_inc[0]->Clone("h_metsum_emu_inc");
  TH1F *h_metsum_ee_tar0 = (TH1F*) v_met_ee_tar0[0]->Clone("h_metsum_ee_tar0");
  TH1F *h_metsum_mumu_tar0 = (TH1F*) v_met_mumu_tar0[0]->Clone("h_metsum_mumu_tar0");
  TH2F *h_metsum_ll_tar0 = (TH2F*) v_met_ll_tar0[0]->Clone("h_metsum_ll_tar0");
  TH1F *h_metsum_emu_tar0 = (TH1F*) v_met_emu_tar0[0]->Clone("h_metsum_emu_tar0");
  TH1F *h_metsum_ee_tar2 = (TH1F*) v_met_ee_tar2[0]->Clone("h_metsum_ee_tar2");
  TH1F *h_metsum_mumu_tar2 = (TH1F*) v_met_mumu_tar2[0]->Clone("h_metsum_mumu_tar2");
  TH2F *h_metsum_ll_tar2 = (TH2F*) v_met_ll_tar2[0]->Clone("h_metsum_ll_tar2");
  TH1F *h_metsum_emu_tar2 = (TH1F*) v_met_emu_tar2[0]->Clone("h_metsum_emu_tar2");
  TH1F *h_nvtxsum_scaled = (TH1F*) v_nvtx_scaled[0]->Clone("h_nvtxsum_scaled");
  TH1F *h_nvtxsum_unscaled = (TH1F*) v_nvtx_unscaled[0]->Clone("h_nvtxsum_unscaled");
  TH1F *h_phisum = (TH1F*) v_phi[0]->Clone("h_phisum");
  TH1F *h_sumetsum_ee_inc = (TH1F*) v_sumet_ee_inc[0]->Clone("h_sumetsum_ee_inc");

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
	  h_metsum_ll_inc->Add(v_met_ll_inc[i]);
	  h_metsum_emu_inc->Add(v_met_emu_inc[i]);
	  h_metsum_ee_tar0->Add(v_met_ee_tar0[i]);
	  h_metsum_mumu_tar0->Add(v_met_mumu_tar0[i]);
	  h_metsum_ll_tar0->Add(v_met_ll_tar0[i]);
	  h_metsum_emu_tar0->Add(v_met_emu_tar0[i]);
	  h_metsum_ee_tar2->Add(v_met_ee_tar2[i]);
	  h_metsum_mumu_tar2->Add(v_met_mumu_tar2[i]);
	  h_metsum_ll_tar2->Add(v_met_ll_tar2[i]);	
	  h_metsum_emu_tar2->Add(v_met_emu_tar2[i]);
	  h_nvtxsum_scaled->Add(v_nvtx_scaled[i]);
	  h_nvtxsum_unscaled->Add(v_nvtx_unscaled[i]);
	  h_phisum->Add(v_phi[i]);	  
	  h_sumetsum_ee_inc->Add(v_sumet_ee_inc[i]);
	}

  cout<<"\n End Summing \n"<<endl;

  //--------------------------------------------------------------------------------//

  //Get rid of stat box
  gStyle->SetOptStat(0);
  // gStyle->SetOptFit(1111);
  //-----------------------

  //-----------------------------Split up X and Y---------------------------------------------//
  TH1D *h_metsum_llx_inc = h_metsum_ll_inc->ProjectionX("h_metsum_llx_inc",0,-1,"e");
  TH1D *h_metsum_lly_inc = h_metsum_ll_inc->ProjectionY("h_metsum_lly_inc",0,-1,"e");
  TH1D *h_metsum_llx_tar0 = h_metsum_ll_tar0->ProjectionX("h_metsum_llx_tar0",0,-1,"e");
  TH1D *h_metsum_lly_tar0 = h_metsum_ll_tar0->ProjectionY("h_metsum_lly_tar0",0,-1,"e");
  TH1D *h_metsum_llx_tar2 = h_metsum_ll_tar2->ProjectionX("h_metsum_llx_tar2",0,-1,"e");
  TH1D *h_metsum_lly_tar2 = h_metsum_ll_tar2->ProjectionY("h_metsum_lly_tar2",0,-1,"e");

  TH1D *h_metdata_llx_inc = v_met_ll_inc[size-1]->ProjectionX("h_metdata_llx_inc",0,-1,"e");
  TH1D *h_metdata_lly_inc = v_met_ll_inc[size-1]->ProjectionY("h_metdata_lly_inc",0,-1,"e");
  TH1D *h_metdata_llx_tar0 = v_met_ll_tar0[size-1]->ProjectionX("h_metdata_llx_tar0",0,-1,"e");
  TH1D *h_metdata_lly_tar0 = v_met_ll_tar0[size-1]->ProjectionY("h_metdata_lly_tar0",0,-1,"e");
  TH1D *h_metdata_llx_tar2 = v_met_ll_tar2[size-1]->ProjectionX("h_metdata_llx_tar2",0,-1,"e");
  TH1D *h_metdata_lly_tar2 = v_met_ll_tar2[size-1]->ProjectionY("h_metdata_lly_tar2",0,-1,"e");
  //------------------------------------------------------------------------------------------//

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
  //don't bother stacking nvtx, or METx,METy yet, or phi
  THStack *hs_sumet_ee_inc = new THStack("hs_sumet_ee_inc","sumet ee inc");

  //Add to Stacked Histograms
  for(int i=size-2; i >= 0; i--)  //Reverse,don't add data
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
	  hs_sumet_ee_inc->Add(v_sumet_ee_inc[i]);
	}
    
  //-----------------------------------------------------------------------------------//
  //----------------------------------Drawing------------------------------------------//
  //-----------------------------------------------------------------------------------//  
  /*
  //-----------------------mll ee inc------------------------- -//
  //Canvas
  TCanvas *c1=new TCanvas("c1","mll ee inc ",800,800);
  TPad *pad_h = new TPad("pad_h","Histo Pad",0., 0, 1., 0.8);
  TPad *pad_r = new TPad("pad_r","Ratio Pad",0., 0.8, 1., 1.);

  pad_h->Draw();
  pad_r->Draw();

  pad_h->SetLogy();

  pad_h->cd();  //change pad

  v_mll_ee_inc[0]->SetTitle(" ee+mumu data vs MC Dilepton Mass (Inclusive)");
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
  */
  /*
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
  */
  /*  
  //-------------------Targeted ee-------------------------------//
  
  //Canvas
  TCanvas *c3=new TCanvas("c3","mll ee targeted",800,800);
  TPad *pad_tarh = new TPad("pad_tarh","Targeted Histo Pad",0., 0, 1., 0.8);
  TPad *pad_tarr = new TPad("pad_tarr","Targeted Ratio Pad",0., 0.8, 1., 1.);

  pad_tarh->Draw();
  pad_tarr->Draw();

  pad_tarh->SetLogy();

  pad_tarh->cd();  //change pad

  v_mll_ee_tar[0]->SetTitle("ee+mumu data vs MC Dilepton Mass (Targeted)");
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
  */
  /*  
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
  */
  
  //---------------------MET ee INC------------------------------//
  //Canvas
  TCanvas *c5=new TCanvas("c5","Met ee inc",800,800);
  TPad *pad_h5 = new TPad("pad_h5","Histo Pad5",0., 0, 1., 0.8);
  TPad *pad_r5 = new TPad("pad_r5","Ratio Pad5",0., 0.8, 1., 1.);

  pad_h5->Draw();
  pad_r5->Draw();

  pad_h5->SetLogy();  //c1->pad_h

  pad_h5->cd();  //change pad

  v_met_ee_inc[0]->SetTitle("ee+mumu data vs MC MET (Inclusive)");
  v_met_ee_inc[0]->Draw();
  hs_met_ee_inc->Draw("histsame");
 
  v_met_ee_inc[size-1]->Draw("same e1 x0");

  pad_h5->RedrawAxis();  //c1->pad_h
  
  double width_mc = find_width(h_metsum_ee_inc);
  double width_data = find_width(v_met_ee_inc[size-1]);
  cout<<"\nINC MC Width = "<<width_mc<<endl;
  cout<<"\nINC data Width = "<<width_data<<endl;
  cout<<"\nINC Width difference = "<<width_data-width_mc<<endl;
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
   
  /*  
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
  */
  
  //---------------------MET ee Tar0-----------------------------//
  //Canvas
  TCanvas *c7=new TCanvas("c7","Met ee tar0",800,800);
  TPad *pad_h7 = new TPad("pad_h7","Histo Pad7",0., 0, 1., 0.8);
  TPad *pad_r7 = new TPad("pad_r7","Ratio Pad7",0., 0.8, 1., 1.);

  pad_h7->Draw();
  pad_r7->Draw();

  pad_h7->SetLogy();  //c1->pad_h

  pad_h7->cd();  //change pad

  v_met_ee_tar0[0]->SetTitle("ee+mumu data vs MC MET (no jets) ");
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
  
  /*
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
  */
    
  //---------------------MET ee Tar2------------------------------//
  //Canvas
  TCanvas *c9=new TCanvas("c9","Met ee tar2",800,800);
  TPad *pad_h9 = new TPad("pad_h9","Histo Pad9",0., 0, 1., 0.8);
  TPad *pad_r9 = new TPad("pad_r9","Ratio Pad9",0., 0.8, 1., 1.);

  pad_h9->Draw();
  pad_r9->Draw();

  pad_h9->SetLogy();  //c1->pad_h

  pad_h9->cd();  //change pad

  v_met_ee_tar2[0]->SetTitle("ee+mumu data vs MC MET (njets >= 2) ");
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
  
  /*
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
  */
  /* 
  //----------------------DONT MAKE SENSE ANYMORE-------------------------------//
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
  */
  /*
 //---------------------nvtx scaled------------------------------//
  //Canvas
  TCanvas *c17=new TCanvas("c17","Number of Vertices (Scaled)",800,800);
  TPad *pad_h17 = new TPad("pad_h17","Histo Pad17",0., 0, 1., 0.8);
  TPad *pad_r17 = new TPad("pad_r17","Ratio Pad17",0., 0.8, 1., 1.);

  pad_h17->Draw();
  pad_r17->Draw();

  // pad_h17->SetLogy();  //c1->pad_h

  pad_h17->cd();  //change pad
  
  h_nvtxsum_scaled->SetLineColor(kRed);
  h_nvtxsum_scaled->SetMarkerColor(kRed);
  h_nvtxsum_scaled->SetMarkerStyle(8);

  h_nvtxsum_scaled->Draw();
  v_nvtx_scaled[size-1]->Draw("same e1");

  pad_h17->RedrawAxis();

  //-------------------------------Legend---------------------------------------//
  TLegend *leg17 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg17->SetLineColor(kWhite);
  leg17->SetTextFont(42); 
  leg17->SetTextSize(0.026);
  leg17->SetShadowColor(kWhite); 
  leg17->SetFillColor(kWhite); 
  leg17->AddEntry("h_nvtxsum_scaled","MC","ep");
  leg17->AddEntry("h_nvtx_scaled_data","data","ep");

  leg17->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r17->cd();

  TH1F *h_nvtx_scaled_data_clone = (TH1F*) v_nvtx_scaled[size-1]->Clone("h_nvtx_scaled_data_clone");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_nvtx_scaled_data_clone->Divide(h_nvtx_scaled_data_clone,h_nvtxsum_scaled);

  h_nvtx_scaled_data_clone->GetXaxis()->SetRangeUser(0,35); 
  h_nvtx_scaled_data_clone->GetYaxis()->SetRangeUser(0,2);
  h_nvtx_scaled_data_clone->GetYaxis()->SetNdivisions(4);
  h_nvtx_scaled_data_clone->GetYaxis()->SetLabelSize(.12);
  h_nvtx_scaled_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r17->SetGridy();

  h_nvtx_scaled_data_clone->Draw();

  //-------------------------------------------------------------//

 //---------------------nvtx unscaled------------------------------//
  //Canvas
  TCanvas *c18=new TCanvas("c18","Number of Vertices (Unscaled)",800,800);
  TPad *pad_h18 = new TPad("pad_h18","Histo Pad18",0., 0, 1., 0.8);
  TPad *pad_r18 = new TPad("pad_r18","Ratio Pad18",0., 0.8, 1., 1.);

  pad_h18->Draw();
  pad_r18->Draw();

  //pad_h18->SetLogy();  //c1->pad_h

  pad_h18->cd();  //change pad
  
  h_nvtxsum_unscaled->SetLineColor(kRed);
  h_nvtxsum_unscaled->SetMarkerColor(kRed);
  h_nvtxsum_unscaled->SetMarkerStyle(8);

  h_nvtxsum_unscaled->Draw();
  v_nvtx_unscaled[size-1]->Draw("same e1");

  pad_h18->RedrawAxis();

  //-------------------------------Legend---------------------------------------//
  TLegend *leg18 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg18->SetLineColor(kWhite);
  leg18->SetTextFont(42); 
  leg18->SetTextSize(0.026);
  leg18->SetShadowColor(kWhite); 
  leg18->SetFillColor(kWhite); 
  leg18->AddEntry("h_nvtxsum_unscaled","MC","ep");
  leg18->AddEntry("h_nvtx_unscaled_data","data","ep");

  leg18->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r18->cd();

  TH1F *h_nvtx_unscaled_data_clone = (TH1F*) v_nvtx_unscaled[size-1]->Clone("h_nvtx_unscaled_data_clone");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_nvtx_unscaled_data_clone->Divide(h_nvtx_unscaled_data_clone,h_nvtxsum_unscaled);

  h_nvtx_unscaled_data_clone->GetXaxis()->SetRangeUser(0,35); 
  h_nvtx_unscaled_data_clone->GetYaxis()->SetRangeUser(0,2);
  h_nvtx_unscaled_data_clone->GetYaxis()->SetNdivisions(4);
  h_nvtx_unscaled_data_clone->GetYaxis()->SetLabelSize(.12);
  h_nvtx_unscaled_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r18->SetGridy();

  h_nvtx_unscaled_data_clone->Draw();

  //-------------------------------------------------------------//
  */
  /*  
  //-----------------------X vs Y MC MET INC-----------------------------//
  //Canvas
  TCanvas *c19=new TCanvas("c19","X vs Y MC Inclusive MET",800,800);
  TPad *pad_h19 = new TPad("pad_h19","Histo Pad19",0., 0, 1., 0.8);
  TPad *pad_r19 = new TPad("pad_r19","Ratio Pad19",0., 0.8, 1., 1.);

  pad_h19->Draw();
  pad_r19->Draw();

  pad_h19->SetLogy();

  pad_h19->cd();
  
  h_metsum_llx_inc->GetXaxis()->SetRangeUser(-200,200);
  h_metsum_lly_inc->GetXaxis()->SetRangeUser(-200,200);
  
  h_metsum_llx_inc->SetMarkerColor(kRed);
  h_metsum_llx_inc->SetMarkerStyle(7);
  h_metsum_llx_inc->SetLineColor(kRed);
  h_metsum_lly_inc->SetMarkerColor(kBlue);
  h_metsum_lly_inc->SetMarkerStyle(7); 
  h_metsum_lly_inc->SetLineColor(kBlue);
  h_metsum_llx_inc->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metsum_llx_inc->GetXaxis()->SetTitle("MET_{ll}");
  h_metsum_llx_inc->SetTitle("MC MET X vs Y (Inclusive)");

  h_metsum_llx_inc->Draw();
  h_metsum_lly_inc->Draw("histsame");//("same e1");
 
  pad_h19->RedrawAxis();
  //--------------------------------
  //all_fitting(h_metsum_llx_inc);
  // all_fitting(h_metsum_lly_inc);
  //--------------------------------
 
  // TFile *metxy_inc_file = new TFile("metxy_inc.root","recreate");
  //metxy_inc_file->cd();
  //h_metsum_llx_inc->Write();
  //h_metsum_lly_inc->Write();

  //-------------------------------Legend---------------------------------------//
  TLegend *leg19 = new TLegend(0.79,0.53,0.87,0.63);//(0.78, 0.63, 0.87, 0.89);
  leg19->SetLineColor(kWhite);
  leg19->SetTextFont(42); 
  leg19->SetTextSize(0.026);
  leg19->SetShadowColor(kWhite); 
  leg19->SetFillColor(kWhite); 
  leg19->AddEntry("h_metsum_llx_inc","MC METx","ep");
  leg19->AddEntry("h_metsum_lly_inc","MC METy","ep");
  // leg19->AddEntry("sum","MC MET Wieghted Average","ep");

  leg19->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r19->cd();

  TH1D *h_metsum_llx_inc_clone = (TH1D*) h_metsum_llx_inc->Clone("h_metsum_llx_inc_clone");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_metsum_llx_inc_clone->Divide(h_metsum_llx_inc_clone,h_metsum_lly_inc);

  h_metsum_llx_inc_clone->GetXaxis()->SetRangeUser(-200,200); 
  h_metsum_llx_inc_clone->GetYaxis()->SetRangeUser(0,2);
  h_metsum_llx_inc_clone->GetYaxis()->SetNdivisions(4);
  h_metsum_llx_inc_clone->GetYaxis()->SetLabelSize(.12);
  h_metsum_llx_inc_clone->GetYaxis()->SetTitle("X/Y");
  h_metsum_llx_inc_clone->SetLineColor(kBlack);
  h_metsum_llx_inc_clone->SetMarkerColor(kBlack);

  pad_r19->SetGridy();

  h_metsum_llx_inc_clone->Draw();

  //-------------------------------------------------------------//
  
  //-----------------------X vs Y MC MET TAR0-----------------------------//
  //Canvas
  TCanvas *c20=new TCanvas("c20","X vs Y MC No Jets MET",800,800);
  TPad *pad_h20 = new TPad("pad_h20","Histo Pad20",0., 0, 1., 0.8);
  TPad *pad_r20 = new TPad("pad_r20","Ratio Pad20",0., 0.8, 1., 1.);

  pad_h20->Draw();
  pad_r20->Draw();

  pad_h20->SetLogy();

  pad_h20->cd();
  
  h_metsum_llx_tar0->GetXaxis()->SetRangeUser(-200,200);
  h_metsum_lly_tar0->GetXaxis()->SetRangeUser(-200,200);
  
  h_metsum_llx_tar0->SetMarkerColor(kRed);
  h_metsum_llx_tar0->SetLineColor(kRed);
  h_metsum_lly_tar0->SetMarkerColor(kBlue);
  h_metsum_lly_tar0->SetMarkerColor(kBlue);
  h_metsum_llx_tar0->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metsum_llx_tar0->GetXaxis()->SetTitle("MET_{ll}");
  h_metsum_llx_tar0->SetTitle("MC MET X vs Y (No Jets)");

  h_metsum_llx_tar0->Draw();
  h_metsum_lly_tar0->Draw("histsame");//("same e1");

  pad_h20->RedrawAxis();

  //-------------------------------Legend---------------------------------------//
  TLegend *leg20 = new TLegend(0.79, 0.53, 0.87, 0.63);
  leg20->SetLineColor(kWhite);
  leg20->SetTextFont(42); 
  leg20->SetTextSize(0.026);
  leg20->SetShadowColor(kWhite); 
  leg20->SetFillColor(kWhite); 
  leg20->AddEntry("h_metsum_llx_tar0","MC METx","ep");
  leg20->AddEntry("h_metsum_lly_tar0","MC METy","ep");

  leg20->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r20->cd();

  TH1D *h_metsum_llx_tar0_clone = (TH1D*) h_metsum_llx_tar0->Clone("h_metsum_llx_tar0_clone");

  //  overflow(h_metsum_ee_tar0);
  
  //Divide
  h_metsum_llx_tar0_clone->Divide(h_metsum_llx_tar0_clone,h_metsum_lly_tar0);

  h_metsum_llx_tar0_clone->GetXaxis()->SetRangeUser(-200,200); 
  h_metsum_llx_tar0_clone->GetYaxis()->SetRangeUser(0,2);
  h_metsum_llx_tar0_clone->GetYaxis()->SetNdivisions(4);
  h_metsum_llx_tar0_clone->GetYaxis()->SetLabelSize(.12);
  h_metsum_llx_tar0_clone->GetYaxis()->SetTitle("X/Y");
  h_metsum_llx_tar0_clone->SetLineColor(kBlack);
  h_metsum_llx_tar0_clone->SetMarkerColor(kBlack);

  pad_r20->SetGridy();

  h_metsum_llx_tar0_clone->Draw();

  //-------------------------------------------------------------//

  //-----------------------X vs Y MC MET TAR2-----------------------------//
  //Canvas
  TCanvas *c21=new TCanvas("c21","X vs Y MC 2 Jets MET",800,800);
  TPad *pad_h21 = new TPad("pad_h21","Histo Pad21",0., 0, 1., 0.8);
  TPad *pad_r21 = new TPad("pad_r21","Ratio Pad21",0., 0.8, 1., 1.);

  pad_h21->Draw();
  pad_r21->Draw();

  pad_h21->SetLogy();

  pad_h21->cd();
  
  h_metsum_llx_tar2->GetXaxis()->SetRangeUser(-200,200);
  h_metsum_lly_tar2->GetXaxis()->SetRangeUser(-200,200);
  
  h_metsum_llx_tar2->SetMarkerColor(kRed);
  h_metsum_llx_tar2->SetLineColor(kRed);
  h_metsum_lly_tar2->SetMarkerColor(kBlue);
  h_metsum_lly_tar2->SetMarkerColor(kBlue);
  h_metsum_llx_tar2->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metsum_llx_tar2->GetXaxis()->SetTitle("MET_{ll}");
  h_metsum_llx_tar2->SetTitle("MC MET X vs Y (2 Jets)");

  h_metsum_llx_tar2->Draw();
  h_metsum_lly_tar2->Draw("histsame");//("same e1");

 
  pad_h21->RedrawAxis();

  //-------------------------------Legend---------------------------------------//
  TLegend *leg21 = new TLegend(0.79, 0.53, 0.87, 0.63);
  leg21->SetLineColor(kWhite);
  leg21->SetTextFont(42); 
  leg21->SetTextSize(0.026);
  leg21->SetShadowColor(kWhite); 
  leg21->SetFillColor(kWhite); 
  leg21->AddEntry("h_metsum_llx_tar2","MC METx","ep");
  leg21->AddEntry("h_metsum_lly_tar2","MC METy","ep");

  leg21->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r21->cd();

  TH1D *h_metsum_llx_tar2_clone = (TH1D*) h_metsum_llx_tar2->Clone("h_metsum_llx_tar2_clone");

  //  overflow(h_metsum_ee_tar2);
  
  //Divide
  h_metsum_llx_tar2_clone->Divide(h_metsum_llx_tar2_clone,h_metsum_lly_tar2);

  h_metsum_llx_tar2_clone->GetXaxis()->SetRangeUser(-200,200); 
  h_metsum_llx_tar2_clone->GetYaxis()->SetRangeUser(0,2);
  h_metsum_llx_tar2_clone->GetYaxis()->SetNdivisions(4);
  h_metsum_llx_tar2_clone->GetYaxis()->SetLabelSize(.12);
  h_metsum_llx_tar2_clone->GetYaxis()->SetTitle("X/Y");
  h_metsum_llx_tar2_clone->SetLineColor(kBlack);
  h_metsum_llx_tar2_clone->SetMarkerColor(kBlack);

  pad_r21->SetGridy();

  h_metsum_llx_tar2_clone->Draw();

  //-------------------------------------------------------------//
  //-----------------------X vs Y DATA MET INC-----------------------------//
  //Canvas
  TCanvas *c22=new TCanvas("c22","X vs Y Data Inclusive MET",800,800);
  TPad *pad_h22 = new TPad("pad_h22","Histo Pad22",0., 0, 1., 0.8);
  TPad *pad_r22 = new TPad("pad_r22","Ratio Pad22",0., 0.8, 1., 1.);

  pad_h22->Draw();
  pad_r22->Draw();

  pad_h22->SetLogy();

  pad_h22->cd();

  h_metdata_llx_inc->GetXaxis()->SetRangeUser(-200,200);
  h_metdata_lly_inc->GetXaxis()->SetRangeUser(-200,200);
  
  h_metdata_llx_inc->SetMarkerColor(kRed);
  h_metdata_llx_inc->SetMarkerStyle(7);
  h_metdata_llx_inc->SetLineColor(kRed);
  h_metdata_lly_inc->SetMarkerColor(kBlue);
  h_metdata_lly_inc->SetMarkerStyle(7); 
  h_metdata_lly_inc->SetLineColor(kBlue);
  h_metdata_llx_inc->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metdata_llx_inc->GetXaxis()->SetTitle("MET_{ll}");
  h_metdata_llx_inc->SetTitle("Data MET X vs Y (Inclusive)");

  h_metdata_llx_inc->Draw();
  h_metdata_lly_inc->Draw("histsame");//("same e1");
 
  pad_h22->RedrawAxis();

  // all_fitting(h_metdata_llx_inc);
  // all_fitting(h_metsum_lly_inc);
  
  // h_metdata_llx_inc->Write();
  //h_metdata_lly_inc->Write();
  //metxy_inc_file->Close();
  //-------------------------------Legend---------------------------------------//
  TLegend *leg22 = new TLegend(0.79, 0.53, 0.87, 0.63);
  leg22->SetLineColor(kWhite);
  leg22->SetTextFont(42); 
  leg22->SetTextSize(0.026);
  leg22->SetShadowColor(kWhite); 
  leg22->SetFillColor(kWhite); 
  leg22->AddEntry("h_metdata_llx_inc","Data METx","ep");
  leg22->AddEntry("h_metdata_lly_inc","Data METy","ep");

  leg22->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r22->cd();

  TH1D *h_metdata_llx_inc_clone = (TH1D*) h_metdata_llx_inc->Clone("h_metdata_llx_inc_clone");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_metdata_llx_inc_clone->Divide(h_metdata_llx_inc_clone,h_metdata_lly_inc);

  h_metdata_llx_inc_clone->GetXaxis()->SetRangeUser(-200,200); 
  h_metdata_llx_inc_clone->GetYaxis()->SetRangeUser(0,2);
  h_metdata_llx_inc_clone->GetYaxis()->SetNdivisions(4);
  h_metdata_llx_inc_clone->GetYaxis()->SetLabelSize(.12);
  h_metdata_llx_inc_clone->GetYaxis()->SetTitle("X/Y");
  h_metdata_llx_inc_clone->SetLineColor(kBlack);
  h_metdata_llx_inc_clone->SetMarkerColor(kBlack);

  pad_r22->SetGridy();

  h_metdata_llx_inc_clone->Draw();

  //-------------------------------------------------------------//
 //-----------------------X vs Y DATA MET TAR0-----------------------------//
  //Canvas
  TCanvas *c23=new TCanvas("c23","X vs Y Data No Jets MET",800,800);
  TPad *pad_h23 = new TPad("pad_h23","Histo Pad23",0., 0, 1., 0.8);
  TPad *pad_r23 = new TPad("pad_r23","Ratio Pad23",0., 0.8, 1., 1.);

  pad_h23->Draw();
  pad_r23->Draw();

  pad_h23->SetLogy();

  pad_h23->cd();

  h_metdata_llx_tar0->GetXaxis()->SetRangeUser(-200,200);
  h_metdata_lly_tar0->GetXaxis()->SetRangeUser(-200,200);
  
  h_metdata_llx_tar0->SetMarkerColor(kRed);
  h_metdata_llx_tar0->SetMarkerStyle(7);
  h_metdata_llx_tar0->SetLineColor(kRed);
  h_metdata_lly_tar0->SetMarkerColor(kBlue);
  h_metdata_lly_tar0->SetMarkerStyle(7); 
  h_metdata_lly_tar0->SetLineColor(kBlue);
  h_metdata_llx_tar0->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metdata_llx_tar0->GetXaxis()->SetTitle("MET_{ll}");
  h_metdata_llx_tar0->SetTitle("Data MET X vs Y (No Jets)");

  h_metdata_llx_tar0->Draw();
  h_metdata_lly_tar0->Draw("histsame");//("same e1");
 
  pad_h23->RedrawAxis();
  
  //-------------------------------Legend---------------------------------------//
  TLegend *leg23 = new TLegend(0.79, 0.53, 0.87, 0.63);
  leg23->SetLineColor(kWhite);
  leg23->SetTextFont(42); 
  leg23->SetTextSize(0.026);
  leg23->SetShadowColor(kWhite); 
  leg23->SetFillColor(kWhite); 
  leg23->AddEntry("h_metdata_llx_tar0","Data METx","ep");
  leg23->AddEntry("h_metdata_lly_tar0","Data METy","ep");

  leg23->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r23->cd();

  TH1D *h_metdata_llx_tar0_clone = (TH1D*) h_metdata_llx_tar0->Clone("h_metdata_llx_tar0_clone");

  //  overflow(h_metsum_ee_tar0);
  
  //Divide
  h_metdata_llx_tar0_clone->Divide(h_metdata_llx_tar0_clone,h_metdata_lly_tar0);

  h_metdata_llx_tar0_clone->GetXaxis()->SetRangeUser(-200,200); 
  h_metdata_llx_tar0_clone->GetYaxis()->SetRangeUser(0,2);
  h_metdata_llx_tar0_clone->GetYaxis()->SetNdivisions(4);
  h_metdata_llx_tar0_clone->GetYaxis()->SetLabelSize(.12);
  h_metdata_llx_tar0_clone->GetYaxis()->SetTitle("X/Y");
  h_metdata_llx_tar0_clone->SetLineColor(kBlack);
  h_metdata_llx_tar0_clone->SetMarkerColor(kBlack);

  pad_r23->SetGridy();

  h_metdata_llx_tar0_clone->Draw();

  //-------------------------------------------------------------//

 //-----------------------X vs Y DATA MET TAR2-----------------------------//
  //Canvas
  TCanvas *c24=new TCanvas("c24","X vs Y Data 2 Jets MET",800,800);
  TPad *pad_h24 = new TPad("pad_h24","Histo Pad24",0., 0, 1., 0.8);
  TPad *pad_r24 = new TPad("pad_r24","Ratio Pad24",0., 0.8, 1., 1.);

  pad_h24->Draw();
  pad_r24->Draw();

  pad_h24->SetLogy();

  pad_h24->cd();

  h_metdata_llx_tar2->GetXaxis()->SetRangeUser(-200,200);
  h_metdata_lly_tar2->GetXaxis()->SetRangeUser(-200,200);
  
  h_metdata_llx_tar2->SetMarkerColor(kRed);
  h_metdata_llx_tar2->SetMarkerStyle(7);
  h_metdata_llx_tar2->SetLineColor(kRed);
  h_metdata_lly_tar2->SetMarkerColor(kBlue);
  h_metdata_lly_tar2->SetMarkerStyle(7); 
  h_metdata_lly_tar2->SetLineColor(kBlue);
  h_metdata_llx_tar2->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metdata_llx_tar2->GetXaxis()->SetTitle("MET_{ll}");
  h_metdata_llx_tar2->SetTitle("Data MET X vs Y (2 Jets)");

  h_metdata_llx_tar2->Draw();
  h_metdata_lly_tar2->Draw("histsame");//("same e1");
 
  pad_h24->RedrawAxis();
  
  //-------------------------------Legend---------------------------------------//
  TLegend *leg24 = new TLegend(0.79, 0.53, 0.87, 0.63);
  leg24->SetLineColor(kWhite);
  leg24->SetTextFont(42); 
  leg24->SetTextSize(0.026);
  leg24->SetShadowColor(kWhite); 
  leg24->SetFillColor(kWhite); 
  leg24->AddEntry("h_metdata_llx_tar2","Data METx","ep");
  leg24->AddEntry("h_metdata_lly_tar2","Data METy","ep");

  leg24->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r24->cd();

  TH1D *h_metdata_llx_tar2_clone = (TH1D*) h_metdata_llx_tar2->Clone("h_metdata_llx_tar2_clone");

  //  overflow(h_metsum_ee_tar2);
  
  //Divide
  h_metdata_llx_tar2_clone->Divide(h_metdata_llx_tar2_clone,h_metdata_lly_tar2);

  h_metdata_llx_tar2_clone->GetXaxis()->SetRangeUser(-200,200); 
  h_metdata_llx_tar2_clone->GetYaxis()->SetRangeUser(0,2);
  h_metdata_llx_tar2_clone->GetYaxis()->SetNdivisions(4);
  h_metdata_llx_tar2_clone->GetYaxis()->SetLabelSize(.12);
  h_metdata_llx_tar2_clone->GetYaxis()->SetTitle("X/Y");
  h_metdata_llx_tar2_clone->SetLineColor(kBlack);
  h_metdata_llx_tar2_clone->SetMarkerColor(kBlack);

  pad_r24->SetGridy();

  h_metdata_llx_tar2_clone->Draw();

  //-------------------------------------------------------------//
  */
  //-------------------------------------------------------------//
  //--------------------------PHI--------------------------------//
  //Canvas
  TCanvas *c25=new TCanvas("c25","MET phi MC vs DataT",800,800);
  TPad *pad_h25 = new TPad("pad_h25","Histo Pad25",0., 0, 1., 0.8);
  TPad *pad_r25 = new TPad("pad_r25","Ratio Pad25",0., 0.8, 1., 1.);

  pad_h25->Draw();
  pad_r25->Draw();

  pad_h25->cd();

  h_phisum->GetXaxis()->SetRangeUser(-4,4);
  h_phisum->GetYaxis()->SetRangeUser(0,6.5e5);
  
  h_phisum->SetMarkerColor(kRed);
  h_phisum->SetMarkerStyle(7);
  h_phisum->SetLineColor(kRed);
  v_phi[size-1]->SetMarkerColor(kBlue);
  v_phi[size-1]->SetMarkerStyle(7); 
  v_phi[size-1]->SetLineColor(kBlue);
  h_phisum->GetYaxis()->SetTitle("Events / 1 GeV");
  h_phisum->GetXaxis()->SetTitle("MET_{ll}");
  h_phisum->SetTitle("MET phi MC vs Data");

  h_phisum->Draw();
  v_phi[size-1]->Draw("histsame");//("same e1");
 
  pad_h25->RedrawAxis();
  

  //-------------------------------Legend---------------------------------------//
  TLegend *leg25 = new TLegend(0.79, 0.53, 0.87, 0.63);
  leg25->SetLineColor(kWhite);
  leg25->SetTextFont(42); 
  leg25->SetTextSize(0.026);
  leg25->SetShadowColor(kWhite); 
  leg25->SetFillColor(kWhite); 
  leg25->AddEntry("h_phisum","MC MET phi","ep");
  //  leg25->AddEntry("v_phi[size-1]","Data MET phi","ep");
  leg25->AddEntry("h_phi_data","Data MET phi","ep");

  leg25->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//
  
  pad_r25->cd();

  TH1D *h_phidata_clone = (TH1D*) v_phi[size-1]->Clone("h_phidata_clone");

  //  overflow(h_metsum_ee_tar2);
  
  //Divide
  h_phidata_clone->Divide(h_phidata_clone,h_phisum);

  h_phidata_clone->GetXaxis()->SetRangeUser(-4,4); 
  h_phidata_clone->GetYaxis()->SetRangeUser(0,2);
  h_phidata_clone->GetYaxis()->SetNdivisions(4);
  h_phidata_clone->GetYaxis()->SetLabelSize(.12);
  h_phidata_clone->GetYaxis()->SetTitle("Data/MC");
  h_phidata_clone->SetLineColor(kBlack);
  h_phidata_clone->SetMarkerColor(kBlack);

  pad_r25->SetGridy();

  h_phidata_clone->Draw();
  /*
  TH1D *h_phidata_clone_scaled = (TH1D*) h_phidata_clone->Clone("h_phidata_clone_scaled");

  double k = h_phisum->Integral(0,-1) / v_phi[size-1]->Integral(0,-1);
  h_phidata_clone_scaled->Scale( k );  
  
  TFile *phiRatioFile = new TFile("phiRatio_t1new.root","recreate");
  phiRatioFile->cd();
  h_phidata_clone->Write();
  h_phidata_clone_scaled->Write();
  phiRatioFile->Close();  
  */
  //---------------------------------------//
  /*
  //-----------------------MC vs Data METx INC-----------------------------//
  //Canvas
  TCanvas *c26=new TCanvas("c26","Data vs MC Inclusive METx",800,800);
  TPad *pad_h26 = new TPad("pad_h26","Histo Pad26",0., 0, 1., 0.8);
  TPad *pad_r26 = new TPad("pad_r26","Ratio Pad26",0., 0.8, 1., 1.);

  pad_h26->Draw();
  pad_r26->Draw();

  pad_h26->SetLogy();

  pad_h26->cd();
  
  h_metsum_llx_inc->GetXaxis()->SetRangeUser(-200,200);
  h_metdata_llx_inc->GetXaxis()->SetRangeUser(-200,200);
  
  h_metsum_llx_inc->SetMarkerColor(kRed);
  h_metsum_llx_inc->SetMarkerStyle(7);
  h_metsum_llx_inc->SetLineColor(kRed);
  h_metdata_llx_inc->SetMarkerColor(kBlue);
  h_metdata_llx_inc->SetMarkerStyle(7); 
  h_metdata_llx_inc->SetLineColor(kBlue);
  h_metsum_llx_inc->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metsum_llx_inc->GetXaxis()->SetTitle("MET_{ll}");
  h_metsum_llx_inc->SetTitle("Data vs MC METx (Inclusive)");

  h_metsum_llx_inc->Draw();
  h_metdata_llx_inc->Draw("histsame");//("same e1");
 
  pad_h26->RedrawAxis();
  //--------------------------------
  //all_fitting(h_metsum_llx_inc);
  // all_fitting(h_metdata_llx_inc);
  //--------------------------------
  
 // TFile *metxy_inc_file = new TFile("metxy_inc.root","recreate");
 // metxy_inc_file->cd();
 // h_metsum_llx_inc->Write();
 // h_metdata_llx_inc->Write();
  
  //-------------------------------Legend---------------------------------------//
  TLegend *leg26 = new TLegend(0.79,0.53,0.87,0.63);//(0.78, 0.63, 0.87, 0.89);
  leg26->SetLineColor(kWhite);
  leg26->SetTextFont(42); 
  leg26->SetTextSize(0.026);
  leg26->SetShadowColor(kWhite); 
  leg26->SetFillColor(kWhite); 
  leg26->AddEntry("h_metsum_llx_inc","MC METx","ep");
  leg26->AddEntry("h_metdata_llx_inc","Data METx","ep");
  // leg26->AddEntry("sum","MC MET Wieghted Average","ep");

  leg26->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r26->cd();

  TH1D *h_metdata_llx_inc_clone2 = (TH1D*) h_metdata_llx_inc->Clone("h_metdata_llx_inc_clone2");

  //  overflow(h_metsum_ee_inc);
  
  //Divide
  h_metdata_llx_inc_clone2->Divide(h_metdata_llx_inc_clone,h_metsum_llx_inc);

  h_metdata_llx_inc_clone2->GetXaxis()->SetRangeUser(-200,200); 
  h_metdata_llx_inc_clone2->GetYaxis()->SetRangeUser(0,2);
  h_metdata_llx_inc_clone2->GetYaxis()->SetNdivisions(4);
  h_metdata_llx_inc_clone2->GetYaxis()->SetLabelSize(.12);
  h_metdata_llx_inc_clone2->GetYaxis()->SetTitle("X/Y");
  h_metdata_llx_inc_clone2->SetLineColor(kBlack);
  h_metdata_llx_inc_clone2->SetMarkerColor(kBlack);

  pad_r26->SetGridy();

  h_metdata_llx_inc_clone2->Draw();

  //-------------------------------------------------------------//

 //-----------------------MC vs Data METx TAR0-----------------------------//
  //Canvas
  TCanvas *c27=new TCanvas("c27","Data vs MC No Jets METx",800,800);
  TPad *pad_h27 = new TPad("pad_h27","Histo Pad27",0., 0, 1., 0.8);
  TPad *pad_r27 = new TPad("pad_r27","Ratio Pad27",0., 0.8, 1., 1.);

  pad_h27->Draw();
  pad_r27->Draw();

  pad_h27->SetLogy();

  pad_h27->cd();
  
  h_metsum_llx_tar0->GetXaxis()->SetRangeUser(-200,200);
  h_metdata_llx_tar0->GetXaxis()->SetRangeUser(-200,200);
  
  h_metsum_llx_tar0->SetMarkerColor(kRed);
  h_metsum_llx_tar0->SetMarkerStyle(7);
  h_metsum_llx_tar0->SetLineColor(kRed);
  h_metdata_llx_tar0->SetMarkerColor(kBlue);
  h_metdata_llx_tar0->SetMarkerStyle(7); 
  h_metdata_llx_tar0->SetLineColor(kBlue);
  h_metsum_llx_tar0->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metsum_llx_tar0->GetXaxis()->SetTitle("MET_{ll}");
  h_metsum_llx_tar0->SetTitle("Data vs MC METx (No Jets)");

  h_metsum_llx_tar0->Draw();
  h_metdata_llx_tar0->Draw("histsame");//("same e1");
 
  pad_h27->RedrawAxis();
  //--------------------------------
  //all_fitting(h_metsum_llx_tar0);
  // all_fitting(h_metdata_llx_tar0);
  //--------------------------------
  
 // TFile *metxy_tar0_file = new TFile("metxy_tar0.root","recreate");
 // metxy_tar0_file->cd();
 // h_metsum_llx_tar0->Write();
 // h_metdata_llx_tar0->Write();
  
  //-------------------------------Legend---------------------------------------//
  TLegend *leg27 = new TLegend(0.79,0.53,0.87,0.63);//(0.78, 0.63, 0.87, 0.89);
  leg27->SetLineColor(kWhite);
  leg27->SetTextFont(42); 
  leg27->SetTextSize(0.027);
  leg27->SetShadowColor(kWhite); 
  leg27->SetFillColor(kWhite); 
  leg27->AddEntry("h_metsum_llx_tar0","MC METx","ep");
  leg27->AddEntry("h_metdata_llx_tar0","Data METx","ep");
  // leg27->AddEntry("sum","MC MET Wieghted Average","ep");

  leg27->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r27->cd();

  TH1D *h_metdata_llx_tar0_clone2 = (TH1D*) h_metdata_llx_tar0->Clone("h_metdata_llx_tar0_clone2");

  //  overflow(h_metsum_ee_tar0);
  
  //Divide
  h_metdata_llx_tar0_clone2->Divide(h_metdata_llx_tar0_clone,h_metsum_llx_tar0);

  h_metdata_llx_tar0_clone2->GetXaxis()->SetRangeUser(-200,200); 
  h_metdata_llx_tar0_clone2->GetYaxis()->SetRangeUser(0,2);
  h_metdata_llx_tar0_clone2->GetYaxis()->SetNdivisions(4);
  h_metdata_llx_tar0_clone2->GetYaxis()->SetLabelSize(.12);
  h_metdata_llx_tar0_clone2->GetYaxis()->SetTitle("X/Y");
  h_metdata_llx_tar0_clone2->SetLineColor(kBlack);
  h_metdata_llx_tar0_clone2->SetMarkerColor(kBlack);

  pad_r27->SetGridy();

  h_metdata_llx_tar0_clone2->Draw();

  //-------------------------------------------------------------//


 //-----------------------MC vs Data METx TAR2-----------------------------//
  //Canvas
  TCanvas *c28=new TCanvas("c28","Data vs MC 2 Jets METx",800,800);
  TPad *pad_h28 = new TPad("pad_h28","Histo Pad28",0., 0, 1., 0.8);
  TPad *pad_r28 = new TPad("pad_r28","Ratio Pad28",0., 0.8, 1., 1.);

  pad_h28->Draw();
  pad_r28->Draw();

  pad_h28->SetLogy();

  pad_h28->cd();
  
  h_metsum_llx_tar2->GetXaxis()->SetRangeUser(-200,200);
  h_metdata_llx_tar2->GetXaxis()->SetRangeUser(-200,200);
  
  h_metsum_llx_tar2->SetMarkerColor(kRed);
  h_metsum_llx_tar2->SetMarkerStyle(7);
  h_metsum_llx_tar2->SetLineColor(kRed);
  h_metdata_llx_tar2->SetMarkerColor(kBlue);
  h_metdata_llx_tar2->SetMarkerStyle(7); 
  h_metdata_llx_tar2->SetLineColor(kBlue);
  h_metsum_llx_tar2->GetYaxis()->SetTitle("Events / 1 GeV");
  h_metsum_llx_tar2->GetXaxis()->SetTitle("MET_{ll}");
  h_metsum_llx_tar2->SetTitle("Data vs MC METx (2 Jets)");

  h_metsum_llx_tar2->Draw();
  h_metdata_llx_tar2->Draw("histsame");//("same e1");
 
  pad_h28->RedrawAxis();
  //--------------------------------
  //all_fitting(h_metsum_llx_tar2);
  // all_fitting(h_metdata_llx_tar2);
  //--------------------------------
  
 // TFile *metxy_tar2_file = new TFile("metxy_tar2.root","recreate");
 // metxy_tar2_file->cd();
 // h_metsum_llx_tar2->Write();
 // h_metdata_llx_tar2->Write();
  
  //-------------------------------Legend---------------------------------------//
  TLegend *leg28 = new TLegend(0.79,0.53,0.87,0.63);//(0.78, 0.63, 0.87, 0.89);
  leg28->SetLineColor(kWhite);
  leg28->SetTextFont(42); 
  leg28->SetTextSize(0.028);
  leg28->SetShadowColor(kWhite); 
  leg28->SetFillColor(kWhite); 
  leg28->AddEntry("h_metsum_llx_tar2","MC METx","ep");
  leg28->AddEntry("h_metdata_llx_tar2","Data METx","ep");
  // leg28->AddEntry("sum","MC MET Wieghted Average","ep");

  leg28->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r28->cd();

  TH1D *h_metdata_llx_tar2_clone2 = (TH1D*) h_metdata_llx_tar2->Clone("h_metdata_llx_tar2_clone2");

  //  overflow(h_metsum_ee_tar2);
  
  //Divide
  h_metdata_llx_tar2_clone2->Divide(h_metdata_llx_tar2_clone,h_metsum_llx_tar2);

  h_metdata_llx_tar2_clone2->GetXaxis()->SetRangeUser(-200,200); 
  h_metdata_llx_tar2_clone2->GetYaxis()->SetRangeUser(0,2);
  h_metdata_llx_tar2_clone2->GetYaxis()->SetNdivisions(4);
  h_metdata_llx_tar2_clone2->GetYaxis()->SetLabelSize(.12);
  h_metdata_llx_tar2_clone2->GetYaxis()->SetTitle("X/Y");
  h_metdata_llx_tar2_clone2->SetLineColor(kBlack);
  h_metdata_llx_tar2_clone2->SetMarkerColor(kBlack);

  pad_r28->SetGridy();

  h_metdata_llx_tar2_clone2->Draw();

  //-------------------------------------------------------------//
  */

  //---------------------SumET ee INC------------------------------//
  //Canvas
  TCanvas *c29=new TCanvas("c29","SumET ee inc",800,800);
  TPad *pad_h29 = new TPad("pad_h29","Histo Pad29",0., 0, 1., 0.8);
  TPad *pad_r29 = new TPad("pad_r29","Ratio Pad29",0., 0.8, 1., 1.);

  pad_h29->Draw();
  pad_r29->Draw();

  pad_h29->SetLogy();  //c1->pad_h

  pad_h29->cd();  //change pad

  v_sumet_ee_inc[0]->SetTitle("ee+mumu data vs MC SumET (Inclusive)");
  v_sumet_ee_inc[0]->Draw();
  hs_sumet_ee_inc->Draw("histsame");
 
  v_sumet_ee_inc[size-1]->Draw("same e1 x0");

  pad_h29->RedrawAxis();  //c1->pad_h
  
  // double width_mc = find_width(h_sumetsum_ee_inc);
  //double width_data = find_width(v_sumet_ee_inc[size-1]);
  // cout<<"\nINC MC Width = "<<width_mc<<endl;
  // cout<<"\nINC data Width = "<<width_data<<endl;
  //cout<<"\nINC Width difference = "<<width_data-width_mc<<endl;
  //-------------------------------Legend---------------------------------------//
  TLegend *leg29 = new TLegend(0.78, 0.63, 0.87, 0.89);
  leg29->SetLineColor(kWhite);
  leg29->SetTextFont(42); 
  leg29->SetTextSize(0.026);
  leg29->SetShadowColor(kWhite); 
  leg29->SetFillColor(kWhite); 
  leg29->AddEntry("h_sumet_ee_inc_data","data","ep");
  leg29->AddEntry("h_sumet_ee_inc_zjets","Z+jets","F");
  leg29->AddEntry("h_sumet_ee_inc_wz","WZ","F");
  leg29->AddEntry("h_sumet_ee_inc_zz","ZZ","F");
  leg29->AddEntry("h_sumet_ee_inc_ttbar","ttbar","F");
  leg29->AddEntry("h_sumet_ee_inc_singlet","single top","F");
  leg29->AddEntry("h_sumet_ee_inc_ww","WW","F");
  leg29->AddEntry("h_sumet_ee_inc_ttv","ttV","F");
  leg29->AddEntry("h_sumet_ee_inc_vvv","VVV","F");

  leg29->Draw();
  //-------------------------------------------------------------//
  //---------------------Ratio Pad-------------------------------//

  pad_r29->cd();

  TH1F *h_sumet_ee_inc_data_clone = (TH1F*) v_sumet_ee_inc[size-1]->Clone("h_sumet_ee_inc_data_clone");

  //  overflow(h_sumetsum_ee_inc);
  
  //Divide
  h_sumet_ee_inc_data_clone->Divide(h_sumet_ee_inc_data_clone,h_sumetsum_ee_inc);

  h_sumet_ee_inc_data_clone->GetXaxis()->SetRangeUser(0,3000); 
  h_sumet_ee_inc_data_clone->GetYaxis()->SetRangeUser(0,3);
  h_sumet_ee_inc_data_clone->GetYaxis()->SetNdivisions(4);
  h_sumet_ee_inc_data_clone->GetYaxis()->SetLabelSize(.12);
  h_sumet_ee_inc_data_clone->GetYaxis()->SetTitle("data/MC");

  pad_r29->SetGridy();

  h_sumet_ee_inc_data_clone->Draw();

  /*
  TH1D *h_sumet_ee_inc_data_clone_scaled = (TH1D*) h_sumet_ee_inc_data_clone->Clone("h_sumet_ee_inc_data_clone_scaled");

  double k_sumet = h_sumetsum_ee_inc->Integral(0,-1) / h_sumet_ee_inc_data->Integral(0,-1);
  h_sumet_ee_inc_data_clone_scaled->Scale( k_sumet );  
  
  TFile *sumetRatioFile = new TFile("sumet_Ratio.root","recreate");
  sumetRatioFile->cd();
  h_sumet_ee_inc_data_clone->Write();
  h_sumet_ee_inc_data_clone_scaled->Write();
  sumetRatioFile->Close();  
  */

  //-------------------------------------------------------------//
  
  /*
  c1->SaveAs("./3lep_pics/mll_inc.jpg");
  //  c2->SaveAs("./3lep_pics/mll_mumu_inc.jpg");
  c3->SaveAs("./3lep_pics/mll_tar.jpg");
  //  c4->SaveAs("./3lep_pics/mll_mumu_inc.jpg");
  c5->SaveAs("./3lep_pics/met_inc.jpg");
  //  c6->SaveAs("./3lep_pics/mumu_inc.jpg");
  c7->SaveAs("./3lep_pics/met_tar0.jpg");
  //  c8->SaveAs("./3lep_pics/mumu_tar0.jpg");
  c9->SaveAs("./3lep_pics/met_tar2.jpg");
  //  c10->SaveAs("./3lep_pics/mumu_tar2.jpg");
  //  c11->SaveAs("./3lep_pics/ee_vs_mumu_tar0_MC.jpg");
  //  c12->SaveAs("./3lep_pics/ee_vs_mumu_tar0_data.jpg");
  //  c13->SaveAs("./3lep_pics/ee_vs_mumu_tar2_MC.jpg");
  //  c14->SaveAs("./3lep_pics/ee_vs_mumu_tar2_data.jpg");
  //  c15->SaveAs("./3lep_pics/ee_vs_mumu_inc_MC.jpg");
  //  c16->SaveAs("./3lep_pics/ee_vs_mumu_inc_data.jpg");   
  c17->SaveAs("./3lep_pics/nvtx_scaled.jpg");
  c18->SaveAs("./3lep_pics/ee_vs_mumu_inc_data.jpg");
  */
  
  return 0;
}
 
