#include "TFile.h" 
#include "TH1F.h" 
#include "TMath.h"

#include <iostream>
#include <iomanip>

using namespace std;

void range(TH1F *histo, double &value, double &error, double min, double max)
{
  value = histo->IntegralAndError(histo->FindBin(min),histo->FindBin(max)-1,error);	  
}

double err_ratio(double data, double data_err, double MC, double MC_err)
{
  double error=TMath::Sqrt(pow(data_err/MC,2) + pow( data*MC_err/pow(MC,2) ,2) );

  return error;
}


int CorrGetYields()
{

  //-------------------------------------------------Get Yields-----------------------------------------------------------------------..

  TFile *InputFile = new TFile("zjets.root","read");
  TH1F *h_mll_zjets = (TH1F*) InputFile->Get("h_mll_ee_inc")->Clone("h_mll_zjets");
  TH1F *h_mll_zjets2 = (TH1F*) InputFile->Get("h_mll_mumu_inc")->Clone("h_mll_zjets2");
  TH1F *h_mll_tarzjets = (TH1F*) InputFile->Get("h_mll_ee_tar")->Clone("h_mll_tarzjets");
  TH1F *h_mll_tarzjets2 = (TH1F*) InputFile->Get("h_mll_mumu_tar")->Clone("h_mll_tarzjets2");
  TH1F *h_met_ee_inc_zjets = (TH1F*) InputFile->Get("h_met_ee_inc")->Clone("h_met_ee_inc_zjets");
  TH1F *h_met_mumu_inc_zjets = (TH1F*) InputFile->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_zjets");
  TH1F *h_met_ee_tar0_zjets = (TH1F*) InputFile->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_zjets");
  TH1F *h_met_mumu_tar0_zjets = (TH1F*) InputFile->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_zjets");
  TH1F *h_met_ee_tar2_zjets = (TH1F*) InputFile->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_zjets");
  TH1F *h_met_mumu_tar2_zjets = (TH1F*) InputFile->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_zjets");
  TH1F *h_met_emu_inc_zjets = (TH1F*) InputFile->Get("h_met_emu_inc")->Clone("h_met_emu_inc_zjets");
  TH1F *h_met_emu_tar0_zjets = (TH1F*) InputFile->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_zjets");
  TH1F *h_met_emu_tar2_zjets = (TH1F*) InputFile->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_zjets");

  double val_zjets=0.;
  double err_zjets=0.;
  double val_zjets2=0.;
  double err_zjets2=0.;
  double val_tarzjets=0.;
  double err_tarzjets=0.;
  double val_tarzjets2=0.;
  double err_tarzjets2=0.;
  /* double err_met_ee_inc_zjets=0;
  double err_met_mumu_inc_zjets=0;
  double err_met_emu_inc_zjets=0;
  double err_met_ee_tar0_zjets=0;
  double err_met_mumu_tar0_zjets=0;
  double err_met_emu_tar0_zjets=0;
  double err_met_ee_tar2_zjets=0;
  double err_met_mumu_tar2_zjets=0;
  double err_met_emu_tar2_zjets=0;*/

  val_zjets = h_mll_zjets->IntegralAndError(h_mll_zjets->FindBin(81.0),h_mll_zjets->FindBin(101.0)-1,err_zjets);
  val_zjets2 = h_mll_zjets2->IntegralAndError(h_mll_zjets2->FindBin(81.0),h_mll_zjets2->FindBin(101.0)-1,err_zjets2);
  val_tarzjets = h_mll_tarzjets->IntegralAndError(h_mll_tarzjets->FindBin(81.0),h_mll_tarzjets->FindBin(101.0)-1,err_tarzjets);
  val_tarzjets2 = h_mll_tarzjets2->IntegralAndError(h_mll_tarzjets2->FindBin(81.0),h_mll_tarzjets2->FindBin(101.0)-1,err_tarzjets2);

  TFile *InputFile_t = new TFile("singlet.root","read");
  TH1F *h_mll_t = (TH1F*) InputFile_t->Get("h_mll_ee_inc")->Clone("h_mll_t");
  TH1F *h_mll_t2 = (TH1F*) InputFile_t->Get("h_mll_mumu_inc")->Clone("h_mll_t2");
  TH1F *h_mll_tart = (TH1F*) InputFile_t->Get("h_mll_ee_tar")->Clone("h_mll_tart");
  TH1F *h_mll_tart2 = (TH1F*) InputFile_t->Get("h_mll_mumu_tar")->Clone("h_mll_tart2");
  TH1F *h_met_ee_inc_singlet = (TH1F*) InputFile_t->Get("h_met_ee_inc")->Clone("h_met_ee_inc_singlet");
  TH1F *h_met_mumu_inc_singlet = (TH1F*) InputFile_t->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_singlet");
  TH1F *h_met_ee_tar0_singlet = (TH1F*) InputFile_t->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_singlet");
  TH1F *h_met_mumu_tar0_singlet = (TH1F*) InputFile_t->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_singlet");
  TH1F *h_met_ee_tar2_singlet = (TH1F*) InputFile_t->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_singlet");
  TH1F *h_met_mumu_tar2_singlet = (TH1F*) InputFile_t->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_singlet");
  TH1F *h_met_emu_inc_singlet = (TH1F*) InputFile_t->Get("h_met_emu_inc")->Clone("h_met_emu_inc_singlet");
  TH1F *h_met_emu_tar0_singlet = (TH1F*) InputFile_t->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_singlet");
  TH1F *h_met_emu_tar2_singlet = (TH1F*) InputFile_t->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_singlet");

  double val_t=0.;
  double err_t=0.;
  double val_t2=0.;
  double err_t2=0.;
  double val_tart=0.;
  double err_tart=0.;
  double val_tart2=0.;
  double err_tart2=0.;

  val_t = h_mll_t->IntegralAndError(h_mll_t->FindBin(81.0),h_mll_t->FindBin(101.0)-1,err_t);
  val_t2 = h_mll_t2->IntegralAndError(h_mll_t2->FindBin(81.0),h_mll_t2->FindBin(101.0)-1,err_t2);
  val_tart = h_mll_tart->IntegralAndError(h_mll_tart->FindBin(81.0),h_mll_tart->FindBin(101.0)-1,err_tart);
  val_tart2 = h_mll_tart2->IntegralAndError(h_mll_tart2->FindBin(81.0),h_mll_tart2->FindBin(101.0)-1,err_tart2);

  TFile *InputFile_wz = new TFile("wz.root","read");
  TH1F *h_mll_wz = (TH1F*) InputFile_wz->Get("h_mll_ee_inc")->Clone("h_mll_wz"); 
  TH1F *h_mll_wz2 = (TH1F*) InputFile_wz->Get("h_mll_mumu_inc")->Clone("h_mll_wz2"); 
  TH1F *h_mll_tarwz = (TH1F*) InputFile_wz->Get("h_mll_ee_tar")->Clone("h_mll_tarwz");
  TH1F *h_mll_tarwz2 = (TH1F*) InputFile_wz->Get("h_mll_mumu_tar")->Clone("h_mll_tarwz2");
  TH1F *h_met_ee_inc_wz = (TH1F*) InputFile_wz->Get("h_met_ee_inc")->Clone("h_met_ee_inc_wz");
  TH1F *h_met_mumu_inc_wz = (TH1F*) InputFile_wz->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_wz");
  TH1F *h_met_ee_tar0_wz = (TH1F*) InputFile_wz->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_wz");
  TH1F *h_met_mumu_tar0_wz = (TH1F*) InputFile_wz->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_wz");
  TH1F *h_met_ee_tar2_wz = (TH1F*) InputFile_wz->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_wz");
  TH1F *h_met_mumu_tar2_wz = (TH1F*) InputFile_wz->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_wz");
  TH1F *h_met_emu_inc_wz = (TH1F*) InputFile_wz->Get("h_met_emu_inc")->Clone("h_met_emu_inc_wz");
  TH1F *h_met_emu_tar0_wz = (TH1F*) InputFile_wz->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_wz");
  TH1F *h_met_emu_tar2_wz = (TH1F*) InputFile_wz->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_wz");

  double val_wz=0.;
  double err_wz=0.;
  double val_wz2=0.;
  double err_wz2=0.;
  double val_tarwz=0.;
  double err_tarwz=0.;
  double val_tarwz2=0.;
  double err_tarwz2=0.;

  val_wz = h_mll_wz->IntegralAndError(h_mll_wz->FindBin(81.0),h_mll_wz->FindBin(101.0)-1,err_wz);
  val_wz2 = h_mll_wz2->IntegralAndError(h_mll_wz2->FindBin(81.0),h_mll_wz2->FindBin(101.0)-1,err_wz2);
  val_tarwz = h_mll_tarwz->IntegralAndError(h_mll_tarwz->FindBin(81.0),h_mll_tarwz->FindBin(101.0)-1,err_tarwz);
  val_tarwz2 = h_mll_tarwz2->IntegralAndError(h_mll_tarwz2->FindBin(81.0),h_mll_tarwz2->FindBin(101.0)-1,err_tarwz2);

  TFile *InputFile_zz = new TFile("zz.root","read");
  TH1F *h_mll_zz = (TH1F*) InputFile_zz->Get("h_mll_ee_inc")->Clone("h_mll_zz"); 
  TH1F *h_mll_zz2 = (TH1F*) InputFile_zz->Get("h_mll_mumu_inc")->Clone("h_mll_zz2"); 
  TH1F *h_mll_tarzz = (TH1F*) InputFile_zz->Get("h_mll_ee_tar")->Clone("h_mll_tarzz");
  TH1F *h_mll_tarzz2 = (TH1F*) InputFile_zz->Get("h_mll_mumu_tar")->Clone("h_mll_tarzz2");
  TH1F *h_met_ee_inc_zz = (TH1F*) InputFile_zz->Get("h_met_ee_inc")->Clone("h_met_ee_inc_zz");
  TH1F *h_met_mumu_inc_zz = (TH1F*) InputFile_zz->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_zz");
  TH1F *h_met_ee_tar0_zz = (TH1F*) InputFile_zz->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_zz");
  TH1F *h_met_mumu_tar0_zz = (TH1F*) InputFile_zz->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_zz");
  TH1F *h_met_ee_tar2_zz = (TH1F*) InputFile_zz->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_zz");
  TH1F *h_met_mumu_tar2_zz = (TH1F*) InputFile_zz->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_zz");
  TH1F *h_met_emu_inc_zz = (TH1F*) InputFile_zz->Get("h_met_emu_inc")->Clone("h_met_emu_inc_zz");
  TH1F *h_met_emu_tar0_zz = (TH1F*) InputFile_zz->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_zz");
  TH1F *h_met_emu_tar2_zz = (TH1F*) InputFile_zz->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_zz");

  double val_zz=0.;
  double err_zz=0.;
  double val_zz2=0.;
  double err_zz2=0.;
  double val_tarzz=0.;
  double err_tarzz=0.;
  double val_tarzz2=0.;
  double err_tarzz2=0.;

  val_zz = h_mll_zz->IntegralAndError(h_mll_zz->FindBin(81.0),h_mll_zz->FindBin(101.0)-1,err_zz);
  val_zz2 = h_mll_zz2->IntegralAndError(h_mll_zz2->FindBin(81.0),h_mll_zz2->FindBin(101.0)-1,err_zz2);
  val_tarzz = h_mll_tarzz->IntegralAndError(h_mll_tarzz->FindBin(81.0),h_mll_tarzz->FindBin(101.0)-1,err_tarzz);
  val_tarzz2 = h_mll_tarzz2->IntegralAndError(h_mll_tarzz2->FindBin(81.0),h_mll_tarzz2->FindBin(101.0)-1,err_tarzz2);

  TFile *InputFile_ttbar = new TFile("ttbar.root","read");
  TH1F *h_mll_ttbar = (TH1F*) InputFile_ttbar->Get("h_mll_ee_inc")->Clone("h_mll_ttbar");
  TH1F *h_mll_ttbar2 = (TH1F*) InputFile_ttbar->Get("h_mll_mumu_inc")->Clone("h_mll_ttbar2");
  TH1F *h_mll_tarttbar = (TH1F*) InputFile_ttbar->Get("h_mll_ee_tar")->Clone("h_mll_tarttbar");
  TH1F *h_mll_tarttbar2 = (TH1F*) InputFile_ttbar->Get("h_mll_mumu_tar")->Clone("h_mll_tarttbar2");
  TH1F *h_met_ee_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_ee_inc")->Clone("h_met_ee_inc_ttbar");
  TH1F *h_met_mumu_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_ttbar");
  TH1F *h_met_ee_tar0_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_ttbar");
  TH1F *h_met_mumu_tar0_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_ttbar");
  TH1F *h_met_ee_tar2_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_ttbar");
  TH1F *h_met_mumu_tar2_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_ttbar");
  TH1F *h_met_emu_inc_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_emu_inc")->Clone("h_met_emu_inc_ttbar");
  TH1F *h_met_emu_tar0_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_ttbar");
  TH1F *h_met_emu_tar2_ttbar = (TH1F*) InputFile_ttbar->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_ttbar");

  double val_ttbar=0.;
  double err_ttbar=0.;
  double val_ttbar2=0.;
  double err_ttbar2=0.;
  double val_tarttbar=0.;
  double err_tarttbar=0.;
  double val_tarttbar2=0.;
  double err_tarttbar2=0.;

  val_ttbar = h_mll_ttbar->IntegralAndError(h_mll_ttbar->FindBin(81.0),h_mll_ttbar->FindBin(101.0)-1,err_ttbar);
  val_ttbar2 = h_mll_ttbar2->IntegralAndError(h_mll_ttbar2->FindBin(81.0),h_mll_ttbar2->FindBin(101.0)-1,err_ttbar2);
  val_tarttbar = h_mll_tarttbar->IntegralAndError(h_mll_tarttbar->FindBin(81.0),h_mll_tarttbar->FindBin(101.0)-1,err_tarttbar);
  val_tarttbar2 = h_mll_tarttbar2->IntegralAndError(h_mll_tarttbar2->FindBin(81.0),h_mll_tarttbar2->FindBin(101.0)-1,err_tarttbar2);

  TFile *InputFile_ww = new TFile("ww.root","read");
  TH1F *h_mll_ww = (TH1F*) InputFile_ww->Get("h_mll_ee_inc")->Clone("h_mll_ww");
  TH1F *h_mll_ww2 = (TH1F*) InputFile_ww->Get("h_mll_mumu_inc")->Clone("h_mll_ww2");
  TH1F *h_mll_tarww = (TH1F*) InputFile_ww->Get("h_mll_ee_tar")->Clone("h_mll_tarww");
  TH1F *h_mll_tarww2 = (TH1F*) InputFile_ww->Get("h_mll_mumu_tar")->Clone("h_mll_tarww2");
  TH1F *h_met_ee_inc_ww = (TH1F*) InputFile_ww->Get("h_met_ee_inc")->Clone("h_met_ee_inc_ww");
  TH1F *h_met_mumu_inc_ww = (TH1F*) InputFile_ww->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_ww");
  TH1F *h_met_ee_tar0_ww = (TH1F*) InputFile_ww->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_ww");
  TH1F *h_met_mumu_tar0_ww = (TH1F*) InputFile_ww->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_ww");
  TH1F *h_met_ee_tar2_ww = (TH1F*) InputFile_ww->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_ww");
  TH1F *h_met_mumu_tar2_ww = (TH1F*) InputFile_ww->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_ww");
  TH1F *h_met_emu_inc_ww = (TH1F*) InputFile_ww->Get("h_met_emu_inc")->Clone("h_met_emu_inc_ww");
  TH1F *h_met_emu_tar0_ww = (TH1F*) InputFile_ww->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_ww");
  TH1F *h_met_emu_tar2_ww = (TH1F*) InputFile_ww->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_ww");

  double val_ww=0.;
  double err_ww=0.;
  double val_ww2=0.;
  double err_ww2=0.;
  double val_tarww=0.;
  double err_tarww=0.;
  double val_tarww2=0.;
  double err_tarww2=0.;

  val_ww = h_mll_ww->IntegralAndError(h_mll_ww->FindBin(81.0),h_mll_ww->FindBin(101.0)-1,err_ww);
  val_ww2 = h_mll_ww2->IntegralAndError(h_mll_ww2->FindBin(81.0),h_mll_ww2->FindBin(101.0)-1,err_ww2);
  val_tarww = h_mll_tarww->IntegralAndError(h_mll_tarww->FindBin(81.0),h_mll_tarww->FindBin(101.0)-1,err_tarww);
  val_tarww2 = h_mll_tarww2->IntegralAndError(h_mll_tarww2->FindBin(81.0),h_mll_tarww2->FindBin(101.0)-1,err_tarww2);

  TFile *InputFile_ttv = new TFile("ttv.root","read");
  TH1F *h_mll_ttv = (TH1F*) InputFile_ttv->Get("h_mll_ee_inc")->Clone("h_mll_ttv");
  TH1F *h_mll_ttv2 = (TH1F*) InputFile_ttv->Get("h_mll_mumu_inc")->Clone("h_mll_ttv2");
  TH1F *h_mll_tarttv = (TH1F*) InputFile_ttv->Get("h_mll_ee_tar")->Clone("h_mll_tarttv");
  TH1F *h_mll_tarttv2 = (TH1F*) InputFile_ttv->Get("h_mll_mumu_tar")->Clone("h_mll_tarttv2");
  TH1F *h_met_ee_inc_ttv = (TH1F*) InputFile_ttv->Get("h_met_ee_inc")->Clone("h_met_ee_inc_ttv");
  TH1F *h_met_mumu_inc_ttv = (TH1F*) InputFile_ttv->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_ttv");
  TH1F *h_met_ee_tar0_ttv = (TH1F*) InputFile_ttv->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_ttv");
  TH1F *h_met_mumu_tar0_ttv = (TH1F*) InputFile_ttv->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_ttv");
  TH1F *h_met_ee_tar2_ttv = (TH1F*) InputFile_ttv->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_ttv");
  TH1F *h_met_mumu_tar2_ttv = (TH1F*) InputFile_ttv->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_ttv");
  TH1F *h_met_emu_inc_ttv = (TH1F*) InputFile_ttv->Get("h_met_emu_inc")->Clone("h_met_emu_inc_ttv");
  TH1F *h_met_emu_tar0_ttv = (TH1F*) InputFile_ttv->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_ttv");
  TH1F *h_met_emu_tar2_ttv = (TH1F*) InputFile_ttv->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_ttv");

  double val_ttv=0.;
  double err_ttv=0.;
  double val_ttv2=0.;
  double err_ttv2=0.;
  double val_tarttv=0.;
  double err_tarttv=0.;
  double val_tarttv2=0.;
  double err_tarttv2=0.;

  val_ttv = h_mll_ttv->IntegralAndError(h_mll_ttv->FindBin(81.0),h_mll_ttv->FindBin(101.0)-1,err_ttv);
  val_ttv2 = h_mll_ttv2->IntegralAndError(h_mll_ttv2->FindBin(81.0),h_mll_ttv2->FindBin(101.0)-1,err_ttv2);
  val_tarttv = h_mll_tarttv->IntegralAndError(h_mll_tarttv->FindBin(81.0),h_mll_tarttv->FindBin(101.0)-1,err_tarttv);
  val_tarttv2 = h_mll_tarttv2->IntegralAndError(h_mll_tarttv2->FindBin(81.0),h_mll_tarttv2->FindBin(101.0)-1,err_tarttv2);
 
  TFile *InputFile_vvv = new TFile("vvv.root","read");
  TH1F *h_mll_vvv = (TH1F*) InputFile_vvv->Get("h_mll_ee_inc")->Clone("h_mll_vvv");
  TH1F *h_mll_vvv2 = (TH1F*) InputFile_vvv->Get("h_mll_mumu_inc")->Clone("h_mll_vvv2");
  TH1F *h_mll_tarvvv = (TH1F*) InputFile_vvv->Get("h_mll_ee_tar")->Clone("h_mll_tarvvv");
  TH1F *h_mll_tarvvv2 = (TH1F*) InputFile_vvv->Get("h_mll_mumu_tar")->Clone("h_mll_tarvvv2");
  TH1F *h_met_ee_inc_vvv = (TH1F*) InputFile_vvv->Get("h_met_ee_inc")->Clone("h_met_ee_inc_vvv");
  TH1F *h_met_mumu_inc_vvv = (TH1F*) InputFile_vvv->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_vvv");
  TH1F *h_met_ee_tar0_vvv = (TH1F*) InputFile_vvv->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_vvv");
  TH1F *h_met_mumu_tar0_vvv = (TH1F*) InputFile_vvv->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_vvv");
  TH1F *h_met_ee_tar2_vvv = (TH1F*) InputFile_vvv->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_vvv");
  TH1F *h_met_mumu_tar2_vvv = (TH1F*) InputFile_vvv->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_vvv");
  TH1F *h_met_emu_inc_vvv = (TH1F*) InputFile_vvv->Get("h_met_emu_inc")->Clone("h_met_emu_inc_vvv");
  TH1F *h_met_emu_tar0_vvv = (TH1F*) InputFile_vvv->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_vvv");
  TH1F *h_met_emu_tar2_vvv = (TH1F*) InputFile_vvv->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_vvv");

  double val_vvv=0.;
  double err_vvv=0.;
  double val_vvv2=0.;
  double err_vvv2=0.;
  double val_tarvvv=0.;
  double err_tarvvv=0.;
  double val_tarvvv2=0.;
  double err_tarvvv2=0.;

  val_vvv = h_mll_vvv->IntegralAndError(h_mll_vvv->FindBin(81.0),h_mll_vvv->FindBin(101.0)-1,err_vvv);
  val_vvv2 = h_mll_vvv2->IntegralAndError(h_mll_vvv2->FindBin(81.0),h_mll_vvv2->FindBin(101.0)-1,err_vvv2);
  val_tarvvv = h_mll_tarvvv->IntegralAndError(h_mll_tarvvv->FindBin(81.0),h_mll_tarvvv->FindBin(101.0)-1,err_tarvvv);
  val_tarvvv2 = h_mll_tarvvv2->IntegralAndError(h_mll_tarvvv2->FindBin(81.0),h_mll_tarvvv2->FindBin(101.0)-1,err_tarvvv2);

  TFile *InputFile_data = new TFile("data.root","read");
  TH1F *h_mll_data = (TH1F*) InputFile_data->Get("h_mll_ee_inc")->Clone("h_mll_data");
  TH1F *h_mll_data2 = (TH1F*) InputFile_data->Get("h_mll_mumu_inc")->Clone("h_mll_data2");
  TH1F *h_mll_tardata = (TH1F*) InputFile_data->Get("h_mll_ee_tar")->Clone("h_mll_tardata");
  TH1F *h_mll_tardata2 = (TH1F*) InputFile_data->Get("h_mll_mumu_tar")->Clone("h_mll_tardata2");
  TH1F *h_met_ee_inc_data = (TH1F*) InputFile_data->Get("h_met_ee_inc")->Clone("h_met_ee_inc_data");
  TH1F *h_met_mumu_inc_data = (TH1F*) InputFile_data->Get("h_met_mumu_inc")->Clone("h_met_mumu_inc_data");
  TH1F *h_met_ee_tar0_data = (TH1F*) InputFile_data->Get("h_met_ee_tar_njets0")->Clone("h_met_ee_tar0_data");
  TH1F *h_met_mumu_tar0_data = (TH1F*) InputFile_data->Get("h_met_mumu_tar_njets0")->Clone("h_met_mumu_tar0_data");
  TH1F *h_met_ee_tar2_data = (TH1F*) InputFile_data->Get("h_met_ee_tar_njets2")->Clone("h_met_ee_tar2_data");
  TH1F *h_met_mumu_tar2_data = (TH1F*) InputFile_data->Get("h_met_mumu_tar_njets2")->Clone("h_met_mumu_tar2_data");
  TH1F *h_met_emu_inc_data = (TH1F*) InputFile_data->Get("h_met_emu_inc")->Clone("h_met_emu_inc_data");
  TH1F *h_met_emu_tar0_data = (TH1F*) InputFile_data->Get("h_met_emu_tar_njets0")->Clone("h_met_emu_tar0_data");
  TH1F *h_met_emu_tar2_data = (TH1F*) InputFile_data->Get("h_met_emu_tar_njets2")->Clone("h_met_emu_tar2_data");

  double val_data=0.;
  double err_data=0.;
  double val_data2=0.;
  double err_data2=0.;
  double val_tardata=0.;
  double err_tardata=0.;
  double val_tardata2=0.;
  double err_tardata2=0.;

  val_data = h_mll_data->IntegralAndError(h_mll_data->FindBin(81.0),h_mll_data->FindBin(101.0)-1,err_data);
  val_data2 = h_mll_data2->IntegralAndError(h_mll_data2->FindBin(81.0),h_mll_data2->FindBin(101.0)-1,err_data2);
  val_tardata = h_mll_tardata->IntegralAndError(h_mll_tardata->FindBin(81.0),h_mll_tardata->FindBin(101.0)-1,err_tardata);
  val_tardata2 = h_mll_tardata2->IntegralAndError(h_mll_tardata2->FindBin(81.0),h_mll_tardata2->FindBin(101.0)-1,err_tardata2);

  //-----------------------------------------------------------------------------------------------------------------//


  //----------------------------------------Make Inclusive Table-----------------------------------------------------//
  
  int  width=23;
  
  double totalsm = 0;
  double totalsmerr = 0;
  double totalsm2 = 0;
  double totalsmerr2 = 0;

  TH1F *sum = new TH1F("sum","Sum",200,0,200);
  sum->Add(h_mll_zjets);
  sum->Add(h_mll_ttbar);
  sum->Add(h_mll_ww);
  sum->Add(h_mll_wz);
  sum->Add(h_mll_zz);
  sum->Add(h_mll_t);
  sum->Add(h_mll_ttv);
  sum->Add(h_mll_vvv);

  TH1F *sum2 = new TH1F("sum2","Sum2",200,0,200);
  sum2->Add(h_mll_zjets2);
  sum2->Add(h_mll_ttbar2);
  sum2->Add(h_mll_ww2);
  sum2->Add(h_mll_wz2);
  sum2->Add(h_mll_zz2);
  sum2->Add(h_mll_t2);
  sum2->Add(h_mll_ttv2);
  sum2->Add(h_mll_vvv2);

  totalsm = sum->IntegralAndError(sum->FindBin(81.0),sum->FindBin(101.0)-1,totalsmerr);
  totalsm2 = sum2->IntegralAndError(sum2->FindBin(81.0),sum2->FindBin(101.0)-1,totalsmerr2);
 
  cout<<"\n"<<setw(2*width)<<"Inclusive"<<endl;
  cout<<"\n"<<setw(width)<<"Sample"<<setw(width)<<"ee"<<setw(width)<<"MuMu"<<endl;

  cout<<"\n"<<setw(width)<<"Z+jets"<<setw(width)<<Form("%.2f +/- %.2f",val_zjets, err_zjets)<<setw(width)<<Form("%.2f +/- %.2f",val_zjets2, err_zjets2)<<endl;
  cout<<"\n"<<setw(width)<<"ttbar"<<setw(width)<<Form("%.2f +/- %.2f",val_ttbar, err_ttbar)<<setw(width)<<Form("%.2f +/- %.2f",val_ttbar2, err_ttbar2)<<endl; 
  cout<<"\n"<<setw(width)<<"WW"<<setw(width)<<Form("%.2f +/- %.2f",val_ww, err_ww)<<setw(width)<<Form("%.2f +/- %.2f",val_ww2, err_ww2)<<endl; 
  cout<<"\n"<<setw(width)<<"WZ"<<setw(width)<<Form("%.2f +/- %.2f",val_wz, err_wz)<<setw(width)<<Form("%.2f +/- %.2f",val_wz2, err_wz2)<<endl; 
  cout<<"\n"<<setw(width)<<"ZZ"<<setw(width)<<Form("%.2f +/- %.2f",val_zz, err_zz)<<setw(width)<<Form("%.2f +/- %.2f",val_zz2, err_zz2)<<endl; 
  cout<<"\n"<<setw(width)<<"single top"<<setw(width)<<Form("%.2f +/- %.2f",val_t, err_t)<<setw(width)<<Form("%.2f +/- %.2f",val_t2, err_t2)<<endl; 
  cout<<"\n"<<setw(width)<<"ttbarV"<<setw(width)<<Form("%.2f +/- %.2f",val_ttv, err_ttv)<<setw(width)<<Form("%.2f +/- %.2f",val_ttv2, err_ttv2)<<endl; 
  cout<<"\n"<<setw(width)<<"VVV"<<setw(width)<<Form("%.2f +/- %.2f",val_vvv, err_vvv)<<setw(width)<<Form("%.2f +/- %.2f",val_vvv2, err_vvv2)<<endl; 
  cout<<"\n"<<setw(width)<<"total SM MC"<<setw(width)<<Form("%.2f +/- %.2f",totalsm, totalsmerr)<<setw(width)<<Form("%.2f +/- %.2f",totalsm2, totalsmerr2)<<endl;
  cout<<"\n"<<setw(width)<<"data"<<setw(width)<<Form("%.2f +/- %.2f",val_data, err_data)<<setw(width)<<Form("%.2f +/- %.2f",val_data2, err_data2)<<endl; 

  //-----------------------------------------------------------------------------------------------------------------//

  //----------------------------------------Make Targeted Table-----------------------------------------------------//
 
  double totalsm_tar = 0;
  double totalsmerr_tar = 0;
  double totalsm_tar2 = 0;
  double totalsmerr_tar2 = 0;

  TH1F *tarsum = new TH1F("tarsum","Targeted Sum",200,0,200);
  tarsum->Add(h_mll_tarzjets);
  tarsum->Add(h_mll_tarttbar);
  tarsum->Add(h_mll_tarww);
  tarsum->Add(h_mll_tarwz);
  tarsum->Add(h_mll_tarzz);
  tarsum->Add(h_mll_tart);
  tarsum->Add(h_mll_tarttv);
  tarsum->Add(h_mll_tarvvv);

  TH1F *tarsum2 = new TH1F("tarsum2","Targeted Sum2",200,0,200);
  tarsum2->Add(h_mll_tarzjets2);
  tarsum2->Add(h_mll_tarttbar2);
  tarsum2->Add(h_mll_tarww2);
  tarsum2->Add(h_mll_tarwz2);
  tarsum2->Add(h_mll_tarzz2);
  tarsum2->Add(h_mll_tart2);
  tarsum2->Add(h_mll_tarttv2);
  tarsum2->Add(h_mll_tarvvv2);

  totalsm_tar = tarsum->IntegralAndError(tarsum->FindBin(81.0),tarsum->FindBin(101.0)-1,totalsmerr_tar);
  totalsm_tar2 = tarsum2->IntegralAndError(tarsum2->FindBin(81.0),tarsum2->FindBin(101.0)-1,totalsmerr_tar2);
 
  cout<<"\n\n"<<setw(2*width)<<"Targeted"<<endl;
  cout<<"\n"<<setw(width)<<"Sample"<<setw(width)<<"ee"<<setw(width)<<"MuMu"<<endl;

  cout<<"\n"<<setw(width)<<"Z+jets"<<setw(width)<<Form("%.2f +/- %.2f",val_tarzjets, err_tarzjets)<<setw(width)<<Form("%.2f +/- %.2f",val_tarzjets2, err_tarzjets2)<<endl;
  cout<<"\n"<<setw(width)<<"ttbar"<<setw(width)<<Form("%.2f +/- %.2f",val_tarttbar, err_tarttbar)<<setw(width)<<Form("%.2f +/- %.2f",val_tarttbar2, err_tarttbar2)<<endl; 
  cout<<"\n"<<setw(width)<<"WW"<<setw(width)<<Form("%.2f +/- %.2f",val_tarww, err_tarww)<<setw(width)<<Form("%.2f +/- %.2f",val_tarww2, err_tarww2)<<endl; 
  cout<<"\n"<<setw(width)<<"WZ"<<setw(width)<<Form("%.2f +/- %.2f",val_tarwz, err_tarwz)<<setw(width)<<Form("%.2f +/- %.2f",val_tarwz2, err_tarwz2)<<endl; 
  cout<<"\n"<<setw(width)<<"ZZ"<<setw(width)<<Form("%.2f +/- %.2f",val_tarzz, err_tarzz)<<setw(width)<<Form("%.2f +/- %.2f",val_tarzz2, err_tarzz2)<<endl; 
  cout<<"\n"<<setw(width)<<"single top"<<setw(width)<<Form("%.2f +/- %.2f",val_tart, err_tart)<<setw(width)<<Form("%.2f +/- %.2f",val_tart2, err_tart2)<<endl; 
  cout<<"\n"<<setw(width)<<"ttbarV"<<setw(width)<<Form("%.2f +/- %.2f",val_tarttv, err_tarttv)<<setw(width)<<Form("%.2f +/- %.2f",val_tarttv2, err_tarttv2)<<endl; 
  cout<<"\n"<<setw(width)<<"VVV"<<setw(width)<<Form("%.2f +/- %.2f",val_tarvvv, err_tarvvv)<<setw(width)<<Form("%.2f +/- %.2f",val_tarvvv2, err_tarvvv2)<<endl; 
  cout<<"\n"<<setw(width)<<"total SM MC"<<setw(width)<<Form("%.2f +/- %.2f",totalsm_tar, totalsmerr_tar)<<setw(width)<<Form("%.2f +/- %.2f",totalsm_tar2, totalsmerr_tar2)<<endl;
  cout<<"\n"<<setw(width)<<"data"<<setw(width)<<Form("%.2f +/- %.2f",val_tardata, err_tardata)<<setw(width)<<Form("%.2f +/- %.2f",val_tardata2, err_tardata2)<<endl; 
  
  //-----------------------------------------------------------------------------------------------------------------//

  //-----------------------------------------------------------------------------------------------------------------//
  //----------------------------------------------------MET----------------------------------------------------------//
  //-----------------------------------------------------------------------------------------------------------------//

  double Ree= 0.440974;
  double Rmumu= 0.524255;

  //------------------------------------------met ee+mumu inc--------------------------------------------------------------//
  h_met_ee_inc_zjets->Add(h_met_emu_inc_zjets,- Ree);
  h_met_ee_inc_ttbar->Add(h_met_emu_inc_ttbar,- Ree);
  h_met_ee_inc_ww->Add(h_met_emu_inc_ww,- Ree);
  h_met_ee_inc_wz->Add(h_met_emu_inc_wz,- Ree);
  h_met_ee_inc_zz->Add(h_met_emu_inc_zz,- Ree);
  h_met_ee_inc_singlet->Add(h_met_emu_inc_singlet,- Ree);
  h_met_ee_inc_ttv->Add(h_met_emu_inc_ttv,- Ree);
  h_met_ee_inc_vvv->Add(h_met_emu_inc_vvv,- Ree);
  h_met_ee_inc_data->Add(h_met_emu_inc_data,- Ree);

  h_met_mumu_inc_zjets->Add(h_met_emu_inc_zjets,- Rmumu);
  h_met_mumu_inc_ttbar->Add(h_met_emu_inc_ttbar,- Rmumu);
  h_met_mumu_inc_ww->Add(h_met_emu_inc_ww,- Rmumu);
  h_met_mumu_inc_wz->Add(h_met_emu_inc_wz,- Rmumu);
  h_met_mumu_inc_zz->Add(h_met_emu_inc_zz,- Rmumu);
  h_met_mumu_inc_singlet->Add(h_met_emu_inc_singlet,- Rmumu);
  h_met_mumu_inc_ttv->Add(h_met_emu_inc_ttv,- Rmumu);
  h_met_mumu_inc_vvv->Add(h_met_emu_inc_vvv,- Rmumu);
  h_met_mumu_inc_data->Add(h_met_emu_inc_data,- Rmumu);
  
  //Merge ee and mumu
  h_met_ee_inc_zjets->Add(h_met_mumu_inc_zjets);
  h_met_ee_inc_ttbar->Add(h_met_mumu_inc_ttbar);
  h_met_ee_inc_ww->Add(h_met_mumu_inc_ww);
  h_met_ee_inc_wz->Add(h_met_mumu_inc_wz);
  h_met_ee_inc_zz->Add(h_met_mumu_inc_zz);
  h_met_ee_inc_singlet->Add(h_met_mumu_inc_singlet);
  h_met_ee_inc_ttv->Add(h_met_mumu_inc_ttv);
  h_met_ee_inc_vvv->Add(h_met_mumu_inc_vvv);
  h_met_ee_inc_data->Add(h_met_mumu_inc_data);

  TH1F *h_metsum_ee_inc = new TH1F("h_metsum_ee_inc","Met ee inc",350,0,350);
  h_metsum_ee_inc->Add(h_met_ee_inc_zjets);
  h_metsum_ee_inc->Add(h_met_ee_inc_ttbar);
  h_metsum_ee_inc->Add(h_met_ee_inc_ww);
  h_metsum_ee_inc->Add(h_met_ee_inc_wz);
  h_metsum_ee_inc->Add(h_met_ee_inc_zz);
  h_metsum_ee_inc->Add(h_met_ee_inc_singlet);
  h_metsum_ee_inc->Add(h_met_ee_inc_ttv);
  h_metsum_ee_inc->Add(h_met_ee_inc_vvv);
  
  double  err1_metsum_ee_inc=0;
  double  val1_metsum_ee_inc=0;
  double  err2_metsum_ee_inc=0;
  double  val2_metsum_ee_inc=0;
  double  err3_metsum_ee_inc=0;
  double  val3_metsum_ee_inc=0;
  double  err4_metsum_ee_inc=0;
  double  val4_metsum_ee_inc=0;
  double  err5_metsum_ee_inc=0;
  double  val5_metsum_ee_inc=0;
  double  err6_metsum_ee_inc=0;
  double  val6_metsum_ee_inc=0;
  double  err7_metsum_ee_inc=0;
  double  val7_metsum_ee_inc=0;
  double  err8_metsum_ee_inc=0;
  double  val8_metsum_ee_inc=0;
 
  range(h_metsum_ee_inc,val1_metsum_ee_inc,err1_metsum_ee_inc,0,20);
  range(h_metsum_ee_inc,val2_metsum_ee_inc,err2_metsum_ee_inc,20,40);
  range(h_metsum_ee_inc,val3_metsum_ee_inc,err3_metsum_ee_inc,40,60);
  range(h_metsum_ee_inc,val4_metsum_ee_inc,err4_metsum_ee_inc,60,80);
  range(h_metsum_ee_inc,val5_metsum_ee_inc,err5_metsum_ee_inc,80,100);
  range(h_metsum_ee_inc,val6_metsum_ee_inc,err6_metsum_ee_inc,100,150);
  range(h_metsum_ee_inc,val7_metsum_ee_inc,err7_metsum_ee_inc,150,200);
  range(h_metsum_ee_inc,val8_metsum_ee_inc,err8_metsum_ee_inc,200,-1);

  double  err1_data_ee_inc=0;
  double  val1_data_ee_inc=0;
  double  err2_data_ee_inc=0;
  double  val2_data_ee_inc=0;
  double  err3_data_ee_inc=0;
  double  val3_data_ee_inc=0;
  double  err4_data_ee_inc=0;
  double  val4_data_ee_inc=0;
  double  err5_data_ee_inc=0;
  double  val5_data_ee_inc=0;
  double  err6_data_ee_inc=0;
  double  val6_data_ee_inc=0;
  double  err7_data_ee_inc=0;
  double  val7_data_ee_inc=0;
  double  err8_data_ee_inc=0;
  double  val8_data_ee_inc=0;

  range(h_met_ee_inc_data,val1_data_ee_inc,err1_data_ee_inc,0,20);
  range(h_met_ee_inc_data,val2_data_ee_inc,err2_data_ee_inc,20,40);
  range(h_met_ee_inc_data,val3_data_ee_inc,err3_data_ee_inc,40,60);
  range(h_met_ee_inc_data,val4_data_ee_inc,err4_data_ee_inc,60,80);
  range(h_met_ee_inc_data,val5_data_ee_inc,err5_data_ee_inc,80,100);
  range(h_met_ee_inc_data,val6_data_ee_inc,err6_data_ee_inc,100,150);
  range(h_met_ee_inc_data,val7_data_ee_inc,err7_data_ee_inc,150,200);
  range(h_met_ee_inc_data,val8_data_ee_inc,err8_data_ee_inc,200,-1);

  double err1_ratio_ee_inc=0;
  double err2_ratio_ee_inc=0;
  double err3_ratio_ee_inc=0;
  double err4_ratio_ee_inc=0;
  double err5_ratio_ee_inc=0;
  double err6_ratio_ee_inc=0;
  double err7_ratio_ee_inc=0;
  double err8_ratio_ee_inc=0;

  err1_ratio_ee_inc=err_ratio(val1_data_ee_inc, err1_data_ee_inc, val1_metsum_ee_inc, err1_metsum_ee_inc);
  err2_ratio_ee_inc=err_ratio(val2_data_ee_inc, err2_data_ee_inc, val2_metsum_ee_inc, err2_metsum_ee_inc);
  err3_ratio_ee_inc=err_ratio(val3_data_ee_inc, err3_data_ee_inc, val3_metsum_ee_inc, err3_metsum_ee_inc);
  err4_ratio_ee_inc=err_ratio(val4_data_ee_inc, err4_data_ee_inc, val4_metsum_ee_inc, err4_metsum_ee_inc);
  err5_ratio_ee_inc=err_ratio(val5_data_ee_inc, err5_data_ee_inc, val5_metsum_ee_inc, err5_metsum_ee_inc);
  err6_ratio_ee_inc=err_ratio(val6_data_ee_inc, err6_data_ee_inc, val6_metsum_ee_inc, err6_metsum_ee_inc);
  err7_ratio_ee_inc=err_ratio(val7_data_ee_inc, err7_data_ee_inc, val7_metsum_ee_inc, err7_metsum_ee_inc);
  err8_ratio_ee_inc=err_ratio(val8_data_ee_inc, err8_data_ee_inc, val8_metsum_ee_inc, err8_metsum_ee_inc);

  //----------------------------------------------------------------------------------------------//
  /*
  //--------------------------------------met mumu inc--------------------------------------------//
  h_met_mumu_inc_zjets->Add(h_met_emu_inc_zjets,- Rmumu);
  h_met_mumu_inc_ttbar->Add(h_met_emu_inc_ttbar,- Rmumu);
  h_met_mumu_inc_ww->Add(h_met_emu_inc_ww,- Rmumu);
  h_met_mumu_inc_wz->Add(h_met_emu_inc_wz,- Rmumu);
  h_met_mumu_inc_zz->Add(h_met_emu_inc_zz,- Rmumu);
  h_met_mumu_inc_singlet->Add(h_met_emu_inc_singlet,- Rmumu);
  h_met_mumu_inc_ttv->Add(h_met_emu_inc_ttv,- Rmumu);
  h_met_mumu_inc_vvv->Add(h_met_emu_inc_vvv,- Rmumu);
  h_met_mumu_inc_data->Add(h_met_emu_inc_data,- Rmumu);
  
  TH1F *h_metsum_mumu_inc = new TH1F("h_metsum_mumu_inc","Met mumu inc",350,0,350);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_zjets);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_ttbar);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_ww);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_wz);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_zz);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_singlet);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_ttv);
  h_metsum_mumu_inc->Add(h_met_mumu_inc_vvv);
  
  double  err1_metsum_mumu_inc=0;
  double  val1_metsum_mumu_inc=0;
  double  err2_metsum_mumu_inc=0;
  double  val2_metsum_mumu_inc=0;
  double  err3_metsum_mumu_inc=0;
  double  val3_metsum_mumu_inc=0;
  double  err4_metsum_mumu_inc=0;
  double  val4_metsum_mumu_inc=0;
  double  err5_metsum_mumu_inc=0;
  double  val5_metsum_mumu_inc=0;
  double  err6_metsum_mumu_inc=0;
  double  val6_metsum_mumu_inc=0;
  double  err7_metsum_mumu_inc=0;
  double  val7_metsum_mumu_inc=0;
  double  err8_metsum_mumu_inc=0;
  double  val8_metsum_mumu_inc=0;
 
  range(h_metsum_mumu_inc,val1_metsum_mumu_inc,err1_metsum_mumu_inc,0,20);
  range(h_metsum_mumu_inc,val2_metsum_mumu_inc,err2_metsum_mumu_inc,20,40);
  range(h_metsum_mumu_inc,val3_metsum_mumu_inc,err3_metsum_mumu_inc,40,60);
  range(h_metsum_mumu_inc,val4_metsum_mumu_inc,err4_metsum_mumu_inc,60,80);
  range(h_metsum_mumu_inc,val5_metsum_mumu_inc,err5_metsum_mumu_inc,80,100);
  range(h_metsum_mumu_inc,val6_metsum_mumu_inc,err6_metsum_mumu_inc,100,150);
  range(h_metsum_mumu_inc,val7_metsum_mumu_inc,err7_metsum_mumu_inc,150,200);
  range(h_metsum_mumu_inc,val8_metsum_mumu_inc,err8_metsum_mumu_inc,200,-1);

  double  err1_data_mumu_inc=0;
  double  val1_data_mumu_inc=0;
  double  err2_data_mumu_inc=0;
  double  val2_data_mumu_inc=0;
  double  err3_data_mumu_inc=0;
  double  val3_data_mumu_inc=0;
  double  err4_data_mumu_inc=0;
  double  val4_data_mumu_inc=0;
  double  err5_data_mumu_inc=0;
  double  val5_data_mumu_inc=0;
  double  err6_data_mumu_inc=0;
  double  val6_data_mumu_inc=0;
  double  err7_data_mumu_inc=0;
  double  val7_data_mumu_inc=0;
  double  err8_data_mumu_inc=0;
  double  val8_data_mumu_inc=0;

  range(h_met_mumu_inc_data,val1_data_mumu_inc,err1_data_mumu_inc,0,20);
  range(h_met_mumu_inc_data,val2_data_mumu_inc,err2_data_mumu_inc,20,40);
  range(h_met_mumu_inc_data,val3_data_mumu_inc,err3_data_mumu_inc,40,60);
  range(h_met_mumu_inc_data,val4_data_mumu_inc,err4_data_mumu_inc,60,80);
  range(h_met_mumu_inc_data,val5_data_mumu_inc,err5_data_mumu_inc,80,100);
  range(h_met_mumu_inc_data,val6_data_mumu_inc,err6_data_mumu_inc,100,150);
  range(h_met_mumu_inc_data,val7_data_mumu_inc,err7_data_mumu_inc,150,200);
  range(h_met_mumu_inc_data,val8_data_mumu_inc,err8_data_mumu_inc,200,-1);

  double err1_ratio_mumu_inc=0;
  double err2_ratio_mumu_inc=0;
  double err3_ratio_mumu_inc=0;
  double err4_ratio_mumu_inc=0;
  double err5_ratio_mumu_inc=0;
  double err6_ratio_mumu_inc=0;
  double err7_ratio_mumu_inc=0;
  double err8_ratio_mumu_inc=0;

  err1_ratio_mumu_inc=err_ratio(val1_data_mumu_inc, err1_data_mumu_inc, val1_metsum_mumu_inc, err1_metsum_mumu_inc);
  err2_ratio_mumu_inc=err_ratio(val2_data_mumu_inc, err2_data_mumu_inc, val2_metsum_mumu_inc, err2_metsum_mumu_inc);
  err3_ratio_mumu_inc=err_ratio(val3_data_mumu_inc, err3_data_mumu_inc, val3_metsum_mumu_inc, err3_metsum_mumu_inc);
  err4_ratio_mumu_inc=err_ratio(val4_data_mumu_inc, err4_data_mumu_inc, val4_metsum_mumu_inc, err4_metsum_mumu_inc);
  err5_ratio_mumu_inc=err_ratio(val5_data_mumu_inc, err5_data_mumu_inc, val5_metsum_mumu_inc, err5_metsum_mumu_inc);
  err6_ratio_mumu_inc=err_ratio(val6_data_mumu_inc, err6_data_mumu_inc, val6_metsum_mumu_inc, err6_metsum_mumu_inc);
  err7_ratio_mumu_inc=err_ratio(val7_data_mumu_inc, err7_data_mumu_inc, val7_metsum_mumu_inc, err7_metsum_mumu_inc);
  err8_ratio_mumu_inc=err_ratio(val8_data_mumu_inc, err8_data_mumu_inc, val8_metsum_mumu_inc, err8_metsum_mumu_inc);
  //----------------------------------------------------------------------------------------------//
  */
  //------------------------------------------met ee+mumu tar0--------------------------------------------------------------//
  h_met_ee_tar0_zjets->Add(h_met_emu_tar0_zjets,- Ree);
  h_met_ee_tar0_ttbar->Add(h_met_emu_tar0_ttbar,- Ree);
  h_met_ee_tar0_ww->Add(h_met_emu_tar0_ww,- Ree);
  h_met_ee_tar0_wz->Add(h_met_emu_tar0_wz,- Ree);
  h_met_ee_tar0_zz->Add(h_met_emu_tar0_zz,- Ree);
  h_met_ee_tar0_singlet->Add(h_met_emu_tar0_singlet,- Ree);
  h_met_ee_tar0_ttv->Add(h_met_emu_tar0_ttv,- Ree);
  h_met_ee_tar0_vvv->Add(h_met_emu_tar0_vvv,- Ree);
  h_met_ee_tar0_data->Add(h_met_emu_tar0_data,- Ree);

  h_met_mumu_tar0_zjets->Add(h_met_emu_tar0_zjets,- Rmumu);
  h_met_mumu_tar0_ttbar->Add(h_met_emu_tar0_ttbar,- Rmumu);
  h_met_mumu_tar0_ww->Add(h_met_emu_tar0_ww,- Rmumu);
  h_met_mumu_tar0_wz->Add(h_met_emu_tar0_wz,- Rmumu);
  h_met_mumu_tar0_zz->Add(h_met_emu_tar0_zz,- Rmumu);
  h_met_mumu_tar0_singlet->Add(h_met_emu_tar0_singlet,- Rmumu);
  h_met_mumu_tar0_ttv->Add(h_met_emu_tar0_ttv,- Rmumu);
  h_met_mumu_tar0_vvv->Add(h_met_emu_tar0_vvv,- Rmumu);
  h_met_mumu_tar0_data->Add(h_met_emu_tar0_data,- Rmumu);

  h_met_ee_tar0_zjets->Add(h_met_mumu_tar0_zjets);
  h_met_ee_tar0_ttbar->Add(h_met_mumu_tar0_ttbar);
  h_met_ee_tar0_ww->Add(h_met_mumu_tar0_ww);
  h_met_ee_tar0_wz->Add(h_met_mumu_tar0_wz);
  h_met_ee_tar0_zz->Add(h_met_mumu_tar0_zz);
  h_met_ee_tar0_singlet->Add(h_met_mumu_tar0_singlet);
  h_met_ee_tar0_ttv->Add(h_met_mumu_tar0_ttv);
  h_met_ee_tar0_vvv->Add(h_met_mumu_tar0_vvv);
  h_met_ee_tar0_data->Add(h_met_mumu_tar0_data);
  
  TH1F *h_metsum_ee_tar0 = new TH1F("h_metsum_ee_tar0","Met ee+mumu tar0",350,0,350);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_zjets);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_ttbar);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_ww);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_wz);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_zz);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_singlet);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_ttv);
  h_metsum_ee_tar0->Add(h_met_ee_tar0_vvv);
  
  double  err1_metsum_ee_tar0=0;
  double  val1_metsum_ee_tar0=0;
  double  err2_metsum_ee_tar0=0;
  double  val2_metsum_ee_tar0=0;
  double  err3_metsum_ee_tar0=0;
  double  val3_metsum_ee_tar0=0;
  double  err4_metsum_ee_tar0=0;
  double  val4_metsum_ee_tar0=0;
  double  err5_metsum_ee_tar0=0;
  double  val5_metsum_ee_tar0=0;
  double  err6_metsum_ee_tar0=0;
  double  val6_metsum_ee_tar0=0;
  double  err7_metsum_ee_tar0=0;
  double  val7_metsum_ee_tar0=0;
  double  err8_metsum_ee_tar0=0;
  double  val8_metsum_ee_tar0=0;
 
  range(h_metsum_ee_tar0,val1_metsum_ee_tar0,err1_metsum_ee_tar0,0,20);
  range(h_metsum_ee_tar0,val2_metsum_ee_tar0,err2_metsum_ee_tar0,20,40);
  range(h_metsum_ee_tar0,val3_metsum_ee_tar0,err3_metsum_ee_tar0,40,60);
  range(h_metsum_ee_tar0,val4_metsum_ee_tar0,err4_metsum_ee_tar0,60,80);
  range(h_metsum_ee_tar0,val5_metsum_ee_tar0,err5_metsum_ee_tar0,80,100);
  range(h_metsum_ee_tar0,val6_metsum_ee_tar0,err6_metsum_ee_tar0,100,150);
  range(h_metsum_ee_tar0,val7_metsum_ee_tar0,err7_metsum_ee_tar0,150,200);
  range(h_metsum_ee_tar0,val8_metsum_ee_tar0,err8_metsum_ee_tar0,200,-1);

  double  err1_data_ee_tar0=0;
  double  val1_data_ee_tar0=0;
  double  err2_data_ee_tar0=0;
  double  val2_data_ee_tar0=0;
  double  err3_data_ee_tar0=0;
  double  val3_data_ee_tar0=0;
  double  err4_data_ee_tar0=0;
  double  val4_data_ee_tar0=0;
  double  err5_data_ee_tar0=0;
  double  val5_data_ee_tar0=0;
  double  err6_data_ee_tar0=0;
  double  val6_data_ee_tar0=0;
  double  err7_data_ee_tar0=0;
  double  val7_data_ee_tar0=0;
  double  err8_data_ee_tar0=0;
  double  val8_data_ee_tar0=0;

  range(h_met_ee_tar0_data,val1_data_ee_tar0,err1_data_ee_tar0,0,20);
  range(h_met_ee_tar0_data,val2_data_ee_tar0,err2_data_ee_tar0,20,40);
  range(h_met_ee_tar0_data,val3_data_ee_tar0,err3_data_ee_tar0,40,60);
  range(h_met_ee_tar0_data,val4_data_ee_tar0,err4_data_ee_tar0,60,80);
  range(h_met_ee_tar0_data,val5_data_ee_tar0,err5_data_ee_tar0,80,100);
  range(h_met_ee_tar0_data,val6_data_ee_tar0,err6_data_ee_tar0,100,150);
  range(h_met_ee_tar0_data,val7_data_ee_tar0,err7_data_ee_tar0,150,200);
  range(h_met_ee_tar0_data,val8_data_ee_tar0,err8_data_ee_tar0,200,-1);

  double err1_ratio_ee_tar0=0;
  double err2_ratio_ee_tar0=0;
  double err3_ratio_ee_tar0=0;
  double err4_ratio_ee_tar0=0;
  double err5_ratio_ee_tar0=0;
  double err6_ratio_ee_tar0=0;
  double err7_ratio_ee_tar0=0;
  double err8_ratio_ee_tar0=0;

  err1_ratio_ee_tar0=err_ratio(val1_data_ee_tar0, err1_data_ee_tar0, val1_metsum_ee_tar0, err1_metsum_ee_tar0);
  err2_ratio_ee_tar0=err_ratio(val2_data_ee_tar0, err2_data_ee_tar0, val2_metsum_ee_tar0, err2_metsum_ee_tar0);
  err3_ratio_ee_tar0=err_ratio(val3_data_ee_tar0, err3_data_ee_tar0, val3_metsum_ee_tar0, err3_metsum_ee_tar0);
  err4_ratio_ee_tar0=err_ratio(val4_data_ee_tar0, err4_data_ee_tar0, val4_metsum_ee_tar0, err4_metsum_ee_tar0);
  err5_ratio_ee_tar0=err_ratio(val5_data_ee_tar0, err5_data_ee_tar0, val5_metsum_ee_tar0, err5_metsum_ee_tar0);
  err6_ratio_ee_tar0=err_ratio(val6_data_ee_tar0, err6_data_ee_tar0, val6_metsum_ee_tar0, err6_metsum_ee_tar0);
  err7_ratio_ee_tar0=err_ratio(val7_data_ee_tar0, err7_data_ee_tar0, val7_metsum_ee_tar0, err7_metsum_ee_tar0);
  err8_ratio_ee_tar0=err_ratio(val8_data_ee_tar0, err8_data_ee_tar0, val8_metsum_ee_tar0, err8_metsum_ee_tar0);
  //----------------------------------------------------------------------------------------------//
  /*
  //--------------------------------------met mumu tar0--------------------------------------------//
  h_met_mumu_tar0_zjets->Add(h_met_emu_tar0_zjets,- Rmumu);
  h_met_mumu_tar0_ttbar->Add(h_met_emu_tar0_ttbar,- Rmumu);
  h_met_mumu_tar0_ww->Add(h_met_emu_tar0_ww,- Rmumu);
  h_met_mumu_tar0_wz->Add(h_met_emu_tar0_wz,- Rmumu);
  h_met_mumu_tar0_zz->Add(h_met_emu_tar0_zz,- Rmumu);
  h_met_mumu_tar0_singlet->Add(h_met_emu_tar0_singlet,- Rmumu);
  h_met_mumu_tar0_ttv->Add(h_met_emu_tar0_ttv,- Rmumu);
  h_met_mumu_tar0_vvv->Add(h_met_emu_tar0_vvv,- Rmumu);
  h_met_mumu_tar0_data->Add(h_met_emu_tar0_data,- Rmumu);
  
  TH1F *h_metsum_mumu_tar0 = new TH1F("h_metsum_mumu_tar0","Met mumu tar0",350,0,350);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_zjets);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_ttbar);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_ww);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_wz);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_zz);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_singlet);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_ttv);
  h_metsum_mumu_tar0->Add(h_met_mumu_tar0_vvv);
  
  double  err1_metsum_mumu_tar0=0;
  double  val1_metsum_mumu_tar0=0;
  double  err2_metsum_mumu_tar0=0;
  double  val2_metsum_mumu_tar0=0;
  double  err3_metsum_mumu_tar0=0;
  double  val3_metsum_mumu_tar0=0;
  double  err4_metsum_mumu_tar0=0;
  double  val4_metsum_mumu_tar0=0;
  double  err5_metsum_mumu_tar0=0;
  double  val5_metsum_mumu_tar0=0;
  double  err6_metsum_mumu_tar0=0;
  double  val6_metsum_mumu_tar0=0;
  double  err7_metsum_mumu_tar0=0;
  double  val7_metsum_mumu_tar0=0;
  double  err8_metsum_mumu_tar0=0;
  double  val8_metsum_mumu_tar0=0;
 
  range(h_metsum_mumu_tar0,val1_metsum_mumu_tar0,err1_metsum_mumu_tar0,0,20);
  range(h_metsum_mumu_tar0,val2_metsum_mumu_tar0,err2_metsum_mumu_tar0,20,40);
  range(h_metsum_mumu_tar0,val3_metsum_mumu_tar0,err3_metsum_mumu_tar0,40,60);
  range(h_metsum_mumu_tar0,val4_metsum_mumu_tar0,err4_metsum_mumu_tar0,60,80);
  range(h_metsum_mumu_tar0,val5_metsum_mumu_tar0,err5_metsum_mumu_tar0,80,100);
  range(h_metsum_mumu_tar0,val6_metsum_mumu_tar0,err6_metsum_mumu_tar0,100,150);
  range(h_metsum_mumu_tar0,val7_metsum_mumu_tar0,err7_metsum_mumu_tar0,150,200);
  range(h_metsum_mumu_tar0,val8_metsum_mumu_tar0,err8_metsum_mumu_tar0,200,-1);

  double  err1_data_mumu_tar0=0;
  double  val1_data_mumu_tar0=0;
  double  err2_data_mumu_tar0=0;
  double  val2_data_mumu_tar0=0;
  double  err3_data_mumu_tar0=0;
  double  val3_data_mumu_tar0=0;
  double  err4_data_mumu_tar0=0;
  double  val4_data_mumu_tar0=0;
  double  err5_data_mumu_tar0=0;
  double  val5_data_mumu_tar0=0;
  double  err6_data_mumu_tar0=0;
  double  val6_data_mumu_tar0=0;
  double  err7_data_mumu_tar0=0;
  double  val7_data_mumu_tar0=0;
  double  err8_data_mumu_tar0=0;
  double  val8_data_mumu_tar0=0;

  range(h_met_mumu_tar0_data,val1_data_mumu_tar0,err1_data_mumu_tar0,0,20);
  range(h_met_mumu_tar0_data,val2_data_mumu_tar0,err2_data_mumu_tar0,20,40);
  range(h_met_mumu_tar0_data,val3_data_mumu_tar0,err3_data_mumu_tar0,40,60);
  range(h_met_mumu_tar0_data,val4_data_mumu_tar0,err4_data_mumu_tar0,60,80);
  range(h_met_mumu_tar0_data,val5_data_mumu_tar0,err5_data_mumu_tar0,80,100);
  range(h_met_mumu_tar0_data,val6_data_mumu_tar0,err6_data_mumu_tar0,100,150);
  range(h_met_mumu_tar0_data,val7_data_mumu_tar0,err7_data_mumu_tar0,150,200);
  range(h_met_mumu_tar0_data,val8_data_mumu_tar0,err8_data_mumu_tar0,200,-1);

  double err1_ratio_mumu_tar0=0;
  double err2_ratio_mumu_tar0=0;
  double err3_ratio_mumu_tar0=0;
  double err4_ratio_mumu_tar0=0;
  double err5_ratio_mumu_tar0=0;
  double err6_ratio_mumu_tar0=0;
  double err7_ratio_mumu_tar0=0;
  double err8_ratio_mumu_tar0=0;

  err1_ratio_mumu_tar0=err_ratio(val1_data_mumu_tar0, err1_data_mumu_tar0, val1_metsum_mumu_tar0, err1_metsum_mumu_tar0);
  err2_ratio_mumu_tar0=err_ratio(val2_data_mumu_tar0, err2_data_mumu_tar0, val2_metsum_mumu_tar0, err2_metsum_mumu_tar0);
  err3_ratio_mumu_tar0=err_ratio(val3_data_mumu_tar0, err3_data_mumu_tar0, val3_metsum_mumu_tar0, err3_metsum_mumu_tar0);
  err4_ratio_mumu_tar0=err_ratio(val4_data_mumu_tar0, err4_data_mumu_tar0, val4_metsum_mumu_tar0, err4_metsum_mumu_tar0);
  err5_ratio_mumu_tar0=err_ratio(val5_data_mumu_tar0, err5_data_mumu_tar0, val5_metsum_mumu_tar0, err5_metsum_mumu_tar0);
  err6_ratio_mumu_tar0=err_ratio(val6_data_mumu_tar0, err6_data_mumu_tar0, val6_metsum_mumu_tar0, err6_metsum_mumu_tar0);
  err7_ratio_mumu_tar0=err_ratio(val7_data_mumu_tar0, err7_data_mumu_tar0, val7_metsum_mumu_tar0, err7_metsum_mumu_tar0);
  err8_ratio_mumu_tar0=err_ratio(val8_data_mumu_tar0, err8_data_mumu_tar0, val8_metsum_mumu_tar0, err8_metsum_mumu_tar0);
  //----------------------------------------------------------------------------------------------//
  */
  //------------------------------------------met ee tar2--------------------------------------------------------------//
  h_met_ee_tar2_zjets->Add(h_met_emu_tar2_zjets,- Ree);
  h_met_ee_tar2_ttbar->Add(h_met_emu_tar2_ttbar,- Ree);
  h_met_ee_tar2_ww->Add(h_met_emu_tar2_ww,- Ree);
  h_met_ee_tar2_wz->Add(h_met_emu_tar2_wz,- Ree);
  h_met_ee_tar2_zz->Add(h_met_emu_tar2_zz,- Ree);
  h_met_ee_tar2_singlet->Add(h_met_emu_tar2_singlet,- Ree);
  h_met_ee_tar2_ttv->Add(h_met_emu_tar2_ttv,- Ree);
  h_met_ee_tar2_vvv->Add(h_met_emu_tar2_vvv,- Ree);
  h_met_ee_tar2_data->Add(h_met_emu_tar2_data,- Ree);

  h_met_mumu_tar2_zjets->Add(h_met_emu_tar2_zjets,- Rmumu);
  h_met_mumu_tar2_ttbar->Add(h_met_emu_tar2_ttbar,- Rmumu);
  h_met_mumu_tar2_ww->Add(h_met_emu_tar2_ww,- Rmumu);
  h_met_mumu_tar2_wz->Add(h_met_emu_tar2_wz,- Rmumu);
  h_met_mumu_tar2_zz->Add(h_met_emu_tar2_zz,- Rmumu);
  h_met_mumu_tar2_singlet->Add(h_met_emu_tar2_singlet,- Rmumu);
  h_met_mumu_tar2_ttv->Add(h_met_emu_tar2_ttv,- Rmumu);
  h_met_mumu_tar2_vvv->Add(h_met_emu_tar2_vvv,- Rmumu);
  h_met_mumu_tar2_data->Add(h_met_emu_tar2_data,- Rmumu);
  //merge ee and mumu
  h_met_ee_tar2_zjets->Add(h_met_mumu_tar2_zjets);
  h_met_ee_tar2_ttbar->Add(h_met_mumu_tar2_ttbar);
  h_met_ee_tar2_ww->Add(h_met_mumu_tar2_ww);
  h_met_ee_tar2_wz->Add(h_met_mumu_tar2_wz);
  h_met_ee_tar2_zz->Add(h_met_mumu_tar2_zz);
  h_met_ee_tar2_singlet->Add(h_met_mumu_tar2_singlet);
  h_met_ee_tar2_ttv->Add(h_met_mumu_tar2_ttv);
  h_met_ee_tar2_vvv->Add(h_met_mumu_tar2_vvv);
  h_met_ee_tar2_data->Add(h_met_mumu_tar2_data);
  
  TH1F *h_metsum_ee_tar2 = new TH1F("h_metsum_ee_tar2","Met ee tar2",350,0,350);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_zjets);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_ttbar);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_ww);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_wz);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_zz);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_singlet);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_ttv);
  h_metsum_ee_tar2->Add(h_met_ee_tar2_vvv);
  
  double  err1_metsum_ee_tar2=0;
  double  val1_metsum_ee_tar2=0;
  double  err2_metsum_ee_tar2=0;
  double  val2_metsum_ee_tar2=0;
  double  err3_metsum_ee_tar2=0;
  double  val3_metsum_ee_tar2=0;
  double  err4_metsum_ee_tar2=0;
  double  val4_metsum_ee_tar2=0;
  double  err5_metsum_ee_tar2=0;
  double  val5_metsum_ee_tar2=0;
  double  err6_metsum_ee_tar2=0;
  double  val6_metsum_ee_tar2=0;
  double  err7_metsum_ee_tar2=0;
  double  val7_metsum_ee_tar2=0;
  double  err8_metsum_ee_tar2=0;
  double  val8_metsum_ee_tar2=0;
 
  range(h_metsum_ee_tar2,val1_metsum_ee_tar2,err1_metsum_ee_tar2,0,20);
  range(h_metsum_ee_tar2,val2_metsum_ee_tar2,err2_metsum_ee_tar2,20,40);
  range(h_metsum_ee_tar2,val3_metsum_ee_tar2,err3_metsum_ee_tar2,40,60);
  range(h_metsum_ee_tar2,val4_metsum_ee_tar2,err4_metsum_ee_tar2,60,80);
  range(h_metsum_ee_tar2,val5_metsum_ee_tar2,err5_metsum_ee_tar2,80,100);
  range(h_metsum_ee_tar2,val6_metsum_ee_tar2,err6_metsum_ee_tar2,100,150);
  range(h_metsum_ee_tar2,val7_metsum_ee_tar2,err7_metsum_ee_tar2,150,200);
  range(h_metsum_ee_tar2,val8_metsum_ee_tar2,err8_metsum_ee_tar2,200,-1);

  double  err1_data_ee_tar2=0;
  double  val1_data_ee_tar2=0;
  double  err2_data_ee_tar2=0;
  double  val2_data_ee_tar2=0;
  double  err3_data_ee_tar2=0;
  double  val3_data_ee_tar2=0;
  double  err4_data_ee_tar2=0;
  double  val4_data_ee_tar2=0;
  double  err5_data_ee_tar2=0;
  double  val5_data_ee_tar2=0;
  double  err6_data_ee_tar2=0;
  double  val6_data_ee_tar2=0;
  double  err7_data_ee_tar2=0;
  double  val7_data_ee_tar2=0;
  double  err8_data_ee_tar2=0;
  double  val8_data_ee_tar2=0;

  range(h_met_ee_tar2_data,val1_data_ee_tar2,err1_data_ee_tar2,0,20);
  range(h_met_ee_tar2_data,val2_data_ee_tar2,err2_data_ee_tar2,20,40);
  range(h_met_ee_tar2_data,val3_data_ee_tar2,err3_data_ee_tar2,40,60);
  range(h_met_ee_tar2_data,val4_data_ee_tar2,err4_data_ee_tar2,60,80);
  range(h_met_ee_tar2_data,val5_data_ee_tar2,err5_data_ee_tar2,80,100);
  range(h_met_ee_tar2_data,val6_data_ee_tar2,err6_data_ee_tar2,100,150);
  range(h_met_ee_tar2_data,val7_data_ee_tar2,err7_data_ee_tar2,150,200);
  range(h_met_ee_tar2_data,val8_data_ee_tar2,err8_data_ee_tar2,200,-1);

  double err1_ratio_ee_tar2=0;
  double err2_ratio_ee_tar2=0;
  double err3_ratio_ee_tar2=0;
  double err4_ratio_ee_tar2=0;
  double err5_ratio_ee_tar2=0;
  double err6_ratio_ee_tar2=0;
  double err7_ratio_ee_tar2=0;
  double err8_ratio_ee_tar2=0;

  err1_ratio_ee_tar2=err_ratio(val1_data_ee_tar2, err1_data_ee_tar2, val1_metsum_ee_tar2, err1_metsum_ee_tar2);
  err2_ratio_ee_tar2=err_ratio(val2_data_ee_tar2, err2_data_ee_tar2, val2_metsum_ee_tar2, err2_metsum_ee_tar2);
  err3_ratio_ee_tar2=err_ratio(val3_data_ee_tar2, err3_data_ee_tar2, val3_metsum_ee_tar2, err3_metsum_ee_tar2);
  err4_ratio_ee_tar2=err_ratio(val4_data_ee_tar2, err4_data_ee_tar2, val4_metsum_ee_tar2, err4_metsum_ee_tar2);
  err5_ratio_ee_tar2=err_ratio(val5_data_ee_tar2, err5_data_ee_tar2, val5_metsum_ee_tar2, err5_metsum_ee_tar2);
  err6_ratio_ee_tar2=err_ratio(val6_data_ee_tar2, err6_data_ee_tar2, val6_metsum_ee_tar2, err6_metsum_ee_tar2);
  err7_ratio_ee_tar2=err_ratio(val7_data_ee_tar2, err7_data_ee_tar2, val7_metsum_ee_tar2, err7_metsum_ee_tar2);
  err8_ratio_ee_tar2=err_ratio(val8_data_ee_tar2, err8_data_ee_tar2, val8_metsum_ee_tar2, err8_metsum_ee_tar2);
  //----------------------------------------------------------------------------------------------//
  /* 
  //--------------------------------------met mumu tar2--------------------------------------------//
  h_met_mumu_tar2_zjets->Add(h_met_emu_tar2_zjets,- Rmumu);
  h_met_mumu_tar2_ttbar->Add(h_met_emu_tar2_ttbar,- Rmumu);
  h_met_mumu_tar2_ww->Add(h_met_emu_tar2_ww,- Rmumu);
  h_met_mumu_tar2_wz->Add(h_met_emu_tar2_wz,- Rmumu);
  h_met_mumu_tar2_zz->Add(h_met_emu_tar2_zz,- Rmumu);
  h_met_mumu_tar2_singlet->Add(h_met_emu_tar2_singlet,- Rmumu);
  h_met_mumu_tar2_ttv->Add(h_met_emu_tar2_ttv,- Rmumu);
  h_met_mumu_tar2_vvv->Add(h_met_emu_tar2_vvv,- Rmumu);
  h_met_mumu_tar2_data->Add(h_met_emu_tar2_data,- Rmumu);
  
  TH1F *h_metsum_mumu_tar2 = new TH1F("h_metsum_mumu_tar2","Met mumu tar2",350,0,350);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_zjets);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_ttbar);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_ww);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_wz);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_zz);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_singlet);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_ttv);
  h_metsum_mumu_tar2->Add(h_met_mumu_tar2_vvv);
  
  double  err1_metsum_mumu_tar2=0;
  double  val1_metsum_mumu_tar2=0;
  double  err2_metsum_mumu_tar2=0;
  double  val2_metsum_mumu_tar2=0;
  double  err3_metsum_mumu_tar2=0;
  double  val3_metsum_mumu_tar2=0;
  double  err4_metsum_mumu_tar2=0;
  double  val4_metsum_mumu_tar2=0;
  double  err5_metsum_mumu_tar2=0;
  double  val5_metsum_mumu_tar2=0;
  double  err6_metsum_mumu_tar2=0;
  double  val6_metsum_mumu_tar2=0;
  double  err7_metsum_mumu_tar2=0;
  double  val7_metsum_mumu_tar2=0;
  double  err8_metsum_mumu_tar2=0;
  double  val8_metsum_mumu_tar2=0;
 
  range(h_metsum_mumu_tar2,val1_metsum_mumu_tar2,err1_metsum_mumu_tar2,0,20);
  range(h_metsum_mumu_tar2,val2_metsum_mumu_tar2,err2_metsum_mumu_tar2,20,40);
  range(h_metsum_mumu_tar2,val3_metsum_mumu_tar2,err3_metsum_mumu_tar2,40,60);
  range(h_metsum_mumu_tar2,val4_metsum_mumu_tar2,err4_metsum_mumu_tar2,60,80);
  range(h_metsum_mumu_tar2,val5_metsum_mumu_tar2,err5_metsum_mumu_tar2,80,100);
  range(h_metsum_mumu_tar2,val6_metsum_mumu_tar2,err6_metsum_mumu_tar2,100,150);
  range(h_metsum_mumu_tar2,val7_metsum_mumu_tar2,err7_metsum_mumu_tar2,150,200);
  range(h_metsum_mumu_tar2,val8_metsum_mumu_tar2,err8_metsum_mumu_tar2,200,-1);

  double  err1_data_mumu_tar2=0;
  double  val1_data_mumu_tar2=0;
  double  err2_data_mumu_tar2=0;
  double  val2_data_mumu_tar2=0;
  double  err3_data_mumu_tar2=0;
  double  val3_data_mumu_tar2=0;
  double  err4_data_mumu_tar2=0;
  double  val4_data_mumu_tar2=0;
  double  err5_data_mumu_tar2=0;
  double  val5_data_mumu_tar2=0;
  double  err6_data_mumu_tar2=0;
  double  val6_data_mumu_tar2=0;
  double  err7_data_mumu_tar2=0;
  double  val7_data_mumu_tar2=0;
  double  err8_data_mumu_tar2=0;
  double  val8_data_mumu_tar2=0;

  range(h_met_mumu_tar2_data,val1_data_mumu_tar2,err1_data_mumu_tar2,0,20);
  range(h_met_mumu_tar2_data,val2_data_mumu_tar2,err2_data_mumu_tar2,20,40);
  range(h_met_mumu_tar2_data,val3_data_mumu_tar2,err3_data_mumu_tar2,40,60);
  range(h_met_mumu_tar2_data,val4_data_mumu_tar2,err4_data_mumu_tar2,60,80);
  range(h_met_mumu_tar2_data,val5_data_mumu_tar2,err5_data_mumu_tar2,80,100);
  range(h_met_mumu_tar2_data,val6_data_mumu_tar2,err6_data_mumu_tar2,100,150);
  range(h_met_mumu_tar2_data,val7_data_mumu_tar2,err7_data_mumu_tar2,150,200);
  range(h_met_mumu_tar2_data,val8_data_mumu_tar2,err8_data_mumu_tar2,200,-1);

  double err1_ratio_mumu_tar2=0;
  double err2_ratio_mumu_tar2=0;
  double err3_ratio_mumu_tar2=0;
  double err4_ratio_mumu_tar2=0;
  double err5_ratio_mumu_tar2=0;
  double err6_ratio_mumu_tar2=0;
  double err7_ratio_mumu_tar2=0;
  double err8_ratio_mumu_tar2=0;

  err1_ratio_mumu_tar2=err_ratio(val1_data_mumu_tar2, err1_data_mumu_tar2, val1_metsum_mumu_tar2, err1_metsum_mumu_tar2);
  err2_ratio_mumu_tar2=err_ratio(val2_data_mumu_tar2, err2_data_mumu_tar2, val2_metsum_mumu_tar2, err2_metsum_mumu_tar2);
  err3_ratio_mumu_tar2=err_ratio(val3_data_mumu_tar2, err3_data_mumu_tar2, val3_metsum_mumu_tar2, err3_metsum_mumu_tar2);
  err4_ratio_mumu_tar2=err_ratio(val4_data_mumu_tar2, err4_data_mumu_tar2, val4_metsum_mumu_tar2, err4_metsum_mumu_tar2);
  err5_ratio_mumu_tar2=err_ratio(val5_data_mumu_tar2, err5_data_mumu_tar2, val5_metsum_mumu_tar2, err5_metsum_mumu_tar2);
  err6_ratio_mumu_tar2=err_ratio(val6_data_mumu_tar2, err6_data_mumu_tar2, val6_metsum_mumu_tar2, err6_metsum_mumu_tar2);
  err7_ratio_mumu_tar2=err_ratio(val7_data_mumu_tar2, err7_data_mumu_tar2, val7_metsum_mumu_tar2, err7_metsum_mumu_tar2);
  err8_ratio_mumu_tar2=err_ratio(val8_data_mumu_tar2, err8_data_mumu_tar2, val8_metsum_mumu_tar2, err8_metsum_mumu_tar2);
  //----------------------------------------------------------------------------------------------//
  */

  //--------------------------------------------------------------------------------------------------------------------------------------//
  //----------------------------------------------Met ee+mumu inc Table-------------------------------------------------------------------//
  //--------------------------------------------------------------------------------------------------------------------------------------//
  cout<<"\n\n"<<setw(70)<<"met ee+mumu inc"<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"0-20 GeV"<<setw(width)<<"20-40 GeV"<<setw(width)<<"40-60 GeV"<<setw(width)<<"60-80 GeV"<<setw(width)<<"80-100 GeV"<<endl;
  //MC
  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val1_metsum_ee_inc, err1_metsum_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val2_metsum_ee_inc, err2_metsum_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val3_metsum_ee_inc, err3_metsum_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val4_metsum_ee_inc, err4_metsum_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val5_metsum_ee_inc, err5_metsum_ee_inc)<<endl;
  //data
  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val1_data_ee_inc, err1_data_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val2_data_ee_inc, err2_data_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val3_data_ee_inc, err3_data_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val4_data_ee_inc, err4_data_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val5_data_ee_inc, err5_data_ee_inc)<<endl;
 
   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val1_data_ee_inc/val1_metsum_ee_inc,err1_ratio_ee_inc)<<setw(width)<<Form("%.2f +/- %.2f", val2_data_ee_inc/val2_metsum_ee_inc, err2_ratio_ee_inc)<<setw(width)<<Form("%.2f +/- %.2f", val3_data_ee_inc/val3_metsum_ee_inc,err3_ratio_ee_inc)<<setw(width)<<Form("%.2f +/- %.2f",  val4_data_ee_inc/val4_metsum_ee_inc,err4_ratio_ee_inc)<<setw(width)<<Form("%.2f +/- %.2f",  val5_data_ee_inc/val5_metsum_ee_inc,err5_ratio_ee_inc)<<endl;
  

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"100-150 GeV"<<setw(width)<<"150-200 GeV"<<setw(width)<<"> 200 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val6_metsum_ee_inc, err6_metsum_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val7_metsum_ee_inc, err7_metsum_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val8_metsum_ee_inc, err8_metsum_ee_inc)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val6_data_ee_inc, err6_data_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val7_data_ee_inc, err7_data_ee_inc)<<setw(width)<<Form("%.1f +/- %.1f", val8_data_ee_inc, err8_data_ee_inc)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val6_data_ee_inc/val6_metsum_ee_inc,err6_ratio_ee_inc)<<setw(width)<<Form("%.2f +/- %.2f", val7_data_ee_inc/val7_metsum_ee_inc, err7_ratio_ee_inc)<<setw(width)<<Form("%.2f +/- %.2f", val8_data_ee_inc/val8_metsum_ee_inc,err8_ratio_ee_inc)<<endl;
  //----------------------------------------------------------------------------------------------------------------------------------------//
  /*
 //----------------------------------------------Met mumu inc Table------------------------------------------------------------------------//
  cout<<"\n\n"<<setw(70)<<"met mumu inc"<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"0-20 GeV"<<setw(width)<<"20-40 GeV"<<setw(width)<<"40-60 GeV"<<setw(width)<<"60-80 GeV"<<setw(width)<<"80-100 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val1_metsum_mumu_inc, err1_metsum_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val2_metsum_mumu_inc, err2_metsum_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val3_metsum_mumu_inc, err3_metsum_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val4_metsum_mumu_inc, err4_metsum_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val5_metsum_mumu_inc, err5_metsum_mumu_inc)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val1_data_mumu_inc, err1_data_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val2_data_mumu_inc, err2_data_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val3_data_mumu_inc, err3_data_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val4_data_mumu_inc, err4_data_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val5_data_mumu_inc, err5_data_mumu_inc)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val1_data_mumu_inc/val1_metsum_mumu_inc,err1_ratio_mumu_inc)<<setw(width)<<Form("%.2f +/- %.2f", val2_data_mumu_inc/val2_metsum_mumu_inc, err2_ratio_mumu_inc)<<setw(width)<<Form("%.2f +/- %.2f", val3_data_mumu_inc/val3_metsum_mumu_inc,err3_ratio_mumu_inc)<<setw(width)<<Form("%.2f +/- %.2f",  val4_data_mumu_inc/val4_metsum_mumu_inc,err4_ratio_mumu_inc)<<setw(width)<<Form("%.2f +/- %.2f",  val5_data_mumu_inc/val5_metsum_mumu_inc,err5_ratio_mumu_inc)<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"100-150 GeV"<<setw(width)<<"150-200 GeV"<<setw(width)<<"> 200 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val6_metsum_mumu_inc, err6_metsum_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val7_metsum_mumu_inc, err7_metsum_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val8_metsum_mumu_inc, err8_metsum_mumu_inc)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val6_data_mumu_inc, err6_data_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val7_data_mumu_inc, err7_data_mumu_inc)<<setw(width)<<Form("%.1f +/- %.1f", val8_data_mumu_inc, err8_data_mumu_inc)<<endl; 

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val6_data_mumu_inc/val6_metsum_mumu_inc,err6_ratio_mumu_inc)<<setw(width)<<Form("%.2f +/- %.2f", val7_data_mumu_inc/val7_metsum_mumu_inc, err7_ratio_mumu_inc)<<setw(width)<<Form("%.2f +/- %.2f", val8_data_mumu_inc/val8_metsum_mumu_inc,err8_ratio_mumu_inc)<<endl;
  //----------------------------------------------------------------------------------------------------------------------------------------//
  */
  //----------------------------------------------Met ee+mumu tar0 Table--------------------------------------------------------------------//

  cout<<"\n\n"<<setw(70)<<"met ee+mumu tar0"<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"0-20 GeV"<<setw(width)<<"20-40 GeV"<<setw(width)<<"40-60 GeV"<<setw(width)<<"60-80 GeV"<<setw(width)<<"80-100 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val1_metsum_ee_tar0, err1_metsum_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val2_metsum_ee_tar0, err2_metsum_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val3_metsum_ee_tar0, err3_metsum_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val4_metsum_ee_tar0, err4_metsum_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val5_metsum_ee_tar0, err5_metsum_ee_tar0)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val1_data_ee_tar0, err1_data_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val2_data_ee_tar0, err2_data_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val3_data_ee_tar0, err3_data_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val4_data_ee_tar0, err4_data_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val5_data_ee_tar0, err5_data_ee_tar0)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val1_data_ee_tar0/val1_metsum_ee_tar0,err1_ratio_ee_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val2_data_ee_tar0/val2_metsum_ee_tar0, err2_ratio_ee_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val3_data_ee_tar0/val3_metsum_ee_tar0,err3_ratio_ee_tar0)<<setw(width)<<Form("%.2f +/- %.2f",  val4_data_ee_tar0/val4_metsum_ee_tar0,err4_ratio_ee_tar0)<<setw(width)<<Form("%.2f +/- %.2f",  val5_data_ee_tar0/val5_metsum_ee_tar0,err5_ratio_ee_tar0)<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"100-150 GeV"<<setw(width)<<"150-200 GeV"<<setw(width)<<"> 200 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val6_metsum_ee_tar0, err6_metsum_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val7_metsum_ee_tar0, err7_metsum_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val8_metsum_ee_tar0, err8_metsum_ee_tar0)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val6_data_ee_tar0, err6_data_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val7_data_ee_tar0, err7_data_ee_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val8_data_ee_tar0, err8_data_ee_tar0)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val6_data_ee_tar0/val6_metsum_ee_tar0,err6_ratio_ee_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val7_data_ee_tar0/val7_metsum_ee_tar0, err7_ratio_ee_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val8_data_ee_tar0/val8_metsum_ee_tar0,err8_ratio_ee_tar0)<<endl;
 
  //----------------------------------------------------------------------------------------------------------------------------------------//
  /*
  //---------------------------------------------Met mumu tar0 Table------------------------------------------------------------------------//
  cout<<"\n\n"<<setw(70)<<"met mumu tar0"<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"0-20 GeV"<<setw(width)<<"20-40 GeV"<<setw(width)<<"40-60 GeV"<<setw(width)<<"60-80 GeV"<<setw(width)<<"80-100 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val1_metsum_mumu_tar0, err1_metsum_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val2_metsum_mumu_tar0, err2_metsum_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val3_metsum_mumu_tar0, err3_metsum_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val4_metsum_mumu_tar0, err4_metsum_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val5_metsum_mumu_tar0, err5_metsum_mumu_tar0)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val1_data_mumu_tar0, err1_data_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val2_data_mumu_tar0, err2_data_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val3_data_mumu_tar0, err3_data_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val4_data_mumu_tar0, err4_data_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val5_data_mumu_tar0, err5_data_mumu_tar0)<<endl;
 
   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val1_data_mumu_tar0/val1_metsum_mumu_tar0,err1_ratio_mumu_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val2_data_mumu_tar0/val2_metsum_mumu_tar0, err2_ratio_mumu_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val3_data_mumu_tar0/val3_metsum_mumu_tar0,err3_ratio_mumu_tar0)<<setw(width)<<Form("%.2f +/- %.2f",  val4_data_mumu_tar0/val4_metsum_mumu_tar0,err4_ratio_mumu_tar0)<<setw(width)<<Form("%.2f +/- %.2f",  val5_data_mumu_tar0/val5_metsum_mumu_tar0,err5_ratio_mumu_tar0)<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"100-150 GeV"<<setw(width)<<"150-200 GeV"<<setw(width)<<"> 200 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val6_metsum_mumu_tar0, err6_metsum_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val7_metsum_mumu_tar0, err7_metsum_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val8_metsum_mumu_tar0, err8_metsum_mumu_tar0)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val6_data_mumu_tar0, err6_data_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val7_data_mumu_tar0, err7_data_mumu_tar0)<<setw(width)<<Form("%.1f +/- %.1f", val8_data_mumu_tar0, err8_data_mumu_tar0)<<endl; 

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val6_data_mumu_tar0/val6_metsum_mumu_tar0,err6_ratio_mumu_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val7_data_mumu_tar0/val7_metsum_mumu_tar0, err7_ratio_mumu_tar0)<<setw(width)<<Form("%.2f +/- %.2f", val8_data_mumu_tar0/val8_metsum_mumu_tar0,err8_ratio_mumu_tar0)<<endl;
  //----------------------------------------------------------------------------------------------------------------------------------------//
  */
  //----------------------------------------------Met ee+mumu tar2 Table------------------------------------------------------------------------//

  cout<<"\n\n"<<setw(70)<<"met ee+mumu tar2"<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"0-20 GeV"<<setw(width)<<"20-40 GeV"<<setw(width)<<"40-60 GeV"<<setw(width)<<"60-80 GeV"<<setw(width)<<"80-100 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val1_metsum_ee_tar2, err1_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val2_metsum_ee_tar2, err2_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val3_metsum_ee_tar2, err3_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val4_metsum_ee_tar2, err4_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val5_metsum_ee_tar2, err5_metsum_ee_tar2)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val1_data_ee_tar2, err1_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val2_data_ee_tar2, err2_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val3_data_ee_tar2, err3_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val4_data_ee_tar2, err4_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val5_data_ee_tar2, err5_data_ee_tar2)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val1_data_ee_tar2/val1_metsum_ee_tar2,err1_ratio_ee_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val2_data_ee_tar2/val2_metsum_ee_tar2, err2_ratio_ee_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val3_data_ee_tar2/val3_metsum_ee_tar2,err3_ratio_ee_tar2)<<setw(width)<<Form("%.2f +/- %.2f",  val4_data_ee_tar2/val4_metsum_ee_tar2,err4_ratio_ee_tar2)<<setw(width)<<Form("%.2f +/- %.2f",  val5_data_ee_tar2/val5_metsum_ee_tar2,err5_ratio_ee_tar2)<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"100-150 GeV"<<setw(width)<<"150-200 GeV"<<setw(width)<<"> 200 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val6_metsum_ee_tar2, err6_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val7_metsum_ee_tar2, err7_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val8_metsum_ee_tar2, err8_metsum_ee_tar2)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val6_data_ee_tar2, err6_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val7_data_ee_tar2, err7_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val8_data_ee_tar2, err8_data_ee_tar2)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val6_data_ee_tar2/val6_metsum_ee_tar2,err6_ratio_ee_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val7_data_ee_tar2/val7_metsum_ee_tar2, err7_ratio_ee_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val8_data_ee_tar2/val8_metsum_ee_tar2,err8_ratio_ee_tar2)<<endl;
 
  //----------------------------------------------------------------------------------------------------------------------------------------//
  /*
 //----------------------------------------------Met mumu tar2 Table------------------------------------------------------------------------//
  cout<<"\n\n"<<setw(70)<<"met mumu tar2"<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"0-20 GeV"<<setw(width)<<"20-40 GeV"<<setw(width)<<"40-60 GeV"<<setw(width)<<"60-80 GeV"<<setw(width)<<"80-100 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val1_metsum_mumu_tar2, err1_metsum_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val2_metsum_mumu_tar2, err2_metsum_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val3_metsum_mumu_tar2, err3_metsum_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val4_metsum_mumu_tar2, err4_metsum_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val5_metsum_mumu_tar2, err5_metsum_mumu_tar2)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val1_data_mumu_tar2, err1_data_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val2_data_mumu_tar2, err2_data_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val3_data_mumu_tar2, err3_data_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val4_data_mumu_tar2, err4_data_mumu_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val5_data_mumu_tar2, err5_data_mumu_tar2)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val1_data_mumu_tar2/val1_metsum_mumu_tar2,err1_ratio_mumu_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val2_data_mumu_tar2/val2_metsum_mumu_tar2, err2_ratio_mumu_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val3_data_mumu_tar2/val3_metsum_mumu_tar2,err3_ratio_mumu_tar2)<<setw(width)<<Form("%.2f +/- %.2f",  val4_data_mumu_tar2/val4_metsum_mumu_tar2,err4_ratio_mumu_tar2)<<setw(width)<<Form("%.2f +/- %.2f",  val5_data_mumu_tar2/val5_metsum_mumu_tar2,err5_ratio_mumu_tar2)<<endl;

  cout<<"\n"<<setw(width)<<" "<<setw(width)<<"100-150 GeV"<<setw(width)<<"150-200 GeV"<<setw(width)<<"> 200 GeV"<<endl;

  cout<<"\n"<<setw(width)<<"Total SM MC"<<setw(width)<<Form("%.1f +/- %.1f", val6_metsum_ee_tar2, err6_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val7_metsum_ee_tar2, err7_metsum_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val8_metsum_ee_tar2, err8_metsum_ee_tar2)<<endl;

  cout<<"\n"<<setw(width)<<"Data-FSBG"<<setw(width)<<Form("%.1f +/- %.1f", val6_data_ee_tar2, err6_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val7_data_ee_tar2, err7_data_ee_tar2)<<setw(width)<<Form("%.1f +/- %.1f", val8_data_ee_tar2, err8_data_ee_tar2)<<endl;

   //ratio
  cout<<"\n"<<setw(width)<<"Data/MC"<<setw(width)<<Form("%.2f +/- %.2f", val6_data_mumu_tar2/val6_metsum_mumu_tar2,err6_ratio_mumu_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val7_data_mumu_tar2/val7_metsum_mumu_tar2, err7_ratio_mumu_tar2)<<setw(width)<<Form("%.2f +/- %.2f", val8_data_mumu_tar2/val8_metsum_mumu_tar2,err8_ratio_mumu_tar2)<<endl;

  //----------------------------------------------------------------------------------------------------------------------------------------//
  */

  return 0;
}
