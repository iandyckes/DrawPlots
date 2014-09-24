// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TCanvas.h"
#include "THStack.h"  
#include "TMap.h"
#include "TLorentzVector.h"
#include "TMath.h"

// ZMET
#include "ZMET.h"

// goodrun (for json)
#include "goodrun.cc"

using namespace std;
using namespace ZMet;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {


  set_goodrun_file("Merged_190456-208686_8TeV_PromptReReco_Collisions12_goodruns.txt");

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  //add hists
  TH1F *h_mll_ee_inc = new TH1F("h_mll_ee_inc", "ee Dilepton Mass", 200,0,200);//200,0,2000
  h_mll_ee_inc->SetDirectory(rootdir);
  h_mll_ee_inc->Sumw2();
  
  TH1F* h_mll_mumu_inc = new TH1F("h_mll_mumu_inc","MuMu Dilepton Mass",200,0,200);
  h_mll_mumu_inc->SetDirectory(rootdir);
  h_mll_mumu_inc->Sumw2();

  TH1F* h_mll_ee_tar = new TH1F("h_mll_ee_tar","ee Dilepton Mass (Targeted)",200,0,200);
  h_mll_ee_tar->SetDirectory(rootdir);
  h_mll_ee_tar->Sumw2();

  TH1F* h_mll_mumu_tar = new TH1F("h_mll_mumu_tar","MuMu Dilepton Mass (Targeted)",200,0,200);
  h_mll_mumu_tar->SetDirectory(rootdir);
  h_mll_mumu_tar->Sumw2();

  TH1F* h_met_ee_inc = new TH1F("h_met_ee_inc","ee MET (Inclusive)",350,0,350);//350,0,350 for all met histos
  h_met_ee_inc->SetDirectory(rootdir);
  h_met_ee_inc->Sumw2();

  TH1F* h_met_mumu_inc = new TH1F("h_met_mumu_inc","mumu MET (Inclusive)",350,0,350);
  h_met_mumu_inc->SetDirectory(rootdir);
  h_met_mumu_inc->Sumw2();

  TH2F* h_met_ll_inc = new TH2F("h_met_ll_inc","ee+mumu MET (Inclusive)",700,-350,350,700,-350,350);//positive and negative met.  2D!
  h_met_ll_inc->SetDirectory(rootdir);
  h_met_ll_inc->Sumw2();

  TH2F* h_met_emu2_inc = new TH2F("h_met_emu2_inc","emu MET (Inclusive)",700,-350,350,700,-350,350);//positive and negative met.  2D!
  h_met_emu2_inc->SetDirectory(rootdir);
  h_met_emu2_inc->Sumw2();

  TH1F* h_met_emu_inc = new TH1F("h_met_emu_inc","emu MET (Inclusive)",350,0,350);
  h_met_emu_inc->SetDirectory(rootdir);
  h_met_emu_inc->Sumw2();

  TH1F* h_met_ee_tar_njets0 = new TH1F("h_met_ee_tar_njets0","ee MET (Targeted, no jets)",350,0,350);
  h_met_ee_tar_njets0->SetDirectory(rootdir);
  h_met_ee_tar_njets0->Sumw2();

  TH1F* h_met_mumu_tar_njets0 = new TH1F("h_met_mumu_tar_njets0","mumu MET (Targeted, no jets)",350,0,350);
  h_met_mumu_tar_njets0->SetDirectory(rootdir);
  h_met_mumu_tar_njets0->Sumw2();

  TH2F* h_met_ll_tar_njets0 = new TH2F("h_met_ll_tar_njets0","ee+mumu MET (Targeted, no jets)",700,-350,350,700,-350,350);//positive and negative met.  2D!
  h_met_ll_tar_njets0->SetDirectory(rootdir);
  h_met_ll_tar_njets0->Sumw2();

  TH2F* h_met_emu2_tar_njets0 = new TH2F("h_met_emu2_tar_njets0","emu MET (Targeted, no jets)",700,-350,350,700,-350,350);//positive and negative met.  2D!
  h_met_emu2_tar_njets0->SetDirectory(rootdir);
  h_met_emu2_tar_njets0->Sumw2();

  TH1F* h_met_emu_tar_njets0 = new TH1F("h_met_emu_tar_njets0","emu MET (Targeted, no jets)",350,0,350);
  h_met_emu_tar_njets0->SetDirectory(rootdir);
  h_met_emu_tar_njets0->Sumw2();

  TH1F* h_met_ee_tar_njets2 = new TH1F("h_met_ee_tar_njets2","ee MET (Targeted, njets >= 2)",350,0,350);
  h_met_ee_tar_njets2->SetDirectory(rootdir);
  h_met_ee_tar_njets2->Sumw2();

  TH1F* h_met_mumu_tar_njets2 = new TH1F("h_met_mumu_tar_njets2","mumu MET (Targeted, njets >= 2)",350,0,350);
  h_met_mumu_tar_njets2->SetDirectory(rootdir);
  h_met_mumu_tar_njets2->Sumw2();

  TH2F* h_met_ll_tar_njets2 = new TH2F("h_met_ll_tar_njets2","ee+mumu MET (Targeted, njets >= 2)",700,-350,350,700,-350,350);//positive and negative met.  2D!
  h_met_ll_tar_njets2->SetDirectory(rootdir);
  h_met_ll_tar_njets2->Sumw2();

  TH2F* h_met_emu2_tar_njets2 = new TH2F("h_met_emu2_tar_njets2","emu MET (Targeted, njets >= 2)",700,-350,350,700,-350,350);//positive and negative met.  2D!
  h_met_emu2_tar_njets2->SetDirectory(rootdir);
  h_met_emu2_tar_njets2->Sumw2();

  TH1F* h_met_emu_tar_njets2 = new TH1F("h_met_emu_tar_njets2","emu MET (Targeted, njets >= 2)",350,0,350);
  h_met_emu_tar_njets2->SetDirectory(rootdir);
  h_met_emu_tar_njets2->Sumw2();

  TH1F* h_nvtx_scaled= new TH1F("h_nvtx_scaled","Number of Vertices (scaled)",35,0,35);
  h_nvtx_scaled->SetDirectory(rootdir);
  h_nvtx_scaled->Sumw2();

  TH1F* h_nvtx_unscaled= new TH1F("h_nvtx_unscaled","Number of Vertices (unscaled)",35,0,35);
  h_nvtx_unscaled->SetDirectory(rootdir);
  h_nvtx_unscaled->Sumw2();

  TH1F* h_phi= new TH1F("h_phi","phi distribution",40,-TMath::Pi(),TMath::Pi());
  h_phi->SetDirectory(rootdir);
  h_phi->Sumw2();

  TFile *InputFile = new TFile("phiRatio.root","read");
  TH1F *h_phi_ratio = (TH1F*) InputFile->Get("h_phidata_clone_scaled")->Clone("h_phi_ratio");

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("T1");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    zmet.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      zmet.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      ZMET::progress( nEventsTotal, nEventsChain );
	  
	  //------------------------- CUT --------------------------------------------
	 
	  //
	  if( (zmet.leptype() != 0) && (zmet.leptype() != 1) && (zmet.leptype() != 2) ){ continue;} // 0 -> e's,  1 -> Mu's, 2->emu


	  float ptcutoff=20.;
	  if( ( zmet.lep1().Pt() <= ptcutoff) || ( zmet.lep2().Pt() <= ptcutoff) ) {continue;}

	  float etacutoff=2.4;
	  if( ( abs( zmet.lep1().Eta()) >= etacutoff) || ( abs( zmet.lep2().Eta()) >= etacutoff) ) {continue;}
	  
	  //barrel/endcap eta cuts
	  if( ( abs( zmet.lep1().Eta()) > 1.4 && abs( zmet.lep1().Eta()) < 1.6) || ( abs( zmet.lep2().Eta()) > 1.4 && abs( zmet.lep2().Eta()) < 1.6) ) 
		{continue;}

	  double weight_mc=zmet.weight()*19.5*zmet.trgeff()*zmet.vtxweight()*( h_phi_ratio->GetBinContent(h_phi_ratio->FindBin(zmet.pfmetphi())) );
	  double weight_data=1.;
	  //cout<< h_phi_ratio->GetBinContent(h_phi_ratio->FindBin(zmet.pfmetphi()))<<endl;
	  
	  //---------------------------------------
	  if (zmet.lep3().Pt() > 10.){continue;}
	  //---------------------------------------
	  //if (zmet.nvtx() > 7.) {continue;}
	  //---------------------------------------
	  //if (zmet.nvtx() <= 7. || zmet.nvtx() >= 20.) {continue;}
	  //---------------------------------------
	  //if (zmet.nvtx() < 20.) {continue;}

	  //Data cuts
	  if(zmet.isdata() && goodrun( zmet.run(), zmet.lumi() )  )
		{
		  if(!(zmet.csc()==0 && 
			   zmet.hbhe()==1 && 
			   zmet.hcallaser()==1 && 
			   //zmet.ecallaser()==1 && 
			   zmet.ecaltp()==1 && 
			   zmet.trkfail()==1 && 
			   zmet.eebadsc()==1 && 
			   zmet.hbhenew()==1))     
			continue;    //met filter variables and trigger variables
		}

	  if(zmet.isdata() && goodrun( zmet.run(), zmet.lumi() )  )
		{
		  if( zmet.leptype()==0 )
			{
			  if( zmet.ee()==0 )  
				{continue;}
			}
		  else if( zmet.leptype()==1 )
			{
			  if( zmet.mm()==0 )              
				{continue;}
			}
		  else if( zmet.leptype()==2 )
			{
			  if( zmet.em()==0 && zmet.me()==0 )
				{continue;}
			}
		}

	  //Fill  histograms
	  if (zmet.isdata() && goodrun( zmet.run(), zmet.lumi() )  )  //If it is data from good run
		{
		  if (zmet.leptype()==0 || zmet.leptype()==1 || zmet.leptype()==2)
			{
			  h_nvtx_scaled->Fill(zmet.nvtx(),weight_data);
			  h_nvtx_unscaled->Fill(zmet.nvtx(),weight_data);
			  h_phi->Fill(zmet.pfmetphi(), weight_data);
			}
		  if (zmet.leptype() == 0)
			{
			  h_mll_ee_inc->Fill(zmet.dilmass(),weight_data);
			 
			  if(zmet.dilmass() > 81. && zmet.dilmass() < 101. /*&& zmet.lep3().Pt() <= 10.*/ )             //Check this.  Only apply to MET histos
				{
				  h_met_ee_inc->Fill(zmet.pfmet(),weight_data);
				  h_met_ll_inc->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);
				 
				  if(zmet.njets()==0)
					{
					  h_met_ee_tar_njets0->Fill(zmet.pfmet(),weight_data);
					  h_met_ll_tar_njets0->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);
					}
				  else if(zmet.njets()>=2)
					{
					  h_met_ee_tar_njets2->Fill(zmet.pfmet(),weight_data);
					  h_met_ll_tar_njets2->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);
					}
				}

			  if ( zmet.lep3().Pt() <= 10.  &&  zmet.njets()>=2  &&  zmet.nbcsvm()==0 && (zmet.mjj() > 70. && zmet.mjj() < 110.) )
				 {h_mll_ee_tar->Fill(zmet.dilmass(),weight_data);}
			}
		  else if (zmet.leptype()==1)
			{
			  h_mll_mumu_inc->Fill(zmet.dilmass(),weight_data);
			 
			  if(zmet.dilmass() > 81. && zmet.dilmass() < 101. /*&& zmet.lep3().Pt() <= 10.*/)             //check this.
				{		  
				  h_met_mumu_inc->Fill(zmet.pfmet(),weight_data);
				  h_met_ll_inc->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);
			
				  if(zmet.njets()==0)
					{
					  h_met_mumu_tar_njets0->Fill(zmet.pfmet(),weight_data);
					  h_met_ll_tar_njets0->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);					
					}
				  else if(zmet.njets()>=2)
					{
					  h_met_mumu_tar_njets2->Fill(zmet.pfmet(),weight_data);
					  h_met_ll_tar_njets2->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);					
					}

				}

			  if ( zmet.lep3().Pt() <= 10.   &&  zmet.njets()>=2  &&  zmet.nbcsvm()==0 && (zmet.mjj() > 70. && zmet.mjj() < 110.) )  
				 {h_mll_mumu_tar->Fill(zmet.dilmass(),weight_data);}
			}
		  else if (zmet.leptype()==2)
			{
			  //h_mll_mumu_inc->Fill(zmet.dilmass(),weight_data);
			 
			  if(zmet.dilmass() > 81. && zmet.dilmass() < 101. /*&& zmet.lep3().Pt() <= 10.*/)             //check this.
				{		  
				  h_met_emu_inc->Fill(zmet.pfmet(),weight_data);
				  h_met_emu2_inc->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);
			
				  if(zmet.njets()==0)
					{
					  h_met_emu_tar_njets0->Fill(zmet.pfmet(),weight_data);
					  h_met_emu2_tar_njets0->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);
					}
				  else if(zmet.njets()>=2)
					{
					  h_met_emu_tar_njets2->Fill(zmet.pfmet(),weight_data);
					  h_met_emu2_tar_njets2->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_data);
					}
				}

			  // if ( zmet.lep3().Pt() <= 10.   &&  zmet.njets()>=2  &&  zmet.nbcsvm()==0 && (zmet.mjj() > 70. && zmet.mjj() < 110.) )  
			  // {h_mll_mumu_tar->Fill(zmet.dilmass(),weight_data);}
			}
		}


	  else if (!zmet.isdata())   //If it is MC
		{
		  if (zmet.leptype()==0 || zmet.leptype()==1 || zmet.leptype()==2)
			{
			  h_nvtx_scaled->Fill(zmet.nvtx(), zmet.weight()*19.5*zmet.trgeff()*zmet.vtxweight() );
			  h_nvtx_unscaled->Fill(zmet.nvtx(), zmet.weight()*19.5*zmet.trgeff() );
			  h_phi->Fill(zmet.pfmetphi(), weight_mc);			}

		  if (zmet.leptype() == 0)
			{
			  h_mll_ee_inc->Fill(zmet.dilmass(),weight_mc);
			  
			  if(zmet.dilmass() > 81. && zmet.dilmass() < 101. /* && zmet.lep3().Pt() <= 10.*/)
				{
				  h_met_ee_inc->Fill(zmet.pfmet(),weight_mc);
				  h_met_ll_inc->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);					
				
				  if(zmet.njets()==0)
					{
					  h_met_ee_tar_njets0->Fill(zmet.pfmet(),weight_mc);
					  h_met_ll_tar_njets0->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);
					}
				  else if(zmet.njets()>=2)
					{
					  h_met_ee_tar_njets2->Fill(zmet.pfmet(),weight_mc);
					  h_met_ll_tar_njets2->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);
					}
				}
			  if (zmet.lep3().Pt() <= 10.  &&  zmet.njets()>=2  &&  zmet.nbcsvm()==0 && (zmet.mjj() > 70. && zmet.mjj() < 110.) )        
   				{h_mll_ee_tar->Fill(zmet.dilmass(),weight_mc);}
			}
		  else if (zmet.leptype()==1)
			{
			  h_mll_mumu_inc->Fill(zmet.dilmass(),weight_mc);
			
			  if(zmet.dilmass() > 81. && zmet.dilmass() < 101. /* && zmet.lep3().Pt() <= 10.*/)
				{
				  h_met_mumu_inc->Fill(zmet.pfmet(),weight_mc);
				  h_met_ll_inc->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);

				  if(zmet.njets()==0)
					{
					  h_met_mumu_tar_njets0->Fill(zmet.pfmet(),weight_mc);
					  h_met_ll_tar_njets0->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);
					}
				  else if(zmet.njets()>=2)
					{
					  h_met_mumu_tar_njets2->Fill(zmet.pfmet(),weight_mc);
					  h_met_ll_tar_njets2->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);
					}
				}

			  if (zmet.lep3().Pt() <= 10.  &&  zmet.njets()>=2  &&  zmet.nbcsvm()==0 && (zmet.mjj() > 70. && zmet.mjj() < 110.) )
   				{h_mll_mumu_tar->Fill(zmet.dilmass(),weight_mc);}
			}
 
		  else if (zmet.leptype()==2)
			{
			  // h_mll_mumu_inc->Fill(zmet.dilmass(),weight_mc);
			
			  if(zmet.dilmass() > 81. && zmet.dilmass() < 101. /* && zmet.lep3().Pt() <= 10.*/)
				{
				  h_met_emu_inc->Fill(zmet.pfmet(),weight_mc);
				  h_met_emu2_inc->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);

				  if(zmet.njets()==0)
					{
					  h_met_emu_tar_njets0->Fill(zmet.pfmet(),weight_mc);
					  h_met_emu2_tar_njets0->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);
					}
				  else if(zmet.njets()>=2)
					{
					  h_met_emu_tar_njets2->Fill(zmet.pfmet(),weight_mc);
					  h_met_emu2_tar_njets2->Fill(zmet.pfmet()*TMath::Cos( zmet.pfmetphi() ), zmet.pfmet()*TMath::Sin( zmet.pfmetphi() ), weight_mc);
					}
				}

			  // if (zmet.lep3().Pt() <= 10.  &&  zmet.njets()>=2  &&  zmet.nbcsvm()==0 && (zmet.mjj() > 70. && zmet.mjj() < 110.) )
			  //	{h_mll_mumu_tar->Fill(zmet.dilmass(),weight_mc);}
			}
	}

	
  

 //---------------------------------------------------------------------------------------------------------
    }    //end events loop.
 
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }


  //Write to file
  TFile *OutputFile = new TFile(Form("%s.root", skimFilePrefix.c_str()),"recreate");
  OutputFile->cd();
  h_mll_ee_inc->Write();
  h_mll_mumu_inc->Write();
  h_mll_ee_tar->Write();
  h_mll_mumu_tar->Write();
  h_met_ee_inc->Write();
  h_met_mumu_inc->Write();
  h_met_ll_inc->Write();
  h_met_emu_inc->Write();
  h_met_emu2_inc->Write();
  h_met_ee_tar_njets0->Write();
  h_met_mumu_tar_njets0->Write();
  h_met_ll_tar_njets0->Write();
  h_met_emu_tar_njets0->Write();
  h_met_emu2_tar_njets0->Write();
  h_met_ee_tar_njets2->Write();
  h_met_mumu_tar_njets2->Write();
  h_met_ll_tar_njets2->Write();
  h_met_emu_tar_njets2->Write();
  h_met_emu2_tar_njets2->Write();
  h_nvtx_scaled->Write();
  h_nvtx_unscaled->Write();
  h_phi->Write();


  OutputFile->Close();


  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;

  return 0;
}
