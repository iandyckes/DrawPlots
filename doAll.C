{

  gROOT->ProcessLine(".L ZMET.cc+");

  gROOT->ProcessLine(".L ScanChain.C+");//remove 2
  
  TChain *ch = new TChain("T1"); 
  // ch->Add("~benhoob/ZMet2012/output/V00-02-16/TChiWZ_madgraph_Zll_baby_oldIso.root");
  ch->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/zjets_53X_slim_baby.root");
  ScanChain(ch,true,-1,"zjets"); 


  //Run ScanChain for single top background
  TChain *ch_t = new TChain("T1");
  ch_t->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/t_53X_slim_baby.root ");
  ScanChain(ch_t,true,-1,"singlet"); 
 
 //Run ScanChain for WZ  background
  TChain *ch_wz = new TChain("T1");
  ch_wz->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/wz2l2q_53X_slim_baby.root ");
  ch_wz->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/wz3lnu_53X_slim_baby.root ");
  ScanChain(ch_wz,true,-1,"wz"); 

 //Run ScanChain for ZZ  background
  TChain *ch_zz = new TChain("T1");
  ch_zz->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/zz2l2nu_53X_slim_baby.root ");
  ch_zz->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/zz2l2q_53X_slim_baby.root ");
  ch_zz->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/zz4l_53X_slim_baby.root ");
  ScanChain(ch_zz,true,-1,"zz"); 

  //Run ScanChain for ttbar  background
  TChain *ch_ttbar = new TChain("T1");
  ch_ttbar->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/ttbar_mgdl_53X_slim_baby.root");
  ch_ttbar->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/ttbar_mgsl_53X_slim_baby.root");
  ch_ttbar->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/ttbar_mghad_53X_slim_baby.root");
  ScanChain(ch_ttbar,true,-1,"ttbar"); 

  //Run ScanChain for WW  background
  TChain *ch_ww = new TChain("T1");
  ch_ww->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/ww_53X_slim_baby.root");
  ScanChain(ch_ww,true,-1,"ww"); 

  //Run ScanChain for ttV  background
  TChain *ch_ttv = new TChain("T1");
  ch_ttv->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/ttw_53X_slim_baby.root ");
  ch_ttv->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/ttz_53X_slim_baby.root ");
  ScanChain(ch_ttv,true,-1,"ttv"); 

  //Run ScanChain for VVV  background
  TChain *ch_vvv = new TChain("T1");
  ch_vvv->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/vvv_53X_slim_baby.root");
  ScanChain(ch_vvv,true,-1,"vvv"); 

  //Run ScanChain for data
  TChain *ch_data = new TChain("T1");
  ch_data->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/data_53X_2012A_baby_nojson.root");
  ch_data->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/data_53X_2012B_baby_nojson.root");
  ch_data->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/data_53X_2012C_baby_nojson.root");
  ch_data->Add("/nfs-6/userdata/cwelke/Zmet2012/output/V00-02-19/data_53X_2012D_baby_nojson.root");
  ScanChain(ch_data,true,-1,"data"); 

}
