#include "LIP/TopTaus/interface/AnalysisMonitoring.hh"


void AnalysisMonitoring::buildStructures(){
 
 
  // if we have already processed info... ////////////////////// 
  if(outFile_){ 
    

    outFile_->Write(); outFile_->Close(); infoFile_.close(); dirsMaps_.clear();

    for(uint i=0;i<keys_.size();i++){
 
      TString key = keys_[i];
      (tauDileptonStepsMons_[key])->getAllMonitors().clear();
      (tauDileptonYieldsMons_[key])->getAllMonitors().clear();
      (tauDileptonDebugMons_[key])->getAllMonitors().clear();  
      (tauDileptonMCYieldsMons_[key])->getAllMonitors().clear();
      (tauDileptonBkgMons_[key])->getAllMonitors().clear();    
      (newPhysicsMons_[key])->getAllMonitors().clear();   

    }
  }
  /////////////////////////////////////////////////////////////


  //create the file ////////////////////////////////////////////////////////////////////
  outFile_ = TFile::Open( outFileName_ ,"RECREATE");// outFile_->SetCompressionLevel( 9 );
  //////////////////////////////////////////////////////////////////////////////////////
  

  //tauDileptonHisto_->
  buildHistograms( outFile_, tauDileptonStepsMons_ , tauDileptonYieldsMons_,tauDileptonMCYieldsMons_,  tauDileptonBkgMons_, tauDileptonDebugMons_,dirsMaps_ , keys_ );
  
  buildTrees(outFile_);
  
}

void AnalysisMonitoring::buildTrees(TFile* f){
  treeMap_.clear();
  
  
  vector<TString> channels;
  channels.push_back("e_tau");
  channels.push_back("m_tau");
  channels.push_back("lep_tau");
  
  vector<TString> dirs;
  dirs.push_back("Selected");
  dirs.push_back("DataDriven");
  
  
  vector<TString> branchNames;
  branchNames.push_back("evtnumb");
  branchNames.push_back("weight");
  
  // isOS
  branchNames.push_back("is_os");
  
  // tau R
  branchNames.push_back("rc_t");

  // transverse masses
  branchNames.push_back("mt_min");
  branchNames.push_back("mt_max");

  //add  weights for MVA (w+jet)
  //add  weights for MVA (qcd multijets)
  //add  average MVA
 

  //tau jet
  branchNames.push_back("pt_t"); branchNames.push_back("eta_t"); branchNames.push_back("radius_t"); 

  //lepton
  branchNames.push_back("pt_l"); branchNames.push_back("eta_l");

  //met
  branchNames.push_back("pt_met");  
  branchNames.push_back("Dphi_lepton_met");  
  branchNames.push_back("Dphi_tau_met");  

  //jets
  branchNames.push_back("multiplicity_j");  
  branchNames.push_back("btagmultiplicity_j");
  //branchNames.push_back("pt_j");    // leading jet filled
  //branchNames.push_back("eta_j");   // leading filled

  
      
  // Build Trees ////////////////////////////////////////////////////////////////////////////  
  for(size_t c=0;c<channels.size();c++){
    TString selectedTreeName          = channels[c]+TString("_Selected");
    TString dataDrivenTreeName        = channels[c]+TString("_DataDriven");
    TString selectedPlusTreeName      = channels[c]+TString("_Selected_plus");
    TString selectedMinusTreeName     = channels[c]+TString("_Selected_minus");

    TString selectedJerPlusTreeName   = channels[c]+TString("_Selected_jerplus");
    TString selectedJerMinusTreeName  = channels[c]+TString("_Selected_jerminus");

    TString selectedBPlusTreeName     = channels[c]+TString("_Selected_bplus");   
    TString selectedBMinusTreeName    = channels[c]+TString("_Selected_bminus");   
    TString selectedUnBPlusTreeName   = channels[c]+TString("_Selected_unbplus");
    TString selectedUnBMinusTreeName  = channels[c]+TString("_Selected_unbminus");      



    TString selectedUncPlusTreeName   = channels[c]+TString("_Selected_uncplus");   
    TString selectedUncMinusTreeName  = channels[c]+TString("_Selected_uncminus");
    

    // no syst
    TTree * selectedTree            = new TTree(selectedTreeName.Data(),selectedTreeName.Data());
    treeMap_[selectedTreeName]      = selectedTree;

    // Tree JES plus
    TTree * selectedPlusTree        = new TTree(selectedPlusTreeName.Data(),selectedPlusTreeName.Data());
    treeMap_[selectedPlusTreeName]  = selectedPlusTree;


    // Tree JES minus
    TTree * selectedMinusTree          = new TTree(selectedMinusTreeName.Data(),selectedMinusTreeName.Data());
    treeMap_[selectedMinusTreeName]    = selectedMinusTree;



    // Tree JER plus
    TTree * selectedJerPlusTree        = new TTree(selectedJerPlusTreeName.Data(),selectedJerPlusTreeName.Data());
    treeMap_[selectedJerPlusTreeName]  = selectedJerPlusTree;

    // Tree JER minus
    TTree * selectedJerMinusTree        = new TTree(selectedJerMinusTreeName.Data(),selectedJerMinusTreeName.Data());
    treeMap_[selectedJerMinusTreeName]  = selectedJerMinusTree;

    
    // Tree b plus
    TTree * selectedBPlusTree          = new TTree(selectedBPlusTreeName.Data(),selectedBPlusTreeName.Data());
    treeMap_[selectedBPlusTreeName]    = selectedBPlusTree;

    // Tree b minus
    TTree * selectedBMinusTree         = new TTree(selectedBMinusTreeName.Data(),selectedBMinusTreeName.Data());
    treeMap_[selectedBMinusTreeName]   = selectedBMinusTree;


    // Tree b plus
    TTree * selectedUnBPlusTree          = new TTree(selectedUnBPlusTreeName.Data(),selectedUnBPlusTreeName.Data());
    treeMap_[selectedUnBPlusTreeName]    = selectedUnBPlusTree;

    // Tree b minus
    TTree * selectedUnBMinusTree         = new TTree(selectedUnBMinusTreeName.Data(),selectedBMinusTreeName.Data());
    treeMap_[selectedUnBMinusTreeName]   = selectedUnBMinusTree;




    // Tree UNC plus
    TTree * selectedUncPlusTree        = new TTree(selectedUncPlusTreeName.Data(),selectedUncPlusTreeName.Data());
    treeMap_[selectedUncPlusTreeName]  = selectedUncPlusTree;

    // Tree UNC minus
    TTree * selectedUncMinusTree       = new TTree(selectedUncMinusTreeName.Data(),selectedUncMinusTreeName.Data());
    treeMap_[selectedUncMinusTreeName] = selectedUncMinusTree;


    // Data Driven
    TTree * dataDrivenTree          = new TTree(dataDrivenTreeName.Data(),dataDrivenTreeName.Data());
    treeMap_[dataDrivenTreeName]    = dataDrivenTree ;



  }
  ///////////////////////////////////////////////////////////////////////////////////////////

  // Build branches on the trees //////////////////////////////////////////////////////////////////////////////
  map<TString,TTree *>::iterator it;
  for( it=treeMap_.begin(); it!=treeMap_.end(); it++ ){
    TString treeName = (*it).first;
    for(size_t b = 0; b < branchNames.size(); b++){
      TString varQualifier  = branchNames[b]+TString("/D");   
      TBranch * myBranch  = ((*it).second)->Branch( (branchNames[b]).Data(), new Double_t ,varQualifier.Data()) ;
    }
   
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////



}


TGraphAsymmErrors*  AnalysisMonitoring::scriptForPlotWithSystematics( TString plot, int level, string plotType="mc"){
  
  
  //////////////////////////////////////////////////////////////////////////
  // WARNING : 
  // NOTE 1 : This script does not take into account systemeatics on tau id 
  //////////////////////////////////////////////////////////////////////////
  // NOTE 2 ; TgraphAsymmErrors should be ploted as   grae->Draw("e2p");
  //////////////////////////////////////////////////////////////////////////


  TString outScript("systematicsPlot.C");

  char numb[100];  sprintf(numb,"_%d",level);

  TString pl   = plot+TString(numb);
  TString pl_p = plot+TString("_plus")+TString(numb);
  TString pl_m = plot+TString("_minus")+TString(numb);

  cout<<endl<<" Processing files for systematics : "<<pl<<"\nplus  : "<<pl_p<<"\nminus : "<<pl_m<<endl;

  // Build h_plus and h_minus //////////////////////////////////////////////////////////////////////
  TH1D *histo, *histo_p, *histo_m;

  initProcessedFiles();

  if(plotType == "data"){
    TH1D *  h   = (TH1D * ) listOfProcessedFiles_[0]->Get(pl);    histo   = (TH1D *) h->Clone();
    TH1D *  h_p = (TH1D * ) listOfProcessedFiles_[0]->Get(pl_p);  histo_p = (TH1D *) h_p->Clone();
    TH1D *  h_m = (TH1D * ) listOfProcessedFiles_[0]->Get(pl_m);  histo_m = (TH1D *) h_m->Clone();
  }
  else {
    TH1D *  h   = (TH1D * ) listOfProcessedFiles_[1]->Get(pl);    histo   = (TH1D *) h->Clone();
    TH1D *  h_p = (TH1D * ) listOfProcessedFiles_[1]->Get(pl_p);  histo_p = (TH1D *) h_p->Clone();
    TH1D *  h_m = (TH1D * ) listOfProcessedFiles_[1]->Get(pl_m);  histo_m = (TH1D *) h_m->Clone();
    // Adding scaled plots
    for( uint i =2; i<listOfProcessedFiles_.size();i++){  // skeep the first plot (data) and cloned plot(1)

      cout<<endl<<" going to process sample i = "<<i<<" "<<endl;

      TH1D *  temp   = (TH1D * ) listOfProcessedFiles_[i]->Get(pl);
      TH1D *  temp_p = (TH1D * ) listOfProcessedFiles_[i]->Get(pl_p);
      TH1D *  temp_m = (TH1D * ) listOfProcessedFiles_[i]->Get(pl_m);
      histo->Add(temp);  histo_p->Add(temp_p); histo_m->Add(temp_m);
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  // build TGraph /////////////////////////////////////////
  TGraphAsymmErrors * graph = new TGraphAsymmErrors(histo);
  graph->SetFillColor(kBlack);
  graph->SetFillStyle(3004);
  /////////////////////////////////////////////////////////


  // Deal with ttbar  ///////////////////////////////////////////////////////////////////////////
  int ttbarSample(1); TH1D * histo_ttbar  = (TH1D * ) listOfProcessedFiles_[ttbarSample]->Get(pl); 
  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Deal with singletop  ////////////////////////////////////////////////////////////////////
  int singleTopSample(2); TH1D * histo_singletop  = (TH1D * ) listOfProcessedFiles_[singleTopSample]->Get(pl); 
  ////////////////////////////////////////////////////////////////////////////////////////////////
    
  // Deal with wjets  ///////////////////////////////////////////////////////////////////////////
  int wjetsSample(3); TH1D * histo_wjets  = (TH1D * ) listOfProcessedFiles_[wjetsSample]->Get(pl); 
  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Deal with zjets  ///////////////////////////////////////////////////////////////////////////
  int zjetsSample(4); TH1D * histo_zjets  = (TH1D * ) listOfProcessedFiles_[zjetsSample]->Get(pl); 
  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Deal with diboson  ////////////////////////////////////////////////////////////////////////////////////
  int dibosonSample(5); TH1D * histo_diboson  = (TH1D * ) listOfProcessedFiles_[dibosonSample]->Get(pl); 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Deal with QCD  /////////////////////////////////////////////////////////////////////
  int qcdSample(6); TH1D * histo_qcd  = (TH1D * ) listOfProcessedFiles_[qcdSample]->Get(pl); 
  ////////////////////////////////////////////////////////////////////////////////////////

  // set systematics errors //////////////////////////////////////////////////////////////////////////////////////////////////
  for( uint i=1; i<=(uint)histo->GetNbinsX(); i++){

    double n_histo_qcd       = histo_qcd->GetBinContent(i);         double xsec_err_qcd      = QCD_CS_ERR*n_histo_qcd;
    double n_histo_ttbar     = histo_ttbar->GetBinContent(i);       double xsec_err_ttbar    = TTBAR_CS_ERR*n_histo_ttbar;
    double n_histo_singletop = histo_singletop->GetBinContent(i);   double xsec_err_sinletop = SINGLETOP_CS_ERR*n_histo_singletop;
    double n_histo_diboson   = histo_diboson->GetBinContent(i);     double xsec_err_diboson  = OTHER_CS_ERR*n_histo_diboson;
    double n_histo_wjets     = histo_wjets->GetBinContent(i);       double xsec_err_wjets    = WJETS_CS_ERR*n_histo_wjets;
    double n_histo_zjets     = histo_zjets->GetBinContent(i);       double xsec_err_zjets    = OTHER_CS_ERR*n_histo_zjets;

    double n_histo           = histo->GetBinContent(i);
    double n_histo_p         = histo_p->GetBinContent(i);
    double n_histo_m         = histo_m->GetBinContent(i);

    double err_lum  = LUM_ERR*(n_histo);
    double err_lep  = LEP_EFF_ERR*(n_histo);

    if(plotType == "data"){
      xsec_err_qcd      = 0;
      xsec_err_ttbar    = 0;
      xsec_err_sinletop = 0;
      xsec_err_diboson  = 0;
      xsec_err_wjets    = 0;
      xsec_err_zjets    = 0;
    }

    vector<double> errors;
    double error_p(0), error_m(0);

    if(n_histo_p > n_histo_m ){ error_p = n_histo_p - n_histo; error_m = n_histo - n_histo_m; }
    else { error_p = n_histo_m - n_histo; error_m = n_histo - n_histo_p; }

    if(error_p<0 ){ error_p=0;}
    if(error_m<0 ){ error_m=0;}


    double temp_err = pow(xsec_err_qcd,2)+pow(xsec_err_ttbar,2)+pow(xsec_err_sinletop,2)+pow(xsec_err_diboson,2)+pow(xsec_err_wjets,2)+pow(xsec_err_zjets,2);
    temp_err += temp_err + pow(err_lum,2) + pow (err_lep,2);

    double error_sist_p = sqrt( temp_err + error_p*error_p);
    double error_sist_m = sqrt( temp_err + error_m*error_m);

    graph->SetPointEYhigh(i-1, error_sist_p);
    graph->SetPointEYlow(i-1,error_sist_m);
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
  Double_t c_w = 600; Double_t c_h = 600;
  TCanvas * c1 = new TCanvas("c", "c", c_w, c_h);
  c1->cd();
  graph->Draw("a2p");
  c1->SaveAs(outScript);
  
  return graph;
}



