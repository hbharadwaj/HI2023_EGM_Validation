#include <TSystem.h>        // needed for gSystem
#include <TStyle.h>         // needed for gStyle 
#include <TChain.h>         // needed for TChain
#include "TROOT.h"          // needed for gROOT
#include <TError.h>         // needed for gErrorIgnoreLevel
#include <TVector2.h>       // needed for TVector2::Phi_mpi_pi        
#include <TH1.h>            // needed fror TH1 and TH2
#include <TH2.h>            // needed fror TH1 and TH2
#include <TLegend.h>        // needed for Legend
#include <TCanvas.h>        // needed for Canvas
#include <TMath.h>          //! needed for floating values in plots for some reason
#include <THStack.h>        // needed for THStack
#include <TLatex.h>         // needed for TLatex
#include <TFile.h>          // needed for TFile

#include <TEfficiency.h>            // needed for TEfficiency
#include <TGraphAsymmErrors.h>      // needed for TGraphAsymmErrors
#include <TLine.h>                  // needed for TLine
#include <TSystemDirectory.h>       // needed for TSystemDirectory
#include <TSystemFile.h>            // needed for TSystemFile

#include <iostream>         // needed for I/O
#include <string>           // needed for string

using namespace std;

TString label="Run_374719_374810";//"2023_Run3_QCDPhoton_Embedded";//"2023_Run3_EmEnrichedDijet30_Embedded"; // "HiMinimumBias2_Data_2018";
TString Tag = "_PromptReco_2023_10_11";
TString output_path = "./";
TFile *fout_2018, *fout_2023, *fout_temp;
// TCanvas c;

bool flag_2018=false;    // To run over 2018 Data
bool flag_MC=false;      // To run over 2023 MC
bool flag_plot=true;    // To plot resulting histograms 
bool flag_ratio=true;   // To plot the ratios. If false only the distributions are shown

// Eta ranges considered for the study {|eta|<1.44, 1.566<|eta|<2.1, 0<|eta|<2.4}  // index ieta
std::vector<float> min_eta = {  0,  1.566,   0};   
std::vector<float> max_eta = {1.44,   2.1, 2.4};
const std::size_t neta = 3;

// PbPb centrality ranges considered for the study // index icent
std::vector<int> min_cent = {  0,  0,  60};   
std::vector<int> max_cent = {180, 60, 180};
const std::size_t ncent = 3;  

// personal histogram plotting script 
void Plot_hist(std::vector<TH1D*>,std::vector<TString> ,TString opt="label",std::vector<TString> eopt={"end"});
void overlay_runs(std::vector<TH1D*>,std::vector<TString> ,TString opt="label",std::vector<TString> eopt={"end"});
void GetFiles(char const* input,std::vector<TString>& files, int file_limit=99999);
void FillChain(TChain& chain,std::vector<TString>& files);
void displayProgress(long current, long max);
double dr(float eta1, float phi1, float eta2, float phi2) {
    float deta = TMath::Abs(eta1 - eta2);
    float dphi = TMath::Abs(phi1 - phi2);
    if (dphi > TMath::Pi()) dphi = TMath::Abs(dphi - 2*TMath::Pi());

    return TMath::Sqrt(dphi*dphi + deta*deta);
}

void Data_fill_pho_ele_2018(std::vector<TString> in_file_path,TString in_label, Bool_t isData=true);
void Data_fill_pho_ele_2023(std::vector<TString> in_file_path,TString in_label, Bool_t isData=true);

void Data_overlay_2018(){

    gROOT->SetBatch();
    gErrorIgnoreLevel = kFatal;
    TString DIR = output_path + "plots_2018/plots_2018_"+label+Tag+"/";
    output_path = DIR;
    TString makedir = "mkdir -p " + DIR;
    const char *mkDIR = makedir.Data();
    
    gSystem->Exec(mkDIR);

    int nfiles = 9999; // Can reduce to test if the script runs
    std::vector<TString> files;
    
    if(flag_2018){
        fout_2018 = new TFile(DIR+"/Output_2018Data"+Tag+".root", "UPDATE");
        fout_2018->cd();

        files.clear();
        TString input_dir_2018  = "/data_CMS/cms/bharikri/2018_Data_EGM_Validation/";
        GetFiles(input_dir_2018.Data(), files,nfiles>100?100:nfiles);
        std::cout<<"Got "<<files.size()<<" files for 2018 Data \n";
        Data_fill_pho_ele_2018(files,"2018_Data",true);
        std::cout<<Form("Output_2018Data_%s%s.root",label.Data(),Tag.Data())<<" file created\n";
        fout_2018->Close();
        delete fout_2018;
    }

    if(flag_MC){
        fout_2023 = new TFile(DIR+"/Output_"+label+Tag+".root", "UPDATE");
        fout_2023->cd();

        files.clear();
        TString input_dir_MC  = "/data_CMS/cms/bharikri/ECAL_study/old/QCDPhoton_30/Run3_Reference/";
        GetFiles(input_dir_MC.Data(), files,nfiles);
        std::cout<<"Got "<<files.size()<<" files for MC \n";
        Data_fill_pho_ele_2023(files,"QCDPhoton_MC_ppNom",false);
        std::cout<<Form("Output_%s%s.root",label.Data(),Tag.Data())<<" file created\n";
        fout_2023->Close();
        delete fout_2023;
    }

    fout_2023 = new TFile(DIR+"/Output_"+label+Tag+".root", "UPDATE");
    fout_2023->cd();

    files.clear();
    TString input_dir_promptreco  = "/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374719_HIPhysicsRawPrime2_2023-10-06/0000/";
    GetFiles(input_dir_promptreco.Data(), files,nfiles);
    input_dir_promptreco  = "/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374719_HIPhysicsRawPrime1_2023-10-06/0000/";
    GetFiles(input_dir_promptreco.Data(), files,nfiles);
    input_dir_promptreco  = "/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374719_HIPhysicsRawPrime0_2023-10-06/0000/";
    GetFiles(input_dir_promptreco.Data(), files,nfiles);

    input_dir_promptreco  = "/eos/cms/store/group/phys_heavyions/mstojano/run3RapidValidation/PbPb2023_run374810_HIPhysicsRawPrime0_Forest_2023-10-09/CRAB_UserFiles/crab_PbPb2023_run374810_HIPhysicsRawPrime0_Forest_2023-10-09/231009_162633/";
    GetFiles(input_dir_promptreco.Data(), files,nfiles);

    std::cout<<"Got "<<files.size()<<" files for run374719_374810_PromptHIPhysicsRawPrime\n";
    Data_fill_pho_ele_2023(files,"run374719_374810_PromptHIPhysicsRawPrime",true);
    std::cout<<"Output_"<<label<<Tag<<".root file created\n";
    fout_2023->Close();
    // delete fout_2023;

    fout_2023 = TFile::Open(DIR+"/Output_"+label+Tag+".root", "UPDATE");

    // Overlay
    if(flag_plot){

        std::vector<std::vector<TString>>dirlist_name = {
            {DIR+"/Output_"+label+Tag+".root","run374719_374810_PromptHIPhysicsRawPrime","Physics Private Reco"},
            {DIR+"/Output_2018Data"+Tag+".root","2018_Data","2018 Data"},
            // {DIR+"/Output_"+label+Tag+".root","QCDPhoton_MC_ppNom","2023 QCDPhoton MC-old ECAL"},
        };

        std::vector<TString>histlist = {
            "h_lead_pho_et",
            "h_lead_pho_eta",
            "h_lead_pho_SIEIE",
            "h_lead_pho_phi",
            "h_lead_pho_HoverE",
            "h_lead_pho_SumCalIso",
            "h_lead_pho_ECALIso",
            "h_lead_pho_HCALIso",
            "h_lead_pho_TRKIso",
            "h_lead_pho_PFPIso",
            "h_lead_pho_PFNIso",
            "h_lead_pho_PFCIso",
            "h_lead_pho_R9",

            "h_lead_ele_Pt",
            "h_lead_ele_eta",
            "h_lead_ele_SIEIE",
            "h_lead_ele_phi",
            "h_lead_ele_HoverE",
            "h_lead_ele_EoverP",
            "h_lead_ele_dEtaSeedAtVtx",
            "h_lead_ele_dPhiAtVtx",
            "h_lead_ele_EoverPInv",

            "h_sublead_ele_Pt",
            "h_sublead_ele_eta",
            "h_sublead_ele_SIEIE",
            "h_sublead_ele_phi",
            "h_sublead_ele_HoverE",
            "h_sublead_ele_EoverP",
            "h_sublead_ele_dEtaSeedAtVtx",
            "h_sublead_ele_dPhiAtVtx",
            "h_sublead_ele_EoverPInv",
        };

        for (std::size_t icent = 0; icent <= 0; ++icent){ // ncent {0-100, 0-30, 30-100} 
       
            for (std::size_t ieta = 0; ieta <=2; ++ieta){ //neta   {|eta|<1.44, 1.566<|eta|<2.1, 0<|eta|<2.4} 
                int temp_hist_index = 1;

                std::vector<TString>sel = {"",Form("Cent. %d-%d%%",min_cent[icent]/2,max_cent[icent]/2),Form("%1.3f<|#eta|<%1.1f",min_eta[ieta],max_eta[ieta]),"p_{T}>40 + HLT40","Loose ID"}; // ,"NOT GEN MATCHED"

                for(TString input_hist:histlist){
                    TString outplotname = input_hist;   

                    if(input_hist.Contains("sublead_ele"))
                        sel = {"",Form("Cent. %d-%d%%",min_cent[icent]/2,max_cent[icent]/2),Form("%1.3f<|#eta|<%1.1f",min_eta[ieta],max_eta[ieta]),"p_{T}>20","Any Double Ele trigger","Loose ID"}; // ,"NOT GEN MATCHED"
                    
                    TString input_histname = Form("%s_%zu_%zu",outplotname.Data(),ieta,icent);
                    // std::cout<<"input hist = "<<input_hist<<"\n";

                    std::vector<TH1D*> h;
                    std::vector<TString>hname;
                    bool flag_ratio_hist = flag_ratio;

                    for(std::size_t i_file = 0; i_file<dirlist_name.size();i_file++){
                        fout_temp = TFile::Open(dirlist_name[i_file][0], "UPDATE");
                        if(!fout_temp){std::cout<<dirlist_name[i_file][0]<<" does not exist\n"; continue;}

                        TH1D *h_temp = (TH1D*)fout_temp->Get(Form("%s/%s",dirlist_name[i_file][1].Data(),input_histname.Data()));
                        if(!h_temp){ std::cout<<input_histname<<" \t histogram doesn't exist\n"; break; }
                        if(h_temp->GetEntries()==0){ flag_ratio_hist = false; continue;}
                        h.push_back((TH1D*)h_temp->Clone());
                        hname.push_back(dirlist_name[i_file][2]);
                    }
                    hname.push_back(Form("canvas_%s",input_histname.Data()));

                    fout_2023->cd();
                    if(dirlist_name.size()>=2 && flag_ratio_hist){
                        overlay_runs(h,hname,"right_norm_width",sel);
                    }
                    else{
                        Plot_hist(h,hname,"right_norm_width",sel);
                    }
                    TCanvas *c_temp = (TCanvas*)fout_2023->Get(Form("canvas_%s",input_histname.Data()));
                    if(!c_temp){ std::cout<<"Combined PDF not saved\n"; break; }
                    c_temp->Draw();
                    if(temp_hist_index==1){
                        // c_temp->Print("Plot_Photon_electron.pdf");
                        c_temp->Print(Form("%sPlot_Photon_electron_2018_%s_%zu_%zu.pdf(",output_path.Data(),label.Data(),ieta,icent),Form("Title:%s",outplotname.Data()));
                    }
                    else if(temp_hist_index==(int)histlist.size()){
                        c_temp->Print(Form("%sPlot_Photon_electron_2018_%s_%zu_%zu.pdf)",output_path.Data(),label.Data(),ieta,icent),Form("Title:%s",outplotname.Data()));
                        std::cout<<Form("Combined Plot_Photon_electron_%s_%zu_%zu.pdf has been saved",label.Data(),ieta,icent)<<"\n";
                        std::cout<<"-------------------------------------------------------------------"<<"\n";
                    }
                    else{                    
                        c_temp->Print(Form("%sPlot_Photon_electron_2018_%s_%zu_%zu.pdf",output_path.Data(),label.Data(),ieta,icent),Form("Title:%s",outplotname.Data()));
                    }

                    temp_hist_index++;

                } // End hist list loop
            }// End eta loop
        } // End centrality loop
    }

    std::cout<<"\n";
    fout_2023->Close();
}

void Data_fill_pho_ele_2018(std::vector<TString> in_file_path,TString in_label, Bool_t isData){
    // in_file_path = Directory name for the specified Data Run
    // in_label = Corresponding legend label

    TChain photonTree("photonTree"), HiTree("HiTree"),skimanalysis("skimanalysis"),hltanalysis("hltanalysis");

    for(TString temp_path:in_file_path){
        // std:cout<<"file path = "<<temp_path+input_pho_tree<<"\n";

        photonTree.Add(temp_path+"/ggHiNtuplizer/EventTree");
        HiTree.Add(temp_path+"/hiEvtAnalyzer/HiTree");
        skimanalysis.Add(temp_path+"/skimanalysis/HltTree");
        hltanalysis.Add(temp_path+"/hltanalysis/HltTree");
    }

    // Centrality weights 
    const float Ncoll[200] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
   
    // ----------------------------------------------------------------------------------------------------------------
    // Variables

        Int_t hiBin=-1;
        Float_t hiHF = -1; 
        Float_t weight = 1; 
        Float_t pthat=0;

        HiTree.LoadTree(0);

        HiTree.SetBranchAddress("hiBin", &hiBin);
        HiTree.SetBranchAddress("hiHF", &hiHF);
        if(HiTree.GetBranch("weight") && !isData)
            HiTree.SetBranchAddress("weight", &weight);
        if(HiTree.GetBranch("pthat") && !isData)
            HiTree.SetBranchAddress("pthat", &pthat);

        // reco::GenParticle
        std::vector<int>*    mcStatus = 0;
        std::vector<int>*    mcPID = 0;
        std::vector<int>*    mcMomPID = 0;
        std::vector<float>*  mcPt = 0;
        std::vector<float>*  mcEta = 0;
        std::vector<float>*  mcPhi = 0;
        std::vector<float>*  mcE = 0;
        std::vector<float>*  mcEt = 0;
        std::vector<float>*  mcMass = 0;
        std::vector<float>*  mcCalIsoDR03 = 0;
        std::vector<float>*  mcCalIsoDR04 = 0;
        std::vector<float>*  mcTrkIsoDR03 = 0;
        std::vector<float>*  mcTrkIsoDR04 = 0;

        std::vector<int>*    pho_genMatchedIndex = 0;

        // reco::Photon
        std::vector<float>*  phoE = 0;
        std::vector<float>*  phoEt = 0;
        std::vector<float>*  phoEta = 0;
        std::vector<float>*  phoPhi = 0;
        std::vector<float>*  phoHoverE = 0;
        std::vector<float>*  pho_ecalClusterIsoR3 = 0;
        std::vector<float>*  pho_hcalRechitIsoR3 = 0;
        std::vector<float>*  pho_trackIsoR3PtCut20 = 0;
        std::vector<float>*  phoSigmaEtaEta_2012 = 0;

        std::vector<float>*  pfcIso3subUE = 0;
        std::vector<float>*  pfnIso3subUE = 0;
        std::vector<float>*  pfpIso3subUE = 0; //! For MiniAOD by default this does not exclude the cone!!

        std::vector<float>*  pfpIso3subSCsubUE = 0; //! For MiniAOD by default this does not exclude the cone!!

        std::vector<float>*  pho_swissCrx = 0;
        std::vector<float>*  pho_seedTime = 0;

        std::vector<float>* phoSCRawE = 0;
        std::vector<float>* phoSCEta = 0;
        std::vector<float>* phoSCPhi = 0;
        std::vector<float>* phoSCEtaWidth = 0;
        std::vector<float>* phoSCPhiWidth = 0;
        std::vector<float>* phoE3x3_2012 = 0;
        std::vector<float>* phoMaxEnergyXtal_2012 = 0;
        std::vector<float>* phoE2nd_2012 = 0;
        std::vector<float>* phoELeft_2012 = 0;
        std::vector<float>* phoERight_2012 = 0;
        std::vector<float>* phoETop_2012 = 0;
        std::vector<float>* phoEBottom_2012 = 0;
        std::vector<float>* phoSigmaIEtaIEta_2012 = 0;
        std::vector<float>* phoSigmaIEtaIPhi_2012 = 0;
        std::vector<float>* phoSigmaIPhiIPhi_2012 = 0;
        Float_t rho = 0;
        std::vector<float>* phoESEn = 0;

        std::vector<int>* phoSCnHits = 0;
        std::vector<float>* phoE5x5_2012 = 0;
        std::vector<float>* phoHadTowerOverEm = 0;
        std::vector<float>* phoHadTowerOverEm1 = 0;
        std::vector<float>* phoHadTowerOverEm2 = 0;
        std::vector<float>* phoR9_2012 = 0;

        // reco::Electron
        std::vector<float>* elePt = 0;
        std::vector<float>* eleEta = 0;
        std::vector<float>* elePhi = 0;
        std::vector<float>* eleEoverP = 0;

        std::vector<float>* eleSigmaIEtaIEta_2012 = 0;
        std::vector<float>* eledEtaSeedAtVtx = 0;
        std::vector<float>* eledPhiAtVtx = 0;
        std::vector<float>* eleHoverE = 0;
        std::vector<float>* eleEoverPInv = 0;
        std::vector<int>* eleMissHits = 0;
        std::vector<int>* eleCharge = 0;
        std::vector<float>* eleIP3D = 0;
        std::vector<float>* eleEn = 0;

        photonTree.LoadTree(0);

        if(photonTree.GetBranch("mcPID") && !isData){

            photonTree.SetBranchAddress("mcStatus",     &mcStatus);
            photonTree.SetBranchAddress("mcPID",        &mcPID);
            photonTree.SetBranchAddress("mcPt",         &mcPt);
            photonTree.SetBranchAddress("mcEta",        &mcEta);
            photonTree.SetBranchAddress("mcPhi",        &mcPhi);
            photonTree.SetBranchAddress("mcE",          &mcE);
            photonTree.SetBranchAddress("mcEt",         &mcEt);   
            photonTree.SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03);
            photonTree.SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04);
            photonTree.SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03);
            photonTree.SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04);
            photonTree.SetBranchAddress("mcMass",       &mcMass);
            photonTree.SetBranchAddress("mcMomPID",     &mcMomPID);
            
            photonTree.SetBranchAddress("pho_genMatchedIndex",   &pho_genMatchedIndex);
        }

        photonTree.SetBranchAddress("phoE",                    &phoE);
        photonTree.SetBranchAddress("phoEt",                   &phoEt);
        photonTree.SetBranchAddress("phoEta",                  &phoEta);
        photonTree.SetBranchAddress("phoPhi",                  &phoPhi);
        photonTree.SetBranchAddress("phoHoverE",               &phoHoverE);
        photonTree.SetBranchAddress("pho_ecalClusterIsoR3",    &pho_ecalClusterIsoR3);
        photonTree.SetBranchAddress("pho_hcalRechitIsoR3",     &pho_hcalRechitIsoR3);
        photonTree.SetBranchAddress("pho_trackIsoR3PtCut20",   &pho_trackIsoR3PtCut20);
        photonTree.SetBranchAddress("phoSigmaEtaEta_2012",     &phoSigmaEtaEta_2012);
        photonTree.SetBranchAddress("pfcIso3subUE",            &pfcIso3subUE);
        photonTree.SetBranchAddress("pfnIso3subUE",            &pfnIso3subUE);
        photonTree.SetBranchAddress("pfpIso3subUE",            &pfpIso3subUE);

        photonTree.SetBranchAddress("pfpIso3subSCsubUE",       &pfpIso3subSCsubUE);

        photonTree.SetBranchAddress("phoSCRawE",               &phoSCRawE);
        photonTree.SetBranchAddress("phoSCEta",                &phoSCEta);
        photonTree.SetBranchAddress("phoSCPhi",                &phoSCPhi);
        photonTree.SetBranchAddress("phoSCEtaWidth",           &phoSCEtaWidth);
        photonTree.SetBranchAddress("phoSCPhiWidth",           &phoSCPhiWidth);
        photonTree.SetBranchAddress("phoE3x3_2012",            &phoE3x3_2012);
        photonTree.SetBranchAddress("phoMaxEnergyXtal_2012",   &phoMaxEnergyXtal_2012);
        photonTree.SetBranchAddress("phoE2nd_2012",            &phoE2nd_2012);
        photonTree.SetBranchAddress("phoELeft_2012",           &phoELeft_2012);
        photonTree.SetBranchAddress("phoERight_2012",          &phoERight_2012);
        photonTree.SetBranchAddress("phoETop_2012",            &phoETop_2012);
        photonTree.SetBranchAddress("phoEBottom_2012",         &phoEBottom_2012);
        photonTree.SetBranchAddress("phoSigmaIEtaIEta_2012",   &phoSigmaIEtaIEta_2012);
        photonTree.SetBranchAddress("phoSigmaIEtaIPhi_2012",   &phoSigmaIEtaIPhi_2012);
        photonTree.SetBranchAddress("phoSigmaIPhiIPhi_2012",   &phoSigmaIPhiIPhi_2012);
        photonTree.SetBranchAddress("rho",                     &rho);
        photonTree.SetBranchAddress("phoESEn",                 &phoESEn);

        photonTree.SetBranchAddress("phoSCnHits",              &phoSCnHits);
        photonTree.SetBranchAddress("phoE5x5_2012",            &phoE5x5_2012);
        photonTree.SetBranchAddress("phoHadTowerOverEm",       &phoHadTowerOverEm);
        photonTree.SetBranchAddress("phoHadTowerOverEm1",      &phoHadTowerOverEm1);
        photonTree.SetBranchAddress("phoHadTowerOverEm2",      &phoHadTowerOverEm2);
        photonTree.SetBranchAddress("phoR9_2012",              &phoR9_2012);

        photonTree.SetBranchAddress("pho_swissCrx",            &pho_swissCrx);
        photonTree.SetBranchAddress("pho_seedTime",            &pho_seedTime);

        photonTree.SetBranchAddress("elePt",                   &elePt);
        photonTree.SetBranchAddress("eleEta",                  &eleEta);
        photonTree.SetBranchAddress("elePhi",                  &elePhi);
        photonTree.SetBranchAddress("eleEoverP",               &eleEoverP);

        photonTree.SetBranchAddress("eleSigmaIEtaIEta_2012",   &eleSigmaIEtaIEta_2012);
        photonTree.SetBranchAddress("eledEtaSeedAtVtx",        &eledEtaSeedAtVtx);
        photonTree.SetBranchAddress("eledPhiAtVtx",            &eledPhiAtVtx);
        photonTree.SetBranchAddress("eleHoverE",               &eleHoverE);
        photonTree.SetBranchAddress("eleEoverPInv",            &eleEoverPInv);
        photonTree.SetBranchAddress("eleMissHits",             &eleMissHits);
        photonTree.SetBranchAddress("eleCharge",               &eleCharge);
        photonTree.SetBranchAddress("eleIP3D",                 &eleIP3D);
        photonTree.SetBranchAddress("eleEn",                   &eleEn);

        int pprimaryVertexFilter = 1;               // Pass if it doesn't exist in the Forest
        int pclusterCompatibilityFilter = 1;
        int pphfCoincFilter2Th4 = 1;

        skimanalysis.SetBranchAddress("pprimaryVertexFilter",  &pprimaryVertexFilter);
        skimanalysis.SetBranchAddress("pclusterCompatibilityFilter",  &pclusterCompatibilityFilter);
        skimanalysis.SetBranchAddress("pphfCoincFilter2Th4",  &pphfCoincFilter2Th4);

        int HLT_HIGEDPhoton20_v1        = 0;
        int HLT_HIGEDPhoton30_v1        = 0;
        int HLT_HIGEDPhoton40_v1        = 0;

        int HLT_HIEle20Gsf_v1           = 0;
        int HLT_HIEle30Gsf_v1           = 0;
        int HLT_HIEle40Gsf_v1           = 0;
        int HLT_HIEle50Gsf_v1           = 0;
        int HLT_HIEle15Ele10Gsf_v1      = 0;
        int HLT_HIEle15Ele10GsfMass50_v1= 0;   

        int HLT_HIDoubleEle10Gsf_v1         = 0;
        int HLT_HIDoubleEle10GsfMass50_v1   = 0;
        int HLT_HIDoubleEle15Gsf_v1         = 0;
        int HLT_HIDoubleEle15GsfMass50_v1   = 0;        

        hltanalysis.SetBranchAddress("HLT_HIGEDPhoton20_v1",           &HLT_HIGEDPhoton20_v1);
        hltanalysis.SetBranchAddress("HLT_HIGEDPhoton30_v1",           &HLT_HIGEDPhoton30_v1);
        hltanalysis.SetBranchAddress("HLT_HIGEDPhoton40_v1",           &HLT_HIGEDPhoton40_v1);

        hltanalysis.SetBranchAddress("HLT_HIEle20Gsf_v1",              &HLT_HIEle20Gsf_v1);
        hltanalysis.SetBranchAddress("HLT_HIEle30Gsf_v1",              &HLT_HIEle30Gsf_v1);
        hltanalysis.SetBranchAddress("HLT_HIEle40Gsf_v1",              &HLT_HIEle40Gsf_v1);
        hltanalysis.SetBranchAddress("HLT_HIEle50Gsf_v1",              &HLT_HIEle50Gsf_v1);
        hltanalysis.SetBranchAddress("HLT_HIEle15Ele10Gsf_v1",         &HLT_HIEle15Ele10Gsf_v1);
        hltanalysis.SetBranchAddress("HLT_HIEle15Ele10GsfMass50_v1",   &HLT_HIEle15Ele10GsfMass50_v1); 

        hltanalysis.SetBranchAddress("HLT_HIDoubleEle10Gsf_v1",        &HLT_HIDoubleEle10Gsf_v1);
        hltanalysis.SetBranchAddress("HLT_HIDoubleEle10GsfMass50_v1",  &HLT_HIDoubleEle10GsfMass50_v1);
        hltanalysis.SetBranchAddress("HLT_HIDoubleEle15Gsf_v1",        &HLT_HIDoubleEle15Gsf_v1);
        hltanalysis.SetBranchAddress("HLT_HIDoubleEle15GsfMass50_v1",  &HLT_HIDoubleEle15GsfMass50_v1);

    // -------- End Tree Variable Declaration
    // ----------------------------------------------------------------------------------------------------------------
    // Histograms
        TH1::SetDefaultSumw2();
        TH2::SetDefaultSumw2();

        TH1D* h_lead_pho_et[neta][ncent];
        TH1D* h_lead_pho_eta[neta][ncent];
        TH1D* h_lead_pho_SIEIE[neta][ncent];
        TH1D* h_lead_pho_phi[neta][ncent];
        TH1D* h_lead_pho_HoverE[neta][ncent];
        TH1D* h_lead_pho_SumCalIso[neta][ncent];
        TH1D* h_lead_pho_ECALIso[neta][ncent];
        TH1D* h_lead_pho_HCALIso[neta][ncent];
        TH1D* h_lead_pho_TRKIso[neta][ncent];
        TH1D* h_lead_pho_PFPIso[neta][ncent];
        TH1D* h_lead_pho_PFNIso[neta][ncent];
        TH1D* h_lead_pho_PFCIso[neta][ncent];
        TH1D* h_lead_pho_R9[neta][ncent];

        TH1D* h_lead_ele_Pt[neta][ncent];
        TH1D* h_lead_ele_eta[neta][ncent];
        TH1D* h_lead_ele_SIEIE[neta][ncent];
        TH1D* h_lead_ele_phi[neta][ncent];
        TH1D* h_lead_ele_HoverE[neta][ncent];
        TH1D* h_lead_ele_EoverP[neta][ncent];
        TH1D* h_lead_ele_dEtaSeedAtVtx[neta][ncent];
        TH1D* h_lead_ele_dPhiAtVtx[neta][ncent];
        TH1D* h_lead_ele_EoverPInv[neta][ncent];

        TH1D* h_sublead_ele_Pt[neta][ncent];
        TH1D* h_sublead_ele_eta[neta][ncent];
        TH1D* h_sublead_ele_SIEIE[neta][ncent];
        TH1D* h_sublead_ele_phi[neta][ncent];
        TH1D* h_sublead_ele_HoverE[neta][ncent];
        TH1D* h_sublead_ele_EoverP[neta][ncent];
        TH1D* h_sublead_ele_dEtaSeedAtVtx[neta][ncent];
        TH1D* h_sublead_ele_dPhiAtVtx[neta][ncent];
        TH1D* h_sublead_ele_EoverPInv[neta][ncent];

        // TH1D* h_ele_Zmass[neta][ncent];    

        const Int_t net_bins = 9;
        Double_t et_edges[net_bins+1] = {20.0, 25, 30, 35, 40, 45, 50, 60, 80, 120}; 

        const Int_t iso_bins = 12;
        Double_t isolation_edges[iso_bins+1] = {-20,-10,-5,-3,-2,-1,1,2,3,5,10,15,30};

        for (std::size_t ieta = 0; ieta < neta; ++ieta){
            for (std::size_t icent = 0; icent < ncent; ++icent){

                // Histograms

                h_lead_pho_et[ieta][icent] = new TH1D(Form("h_lead_pho_et_%zu_%zu",ieta,icent),Form("h_lead_pho_et_%zu_%zu;Leading #gamma p_{T} (GeV);Norm. Events",ieta,icent),net_bins,et_edges); // 30,0,150
                if(ieta==0){
                    h_lead_pho_eta[ieta][icent] = new TH1D(Form("h_lead_pho_eta_%zu_%zu",ieta,icent),Form("h_lead_pho_eta_%zu_%zu;Leading #gamma SC |#eta|;Norm. Events",ieta,icent),10,0,1.44);
                    h_lead_pho_SIEIE[ieta][icent] = new TH1D(Form("h_lead_pho_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_pho_SIEIE_%zu_%zu;Leading #gamma #sigma_{#eta #eta};Norm. Events",ieta,icent),15,0,0.03);
                }
                if(ieta==1){
                    h_lead_pho_eta[ieta][icent] = new TH1D(Form("h_lead_pho_eta_%zu_%zu",ieta,icent),Form("h_lead_pho_eta_%zu_%zu;Leading #gamma SC |#eta|;Norm. Events",ieta,icent),10,1.566,2.1);
                    h_lead_pho_SIEIE[ieta][icent] = new TH1D(Form("h_lead_pho_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_pho_SIEIE_%zu_%zu;Leading #gamma #sigma_{#eta #eta};Norm. Events",ieta,icent),30,0,0.1);
                }
                if(ieta==2){
                    h_lead_pho_eta[ieta][icent] = new TH1D(Form("h_lead_pho_eta_%zu_%zu",ieta,icent),Form("h_lead_pho_eta_%zu_%zu;Leading #gamma SC |#eta|;Norm. Events",ieta,icent),10,0,2.1);
                    h_lead_pho_SIEIE[ieta][icent] = new TH1D(Form("h_lead_pho_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_pho_SIEIE_%zu_%zu;Leading #gamma #sigma_{#eta #eta};Norm. Events",ieta,icent),30,0,0.1);
                }                
                h_lead_pho_phi[ieta][icent] = new TH1D(Form("h_lead_pho_phi_%zu_%zu",ieta,icent),Form("h_lead_pho_phi_%zu_%zu;Leading #gamma SC #phi;Norm. Events",ieta,icent),10,-3.14,3.14);
                h_lead_pho_HoverE[ieta][icent] = new TH1D(Form("h_lead_pho_HoverE_%zu_%zu",ieta,icent),Form("h_lead_pho_HoverE_%zu_%zu;Leading #gamma H/E;Norm. Events",ieta,icent),25,0,0.5);
                h_lead_pho_SumCalIso[ieta][icent] = new TH1D(Form("h_lead_pho_SumCalIso_%zu_%zu",ieta,icent),Form("h_lead_pho_SumCalIso_%zu_%zu;Leading #gamma #Sigma Det-Iso;Norm. Events",ieta,icent),30,-20,50);
                h_lead_pho_ECALIso[ieta][icent] = new TH1D(Form("h_lead_pho_ECALIso_%zu_%zu",ieta,icent),Form("h_lead_pho_ECALIso_%zu_%zu;Leading #gamma ECAL Cluster IsoR3;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_HCALIso[ieta][icent] = new TH1D(Form("h_lead_pho_HCALIso_%zu_%zu",ieta,icent),Form("h_lead_pho_HCALIso_%zu_%zu;Leading #gamma HCAL Rechit IsoR3;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_TRKIso[ieta][icent] = new TH1D(Form("h_lead_pho_TRKIso_%zu_%zu",ieta,icent),Form("h_lead_pho_TRKIso_%zu_%zu;Leading #gamma trackIsoR3PtCut20;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_PFPIso[ieta][icent] = new TH1D(Form("h_lead_pho_PFPIso_%zu_%zu",ieta,icent),Form("h_lead_pho_PFPIso_%zu_%zu;Leading #gamma pfpIso3subUE;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_PFNIso[ieta][icent] = new TH1D(Form("h_lead_pho_PFNIso_%zu_%zu",ieta,icent),Form("h_lead_pho_PFNIso_%zu_%zu;Leading #gamma pfnIso3subUE;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_PFCIso[ieta][icent] = new TH1D(Form("h_lead_pho_PFCIso_%zu_%zu",ieta,icent),Form("h_lead_pho_PFCIso_%zu_%zu;Leading #gamma pfcIso3subUE;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_R9[ieta][icent] = new TH1D(Form("h_lead_pho_R9_%zu_%zu",ieta,icent),Form("h_lead_pho_R9_%zu_%zu;Leading #gamma R9;Norm. Events",ieta,icent),40,0,1.1);

                h_lead_ele_Pt[ieta][icent] = new TH1D(Form("h_lead_ele_Pt_%zu_%zu",ieta,icent),Form("h_lead_ele_Pt_%zu_%zu;Leading Electron p_{T} (GeV);Norm. Events",ieta,icent),net_bins,et_edges); // 30,0,150
                h_lead_ele_eta[ieta][icent] = new TH1D(Form("h_lead_ele_eta_%zu_%zu",ieta,icent),Form("h_lead_ele_eta_%zu_%zu;Leading Electron SC #eta;Norm. Events",ieta,icent),10,-2.4,2.4);
                h_lead_ele_phi[ieta][icent] = new TH1D(Form("h_lead_ele_phi_%zu_%zu",ieta,icent),Form("h_lead_ele_phi_%zu_%zu;Leading Electron SC #phi;Norm. Events",ieta,icent),10,-3.14,3.14);
                h_lead_ele_HoverE[ieta][icent] = new TH1D(Form("h_lead_ele_HoverE_%zu_%zu",ieta,icent),Form("h_lead_ele_HoverE_%zu_%zu;Leading Electron H/E;Norm. Events",ieta,icent),20,0,0.1);
                h_lead_ele_SIEIE[ieta][icent] = new TH1D(Form("h_lead_ele_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_ele_SIEIE_%zu_%zu;Leading Electron #sigma_{#eta #eta};Norm. Events",ieta,icent),30,0,0.1);
                h_lead_ele_EoverP[ieta][icent] = new TH1D(Form("h_lead_ele_EoverP_%zu_%zu",ieta,icent),Form("h_lead_ele_EoverP_%zu_%zu;Leading Electron E/P;Norm. Events",ieta,icent),20,0,20);
                h_lead_ele_dEtaSeedAtVtx[ieta][icent] = new TH1D(Form("h_lead_ele_dEtaSeedAtVtx_%zu_%zu",ieta,icent),Form("h_lead_ele_dEtaSeedAtVtx_%zu_%zu;Leading Electron #Delta #eta_{seed, Track};Norm. Events",ieta,icent),20,-0.02,0.02);
                h_lead_ele_dPhiAtVtx[ieta][icent] = new TH1D(Form("h_lead_ele_dPhiAtVtx_%zu_%zu",ieta,icent),Form("h_lead_ele_dPhiAtVtx_%zu_%zu;Leading Electron #Delta #phi_{SC, Track};Norm. Events",ieta,icent),20,-0.2,0.2);
                h_lead_ele_EoverPInv[ieta][icent] = new TH1D(Form("h_lead_ele_EoverPInv_%zu_%zu",ieta,icent),Form("h_lead_ele_EoverPInv_%zu_%zu;Leading Electron EoverPInv;Norm. Events",ieta,icent),20,-0.2,0.2);

                h_sublead_ele_Pt[ieta][icent] = new TH1D(Form("h_sublead_ele_Pt_%zu_%zu",ieta,icent),Form("h_sublead_ele_Pt_%zu_%zu;Sub-Leading Electron p_{T} (GeV);Norm. Events",ieta,icent),net_bins,et_edges); // 30,0,150
                h_sublead_ele_eta[ieta][icent] = new TH1D(Form("h_sublead_ele_eta_%zu_%zu",ieta,icent),Form("h_sublead_ele_eta_%zu_%zu;Sub-Leading Electron SC #eta;Norm. Events",ieta,icent),10,-2.4,2.4);
                h_sublead_ele_phi[ieta][icent] = new TH1D(Form("h_sublead_ele_phi_%zu_%zu",ieta,icent),Form("h_sublead_ele_phi_%zu_%zu;Sub-Leading Electron SC #phi;Norm. Events",ieta,icent),10,-3.14,3.14);
                h_sublead_ele_HoverE[ieta][icent] = new TH1D(Form("h_sublead_ele_HoverE_%zu_%zu",ieta,icent),Form("h_sublead_ele_HoverE_%zu_%zu;Sub-Leading Electron H/E;Norm. Events",ieta,icent),20,0,0.1);
                h_sublead_ele_SIEIE[ieta][icent] = new TH1D(Form("h_sublead_ele_SIEIE_%zu_%zu",ieta,icent),Form("h_sublead_ele_SIEIE_%zu_%zu;Sub-Leading Electron #sigma_{#eta #eta} ;Norm. Events",ieta,icent),30,0,0.1);
                h_sublead_ele_EoverP[ieta][icent] = new TH1D(Form("h_sublead_ele_EoverP_%zu_%zu",ieta,icent),Form("h_sublead_ele_EoverP_%zu_%zu;Sub-Leading Electron E/P;Norm. Events",ieta,icent),20,0,20);
                h_sublead_ele_dEtaSeedAtVtx[ieta][icent] = new TH1D(Form("h_sublead_ele_dEtaSeedAtVtx_%zu_%zu",ieta,icent),Form("h_sublead_ele_dEtaSeedAtVtx_%zu_%zu;Sub-Leading Electron #Delta #eta_{seed, Track};Norm. Events",ieta,icent),20,-0.02,0.02);
                h_sublead_ele_dPhiAtVtx[ieta][icent] = new TH1D(Form("h_sublead_ele_dPhiAtVtx_%zu_%zu",ieta,icent),Form("h_sublead_ele_dPhiAtVtx_%zu_%zu;Sub-Leading Electron #Delta #phi_{SC, Track};Norm. Events",ieta,icent),20,-0.2,0.2);
                h_sublead_ele_EoverPInv[ieta][icent] = new TH1D(Form("h_sublead_ele_EoverPInv_%zu_%zu",ieta,icent),Form("h_sublead_ele_EoverPInv_%zu_%zu;Sub-Leading Electron EoverPInv;Norm. Events",ieta,icent),20,-0.2,0.2);

                // h_ele_Zmass[ieta][icent] = new TH1D(Form("h_ele_Zmass_%zu_%zu",ieta,icent),Form("h_ele_Zmass_%zu_%zu;Invariant Mass;Norm. Events",ieta,icent),60,0,180);
            }
        }
    // ----------------------------------------------------------------------------------------------------------------
    //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
    // ----------------------------------------------------------------------------------------------------------------

    Int_t nEv=photonTree.GetEntries();
    std::cout<<in_label<<" \t => "<<nEv<<"\n";

    for(int iEntry=0; iEntry< nEv; iEntry++){
        displayProgress(iEntry,nEv);
        photonTree.GetEntry(iEntry);
        HiTree.GetEntry(iEntry);
        hltanalysis.GetEntry(iEntry);
        skimanalysis.GetEntry(iEntry);
        // std::cout<<"Inside event loop = "<<iEntry;
        // printf("\n");
        float scale = 1;
        if(!isData) scale*=weight*Ncoll[hiBin];

        int pho_index=-1;
        float Etmax=-1;
        // if(pprimaryVertexFilter)
        // std::cout<<"Event = "<<iEntry<<"    photon = "<<phoEt->size()<<"\t pprimaryVertexFilter ="<<pprimaryVertexFilter<<"\n";
        // if(pclusterCompatibilityFilter)
        // std::cout<<"Event = "<<iEntry<<"    photon = "<<phoEt->size()<<"\t pclusterCompatibilityFilter ="<<pclusterCompatibilityFilter<<"\n";
        // if(pphfCoincFilter2Th4)
        // std::cout<<"Event = "<<iEntry<<"    photon = "<<phoEt->size()<<"\t pphfCoincFilter2Th4 ="<<pphfCoincFilter2Th4<<"\n";

        if(isData && pprimaryVertexFilter<=0) continue;
        if(isData && pclusterCompatibilityFilter<=0) continue;
        if(isData && pphfCoincFilter2Th4<=0) continue;
        if(isData && hiHF>6000) continue;
        
        // Leading Photon
        for(int ipho=0; ipho<(int)phoEt->size(); ipho++){
            if(!isData){
                if(pho_genMatchedIndex->at(ipho)==-1) continue;
                if(mcPID->at(pho_genMatchedIndex->at(ipho)) != 22) continue;
                if(!(abs(mcMomPID->at(pho_genMatchedIndex->at(ipho))) <= 22 || mcMomPID->at(pho_genMatchedIndex->at(ipho)) == -999) ) continue;
                if(!(mcCalIsoDR04->at(pho_genMatchedIndex->at(ipho)) < 5)) continue;
            }
            if(Etmax<phoEt->at(ipho)){
                Etmax = phoEt->at(ipho);
                pho_index = ipho;
            }
        }

        // if((int)phoEt->size()>0) std::cout<<"Event = "<<iEntry<<"    photon index = "<<pho_index<<std::endl;

        // Photon Histograms
        if(pho_index!=-1){
            for (std::size_t ieta = 0; ieta < neta; ++ieta){
                if(!(fabs(phoSCEta->at(pho_index))>=min_eta[ieta] && fabs(phoSCEta->at(pho_index))<max_eta[ieta])) continue;   

                for (std::size_t icent = 0; icent < ncent; ++icent){
                    if(!(hiBin>=min_cent[icent] && hiBin<max_cent[icent])) continue;  

                    // std::cout<<iEntry<<"  Inside Photon Histogram loop \n";
                    
                    if(phoEt->at(pho_index)<40) continue;
                    if(phoSigmaEtaEta_2012->at(pho_index)<=0.002) continue;        // Spike Cuts
                    if(fabs(pho_seedTime->at(pho_index))>=3 || pho_swissCrx->at(pho_index)>=0.9) continue;

                    // if(!(HLT_HIGEDPhoton20_v1 || HLT_HIGEDPhoton30_v1 || HLT_HIGEDPhoton40_v1)) continue;
                    if(!(HLT_HIGEDPhoton40_v1)) continue;

                    double isolation=pho_ecalClusterIsoR3->at(pho_index)+pho_hcalRechitIsoR3->at(pho_index)+pho_trackIsoR3PtCut20->at(pho_index);
                    // if(!(phoHoverE->at(pho_index)<=0.2)) continue;
                    // if(!(isolation<5)) continue;

                    // Loose Photon ID
                    if(!(phoHoverE->at(pho_index)<=0.247995)) continue;
                    if(!(phoSigmaIEtaIEta_2012->at(pho_index)<=0.012186)) continue;
                    if(!(isolation<11.697505)) continue;

                    h_lead_pho_et[ieta][icent]->Fill(         phoEt->at(pho_index), scale);
                    h_lead_pho_eta[ieta][icent]->Fill(        phoSCEta->at(pho_index), scale);
                    h_lead_pho_phi[ieta][icent]->Fill(        phoSCPhi->at(pho_index), scale);
                    h_lead_pho_HoverE[ieta][icent]->Fill(     phoHoverE->at(pho_index), scale);
                    h_lead_pho_SIEIE[ieta][icent]->Fill(      phoSigmaIEtaIEta_2012->at(pho_index), scale);
                    h_lead_pho_SumCalIso[ieta][icent]->Fill(  pho_ecalClusterIsoR3->at(pho_index)+pho_hcalRechitIsoR3->at(pho_index)+pho_trackIsoR3PtCut20->at(pho_index), scale);
                    h_lead_pho_ECALIso[ieta][icent]->Fill(    pho_ecalClusterIsoR3->at(pho_index), scale);
                    h_lead_pho_HCALIso[ieta][icent]->Fill(    pho_hcalRechitIsoR3->at(pho_index), scale);
                    h_lead_pho_TRKIso[ieta][icent]->Fill(     pho_trackIsoR3PtCut20->at(pho_index), scale);
                    h_lead_pho_PFPIso[ieta][icent]->Fill(     pfpIso3subUE->at(pho_index), scale);
                    h_lead_pho_PFNIso[ieta][icent]->Fill(     pfnIso3subUE->at(pho_index), scale);
                    h_lead_pho_PFCIso[ieta][icent]->Fill(     pfcIso3subUE->at(pho_index), scale);
                    h_lead_pho_R9[ieta][icent]->Fill(         phoR9_2012->at(pho_index), scale);
                }
            }
        }

        int lead_ele_index=-1,sublead_ele_index=-1;
        float ele_Etmax=-1,ele_Etsub=-1;

        // Electron Gen Matching for MC
        std::vector<int> ele_genMatchedIndex;
        if(!isData){
            for(int iele=0;iele<(int)elePt->size();iele++){
                int matchedIndex = -1;
                float minDR = 0.3;
                for (unsigned igen = 0; igen < mcEt->size(); ++igen) {
                    if ( abs(mcPID->at(igen)) != 11) 
                    continue;                    
                    if (dr(eleEta->at(iele), elePhi->at(iele), mcEta->at(igen), mcPhi->at(igen)) < minDR){ // delta2 && mcPt->at(igen) > currentMaxPt) {
                    minDR = dr(eleEta->at(iele), elePhi->at(iele), mcEta->at(igen), mcPhi->at(igen));
                    matchedIndex = igen;
                    }
                }
                ele_genMatchedIndex.push_back(matchedIndex);
            }
        }
        // Leading Electron
        for(int iele=0;iele<(int)elePt->size();iele++){
            if(!isData){
                if(ele_genMatchedIndex.at(iele)==-1) continue;
                if(abs(mcPID->at(ele_genMatchedIndex.at(iele))) != 11) continue;       // True electron
                if(!(abs(mcMomPID->at(ele_genMatchedIndex.at(iele))) == 23)) continue; // Parent Z
            }
            if(ele_Etmax<elePt->at(iele)){
                ele_Etmax = elePt->at(iele);
                lead_ele_index = iele;
            }
        }
        // Sub Leading Electron
        for(int iele=0;iele<(int)elePt->size();iele++){
            if(iele==lead_ele_index) continue;
            if(!isData){
                if(ele_genMatchedIndex.at(iele)==-1) continue;
                if(abs(mcPID->at(ele_genMatchedIndex.at(iele))) != 11) continue;       // True electron
                if(!(abs(mcMomPID->at(ele_genMatchedIndex.at(iele))) == 23)) continue; // Parent Z
            }
            if(ele_Etsub<elePt->at(iele)){
                ele_Etsub = elePt->at(iele);
                sublead_ele_index = iele;
            }
        }

        if(lead_ele_index!=-1){
            for (std::size_t ieta = 0; ieta < neta; ++ieta){
                if(!(fabs(eleEta->at(lead_ele_index))>=min_eta[ieta] && fabs(eleEta->at(lead_ele_index))<max_eta[ieta])) continue;   

                for (std::size_t icent = 0; icent < ncent; ++icent){
                    if(!(hiBin>=min_cent[icent] && hiBin<max_cent[icent])) continue;  

                    // if(!(HLT_HIEle20Gsf_v1 || HLT_HIEle30Gsf_v1 || HLT_HIEle40Gsf_v1 || HLT_HIEle50Gsf_v1)) continue;
                    if(!(HLT_HIEle40Gsf_v1)) continue;

                    if(elePt->at(lead_ele_index)<40) continue;
                    Bool_t flagLooseEle = false;
                    if(fabs(eleEta->at(lead_ele_index))<1.442) {
                        if(eleMissHits->at(lead_ele_index)>1) continue;
                        if(eleIP3D->at(lead_ele_index)>=0.03) continue;
                    }
                    // Electron Veto ID
                    if((hiBin)>=0 && (hiBin)<60){
                        if(eleSigmaIEtaIEta_2012->at(lead_ele_index)>=0.0147) continue;
                        if(eledEtaSeedAtVtx->at(lead_ele_index)>=0.0041) continue;
                        if(eledPhiAtVtx->at(lead_ele_index)>=0.0853) continue;
                        if(eleHoverE->at(lead_ele_index)>=0.2733) continue;
                        if(eleEoverPInv->at(lead_ele_index)>=0.0367) continue;
                        flagLooseEle=true;
                    }
                    else{
                        if(eleSigmaIEtaIEta_2012->at(lead_ele_index)>=0.0113) continue;
                        if(eledEtaSeedAtVtx->at(lead_ele_index)>=0.0037) continue;
                        if(eledPhiAtVtx->at(lead_ele_index)>=0.1280) continue;
                        if(eleHoverE->at(lead_ele_index)>=0.1814) continue;
                        if(eleEoverPInv->at(lead_ele_index)>=0.1065) continue;
                        flagLooseEle=true;
                    }

                    if(!flagLooseEle) continue;

                    h_lead_ele_Pt[ieta][icent]->Fill(                elePt->at(lead_ele_index) ,scale);
                    h_lead_ele_eta[ieta][icent]->Fill(               eleEta->at(lead_ele_index) ,scale);
                    h_lead_ele_SIEIE[ieta][icent]->Fill(             eleSigmaIEtaIEta_2012->at(lead_ele_index) ,scale);
                    h_lead_ele_phi[ieta][icent]->Fill(               elePhi->at(lead_ele_index) ,scale);
                    h_lead_ele_HoverE[ieta][icent]->Fill(            eleHoverE->at(lead_ele_index) ,scale);
                    h_lead_ele_EoverP[ieta][icent]->Fill(            eleEoverP->at(lead_ele_index) ,scale);
                    h_lead_ele_dEtaSeedAtVtx[ieta][icent]->Fill(     eledEtaSeedAtVtx->at(lead_ele_index) ,scale);
                    h_lead_ele_dPhiAtVtx[ieta][icent]->Fill(         eledPhiAtVtx->at(lead_ele_index) ,scale);
                    h_lead_ele_EoverPInv[ieta][icent]->Fill(         eleEoverPInv->at(lead_ele_index) ,scale);
                    
                }
            }
        }

        if(sublead_ele_index!=-1){
            for (std::size_t ieta = 0; ieta < neta; ++ieta){
                if(!(fabs(eleEta->at(sublead_ele_index))>=min_eta[ieta] && fabs(eleEta->at(sublead_ele_index))<max_eta[ieta])) continue;   

                for (std::size_t icent = 0; icent < ncent; ++icent){
                    if(!(hiBin>=min_cent[icent] && hiBin<max_cent[icent])) continue; 
                    
                    if(!(HLT_HIEle15Ele10Gsf_v1 || HLT_HIEle15Ele10GsfMass50_v1 || HLT_HIDoubleEle10Gsf_v1 || HLT_HIDoubleEle10GsfMass50_v1 || HLT_HIDoubleEle15Gsf_v1 || HLT_HIDoubleEle15GsfMass50_v1)) continue;

                    if(elePt->at(sublead_ele_index)<20) continue;
                    Bool_t flagLooseEle = false;
                    if(fabs(eleEta->at(sublead_ele_index))<1.442) {
                        if(eleMissHits->at(sublead_ele_index)>1) continue;
                        if(eleIP3D->at(sublead_ele_index)>=0.03) continue;
                    }
                    // Electron Veto ID
                    if((hiBin)>=0 && (hiBin)<60){
                        if(eleSigmaIEtaIEta_2012->at(sublead_ele_index)>=0.0147) continue;
                        if(eledEtaSeedAtVtx->at(sublead_ele_index)>=0.0041) continue;
                        if(eledPhiAtVtx->at(sublead_ele_index)>=0.0853) continue;
                        if(eleHoverE->at(sublead_ele_index)>=0.2733) continue;
                        if(eleEoverPInv->at(sublead_ele_index)>=0.0367) continue;
                        flagLooseEle=true;
                    }
                    else{
                        if(eleSigmaIEtaIEta_2012->at(sublead_ele_index)>=0.0113) continue;
                        if(eledEtaSeedAtVtx->at(sublead_ele_index)>=0.0037) continue;
                        if(eledPhiAtVtx->at(sublead_ele_index)>=0.1280) continue;
                        if(eleHoverE->at(sublead_ele_index)>=0.1814) continue;
                        if(eleEoverPInv->at(sublead_ele_index)>=0.1065) continue;
                        flagLooseEle=true;
                    }

                    if(!flagLooseEle) continue;

                    h_sublead_ele_Pt[ieta][icent]->Fill(                elePt->at(sublead_ele_index) ,scale);
                    h_sublead_ele_eta[ieta][icent]->Fill(               eleEta->at(sublead_ele_index) ,scale);
                    h_sublead_ele_SIEIE[ieta][icent]->Fill(             eleSigmaIEtaIEta_2012->at(sublead_ele_index) ,scale);
                    h_sublead_ele_phi[ieta][icent]->Fill(               elePhi->at(sublead_ele_index) ,scale);
                    h_sublead_ele_HoverE[ieta][icent]->Fill(            eleHoverE->at(sublead_ele_index) ,scale);
                    h_sublead_ele_EoverP[ieta][icent]->Fill(            eleEoverP->at(sublead_ele_index) ,scale);
                    h_sublead_ele_dEtaSeedAtVtx[ieta][icent]->Fill(     eledEtaSeedAtVtx->at(sublead_ele_index) ,scale);
                    h_sublead_ele_dPhiAtVtx[ieta][icent]->Fill(         eledPhiAtVtx->at(sublead_ele_index) ,scale);
                    h_sublead_ele_EoverPInv[ieta][icent]->Fill(         eleEoverPInv->at(sublead_ele_index) ,scale);
                }
            }
        }        

    } // End event loop

    std::cout<<"\nNumber of selected photons = "<<h_lead_pho_et[0][0]->GetEntries()<<"\n";

    if(h_lead_pho_et[0][0]->GetEntries()==0){ std::cout<<"Return from function 0 entries => Not writing to file\n"; return;}
    
    // ------ Save histograms to File
    fout_2018->cd();
    gDirectory->mkdir(in_label);
    fout_2018->cd(in_label);

    for (std::size_t ieta = 0; ieta < neta; ++ieta){
        for (std::size_t icent = 0; icent < ncent; ++icent){
            
            h_lead_pho_et[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_eta[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_phi[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_HoverE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_SIEIE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_SumCalIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_ECALIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_HCALIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_TRKIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_PFPIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_PFNIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_PFCIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_R9[ieta][icent]->Write("",TObject::kOverwrite);

            h_lead_ele_Pt[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_eta[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_SIEIE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_phi[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_HoverE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_EoverP[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_dEtaSeedAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_dPhiAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_EoverPInv[ieta][icent]->Write("",TObject::kOverwrite);

            h_sublead_ele_Pt[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_eta[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_SIEIE[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_phi[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_HoverE[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_EoverP[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_dEtaSeedAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_dPhiAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_EoverPInv[ieta][icent]->Write("",TObject::kOverwrite);
        }
    }

    fout_2018->cd();
}

void Data_fill_pho_ele_2023(std::vector<TString> in_file_path,TString in_label, Bool_t isData){
    // in_file_path = Directory name for the specified Data Run
    // in_label = Corresponding legend label

    TChain photonTree("photonTree"), HiTree("HiTree"),skimanalysis("skimanalysis"),hltanalysis("hltanalysis");
    
    for(TString temp_path:in_file_path){
        // std:cout<<"file path = "<<temp_path+input_pho_tree<<"\n";

        photonTree.Add(temp_path+"/ggHiNtuplizer/EventTree");
        HiTree.Add(temp_path+"/hiEvtAnalyzer/HiTree");
        skimanalysis.Add(temp_path+"/skimanalysis/HltTree");
        hltanalysis.Add(temp_path+"/hltanalysis/HltTree");
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Variables

        Int_t hiBin=-1;
        Float_t hiHF = -1; 
        Float_t weight = 1; 
        Float_t pthat=0;

        HiTree.LoadTree(0);

        HiTree.SetBranchAddress("hiBin", &hiBin);
        HiTree.SetBranchAddress("hiHF", &hiHF);
        if(HiTree.GetBranch("weight") && !isData)
            HiTree.SetBranchAddress("weight", &weight);
        if(HiTree.GetBranch("pthat") && !isData)
            HiTree.SetBranchAddress("pthat", &pthat);

        // reco::GenParticle
        std::vector<int>*    mcStatus = 0;
        std::vector<int>*    mcPID = 0;
        std::vector<int>*    mcMomPID = 0;
        std::vector<float>*  mcPt = 0;
        std::vector<float>*  mcEta = 0;
        std::vector<float>*  mcPhi = 0;
        std::vector<float>*  mcE = 0;
        std::vector<float>*  mcEt = 0;
        std::vector<float>*  mcMass = 0;
        std::vector<float>*  mcCalIsoDR03 = 0;
        std::vector<float>*  mcCalIsoDR04 = 0;
        std::vector<float>*  mcTrkIsoDR03 = 0;
        std::vector<float>*  mcTrkIsoDR04 = 0;

        std::vector<int>*    pho_genMatchedIndex = 0;

        // reco::Photon
        std::vector<float>*  phoE = 0;
        std::vector<float>*  phoEt = 0;
        std::vector<float>*  phoEta = 0;
        std::vector<float>*  phoPhi = 0;
        std::vector<float>*  phoHoverE = 0;
        std::vector<float>*  pho_ecalClusterIsoR3 = 0;
        std::vector<float>*  pho_hcalRechitIsoR3 = 0;
        std::vector<float>*  pho_trackIsoR3PtCut20 = 0;
        std::vector<float>*  phoSigmaEtaEta_2012 = 0;

        std::vector<float>*  pfcIso3subUEec = 0;
        std::vector<float>*  pfnIso3subUEec = 0;
        std::vector<float>*  pfpIso3subUEec = 0; //! For MiniAOD by default this does not exclude the cone!!

        std::vector<float>*  pfpIso3subSCsubUE = 0; //! For MiniAOD by default this does not exclude the cone!!

        std::vector<float>*  pho_swissCrx = 0;
        std::vector<float>*  pho_seedTime = 0;

        std::vector<float>* phoSCRawE = 0;
        std::vector<float>* phoSCEta = 0;
        std::vector<float>* phoSCPhi = 0;
        std::vector<float>* phoSCEtaWidth = 0;
        std::vector<float>* phoSCPhiWidth = 0;
        std::vector<float>* phoE3x3_2012 = 0;
        std::vector<float>* phoMaxEnergyXtal_2012 = 0;
        std::vector<float>* phoE2nd_2012 = 0;
        std::vector<float>* phoELeft_2012 = 0;
        std::vector<float>* phoERight_2012 = 0;
        std::vector<float>* phoETop_2012 = 0;
        std::vector<float>* phoEBottom_2012 = 0;
        std::vector<float>* phoSigmaIEtaIEta_2012 = 0;
        std::vector<float>* phoSigmaIEtaIPhi_2012 = 0;
        std::vector<float>* phoSigmaIPhiIPhi_2012 = 0;
        Float_t rho = 0;
        std::vector<float>* phoESEn = 0;

        std::vector<int>* phoSCnHits = 0;
        std::vector<float>* phoE5x5_2012 = 0;
        std::vector<float>* phoHadTowerOverEm = 0;
        std::vector<float>* phoHadTowerOverEm1 = 0;
        std::vector<float>* phoHadTowerOverEm2 = 0;
        std::vector<float>* phoR9_2012 = 0;

        // reco::Electron
        std::vector<float>* elePt = 0;
        std::vector<float>* eleEta = 0;
        std::vector<float>* elePhi = 0;
        std::vector<float>* eleEoverP = 0;

        std::vector<float>* eleSigmaIEtaIEta_2012 = 0;
        std::vector<float>* eledEtaSeedAtVtx = 0;
        std::vector<float>* eledPhiAtVtx = 0;
        std::vector<float>* eleHoverE = 0;
        std::vector<float>* eleEoverPInv = 0;
        std::vector<int>* eleMissHits = 0;
        std::vector<int>* eleCharge = 0;
        std::vector<float>* eleIP3D = 0;
        std::vector<float>* eleEn = 0;

        if(photonTree.GetBranch("mcPID") && !isData){

            photonTree.SetBranchAddress("mcStatus",     &mcStatus);
            photonTree.SetBranchAddress("mcPID",        &mcPID);
            photonTree.SetBranchAddress("mcPt",         &mcPt);
            photonTree.SetBranchAddress("mcEta",        &mcEta);
            photonTree.SetBranchAddress("mcPhi",        &mcPhi);
            photonTree.SetBranchAddress("mcE",          &mcE);
            photonTree.SetBranchAddress("mcEt",         &mcEt);   
            photonTree.SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03);
            photonTree.SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04);
            photonTree.SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03);
            photonTree.SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04);
            photonTree.SetBranchAddress("mcMass",       &mcMass);
            photonTree.SetBranchAddress("mcMomPID",     &mcMomPID);
            
            photonTree.SetBranchAddress("pho_genMatchedIndex",   &pho_genMatchedIndex);
        }

        photonTree.SetBranchAddress("phoE",                    &phoE);
        photonTree.SetBranchAddress("phoEt",                   &phoEt);
        photonTree.SetBranchAddress("phoEta",                  &phoEta);
        photonTree.SetBranchAddress("phoPhi",                  &phoPhi);
        photonTree.SetBranchAddress("phoHoverE",               &phoHoverE);
        photonTree.SetBranchAddress("pho_ecalClusterIsoR3",    &pho_ecalClusterIsoR3);
        photonTree.SetBranchAddress("pho_hcalRechitIsoR3",     &pho_hcalRechitIsoR3);
        photonTree.SetBranchAddress("pho_trackIsoR3PtCut20",   &pho_trackIsoR3PtCut20);
        photonTree.SetBranchAddress("phoSigmaEtaEta_2012",     &phoSigmaEtaEta_2012);
        photonTree.SetBranchAddress("pfcIso3subUEec",          &pfcIso3subUEec);
        photonTree.SetBranchAddress("pfnIso3subUEec",          &pfnIso3subUEec);
        photonTree.SetBranchAddress("pfpIso3subUEec",          &pfpIso3subUEec);

        photonTree.SetBranchAddress("pfpIso3subSCsubUE",       &pfpIso3subSCsubUE);

        photonTree.SetBranchAddress("phoSCRawE",               &phoSCRawE);
        photonTree.SetBranchAddress("phoSCEta",                &phoSCEta);
        photonTree.SetBranchAddress("phoSCPhi",                &phoSCPhi);
        photonTree.SetBranchAddress("phoSCEtaWidth",           &phoSCEtaWidth);
        photonTree.SetBranchAddress("phoSCPhiWidth",           &phoSCPhiWidth);
        photonTree.SetBranchAddress("phoE3x3_2012",            &phoE3x3_2012);
        photonTree.SetBranchAddress("phoMaxEnergyXtal_2012",   &phoMaxEnergyXtal_2012);
        photonTree.SetBranchAddress("phoE2nd_2012",            &phoE2nd_2012);
        photonTree.SetBranchAddress("phoELeft_2012",           &phoELeft_2012);
        photonTree.SetBranchAddress("phoERight_2012",          &phoERight_2012);
        photonTree.SetBranchAddress("phoETop_2012",            &phoETop_2012);
        photonTree.SetBranchAddress("phoEBottom_2012",         &phoEBottom_2012);
        photonTree.SetBranchAddress("phoSigmaIEtaIEta_2012",   &phoSigmaIEtaIEta_2012);
        photonTree.SetBranchAddress("phoSigmaIEtaIPhi_2012",   &phoSigmaIEtaIPhi_2012);
        photonTree.SetBranchAddress("phoSigmaIPhiIPhi_2012",   &phoSigmaIPhiIPhi_2012);
        photonTree.SetBranchAddress("rho",                     &rho);
        photonTree.SetBranchAddress("phoESEn",                 &phoESEn);

        photonTree.SetBranchAddress("phoSCnHits",              &phoSCnHits);
        photonTree.SetBranchAddress("phoE5x5_2012",            &phoE5x5_2012);
        photonTree.SetBranchAddress("phoHadTowerOverEm",       &phoHadTowerOverEm);
        photonTree.SetBranchAddress("phoHadTowerOverEm1",      &phoHadTowerOverEm1);
        photonTree.SetBranchAddress("phoHadTowerOverEm2",      &phoHadTowerOverEm2);
        photonTree.SetBranchAddress("phoR9_2012",              &phoR9_2012);

        photonTree.SetBranchAddress("pho_swissCrx",            &pho_swissCrx);
        photonTree.SetBranchAddress("pho_seedTime",            &pho_seedTime);

        photonTree.SetBranchAddress("elePt",                   &elePt);
        photonTree.SetBranchAddress("eleEta",                  &eleEta);
        photonTree.SetBranchAddress("elePhi",                  &elePhi);
        photonTree.SetBranchAddress("eleEoverP",               &eleEoverP);

        photonTree.SetBranchAddress("eleSigmaIEtaIEta_2012",   &eleSigmaIEtaIEta_2012);
        photonTree.SetBranchAddress("eledEtaSeedAtVtx",        &eledEtaSeedAtVtx);
        photonTree.SetBranchAddress("eledPhiAtVtx",            &eledPhiAtVtx);
        photonTree.SetBranchAddress("eleHoverE",               &eleHoverE);
        photonTree.SetBranchAddress("eleEoverPInv",            &eleEoverPInv);
        photonTree.SetBranchAddress("eleMissHits",             &eleMissHits);
        photonTree.SetBranchAddress("eleCharge",               &eleCharge);
        photonTree.SetBranchAddress("eleIP3D",                 &eleIP3D);
        photonTree.SetBranchAddress("eleEn",                   &eleEn);

        int pprimaryVertexFilter = 1;           // Pass if it doesn't exist in the Forest
        int pclusterCompatibilityFilter = 1;
        int pphfCoincFilter2Th4 = 1;            // Different name from AOD

        if(skimanalysis.GetBranch("pprimaryVertexFilter")) skimanalysis.SetBranchAddress("pprimaryVertexFilter",  &pprimaryVertexFilter);
        if(skimanalysis.GetBranch("pclusterCompatibilityFilter")) skimanalysis.SetBranchAddress("pclusterCompatibilityFilter",  &pclusterCompatibilityFilter);
        if(skimanalysis.GetBranch("pphfCoincFilter2Th4")) skimanalysis.SetBranchAddress("pphfCoincFilter2Th4",  &pphfCoincFilter2Th4);

        int HLT_HIGEDPhoton20_v10        = 1;       // Pass if it doesn't exist in the Forest
        int HLT_HIGEDPhoton30_v10        = 1;
        int HLT_HIGEDPhoton40_v10        = 1;

        int HLT_HIEle20Gsf_v10           = 1;
        int HLT_HIEle30Gsf_v10           = 1;
        int HLT_HIEle40Gsf_v10           = 1;
        int HLT_HIEle50Gsf_v10           = 1;
        int HLT_HIEle15Ele10Gsf_v10      = 1;
        int HLT_HIEle15Ele10GsfMass50_v10= 1;   

        int HLT_HIDoubleEle10Gsf_v10         = 1;
        int HLT_HIDoubleEle10GsfMass50_v10   = 1;
        int HLT_HIDoubleEle15Gsf_v10         = 1;
        int HLT_HIDoubleEle15GsfMass50_v10   = 1;  

        if(hltanalysis.GetBranch("HLT_HIGEDPhoton40_v10")) std::cout<<"HLT_HIGEDPhoton40_v10 branch doesn't exist\n"; 

        hltanalysis.SetBranchAddress("HLT_HIGEDPhoton20_v10",           &HLT_HIGEDPhoton20_v10);
        hltanalysis.SetBranchAddress("HLT_HIGEDPhoton30_v10",           &HLT_HIGEDPhoton30_v10);
        hltanalysis.SetBranchAddress("HLT_HIGEDPhoton40_v10",           &HLT_HIGEDPhoton40_v10);

        hltanalysis.SetBranchAddress("HLT_HIEle20Gsf_v10",              &HLT_HIEle20Gsf_v10);
        hltanalysis.SetBranchAddress("HLT_HIEle30Gsf_v10",              &HLT_HIEle30Gsf_v10);
        hltanalysis.SetBranchAddress("HLT_HIEle40Gsf_v10",              &HLT_HIEle40Gsf_v10);
        hltanalysis.SetBranchAddress("HLT_HIEle50Gsf_v10",              &HLT_HIEle50Gsf_v10);
        hltanalysis.SetBranchAddress("HLT_HIEle15Ele10Gsf_v10",         &HLT_HIEle15Ele10Gsf_v10);
        hltanalysis.SetBranchAddress("HLT_HIEle15Ele10GsfMass50_v10",   &HLT_HIEle15Ele10GsfMass50_v10); 

        hltanalysis.SetBranchAddress("HLT_HIDoubleEle10Gsf_v10",        &HLT_HIDoubleEle10Gsf_v10);
        hltanalysis.SetBranchAddress("HLT_HIDoubleEle10GsfMass50_v10",  &HLT_HIDoubleEle10GsfMass50_v10);
        hltanalysis.SetBranchAddress("HLT_HIDoubleEle15Gsf_v10",        &HLT_HIDoubleEle15Gsf_v10);
        hltanalysis.SetBranchAddress("HLT_HIDoubleEle15GsfMass50_v10",  &HLT_HIDoubleEle15GsfMass50_v10);

    // -------- End Tree Variable Declaration
    // ----------------------------------------------------------------------------------------------------------------
    // Histograms
        TH1::SetDefaultSumw2();
        TH2::SetDefaultSumw2();

        TH1D* h_lead_pho_et[neta][ncent];
        TH1D* h_lead_pho_eta[neta][ncent];
        TH1D* h_lead_pho_SIEIE[neta][ncent];
        TH1D* h_lead_pho_phi[neta][ncent];
        TH1D* h_lead_pho_HoverE[neta][ncent];
        TH1D* h_lead_pho_SumCalIso[neta][ncent];
        TH1D* h_lead_pho_ECALIso[neta][ncent];
        TH1D* h_lead_pho_HCALIso[neta][ncent];
        TH1D* h_lead_pho_TRKIso[neta][ncent];
        TH1D* h_lead_pho_PFPIso[neta][ncent];
        TH1D* h_lead_pho_PFNIso[neta][ncent];
        TH1D* h_lead_pho_PFCIso[neta][ncent];
        TH1D* h_lead_pho_R9[neta][ncent];

        TH1D* h_lead_ele_Pt[neta][ncent];
        TH1D* h_lead_ele_eta[neta][ncent];
        TH1D* h_lead_ele_SIEIE[neta][ncent];
        TH1D* h_lead_ele_phi[neta][ncent];
        TH1D* h_lead_ele_HoverE[neta][ncent];
        TH1D* h_lead_ele_EoverP[neta][ncent];
        TH1D* h_lead_ele_dEtaSeedAtVtx[neta][ncent];
        TH1D* h_lead_ele_dPhiAtVtx[neta][ncent];
        TH1D* h_lead_ele_EoverPInv[neta][ncent];

        TH1D* h_sublead_ele_Pt[neta][ncent];
        TH1D* h_sublead_ele_eta[neta][ncent];
        TH1D* h_sublead_ele_SIEIE[neta][ncent];
        TH1D* h_sublead_ele_phi[neta][ncent];
        TH1D* h_sublead_ele_HoverE[neta][ncent];
        TH1D* h_sublead_ele_EoverP[neta][ncent];
        TH1D* h_sublead_ele_dEtaSeedAtVtx[neta][ncent];
        TH1D* h_sublead_ele_dPhiAtVtx[neta][ncent];
        TH1D* h_sublead_ele_EoverPInv[neta][ncent];

        // TH1D* h_ele_Zmass[neta][ncent];    

        const Int_t net_bins = 9;
        Double_t et_edges[net_bins+1] = {20.0, 25, 30, 35, 40, 45, 50, 60, 80, 120}; 

        const Int_t iso_bins = 12;
        Double_t isolation_edges[iso_bins+1] = {-20,-10,-5,-3,-2,-1,1,2,3,5,10,15,30};

        for (std::size_t ieta = 0; ieta < neta; ++ieta){
            for (std::size_t icent = 0; icent < ncent; ++icent){
                // Histograms

                h_lead_pho_et[ieta][icent] = new TH1D(Form("h_lead_pho_et_%zu_%zu",ieta,icent),Form("h_lead_pho_et_%zu_%zu;Leading #gamma p_{T} (GeV);Norm. Events",ieta,icent),net_bins,et_edges); // 30,0,150
                if(ieta==0){
                    h_lead_pho_eta[ieta][icent] = new TH1D(Form("h_lead_pho_eta_%zu_%zu",ieta,icent),Form("h_lead_pho_eta_%zu_%zu;Leading #gamma SC |#eta|;Norm. Events",ieta,icent),10,0,1.44);
                    h_lead_pho_SIEIE[ieta][icent] = new TH1D(Form("h_lead_pho_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_pho_SIEIE_%zu_%zu;Leading #gamma #sigma_{#eta #eta};Norm. Events",ieta,icent),15,0,0.03);
                }
                if(ieta==1){
                    h_lead_pho_eta[ieta][icent] = new TH1D(Form("h_lead_pho_eta_%zu_%zu",ieta,icent),Form("h_lead_pho_eta_%zu_%zu;Leading #gamma SC |#eta|;Norm. Events",ieta,icent),10,1.566,2.1);
                    h_lead_pho_SIEIE[ieta][icent] = new TH1D(Form("h_lead_pho_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_pho_SIEIE_%zu_%zu;Leading #gamma #sigma_{#eta #eta};Norm. Events",ieta,icent),30,0,0.1);
                }
                if(ieta==2){
                    h_lead_pho_eta[ieta][icent] = new TH1D(Form("h_lead_pho_eta_%zu_%zu",ieta,icent),Form("h_lead_pho_eta_%zu_%zu;Leading #gamma SC |#eta|;Norm. Events",ieta,icent),10,0,2.1);
                    h_lead_pho_SIEIE[ieta][icent] = new TH1D(Form("h_lead_pho_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_pho_SIEIE_%zu_%zu;Leading #gamma #sigma_{#eta #eta};Norm. Events",ieta,icent),30,0,0.1);
                }                
                h_lead_pho_phi[ieta][icent] = new TH1D(Form("h_lead_pho_phi_%zu_%zu",ieta,icent),Form("h_lead_pho_phi_%zu_%zu;Leading #gamma SC #phi;Norm. Events",ieta,icent),10,-3.14,3.14);
                h_lead_pho_HoverE[ieta][icent] = new TH1D(Form("h_lead_pho_HoverE_%zu_%zu",ieta,icent),Form("h_lead_pho_HoverE_%zu_%zu;Leading #gamma H/E;Norm. Events",ieta,icent),25,0,0.5);
                h_lead_pho_SumCalIso[ieta][icent] = new TH1D(Form("h_lead_pho_SumCalIso_%zu_%zu",ieta,icent),Form("h_lead_pho_SumCalIso_%zu_%zu;Leading #gamma #Sigma Det-Iso;Norm. Events",ieta,icent),30,-20,50);
                h_lead_pho_ECALIso[ieta][icent] = new TH1D(Form("h_lead_pho_ECALIso_%zu_%zu",ieta,icent),Form("h_lead_pho_ECALIso_%zu_%zu;Leading #gamma ECAL Cluster IsoR3;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_HCALIso[ieta][icent] = new TH1D(Form("h_lead_pho_HCALIso_%zu_%zu",ieta,icent),Form("h_lead_pho_HCALIso_%zu_%zu;Leading #gamma HCAL Rechit IsoR3;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_TRKIso[ieta][icent] = new TH1D(Form("h_lead_pho_TRKIso_%zu_%zu",ieta,icent),Form("h_lead_pho_TRKIso_%zu_%zu;Leading #gamma trackIsoR3PtCut20;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_PFPIso[ieta][icent] = new TH1D(Form("h_lead_pho_PFPIso_%zu_%zu",ieta,icent),Form("h_lead_pho_PFPIso_%zu_%zu;Leading #gamma pfpIso3subUE;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_PFNIso[ieta][icent] = new TH1D(Form("h_lead_pho_PFNIso_%zu_%zu",ieta,icent),Form("h_lead_pho_PFNIso_%zu_%zu;Leading #gamma pfnIso3subUE;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_PFCIso[ieta][icent] = new TH1D(Form("h_lead_pho_PFCIso_%zu_%zu",ieta,icent),Form("h_lead_pho_PFCIso_%zu_%zu;Leading #gamma pfcIso3subUE;Norm. Events",ieta,icent),iso_bins, isolation_edges);
                h_lead_pho_R9[ieta][icent] = new TH1D(Form("h_lead_pho_R9_%zu_%zu",ieta,icent),Form("h_lead_pho_R9_%zu_%zu;Leading #gamma R9;Norm. Events",ieta,icent),40,0,1.1);

                h_lead_ele_Pt[ieta][icent] = new TH1D(Form("h_lead_ele_Pt_%zu_%zu",ieta,icent),Form("h_lead_ele_Pt_%zu_%zu;Leading Electron p_{T} (GeV);Norm. Events",ieta,icent),net_bins,et_edges); // 30,0,150
                h_lead_ele_eta[ieta][icent] = new TH1D(Form("h_lead_ele_eta_%zu_%zu",ieta,icent),Form("h_lead_ele_eta_%zu_%zu;Leading Electron SC #eta;Norm. Events",ieta,icent),10,-2.4,2.4);
                h_lead_ele_phi[ieta][icent] = new TH1D(Form("h_lead_ele_phi_%zu_%zu",ieta,icent),Form("h_lead_ele_phi_%zu_%zu;Leading Electron SC #phi;Norm. Events",ieta,icent),10,-3.14,3.14);
                h_lead_ele_HoverE[ieta][icent] = new TH1D(Form("h_lead_ele_HoverE_%zu_%zu",ieta,icent),Form("h_lead_ele_HoverE_%zu_%zu;Leading Electron H/E;Norm. Events",ieta,icent),20,0,0.1);
                h_lead_ele_SIEIE[ieta][icent] = new TH1D(Form("h_lead_ele_SIEIE_%zu_%zu",ieta,icent),Form("h_lead_ele_SIEIE_%zu_%zu;Leading Electron #sigma_{#eta #eta};Norm. Events",ieta,icent),30,0,0.1);
                h_lead_ele_EoverP[ieta][icent] = new TH1D(Form("h_lead_ele_EoverP_%zu_%zu",ieta,icent),Form("h_lead_ele_EoverP_%zu_%zu;Leading Electron E/P;Norm. Events",ieta,icent),20,0,20);
                h_lead_ele_dEtaSeedAtVtx[ieta][icent] = new TH1D(Form("h_lead_ele_dEtaSeedAtVtx_%zu_%zu",ieta,icent),Form("h_lead_ele_dEtaSeedAtVtx_%zu_%zu;Leading Electron #Delta #eta_{seed, Track};Norm. Events",ieta,icent),20,-0.02,0.02);
                h_lead_ele_dPhiAtVtx[ieta][icent] = new TH1D(Form("h_lead_ele_dPhiAtVtx_%zu_%zu",ieta,icent),Form("h_lead_ele_dPhiAtVtx_%zu_%zu;Leading Electron #Delta #phi_{SC, Track};Norm. Events",ieta,icent),20,-0.2,0.2);
                h_lead_ele_EoverPInv[ieta][icent] = new TH1D(Form("h_lead_ele_EoverPInv_%zu_%zu",ieta,icent),Form("h_lead_ele_EoverPInv_%zu_%zu;Leading Electron EoverPInv;Norm. Events",ieta,icent),20,-0.2,0.2);

                h_sublead_ele_Pt[ieta][icent] = new TH1D(Form("h_sublead_ele_Pt_%zu_%zu",ieta,icent),Form("h_sublead_ele_Pt_%zu_%zu;Sub-Leading Electron p_{T} (GeV);Norm. Events",ieta,icent),net_bins,et_edges); // 30,0,150
                h_sublead_ele_eta[ieta][icent] = new TH1D(Form("h_sublead_ele_eta_%zu_%zu",ieta,icent),Form("h_sublead_ele_eta_%zu_%zu;Sub-Leading Electron SC #eta;Norm. Events",ieta,icent),10,-2.4,2.4);
                h_sublead_ele_phi[ieta][icent] = new TH1D(Form("h_sublead_ele_phi_%zu_%zu",ieta,icent),Form("h_sublead_ele_phi_%zu_%zu;Sub-Leading Electron SC #phi;Norm. Events",ieta,icent),10,-3.14,3.14);
                h_sublead_ele_HoverE[ieta][icent] = new TH1D(Form("h_sublead_ele_HoverE_%zu_%zu",ieta,icent),Form("h_sublead_ele_HoverE_%zu_%zu;Sub-Leading Electron H/E;Norm. Events",ieta,icent),20,0,0.1);
                h_sublead_ele_SIEIE[ieta][icent] = new TH1D(Form("h_sublead_ele_SIEIE_%zu_%zu",ieta,icent),Form("h_sublead_ele_SIEIE_%zu_%zu;Sub-Leading Electron #sigma_{#eta #eta} ;Norm. Events",ieta,icent),30,0,0.1);
                h_sublead_ele_EoverP[ieta][icent] = new TH1D(Form("h_sublead_ele_EoverP_%zu_%zu",ieta,icent),Form("h_sublead_ele_EoverP_%zu_%zu;Sub-Leading Electron E/P;Norm. Events",ieta,icent),20,0,20);
                h_sublead_ele_dEtaSeedAtVtx[ieta][icent] = new TH1D(Form("h_sublead_ele_dEtaSeedAtVtx_%zu_%zu",ieta,icent),Form("h_sublead_ele_dEtaSeedAtVtx_%zu_%zu;Sub-Leading Electron #Delta #eta_{seed, Track};Norm. Events",ieta,icent),20,-0.02,0.02);
                h_sublead_ele_dPhiAtVtx[ieta][icent] = new TH1D(Form("h_sublead_ele_dPhiAtVtx_%zu_%zu",ieta,icent),Form("h_sublead_ele_dPhiAtVtx_%zu_%zu;Sub-Leading Electron #Delta #phi_{SC, Track};Norm. Events",ieta,icent),20,-0.2,0.2);
                h_sublead_ele_EoverPInv[ieta][icent] = new TH1D(Form("h_sublead_ele_EoverPInv_%zu_%zu",ieta,icent),Form("h_sublead_ele_EoverPInv_%zu_%zu;Sub-Leading Electron EoverPInv;Norm. Events",ieta,icent),20,-0.2,0.2);

                // h_ele_Zmass[ieta][icent] = new TH1D(Form("h_ele_Zmass_%zu_%zu",ieta,icent),Form("h_ele_Zmass_%zu_%zu;Invariant Mass;Norm. Events",ieta,icent),60,0,180);
            }
        }
    // ----------------------------------------------------------------------------------------------------------------
    //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
    // ----------------------------------------------------------------------------------------------------------------

    Int_t nEv=photonTree.GetEntries();
    std::cout<<in_label<<" \t => "<<nEv<<"\n";

    for(int iEntry=0; iEntry< nEv; iEntry++){
        displayProgress(iEntry,nEv);
        photonTree.GetEntry(iEntry);
        HiTree.GetEntry(iEntry);
        hltanalysis.GetEntry(iEntry);
        skimanalysis.GetEntry(iEntry);
        // std::cout<<"Inside event loop = "<<iEntry;
        // printf("\n");
        float scale = 1;
        if(!isData) scale*=weight; //*Ncoll[hiBin];

        int pho_index=-1;
        float Etmax=-1;
        // if(pprimaryVertexFilter)
        // std::cout<<"Event = "<<iEntry<<"    photon = "<<phoEt->size()<<"\t pprimaryVertexFilter ="<<pprimaryVertexFilter<<"\n";
        // if(pclusterCompatibilityFilter)
        // std::cout<<"Event = "<<iEntry<<"    photon = "<<phoEt->size()<<"\t pclusterCompatibilityFilter ="<<pclusterCompatibilityFilter<<"\n";
        // if(pphfCoincFilter2Th4)
        // std::cout<<"Event = "<<iEntry<<"    photon = "<<phoEt->size()<<"\t pphfCoincFilter2Th4 ="<<pphfCoincFilter2Th4<<"\n";

        if(isData && pprimaryVertexFilter<=0 && !in_label.Contains("Express")) continue;
        if(isData && pclusterCompatibilityFilter<=0 && !in_label.Contains("Express")) continue;
        if(isData && pphfCoincFilter2Th4<=0 && !in_label.Contains("Express")) continue;
        if(isData && hiHF>6000) continue;
        
        // Leading Photon
        for(int ipho=0; ipho<(int)phoEt->size(); ipho++){
            if(!isData){
                if(pho_genMatchedIndex->at(ipho)==-1) continue;
                if(mcPID->at(pho_genMatchedIndex->at(ipho)) != 22) continue;
                if(!(abs(mcMomPID->at(pho_genMatchedIndex->at(ipho))) <= 22 || mcMomPID->at(pho_genMatchedIndex->at(ipho)) == -999) ) continue;
                if(!(mcCalIsoDR04->at(pho_genMatchedIndex->at(ipho)) < 5)) continue;
            }
            if(Etmax<phoEt->at(ipho)){
                Etmax = phoEt->at(ipho);
                pho_index = ipho;
            }
        }

        // Photon Histograms
        if(pho_index!=-1){
            for (std::size_t ieta = 0; ieta < neta; ++ieta){
                if(!(fabs(phoSCEta->at(pho_index))>=min_eta[ieta] && fabs(phoSCEta->at(pho_index))<max_eta[ieta])) continue;   

                for (std::size_t icent = 0; icent < ncent; ++icent){
                    if(!(hiBin>=min_cent[icent] && hiBin<max_cent[icent])) continue;  

                    if(phoEt->at(pho_index)<40) continue;
                    if(phoSigmaEtaEta_2012->at(pho_index)<=0.002) continue;        // Spike Cuts
                    if(fabs(pho_seedTime->at(pho_index))>=3 || pho_swissCrx->at(pho_index)>=0.9) continue;

                    // if(!(HLT_HIGEDPhoton20_v10 || HLT_HIGEDPhoton30_v10 || HLT_HIGEDPhoton40_v10)) continue;
                    if(!(HLT_HIGEDPhoton40_v10)) continue;

                    double isolation=pho_ecalClusterIsoR3->at(pho_index)+pho_hcalRechitIsoR3->at(pho_index)+pho_trackIsoR3PtCut20->at(pho_index);
                    // if(!(phoHoverE->at(pho_index)<=0.2)) continue;
                    // if(!(isolation<5)) continue;

                    // Loose Photon ID
                    if(!(phoHoverE->at(pho_index)<=0.247995)) continue;
                    if(!(phoSigmaIEtaIEta_2012->at(pho_index)<=0.012186)) continue;
                    if(!(isolation<11.697505)) continue;

                    h_lead_pho_et[ieta][icent]->Fill(         phoEt->at(pho_index), scale);
                    h_lead_pho_eta[ieta][icent]->Fill(        phoSCEta->at(pho_index), scale);
                    h_lead_pho_phi[ieta][icent]->Fill(        phoSCPhi->at(pho_index), scale);
                    h_lead_pho_HoverE[ieta][icent]->Fill(     phoHoverE->at(pho_index), scale);
                    h_lead_pho_SIEIE[ieta][icent]->Fill(      phoSigmaIEtaIEta_2012->at(pho_index), scale);
                    h_lead_pho_SumCalIso[ieta][icent]->Fill(  pho_ecalClusterIsoR3->at(pho_index)+pho_hcalRechitIsoR3->at(pho_index)+pho_trackIsoR3PtCut20->at(pho_index), scale);
                    h_lead_pho_ECALIso[ieta][icent]->Fill(    pho_ecalClusterIsoR3->at(pho_index), scale);
                    h_lead_pho_HCALIso[ieta][icent]->Fill(    pho_hcalRechitIsoR3->at(pho_index), scale);
                    h_lead_pho_TRKIso[ieta][icent]->Fill(     pho_trackIsoR3PtCut20->at(pho_index), scale);
                    h_lead_pho_PFPIso[ieta][icent]->Fill(     pfpIso3subUEec->at(pho_index), scale);
                    h_lead_pho_PFNIso[ieta][icent]->Fill(     pfnIso3subUEec->at(pho_index), scale);
                    h_lead_pho_PFCIso[ieta][icent]->Fill(     pfcIso3subUEec->at(pho_index), scale);
                    h_lead_pho_R9[ieta][icent]->Fill(         phoR9_2012->at(pho_index), scale);
                }
            }
        }

        int lead_ele_index=-1,sublead_ele_index=-1;
        float ele_Etmax=-1,ele_Etsub=-1;

        // Electron Gen Matching for MC
        std::vector<int> ele_genMatchedIndex;
        if(!isData){
            for(int iele=0;iele<(int)elePt->size();iele++){
                int matchedIndex = -1;
                float minDR = 0.3;
                for (unsigned igen = 0; igen < mcEt->size(); ++igen) {
                    if ( abs(mcPID->at(igen)) != 11) 
                    continue;                    
                    if (dr(eleEta->at(iele), elePhi->at(iele), mcEta->at(igen), mcPhi->at(igen)) < minDR){ // delta2 && mcPt->at(igen) > currentMaxPt) {
                    minDR = dr(eleEta->at(iele), elePhi->at(iele), mcEta->at(igen), mcPhi->at(igen));
                    matchedIndex = igen;
                    }
                }
                ele_genMatchedIndex.push_back(matchedIndex);
            }
        }
        // Leading Electron
        for(int iele=0;iele<(int)elePt->size();iele++){
            if(!isData){
                if(ele_genMatchedIndex.at(iele)==-1) continue;
                if(abs(mcPID->at(ele_genMatchedIndex.at(iele))) != 11) continue;       // True electron
                // if(!(abs(mcMomPID->at(ele_genMatchedIndex.at(iele))) == 23)) continue; // Parent Z
            }
            if(ele_Etmax<elePt->at(iele)){
                ele_Etmax = elePt->at(iele);
                lead_ele_index = iele;
            }
        }
        // Sub Leading Electron
        for(int iele=0;iele<(int)elePt->size();iele++){
            if(iele==lead_ele_index) continue;
            if(!isData){
                if(ele_genMatchedIndex.at(iele)==-1) continue;
                if(abs(mcPID->at(ele_genMatchedIndex.at(iele))) != 11) continue;       // True electron
                // if(!(abs(mcMomPID->at(ele_genMatchedIndex.at(iele))) == 23)) continue; // Parent Z
            }
            if(ele_Etsub<elePt->at(iele)){
                ele_Etsub = elePt->at(iele);
                sublead_ele_index = iele;
            }
        }

        if(lead_ele_index!=-1){
            for (std::size_t ieta = 0; ieta < neta; ++ieta){
                if(!(fabs(eleEta->at(lead_ele_index))>=min_eta[ieta] && fabs(eleEta->at(lead_ele_index))<max_eta[ieta])) continue;   

                for (std::size_t icent = 0; icent < ncent; ++icent){
                    if(!(hiBin>=min_cent[icent] && hiBin<max_cent[icent])) continue;  
                    
                    // if(!(HLT_HIEle20Gsf_v10 || HLT_HIEle30Gsf_v10 || HLT_HIEle40Gsf_v10 || HLT_HIEle50Gsf_v10)) continue;
                    if(!(HLT_HIEle40Gsf_v10)) continue;
                    
                    if(elePt->at(lead_ele_index)<40) continue;
                    Bool_t flagLooseEle = false;
                    if(fabs(eleEta->at(lead_ele_index))<1.442) {
                        if(eleMissHits->at(lead_ele_index)>1) continue;
                        if(eleIP3D->at(lead_ele_index)>=0.03) continue;
                    }
                    // Electron Veto ID
                    if((hiBin)>=0 && (hiBin)<60){
                        if(eleSigmaIEtaIEta_2012->at(lead_ele_index)>=0.0147) continue;
                        if(eledEtaSeedAtVtx->at(lead_ele_index)>=0.0041) continue;
                        if(eledPhiAtVtx->at(lead_ele_index)>=0.0853) continue;
                        if(eleHoverE->at(lead_ele_index)>=0.2733) continue;
                        if(eleEoverPInv->at(lead_ele_index)>=0.0367) continue;
                        flagLooseEle=true;
                    }
                    else{
                        if(eleSigmaIEtaIEta_2012->at(lead_ele_index)>=0.0113) continue;
                        if(eledEtaSeedAtVtx->at(lead_ele_index)>=0.0037) continue;
                        if(eledPhiAtVtx->at(lead_ele_index)>=0.1280) continue;
                        if(eleHoverE->at(lead_ele_index)>=0.1814) continue;
                        if(eleEoverPInv->at(lead_ele_index)>=0.1065) continue;
                        flagLooseEle=true;
                    }

                    if(!flagLooseEle) continue;

                    h_lead_ele_Pt[ieta][icent]->Fill(                elePt->at(lead_ele_index) ,scale);
                    h_lead_ele_eta[ieta][icent]->Fill(               eleEta->at(lead_ele_index) ,scale);
                    h_lead_ele_SIEIE[ieta][icent]->Fill(             eleSigmaIEtaIEta_2012->at(lead_ele_index) ,scale);
                    h_lead_ele_phi[ieta][icent]->Fill(               elePhi->at(lead_ele_index) ,scale);
                    h_lead_ele_HoverE[ieta][icent]->Fill(            eleHoverE->at(lead_ele_index) ,scale);
                    h_lead_ele_EoverP[ieta][icent]->Fill(            eleEoverP->at(lead_ele_index) ,scale);
                    h_lead_ele_dEtaSeedAtVtx[ieta][icent]->Fill(     eledEtaSeedAtVtx->at(lead_ele_index) ,scale);
                    h_lead_ele_dPhiAtVtx[ieta][icent]->Fill(         eledPhiAtVtx->at(lead_ele_index) ,scale);
                    h_lead_ele_EoverPInv[ieta][icent]->Fill(         eleEoverPInv->at(lead_ele_index) ,scale);
                    
                }
            }
        }

        if(sublead_ele_index!=-1){
            for (std::size_t ieta = 0; ieta < neta; ++ieta){
                if(!(fabs(eleEta->at(sublead_ele_index))>=min_eta[ieta] && fabs(eleEta->at(sublead_ele_index))<max_eta[ieta])) continue;   

                for (std::size_t icent = 0; icent < ncent; ++icent){
                    if(!(hiBin>=min_cent[icent] && hiBin<max_cent[icent])) continue;  
                    
                    if(!(HLT_HIEle15Ele10Gsf_v10 || HLT_HIEle15Ele10GsfMass50_v10 || HLT_HIDoubleEle10Gsf_v10 || HLT_HIDoubleEle10GsfMass50_v10 || HLT_HIDoubleEle15Gsf_v10 || HLT_HIDoubleEle15GsfMass50_v10)) continue;

                    if(elePt->at(sublead_ele_index)<20) continue;
                    Bool_t flagLooseEle = false;
                    if(fabs(eleEta->at(sublead_ele_index))<1.442) {
                        if(eleMissHits->at(sublead_ele_index)>1) continue;
                        if(eleIP3D->at(sublead_ele_index)>=0.03) continue;
                    }
                    // Electron Veto ID
                    if((hiBin)>=0 && (hiBin)<60){
                        if(eleSigmaIEtaIEta_2012->at(sublead_ele_index)>=0.0147) continue;
                        if(eledEtaSeedAtVtx->at(sublead_ele_index)>=0.0041) continue;
                        if(eledPhiAtVtx->at(sublead_ele_index)>=0.0853) continue;
                        if(eleHoverE->at(sublead_ele_index)>=0.2733) continue;
                        if(eleEoverPInv->at(sublead_ele_index)>=0.0367) continue;
                        flagLooseEle=true;
                    }
                    else{
                        if(eleSigmaIEtaIEta_2012->at(sublead_ele_index)>=0.0113) continue;
                        if(eledEtaSeedAtVtx->at(sublead_ele_index)>=0.0037) continue;
                        if(eledPhiAtVtx->at(sublead_ele_index)>=0.1280) continue;
                        if(eleHoverE->at(sublead_ele_index)>=0.1814) continue;
                        if(eleEoverPInv->at(sublead_ele_index)>=0.1065) continue;
                        flagLooseEle=true;
                    }

                    if(!flagLooseEle) continue;

                    h_sublead_ele_Pt[ieta][icent]->Fill(                elePt->at(sublead_ele_index) ,scale);
                    h_sublead_ele_eta[ieta][icent]->Fill(               eleEta->at(sublead_ele_index) ,scale);
                    h_sublead_ele_SIEIE[ieta][icent]->Fill(             eleSigmaIEtaIEta_2012->at(sublead_ele_index) ,scale);
                    h_sublead_ele_phi[ieta][icent]->Fill(               elePhi->at(sublead_ele_index) ,scale);
                    h_sublead_ele_HoverE[ieta][icent]->Fill(            eleHoverE->at(sublead_ele_index) ,scale);
                    h_sublead_ele_EoverP[ieta][icent]->Fill(            eleEoverP->at(sublead_ele_index) ,scale);
                    h_sublead_ele_dEtaSeedAtVtx[ieta][icent]->Fill(     eledEtaSeedAtVtx->at(sublead_ele_index) ,scale);
                    h_sublead_ele_dPhiAtVtx[ieta][icent]->Fill(         eledPhiAtVtx->at(sublead_ele_index) ,scale);
                    h_sublead_ele_EoverPInv[ieta][icent]->Fill(         eleEoverPInv->at(sublead_ele_index) ,scale);
                }
            }
        }        

    } // End event loop

    std::cout<<"\nNumber of selected photons = "<<h_lead_pho_et[0][0]->GetEntries()<<"\n";

    if(h_lead_pho_et[0][0]->GetEntries()==0){ std::cout<<"Return from function 0 entries => Not writing to file\n"; return;}
    
    // ------ Save histograms to File
    fout_2023->cd();
    gDirectory->mkdir(in_label);
    fout_2023->cd(in_label);

    for (std::size_t ieta = 0; ieta < neta; ++ieta){
        for (std::size_t icent = 0; icent < ncent; ++icent){
            
            h_lead_pho_et[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_eta[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_phi[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_HoverE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_SIEIE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_SumCalIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_ECALIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_HCALIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_TRKIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_PFPIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_PFNIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_PFCIso[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_pho_R9[ieta][icent]->Write("",TObject::kOverwrite);

            h_lead_ele_Pt[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_eta[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_SIEIE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_phi[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_HoverE[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_EoverP[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_dEtaSeedAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_dPhiAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_lead_ele_EoverPInv[ieta][icent]->Write("",TObject::kOverwrite);

            h_sublead_ele_Pt[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_eta[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_SIEIE[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_phi[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_HoverE[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_EoverP[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_dEtaSeedAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_dPhiAtVtx[ieta][icent]->Write("",TObject::kOverwrite);
            h_sublead_ele_EoverPInv[ieta][icent]->Write("",TObject::kOverwrite);
        }
    }

    fout_2023->cd();
}

void overlay_runs(std::vector<TH1D*> hist,std::vector<TString> histname,TString opt,std::vector<TString> eopt){
    // Based on the old Ratio Plot script from https://root.cern/doc/master/ratioplotOld_8C.html

    // opt contains options
    // "eff" = Efficiency Plot. Divide by 0th index
    // "left","right", "bcenter" = Legend location 
    // "label" = 2 before the last entries of histname is the X and Y label otherwise use the default 
    // "log" = Set log scale
    // "width" = Divide by Bin Width
    // "opt" = drawopt is last element of eopt

    if(hist.size()<2){
        std::cout<<"Not Enough Histograms"<<std::endl;
        return;
    }

    gStyle->SetOptStat(0);      // No Stat Box
    gStyle->SetOptTitle(0);     // No Title
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetHistTopMargin(0);
    if(opt.Contains("log")) gStyle->SetOptLogy(1);
    
    TString drawopt = "E1][P0"; 
    if(opt.Contains("opt"))
        drawopt = eopt.back();
    gStyle->SetPalette(1);
    // TColor *pal = new TColor();
    // good for primary marker colors      

    /*                                                                                       
    Int_t kmagenta = pal->GetColor(124,  0,124);
    Int_t kviolet  = pal->GetColor( 72,  0,190);
    Int_t kblue    = pal->GetColor(  9,  0,200);
    Int_t kazure   = pal->GetColor(  0, 48, 97);
    Int_t kcyan    = pal->GetColor(  0, 83, 98);
    Int_t kteal    = pal->GetColor(  0, 92, 46);
    Int_t kgreen   = pal->GetColor( 15, 85, 15);
    Int_t kspring  = pal->GetColor( 75, 97, 53);
    Int_t kyellow  = pal->GetColor(117,118,  0);
    Int_t korange  = pal->GetColor(101, 42,  0);
    Int_t kred     = pal->GetColor(190,  0,  3);
    Int_t kpink    = pal->GetColor(180, 35,145);

    Int_t kmagentaLight = pal->GetColor(215,165,215);
    Int_t kvioletLight  = pal->GetColor(200,160,255);
    Int_t kblueLight    = pal->GetColor(178,185,254);
    Int_t kazureLight   = pal->GetColor(153,195,225);
    Int_t kcyanLight    = pal->GetColor(140,209,224);
    Int_t ktealLight    = pal->GetColor( 92,217,141);
    Int_t kgreenLight   = pal->GetColor(135,222,135);
    Int_t kspringLight  = pal->GetColor(151,207,116);
    Int_t kyellowLight  = pal->GetColor(225,225,100);
    Int_t korangeLight  = pal->GetColor(255,168,104);
    Int_t kredLight     = pal->GetColor(253,169,179);
    Int_t kpinkLight    = pal->GetColor(255,192,224);
    */
    //const std::vector<int> colarray  = {kred,kblue,kmagenta,kgreen,kviolet,kazure,korange,kcyan,kteal,kspring,kpink,kyellow,kmagenta,kviolet,kblue};
    //const std::vector<int> markarray = {25, 22, 32, 29, 28, 39, 40,
    //                                    24, 21, 26, 23, 30, 34, 37, 41};
    //
    const std::vector<int> colarray  = {1,632,600,kMagenta+2,419,  kOrange+3, kViolet+6,898,
                                       922,910,851,877,811,804,434,606,
                                       1,632,600,616,419,804,425,898,
                                       922,910,851,877,811,804,434,606};// { 1, 2, 4, 6, 8,20,28};
    const std::vector<int> markarray = {20, 25, 22, 32, 29, 28, 39, 40,
                                        24, 21, 26, 23, 30, 34, 37, 41,
                                        20, 25, 22, 32, 29, 28, 39, 40,
                                        24, 21, 26, 23, 30, 34, 37, 41};

    
    TCanvas c;
    c.cd();
    TLegend *l;
    float leg_x1 = 0.6;
    float leg_y1 = 0.57;
    float leg_x2 = 0.85;
    float leg_y2 = 0.82;
    if(opt.Contains("bcenter")){
        leg_x1 = 0.45;
        leg_y1 = 0.13;
        leg_x2 = 0.55;
        leg_y2 = 0.42;
    }
    else if(opt.Contains("left")){
        leg_x1 = 0.15;
        leg_y1 = 0.57;
        leg_x2 = 0.25;
        leg_y2 = 0.82;
    }
    l = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2,"","brNDC");
    l->SetFillStyle(0);
    l->SetFillColor(0);
    l->SetLineColor(0);
    l->SetTextSize(0.035);
    l->SetTextFont(42);
        
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetTopMargin(0.15);
    pad1->SetRightMargin(0.05); 
    pad1->SetLeftMargin(0.12); 
    pad1->SetBottomMargin(0.05); 
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    float ymax = -99999;
    for(std::size_t ihist=0; ihist<hist.size();ihist++){

        hist[ihist]->SetLineColor(colarray[ihist]);
        hist[ihist]->SetMarkerColor(colarray[ihist]);
        hist[ihist]->SetMarkerStyle(markarray[ihist]);        
        if(opt.Contains("norm_width")){
            hist[ihist]->Scale(1.0/hist[ihist]->Integral(1,hist[ihist]->GetNbinsX()+1),"width");
        }
        else if(opt.Contains("norm")){
            hist[ihist]->Scale(1.0/hist[ihist]->Integral(0,hist[ihist]->GetNbinsX()+2));
        }
        if(opt.Contains("flow"))
            hist[ihist]->GetXaxis()->SetRange(1,hist[ihist]->GetNbinsX()+1);   
        
        if(hist[ihist]->GetMaximum()>ymax) ymax = hist[ihist]->GetMaximum();
        hist[ihist]->Draw(drawopt);
        hist[ihist]->GetXaxis()->SetLabelSize(0);
        hist[ihist]->GetXaxis()->SetTitleOffset(999999);
        l->AddEntry(hist[ihist], histname[ihist], "lep");
        if(ihist==0) drawopt+="SAME";
    }
    hist[0]->SetMinimum(0.0);
    if(!opt.Contains("log"))
        hist[0]->SetMaximum(ymax*1.1);
    else
        hist[0]->SetMaximum(ymax*2);
    l->Draw();

    TLatex latex;
    latex.SetTextSize(0.05);
    latex.DrawLatexNDC(0.12,0.9,"CMS #it{#bf{Internal}} Data Validation "+label); // PbPb 1.69nb^{-1}, pp 300.6pb^{-1} (5.02 TeV)

    if(eopt[0].Contains("Cent")){
        latex.DrawLatexNDC(0.78,0.9,eopt[0]);
    }
    else if(eopt.size()>=2){
        latex.DrawLatexNDC(0.55,0.9,eopt[0]);
        latex.DrawLatexNDC(0.82,0.9,eopt[1]);
    }    
    latex.SetTextSize(0.035);
    if(eopt.size()>2){
        for(std::size_t ind=2;ind<eopt.size() && !eopt[ind].Contains("end"); ind++){
            if(opt.Contains("bcenter")){
                leg_y2+=0.07;
                latex.DrawLatexNDC(leg_x1,leg_y2,eopt[ind]);
            }
            else{
                leg_y1-=0.07;
                latex.DrawLatexNDC(leg_x1,leg_y1,eopt[ind]);
            }
        }
    }
    
    // lower plot will be in pad
    c.cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetRightMargin(0.05); 
    pad2->SetLeftMargin(0.12); 
    pad2->SetBottomMargin(0.3); 
    pad2->SetGrid();
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    // Define the ratio plot
    TH1D *hratio = (TH1D*)hist[0]->Clone(histname[0]);
    for(std::size_t ihist=1; ihist<hist.size();ihist++){
        hratio = (TH1D*)hist[ihist]->Clone(histname[ihist]);
        hratio->SetLineColor(colarray[ihist]);
        hratio->SetMarkerColor(colarray[ihist]);
        hratio->SetMarkerStyle(markarray[ihist]);   

        if(!opt.Contains("log")){
            hratio->SetMinimum(0.);  // Define Y ..
            hratio->SetMaximum(2.0); // .. range
        }

        hratio->Sumw2();
        if(hratio)
        hratio->Divide(hist[0]);
        // hratio->Write("Ratio_"+(TString)hist[ihist]->GetName()+"_"+ihist,TObject::kOverwrite);
        hratio->Draw(drawopt);

        // To Get histogram entries in an array
        // double *bins = new double[hratio->GetSize()];
        // std::cout<<"\n{";
        // for (Int_t ieta=1;ieta<hratio->GetSize();ieta++) bins[ieta] = hratio->GetBinContent(ieta);
        // for (Int_t ieta=1;ieta<hratio->GetSize();ieta++)printf("%.7f ,",bins[ieta]); // <<" ,"
        // std::cout<<"}\n";        

        // Y axis ratio plot settings
        hratio->GetYaxis()->SetTitle("Ratio   ");
        hratio->GetYaxis()->SetNdivisions(10);
        hratio->GetYaxis()->SetTitleSize(0.12);
        hratio->GetYaxis()->SetTitleOffset(0.32);
        hratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        hratio->GetYaxis()->SetLabelSize(15);
        
        // X axis ratio plot settings
        hratio->GetXaxis()->SetTitleSize(0.12);
        hratio->GetXaxis()->SetTitleOffset(1.0);
        hratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        hratio->GetXaxis()->SetLabelSize(15);
    }
    c.cd();  

    if(opt.Contains("label")){
        hratio->GetXaxis()->SetTitle(histname.at(histname.size()-3));
        hist[0]->GetYaxis()->SetTitle(histname.at(histname.size()-2));
    }
    else{
        hratio->GetXaxis()->SetTitle(hist[0]->GetXaxis()->GetTitle());
        hist[0]->GetYaxis()->SetTitle(hist[0]->GetYaxis()->GetTitle());
    }

    // Y axis upper plot settings
        hist[0]->GetYaxis()->SetTitleSize(15);
        hist[0]->GetYaxis()->SetTitleFont(43);
        hist[0]->GetYaxis()->SetTitleOffset(2.0);
        hist[0]->GetYaxis()->SetLabelSize(0.05);

    
    gPad->Update();
    c.SaveAs(output_path +histname.back()+".png");
    c.Write(histname.back(),TObject::kOverwrite);
    std::cout<<histname.back()<<" has been saved"<<std::endl;
    delete l;
    if(opt.Contains("log")) gStyle->SetOptLogy(0);
}

void Plot_hist(std::vector<TH1D*> hist,std::vector<TString> histname,TString opt,std::vector<TString> eopt){
    // opt contains options
    // "eff" = Efficiency Plot. Divide by 0th index
    // "left","right", "bcenter" = Legend location 
    // "label" = 2 before the last entries of histname is the X and Y label otherwise use the default 
    // "log" = Set log scale
    // "width" = Divide by Bin Width
    // "OBJ" = optional stuff for later? 

    gStyle->SetOptStat(0);      // No Stat Box
    gStyle->SetOptTitle(0);     // No Title
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    if(opt.Contains("log")) gStyle->SetOptLogy(1);
    
    TString drawopt = "nostackE1][P0"; 
    const std::vector<int> colarray  = { 1,632,600,616,419,800,427,898,
                                       922,910,851,877,811,804,434,606,
                                       1,632,600,616,419,800,427,898,
                                       922,910,851,877,811,804,434,606};// { 1, 2, 4, 6, 8,20,28};
    const std::vector<int> markarray = {20, 25, 22, 32, 29, 28, 39, 40,
                                        24, 21, 26, 23, 30, 34, 37, 41,
                                        20, 25, 22, 32, 29, 28, 39, 40,
                                        24, 21, 26, 23, 30, 34, 37, 41};
    TCanvas c;
    c.cd();
    TLegend *l;
    float leg_x1 = 0.7;
    float leg_y1 = 0.7;
    float leg_x2 = 0.8;
    float leg_y2 = 0.85;
    if(opt.Contains("bcenter")){
        leg_x1 = 0.45;
        leg_y1 = 0.15;
        leg_x2 = 0.55;
        leg_y2 = 0.3;
    }
    else if(opt.Contains("left")){
        leg_x1 = 0.15;
        leg_y1 = 0.7;
        leg_x2 = 0.25;
        leg_y2 = 0.85;
    }
    l = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2,"","brNDC");
    l->SetFillStyle(0);
    l->SetFillColor(0);
    l->SetLineColor(0);
    l->SetTextSize(0.025);
    l->SetTextFont(42);

    THStack hs("hs","hs");

    for(std::size_t ihist=0; ihist<hist.size();){
        if(hist[ihist]->GetEntries()==0) continue;
        if(opt.Contains("eff")){
            if((ihist+1)>=hist.size()) break;
            hist[ihist+1]->SetLineColor(colarray[ihist]);
            hist[ihist+1]->SetMarkerColor(colarray[ihist]);
            hist[ihist+1]->SetMarkerStyle(markarray[ihist]);
            hist[ihist+1]->Divide(hist[ihist+1],hist[0],1,1,"B");
            ihist++;
        }
        else{
            hist[ihist]->SetLineColor(colarray[ihist]);
            hist[ihist]->SetMarkerColor(colarray[ihist]);
            hist[ihist]->SetMarkerStyle(markarray[ihist]);
        }
        if(opt.Contains("norm")){
            hist[ihist]->Scale(1.0/hist[ihist]->Integral(0,hist[ihist]->GetNbinsX()+2));
        }
        if(opt.Contains("width")){
            hist[ihist]->Scale(1.0,"width");
        }
        if(opt.Contains("flow"))
            hist[ihist]->GetXaxis()->SetRange(0,hist[ihist]->GetNbinsX()+2);
        hs.Add(hist[ihist]);      
        l->AddEntry(hist[ihist], histname[ihist], "lep");
        if(!opt.Contains("eff")) ihist++;
    }
    if(hs.GetNhists()==0){ std::cout<<"\n---------------------------\n"<<histname.back()<<" NOT SAVED. HISTOGRAMS ARE EMPTY\n-----------------------------\n"<<std::endl; return;}
    hs.Draw(drawopt);
    if(opt.Contains("eff")){
        hs.SetMaximum(1.1);
        hs.SetMinimum(0);
    }
    if(opt.Contains("label")){
        hs.GetXaxis()->SetTitle(histname.at(histname.size()-3));
        hs.GetYaxis()->SetTitle(histname.at(histname.size()-2));
    }
    else{
        hs.GetXaxis()->SetTitle(hist[0]->GetXaxis()->GetTitle());
        hs.GetYaxis()->SetTitle(hist[0]->GetYaxis()->GetTitle());
    }
    l->Draw();


    TLatex latex;
    latex.SetTextSize(0.035);
    latex.DrawLatexNDC(0.12,0.92,"CMS #it{#bf{Internal}} Data Validation "+label);

    if(eopt[0].Contains("Cent")){
        latex.DrawLatexNDC(0.78,0.92,eopt[0]);
    }
    else if(eopt.size()>=2){
        latex.DrawLatexNDC(0.6,0.92,eopt[0]);
        latex.DrawLatexNDC(0.78,0.92,eopt[1]);
    }    
    latex.SetTextSize(0.025);
    if(eopt.size()>2){
        for(std::size_t ind=2;ind<eopt.size() && !eopt[ind].Contains("end"); ind++){
            if(opt.Contains("bcenter")){
                leg_y2+=0.05;
                latex.DrawLatexNDC(leg_x1,leg_y2,eopt[ind]);
            }
            else{
                leg_y1-=0.05;
                latex.DrawLatexNDC(leg_x1,leg_y1,eopt[ind]);
            }
        }
    }
    gPad->Update();
    c.SaveAs(output_path +"/"+histname.back()+".png");
    c.Write(histname.back(),TObject::kWriteDelete);
    std::cout<<histname.back()<<" has been saved"<<std::endl;
    delete l;
    if(opt.Contains("log")) gStyle->SetOptLogy(0);
}


void GetFiles(char const* input, std::vector<TString>& files, int file_limit) {
    TSystemDirectory dir(input, input);
    TList *list = dir.GetListOfFiles();

    if (list) {
        TSystemFile *file;
        std::string fname;
        TIter next(list);
        while ((file = (TSystemFile*) next())) {
            fname = file->GetName();

            if (file->IsDirectory() && (fname.find(".") == std::string::npos)) {
                std::string newDir = std::string(input) + fname + "/";
                GetFiles(newDir.c_str(), files, file_limit);
            }
            else if ((fname.find(".root") != std::string::npos)) {
                files.push_back(std::string(input) + fname);
                // cout << files.back() << endl;
                if((int)files.size()>file_limit) return;
            }
        }
    }

    return;
}

void FillChain(TChain& chain, std::vector<TString>& files) {
    for (auto file : files) {
        chain.Add(file.Data());
    }
}

void displayProgress(long current, long max){
   using std::cerr;
   if (max < 2500){
      cerr <<Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0 * current / max);
      return;
   }
   if (current % (max / 2500) != 0 && current < max - 1)
      return;

   int width = 52; // Hope the terminal is at least that wide.
   int barWidth = width - 2;
   cerr << "\x1B[2K";    // Clear line
   cerr << "\x1B[2000D"; // Cursor left
   cerr << '[';
   for (int ieta = 0; ieta < barWidth; ++ieta)
   {
      if (ieta < barWidth * current / max)
      {
         cerr << '=';
      }
      else
      {
         cerr << ' ';
      }
   }
   cerr << ']';
   cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0 * current / max);
   cerr.flush();
}

int main(){ // int argc, char* argv[]
    // Run with 
    // ./Data_plot_pho_ele

    Data_overlay_2018();
    return 0;

}
