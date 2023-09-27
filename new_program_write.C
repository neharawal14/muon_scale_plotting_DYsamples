// this file is just to analyze those produced histograms by the DY smaples
//
#include<iostream>
using namespace std;
#include <TROOT.h>  
#include <TChain.h>
#include <TFile.h> 
// Header file for the classes stored in the TTree if any.                                                                                                                          
#include "string" 
#include "vector" 
#include "RooRealVar.h"   
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h" 
#include "RooHistPdf.h" 
#include "RooPolynomial.h"
#include "RooAbsArg.h"   
#include "RooPlot.h"    
//#include "TRatioPlot.h" 
#include "RooAddPdf.h" 
#include "RooFitResult.h"
#include "TAxis.h"      
#include "TH1.h"
using namespace RooFit ;
const int pT_bins=13;
const int eta_bins=14;
const TString corrections="without";
// base class having all the functions to analyze the histograms
// The functions defined in base class are called by derived class and used in that way
class class_reading{
	public:
   TFile *file = TFile::Open("/afs/cern.ch/user/n/nrawal/work/UL_sample_generation_with_without_rochester/analysis_DYsamples_with_different_variables/inclusive_2D_plots/pT_eta_bins/DY_ntuple_"+corrections+"_rochester_correction_resolution_etabins_with_Zmass_withFSR.root"); 
   TString saving_path = "/afs/cern.ch/user/n/nrawal/work/UL_sample_generation_with_without_rochester/analysis_DYsamples_with_different_variables/inclusive_2D_plots/pT_eta_bins/results_"+corrections+"_rochester_correction_Zmass/";    
   float pt_list[14] = {0,10,15,20, 25,31,35,40,43,46,50,60,90,200};
   float eta_list[15] = {-2.4, -2.1, -1.85, -1.6, -1.2, -0.8, -0.4, 0,0.4, 0.8, 1.2, 1.6, 1.85, 2.1, 2.4};
   TString pt_list_symbol[14] ={"0","10","15", "20", "25","31","35","40","43","46","50","60","90","200"};  
   TString eta_list_symbol[15] ={"-2.4", "-2.1", "-1.85", "-1.6", "-1.2",  "-0.8", "-0.4", "0","0.4", "0.8", "1.2", "1.6", "1.85", "2.1","2.4"};  
   TString bin_number_pt[13] = {"1st","2nd","3rd","4th","5th","6th","7th","8th","9th","10th","11th","12th","13th"};
   TString bin_number_eta[14] = {"1st","2nd","3rd","4th","5th", "6th", "7th", "8th", "9th", "10th","11th", "12th", "13th", "14th"};

//   float pt_list_mean[10] = {12.5,28,33, 38.5,42.5, 44.5,48,55,75,145};
   // float pt_list_mean[10] = {7.5,22.5,35,45,55,65,82.5,100,120,165};
   float mean_positive_delta_pT[pT_bins][eta_bins]={0};
   float mean_positive_error_delta_pT[pT_bins][eta_bins]={0};
   float sigma_positive_delta_pT[pT_bins][eta_bins]={0};
   float sigma_positive_error_delta_pT[pT_bins][eta_bins]={0};

   float mean_negative_delta_pT[pT_bins][eta_bins]={0};
   float mean_negative_error_delta_pT[pT_bins][eta_bins]={0};
   float sigma_negative_delta_pT[pT_bins][eta_bins]={0};
   float sigma_negative_error_delta_pT[pT_bins][eta_bins]={0};
   float chi2_ndf_positive_delta_pT[pT_bins][eta_bins] = {0};
   float chi2_ndf_negative_delta_pT[pT_bins][eta_bins] = {0};

   ofstream new_file;

    // Z mass mean values 
   float mean_positive_gen_Zmass[pT_bins][eta_bins]={0};
   float mean_negative_gen_Zmass[pT_bins][eta_bins]={0};
   float mean_positive_error_gen_Zmass[pT_bins][eta_bins]={0};
   float mean_negative_error_gen_Zmass[pT_bins][eta_bins]={0};

   float mean_positive_gen_Zmass_with_pt[pT_bins]={0};
   float mean_negative_gen_Zmass_with_pt[pT_bins]={0};
   float mean_positive_error_gen_Zmass_with_pt[pT_bins]={0};
   float mean_negative_error_gen_Zmass_with_pt[pT_bins]={0};

   float mean_positive_Zmass_with_pt[pT_bins]={0};
   float mean_negative_Zmass_with_pt[pT_bins]={0};
   float mean_positive_error_Zmass_with_pt[pT_bins]={0};
   float mean_negative_error_Zmass_with_pt[pT_bins]={0};

   float mean_positive_Zmass_FSR_with_pt[pT_bins]={0};
   float mean_negative_Zmass_FSR_with_pt[pT_bins]={0};
   float mean_positive_error_Zmass_FSR_with_pt[pT_bins]={0};
   float mean_negative_error_Zmass_FSR_with_pt[pT_bins]={0};

   float mean_positive_gen_Zmass_with_eta[eta_bins]={0};
   float mean_negative_gen_Zmass_with_eta[eta_bins]={0};
   float mean_positive_error_gen_Zmass_with_eta[eta_bins]={0};
   float mean_negative_error_gen_Zmass_with_eta[eta_bins]={0};

   float mean_positive_Zmass_with_eta[eta_bins]={0};
   float mean_negative_Zmass_with_eta[eta_bins]={0};
   float mean_positive_error_Zmass_with_eta[eta_bins]={0};
   float mean_negative_error_Zmass_with_eta[eta_bins]={0};

   float mean_positive_Zmass_FSR_with_eta[eta_bins]={0};
   float mean_negative_Zmass_FSR_with_eta[eta_bins]={0};
   float mean_positive_error_Zmass_FSR_with_eta[eta_bins]={0};
   float mean_negative_error_Zmass_FSR_with_eta[eta_bins]={0};

   float sigma_positive_gen_Zmass_with_eta[eta_bins]={0};
   float sigma_negative_gen_Zmass_with_eta[eta_bins]={0};
   float sigma_positive_error_gen_Zmass_with_eta[eta_bins]={0};
   float sigma_negative_error_gen_Zmass_with_eta[eta_bins]={0};

   float sigma_positive_Zmass_with_eta[eta_bins]={0};
   float sigma_negative_Zmass_with_eta[eta_bins]={0};
   float sigma_positive_error_Zmass_with_eta[eta_bins]={0};
   float sigma_negative_error_Zmass_with_eta[eta_bins]={0};

   float sigma_positive_Zmass_FSR_with_eta[eta_bins]={0};
   float sigma_negative_Zmass_FSR_with_eta[eta_bins]={0};
   float sigma_positive_error_Zmass_FSR_with_eta[eta_bins]={0};
   float sigma_negative_error_Zmass_FSR_with_eta[eta_bins]={0};

   float sigma_positive_gen_Zmass_with_pt[eta_bins]={0};
   float sigma_negative_gen_Zmass_with_pt[eta_bins]={0};
   float sigma_positive_error_gen_Zmass_with_pt[eta_bins]={0};
   float sigma_negative_error_gen_Zmass_with_pt[eta_bins]={0};

   float sigma_positive_Zmass_with_pt[eta_bins]={0};
   float sigma_negative_Zmass_with_pt[eta_bins]={0};
   float sigma_positive_error_Zmass_with_pt[eta_bins]={0};
   float sigma_negative_error_Zmass_with_pt[eta_bins]={0};

   float sigma_positive_Zmass_FSR_with_pt[eta_bins]={0};
   float sigma_negative_Zmass_FSR_with_pt[eta_bins]={0};
   float sigma_positive_error_Zmass_FSR_with_pt[eta_bins]={0};
   float sigma_negative_error_Zmass_FSR_with_pt[eta_bins]={0};



   float mean_positive_Zmass[pT_bins][eta_bins]={0};
   float mean_negative_Zmass[pT_bins][eta_bins]={0};
   float mean_positive_error_Zmass[pT_bins][eta_bins]={0};
   float mean_negative_error_Zmass[pT_bins][eta_bins]={0};

   float mean_positive_Zmass_FSR[pT_bins][eta_bins]={0};
   float mean_negative_Zmass_FSR[pT_bins][eta_bins]={0};
   float mean_positive_error_Zmass_FSR[pT_bins][eta_bins]={0};
   float mean_negative_error_Zmass_FSR[pT_bins][eta_bins]={0};

   float sigma_positive_Zmass_FSR[pT_bins][eta_bins]={0};
   float sigma_negative_Zmass_FSR[pT_bins][eta_bins]={0};
   float sigma_positive_error_Zmass_FSR[pT_bins][eta_bins]={0};
   float sigma_negative_error_Zmass_FSR[pT_bins][eta_bins]={0};
   float chi2_ndf_positive_Zmass_FSR[pT_bins][eta_bins] = {0};
   float chi2_ndf_negative_Zmass_FSR[pT_bins][eta_bins] = {0};

   float sigma_positive_gen_Zmass[pT_bins][eta_bins]={0};
   float sigma_negative_gen_Zmass[pT_bins][eta_bins]={0};
   float sigma_positive_error_gen_Zmass[pT_bins][eta_bins]={0};
   float sigma_negative_error_gen_Zmass[pT_bins][eta_bins]={0};
   float chi2_ndf_positive_gen_Zmass[pT_bins][eta_bins] = {0};
   float chi2_ndf_negative_gen_Zmass[pT_bins][eta_bins] = {0};

   float sigma_positive_Zmass[pT_bins][eta_bins]={0};
   float sigma_negative_Zmass[pT_bins][eta_bins]={0};
   float sigma_positive_error_Zmass[pT_bins][eta_bins]={0};
   float sigma_negative_error_Zmass[pT_bins][eta_bins]={0};
   float chi2_ndf_positive_Zmass[pT_bins][eta_bins] = {0};
   float chi2_ndf_negative_Zmass[pT_bins][eta_bins] = {0};
  
	 // hisotgram declaration to read the histograms from the root file
    TH1F * histogram_gen_Zmass; 
    TH1F * histogram_reco_Zmass; 
    TH1F * histogram_reco_Zmass_FSR; 
    TH1F * histogram_gen_positive_mu ;
    TH1F * histogram_gen_positive_eta ;
    TH1F * histogram_gen_positive_phi ;
    TH1F * histogram_id_positive_mu; 
    TH1F * histogram_gen_negative_mu; 
    TH1F * histogram_gen_negative_eta ;
    TH1F * histogram_gen_negative_phi ;

    TH1F * histogram_id_negative_mu; 
//		TH1F * histogram_eta_gen_positive_mu; 
//		TH1F * histogram_eta_gen_negative_mu; 
//		TH1F * histogram_phi_gen_positive_mu; 
//		TH1F * histogram_phi_gen_negative_mu; 
//    TH1F * histogram_pt_gen_Z; 
//    TH1F * histogram_eta_gen_Z; 
//    TH1F * histogram_phi_gen_Z; 

     // this histograms are defined in a loop since its an array  
  
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
   TH2F *histogram_delta_pT_positive_bin_with_pT;  
   TH2F *histogram_delta_pT_negative_bin_with_pT;  
 
    TH1F *histogram_pT_positive_bin[pT_bins][eta_bins] ; 
    TH1F *histogram_pT_negative_bin[pT_bins][eta_bins] ; 
    TH1F *histogram_pT_positive_bin_reco[pT_bins][eta_bins] ; 
    TH1F *histogram_pT_negative_bin_reco[pT_bins][eta_bins] ; 

     TH1F *histogram_delta_pT_positive_bin_reco[pT_bins][eta_bins]; 
     TH1F *histogram_delta_pT_negative_bin_reco[pT_bins][eta_bins]; 

     TH1F *histogram_positive_gen_Zmass[pT_bins][eta_bins] ; 
     TH1F *histogram_negative_gen_Zmass[pT_bins][eta_bins] ; 
     TH1F *histogram_positive_reco_Zmass[pT_bins][eta_bins] ; 
     TH1F *histogram_positive_reco_Zmass_FSR[pT_bins][eta_bins] ; 
     TH1F *histogram_negative_reco_Zmass[pT_bins][eta_bins] ; 
     TH1F *histogram_negative_reco_Zmass_FSR[pT_bins][eta_bins] ; 

     TH1F *histogram_positive_gen_Zmass_with_pt[pT_bins] ; 
     TH1F *histogram_negative_gen_Zmass_with_pt[pT_bins] ; 
     TH1F *histogram_positive_Zmass_with_pt[pT_bins] ; 
     TH1F *histogram_positive_Zmass_FSR_with_pt[pT_bins]; 
     TH1F *histogram_negative_Zmass_with_pt[pT_bins] ; 
     TH1F *histogram_negative_Zmass_FSR_with_pt[pT_bins] ; 

     TH1F *histogram_positive_gen_Zmass_with_eta[eta_bins] ; 
     TH1F *histogram_negative_gen_Zmass_with_eta[eta_bins] ; 
     TH1F *histogram_positive_Zmass_with_eta[eta_bins] ; 
     TH1F *histogram_positive_Zmass_FSR_with_eta[eta_bins]; 
     TH1F *histogram_negative_Zmass_with_eta[eta_bins] ; 
     TH1F *histogram_negative_Zmass_FSR_with_eta[eta_bins] ; 


     float pt_error_positive[pT_bins][eta_bins] = {0};
     float pt_error_negative[pT_bins][eta_bins] = {0};
     float pt_measured_positive[pT_bins][eta_bins] = {0}; 
     float pt_measured_negative[pT_bins][eta_bins] = {0}; 
//    TH1F *histogram_pT_negative_corresponding_positive_bin[9] ; 
//    TH1F *histogram_pT_positive_corresponding_negative_bin[9] ; 
//    TH1F *histogram_pT_negative_corresponding_positive_bin_smear[9] ; 
//    TH1F *histogram_pT_positive_corresponding_negative_bin_smear[9] ; 

		// delcared variables to analyze the mean and sigma of the distribution

   // to make a final plot for deltapT mean and pT 

	 // constructor to open the root file
	 void initializing();

  ofstream output_file;
   TFile * output_root_file;
  	// to plot histogram for pT distribution and mass distribution individually 
//   void plotting_hist(TH1F* hist_draw, TString title, TString saving_name, TString title_name_axis);
   void plotting_hist_2D(TH2F* hist_draw, TString title, TString saving_name, TString title_name_axis);
	 	// to count total number of entries in each pT bin and also above 200 GeV and print them in a text file myfile_n_entries                                                          
  	// these 3 functions are for sanity check
    // void plot_entries();   
    // void total_entries();   
    // void evaluating_mean_pT(); 
    // plotting mean and sigma distribution 		

    // to save all the histograms in pdf plots	 
    void saving_histogram_pT();
    void saving_histogram_delta_pT();
    void fitting_delta_pT_distribution();
		// graphs to plot mean and sigma, and difference due to smearing, and  also print them in  text files
    void plotting_mean_delta_pT(TString);
    void plotting_mean_with_pt(TString);
    void plotting_mean_with_eta(TString);
    void temperature_plot_delta_pT();
    void fitting_DSCB(TH1F * , TString saving_name, TString title_name, int, int, TString);
    void fitting_DSCB_Zmass(TH1F * , TString saving_name, TString title_name, int, int, TString,TString );
    void fitting_DSCB_Zmass_with_pt(TH1F * , TString saving_name, TString title_name, int, TString,TString );
    void fitting_DSCB_Zmass_with_eta(TH1F * , TString saving_name, TString title_name, int, TString,TString );
    void fitting_BW(TH1F * hist_fit, TString saving_name, TString title_name );

    void fitting_DSCB_Zmass_alone(TH1F *hist_fit, TString saving_name, TString title_name);
    void fitting_Zmass_distribution();
    void fitting_Zmass_FSR_distribution();
};

// declaring other class and its variables - the one which is actually used -=> inherited class from class reading
/*class derived_class_reading : public class_reading{
	public:

   //TString pt_list_symbol[11] ={"0","25","31","35","40","43","46","50","60","90","200"};  

   //TString eta_list_symbol[5] ={"0","0.9","1.2","2.0","2.4"};  
//  std::pair<float,float>mean_sigma_Zmass_gauss; 
//  std::pair<float,float> mean_sigma_Zmass_gauss_smear; 
  
  void fitting_delta_pT();
  void saving_histograms();
};

void derived_class_reading :: saving_histograms(){
    saving_histogram_pT();
    saving_histogram_delta_pT();
}

void derived_class_reading :: fitting_delta_pT(){
   fitting_delta_pT_distribution();
   plotting_mean_delta_pT();
   temperature_plot_delta_pT();
}
*/
  // declaring base class functions
	//
	// to plot mean Jpsi after smearing
  
// this is to plot mean Jpsi mass before smearing
 /*  void class_reading :: plotting_hist(TH1F * hist_draw, TString title, TString saving_name, TString title_name_axis){
    
    gStyle->SetOptStat(000001111);
    gStyle->SetStatFontSize(0.02); 
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.1); 

    TCanvas * c = new TCanvas("c",title,900,600);
    c->cd();
    hist_draw->GetXaxis()->SetTitle(title_name_axis+"");
    hist_draw->Draw();
    c->SaveAs(saving_path+saving_name+".pdf");
   }
*/
  void class_reading :: plotting_hist_2D(TH2F * hist_draw, TString title, TString saving_name, TString title_name_axis){
    hist_draw->SetStats(0); 
//    gStyle->SetOptStat(000001111);
//    gStyle->SetStatFontSize(0.02); 
//    gStyle->SetStatW(0.2);
//    gStyle->SetStatH(0.1); 

    TCanvas * c = new TCanvas("c",title,900,600);
    c->cd();
    c->SetRightMargin(0.2);
    //hist_draw->GetXaxis()->SetTitle(title_name_axis+"");
    hist_draw->Draw("colZ");
    hist_draw->SetTitle(title+ " : "+corrections+" rochester corrections");
    hist_draw->GetXaxis()->SetTitle("pT (GeV)");
    hist_draw->GetYaxis()->SetTitle("#eta");
    c->SaveAs(saving_path+saving_name+".pdf");
   }

void class_reading :: saving_histogram_delta_pT(){
   TString string_delta_pT_distribution_negative_reco, string_delta_pT_distribution_positive_reco;  
  // TString string_delta_pT_distribution_negative_smear, string_delta_pT_distribution_positive_smear;  
   // histograms and corresponding other pT distribution
    gStyle->SetOptStat(000001111);
    gStyle->SetStatFontSize(0.02); 
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.1); 

  TCanvas *canvas_pT_positive_reco;
  TCanvas *canvas_pT_negative_reco;
  for(int i=0; i<pT_bins; i++){
  for(int j=0; j<eta_bins; j++){
  canvas_pT_positive_reco = new TCanvas("canvas_pT_positive_reco", "reco #frac{#Delta pT}{pT} distribution in each bin #mu^{+}", 900,600);
  canvas_pT_positive_reco->cd();
  histogram_delta_pT_positive_bin_reco[i][j]->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
  histogram_delta_pT_positive_bin_reco[i][j]->Draw();
  histogram_delta_pT_positive_bin_reco[i][j]->GetXaxis()->SetTitleOffset(1.2);
  histogram_delta_pT_positive_bin_reco[i][j]->GetXaxis()->CenterTitle(true);
  string_delta_pT_distribution_positive_reco ="delta_pT_positive_" + bin_number_pt[i] + "_eta_"+bin_number_eta[j]+"_bin_reco";
  canvas_pT_positive_reco->SaveAs(saving_path+"pT_plots/"+string_delta_pT_distribution_positive_reco+".pdf");
  canvas_pT_positive_reco->Close();
 
  canvas_pT_negative_reco = new TCanvas("canvas_pT_negative_reco", "reco #frac{#Delta pT}{pT}  distribution in each bin #mu^{+}", 900,600);
  canvas_pT_negative_reco->cd();
  histogram_delta_pT_negative_bin_reco[i][j]->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
  histogram_delta_pT_negative_bin_reco[i][j]->Draw();
  histogram_delta_pT_negative_bin_reco[i][j]->GetXaxis()->SetTitleOffset(1.2);
  histogram_delta_pT_negative_bin_reco[i][j]->GetXaxis()->CenterTitle(true);
  string_delta_pT_distribution_negative_reco ="delta_pT_negative_" + bin_number_pt[i] + "_eta_"+bin_number_eta[j]+"_bin_reco";
  canvas_pT_negative_reco->SaveAs(saving_path+"pT_plots/"+string_delta_pT_distribution_negative_reco+".pdf");
  canvas_pT_negative_reco->Close();
  }
 }

}
 
void class_reading :: saving_histogram_pT(){
    TString string_pT_distribution_negative, string_pT_distribution_positive;  
    TString string_pT_distribution_negative_reco, string_pT_distribution_positive_reco;  
    gStyle->SetOptStat(000001111);
    gStyle->SetStatFontSize(0.02); 
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.1); 
  //  double range_yaxis[10] = {100000, 1000000, 1000000, 1000000, 500000, 80000, 40000, 20000, 15000, 10000}; 
  // histograms and corresponding other pT distribution
  TCanvas *canvas_pT_positive;
  TCanvas *canvas_pT_negative;
  std::cout<<"entered into saving histogram pT "<<std::endl;
  
  for(int i=0; i<pT_bins; i++){
  for(int j=0; j<eta_bins; j++){
    std::cout<<"entered into saving histogram pT : i  "<<i<<" j "<<j<<std::endl;
  canvas_pT_positive = new TCanvas("canvas_pT_positive", "pT distribution in each bin #mu^{+}", 900,600);
  canvas_pT_positive->cd();

  histogram_pT_positive_bin[i][j]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_positive_bin[i][j]->Draw();
  string_pT_distribution_positive ="pT_positive_" + bin_number_pt[i] + "_eta " + bin_number_eta[j] +"bin";
  canvas_pT_positive->SaveAs(saving_path+"pT_plots/"+string_pT_distribution_positive+".pdf");
  canvas_pT_positive->Close();
  
    std::cout<<"could not finish into saving histogram pT "<<std::endl;
  pt_measured_positive[i][j] =  histogram_pT_positive_bin[i][j]->GetMean();
  std::cout<<"pt measured value : positive "<<pt_measured_positive[i][j]<<std::endl;
//  pt_error_positive[i][j] = histogram_pT_positive_bin[i][j]->GetRMS()/histogram_pT_positive_bin[i][j]->Integral();
  pt_error_positive[i][j] = histogram_pT_positive_bin[i][j]->GetMeanError();

  canvas_pT_negative = new TCanvas("canvas_pT_negative", "pT distribution in each bin #mu^{-}", 900,600);
  canvas_pT_negative->cd();
  histogram_pT_negative_bin[i][j]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_negative_bin[i][j]->Draw();
  string_pT_distribution_negative ="pT_negative_" + bin_number_pt[i] + "_eta " + bin_number_eta[j] +"bin";
  canvas_pT_negative->SaveAs(saving_path+"pT_plots/"+string_pT_distribution_negative+".pdf");
  canvas_pT_negative->Close();

  pt_measured_negative[i][j] =  histogram_pT_negative_bin[i][j]->GetMean();
  //pt_error_negative[i][j] = histogram_pT_negative_bin[i][j]->GetRMS()/histogram_pT_negative_bin[i][j]->Integral();
  pt_error_negative[i][j] = histogram_pT_negative_bin[i][j]->GetMeanError();
  std::cout<<" pt values "<<pt_measured_negative[i][j]<<" error "<<pt_error_negative[i][j]<<std::endl;

  }
} 
} //end of saving histogram

void class_reading :: initializing(){

      new_file.open("new_file.txt"); 
    output_file.open("mean_sigma.txt");
//initialize all the 9 histograms positive and negative one before and afer smearing
    std::cout<<"started initialization "<<std::endl;
    output_root_file = new TFile(corrections+"_rochester_histograms_delta_pT_with_fit.root","RECREATE");
    histogram_gen_Zmass = (TH1F*) file->Get("hist_gen_Zmass"); 
    histogram_reco_Zmass = (TH1F*) file->Get("hist_reco_Zmass"); 
    histogram_reco_Zmass_FSR = (TH1F*) file->Get("hist_reco_Zmass_FSR"); 
   histogram_gen_positive_mu= (TH1F*) file->Get("hist_gen_positive_mu") ;
   histogram_id_positive_mu= (TH1F*) file->Get("hist_id_positive_mu"); 
   histogram_gen_negative_mu= (TH1F*) file->Get("hist_gen_negative_mu"); 
   histogram_id_negative_mu= (TH1F*) file->Get("hist_id_negative_mu"); 

   histogram_gen_positive_eta= (TH1F*) file->Get("hist_gen_positive_eta") ;
   histogram_gen_negative_eta= (TH1F*) file->Get("hist_gen_negative_eta") ;
   histogram_gen_positive_phi= (TH1F*) file->Get("hist_gen_positive_phi") ;
   histogram_gen_negative_phi= (TH1F*) file->Get("hist_gen_negative_phi") ;
    
   histogram_delta_pT_positive_bin_with_pT = (TH2F*) file->Get("delta_pT_positive_bin_with_pT") ; 
   histogram_delta_pT_negative_bin_with_pT = (TH2F*) file->Get("delta_pT_negative_bin_with_pT") ; 
//	 histogram_eta_gen_positive_mu= (TH1F*) file->Get("hist_eta_gen_positive_mu"); 
//	 histogram_eta_gen_negative_mu= (TH1F*) file->Get("hist_eta_gen_negative_mu"); 
//	 histogram_phi_gen_positive_mu= (TH1F*) file->Get("hist_phi_gen_positive_mu"); 
//	 histogram_phi_gen_negative_mu= (TH1F*) file->Get("hist_phi_gen_negative_mu"); 
//
    std::cout<<"continuing initialization "<<std::endl;
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution

   for(int i=0; i<pT_bins; i++){
   for(int j=0; j<eta_bins; j++){
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
   
		TString title_pT_positive_bin = TString::Format("pT_positive_bin[%d][%d]",i,j);
		TString title_pT_negative_bin = TString::Format("pT_negative_bin[%d][%d]",i,j);
    TString title_pT_positive_bin_reco = TString::Format("pT_positive_bin_reco[%d][%d]",i,j);
		TString title_pT_negative_bin_reco = TString::Format("pT_negative_bin_reco[%d][%d]",i,j);

		TString title_pT_positive_gen_Zmass =  TString::Format("hist_positive_gen_Zmass[%d][%d]",i,j);
		TString title_pT_negative_gen_Zmass =  TString::Format("hist_negative_gen_Zmass[%d][%d]",i,j);
		TString title_pT_positive_Zmass =  TString::Format("hist_positive_reco_Zmass[%d][%d]",i,j);
		TString title_pT_negative_Zmass =  TString::Format("hist_negative_reco_Zmass[%d][%d]",i,j);
		TString title_pT_positive_Zmass_FSR =  TString::Format("hist_positive_reco_Zmass_FSR[%d][%d]",i,j);
		TString title_pT_negative_Zmass_FSR =  TString::Format("hist_negative_reco_Zmass_FSR[%d][%d]",i,j);

    std::cout<<"print title "<<title_pT_positive_bin<<std::endl;
    // delta pT distribution
    TString title_delta_pT_positive_bin_reco = TString::Format("delta_pT_positive_bin_reco[%d][%d]",i,j);
		TString title_delta_pT_negative_bin_reco = TString::Format("delta_pT_negative_bin_reco[%d][%d]",i,j);

    // getting Z mass from root file
    histogram_positive_gen_Zmass[i][j] = (TH1F*) file->Get(title_pT_positive_gen_Zmass) ; 
    histogram_positive_reco_Zmass[i][j] = (TH1F*) file->Get(title_pT_positive_Zmass) ; 
    histogram_positive_reco_Zmass_FSR[i][j] = (TH1F*) file->Get(title_pT_positive_Zmass_FSR) ; 
    histogram_negative_gen_Zmass[i][j] = (TH1F*) file->Get(title_pT_negative_gen_Zmass) ; 
    histogram_negative_reco_Zmass[i][j] = (TH1F*) file->Get(title_pT_negative_Zmass) ; 
    histogram_negative_reco_Zmass_FSR[i][j] = (TH1F*) file->Get(title_pT_negative_Zmass_FSR) ; 

    histogram_pT_positive_bin[i][j] = (TH1F*) file->Get(title_pT_positive_bin) ; 
    histogram_pT_negative_bin[i][j] = (TH1F*) file->Get(title_pT_negative_bin) ; 
    histogram_pT_positive_bin_reco[i][j] = (TH1F*) file->Get(title_pT_positive_bin_reco) ; 
    histogram_pT_negative_bin_reco[i][j] = (TH1F*) file->Get(title_pT_negative_bin_reco) ; 
    histogram_delta_pT_positive_bin_reco[i][j] = (TH1F*) file->Get(title_delta_pT_positive_bin_reco) ; 
    histogram_delta_pT_negative_bin_reco[i][j] = (TH1F*) file->Get(title_delta_pT_negative_bin_reco) ; 
	}
}
   for(int i=0; i<pT_bins; i++){
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
   
		TString title_pT_positive_gen_Zmass =  TString::Format("hist_positive_gen_Zmass_with_pt[%d]",i);
		TString title_pT_negative_gen_Zmass =  TString::Format("hist_negative_gen_Zmass_with_pt[%d]",i);
		TString title_pT_positive_Zmass =  TString::Format("hist_positive_reco_Zmass_with_pt[%d]",i);
		TString title_pT_negative_Zmass =  TString::Format("hist_negative_reco_Zmass_with_pt[%d]",i);
		TString title_pT_positive_Zmass_FSR =  TString::Format("hist_positive_reco_Zmass_FSR_with_pt[%d]",i);
		TString title_pT_negative_Zmass_FSR =  TString::Format("hist_negative_reco_Zmass_FSR_with_pt[%d]",i);

    // getting Z mass from root file
    histogram_positive_gen_Zmass_with_pt[i] = (TH1F*) file->Get(title_pT_positive_gen_Zmass) ; 
    histogram_positive_Zmass_with_pt[i] = (TH1F*) file->Get(title_pT_positive_Zmass) ; 
    histogram_positive_Zmass_FSR_with_pt[i] = (TH1F*) file->Get(title_pT_positive_Zmass_FSR) ; 
    histogram_negative_gen_Zmass_with_pt[i] = (TH1F*) file->Get(title_pT_negative_gen_Zmass) ; 
    histogram_negative_Zmass_with_pt[i] = (TH1F*) file->Get(title_pT_negative_Zmass) ; 
    histogram_negative_Zmass_FSR_with_pt[i] = (TH1F*) file->Get(title_pT_negative_Zmass_FSR) ; 
	}
   for(int i=0; i<eta_bins; i++){
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
   
		TString title_pT_positive_gen_Zmass =  TString::Format("hist_positive_gen_Zmass_with_eta[%d]",i);
		TString title_pT_negative_gen_Zmass =  TString::Format("hist_negative_gen_Zmass_with_eta[%d]",i);
		TString title_pT_positive_Zmass =  TString::Format("hist_positive_reco_Zmass_with_eta[%d]",i);
		TString title_pT_negative_Zmass =  TString::Format("hist_negative_reco_Zmass_with_eta[%d]",i);
		TString title_pT_positive_Zmass_FSR =  TString::Format("hist_positive_reco_Zmass_FSR_with_eta[%d]",i);
		TString title_pT_negative_Zmass_FSR =  TString::Format("hist_negative_reco_Zmass_FSR_with_eta[%d]",i);

    // getting Z mass from root file
    histogram_positive_gen_Zmass_with_eta[i] = (TH1F*) file->Get(title_pT_positive_gen_Zmass) ; 
    histogram_positive_Zmass_with_eta[i] = (TH1F*) file->Get(title_pT_positive_Zmass) ; 
    histogram_positive_Zmass_FSR_with_eta[i] = (TH1F*) file->Get(title_pT_positive_Zmass_FSR) ; 
    histogram_negative_gen_Zmass_with_eta[i] = (TH1F*) file->Get(title_pT_negative_gen_Zmass) ; 
    histogram_negative_Zmass_with_eta[i] = (TH1F*) file->Get(title_pT_negative_Zmass) ; 
    histogram_negative_Zmass_FSR_with_eta[i] = (TH1F*) file->Get(title_pT_negative_Zmass_FSR) ; 
    }
std::cout<<"finished initialization "<<std::endl;
}
  
void class_reading::fitting_delta_pT_distribution(){
    //fitting both deltapT distributions with a gaussian functions
    gStyle->SetOptFit();
    std::cout<<"declaring canvas 1 fitting deltapTfunction"<<std::endl;
    std::cout<<"inside class reading fitting delta pTfunction"<<std::endl; 
     TDirectory *dir_fit = (TDirectory*) output_root_file->mkdir("fits");
     dir_fit->cd();

    std::pair<float , float> mean_sigma_positive_delta_pT[pT_bins][eta_bins]; 
    output_file<<"positive muon "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "delta_pT_plots/delta_pT_positive_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "#frac{#Delta pT}{pT} (#mu^{+}) : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        //mean_sigma_positive_delta_pT[i][j] = fitting_DSCB(histogram_delta_pT_positive_bin_reco[i][j], saving_string, title_name);
        fitting_DSCB(histogram_delta_pT_positive_bin_reco[i][j], saving_string, title_name, i, j, "positive");
        output_file<<"mean : "<<mean_positive_delta_pT[i][j]<<" sigma "<<sigma_positive_delta_pT[i][j]<<"mean erro : "<<mean_positive_error_delta_pT[i][j]<<" sigma "<<sigma_positive_error_delta_pT[i][j]<<std::endl;
//        std::cout<<"mean sigma pair "<<mean_positive_delta_pT[i][j]<<" second "<<sigma_positive_delta_pT[i][j]<<std::endl;
      }
    }

    output_file<<"negative muon "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "delta_pT_plots/delta_pT_negative_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "#frac{#Delta pT}{pT} (#mu^{-}) : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        //mean_sigma_positive_delta_pT[i][j] = fitting_DSCB(histogram_delta_pT_positive_bin_reco[i][j], saving_string, title_name);
        fitting_DSCB(histogram_delta_pT_negative_bin_reco[i][j], saving_string, title_name, i, j,"negative");
        output_file<<"mean : "<<mean_negative_delta_pT[i][j]<<" sigma "<<sigma_negative_delta_pT[i][j]<<"mean erro : "<<mean_negative_error_delta_pT[i][j]<<" sigma "<<sigma_negative_error_delta_pT[i][j]<<std::endl;
//        std::cout<<"mean sigma pair "<<mean_positive_delta_pT[i][j]<<" second "<<sigma_positive_delta_pT[i][j]<<std::endl;
      }
    }


    std::cout<<"will save 1st deltapT postiive"<<std::endl;

/*
    gStyle->SetOptFit();
    std::cout<<"declaring canvas 1 fitting deltapTfunction"<<std::endl;
    TCanvas *canvas_1;   
    TCanvas *canvas_2;
    int j;
    std::cout<<"inside class reading fitting delta pTfunction"<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        canvas_1 = new TCanvas("canvas_1","positive muon deltapT",1000,1000);
        TF1 * f1 = new TF1("f1","gaus",-0.1,0.1);
        canvas_1->cd(); 
        histogram_delta_pT_positive_bin_reco[i][j]->Fit(f1,"R");
        gStyle->SetOptFit(00111111);
        mean_positive_delta_pT[i][j] = f1->GetParameter(1);
        std::cout<<"mean positive "<<mean_positive_delta_pT[i][j]<<std::endl;
        mean_positive_error_delta_pT[i][j] = f1->GetParError(1);
        sigma_positive_delta_pT[i][j] = f1->GetParameter(2);
        sigma_positive_error_delta_pT[i][j] = f1->GetParError(2);
       canvas_1->SaveAs(saving_path+"delta_pT_plots/delta_pT_positive_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j]+".pdf"); 
       canvas_1->Close();
      }
    }
    std::cout<<"will save 1st deltapT postiive"<<std::endl;
   canvas_1->Delete();

    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
           canvas_2 = new TCanvas("canvas_2","negative muon deltapT",1000,1000);

        TF1 * f2 = new TF1("f2","gaus",-0.1,0.1);
        canvas_2->cd();
        histogram_delta_pT_negative_bin_reco[i][j]->Fit(f2,"R");
        gStyle->SetOptFit(00111111);
        mean_negative_delta_pT[i][j] = f2->GetParameter(1);
        std::cout<<"mean negative "<<mean_negative_delta_pT[i][j]<<std::endl;
        mean_negative_error_delta_pT[i][j] = f2->GetParError(1);
        sigma_negative_delta_pT[i][j] = f2->GetParameter(2);
        sigma_negative_error_delta_pT[i][j] = f2->GetParError(2);
       canvas_2->SaveAs(saving_path+"delta_pT_plots/delta_pT_negative_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j]+".pdf"); 
       canvas_2->Close();


        std::cout<<" mean values "<<mean_negative_delta_pT[i][j]<<" error "<<mean_negative_error_delta_pT[i][j]<<std::endl;
     }
    }
    std::cout<<"will save 2nd deltapT postiive"<<std::endl;
*/
    } // end of fitting delta pT

  void class_reading::fitting_Zmass_distribution(){
    //fitting both deltapT distributions with a gaussian functions
    gStyle->SetOptFit();


        fitting_BW(histogram_gen_Zmass, "inclusive_gen_Zmass", " GEN Z mass : inclusive");
        fitting_DSCB_Zmass_alone(histogram_gen_Zmass, "inclusive_gen_Zmass_with_DSCB", " Gen Z mass : inclusive");
        fitting_DSCB_Zmass_alone(histogram_reco_Zmass, "inclusive_reco_Zmass", " RECO Z mass : inclusive");
        fitting_DSCB_Zmass_alone(histogram_reco_Zmass_FSR, "inclusive_reco_Zmass_FSR", " RECO Z mass with FSR : inclusive");
        //fitting_DSCB_Zmass_alone(histogram_reco_Zmass_FSR, "inclusive_reco_Zmass_FSR", " RECO Z mass : inclusive");
        
     TDirectory *dir_fit = (TDirectory*) output_root_file->mkdir("fits_Zmass");
     dir_fit->cd();
     // gen Z  mass fitting
    output_file<<"positive GENZ  "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "Zmass_plots/gen_Zmass_positive_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "GEN Z (GeV) : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        //fitting_DSCB_Zmass(histogram_positive_gen_Zmass[i][j], saving_string, title_name, i, j, "gen","positive");
        fitting_DSCB_Zmass(histogram_positive_gen_Zmass[i][j], saving_string, title_name, i, j,"Zmass_gen","positive");
        output_file<<"mean : "<<mean_positive_gen_Zmass[i][j]<<" sigma "<<sigma_positive_gen_Zmass[i][j]<<"mean erro : "<<mean_positive_error_gen_Zmass[i][j]<<" sigma "<<sigma_positive_error_gen_Zmass[i][j]<<std::endl;
      }
    }
    output_file<<"negative GENZ  "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "Zmass_plots/gen_Zmass_negative_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "GEN Z (GeV) : #mu^{-} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
//        fitting_DSCB_Zmass(histogram_negative_gen_Zmass[i][j], saving_string, title_name, i, j, "gen","negative");
        fitting_DSCB_Zmass(histogram_negative_gen_Zmass[i][j], saving_string, title_name, i, j,"Zmass_gen", "negative");
        output_file<<"mean : "<<mean_negative_gen_Zmass[i][j]<<" sigma "<<sigma_negative_gen_Zmass[i][j]<<"mean erro : "<<mean_negative_error_gen_Zmass[i][j]<<" sigma "<<sigma_negative_error_gen_Zmass[i][j]<<std::endl;
      }
    }


     // reco muons fitting
    output_file<<"positive muon "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "Zmass_plots/reco_Zmass_positive_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "m(#mu^{+} #mu^{-}) GeV : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        fitting_DSCB_Zmass(histogram_positive_reco_Zmass[i][j], saving_string, title_name, i, j, "Zmass","positive");
        output_file<<"mean : "<<mean_positive_Zmass[i][j]<<" sigma "<<sigma_positive_Zmass[i][j]<<"mean erro : "<<mean_positive_error_Zmass[i][j]<<" sigma "<<sigma_positive_error_Zmass[i][j]<<std::endl;
      }
    }

    output_file<<"negative muon "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "Zmass_plots/reco_Zmass_negative_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "m(#mu^{+} #mu^{-}) GeV : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        fitting_DSCB_Zmass(histogram_negative_reco_Zmass[i][j], saving_string, title_name, i, j, "Zmass","negative");
        output_file<<"mean : "<<mean_negative_Zmass[i][j]<<" sigma "<<sigma_negative_Zmass[i][j]<<"mean erro : "<<mean_negative_error_Zmass[i][j]<<" sigma "<<sigma_negative_error_Zmass[i][j]<<std::endl;
      }
    }

    TString saving_string; 
    TString title_name;
    for(int i=0; i<pT_bins; i++){
        saving_string = "Zmass_plots_inclusive/gen_Zmass_positive_mu_"+bin_number_pt[i] +"_with_pt";
        title_name = "GEN Z (GeV) : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV  ";
        fitting_DSCB_Zmass_with_pt(histogram_positive_gen_Zmass_with_pt[i], saving_string, title_name, i,"Zmass_gen_with_pt","positive");
        
        saving_string = "Zmass_plots_inclusive/gen_Zmass_negative_mu_"+bin_number_pt[i] +"_with_pt";
        title_name = "GEN Z (GeV) : #mu^{-} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV  ";
        fitting_DSCB_Zmass_with_pt(histogram_negative_gen_Zmass_with_pt[i], saving_string, title_name, i,"Zmass_gen_with_pt","negative");

        saving_string = "Zmass_plots_inclusive/reco_Zmass_positive_mu_"+bin_number_pt[i] + "_with_pt";
        title_name = "Reco Z (GeV) : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV  ";
        fitting_DSCB_Zmass_with_pt(histogram_positive_Zmass_with_pt[i], saving_string, title_name, i,"Zmass_reco_with_pt","positive");

        saving_string = "Zmass_plots_inclusive/reco_Zmass_negative_mu_"+bin_number_pt[i] + "_with_pt";
        title_name = "Reco Z (GeV) : #mu^{-} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV  ";
        fitting_DSCB_Zmass_with_pt(histogram_negative_Zmass_with_pt[i], saving_string, title_name, i,"Zmass_reco_with_pt","negative");

      }
    // fitting eta bins 
    for(int i=0; i<eta_bins; i++){
        saving_string = "Zmass_plots_inclusive/gen_Zmass_positive_mu_"+bin_number_eta[i] +"_with_eta";
        title_name = "GEN Z (GeV) : #mu^{+} : "+eta_list_symbol[i] + " #leq #eta < "+ eta_list_symbol[i+1];
        fitting_DSCB_Zmass_with_eta(histogram_positive_gen_Zmass_with_eta[i], saving_string, title_name, i,"Zmass_gen_with_eta","positive");
        
        saving_string = "Zmass_plots_inclusive/gen_Zmass_negative_mu_"+bin_number_eta[i] +"_with_eta";
        title_name = "GEN Z (GeV) : #mu^{-} : "+eta_list_symbol[i] + " #leq #eta < "+ eta_list_symbol[i+1] ;
        fitting_DSCB_Zmass_with_eta(histogram_negative_gen_Zmass_with_eta[i], saving_string, title_name, i,"Zmass_gen_with_eta","negative");

        saving_string = "Zmass_plots_inclusive/reco_Zmass_positive_mu_"+bin_number_eta[i] + "_with_eta";
        title_name = "Reco Z (GeV) : #mu^{+} : "+eta_list_symbol[i] + " #leq #eta < "+ eta_list_symbol[i+1] ;
        fitting_DSCB_Zmass_with_eta(histogram_positive_Zmass_with_eta[i], saving_string, title_name, i,"Zmass_reco_with_eta","positive");

        saving_string = "Zmass_plots_inclusive/reco_Zmass_negative_mu_"+bin_number_eta[i] + "_with_eta";
        title_name = "Reco Z (GeV) : #mu^{-} : "+eta_list_symbol[i] + " #leq #eta < "+ eta_list_symbol[i+1];
        fitting_DSCB_Zmass_with_eta(histogram_negative_Zmass_with_eta[i], saving_string, title_name, i,"Zmass_reco_with_eta","negative");
      }
 
  std::cout<<"saved all Z mass fitting "<<std::endl;
   }


  void class_reading::fitting_Zmass_FSR_distribution(){
    //fitting both deltapT distributions with a gaussian functions
    gStyle->SetOptFit();
     TDirectory *dir_fit = (TDirectory*) output_root_file->mkdir("fits_Zmass_FSR");
     dir_fit->cd();

     // reco muons fitting
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "Zmass_plots/reco_Zmass_FSR_positive_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "m(#mu^{+} #mu^{-}) GeV : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        fitting_DSCB_Zmass(histogram_positive_reco_Zmass_FSR[i][j], saving_string, title_name, i, j, "Zmass_FSR","positive");
      }
    }

    output_file<<"negative muon "<<std::endl;
    for(int i=0; i<pT_bins; i++){
    for(int j=0; j<eta_bins; j++){
        TString saving_string = "Zmass_plots/reco_Zmass_FSR_negative_mu_"+bin_number_pt[i] + "_eta_"+bin_number_eta[j];
        TString title_name = "m(#mu^{+} #mu^{-}) GeV : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV , "+eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
        fitting_DSCB_Zmass(histogram_negative_reco_Zmass_FSR[i][j], saving_string, title_name, i, j, "Zmass_FSR","negative");
      }
    }

    TString saving_string, title_name;
    for(int i=0; i<pT_bins; i++){
        saving_string = "Zmass_plots_inclusive/reco_Zmass_FSR_positive_mu_"+bin_number_pt[i] + "_with_pt";
        title_name = "Reco Z (GeV) with FSR : #mu^{+} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV  ";
        fitting_DSCB_Zmass_with_pt(histogram_positive_Zmass_FSR_with_pt[i], saving_string, title_name, i,"Zmass_FSR_reco_with_pt","positive");

        saving_string = "Zmass_plots_inclusive/reco_Zmass_FSR_negative_mu_"+bin_number_pt[i] + "_with_pt";
        title_name = "Reco Z (GeV) with FSR : #mu^{-} : "+pt_list_symbol[i] + " #leq pT(GEN) < "+ pt_list_symbol[i+1] + " GeV  ";
        fitting_DSCB_Zmass_with_pt(histogram_negative_Zmass_FSR_with_pt[i], saving_string, title_name, i,"Zmass_FSR_reco_with_pt","negative");
      }
    // fitting eta bins 
    for(int i=0; i<eta_bins; i++){

        saving_string = "Zmass_plots_inclusive/reco_Zmass_FSR_positive_mu_"+bin_number_eta[i] + "_with_eta";
        title_name = "Reco Z (GeV) with FSR : #mu^{+} : "+eta_list_symbol[i] + " #leq #eta < "+ eta_list_symbol[i+1] ;
        fitting_DSCB_Zmass_with_eta(histogram_positive_Zmass_FSR_with_eta[i], saving_string, title_name, i,"Zmass_FSR_reco_with_eta","positive");

        saving_string = "Zmass_plots_inclusive/reco_Zmass_FSR_negative_mu_"+bin_number_eta[i] + "_with_eta";
        title_name = "Reco Z (GeV) with FSR : #mu^{-} : "+eta_list_symbol[i] + " #leq #eta < "+ eta_list_symbol[i+1];
        fitting_DSCB_Zmass_with_eta(histogram_negative_Zmass_FSR_with_eta[i], saving_string, title_name, i,"Zmass_FSR_reco_with_eta","negative");
      }
 

  std::cout<<"saved all Z mass fitting "<<std::endl;
   }

 void class_reading::temperature_plot_delta_pT(){

   TH2F *temp_mean_positive = new TH2F("temp_mean_positive", "Mean (#frac{#Delta pT}{pT}) : #mu^{+}", pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_mean_negative = new TH2F("temp_mean_negative", "Mean (#frac{#Delta pT}{pT}) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_positive = new TH2F("temp_sigma_positive", "#sigma (#frac{#Delta pT}{pT}) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_negative = new TH2F("temp_sigma_negative", "#sigma (#frac{#Delta pT}{pT}) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_positive = new TH2F("chi2_fit_positive", "#chi^2/ndf (#frac{#Delta pT}{pT}) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_negative = new TH2F("chi2_fit_negative", "#chi^2/ndf (#frac{#Delta pT}{pT}) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);

   TH2F *temp_mean_positive_Zmass_FSR = new TH2F("temp_mean_positive_Zmass_FSR", "Mean (Z mass) with FSR : #mu^{+}", pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_mean_negative_Zmass_FSR = new TH2F("temp_mean_negative_Zmass_FSR", "Mean (Z mass) with FSR : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_positive_Zmass_FSR = new TH2F("temp_sigma_positive_Zmass_FSR", "#sigma (Z mass) with FSR : #mu^{+}", pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_negative_Zmass_FSR = new TH2F("temp_sigma_negative_Zmass_FSR", "#sigma (Z mass) with FSR : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);

   TH2F *temp_mean_positive_Zmass = new TH2F("temp_mean_positive_Zmass", "Mean (Z mass) : #mu^{+}", pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_mean_negative_Zmass = new TH2F("temp_mean_negative_Zmass", "Mean (Z mass) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_positive_Zmass = new TH2F("temp_sigma_positive_Zmass", "#sigma (Z mass)  : #mu^{+}", pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_sigma_negative_Zmass = new TH2F("temp_sigma_negative_Zmass", "#sigma (Z mass)  : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);


//   TH2F *temp_sigma_positive_Zmass = new TH2F("temp_sigma_positive_Zmass", "#sigma (Z mass) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
//   TH2F *temp_sigma_negative_Zmass = new TH2F("temp_sigma_negative_Zmass", "#sigma (Z mass) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_positive_Zmass = new TH2F("chi2_fit_positive_Zmass", "#chi^2/ndf (Z mass) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_negative_Zmass = new TH2F("chi2_fit_negative_Zmass", "#chi^2/ndf (Z mass) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_positive_Zmass_FSR = new TH2F("chi2_fit_positive_Zmass_FSR", "#chi^2/ndf (Z mass) with FSR : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_negative_Zmass_FSR = new TH2F("chi2_fit_negative_Zmass_FSR", "#chi^2/ndf (Z mass) with FSR : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);


   TH2F *temp_mean_positive_gen_Zmass = new TH2F("temp_mean_positive_gen_Zmass", "Mean (GEN Z mass) : #mu^{+}", pT_bins, pt_list, eta_bins, eta_list);
   TH2F *temp_mean_negative_gen_Zmass = new TH2F("temp_mean_negative_gen_Zmass", "Mean (GEN Z mass) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
//   TH2F *temp_sigma_positive_gen_Zmass = new TH2F("temp_sigma_positive_gen_Zmass", "#sigma (GEN Z mass) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
//   TH2F *temp_sigma_negative_gen_Zmass = new TH2F("temp_sigma_negative_gen_Zmass", "#sigma (GEN Z mass) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_positive_gen_Zmass = new TH2F("chi2_fit_positive_gen_Zmass", "#chi^2/ndf (GEN Z mass) : #mu^{+}",pT_bins, pt_list, eta_bins, eta_list);
   TH2F *chi2_fit_negative_gen_Zmass = new TH2F("chi2_fit_negative_gen_Zmass", "#chi^2/ndf (GEN Z mass) : #mu^{-}",pT_bins, pt_list, eta_bins, eta_list);

   for(int i=0; i< pT_bins; i++){ 
    for(int j=0; j<eta_bins;j++){
      temp_mean_positive->SetBinContent(i+1,j+1, mean_positive_delta_pT[i][j]);
      temp_sigma_positive->SetBinContent(i+1,j+1, sigma_positive_delta_pT[i][j]);
      temp_mean_negative->SetBinContent(i+1,j+1, mean_negative_delta_pT[i][j]);
      temp_sigma_negative->SetBinContent(i+1,j+1, sigma_negative_delta_pT[i][j]);
      chi2_fit_positive->SetBinContent(i+1, j+1, chi2_ndf_positive_delta_pT[i][j]);
      chi2_fit_negative->SetBinContent(i+1, j+1, chi2_ndf_negative_delta_pT[i][j]);

      temp_mean_positive_Zmass_FSR->SetBinContent(i+1,j+1, mean_positive_Zmass_FSR[i][j]);
      //temp_sigma_positive_Zmass->SetBinContent(i+1,j+1, sigma_positive_Zmass[i][j]);
      temp_mean_negative_Zmass_FSR->SetBinContent(i+1,j+1, mean_negative_Zmass_FSR[i][j]);
      chi2_fit_positive_Zmass_FSR->SetBinContent(i+1, j+1, chi2_ndf_positive_Zmass_FSR[i][j]);
      chi2_fit_negative_Zmass_FSR->SetBinContent(i+1, j+1, chi2_ndf_negative_Zmass_FSR[i][j]);

      temp_mean_positive_Zmass->SetBinContent(i+1,j+1, mean_positive_Zmass[i][j]);
      //temp_sigma_positive_Zmass->SetBinContent(i+1,j+1, sigma_positive_Zmass[i][j]);
      temp_mean_negative_Zmass->SetBinContent(i+1,j+1, mean_negative_Zmass[i][j]);
      //temp_sigma_negative_Zmass->SetBinContent(i+1,j+1, sigma_negative_Zmass[i][j]);
      chi2_fit_positive_Zmass->SetBinContent(i+1, j+1, chi2_ndf_positive_Zmass[i][j]);
      chi2_fit_negative_Zmass->SetBinContent(i+1, j+1, chi2_ndf_negative_Zmass[i][j]);

      temp_mean_positive_gen_Zmass->SetBinContent(i+1,j+1, mean_positive_gen_Zmass[i][j]);
      //temp_sigma_positive_gen_Zmass->SetBinContent(i+1,j+1, sigma_positive_gen_Zmass[i][j]);
      temp_mean_negative_gen_Zmass->SetBinContent(i+1,j+1, mean_negative_gen_Zmass[i][j]);
      //temp_sigma_negative_gen_Zmass->SetBinContent(i+1,j+1, sigma_negative_gen_Zmass[i][j]);
      chi2_fit_positive_gen_Zmass->SetBinContent(i+1, j+1, chi2_ndf_positive_gen_Zmass[i][j]);
      chi2_fit_negative_gen_Zmass->SetBinContent(i+1, j+1, chi2_ndf_negative_gen_Zmass[i][j]);

    }
   }
  // plotting this temperature plots 
  
    temp_mean_positive->GetZaxis()->SetRangeUser(-0.01,0.004);
    temp_sigma_positive->GetZaxis()->SetRangeUser(0,0.07);
    temp_mean_negative->GetZaxis()->SetRangeUser(-0.01,0.004);
    temp_sigma_negative->GetZaxis()->SetRangeUser(0,0.07);

		class_reading::plotting_hist_2D(temp_mean_positive," Mean #frac{#Delta pT}{pT} : #mu^{+} ","temperature_plots/temperature_mean_delta_pT_over_pT_2D_positive_distribution", " ");
		class_reading::plotting_hist_2D(temp_mean_negative," Mean #frac{#Delta pT}{pT} : #mu^{-} ","temperature_plots/temperature_mean_delta_pT_over_pT_2D_negative_distribution", " ");
//		class_reading::plotting_hist_2D(temp_sigma_positive," #sigma #frac{#Delta pT}{pT} : #mu^{+} ","temperature_plots/temperature_sigma_delta_pT_over_pT_2D_positive_distribution", " ");
//		class_reading::plotting_hist_2D(temp_sigma_negative," #sigma #frac{#Delta pT}{pT} : #mu^{-} ","temperature_plots/temperature_sigma_delta_pT_over_pT_2D_negative_distribution", " ");
		class_reading::plotting_hist_2D(chi2_fit_positive," #chi^2/ndf (#frac{#Delta pT}{pT}) : #mu^{+} ","chi2_plots/chi2_fit_positive", " ");
		class_reading::plotting_hist_2D(chi2_fit_negative," #chi^2/ndf (#frac{#Delta pT}{pT}) : #mu^{-} ","chi2_plots/chi2_fit_negative", " ");

    // Z mass plots
    class_reading::plotting_hist_2D(temp_mean_positive_Zmass," Mean (reco Zmass) : #mu^{+} ","temperature_plots/temperature_mean_Zmass_2D_positive_distribution", " ");
		class_reading::plotting_hist_2D(temp_mean_negative_Zmass," Mean (reco Zmass) : #mu^{-} ","temperature_plots/temperature_mean_Zmass_2D_negative_distribution", " ");
//		class_reading::plotting_hist_2D(temp_sigma_positive_Zmass," #sigma (reco Zmass) : #mu^{+} ","temperature_plots/temperature_sigma_Zmass_2D_positive_distribution", " ");
//		class_reading::plotting_hist_2D(temp_sigma_negative_Zmass," #sigma (reco Zmass) : #mu^{-} ","temperature_plots/temperature_sigma_Zmass_2D_negative_distribution", " ");
		class_reading::plotting_hist_2D(chi2_fit_positive_Zmass," #chi^2/ndf (reco Zmass) : #mu^{+} ","chi2_plots/chi2_fit_positive_Zmass", " ");
		class_reading::plotting_hist_2D(chi2_fit_negative_Zmass," #chi^2/ndf (reco Zmass) : #mu^{-} ","chi2_plots/chi2_fit_negative_Zmass", " ");

    // Z massgen  plots
    class_reading::plotting_hist_2D(temp_mean_positive_gen_Zmass," Mean (gen Zmass) : #mu^{+} ","temperature_plots/temperature_gen_Zmass_2D_positive_distribution", " ");
		class_reading::plotting_hist_2D(temp_mean_negative_gen_Zmass," Mean (gen Zmass) : #mu^{-} ","temperature_plots/temperature_gen_mass_2D_negative_distribution", " ");
//		class_reading::plotting_hist_2D(temp_sigma_positive_gen_Zmass," #sigma (gen Zmass) : #mu^{+} ","temperature_plots/temperature_sigma_gen_Zmass_positive_distribution", " ");
//		class_reading::plotting_hist_2D(temp_sigma_negative_gen_Zmass," #sigma (gen Zmass) : #mu^{-} ","temperature_plots/temperature_sigma_gen_Zmass_negative_distribution", " ");
		class_reading::plotting_hist_2D(chi2_fit_positive_gen_Zmass," #chi^2/ndf (gen Zmass) : #mu^{+} ","chi2_plots/chi2_fit_positive_genZmass", " ");
		class_reading::plotting_hist_2D(chi2_fit_negative_gen_Zmass," #chi^2/ndf (gen Zmass) : #mu^{-} ","chi2_plots/chi2_fit_negative_genZmass", " ");
 
    // Z mass reco with FSR plots

    class_reading::plotting_hist_2D(temp_mean_positive_Zmass_FSR," Mean (reco Zmass) with FSR : #mu^{+} ","temperature_plots/temperature_mean_Zmass_FSR_2D_positive_distribution", " ");
		class_reading::plotting_hist_2D(temp_mean_negative_Zmass_FSR," Mean (reco Zmass) with FSR : #mu^{-} ","temperature_plots/temperature_mean_Zmass_FSR_2D_negative_distribution", " ");
//		class_reading::plotting_hist_2D(temp_sigma_positive_Zmass," #sigma (reco Zmass) : #mu^{+} ","temperature_plots/temperature_sigma_Zmass_2D_positive_FSR_distribution", " ");
//	class_reading::plotting_hist_2D(temp_sigma_negative_Zmass," #sigma (reco Zmass) : #mu^{-} ","temperature_plots/temperature_sigma_Zmass_2D_negative_FSR_distribution", " ");
		class_reading::plotting_hist_2D(chi2_fit_positive_Zmass_FSR," #chi^2/ndf (reco Zmass) with FSR : #mu^{+} ","chi2_plots/chi2_fit_positive_Zmass_FSR", " ");
		class_reading::plotting_hist_2D(chi2_fit_negative_Zmass_FSR," #chi^2/ndf (reco Zmass) with FSR : #mu^{-} ","chi2_plots/chi2_fit_negative_Zmass_FSR", " ");


    // done with temperature plots
 }

// plot new with eta things  
void class_reading::plotting_mean_with_eta(TString type){
     TString graph_name_positive = "";
     TString saving_name_positive = "";

     TString graph_name_negative = "";
     TString saving_name_negative = "";
     TString graph_name_positive_sigma = "";
     TString graph_name_negative_sigma = "";

     TString saving_name_positive_sigma = "";
     TString saving_name_negative_sigma = "";

     TString string_title_positive = "";
     TString string_title_negative = "";
     TString string_title_positive_sigma = "";
     TString string_title_negative_sigma = "";

   //  TString string_title_positive = "";
   //  TString string_title_negative = "";

      float range_low ;
      float range_high ; 
      float range_high_sigma ; 

     if(type=="Zmass_gen_with_eta"){
       graph_name_positive = "mean_gen_Zmass_positive_with_eta";
       graph_name_negative = "mean_gen_Zmass_negative_with_eta";
       graph_name_positive_sigma = "sigam_gen_Zmass_positive_with_eta";
       graph_name_negative_sigma = "sigam_gen_Zmass_negative_with_eta";

       saving_name_negative = "mean_sigma_Zmass_plots_with_eta/graph_resolution_gen_Zmass_mean_negative_with_eta";
       saving_name_positive = "mean_sigma_Zmass_plots_with_eta/graph_resolution_gen_Zmass_mean_positive_with_eta";
       saving_name_positive_sigma =  "mean_sigma_plots_with_eta/graph_resolution_gen_Zmass_sigma_positive_with_eta";
       saving_name_negative_sigma =  "mean_sigma_plots_with_eta/graph_resolution_gen_Zmass_sigma_negative_with_eta";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;
    } 

     if(type=="Zmass_reco_with_eta"){
       graph_name_positive = "mean_reco_Zmass_positive_with_eta";
       graph_name_negative = "mean_reco_Zmass_negative_with_eta";
       graph_name_positive_sigma = "sigam_reco_Zmass_positive_with_eta";
       graph_name_negative_sigma = "sigam_reco_Zmass_negative_with_eta";

       saving_name_negative = "mean_sigma_Zmass_plots_with_eta/graph_resolution_reco_Zmass_mean_negative_with_eta";
       saving_name_positive = "mean_sigma_Zmass_plots_with_eta/graph_resolution_reco_Zmass_mean_positive_with_eta";
       saving_name_positive_sigma =  "mean_sigma_plots_with_eta/graph_resolution_reco_Zmass_sigma_positive_with_eta";
       saving_name_negative_sigma =  "mean_sigma_plots_with_eta/graph_resolution_reco_Zmass_sigma_negative_with_eta";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;
    }


    if(type=="Zmass_reco_FSR_with_eta"){
       graph_name_positive = "mean_reco_Zmass_FSR_positive_with_eta";
       graph_name_negative = "mean_reco_Zmass_FSR_negative_with_eta";
       graph_name_positive_sigma = "sigam_reco_Zmass_FSR_positive_with_eta";
       graph_name_negative_sigma = "sigam_reco_Zmass_FSR_negative_with_eta";

       saving_name_negative = "mean_sigma_Zmass_plots_with_eta/graph_resolution_reco_Zmass_FSR_mean_negative_with_eta";
       saving_name_positive = "mean_sigma_Zmass_plots_with_eta/graph_resolution_reco_Zmass_FSR_mean_positive_with_eta";
       saving_name_positive_sigma =  "mean_sigma_plots_with_eta/graph_resolution_reco_Zmass_FSR_sigma_positive_with_eta";
       saving_name_negative_sigma =  "mean_sigma_plots_with_eta/graph_resolution_reco_Zmass_FSR_sigma_negative_with_eta";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;
    }


    TGraphErrors *gr_mean1; 
    TGraphErrors *gr_mean2; 
    TGraphErrors *gr1; 
    TGraphErrors *gr2; 

    TCanvas *c1_mean_graph;
    TCanvas *c2_mean_graph;
    TCanvas *c1_graph;
    TCanvas *c2_graph;
 
    TString string_title;

     if(type=="Zmass_gen_with_eta") {
       string_title_positive = TString::Format("Mean Zmass GEN (#mu^{+})" );
       string_title_negative = TString::Format("Mean Zmass GEN (#mu^{-})" );
       string_title_positive_sigma = TString::Format("#sigma : Z mass GEN (#mu^{+})");
       string_title_negative_sigma = TString::Format("#sigma : Z mass GEN (#mu^{-})");
     }


     if(type=="Zmass_reco_with_eta") {
       string_title_positive = TString::Format("Mean Zmass reco (#mu^{+})");
       string_title_negative = TString::Format("Mean Zmass reco  (#mu^{-})");
       string_title_positive_sigma = TString::Format("#sigma : Z mass reco  (#mu^{+})");
       string_title_negative_sigma = TString::Format("#sigma : Z mass reco  (#mu^{-})");
     }
     if(type=="Zmass_reco_FSR_with_eta") {
       string_title_positive = TString::Format("Mean Zmass reco with FSR (#mu^{+})");
       string_title_negative = TString::Format("Mean Zmass reco  with FSR (#mu^{-})");
       string_title_positive_sigma = TString::Format("#sigma : Z mass reco with FSR (#mu^{+})");
       string_title_negative_sigma = TString::Format("#sigma : Z mass reco with FSR (#mu^{-})");
     }

   
       std::cout<<"started graph plotting "<<std::endl; 


       if(type=="Zmass_gen_with_eta") {
        gr_mean1 = new TGraphErrors(eta_bins, eta_list, mean_positive_gen_Zmass_with_eta, 0, mean_positive_error_gen_Zmass_with_eta);
        gr_mean2 = new TGraphErrors(eta_bins, eta_list, mean_negative_gen_Zmass_with_eta, 0, mean_negative_error_gen_Zmass_with_eta);
       }

       if(type=="Zmass_reco_FSR_with_eta") {
        gr_mean1 = new TGraphErrors(eta_bins, eta_list, mean_positive_Zmass_FSR_with_eta, 0, mean_positive_error_Zmass_FSR_with_eta);
        gr_mean2 = new TGraphErrors(eta_bins, eta_list, mean_negative_Zmass_FSR_with_eta, 0, mean_negative_error_Zmass_FSR_with_eta);
       }
       if(type=="Zmass_reco_with_eta") {
        gr_mean1 = new TGraphErrors(eta_bins, eta_list, mean_positive_Zmass_with_eta,0 , mean_positive_error_Zmass_with_eta);
        gr_mean2 = new TGraphErrors(eta_bins, eta_list, mean_negative_Zmass_with_eta,0, mean_negative_error_Zmass_with_eta);
       }

    c1_mean_graph = new TCanvas("c1_mean_graph","canvas",1000,1000);
    c1_mean_graph->cd();
    c1_mean_graph->SetLeftMargin(0.15);
    c1_mean_graph->SetTopMargin(0.2);
    gr_mean1->SetTitle(string_title_positive);
    gr_mean1->GetXaxis()->SetTitleOffset(1.2);
    gr_mean1->GetYaxis()->SetTitleOffset(1.2);
    gr_mean1->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean1->GetYaxis()->SetRangeUser(range_low, range_high);
    gr_mean1->Draw("AP*");
    gr_mean1->SetName(graph_name_positive);
    gr_mean1->GetXaxis()->SetRangeUser(-2.5, 2.5);
    c1_mean_graph->Update();
   // c1_mean_graph->SaveAs(saving_path+"mean_sigma_plots/graph_resolution_pT_mean_positive_eta_"+bin_number_eta[j]+".pdf");
    c1_mean_graph->SaveAs(saving_path+saving_name_positive+".pdf");
    c1_mean_graph->Close();

    // making an overlapping TGraph for comparing variation in eta bins 


    //c2_mean_graph = new TCanvas("c2_mean_graph","Mean #frac{#Delta pT}{pT} mu^{-}",1000,1000);
    c2_mean_graph = new TCanvas("c2_mean_graph","canvas",1000,1000);
    c2_mean_graph->cd();
    c2_mean_graph->SetLeftMargin(0.15);
    c2_mean_graph->SetTopMargin(0.2);
    gr_mean2->SetTitle(string_title_negative);
    gr_mean2->GetXaxis()->SetTitleOffset(1.2);
    gr_mean2->GetYaxis()->SetTitleOffset(1.2);
    gr_mean2->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean2->GetYaxis()->SetRangeUser(range_low,range_high);
    gr_mean2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    //gr_mean2->SetName("mean_delta_pT_negative_"+bin_number_eta[j]);
    gr_mean2->SetName(graph_name_negative);
    gr_mean2->Draw("AP*");
    c2_mean_graph->Update();
    //c2_mean_graph->SaveAs(saving_path+"mean_sigma_plots/graph_resolution_pT_mean_negative_eta_"+bin_number_eta[j]+".pdf");
    c2_mean_graph->SaveAs(saving_path+saving_name_negative+".pdf");
    c2_mean_graph->Close();

    c1_graph = new TCanvas("c1_graph","canvas",1000,1000);
    c1_graph->cd();
    c1_graph->SetLeftMargin(0.15);
    c1_graph->SetTopMargin(0.2);
    if(type=="Zmass_gen_with_eta") { 
    gr1 = new TGraphErrors(eta_bins, eta_list, sigma_positive_gen_Zmass_with_eta, 0, sigma_positive_error_gen_Zmass_with_eta);
    gr2 = new TGraphErrors(eta_bins, eta_list, sigma_negative_gen_Zmass_with_eta, 0,  sigma_negative_error_gen_Zmass_with_eta);
    }

    if(type=="Zmass_reco_with_eta") { 
    gr1 = new TGraphErrors(eta_bins, eta_list, sigma_positive_Zmass_with_eta, 0, sigma_positive_error_Zmass_with_eta);
    gr2 = new TGraphErrors(eta_bins, eta_list, sigma_negative_Zmass_with_eta, 0 ,  sigma_negative_error_Zmass_with_eta);
    }
    if(type=="Zmass_reco_FSR_with_eta") {
    gr1 = new TGraphErrors(eta_bins, eta_list, sigma_positive_Zmass_FSR_with_eta, 0, sigma_positive_error_Zmass_FSR_with_eta);
    gr2 = new TGraphErrors(eta_bins, eta_list, sigma_negative_Zmass_FSR_with_eta, 0,  sigma_negative_error_Zmass_FSR_with_eta);
    }
    
    gr1->SetTitle(string_title_positive_sigma);
    gr1->GetXaxis()->SetTitleOffset(1.2);
    gr1->GetYaxis()->SetTitleOffset(1.2);
    gr1->GetXaxis()->SetTitle("#eta");
    gr1->GetYaxis()->SetRangeUser(0, range_high_sigma);
    gr1->GetXaxis()->SetRangeUser(-2.5,2.5);
    gr1->SetName(graph_name_positive_sigma);
    gr1->Draw("AP*");
    c1_graph->Update();
    c1_graph->SaveAs(saving_path+saving_name_positive_sigma+".pdf");
    c1_graph->Close();

    c2_graph = new TCanvas("c2_graph","canvas",1000,1000);
    c2_graph->cd();
    c2_graph->SetLeftMargin(0.15);
    c2_graph->SetTopMargin(0.2);
    gr2->SetTitle(string_title_negative_sigma);
    gr2->GetXaxis()->SetTitleOffset(1.2);
    gr2->GetYaxis()->SetTitleOffset(1.2);
    gr2->GetXaxis()->SetTitle("#eta");
    gr2->GetYaxis()->SetRangeUser(0,range_high_sigma);
    gr2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    gr2->Draw("AP*");
    gr2->SetName(graph_name_negative_sigma);
    c2_graph->Update();
    c2_graph->SaveAs(saving_path+saving_name_negative_sigma+".pdf");
    c2_graph->Close(); 
  } // end of plotting mean delta with pT



// plot new with pT things  
void class_reading::plotting_mean_with_pt(TString type){
     TString graph_name_positive = "";
     TString saving_name_positive = "";

     TString graph_name_negative = "";
     TString saving_name_negative = "";
     TString graph_name_positive_sigma = "";
     TString graph_name_negative_sigma = "";

     TString saving_name_positive_sigma = "";
     TString saving_name_negative_sigma = "";

     TString string_title_positive = "";
     TString string_title_negative = "";
     TString string_title_positive_sigma = "";
     TString string_title_negative_sigma = "";

   //  TString string_title_positive = "";
   //  TString string_title_negative = "";

      float range_low ;
      float range_high ; 
      float range_high_sigma ; 
     if(type=="Zmass_gen_with_pt"){
       graph_name_positive = "mean_gen_Zmass_positive_with_pt";
       graph_name_negative = "mean_gen_Zmass_negative_with_pt";
       graph_name_positive_sigma = "sigam_gen_Zmass_positive_with_pt";
       graph_name_negative_sigma = "sigam_gen_Zmass_negative_with_pt";

       saving_name_negative = "mean_sigma_Zmass_plots_with_pt/graph_resolution_gen_Zmass_mean_negative_with_pt";
       saving_name_positive = "mean_sigma_Zmass_plots_with_pt/graph_resolution_gen_Zmass_mean_positive_with_pt";
       saving_name_positive_sigma =  "mean_sigma_plots_with_pt/graph_resolution_gen_Zmass_sigma_positive_with_pt";
       saving_name_negative_sigma =  "mean_sigma_plots_with_pt/graph_resolution_gen_Zmass_sigma_negative_with_pt";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;
    }
     if(type=="Zmass_reco_with_pt"){
       graph_name_positive = "mean_reco_Zmass_positive_with_pt";
       graph_name_negative = "mean_reco_Zmass_negative_with_pt";
       graph_name_positive_sigma = "sigam_reco_Zmass_positive_with_pt";
       graph_name_negative_sigma = "sigam_reco_Zmass_negative_with_pt";

       saving_name_negative = "mean_sigma_Zmass_plots_with_pt/graph_resolution_reco_Zmass_mean_negative_with_pt";
       saving_name_positive = "mean_sigma_Zmass_plots_with_pt/graph_resolution_reco_Zmass_mean_positive_with_pt";
       saving_name_positive_sigma =  "mean_sigma_plots_with_pt/graph_resolution_reco_Zmass_sigma_positive_with_pt";
       saving_name_negative_sigma =  "mean_sigma_plots_with_pt/graph_resolution_reco_Zmass_sigma_negative_with_pt";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;
    }
     if(type=="Zmass_reco_FSR_with_pt"){
       graph_name_positive = "mean_reco_Zmass_FSR_positive_with_pt";
       graph_name_negative = "mean_reco_Zmass_FSR_negative_with_pt";
       graph_name_positive_sigma = "sigam_reco_Zmass_FSR_positive_with_pt";
       graph_name_negative_sigma = "sigam_reco_Zmass_FSR_negative_with_pt";

       saving_name_negative = "mean_sigma_Zmass_plots_with_pt/graph_resolution_reco_Zmass_FSR_mean_negative_with_pt";
       saving_name_positive = "mean_sigma_Zmass_plots_with_pt/graph_resolution_reco_Zmass_FSR_mean_positive_with_pt";
       saving_name_positive_sigma =  "mean_sigma_plots_with_pt/graph_resolution_reco_Zmass_FSR_sigma_positive_with_pt";
       saving_name_negative_sigma =  "mean_sigma_plots_with_pt/graph_resolution_reco_Zmass_FSR_sigma_negative_with_pt";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;
    }

    TGraphErrors *gr_mean1; 
    TGraphErrors *gr_mean2; 
    TGraphErrors *gr1; 
    TGraphErrors *gr2; 

    TCanvas *c1_mean_graph;
    TCanvas *c2_mean_graph;
    TCanvas *c1_graph;
    TCanvas *c2_graph;
 
    TString string_title;

     if(type=="Zmass_gen_with_pt") {
       string_title_positive = TString::Format("Mean Zmass GEN (#mu^{+})" );
       string_title_negative = TString::Format("Mean Zmass GEN (#mu^{-})" );
       string_title_positive_sigma = TString::Format("#sigma : Z mass GEN (#mu^{+})");
       string_title_negative_sigma = TString::Format("#sigma : Z mass GEN (#mu^{-})");
     }


     if(type=="Zmass_reco_with_pt") {
       string_title_positive = TString::Format("Mean Zmass reco (#mu^{+})");
       string_title_negative = TString::Format("Mean Zmass reco  (#mu^{-})");
       string_title_positive_sigma = TString::Format("#sigma : Z mass reco  (#mu^{+})");
       string_title_negative_sigma = TString::Format("#sigma : Z mass reco  (#mu^{-})");
     }
     if(type=="Zmass_reco_FSR_with_pt") {
       string_title_positive = TString::Format("Mean Zmass reco with FSR (#mu^{+})");
       string_title_negative = TString::Format("Mean Zmass reco  with FSR (#mu^{-})");
       string_title_positive_sigma = TString::Format("#sigma : Z mass reco with FSR (#mu^{+})");
       string_title_negative_sigma = TString::Format("#sigma : Z mass reco with FSR (#mu^{-})");
     }

   
       std::cout<<"started graph plotting "<<std::endl; 

       float new_pt_measured_positive[pT_bins] ={0};
       float new_pt_measured_negative[pT_bins] ={0};
       float new_pt_error_positive[pT_bins] ={0};
       float new_pt_error_negative[pT_bins] ={0};
       for(int j=0; j<eta_bins; j++){
       for(int i=0; i<pT_bins; i++){
        new_pt_measured_positive[i]= pt_measured_positive[i][j];
        new_pt_error_positive[i]= pt_error_positive[i][j];
        new_pt_measured_negative[i]= pt_measured_negative[i][j];
        new_pt_error_negative[i]= pt_error_negative[i][j];
        }
       }


       if(type=="Zmass_gen_with_pt") {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, mean_positive_gen_Zmass_with_pt, new_pt_error_positive,mean_positive_error_gen_Zmass_with_pt);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, mean_negative_gen_Zmass_with_pt, new_pt_error_negative,mean_negative_error_gen_Zmass_with_pt);
       }

       if(type=="Zmass_reco_FSR_with_pt") {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, mean_positive_Zmass_FSR_with_pt, new_pt_error_positive, mean_positive_error_Zmass_FSR_with_pt);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, mean_negative_Zmass_FSR_with_pt, new_pt_error_negative, mean_negative_error_Zmass_FSR_with_pt);
       }
       if(type=="Zmass_reco_with_pt") {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, mean_positive_Zmass_with_pt, new_pt_error_positive, mean_positive_error_Zmass_with_pt);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, mean_negative_Zmass_with_pt, new_pt_error_negative, mean_negative_error_Zmass_with_pt);
       }

    c1_mean_graph = new TCanvas("c1_mean_graph","canvas",1000,1000);
    c1_mean_graph->cd();
    c1_mean_graph->SetLeftMargin(0.15);
    c1_mean_graph->SetTopMargin(0.2);
    gr_mean1->SetTitle(string_title_positive);
    gr_mean1->GetXaxis()->SetTitleOffset(1.2);
    gr_mean1->GetYaxis()->SetTitleOffset(1.2);
    gr_mean1->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean1->GetYaxis()->SetRangeUser(range_low, range_high);
    gr_mean1->Draw("AP*");
    gr_mean1->SetName(graph_name_positive);
    gr_mean1->GetXaxis()->SetRangeUser(0,200);
    c1_mean_graph->Update();
   // c1_mean_graph->SaveAs(saving_path+"mean_sigma_plots/graph_resolution_pT_mean_positive_eta_"+bin_number_eta[j]+".pdf");
    c1_mean_graph->SaveAs(saving_path+saving_name_positive+".pdf");
    c1_mean_graph->Close();

    // making an overlapping TGraph for comparing variation in eta bins 


    //c2_mean_graph = new TCanvas("c2_mean_graph","Mean #frac{#Delta pT}{pT} mu^{-}",1000,1000);
    c2_mean_graph = new TCanvas("c2_mean_graph","canvas",1000,1000);
    c2_mean_graph->cd();
    c2_mean_graph->SetLeftMargin(0.15);
    c2_mean_graph->SetTopMargin(0.2);
    gr_mean2->SetTitle(string_title_negative);
    gr_mean2->GetXaxis()->SetTitleOffset(1.2);
    gr_mean2->GetYaxis()->SetTitleOffset(1.2);
    gr_mean2->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean2->GetYaxis()->SetRangeUser(range_low,range_high);
    gr_mean2->GetXaxis()->SetRangeUser(0,200);
    //gr_mean2->SetName("mean_delta_pT_negative_"+bin_number_eta[j]);
    gr_mean2->SetName(graph_name_negative);
    gr_mean2->Draw("AP*");
    c2_mean_graph->Update();
    //c2_mean_graph->SaveAs(saving_path+"mean_sigma_plots/graph_resolution_pT_mean_negative_eta_"+bin_number_eta[j]+".pdf");
    c2_mean_graph->SaveAs(saving_path+saving_name_negative+".pdf");
    c2_mean_graph->Close();

    c1_graph = new TCanvas("c1_graph","canvas",1000,1000);
    c1_graph->cd();
    c1_graph->SetLeftMargin(0.15);
    c1_graph->SetTopMargin(0.2);
    if(type=="Zmass_gen_with_pt") { 
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, sigma_positive_gen_Zmass_with_pt,new_pt_error_positive, sigma_positive_error_gen_Zmass_with_pt);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, sigma_negative_gen_Zmass_with_pt, new_pt_error_negative,  sigma_negative_error_gen_Zmass_with_pt);
    }

    if(type=="Zmass_reco_with_pt") { 
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, sigma_positive_Zmass_with_pt,new_pt_error_positive, sigma_positive_error_Zmass_with_pt);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, sigma_negative_Zmass_with_pt, new_pt_error_negative,  sigma_negative_error_Zmass_with_pt);
    }
    if(type=="Zmass_reco_FSR_with_pt") { 
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, sigma_positive_Zmass_FSR_with_pt,new_pt_error_positive, sigma_positive_error_Zmass_FSR_with_pt);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, sigma_negative_Zmass_FSR_with_pt, new_pt_error_negative,  sigma_negative_error_Zmass_FSR_with_pt);
    }
    
    gr1->SetTitle(string_title_positive_sigma);
    gr1->GetXaxis()->SetTitleOffset(1.2);
    gr1->GetYaxis()->SetTitleOffset(1.2);
    gr1->GetXaxis()->SetTitle("pT(GEN)");
    gr1->GetYaxis()->SetRangeUser(0, range_high_sigma);
    gr1->GetXaxis()->SetRangeUser(0,200);
    gr1->SetName(graph_name_positive_sigma);
    gr1->Draw("AP*");
    c1_graph->Update();
    c1_graph->SaveAs(saving_path+saving_name_positive_sigma+".pdf");
    c1_graph->Close();

    c2_graph = new TCanvas("c2_graph","canvas",1000,1000);
    c2_graph->cd();
    c2_graph->SetLeftMargin(0.15);
    c2_graph->SetTopMargin(0.2);
    gr2->SetTitle(string_title_negative_sigma);
    gr2->GetXaxis()->SetTitleOffset(1.2);
    gr2->GetYaxis()->SetTitleOffset(1.2);
    gr2->GetXaxis()->SetTitle("pT(GEN)");
    gr2->GetYaxis()->SetRangeUser(0,range_high_sigma);
    gr2->GetXaxis()->SetRangeUser(0,200);
    gr2->Draw("AP*");
    gr2->SetName(graph_name_negative_sigma);
    c2_graph->Update();
    c2_graph->SaveAs(saving_path+saving_name_negative_sigma+".pdf");
    c2_graph->Close(); 
  } // end of plotting mean delta with pT

 
   void class_reading::plotting_mean_delta_pT(TString type){
     bool diff=false;
     TString dir_name = "";
     TString graph_name_positive = "";
     TString saving_name_positive = "";

     TString graph_name_negative = "";
     TString saving_name_negative = "";
     TString graph_name_positive_sigma = "";

     TString graph_name_negative_sigma = "";

     TString saving_name_positive_sigma = "";
     TString saving_name_negative_sigma = "";

     TString text_file_name = ""; 
     TString string_title_positive = "";
     TString string_title_negative = "";
     TString string_title_positive_sigma = "";
     TString string_title_negative_sigma = "";

     TString title_multigraph_positive = "";
     TString title_multigraph_negative = "";
   //  TString string_title_positive = "";
   //  TString string_title_negative = "";

    TString saving_name_multigraph_positive_muon_positive_eta = "" ;
    TString saving_name_multigraph_positive_muon_negative_eta = ""  ;
    TString saving_name_multigraph_negative_muon_positive_eta = "" ; 
    TString saving_name_multigraph_negative_muon_negative_eta = "" ; 
      float range_low ;
      float range_high ; 
      float range_high_sigma ; 
      if(type=="delta_pT"){
       dir_name = "mean_resolution_delta_pT";
       graph_name_positive = "mean_delta_pT_positive_";
       graph_name_negative = "mean_delta_pT_negative_";
       graph_name_positive_sigma = "sigam_delta_pT_positive_";
       graph_name_negative_sigma = "sigam_delta_pT_negative_";

       saving_name_negative = "mean_sigma_delta_pT_plots/graph_resolution_delta_pT_mean_negative_eta_";
       saving_name_positive = "mean_sigma_delta_pT_plots/graph_resolution_delta_pT_mean_positive_eta_";
       saving_name_positive_sigma =  "mean_sigma_plots/graph_resolution_pT_sigma_positive_eta_";
       saving_name_negative_sigma =  "mean_sigma_plots/graph_resolution_pT_sigma_negative_eta_";

       text_file_name = "mean_sigma_delta_pT.txt";

//       string_title_mulitgraph_positive = "Mean #frac{#Delta pT}{pT} mu^{+}";
//       string_title_mulitgraph_negative = "Mean #frac{#Delta pT}{pT} mu^{-}";

       range_low = -0.01;
       range_high = 0.004;
       range_high_sigma = 0.07;

      title_multigraph_positive = "Mean #frac{#Delta pT}{pT} : #mu^{+}";
      title_multigraph_negative = "Mean #frac{#Delta pT}{pT} : #mu^{-}";

      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_delta_pT_comparison_plots/mean_delta_pT_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_delta_pT_comparison_plots/mean_delta_pT_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_delta_pT_comparison_plots/mean_delta_pT_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_delta_pT_comparison_plots/mean_delta_pT_negative_muon_negative_eta.pdf";
     
     }
      if(type=="Zmass_gen_with_pt"){
       dir_name = "mean_resolution_gen_Zmass_with_pt";
       graph_name_positive = "mean_gen_Zmass_positive_with_pt";
       graph_name_negative = "mean_gen_Zmass_negative_with_pt";
       graph_name_positive_sigma = "sigam_gen_Zmass_positive_with_pt";
       graph_name_negative_sigma = "sigam_gen_Zmass_negative_with_pt";

       saving_name_negative = "mean_sigma_Zmass_plots_with_pt/graph_resolution_gen_Zmass_mean_negative_with_pt_";
       saving_name_positive = "mean_sigma_Zmass_plots_with_pt/graph_resolution_gen_Zmass_mean_positive_with_pt";
       saving_name_positive_sigma =  "mean_sigma_plots_with_pt/graph_resolution_gen_Zmass_sigma_positive_with_pt";
       saving_name_negative_sigma =  "mean_sigma_plots_with_pt/graph_resolution_gen_Zmass_sigma_negative_with_pt";

       text_file_name = "mean_sigma_Zmass.txt";
 //      string_title_mulitgraph_positive = "Mean Z mass mu^{+}";
 //      string_title_mulitgraph_negative = "Mean Z mass mu^{-}";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;

      title_multigraph_positive = "Mean Z mass GEN : #mu^{+}";
      title_multigraph_negative = "Mean Z mass GEN : #mu^{-}";
    
      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_negative_muon_negative_eta.pdf";

     }

      if(type=="Zmass_gen"){
       dir_name = "mean_resolution_gen_Zmass";
       graph_name_positive = "mean_gen_Zmass_positive_";
       graph_name_negative = "mean_gen_Zmass_negative_";
       graph_name_positive_sigma = "sigam_gen_Zmass_positive_";
       graph_name_negative_sigma = "sigam_gen_Zmass_negative_";

       saving_name_negative = "mean_sigma_Zmass_plots/graph_resolution_gen_Zmass_mean_negative_eta_";
       saving_name_positive = "mean_sigma_Zmass_plots/graph_resolution_gen_Zmass_mean_positive_eta_";
       saving_name_positive_sigma =  "mean_sigma_plots/graph_resolution_gen_Zmass_sigma_positive_eta_";
       saving_name_negative_sigma =  "mean_sigma_plots/graph_resolution_gen_Zmass_sigma_negative_eta_";

       text_file_name = "mean_sigma_Zmass.txt";
 //      string_title_mulitgraph_positive = "Mean Z mass mu^{+}";
 //      string_title_mulitgraph_negative = "Mean Z mass mu^{-}";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;

      title_multigraph_positive = "Mean Z mass GEN : #mu^{+}";
      title_multigraph_negative = "Mean Z mass GEN : #mu^{-}";
    
      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_gen_Zmass_negative_muon_negative_eta.pdf";

     }


      if(type=="Zmass" && diff==false){
       dir_name = "mean_resolution_Zmass";
       graph_name_positive = "mean_Zmass_positive_";
       graph_name_negative = "mean_Zmass_negative_";
       graph_name_positive_sigma = "sigam_Zmass_positive_";
       graph_name_negative_sigma = "sigam_Zmass_negative_";

       saving_name_negative = "mean_sigma_Zmass_plots/graph_resolution_Zmass_mean_negative_eta_";
       saving_name_positive = "mean_sigma_Zmass_plots/graph_resolution_Zmass_mean_positive_eta_";
       saving_name_positive_sigma =  "mean_sigma_plots/graph_resolution_Zmass_sigma_positive_eta_";
       saving_name_negative_sigma =  "mean_sigma_plots/graph_resolution_Zmass_sigma_negative_eta_";

       text_file_name = "mean_sigma_Zmass.txt";
 //      string_title_mulitgraph_positive = "Mean Z mass mu^{+}";
 //      string_title_mulitgraph_negative = "Mean Z mass mu^{-}";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;

      title_multigraph_positive = "Mean Z mass : #mu^{+}";
      title_multigraph_negative = "Mean Z mass : #mu^{-}";
    
      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_negative_muon_negative_eta.pdf";

     }

      if(type=="Zmass_FSR" && diff==false){
       dir_name = "mean_resolution_Zmass_FSR";
       graph_name_positive = "mean_Zmass_FSR_positive_";
       graph_name_negative = "mean_Zmass_FSR_negative_";
       graph_name_positive_sigma = "sigam_Zmass_FSR_positive_";
       graph_name_negative_sigma = "sigam_Zmass_FSR_negative_";

       saving_name_negative = "mean_sigma_Zmass_FSR_plots/graph_resolution_Zmass_FSR_mean_negative_eta_";
       saving_name_positive = "mean_sigma_Zmass_FSR_plots/graph_resolution_Zmass_FSR_mean_positive_eta_";
       saving_name_positive_sigma =  "mean_sigma_plots/graph_resolution_Zmass_FSR_sigma_positive_eta_";
       saving_name_negative_sigma =  "mean_sigma_plots/graph_resolution_Zmass_FSR_sigma_negative_eta_";

       text_file_name = "mean_sigma_Zmass_FSR.txt";
 //      string_title_mulitgraph_positive = "Mean Z mass mu^{+}";
 //      string_title_mulitgraph_negative = "Mean Z mass mu^{-}";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 0.07;

      title_multigraph_positive = "Mean Z mass with FSR : #mu^{+}";
      title_multigraph_negative = "Mean Z mass with FSR : #mu^{-}";
    
      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_FSR_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_FSR_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_FSR_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_FSR_negative_muon_negative_eta.pdf";

     }
      if(type=="Zmass_FSR" && diff==true){
       dir_name = "mean_diff_resolution_Zmass_FSR";
       graph_name_positive = "mean_diff_Zmass_FSR_positive_";
       graph_name_negative = "mean_diff_Zmass_FSR_negative_";

       saving_name_negative = "mean_sigma_diff_Zmass_FSR_plots/graph_resolution_Zmass_FSR_mean_negative_eta_";
       saving_name_positive = "mean_sigma_diff_Zmass_FSR_plots/graph_resolution_Zmass_FSR_mean_positive_eta_";

       text_file_name = "mean_sigma_diff_Zmass_FSR.txt";

       range_high = 0.0005;
       range_low =-0.003;
       range_high_sigma = 0.07;

      title_multigraph_positive = "Mean Z mass with FSR : #mu^{+}";
      title_multigraph_negative = "Mean Z mass with FSR : #mu^{-}";
    
      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_FSR_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_FSR_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_FSR_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_FSR_negative_muon_negative_eta.pdf";

     }

      if(type=="Zmass" && diff==true){
       dir_name = "mean_diff_resolution_Zmass";
       graph_name_positive = "mean_diff_Zmass_positive_";
       graph_name_negative = "mean_diff_Zmass_negative_";

       saving_name_negative = "mean_sigma_diff_Zmass_plots/graph_resolution_Zmass_mean_negative_eta_";
       saving_name_positive = "mean_sigma_diff_mass_plots/graph_resolution_Zmass_mean_positive_eta_";

       text_file_name = "mean_sigma_diff_Zmass.txt";

       range_high = 0.0005;
       range_low = -.003;
       range_high_sigma = 0.07;

      title_multigraph_positive = "Mean Z mass : #mu^{+}";
      title_multigraph_negative = "Mean Z mass : #mu^{-}";
    
      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_diff_Zmass_comparison_plots/mean_Zmass_negative_muon_negative_eta.pdf";

     }

    ofstream sigma_text_file_positive;    
    ofstream sigma_text_file_negative;    
    sigma_text_file_positive.open("sigma_delta_pT_resolution_positive.txt");
    sigma_text_file_negative.open("sigma_delta_pT_resolution_negative.txt");
    TDirectory *dir_graph = (TDirectory*) output_root_file->mkdir(dir_name);
    dir_graph->cd();
    TGraphErrors *gr_mean1; 
    TGraphErrors *gr_mean2; 
    TGraphErrors *gr1; 
    TGraphErrors *gr2; 

    TCanvas *c1_mean_graph;
    TCanvas *c2_mean_graph;
    TCanvas *c1_graph;
    TCanvas *c2_graph;
 
    TString string_title;

    ofstream text_file;
    text_file.open("difference_Zmass.txt");
    std::cout<<"strating function plotting mean delta pT"<<std::endl; 
    // mean and sigma plots for different eta bins
    ofstream new_text;
    new_text.open(text_file_name);
   
/*    TMultiGraph *mg_positive =  new TMultiGraph();  
    TMultiGraph *mg_negative =  new TMultiGraph();  
    TMultiGraph *mg_sigma_positive =  new TMultiGraph();  
    TMultiGraph *mg_sigma_negative =  new TMultiGraph();  
*/
    TMultiGraph * mg_positive_muon_positive_eta = new TMultiGraph();
    TMultiGraph * mg_positive_muon_negative_eta = new TMultiGraph();
    TMultiGraph * mg_negative_muon_positive_eta = new TMultiGraph();
    TMultiGraph * mg_negative_muon_negative_eta = new TMultiGraph();

 //   TMultiGraph * mg_sigma_positive_muon_positive_eta = new TMultiGraph();
 //   TMultiGraph * mg_sigma_positive_muon_negative_eta = new TMultiGraph();
 //   TMultiGraph * mg_sigma_negative_muon_positive_eta = new TMultiGraph();
 //   TMultiGraph * mg_sigma_negative_muon_negative_eta = new TMultiGraph();
//     auto  * legend1_positive = new TLegend(0.15,0.15,0.30,0.35);
//     auto  * legend2_positive = new TLegend(0.15,0.15,0.30,0.35);

     auto  * legend1_positive = new TLegend(0.10,0.70,0.25,0.90);
     auto  * legend2_positive = new TLegend(0.10,0.70,0.25,0.90);
//     auto  * legend3_positive = new TLegend(0.13,0.65,0.25,0.85);
//     auto  * legend4_positive = new TLegend(0.13,0.65,0.25,0.85);

     auto  * legend1_negative = new TLegend(0.10,0.70,0.25,0.90);
     auto  * legend2_negative = new TLegend(0.10,0.70,0.25,0.90);

     auto  * legend_1_positive = new TLegend(0.10,0.10,0.25,0.30);
     auto  * legend_2_positive = new TLegend(0.10,0.10,0.25,0.30);

     auto  * legend_1_negative = new TLegend(0.10,0.10,0.25,0.30);
     auto  * legend_2_negative = new TLegend(0.10,0.10,0.25,0.30);

 //    auto  * legend3_negative = new TLegend(0.13,0.65,0.25,0.85);
 //    auto  * legend4_negative = new TLegend(0.13,0.65,0.25,0.85);

     TString title_legend; 
    
     int colour[7] = {1, 2, 3, 4, 6, 7, 28}; 
     int colour_2[7] = {28, 7, 6, 4, 3, 2, 1}; 

     sigma_text_file_positive<<"positive"<<std::endl;
     sigma_text_file_positive<<"row : eta bins"<<std::endl;
     sigma_text_file_positive<<"column : pT bins"<<std::endl;

     sigma_text_file_negative<<"positive"<<std::endl;
     sigma_text_file_negative<<"row : pT bins"<<std::endl;
     sigma_text_file_negative<<"column : eta bins"<<std::endl;

     // for separting each eta bins
      sigma_text_file_positive<<" pT : \t";
      sigma_text_file_negative<<" pT : \t";

     for(int i=0; i<pT_bins ; i++){
      sigma_text_file_positive<<" [ "<<pt_list[i]<<" - "<<pt_list[i+1]<<"] \t";
      sigma_text_file_negative<<" [ "<<pt_list[i]<<" - "<<pt_list[i+1]<<"] \t";
     }
       sigma_text_file_positive<<std::endl;
       sigma_text_file_negative<<std::endl;

     for(int j=0 ; j< eta_bins; j++){
      
      sigma_text_file_positive<<" eta : [ "<<eta_list[j]<<" - "<<eta_list[j+1]<<"] \t";
      sigma_text_file_negative<<" eta : [ "<<eta_list[j]<<" - "<<eta_list[j+1]<<"] \t";
      float new_pt_measured_positive[pT_bins] = {0}; 
      float new_pt_error_positive[pT_bins] = {0}; 
      float new_pt_measured_negative[pT_bins] = {0}; 
      float new_pt_error_negative[pT_bins] = {0}; 

      float diff_mean_positive_Zmass[pT_bins] = {0}; 
      float diff_mean_negative_Zmass[pT_bins] = {0}; 
      float diff_mean_error_positive_Zmass[pT_bins] = {0}; 
      float diff_mean_error_negative_Zmass[pT_bins] = {0}; 

      float diff_mean_positive_Zmass_FSR[pT_bins] = {0}; 
      float diff_mean_negative_Zmass_FSR[pT_bins] = {0}; 
      float diff_mean_error_positive_Zmass_FSR[pT_bins] = {0}; 
      float diff_mean_error_negative_Zmass_FSR[pT_bins] = {0}; 

      float new_mean_positive_delta_pT[pT_bins] = {0}; 
      float new_mean_negative_delta_pT[pT_bins] = {0}; 
      float new_sigma_positive_delta_pT[pT_bins] = {0}; 
      float new_sigma_negative_delta_pT[pT_bins] = {0}; 

      float new_mean_positive_error_delta_pT[pT_bins] = {0}; 
      float new_mean_negative_error_delta_pT[pT_bins] = {0}; 
      float new_sigma_positive_error_delta_pT[pT_bins] = {0}; 
      float new_sigma_negative_error_delta_pT[pT_bins] = {0}; 

      float new_mean_positive_gen_Zmass[pT_bins] = {0}; 
      float new_mean_negative_gen_Zmass[pT_bins] = {0}; 
      float new_sigma_positive_gen_Zmass[pT_bins] = {0}; 
      float new_sigma_negative_gen_Zmass[pT_bins] = {0}; 

      float new_mean_positive_error_gen_Zmass[pT_bins] = {0}; 
      float new_mean_negative_error_gen_Zmass[pT_bins] = {0}; 
      float new_sigma_positive_error_gen_Zmass[pT_bins] = {0}; 
      float new_sigma_negative_error_gen_Zmass[pT_bins] = {0}; 

      float new_mean_positive_Zmass[pT_bins] = {0}; 
      float new_mean_negative_Zmass[pT_bins] = {0}; 
      float new_sigma_positive_Zmass[pT_bins] = {0}; 
      float new_sigma_negative_Zmass[pT_bins] = {0}; 

      float new_mean_positive_error_Zmass[pT_bins] = {0}; 
      float new_mean_negative_error_Zmass[pT_bins] = {0}; 
      float new_sigma_positive_error_Zmass[pT_bins] = {0}; 
      float new_sigma_negative_error_Zmass[pT_bins] = {0}; 

      float new_mean_positive_Zmass_FSR[pT_bins] = {0}; 
      float new_mean_negative_Zmass_FSR[pT_bins] = {0}; 
      float new_sigma_positive_Zmass_FSR[pT_bins] = {0}; 
      float new_sigma_negative_Zmass_FSR[pT_bins] = {0}; 

      float new_mean_positive_error_Zmass_FSR[pT_bins] = {0}; 
      float new_mean_negative_error_Zmass_FSR[pT_bins] = {0}; 
      float new_sigma_positive_error_Zmass_FSR[pT_bins] = {0}; 
      float new_sigma_negative_error_Zmass_FSR[pT_bins] = {0}; 


     if(type=="delta_pT") {
      string_title_positive = TString::Format("Mean #frac{#Delta pT}{pT} (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
      string_title_negative = TString::Format("Mean #frac{#Delta pT}{pT} (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
      string_title_positive_sigma = TString::Format("#sigma : #frac{#Delta pT}{pT}(#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
      string_title_positive_sigma = TString::Format("#sigma : #frac{#Delta pT}{pT} (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
     }

     if(type=="Zmass_gen") {
       string_title_positive = TString::Format("Mean Zmass GEN (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_negative = TString::Format("Mean Zmass GEN (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_positive_sigma = TString::Format("#sigma : Z mass GEN (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_positive_sigma = TString::Format("#sigma : Z mass GEN (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);

       string_title_negative_sigma = TString::Format("#sigma : Z mass GEN (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_negative_sigma = TString::Format("#sigma : Z mass GEN (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
     }


     if(type=="Zmass") {
       string_title_positive = TString::Format("Mean Zmass (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_negative = TString::Format("Mean Zmass (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_positive_sigma = TString::Format("#sigma : Z mass (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_positive_sigma = TString::Format("#sigma : Z mass (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);

       string_title_negative_sigma = TString::Format("#sigma : Z mass (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_negative_sigma = TString::Format("#sigma : Z mass (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
     }
     if(type=="Zmass_FSR") {
       string_title_positive = TString::Format("Mean Zmass with FSR (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_negative = TString::Format("Mean Zmass with FSR (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_positive_sigma = TString::Format("#sigma : Z mass with FSR (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_positive_sigma = TString::Format("#sigma : Z mass with FSR (#mu^{+}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);

       string_title_negative_sigma = TString::Format("#sigma : Z mass with FSR (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
       string_title_negative_sigma = TString::Format("#sigma : Z mass with FSR (#mu^{-}) : %.2f  #leq #eta < %.2f", eta_list[j], eta_list[j+1]);
     }

  
       for(int i=0; i<pT_bins; i++){
        std::cout<<"first bin of eta : pT mean value  "<<pt_measured_positive[i][j]<<std::endl; 
        new_pt_measured_positive[i]= pt_measured_positive[i][j];
        new_pt_error_positive[i]= pt_error_positive[i][j];
        new_pt_measured_negative[i]= pt_measured_negative[i][j];
        new_pt_error_negative[i]= pt_error_negative[i][j];
//        std::cout<<"first bin of eta : deltapT mean "<<mean_positive_delta_pT[i][j]<<std::endl; 

       if(type=="delta_pT"){
        new_mean_positive_delta_pT[i]= mean_positive_delta_pT[i][j];
        new_mean_negative_delta_pT[i]= mean_negative_delta_pT[i][j];
        new_mean_positive_error_delta_pT[i]= mean_positive_error_delta_pT[i][j];
        new_mean_negative_error_delta_pT[i]= mean_negative_error_delta_pT[i][j];
        new_sigma_positive_delta_pT[i]= sigma_positive_delta_pT[i][j];
        new_sigma_negative_delta_pT[i]= sigma_negative_delta_pT[i][j];
        new_sigma_positive_error_delta_pT[i]= sigma_positive_error_delta_pT[i][j];
        new_sigma_negative_error_delta_pT[i]= sigma_negative_error_delta_pT[i][j];

        new_text<<"positive muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_positive_delta_pT[i]<<" pt "<<new_pt_measured_positive[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_positive_delta_pT[i]<<" pt error "<<new_pt_error_positive[i]<<std::endl;

        new_text<<"negative muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_negative_delta_pT[i]<<" pt "<<new_pt_measured_negative[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_negative_delta_pT[i]<<" pt error "<<new_pt_error_negative[i]<<std::endl;
         
        sigma_text_file_positive<<sigma_positive_delta_pT[i][j]<<"\t";
        sigma_text_file_negative<<sigma_negative_delta_pT[i][j]<<"\t";

       }
       if(type=="Zmass_gen" && diff==false){ 
        new_mean_positive_gen_Zmass[i]= mean_positive_gen_Zmass[i][j];
        new_mean_negative_gen_Zmass[i]= mean_negative_gen_Zmass[i][j];
        new_mean_positive_error_gen_Zmass[i]= mean_positive_error_gen_Zmass[i][j];
        new_mean_negative_error_gen_Zmass[i]= mean_negative_error_gen_Zmass[i][j];
        new_sigma_positive_gen_Zmass[i]= sigma_positive_gen_Zmass[i][j];
        new_sigma_negative_gen_Zmass[i]= sigma_negative_gen_Zmass[i][j];
        new_sigma_positive_error_gen_Zmass[i]= sigma_positive_error_gen_Zmass[i][j];
        new_sigma_negative_error_gen_Zmass[i]= sigma_negative_error_gen_Zmass[i][j];

        new_text<<"positive muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_positive_gen_Zmass[i]<<" pt "<<new_pt_measured_positive[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_positive_gen_Zmass[i]<<" pt error "<<new_pt_error_positive[i]<<std::endl;

        new_text<<"negative muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_negative_gen_Zmass[i]<<" pt "<<new_pt_measured_negative[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_negative_gen_Zmass[i]<<" pt error "<<new_pt_error_negative[i]<<std::endl;
       }


       if(type=="Zmass" && diff==false){ 
        new_mean_positive_Zmass[i]= mean_positive_Zmass[i][j];
        new_mean_negative_Zmass[i]= mean_negative_Zmass[i][j];
        new_mean_positive_error_Zmass[i]= mean_positive_error_Zmass[i][j];
        new_mean_negative_error_Zmass[i]= mean_negative_error_Zmass[i][j];
        new_sigma_positive_Zmass[i]= sigma_positive_Zmass[i][j];
        new_sigma_negative_Zmass[i]= sigma_negative_Zmass[i][j];
        new_sigma_positive_error_Zmass[i]= sigma_positive_error_Zmass[i][j];
        new_sigma_negative_error_Zmass[i]= sigma_negative_error_Zmass[i][j];

        new_text<<"positive muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_positive_Zmass[i]<<" pt "<<new_pt_measured_positive[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_positive_Zmass[i]<<" pt error "<<new_pt_error_positive[i]<<std::endl;

        new_text<<"negative muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_negative_Zmass[i]<<" pt "<<new_pt_measured_negative[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_negative_Zmass[i]<<" pt error "<<new_pt_error_negative[i]<<std::endl;
       }


       if(type=="Zmass_FSR" && diff==false){ 
        new_mean_positive_Zmass_FSR[i]= mean_positive_Zmass_FSR[i][j];
        new_mean_negative_Zmass_FSR[i]= mean_negative_Zmass_FSR[i][j];
        new_mean_positive_error_Zmass_FSR[i]= mean_positive_error_Zmass_FSR[i][j];
        new_mean_negative_error_Zmass_FSR[i]= mean_negative_error_Zmass_FSR[i][j];
        new_sigma_positive_Zmass_FSR[i]= sigma_positive_Zmass_FSR[i][j];
        new_sigma_negative_Zmass_FSR[i]= sigma_negative_Zmass_FSR[i][j];
        new_sigma_positive_error_Zmass_FSR[i]= sigma_positive_error_Zmass_FSR[i][j];
        new_sigma_negative_error_Zmass_FSR[i]= sigma_negative_error_Zmass_FSR[i][j];

        new_text<<"positive muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_positive_Zmass_FSR[i]<<" pt "<<new_pt_measured_positive[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_positive_Zmass_FSR[i]<<" pt error "<<new_pt_error_positive[i]<<std::endl;

        new_text<<"negative muon : etabin :  "<<j<<std::endl;
        new_text<<"new mean "<<new_mean_negative_Zmass_FSR[i]<<" pt "<<new_pt_measured_negative[i]<<std::endl;
        new_text<<"new sigma "<<new_sigma_negative_Zmass_FSR[i]<<" pt error "<<new_pt_error_negative[i]<<std::endl;
       }
        
       // evaluating difference here
       if(type=="Zmass" && diff==true){
        new_mean_positive_Zmass[i]= mean_positive_Zmass[i][j];
        new_mean_negative_Zmass[i]= mean_negative_Zmass[i][j];
        new_mean_positive_error_Zmass[i]= mean_positive_error_Zmass[i][j];
        new_mean_negative_error_Zmass[i]= mean_negative_error_Zmass[i][j];
        new_sigma_positive_Zmass[i]= sigma_positive_Zmass[i][j];
        new_sigma_negative_Zmass[i]= sigma_negative_Zmass[i][j];
        new_sigma_positive_error_Zmass[i]= sigma_positive_error_Zmass[i][j];
        new_sigma_negative_error_Zmass[i]= sigma_negative_error_Zmass[i][j];

        new_mean_positive_gen_Zmass[i]= mean_positive_gen_Zmass[i][j];
        new_mean_negative_gen_Zmass[i]= mean_negative_gen_Zmass[i][j];
        new_mean_positive_error_gen_Zmass[i]= mean_positive_error_gen_Zmass[i][j];
        new_mean_negative_error_gen_Zmass[i]= mean_negative_error_gen_Zmass[i][j];
        new_sigma_positive_gen_Zmass[i]= sigma_positive_gen_Zmass[i][j];
        new_sigma_negative_gen_Zmass[i]= sigma_negative_gen_Zmass[i][j];
        new_sigma_positive_error_gen_Zmass[i]= sigma_positive_error_gen_Zmass[i][j];
        new_sigma_negative_error_gen_Zmass[i]= sigma_negative_error_gen_Zmass[i][j];

        text_file<<" mean  Z mass posiitive"<<new_mean_positive_Zmass[i]<<" gen "<<new_mean_positive_gen_Zmass[i]<<std::endl;
        text_file<<" mean  Z mass negative"<<new_mean_negative_Zmass[i]<<" gen "<<new_mean_negative_gen_Zmass[i]<<std::endl;
        diff_mean_positive_Zmass[i]= ((new_mean_positive_Zmass[i]- new_mean_positive_gen_Zmass[i])/new_mean_positive_gen_Zmass[i]);
        diff_mean_negative_Zmass[i]= ((new_mean_negative_Zmass[i]- new_mean_negative_gen_Zmass[i])/new_mean_negative_gen_Zmass[i]);

        text_file<<" diff Z mass posiitive"<<diff_mean_positive_Zmass[i]<<std::endl;
        text_file<<" diff Z mass negative"<<diff_mean_negative_Zmass[i]<<std::endl;

        diff_mean_error_positive_Zmass[i]= diff_mean_positive_Zmass[i] * TMath::Sqrt( pow((new_mean_positive_error_Zmass[i]/new_mean_positive_Zmass[i]) ,2) + pow(( new_mean_positive_error_gen_Zmass[i]/new_mean_positive_gen_Zmass[i]),2) );
        diff_mean_error_negative_Zmass[i]= diff_mean_negative_Zmass[i] * TMath::Sqrt( pow((new_mean_negative_error_Zmass[i]/new_mean_negative_Zmass[i]) ,2) + pow(( new_mean_negative_error_gen_Zmass[i]/new_mean_negative_gen_Zmass[i]),2) );
       }
       if(type=="Zmass_FSR" && diff==true){
        new_mean_positive_Zmass_FSR[i]= mean_positive_Zmass_FSR[i][j];
        new_mean_negative_Zmass_FSR[i]= mean_negative_Zmass_FSR[i][j];
        new_mean_positive_error_Zmass_FSR[i]= mean_positive_error_Zmass_FSR[i][j];
        new_mean_negative_error_Zmass_FSR[i]= mean_negative_error_Zmass_FSR[i][j];
        new_sigma_positive_Zmass_FSR[i]= sigma_positive_Zmass_FSR[i][j];
        new_sigma_negative_Zmass_FSR[i]= sigma_negative_Zmass_FSR[i][j];
        new_sigma_positive_error_Zmass_FSR[i]= sigma_positive_error_Zmass_FSR[i][j];
        new_sigma_negative_error_Zmass_FSR[i]= sigma_negative_error_Zmass_FSR[i][j];

        new_mean_positive_gen_Zmass[i]= mean_positive_gen_Zmass[i][j];
        new_mean_negative_gen_Zmass[i]= mean_negative_gen_Zmass[i][j];
        new_mean_positive_error_gen_Zmass[i]= mean_positive_error_gen_Zmass[i][j];
        new_mean_negative_error_gen_Zmass[i]= mean_negative_error_gen_Zmass[i][j];
        new_sigma_positive_gen_Zmass[i]= sigma_positive_gen_Zmass[i][j];
        new_sigma_negative_gen_Zmass[i]= sigma_negative_gen_Zmass[i][j];
        new_sigma_positive_error_gen_Zmass[i]= sigma_positive_error_gen_Zmass[i][j];
        new_sigma_negative_error_gen_Zmass[i]= sigma_negative_error_gen_Zmass[i][j];


        diff_mean_positive_Zmass_FSR[i]= ((new_mean_positive_Zmass_FSR[i]- new_mean_positive_gen_Zmass[i])/new_mean_positive_gen_Zmass[i]);
        diff_mean_negative_Zmass_FSR[i]= ((new_mean_negative_Zmass_FSR[i]- new_mean_negative_gen_Zmass[i])/new_mean_negative_gen_Zmass[i]);

        text_file<<" diff Z mass FSR posiitive"<<diff_mean_positive_Zmass[i];
        text_file<<" diff Z mass FSR negative"<<diff_mean_negative_Zmass[i];
       diff_mean_error_positive_Zmass_FSR[i]= diff_mean_positive_Zmass_FSR[i] * TMath::Sqrt( pow((new_mean_positive_error_Zmass_FSR[i]/new_mean_positive_Zmass_FSR[i]) ,2) + pow(( new_mean_positive_error_gen_Zmass[i]/new_mean_positive_gen_Zmass[i]),2) );
        diff_mean_error_negative_Zmass_FSR[i]= diff_mean_negative_Zmass_FSR[i] * TMath::Sqrt( pow((new_mean_negative_error_Zmass_FSR[i]/new_mean_negative_Zmass_FSR[i]) ,2) + pow(( new_mean_negative_error_gen_Zmass[i]/new_mean_negative_gen_Zmass[i]),2) );
       }

      } // end of filling the mean and sigma values

      // time to plot these mean and sigma values
    
       std::cout<<"started graph plotting "<<std::endl; 

       if(type=="delta_pT") {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_mean_positive_delta_pT, new_pt_error_positive, new_mean_positive_error_delta_pT);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_mean_negative_delta_pT, new_pt_error_negative, new_mean_negative_error_delta_pT);
       sigma_text_file_positive<<std::endl;
       sigma_text_file_negative<<std::endl;

       }
       if(type=="Zmass_gen") {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_mean_positive_gen_Zmass,new_pt_error_positive, new_mean_positive_error_gen_Zmass);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_mean_negative_gen_Zmass, new_pt_error_negative, new_mean_negative_error_gen_Zmass);
       }
       if(type=="Zmass") {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_mean_positive_Zmass,new_pt_error_positive, new_mean_positive_error_Zmass);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_mean_negative_Zmass, new_pt_error_negative, new_mean_negative_error_Zmass);
       }
       if(type=="Zmass_FSR") {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_mean_positive_Zmass_FSR,new_pt_error_positive, new_mean_positive_error_Zmass_FSR);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_mean_negative_Zmass_FSR, new_pt_error_negative, new_mean_negative_error_Zmass_FSR);
       } 
       if(type=="Zmass_FSR" && diff==true) {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, diff_mean_positive_Zmass_FSR,new_pt_error_positive, diff_mean_error_positive_Zmass_FSR);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, diff_mean_negative_Zmass_FSR, new_pt_error_negative, diff_mean_error_negative_Zmass_FSR);
       }
       if(type=="Zmass" && diff==true) {
        gr_mean1 = new TGraphErrors(pT_bins, new_pt_measured_positive, diff_mean_positive_Zmass,new_pt_error_positive, diff_mean_error_positive_Zmass);
        gr_mean2 = new TGraphErrors(pT_bins, new_pt_measured_negative, diff_mean_negative_Zmass, new_pt_error_negative, diff_mean_error_negative_Zmass);
       }

    c1_mean_graph = new TCanvas("c1_mean_graph","canvas",1000,1000);
    c1_mean_graph->cd();
    c1_mean_graph->SetLeftMargin(0.15);
    c1_mean_graph->SetTopMargin(0.2);
    gr_mean1->SetTitle(string_title_positive);
    gr_mean1->GetXaxis()->SetTitleOffset(1.2);
    gr_mean1->GetYaxis()->SetTitleOffset(1.2);
    gr_mean1->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean1->GetYaxis()->SetRangeUser(range_low, range_high);
    gr_mean1->Draw("AP*");
    gr_mean1->SetName(graph_name_positive+bin_number_eta[j]);
    gr_mean1->GetXaxis()->SetRangeUser(0,200);
    c1_mean_graph->Update();
   // c1_mean_graph->SaveAs(saving_path+"mean_sigma_plots/graph_resolution_pT_mean_positive_eta_"+bin_number_eta[j]+".pdf");
    c1_mean_graph->SaveAs(saving_path+saving_name_positive+bin_number_eta[j]+".pdf");
    c1_mean_graph->Close();

    // making an overlapping TGraph for comparing variation in eta bins 


    std::cout<<"ended 1st  graph plotting "<<std::endl; 
    //c2_mean_graph = new TCanvas("c2_mean_graph","Mean #frac{#Delta pT}{pT} mu^{-}",1000,1000);
    c2_mean_graph = new TCanvas("c2_mean_graph","canvas",1000,1000);
    c2_mean_graph->cd();
    c2_mean_graph->SetLeftMargin(0.15);
    c2_mean_graph->SetTopMargin(0.2);
    gr_mean2->SetTitle(string_title_negative);
    gr_mean2->GetXaxis()->SetTitleOffset(1.2);
    gr_mean2->GetYaxis()->SetTitleOffset(1.2);
    gr_mean2->GetXaxis()->SetTitle("pT(GEN)");
    gr_mean2->GetYaxis()->SetRangeUser(range_low,range_high);
    gr_mean2->GetXaxis()->SetRangeUser(0,200);
    //gr_mean2->SetName("mean_delta_pT_negative_"+bin_number_eta[j]);
    gr_mean2->SetName(graph_name_negative+bin_number_eta[j]);
    gr_mean2->Draw("AP*");
    c2_mean_graph->Update();
    //c2_mean_graph->SaveAs(saving_path+"mean_sigma_plots/graph_resolution_pT_mean_negative_eta_"+bin_number_eta[j]+".pdf");
    c2_mean_graph->SaveAs(saving_path+saving_name_negative+bin_number_eta[j]+".pdf");
    c2_mean_graph->Close();

    if(j<7){
    gr_mean1->SetMarkerColor(colour[j]);
    gr_mean2->SetMarkerColor(colour[j]);
    }

    else { 
    gr_mean1->SetMarkerColor(colour_2[j-7]);
    gr_mean2->SetMarkerColor(colour_2[j-7]);
    }

    //mg_positive->Add(gr_mean1);
    //mg_negative->Add(gr_mean2);

    //title_legend = TString::Format("%0.0f #leq #eta < %0.0f", eta_list_symbol[j],eta_list_symbol[j+1]);
if(diff==false){ 
    c1_graph = new TCanvas("c1_graph","canvas",1000,1000);
    c1_graph->cd();
    c1_graph->SetLeftMargin(0.15);
    c1_graph->SetTopMargin(0.2);
    if(type=="delta_pT") { 
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_sigma_positive_delta_pT,new_pt_error_positive, new_sigma_positive_error_delta_pT);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_sigma_negative_delta_pT, new_pt_error_negative,  new_sigma_negative_error_delta_pT);
    }
    if(type=="Zmass_gen") { 
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_sigma_positive_gen_Zmass,new_pt_error_positive, new_sigma_positive_error_gen_Zmass);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_sigma_negative_gen_Zmass, new_pt_error_negative,  new_sigma_negative_error_gen_Zmass);
    }

    if(type=="Zmass") { 
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_sigma_positive_Zmass,new_pt_error_positive, new_sigma_positive_error_Zmass);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_sigma_negative_Zmass, new_pt_error_negative,  new_sigma_negative_error_Zmass);
    }
    if(type=="Zmass_FSR") { 
    gr1 = new TGraphErrors(pT_bins, new_pt_measured_positive, new_sigma_positive_Zmass_FSR,new_pt_error_positive, new_sigma_positive_error_Zmass_FSR);
    gr2 = new TGraphErrors(pT_bins, new_pt_measured_negative, new_sigma_negative_Zmass_FSR, new_pt_error_negative,  new_sigma_negative_error_Zmass_FSR);
    }
    
    gr1->SetTitle(string_title_positive_sigma);
    gr1->GetXaxis()->SetTitleOffset(1.2);
    gr1->GetYaxis()->SetTitleOffset(1.2);
    gr1->GetXaxis()->SetTitle("pT(GEN)");
    gr1->GetYaxis()->SetRangeUser(0, range_high_sigma);
    gr1->GetXaxis()->SetRangeUser(0,200);
    gr1->SetName(graph_name_positive_sigma+bin_number_eta[j]);
    gr1->Draw("AP*");
    c1_graph->Update();
    c1_graph->SaveAs(saving_path+saving_name_positive_sigma+bin_number_eta[j]+".pdf");
    c1_graph->Close();

    c2_graph = new TCanvas("c2_graph","canvas",1000,1000);
    c2_graph->cd();
    c2_graph->SetLeftMargin(0.15);
    c2_graph->SetTopMargin(0.2);
    gr2->SetTitle(string_title_negative_sigma);
    gr2->GetXaxis()->SetTitleOffset(1.2);
    gr2->GetYaxis()->SetTitleOffset(1.2);
    gr2->GetXaxis()->SetTitle("pT(GEN)");
    gr2->GetYaxis()->SetRangeUser(0,range_high_sigma);
    gr2->GetXaxis()->SetRangeUser(0,200);
    gr2->Draw("AP*");
    gr2->SetName(graph_name_negative_sigma+bin_number_eta[j]);
    c2_graph->Update();
    c2_graph->SaveAs(saving_path+saving_name_negative_sigma+bin_number_eta[j]+".pdf");
    c2_graph->Close(); 

    if(j<7){
    gr1->SetMarkerColor(colour[j]);
    gr2->SetMarkerColor(colour[j]);
    }
    else {
    gr1->SetMarkerColor(colour_2[j-7]);
    gr2->SetMarkerColor(colour_2[j-7]);
    }
   
    gr1->Write();
    gr2->Write();
}
    gr_mean1->Write();
    gr_mean2->Write();
    
    title_legend =   eta_list_symbol[j] + " #leq #eta < "+eta_list_symbol[j+1];
    if(j<7) {  
       mg_positive_muon_negative_eta->Add(gr_mean1);
       mg_negative_muon_negative_eta->Add(gr_mean2);
//       mg_sigma_positive_muon_negative_eta->Add(gr1);
//       mg_sigma_negative_muon_negative_eta->Add(gr2);

       legend1_negative->AddEntry(gr_mean1,title_legend,"P");
       legend2_negative->AddEntry(gr_mean2,title_legend,"P");
 //      legend3_negative->AddEntry(gr1,title_legend,"P");
 //      legend4_negative->AddEntry(gr2,title_legend,"P");
    }

    else  {     
      mg_positive_muon_positive_eta->Add(gr_mean1);
      mg_negative_muon_positive_eta->Add(gr_mean2);
  //    mg_sigma_positive_muon_positive_eta->Add(gr1);
  //    mg_sigma_negative_muon_positive_eta->Add(gr2);

      legend1_positive->AddEntry(gr_mean1,title_legend,"P");
      legend2_positive->AddEntry(gr_mean2,title_legend,"P");
   //   legend3_positive->AddEntry(gr1,title_legend,"P");
   //   legend4_positive->AddEntry(gr2,title_legend,"P");
    }

   
	} // end of eta bins and plots
 
   TCanvas *mg_first_positive_1 = new TCanvas("mg_first_positive_1","canvas" ,1000,1000); 
   TCanvas *mg_first_negative_1 = new TCanvas("mg_first_negative_1","canvas",1000,1000); 
 //  TCanvas *mg_second_positive_1 = new TCanvas("mg_second_positive_1", "#sigma #frac{#Delta pT}{pT} mu^{+}",1000,1000); 
//   TCanvas *mg_second_negative_1 = new TCanvas("mg_second_negative_1", "#sigma #frac{#Delta pT}{pT} mu^{-}",1000,1000); 

   TCanvas *mg_first_positive_2 = new TCanvas("mg_first_positive_2", "canvas",1000,1000); 
   TCanvas *mg_first_negative_2 = new TCanvas("mg_first_negative_2", "canvas",1000,1000); 
  // TCanvas *mg_second_positive_2 = new TCanvas("mg_second_positive_2", "#sigma #frac{#Delta pT}{pT} mu^{+}",1000,1000); 
   //TCanvas *mg_second_negative_2 = new TCanvas("mg_second_negative_2", "#sigma #frac{#Delta pT}{pT} mu^{-}",1000,1000); 

   mg_first_positive_1->cd();
   mg_positive_muon_positive_eta->Draw("AP");
   //mg_positive_muon_positive_eta->SetTitle("Mean (#frac{#Delta pT}{pT}) : #mu^{+} ("+corrections+" rochester corrections)");
   mg_positive_muon_positive_eta->SetTitle(title_multigraph_positive+"("+corrections+" rochester corrections)");
   mg_positive_muon_positive_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_positive_muon_positive_eta->GetXaxis()->SetRangeUser(0,200);
   mg_positive_muon_positive_eta->GetYaxis()->SetRangeUser(range_low,range_high);
   if(type=="Zmass_gen") legend_1_positive->Draw();
   else legend1_positive->Draw();
   mg_first_positive_1->SaveAs(saving_path+saving_name_multigraph_positive_muon_positive_eta+".pdf");

   mg_first_negative_1->cd();
   mg_negative_muon_positive_eta->Draw("AP");
   //mg_negative_muon_positive_eta->SetTitle("Mean (#frac{#Delta pT}{pT}) : #mu^{-} ("+corrections+" rochester corrections)");
   mg_negative_muon_positive_eta->SetTitle(title_multigraph_negative+" ("+corrections+" rochester corrections)");
   mg_negative_muon_positive_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_negative_muon_positive_eta->GetXaxis()->SetRangeUser(0,200);
   mg_negative_muon_positive_eta->GetYaxis()->SetRangeUser(range_low,range_high);
   if(type=="Zmass_gen") legend_2_positive->Draw();
   else legend2_positive->Draw();
   mg_first_negative_1->SaveAs(saving_path+saving_name_multigraph_negative_muon_positive_eta+".pdf");

/*   mg_second_positive_1->cd();
   mg_sigma_positive_muon_positive_eta->Draw("AP");
   mg_sigma_positive_muon_positive_eta->SetTitle("#sigma (#frac{#Delta pT}{pT}) : #mu^{+} ("+corrections+" rochester corrections)");
   mg_sigma_positive_muon_positive_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_sigma_positive_muon_positive_eta->GetXaxis()->SetRangeUser(0,200);
   mg_sigma_positive_muon_positive_eta->GetYaxis()->SetRangeUser(0,0.07);
   legend3_positive->Draw();
   mg_second_positive_1->SaveAs(saving_path+"mean_sigma_comparison_plots/sigma_delta_pT_positive_muon_positive_eta.pdf");

   mg_second_negative_1->cd();
   mg_sigma_negative_muon_positive_eta->Draw("AP");
   mg_sigma_negative_muon_positive_eta->SetTitle("#sigma (#frac{#Delta pT}{pT}) : #mu^{-} ("+corrections+" rochester corrections)");
   mg_sigma_negative_muon_positive_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_sigma_negative_muon_positive_eta->GetXaxis()->SetRangeUser(0,200);
   mg_sigma_negative_muon_positive_eta->GetYaxis()->SetRangeUser(0.0,0.07);
   legend4_positive->Draw();
   mg_second_negative_1->SaveAs(saving_path+"mean_sigma_comparison_plots/sigma_delta_pT_negative_muon_positive_eta.pdf");
 */ 
   //negative eta bins
   mg_first_positive_2->cd();
   mg_positive_muon_negative_eta->Draw("AP");
   mg_positive_muon_negative_eta->SetTitle(title_multigraph_positive+" ("+corrections+" rochester corrections)");
   mg_positive_muon_negative_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_positive_muon_negative_eta->GetXaxis()->SetRangeUser(0,200);
   mg_positive_muon_negative_eta->GetYaxis()->SetRangeUser(range_low,range_high);
   if(type=="Zmass_gen") legend_1_negative->Draw();
   else  legend1_negative->Draw();

   mg_first_positive_2->SaveAs(saving_path+saving_name_multigraph_positive_muon_negative_eta+".pdf");

   mg_first_negative_2->cd();
   mg_negative_muon_negative_eta->Draw("AP");
   mg_negative_muon_negative_eta->SetTitle(title_multigraph_negative+corrections+" rochester corrections)");
   mg_negative_muon_negative_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_negative_muon_negative_eta->GetXaxis()->SetRangeUser(0,200);
   mg_negative_muon_negative_eta->GetYaxis()->SetRangeUser(range_low,range_high);
   if(type=="Zmass_gen") legend_2_negative->Draw();
   else legend2_negative->Draw();
   mg_first_negative_2->SaveAs(saving_path+ saving_name_multigraph_negative_muon_negative_eta +".pdf");

/*   mg_second_positive_2->cd();
   mg_sigma_positive_muon_negative_eta->Draw("AP");
   mg_sigma_positive_muon_negative_eta->SetTitle("#sigma (#frac{#Delta pT}{pT}) : #mu^{+} ("+corrections+" rochester corrections)");
   mg_sigma_positive_muon_negative_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_sigma_positive_muon_negative_eta->GetXaxis()->SetRangeUser(0,200);
   mg_sigma_positive_muon_negative_eta->GetYaxis()->SetRangeUser(0,0.07);
   legend3_negative->Draw();
   mg_second_positive_2->SaveAs(saving_path+"mean_sigma_comparison_plots/sigma_delta_pT_positive_muon_negative_eta.pdf");

   mg_second_negative_2->cd();
   mg_sigma_negative_muon_negative_eta->Draw("AP");
   mg_sigma_negative_muon_negative_eta->SetTitle("#sigma (#frac{#Delta pT}{pT}) : #mu^{-} ("+corrections+" rochester corrections)");
   mg_sigma_negative_muon_negative_eta->GetXaxis()->SetTitle("pT (GeV) ");
   mg_sigma_negative_muon_negative_eta->GetXaxis()->SetRangeUser(0,200);
   mg_sigma_negative_muon_negative_eta->GetYaxis()->SetRangeUser(0.0,0.07);
   legend4_negative->Draw();
   mg_second_negative_2->SaveAs(saving_path+"mean_sigma_comparison_plots/sigma_delta_pT_negative_muon_negative_eta.pdf");
*/
  }


    void class_reading :: fitting_BW(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into Breit Wigner
     RooRealVar mass_var("mass_var","mass_var",65,110);
     RooDataHist histo("histo","mass dataset",mass_var,hist_fit);
     RooRealVar mean_mass("mean_mass","mean of Z mass",91,87,95);
     RooRealVar width("width","width of Z mass",2.49,1,5);
   
     RooPlot *xframe=mass_var.frame(Title(title_name));
     histo.plotOn(xframe);
//     histo.statOn(xframe);
   
     RooBreitWigner BW("BW","Breit Wigner fit",mass_var, mean_mass,width);
     BW.fitTo(histo,Range(70,110));
    
     BW.plotOn(xframe,RooFit::LineColor(kRed+2),Name("BW_sig"));
    
     BW.paramOn(xframe,RooFit::Layout(0.6,0.9,0.7));
     
     TCanvas *tmp = new TCanvas("tmp","Gen Z mass", 900,600);
     tmp->cd();
     gPad->SetLeftMargin(0.15);
     xframe->getAttText()->SetTextSize(0.025);
     xframe->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
   
     xframe->GetYaxis()->SetTitle("N/0.5 (GeV)");
     xframe->Draw();
     float chi_square_value = xframe->chiSquare();
     TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen Z", "EP");
     leg2->AddEntry("BW_sig","BW fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       

  
     gStyle->SetOptStat();
     tmp->SaveAs(saving_path+ saving_name+".pdf");
     tmp->Close(); 
     std::cout<<" Z mass mean : "<<mean_mass.getVal()<<" width : "<<width.getVal()<<std::endl;
    
/*    if(type=="positive"){ 
     mean_positive_gen_Zmass[pt_bin_number][eta_bin_number] = mean_mass.getVal();
     mean_positive_error_gen_Zmass[pt_bin_number][eta_bin_number] = mean_mass.getError();
     chi2_ndf_positive_gen_Zmass[pt_bin_number][eta_bin_number] = chi_square_value ;
     }
    if(type=="negative"){ 
     mean_negative_gen_Zmass[pt_bin_number][eta_bin_number] = mean_mass.getVal();
     mean_negative_error_gen_Zmass[pt_bin_number][eta_bin_number] = mean_mass.getError();
     chi2_ndf_negative_gen_Zmass[pt_bin_number][eta_bin_number] =  chi_square_value;
     }
*/
  }

 
    // implementing DSCB fit
    void class_reading :: fitting_DSCB(TH1F * hist_fit, TString saving_name, TString title_name, int pt_bin_number, int eta_bin_number, TString type){
      TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
      TF1* f1 = new TF1("f1", "gaus", -0.02,0.02);
      hist_fit_clone->Fit(f1,"R+");

      float mean_gaus =0; float sigma_gaus = 0; 
      mean_gaus = f1->GetParameter(1);
      sigma_gaus = f1->GetParameter(2);
 
     RooRealVar delta_pT("delta_pT", "deltapT", -0.2,0.2, "");
     RooDataHist histo("histo","dataset with var",delta_pT,hist_fit);
     RooRealVar Mean("Mean", "Mean",mean_gaus, -0.02, 0.02);
     RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.0001,0.1);//sigma[decay]);
     RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.0001, 500);//alphaL[decay]);
     RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 500);//expL[decay]);
     RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.0001, 500);//alphaR[decay]);
     RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 500);//expR[decay]);
     RooMyPDF_DSCB DSCB("DSCB", "DSCB", delta_pT, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
 
      output_file<<"fitting positive muon : mean "<<mean_gaus<<" sigma : "<<sigma_gaus<<std::endl;
     // Fit the mass into DSCB

     TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
     c_MC->SetFrameFillColor(0);
 
     RooPlot* xframe = delta_pT.frame(RooFit::Title(title_name));
     histo.plotOn(xframe);
 
     Int_t color = kRed+2;
     Double_t size_text = 0.020;
     DSCB.fitTo(histo, Range(-0.15,0.15));
     DSCB.plotOn(xframe, RooFit::LineColor(color),Name("#frac{#delta p_T}{p_T}"));
     DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
     c_MC->cd();
     int pt_number = pt_bin_number+1;
     int eta_number = eta_bin_number+1;
     TString title_s;
     if(type=="positive") {
     title_s =TString::Format("delta_pT_fit_positive_mu_pt_%d_eta_%d",pt_number, eta_number); 
     }
     if(type=="negative") {
     title_s =TString::Format("delta_pT_fit_negative_mu_pt_%d_eta_%d",pt_number, eta_number); 
     }
     c_MC->SetName(title_s);
     c_MC->SetTopMargin(1.5);
     xframe->getAttText()->SetTextSize(size_text);
     xframe->getAttText()->SetTextColor(color);
//     xframe->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
//     xframe->GetXaxis()->CenterTitle(true);
     xframe->GetYaxis()->SetTitle("");
     xframe->GetXaxis()->SetTitle("");
//     xframe->GetXaxis()->SetTitleOffset(1.4);
     xframe->chiSquare();

     float chi_square_value = xframe->chiSquare();
     //xframe->SetStats();
     xframe->Draw();
    
     
     //const TH1* histogram = histo.createHistogram("delta_pT");

    TLegend* leg2 = new TLegend(0.68, 0.68, 0.89, 0.89);
      leg2->SetFillColor(kWhite);
      leg2->SetLineColor(kBlack);
      leg2->AddEntry("histo","#frac{#Delta pT}{pT}", "EP");
      leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
      leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
//      leg2->AddEntry("histogram->GetBinContent(0)",Form("Underflow= %.0f",histogram->GetBinContent(0)),"");
//      leg2->AddEntry("histogram->GetBinContent(histogram->GetNbinsX() + 1)",Form("Overflow= %.0f",histogram->GetBinContent(histogram->GetNbinsX() + 1)),"");
      leg2->Draw("same");
 
//    std::cout<<"Underflow= "<<histogram->GetBinContent(0)<<std::endl;
//    std::cout<<"Overflow= "<<histogram->GetBinContent(histogram->GetNbinsX()+1)<<std::endl;
     gStyle->SetOptStat();
     c_MC->SaveAs((saving_path+ saving_name + ".pdf"));// + ".pdf");
     c_MC->Write(); 
     c_MC->Close();
     std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;

     if(type=="positive"){ 
     mean_positive_delta_pT[pt_bin_number][eta_bin_number] = Mean.getVal();
     mean_positive_error_delta_pT[pt_bin_number][eta_bin_number] = Mean.getError();

     sigma_positive_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getVal();
     sigma_positive_error_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getError();

     chi2_ndf_positive_delta_pT[pt_bin_number][eta_bin_number] =  chi_square_value ;
     }
    if(type=="negative"){ 
     mean_negative_delta_pT[pt_bin_number][eta_bin_number] = Mean.getVal();
     mean_negative_error_delta_pT[pt_bin_number][eta_bin_number] = Mean.getError();

     sigma_negative_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getVal();
     sigma_negative_error_delta_pT[pt_bin_number][eta_bin_number] = Sigma.getError();

     chi2_ndf_negative_delta_pT[pt_bin_number][eta_bin_number] = chi_square_value;
     }

 }

// DSCB eta bins 
    void class_reading :: fitting_DSCB_Zmass_alone(TH1F * hist_fit, TString saving_name,  TString title_name){
     TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
     TF1* f1 = new TF1("f1", "gaus", 87,96);
     hist_fit_clone->Fit(f1,"R+");
     float mean_gaus =91.19; float sigma_gaus = 2.49;
     mean_gaus = f1->GetParameter(1);
     sigma_gaus = f1->GetParameter(2);

       RooRealVar Z_mass("Z_mass", "Z mass", 65,115, "");
     RooDataHist histo("histo","dataset with var",Z_mass,hist_fit);
     RooRealVar Mean("Mean", "Mean",mean_gaus, 88, 95);
     RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.1,4);//sigma[decay]);
     RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.1, 50);//alphaL[decay]);
     RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 200);//expL[decay]);
     RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.1, 50);//alphaR[decay]);
     RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 200);//expR[decay]);
     RooMyPDF_DSCB DSCB("DSCB", "DSCB", Z_mass, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    // Fit the mass into DSCB
 
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);

      RooPlot* xframe = Z_mass.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);

    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(70,110));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("m(#mu^{+} #mu^{-})"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    TString title_s =title_name;
    c_MC->SetName(title_s);
    c_MC->SetTopMargin(1.5);
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
    xframe->GetYaxis()->SetTitle("");
    xframe->GetXaxis()->SetTitle("GeV");
//     xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    //xframe->SetStats();
    xframe->Draw();

      double chi_square_value = xframe->chiSquare();
    //const TH1* histogram = histo.createHistogram("delta_pT");

     TLegend* leg2 = new TLegend(0.68, 0.68, 0.89, 0.89);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Z_mass", "EP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");
     gStyle->SetOptStat();
     c_MC->SaveAs((saving_path+ saving_name + ".pdf"));// + ".pdf");
     c_MC->Write();
     c_MC->Close();
 } // end of fitting DSCB Z mass function


// DSCB eta bins 
    void class_reading :: fitting_DSCB_Zmass_with_eta(TH1F * hist_fit, TString saving_name, TString title_name, int eta_bin_number, TString thecandidate, TString type){
     TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
     TF1* f1 = new TF1("f1", "gaus", 87,96);
     hist_fit_clone->Fit(f1,"R+");
     float mean_gaus =91.19; float sigma_gaus = 2.49;
     mean_gaus = f1->GetParameter(1);
     sigma_gaus = f1->GetParameter(2);

     RooRealVar Z_mass("Z_mass", "Z mass", 65,115, "");
     RooDataHist histo("histo","dataset with var",Z_mass,hist_fit);
     RooRealVar Mean("Mean", "Mean",mean_gaus, 88, 95);
     RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.1,4);//sigma[decay]);
     RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.1, 50);//alphaL[decay]);
     RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 200);//expL[decay]);
     RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.1, 50);//alphaR[decay]);
     RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 200);//expR[decay]);
     RooMyPDF_DSCB DSCB("DSCB", "DSCB", Z_mass, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    // Fit the mass into DSCB
 
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);

      RooPlot* xframe = Z_mass.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);

    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(70,110));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("m(#mu^{+} #mu^{-})"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    int eta_number = eta_bin_number+1;
    TString title_s =TString::Format("Zmass_fit_positive_mu_eta_%d",eta_number);
    c_MC->SetName(title_s);
    c_MC->SetTopMargin(1.5);
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
//     xframe->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
//     xframe->GetXaxis()->CenterTitle(true);
    xframe->GetYaxis()->SetTitle("");
    xframe->GetXaxis()->SetTitle("GeV");
//     xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    //xframe->SetStats();
    xframe->Draw();

      double chi_square_value = xframe->chiSquare();
    //const TH1* histogram = histo.createHistogram("delta_pT");

     TLegend* leg2 = new TLegend(0.68, 0.68, 0.89, 0.89);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Z_mass", "EP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
       leg2->Draw("same");
            gStyle->SetOptStat();
            c_MC->SaveAs((saving_path+ saving_name + ".pdf"));// + ".pdf");
            c_MC->Write();
            c_MC->Close();

   
    if(thecandidate=="Zmass_gen_with_eta"){
      if(type=="positive") {
       mean_positive_gen_Zmass_with_eta[eta_bin_number] = Mean.getVal();
       mean_positive_error_gen_Zmass_with_eta[eta_bin_number] = Mean.getError();
       sigma_positive_gen_Zmass_with_eta[eta_bin_number] = Sigma.getVal();
       sigma_positive_error_gen_Zmass_with_eta[eta_bin_number] = Sigma.getError();
//chi2_ndf_positive_gen_Zmass_with_eta[eta_bin_number] = chi_square_value;
      }
      if(type=="negative") {
       mean_negative_gen_Zmass_with_eta[eta_bin_number]= Mean.getVal();
       mean_negative_error_gen_Zmass_with_eta[eta_bin_number]= Mean.getError();
       sigma_negative_gen_Zmass_with_eta[eta_bin_number]= Sigma.getVal();
       sigma_negative_error_gen_Zmass_with_eta[eta_bin_number]= Sigma.getError();
//      chi2_ndf_negative_gen_Zmass_with_eta[eta_bin_number]= chi_square_value ;
      }
    }


    if(thecandidate=="Zmass_reco_with_eta"){
      if(type=="positive") {
       mean_positive_Zmass_with_eta[eta_bin_number] = Mean.getVal();
       mean_positive_error_Zmass_with_eta[eta_bin_number] = Mean.getError();
       sigma_positive_Zmass_with_eta[eta_bin_number] = Sigma.getVal();
       sigma_positive_error_Zmass_with_eta[eta_bin_number] = Sigma.getError();
//      chi2_ndf_positive_Zmass_with_eta[eta_bin_number] =  chi_square_value;
      }
      if(type=="negative") {
       mean_negative_Zmass_with_eta[eta_bin_number] = Mean.getVal();
       mean_negative_error_Zmass_with_eta[eta_bin_number] = Mean.getError();
       sigma_negative_Zmass_with_eta[eta_bin_number] = Sigma.getVal();
       sigma_negative_error_Zmass_with_eta[eta_bin_number] = Sigma.getError();
//      chi2_ndf_negative_Zmass_with_eta[eta_bin_number] = chi_square_value ;
      }
    }

    if(thecandidate=="Zmass_FSR_reco_with_eta"){
      if(type=="positive"){
       mean_positive_Zmass_FSR_with_eta[eta_bin_number] = Mean.getVal();
       mean_positive_error_Zmass_FSR_with_eta[eta_bin_number] = Mean.getError();
       sigma_positive_Zmass_FSR_with_eta[eta_bin_number] = Sigma.getVal();
       sigma_positive_error_Zmass_FSR_with_eta[eta_bin_number] = Sigma.getError();
//      chi2_ndf_positive_Zmass_FSR_with_eta[eta_bin_number] = chi_square_value;
      }
      if(type=="negative"){
       mean_negative_Zmass_FSR_with_eta[eta_bin_number] = Mean.getVal();
       mean_negative_error_Zmass_FSR_with_eta[eta_bin_number] = Mean.getError();
       sigma_negative_Zmass_FSR_with_eta[eta_bin_number] = Sigma.getVal();
       sigma_negative_error_Zmass_FSR_with_eta[eta_bin_number] = Sigma.getError();
//      chi2_ndf_negative_Zmass_FSR_with_eta[eta_bin_number] =  chi_square_value;
      }
    } // end of else if lop
 } // end of fitting DSCB Z mass function

// DSCB pt bins 
    void class_reading :: fitting_DSCB_Zmass_with_pt(TH1F * hist_fit, TString saving_name, TString title_name, int pt_bin_number, TString thecandidate, TString type){
     TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
     TF1* f1 = new TF1("f1", "gaus", 87,96);
     hist_fit_clone->Fit(f1,"R+");
     float mean_gaus =91.19; float sigma_gaus = 2.49;
     mean_gaus = f1->GetParameter(1);
     sigma_gaus = f1->GetParameter(2);

       RooRealVar Z_mass("Z_mass", "Z mass", 65,115, "");
     RooDataHist histo("histo","dataset with var",Z_mass,hist_fit);
     RooRealVar Mean("Mean", "Mean",mean_gaus, 88, 95);
     RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.1,4);//sigma[decay]);
     RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.1, 50);//alphaL[decay]);
     RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 200);//expL[decay]);
     RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.1, 50);//alphaR[decay]);
     RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 200);//expR[decay]);
     RooMyPDF_DSCB DSCB("DSCB", "DSCB", Z_mass, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    // Fit the mass into DSCB
 
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);

      RooPlot* xframe = Z_mass.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);

    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(70,110));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("m(#mu^{+} #mu^{-})"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    int pt_number = pt_bin_number+1;
    TString title_s =TString::Format("Zmass_fit_positive_mu_pt_%d",pt_number);
    c_MC->SetName(title_s);
    c_MC->SetTopMargin(1.5);
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
//     xframe->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
//     xframe->GetXaxis()->CenterTitle(true);
    xframe->GetYaxis()->SetTitle("");
    xframe->GetXaxis()->SetTitle("GeV");
//     xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    //xframe->SetStats();
    xframe->Draw();

      double chi_square_value = xframe->chiSquare();
    //const TH1* histogram = histo.createHistogram("delta_pT");

     TLegend* leg2 = new TLegend(0.68, 0.68, 0.89, 0.89);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Z_mass", "EP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
       leg2->Draw("same");
            gStyle->SetOptStat();
            c_MC->SaveAs((saving_path+ saving_name + ".pdf"));// + ".pdf");
            c_MC->Write();
            c_MC->Close();

   
    if(thecandidate=="Zmass_gen_with_pt"){
      if(type=="positive") {
       mean_positive_gen_Zmass_with_pt[pt_bin_number] = Mean.getVal();
       mean_positive_error_gen_Zmass_with_pt[pt_bin_number] = Mean.getError();
       sigma_positive_gen_Zmass_with_pt[pt_bin_number] = Sigma.getVal();
       sigma_positive_error_gen_Zmass_with_pt[pt_bin_number] = Sigma.getError();
//      chi2_ndf_positive_gen_Zmass_with_pt[pt_bin_number] = chi_square_value;
      }
      if(type=="negative") {
       mean_negative_gen_Zmass_with_pt[pt_bin_number]= Mean.getVal();
       mean_negative_error_gen_Zmass_with_pt[pt_bin_number]= Mean.getError();
       sigma_negative_gen_Zmass_with_pt[pt_bin_number]= Sigma.getVal();
       sigma_negative_error_gen_Zmass_with_pt[pt_bin_number]= Sigma.getError();
//      chi2_ndf_negative_gen_Zmass_with_pt[pt_bin_number]= chi_square_value ;
      }
    }


    if(thecandidate=="Zmass_reco_with_pt"){
      if(type=="positive") {
       mean_positive_Zmass_with_pt[pt_bin_number] = Mean.getVal();
       mean_positive_error_Zmass_with_pt[pt_bin_number] = Mean.getError();
       sigma_positive_Zmass_with_pt[pt_bin_number] = Sigma.getVal();
       sigma_positive_error_Zmass_with_pt[pt_bin_number] = Sigma.getError();
//      chi2_ndf_positive_Zmass_with_pt[pt_bin_number] =  chi_square_value;
      }
      if(type=="negative") {
       mean_negative_Zmass_with_pt[pt_bin_number] = Mean.getVal();
       mean_negative_error_Zmass_with_pt[pt_bin_number] = Mean.getError();
       sigma_negative_Zmass_with_pt[pt_bin_number] = Sigma.getVal();
       sigma_negative_error_Zmass_with_pt[pt_bin_number] = Sigma.getError();
//      chi2_ndf_negative_Zmass_with_pt[pt_bin_number] = chi_square_value ;
      }
    }

    if(thecandidate=="Zmass_FSR_reco_with_pt"){
      if(type=="positive"){
       mean_positive_Zmass_FSR_with_pt[pt_bin_number] = Mean.getVal();
       mean_positive_error_Zmass_FSR_with_pt[pt_bin_number] = Mean.getError();
       sigma_positive_Zmass_FSR_with_pt[pt_bin_number] = Sigma.getVal();
       sigma_positive_error_Zmass_FSR_with_pt[pt_bin_number] = Sigma.getError();
//      chi2_ndf_positive_Zmass_FSR_with_pt[pt_bin_number] = chi_square_value;
      }
      if(type=="negative"){
       mean_negative_Zmass_FSR_with_pt[pt_bin_number] = Mean.getVal();
       mean_negative_error_Zmass_FSR_with_pt[pt_bin_number] = Mean.getError();
       sigma_negative_Zmass_FSR_with_pt[pt_bin_number] = Sigma.getVal();
       sigma_negative_error_Zmass_FSR_with_pt[pt_bin_number] = Sigma.getError();
//      chi2_ndf_negative_Zmass_FSR_with_pt[pt_bin_number] =  chi_square_value;
      }
    } // end of else if lop
 } // end of fitting DSCB Z mass function

// just a seaprate DSCB function to fit Z mass also, I know its bad coding, but I have to submit the plots tomorrow 
   // implementing DSCB fit for Z mass
    void class_reading :: fitting_DSCB_Zmass(TH1F * hist_fit, TString saving_name, TString title_name, int pt_bin_number, int eta_bin_number, TString thecandidate, TString type){
     TH1F * hist_fit_clone = (TH1F*) hist_fit->Clone();
     TF1* f1 = new TF1("f1", "gaus", 87,96);
     hist_fit_clone->Fit(f1,"R+");
     float mean_gaus =91.19; float sigma_gaus = 2.49;
     mean_gaus = f1->GetParameter(1);
     sigma_gaus = f1->GetParameter(2);

       RooRealVar Z_mass("Z_mass", "Z mass", 65,115, "");
     RooDataHist histo("histo","dataset with var",Z_mass,hist_fit);
     RooRealVar Mean("Mean", "Mean",mean_gaus, 88, 95);
     RooRealVar Sigma("#sigma", "#sigma", sigma_gaus, 0.1,4);//sigma[decay]);
     RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.1, 50);//alphaL[decay]);
     RooRealVar ExpL("n_{L}", "n_{L}", 5, 0.1, 200);//expL[decay]);
     RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.3, 0.1, 50);//alphaR[decay]);
     RooRealVar ExpR("n_{R}", "n_{R}", 12, 0.1, 200);//expR[decay]);
     RooMyPDF_DSCB DSCB("DSCB", "DSCB", Z_mass, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
      output_file<<"fitting positive muon : mean "<<mean_gaus<<" sigma : "<<sigma_gaus<<std::endl;
    // Fit the mass into DSCB
 
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);

      RooPlot* xframe = Z_mass.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);

    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(70,110));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("m(#mu^{+} #mu^{-})"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    int pt_number = pt_bin_number+1;
    int eta_number = eta_bin_number+1;
    TString title_s =TString::Format("Zmass_fit_positive_mu_pt_%d_eta_%d",pt_number, eta_number);
    c_MC->SetName(title_s);
    c_MC->SetTopMargin(1.5);
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
//     xframe->GetXaxis()->SetTitle("#frac{#Delta pT}{pT}");
//     xframe->GetXaxis()->CenterTitle(true);
    xframe->GetYaxis()->SetTitle("");
    xframe->GetXaxis()->SetTitle("GeV");
//     xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    //xframe->SetStats();
    xframe->Draw();

      double chi_square_value = xframe->chiSquare();
    //const TH1* histogram = histo.createHistogram("delta_pT");

     TLegend* leg2 = new TLegend(0.68, 0.68, 0.89, 0.89);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Z_mass", "EP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
       leg2->Draw("same");
            gStyle->SetOptStat();
            c_MC->SaveAs((saving_path+ saving_name + ".pdf"));// + ".pdf");
            c_MC->Write();
            c_MC->Close();
            std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;
             new_file<<" candidate "<<thecandidate<<" \t type "<<type<<std::endl;
       
              new_file<<"Zmass mean "<<Mean.getVal()<<std::endl;

   
    if(thecandidate=="Zmass_gen"){
      if(type=="positive") {
       new_file<<"Zmass inside test "<<Mean.getVal()<<std::endl;
       mean_positive_gen_Zmass[pt_bin_number][eta_bin_number] = Mean.getVal();
       new_file<<" printing Zmass inside test "<<mean_positive_Zmass[pt_bin_number][eta_bin_number]<<std::endl;
       mean_positive_error_gen_Zmass[pt_bin_number][eta_bin_number] = Mean.getError();
       sigma_positive_gen_Zmass[pt_bin_number][eta_bin_number] = Sigma.getVal();
       sigma_positive_error_gen_Zmass[pt_bin_number][eta_bin_number] = Sigma.getError();
       chi2_ndf_positive_gen_Zmass[pt_bin_number][eta_bin_number] = chi_square_value;
      }
      if(type=="negative") {
       mean_negative_gen_Zmass[pt_bin_number][eta_bin_number] = Mean.getVal();
       mean_negative_error_gen_Zmass[pt_bin_number][eta_bin_number] = Mean.getError();
       sigma_negative_gen_Zmass[pt_bin_number][eta_bin_number] = Sigma.getVal();
       sigma_negative_error_gen_Zmass[pt_bin_number][eta_bin_number] = Sigma.getError();
       chi2_ndf_negative_gen_Zmass[pt_bin_number][eta_bin_number] = chi_square_value ;
      }
    }


    if(thecandidate=="Zmass"){
      if(type=="positive") {
       new_file<<"Zmass inside test "<<Mean.getVal()<<std::endl;
       mean_positive_Zmass[pt_bin_number][eta_bin_number] = Mean.getVal();
       new_file<<" printing Zmass inside test "<<mean_positive_Zmass[pt_bin_number][eta_bin_number]<<std::endl;
       mean_positive_error_Zmass[pt_bin_number][eta_bin_number] = Mean.getError();
       sigma_positive_Zmass[pt_bin_number][eta_bin_number] = Sigma.getVal();
       sigma_positive_error_Zmass[pt_bin_number][eta_bin_number] = Sigma.getError();
       chi2_ndf_positive_Zmass[pt_bin_number][eta_bin_number] =  chi_square_value;
      }
      if(type=="negative") {
       mean_negative_Zmass[pt_bin_number][eta_bin_number] = Mean.getVal();
       mean_negative_error_Zmass[pt_bin_number][eta_bin_number] = Mean.getError();
       sigma_negative_Zmass[pt_bin_number][eta_bin_number] = Sigma.getVal();
       sigma_negative_error_Zmass[pt_bin_number][eta_bin_number] = Sigma.getError();
       chi2_ndf_negative_Zmass[pt_bin_number][eta_bin_number] = chi_square_value ;
      }
    }

    if(thecandidate=="Zmass_FSR"){
      if(type=="positive"){
       mean_positive_Zmass_FSR[pt_bin_number][eta_bin_number] = Mean.getVal();
       mean_positive_error_Zmass_FSR[pt_bin_number][eta_bin_number] = Mean.getError();
       sigma_positive_Zmass_FSR[pt_bin_number][eta_bin_number] = Sigma.getVal();
       sigma_positive_error_Zmass_FSR[pt_bin_number][eta_bin_number] = Sigma.getError();
       chi2_ndf_positive_Zmass_FSR[pt_bin_number][eta_bin_number] = chi_square_value;
      }
      if(type=="negative"){
       mean_negative_Zmass_FSR[pt_bin_number][eta_bin_number] = Mean.getVal();
       mean_negative_error_Zmass_FSR[pt_bin_number][eta_bin_number] = Mean.getError();
       sigma_negative_Zmass_FSR[pt_bin_number][eta_bin_number] = Sigma.getVal();
       sigma_negative_error_Zmass_FSR[pt_bin_number][eta_bin_number] = Sigma.getError();
       chi2_ndf_negative_Zmass_FSR[pt_bin_number][eta_bin_number] =  chi_square_value;
      }
    } // end of else if lop
 } // end of fitting DSCB Z mass function

/*void class_reading :: diff_gen_reco(TString type){
      if(type=="Zmass_FSR"){
       dir_name = "mean_resolution_Zmass_diff";
       graph_name_positive = "mean_Zmass_positive_";
       graph_name_negative = "mean_Zmass_negative_";
       graph_name_positive_sigma = "sigam_Zmass_positive_";
       graph_name_negative_sigma = "sigam_Zmass_negative_";

       saving_name_negative = "mean_sigma_Zmass_plots/graph_resolution_Zmass_mean_negative_eta_";
       saving_name_positive = "mean_sigma_Zmass_plots/graph_resolution_Zmass_mean_positive_eta_";
       saving_name_positive_sigma =  "mean_sigma_plots/graph_resolution_Zmass_sigma_positive_eta_";
       saving_name_negative_sigma =  "mean_sigma_plots/graph_resolution_Zmass_sigma_negative_eta_";

       text_file_name = "mean_sigma_Zmass.txt";
 //      string_title_mulitgraph_positive = "Mean Z mass mu^{+}";
 //      string_title_mulitgraph_negative = "Mean Z mass mu^{-}";

       range_low = 90.5;
       range_high = 91.3;
       range_high_sigma = 3.5;

      title_multigraph_positive = "Mean Z mass : #mu^{+}";
      title_multigraph_negative = "Mean Z mass : #mu^{-}";
    
      saving_name_multigraph_positive_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_positive_muon_positive_eta.pdf";
      saving_name_multigraph_positive_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_positive_muon_negative_eta.pdf";
      saving_name_multigraph_negative_muon_positive_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_negative_muon_positive_eta.pdf";
      saving_name_multigraph_negative_muon_negative_eta = "mean_sigma_Zmass_comparison_plots/mean_Zmass_negative_muon_negative_eta.pdf";

     }


  
} */
 void new_program_write(){ 
 
   // this code is to analyze all the histograms produced for the Jpsi sample : histograms ar     e Z mass histogram in all the pT categories, for both mu^+ and mu^-   
   // derived class reading is the class where I am performing operation to fit the histogram     s , and dervied class reading is using functions from base class reading to perform this oper     ations   
   gROOT->SetBatch(kTRUE);    
   std::cout<<" initiated a class obj"<<std::endl;   
   /*derived_class_reading obj; 
   // initializing is to open the root file and to read all the histograms 
   obj.initializing();  
   std::cout<<" class obj declared successfully"<<std::endl;   
   obj.saving_histograms(); 
   std::cout<<" saving histogram delta pT"<<std::endl;   
   obj.fitting_delta_pT(); 
  */
   // plotting every histogram on canvas with the plot_histograms function   
   //   obj.plot_histograms();  
   // plotting inclusive and binned distribution on same canvas  
   //obj.plotting_inclusive_bin(); 
   // to do BW and DSCB fits and find mean and sigma values  and this mean and sigma values a     re further stored and used to make final plots in the functio nmean_sigma_calculation 
   //   obj.fitting_histograms(); 
   // saving text file also has the calculation for mean and sigma which will be used further     e 
   //   obj.saving_text_file(); 
   // graph to plot mean z mass with pT 
  //   obj.mean_sigma_calculation(); 
   // to save the pT distributions for individual histograms 
   // calling a derived class functions available with these functions 
 

   // writing in terms of objection of normal class
   std::cout<<"starting the program"<<std::endl;
   class_reading obj_new; 
   obj_new.initializing();
   obj_new.saving_histogram_pT();
   obj_new.saving_histogram_delta_pT();
   obj_new.fitting_Zmass_distribution();
   obj_new.fitting_Zmass_FSR_distribution();
   // obj_new.fitting_delta_pT_distribution();

   // obj_new.plotting_mean_delta_pT("delta_pT");
   obj_new.plotting_mean_delta_pT("Zmass");
   obj_new.plotting_mean_with_pt("Zmass_gen_with_pt");
   obj_new.plotting_mean_with_pt("Zmass_reco_with_pt");
   obj_new.plotting_mean_with_pt("Zmass_reco_FSR_with_pt");
   obj_new.plotting_mean_with_eta("Zmass_gen_with_eta");
   obj_new.plotting_mean_with_eta("Zmass_reco_with_eta");
   obj_new.plotting_mean_with_eta("Zmass_reco_FSR_with_eta");

//   obj_new.plotting_mean_delta_pT("Zmass_gen");
//   obj_new.plotting_mean_delta_pT("Zmass_FSR");
   //obj_new.temperature_plot_delta_pT();

//   obj_new.diff_gen_reco("Zmass_FSR");
  } 




