#define DY_2018_macro_cxx 
#include "DY_2018_macro.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//#include "RooMyPDF_DSCB.h"
//#include "RooMyPDF_DSCB.cxx"

void DY_2018_macro::Loop()
{
     gROOT->SetBatch(kTRUE);
//     double pi = TMath::Pi();
//   In a ROOT session, you can do:
//      root> .L DY_2018_macro.C
//      root> DY_2018_macro t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//   for (Long64_t jentry=0; jentry<10000;jentry++){
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) { 
      std::cout<<"entry is breaking"<<std::endl;
      std::cout<<"entry number "<<jentry<<std::endl; 
      break;
     }
     nb = fChain->GetEntry(jentry);   
     nbytes += nb;
     if(jentry%100000==0) std::cout<<" entry "<<jentry<<std::endl;

//  if(pT1_genFromReco <2) continue;
//  if(pT1_genFromReco <2) continue;
     if(GENmass2l <70 || GENmass2l >110) continue;
       // filling all the inclusive gen and reco distributions : genLep - gen, pT1 - reco
       // filling pT, id and mass
       //
       counter++; 
       hist_gen_pT1 ->Fill(pT1_genFromReco);      
       hist_gen_pT2 ->Fill(pT2_genFromReco);      
       hist_gen_eta1 ->Fill(genLep_eta1);      
       hist_gen_eta2 ->Fill(genLep_eta2);      
       hist_gen_phi1 ->Fill(genLep_phi1);      
       hist_gen_phi2 ->Fill(genLep_phi2);      

       hist_reco_eta1 ->Fill(eta1);      
       hist_reco_eta2 ->Fill(eta2);      
       hist_reco_pT1 ->Fill(pT_FSR1);      
       hist_reco_pT2 ->Fill(pT_FSR2);      
       hist_reco_phi1 ->Fill(phi1);      
       hist_reco_phi2 ->Fill(phi2);      

       //  hist_gen_pT1 ->Fill(pT1_genFromReco);      
       //  hist_gen_pT2 ->Fill(pT2_genFromReco);      
       //  hist_reco_pT1 ->Fill(pT1);      
       //  hist_reco_pT2 ->Fill(pT2);      
       hist_id1 ->Fill(Id1);      
       hist_id2 ->Fill(Id2);      
       hist_gen_Zmass ->Fill(GENmass2l);      
       hist_reco_Zmass ->Fill(massZ);      
       hist_reco_Zmass_FSR->Fill(massZ_FSR);      

       // when positive muon has higher pT than negative one
      if(Id1==-13)    
      {
         // first muon is positive and second is negative
         // dividing in pT bins
         counter_first++; 
         // pT and Id of positive and negative muon in this event
         hist_gen_positive_mu->Fill(pT1_genFromReco);
         hist_id_positive_mu->Fill(Id1);        
         hist_gen_negative_mu->Fill(pT2_genFromReco);
         hist_id_negative_mu->Fill(Id2); 

         hist_gen_positive_eta->Fill(genLep_eta1);
         hist_gen_positive_phi->Fill(genLep_phi1);

         hist_gen_negative_eta->Fill(genLep_eta2);
         hist_gen_negative_phi->Fill(genLep_phi2);
         // reco distribution of gen leptons
         // smear pT of individual lepton and then reconstruct Z mass

         // analyze pT meaning we are putting events into bins of pT and eta
         
//          std::cout<<" id : entering checking loop "<<Id1<<" :  "<<pT1_genFromReco<<" second "<<pT2_genFromReco<<" pT1 "<<pT1<<std::endl;
          for(int i=0; i<pT_bins; i++){
            for(int j=0; j<eta_bins; j++){
             // falling in ith bin for positive muon
             if( (pt_list[i] <= pT1_genFromReco && pT1_genFromReco < pt_list[i+1]) && (eta_list[j]<= genLep_eta1 && genLep_eta1 <eta_list[j+1]) ){
 //            std::cout<<"positive first :  bin number "<<i<<" j "<<j<<std::endl;
              //std::cout<<"positive first :  bin number "<<i<<" j "<<j<<std::endl;
                // falling in eta bins for positive muon
                pT_positive_bin[i][j]->Fill(pT1_genFromReco);
                // filling reco pT depending on the genpT bin: 
                pT_positive_bin_reco[i][j]->Fill(pT_FSR1);
                // first since first muon is positive
                delta_pT_positive_bin_reco_first(i,j);
                hist_positive_gen_Zmass[i][j]->Fill(GENmass2l);
                hist_positive_reco_Zmass[i][j]->Fill(massZ);
                hist_positive_reco_Zmass_FSR[i][j]->Fill(massZ_FSR);
            }   // ended the binning thing
           } //loop over eta
         } // loop over pT

         // in eta bins 
          for(int i=0; i<eta_bins; i++){
             if(  eta_list[i]<= genLep_eta1 && genLep_eta1 <eta_list[i+1] ){
                hist_positive_gen_Zmass_with_eta[i]->Fill(GENmass2l);
                hist_positive_reco_Zmass_with_eta[i]->Fill(massZ);
                hist_positive_reco_Zmass_FSR_with_eta[i]->Fill(massZ_FSR);
             }
          }
         // in pt bins 
          for(int i=0; i<pT_bins; i++){
             if( (pt_list[i] <= pT1_genFromReco && pT1_genFromReco < pt_list[i+1])){
                hist_positive_gen_Zmass_with_pt[i]->Fill(GENmass2l);
                hist_positive_reco_Zmass_with_pt[i]->Fill(massZ);
                hist_positive_reco_Zmass_FSR_with_pt[i]->Fill(massZ_FSR);
             }
          }


    
       // Loop for second muon ie negative muon 
          for(int i=0; i<pT_bins; i++){
            for(int j=0; j<eta_bins; j++){
             // falling in ith bin for positive muon
             if( (pt_list[i] <= pT2_genFromReco && pT2_genFromReco < pt_list[i+1]) && (eta_list[j]<= genLep_eta2 && genLep_eta2 <eta_list[j+1]) ){
              //std::cout<<"negative first :  bin number "<<i<<" j "<<j<<std::endl;
                // falling in eta bins for positive muon
                if(debug_pT == true) pT_negative_bin[i][j]->Fill(pT2_genFromReco);
                // filling reco pT depending on the genpT bin: 
                if(debug_pT == true) pT_negative_bin_reco[i][j] ->Fill(pT_FSR2);
                // first since first muon is positive
                delta_pT_negative_bin_reco_first(i,j);
                hist_negative_gen_Zmass[i][j]->Fill(GENmass2l);
                hist_negative_reco_Zmass[i][j]->Fill(massZ);
                hist_negative_reco_Zmass_FSR[i][j]->Fill(massZ_FSR);

            }   // ended the binning thing
           } //loop over eta
         } // loop over pT

          // in eta bins
          for(int i=0; i<eta_bins; i++){
             if(  eta_list[i]<= genLep_eta2 && genLep_eta2 <eta_list[i+1] ){
                hist_negative_gen_Zmass_with_eta[i]->Fill(GENmass2l);
                hist_negative_reco_Zmass_with_eta[i]->Fill(massZ);
                hist_negative_reco_Zmass_FSR_with_eta[i]->Fill(massZ_FSR);
             }
          }
         // in pT bins 
          for(int i=0; i<pT_bins; i++){
             if( (pt_list[i] <= pT2_genFromReco && pT2_genFromReco < pt_list[i+1])){
                hist_negative_gen_Zmass_with_pt[i]->Fill(GENmass2l);
                hist_negative_reco_Zmass_with_pt[i]->Fill(massZ);
                hist_negative_reco_Zmass_FSR_with_pt[i]->Fill(massZ_FSR);
             }
          }

          if(debug_count == true){ 
             if(pT2_genFromReco >=200 || pT2_genFromReco <=0){
             count_negative_mu++;
             }
            } 
     } // end of Id checking for first muon as positive one

       // Storing positive muons with pT less than negative muons
       if(Id2==-13)    
       {
         counter_second++; 
         hist_gen_positive_mu->Fill(pT2_genFromReco);
         hist_id_positive_mu->Fill(Id2);        
         hist_gen_negative_mu->Fill(pT1_genFromReco);
         hist_id_negative_mu->Fill(Id1); 
         
         hist_gen_positive_eta->Fill(genLep_eta2);
         hist_gen_positive_phi->Fill(genLep_phi2);
         hist_gen_negative_eta->Fill(genLep_eta1);
         hist_gen_negative_phi->Fill(genLep_phi1);
//          std::cout<<" id : entering checking loop "<<Id1<<" :  "<<pT1_genFromReco<<" second "<<pT2_genFromReco<<" pT1 "<<pT1<<std::endl;

         for(int i=0; i<pT_bins; i++){
          for(int j=0; j<eta_bins; j++){
           if( (pt_list[i] <= pT2_genFromReco && pT2_genFromReco < pt_list[i+1]) && (eta_list[j]<= genLep_eta2 && genLep_eta2 <eta_list[j+1]) ){
//             std::cout<<"positive second :  bin number "<<i<<" j "<<j<<std::endl;
           pT_positive_bin[i][j] ->Fill(pT2_genFromReco);
           pT_positive_bin_reco[i][j] ->Fill(pT_FSR2);
           delta_pT_positive_bin_reco_second(i,j);
           hist_positive_gen_Zmass[i][j]->Fill(GENmass2l);
           hist_positive_reco_Zmass[i][j]->Fill(massZ);
           hist_positive_reco_Zmass_FSR[i][j]->Fill(massZ_FSR);

          }
         }
        }

         // in eta bins 
          for(int i=0; i<eta_bins; i++){
             if(  eta_list[i]<= genLep_eta2 && genLep_eta2 <eta_list[i+1] ){
                hist_positive_gen_Zmass_with_eta[i]->Fill(GENmass2l);
                hist_positive_reco_Zmass_with_eta[i]->Fill(massZ);
                hist_positive_reco_Zmass_FSR_with_eta[i]->Fill(massZ_FSR);
             }
          }
         // in pt bins 
          for(int i=0; i<pT_bins; i++){
             if( (pt_list[i] <= pT2_genFromReco && pT2_genFromReco < pt_list[i+1])){
                hist_positive_gen_Zmass_with_pt[i]->Fill(GENmass2l);
                hist_positive_reco_Zmass_with_pt[i]->Fill(massZ);
                hist_positive_reco_Zmass_FSR_with_pt[i]->Fill(massZ_FSR);
             }
          }



         for(int i=0; i<pT_bins; i++){
           for(int j=0; j<eta_bins; j++){
           // falling in ith bin for positive muon
           if( (pt_list[i] <= pT1_genFromReco && pT1_genFromReco < pt_list[i+1]) && (eta_list[j]<= genLep_eta1 && genLep_eta1 <eta_list[j+1]) ){
              //std::cout<<"negative second :  bin number "<<i<<" j "<<j<<std::endl;
              // falling in eta bins for positive muon
              if(debug_pT == true) pT_negative_bin[i][j]->Fill(pT1_genFromReco);
              // filling reco pT depending on the genpT bin: 
              if(debug_pT == true) pT_negative_bin_reco[i][j] ->Fill(pT_FSR1);
              // first since first muon is positive
              delta_pT_negative_bin_reco_second(i,j);
              hist_negative_gen_Zmass[i][j]->Fill(GENmass2l);
              hist_negative_reco_Zmass[i][j]->Fill(massZ);
              hist_negative_reco_Zmass_FSR[i][j]->Fill(massZ_FSR);

            }   // ended the binning thing
           } //loop over eta
         } // loop over pT

         for(int i=0; i<eta_bins; i++){
             if(  eta_list[i]<= genLep_eta1 && genLep_eta1 <eta_list[i+1] ){
                hist_negative_gen_Zmass_with_eta[i]->Fill(GENmass2l);
                hist_negative_reco_Zmass_with_eta[i]->Fill(massZ);
                hist_negative_reco_Zmass_FSR_with_eta[i]->Fill(massZ_FSR);
             }
         }
         // in pT bins 
          for(int i=0; i<pT_bins; i++){
             if( (pt_list[i] <= pT1_genFromReco && pT1_genFromReco < pt_list[i+1])){
                hist_negative_gen_Zmass_with_pt[i]->Fill(GENmass2l);
                hist_negative_reco_Zmass_with_pt[i]->Fill(massZ);
                hist_negative_reco_Zmass_FSR_with_pt[i]->Fill(massZ_FSR);
             }
          }


        if(debug_count == true){ 
           if(pT2_genFromReco >=200 ||pT2_genFromReco <=0){
           count_positive_mu++;
          }
        }
 
        if(debug_count == true){ 
         if(pT1_genFromReco >=200 || pT1_genFromReco <=0){
          count_negative_mu++;
          }
        } 
    } //end of second iD loop for first muon as negative muon

  }  // end of entries for statement
   if(debug_base == true){
   std::cout<<"number of positive muons ignored in analyzing"<<count_positive_mu<<std::endl;
   std::cout<<"number of negative muons ignored in analyzing"<<count_negative_mu<<std::endl;
   std::cout<<"total events "<<counter<<std::endl;
   std::cout<<"total events first poisitive "<<counter_first<<std::endl;
   std::cout<<"total events second poisitive "<<counter_second<<std::endl;
   } 
   
   // to write every histogram to an output root file 
   saving_histogram();
   total_entries();
   // end of program
  
}
