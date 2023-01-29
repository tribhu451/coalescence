#include "observables.h"

observables::observables(input_paramters &iparam_, read_input_file* rif_) : iparam(iparam_), rif(rif_) {
  deutron_mass = 1.87561 ; // in GeV
  proton_mass  = 0.93800 ; // in GeV
  neutron_mass = 0.93800 ; // in GeV
  deutron_sigma_rho = 2.473 ; // in fm
  triton_sigma_lambda = 1.759 ;  // in fm
  triton_sigma_rho = 1.759 ;  // in fm 
  He3_sigma_lambda =  1.390; // in fm .  Reference : Phys. Rev. C 95, 054913
  He3_sigma_rho = 1.390 ;    // in fm
  velocity_failure_counter = 0 ; 

  rand=new TRandom3();
  rand->SetSeed(0);
}


void observables::calculate_dnchdeta_eta( double pT_min, double pT_max ){
    std::cout << "Calculating dN/deta within pT range " 
              << pT_min << "-" << pT_max << " ... "
              << std::endl ; 


  const int _N_HISTOGRAMS_DNCHDETA_ETA = 9 ; 

  double Eta_min   = -8. ; 
  double Eta_max   =  8. ; 
  int    Eta_bins  =  48 ; 

  TH1D*                 H1D_DNCHDETA_ETA[_N_HISTOGRAMS_DNCHDETA_ETA] ; 
  H1D_DNCHDETA_ETA[0] = new TH1D("H0", "charged_particle_eta_differential_yield", Eta_bins, Eta_min, Eta_max);
  H1D_DNCHDETA_ETA[1] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H1");  H1D_DNCHDETA_ETA[1]->SetTitle("pion_plus");
  H1D_DNCHDETA_ETA[2] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H2");  H1D_DNCHDETA_ETA[2]->SetTitle("pion_minus");
  H1D_DNCHDETA_ETA[3] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H3");  H1D_DNCHDETA_ETA[3]->SetTitle("kaon_plus");
  H1D_DNCHDETA_ETA[4] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H4");  H1D_DNCHDETA_ETA[4]->SetTitle("kaon_minus");
  H1D_DNCHDETA_ETA[5] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H5");  H1D_DNCHDETA_ETA[5]->SetTitle("proton");
  H1D_DNCHDETA_ETA[6] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H6");  H1D_DNCHDETA_ETA[6]->SetTitle("anti_proton");
  H1D_DNCHDETA_ETA[7] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H7");  H1D_DNCHDETA_ETA[7]->SetTitle("lambda");
  H1D_DNCHDETA_ETA[8] = (TH1D*) H1D_DNCHDETA_ETA[0]->Clone("H8");  H1D_DNCHDETA_ETA[8]->SetTitle("anti_lambda");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      //double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Eta = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Pt > pT_max || Pt < pT_min)
        continue ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        H1D_DNCHDETA_ETA[0]->Fill(Eta,1.);
      }

      if(PID == 211){
        H1D_DNCHDETA_ETA[1]->Fill(Eta,1.);
      }
      if(PID == -211){
        H1D_DNCHDETA_ETA[2]->Fill(Eta,1.);
      }
      if(PID == 321){
        H1D_DNCHDETA_ETA[3]->Fill(Eta,1.);
      }
      if(PID == -321){
        H1D_DNCHDETA_ETA[4]->Fill(Eta,1.);
      }
      if(PID == 2212){
        H1D_DNCHDETA_ETA[5]->Fill(Eta,1.);
      }
      if(PID == -2212){
        H1D_DNCHDETA_ETA[6]->Fill(Eta,1.);
      }
      if(PID == 3122){
        H1D_DNCHDETA_ETA[7]->Fill(Eta,1.);
      }
      if(PID == -3122){
        H1D_DNCHDETA_ETA[8]->Fill(Eta,1.);
      }
     } // particle loop
    } // event loop


  double dX = ( Eta_max - Eta_min ) / Eta_bins ; 
  for(int ii=0; ii<_N_HISTOGRAMS_DNCHDETA_ETA; ii++)
        H1D_DNCHDETA_ETA[ii]->Scale(1.0 / (nEvents * dX));

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_DNCHDETA_ETA] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122" };
  int        hadron_index[_N_HISTOGRAMS_DNCHDETA_ETA] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8   };


 for(int ix =0; ix < _N_HISTOGRAMS_DNCHDETA_ETA; ix++){
   output_filename.str("");
   output_filename << "results/dnchdeta_eta-" << hadron_name[ix];
   output_filename << "_pT_" << pT_min << "_" << pT_max ;
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= H1D_DNCHDETA_ETA[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << H1D_DNCHDETA_ETA[hadron_index[ix]]->GetBinCenter(i) << "\t" << H1D_DNCHDETA_ETA[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << H1D_DNCHDETA_ETA[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_DNCHDETA_ETA ; ixx++){
     H1D_DNCHDETA_ETA[ixx]->Clear(); 
   }


delete H1D_DNCHDETA_ETA[_N_HISTOGRAMS_DNCHDETA_ETA] ; 

}



void observables::calculate_dndy_y( double pT_min, double pT_max ){
    std::cout << "Calculating dN/dy within pT range " 
              << pT_min << "-" << pT_max << " ... "
              << std::endl ; 

  const int _N_HISTOGRAMS_DNDY_Y = 13 ; 

  double Y_min   = -8. ; 
  double Y_max   =  8. ; 
  int    Y_bins  =  48 ; 

  TH1D*                 H1D_DNDY_Y[_N_HISTOGRAMS_DNDY_Y] ; 
  H1D_DNDY_Y[0] = new TH1D("YH0", "charged_particle_eta_differential_yield", Y_bins, Y_min, Y_max);
  H1D_DNDY_Y[1] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH1");  H1D_DNDY_Y[1]->SetTitle("pion_plus");
  H1D_DNDY_Y[2] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH2");  H1D_DNDY_Y[2]->SetTitle("pion_minus");
  H1D_DNDY_Y[3] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH3");  H1D_DNDY_Y[3]->SetTitle("kaon_plus");
  H1D_DNDY_Y[4] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH4");  H1D_DNDY_Y[4]->SetTitle("kaon_minus");
  H1D_DNDY_Y[5] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH5");  H1D_DNDY_Y[5]->SetTitle("proton");
  H1D_DNDY_Y[6] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH6");  H1D_DNDY_Y[6]->SetTitle("anti_proton");
  H1D_DNDY_Y[7] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH7");  H1D_DNDY_Y[7]->SetTitle("lambda");
  H1D_DNDY_Y[8] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH8");  H1D_DNDY_Y[8]->SetTitle("anti_lambda");
  H1D_DNDY_Y[9] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH9");  H1D_DNDY_Y[9]->SetTitle("deuteron");
  H1D_DNDY_Y[10] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH10");  H1D_DNDY_Y[10]->SetTitle("anti_deuteron");
  H1D_DNDY_Y[11] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH11");  H1D_DNDY_Y[11]->SetTitle("triton");
  H1D_DNDY_Y[12] = (TH1D*) H1D_DNDY_Y[0]->Clone("YH12");  H1D_DNDY_Y[12]->SetTitle("He3");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      //double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Pt > pT_max || Pt < pT_min)
        continue ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        H1D_DNDY_Y[0]->Fill(Rap,1.);
      }

      if(PID == 211){
        H1D_DNDY_Y[1]->Fill(Rap,1.);
      }
      if(PID == -211){
        H1D_DNDY_Y[2]->Fill(Rap,1.);
      }
      if(PID == 321){
        H1D_DNDY_Y[3]->Fill(Rap,1.);
      }
      if(PID == -321){
        H1D_DNDY_Y[4]->Fill(Rap,1.);
      }
      if(PID == 2212){
        H1D_DNDY_Y[5]->Fill(Rap,1.);
      }
      if(PID == -2212){
        H1D_DNDY_Y[6]->Fill(Rap,1.);
      }
      if(PID == 3122){
        H1D_DNDY_Y[7]->Fill(Rap,1.);
      }
      if(PID == -3122){
        H1D_DNDY_Y[8]->Fill(Rap,1.);
      }
      if(PID == 11111){
        H1D_DNDY_Y[9]->Fill(Rap,1.);
      }
      if(PID == -11111){
        H1D_DNDY_Y[10]->Fill(Rap,1.);
      }
      if(PID == 33333){
        H1D_DNDY_Y[11]->Fill(Rap,1.);
      }
      if(PID == 44444){
        H1D_DNDY_Y[12]->Fill(Rap,1.);
      }

     } // particle loop
    } // event loop


  double dX = ( Y_max - Y_min ) / Y_bins ; 
  for(int ii=0; ii<_N_HISTOGRAMS_DNDY_Y; ii++)
        H1D_DNDY_Y[ii]->Scale(1.0 / (nEvents * dX));

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_DNDY_Y] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "11111", "-11111", "33333", "44444" };
  int        hadron_index[_N_HISTOGRAMS_DNDY_Y] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,     9  ,     10,      11,       12  };


 for(int ix =0; ix < _N_HISTOGRAMS_DNDY_Y; ix++){
   output_filename.str("");
   output_filename << "results/dndy_y-" << hadron_name[ix];
   output_filename << "_pT_" << pT_min << "_" << pT_max ;
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= H1D_DNDY_Y[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << H1D_DNDY_Y[hadron_index[ix]]->GetBinCenter(i) << "\t" << H1D_DNDY_Y[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << H1D_DNDY_Y[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_DNDY_Y ; ixx++){
     H1D_DNDY_Y[ixx]->Clear(); 
   }


 delete H1D_DNDY_Y[_N_HISTOGRAMS_DNDY_Y] ; 

}


void observables::calculate_invariant_yield_vs_pt( int yflag, double Rap_min, double Rap_max ){
    if(yflag > 1){
      std::cout << "Calculating dN/pTdpTdy within y range " 
                << Rap_min << "-" << Rap_max << " ... "
                << std::endl ;
    } 
    else{
      std::cout << "Calculating dN/pTdpTdeta within eta range " 
                << Rap_min << "-" << Rap_max << " ... "
                << std::endl ;
    }

  const int _N_HISTOGRAMS_INVYLD_PT_ = 12 ; 

  double pT_min  = 0.15 ; 
  double pT_max  = 4.5  ;
  int    pT_bins = 24   ;
  double dpT = ( pT_max - pT_min ) / pT_bins ; 
  double dY = ( Rap_max - Rap_min )  ; 

  TH1D*                 H1D_INVYLD_PT[_N_HISTOGRAMS_INVYLD_PT_] ; 
  H1D_INVYLD_PT[0] = new TH1D("PT0", "inavriant_yield_pion_plus", pT_bins, pT_min, pT_max);
  H1D_INVYLD_PT[1] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT1");  H1D_INVYLD_PT[1]->SetTitle("pion_minus");
  H1D_INVYLD_PT[2] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT2");  H1D_INVYLD_PT[2]->SetTitle("kaon_plus");
  H1D_INVYLD_PT[3] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT3");  H1D_INVYLD_PT[3]->SetTitle("kaon_minus");
  H1D_INVYLD_PT[4] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT4");  H1D_INVYLD_PT[4]->SetTitle("proton");
  H1D_INVYLD_PT[5] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT5");  H1D_INVYLD_PT[5]->SetTitle("anti_proton");
  H1D_INVYLD_PT[6] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT6");  H1D_INVYLD_PT[6]->SetTitle("Lambda");
  H1D_INVYLD_PT[7] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT7");  H1D_INVYLD_PT[7]->SetTitle("anti_Lambda");
  H1D_INVYLD_PT[8] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT8");  H1D_INVYLD_PT[8]->SetTitle("deutron");
  H1D_INVYLD_PT[9] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT9");  H1D_INVYLD_PT[9]->SetTitle("anti_deutron");
  H1D_INVYLD_PT[10] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT10");  H1D_INVYLD_PT[10]->SetTitle("triton");
  H1D_INVYLD_PT[11] = (TH1D*) H1D_INVYLD_PT[0]->Clone("PT11");  H1D_INVYLD_PT[11]->SetTitle("He3");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double W   = Event->get_particle(jj)->get_weight()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 

      if(PID == 211){
        H1D_INVYLD_PT[0]->Fill(Pt,1.0/Pt);
      }
      if(PID == -211){
        H1D_INVYLD_PT[1]->Fill(Pt,1.0/Pt);
      }
      if(PID == 321){
        H1D_INVYLD_PT[2]->Fill(Pt,1.0/Pt);
      }
      if(PID == -321){
        H1D_INVYLD_PT[3]->Fill(Pt,1.0/Pt);
      }
      if(PID == 2212){
        H1D_INVYLD_PT[4]->Fill(Pt,1.0/Pt);
      }
      if(PID == -2212){
        H1D_INVYLD_PT[5]->Fill(Pt,1.0/Pt);
      }
      if(PID == 3122){
        H1D_INVYLD_PT[6]->Fill(Pt,1.0/Pt);
      }
      if(PID == -3122){
        H1D_INVYLD_PT[7]->Fill(Pt,1.0/Pt);
      }
      if(PID == 11111){
        H1D_INVYLD_PT[8]->Fill(Pt,W/Pt);
      }
      if(PID == -11111){
        H1D_INVYLD_PT[9]->Fill(Pt,W/Pt);
      }
      if(PID == 33333){
        H1D_INVYLD_PT[10]->Fill(Pt,W/Pt);
      }
      if(PID == 44444){
        H1D_INVYLD_PT[11]->Fill(Pt,W/Pt);
      }


     } // particle loop
    } // event loop


  H1D_INVYLD_PT[0]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[1]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[2]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[3]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[4]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[5]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[6]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[7]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[8]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[9]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[10]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );
  H1D_INVYLD_PT[11]->Scale( 1.0 / ( nEvents  * 2.0 * TMath::Pi() * dY * dpT ) );

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_INVYLD_PT_] = { "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "11111", "-11111", "33333", "44444" };
  int        hadron_index[_N_HISTOGRAMS_INVYLD_PT_] = {   0,      1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7  ,     8   ,     9,       10,       11  };


 for(int ix =0; ix < _N_HISTOGRAMS_INVYLD_PT_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/dnptdptdy_pt-" << hadron_name[ix];
     output_filename << "_y_" << Rap_min << "_" << Rap_max ;
   }
  else{
     output_filename << "results/dnptdptdy_pt-" << hadron_name[ix];
     output_filename << "_eta_" << Rap_min << "_" << Rap_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= H1D_INVYLD_PT[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << H1D_INVYLD_PT[hadron_index[ix]]->GetBinCenter(i) << "\t" << H1D_INVYLD_PT[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << H1D_INVYLD_PT[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_INVYLD_PT_ ; ixx++){
     H1D_INVYLD_PT[ixx]->Clear(); 
   }


}


void observables::calculate_v1_vs_y_or_eta(int yflag, double psi1,  double pT_min, double pT_max ){
    if(yflag > 1){
      std::cout << "Calculating v1-y within pT range " 
                << pT_min << "-" << pT_max << " ... "
                << std::endl ;
    } 
    else{
      std::cout << "Calculating v1-eta within pT range " 
                << pT_min << "-" << pT_max << " ... "
                << std::endl ;
    }

   double Y_min   = -8. ; 
   double Y_max   =  8. ; 
   int    Y_bins  =  47 ; 
 
   const int   _N_HISTOGRAMS_V1_Y_ = 13 ; 
   TProfile*   PROFILE_V1_Y[_N_HISTOGRAMS_V1_Y_] ; 
   PROFILE_V1_Y[0] = new TProfile("PV0", "v1_y_charged_particle", Y_bins, Y_min, Y_max);
   PROFILE_V1_Y[1] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV1");  PROFILE_V1_Y[1]->SetTitle("pion_plus");
   PROFILE_V1_Y[2] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV2");  PROFILE_V1_Y[2]->SetTitle("pion_minus");
   PROFILE_V1_Y[3] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV3");  PROFILE_V1_Y[3]->SetTitle("kaon_plus");
   PROFILE_V1_Y[4] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV4");  PROFILE_V1_Y[4]->SetTitle("kaon_minus");
   PROFILE_V1_Y[5] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV5");  PROFILE_V1_Y[5]->SetTitle("proton");
   PROFILE_V1_Y[6] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV6");  PROFILE_V1_Y[6]->SetTitle("anti_proton");
   PROFILE_V1_Y[7] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV7");  PROFILE_V1_Y[7]->SetTitle("lambda");
   PROFILE_V1_Y[8] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV8");  PROFILE_V1_Y[8]->SetTitle("anti_lambda");
   PROFILE_V1_Y[9] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV9");  PROFILE_V1_Y[9]->SetTitle("deuteron");
   PROFILE_V1_Y[10] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV10");  PROFILE_V1_Y[10]->SetTitle("anti_deuteron");
   PROFILE_V1_Y[11] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV11");  PROFILE_V1_Y[11]->SetTitle("triton");
   PROFILE_V1_Y[12] = (TProfile*) PROFILE_V1_Y[0]->Clone("PV12");  PROFILE_V1_Y[12]->SetTitle("He3");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ;
      double W   = Event->get_particle(jj)->get_weight()   ;  
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Pt = sqrt( Px * Px + Py * Py )  ;

      if(Pt > pT_max || Pt < pT_min)
        continue ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 

      double Rap ;


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

     double v1 = Px / Pt ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_V1_Y[0]->Fill(Rap,v1);
      }

      if(PID == 211){
        PROFILE_V1_Y[1]->Fill(Rap,v1);
      }
      if(PID == -211){
        PROFILE_V1_Y[2]->Fill(Rap,v1);
      }
      if(PID == 321){
        PROFILE_V1_Y[3]->Fill(Rap,v1);
      }
      if(PID == -321){
        PROFILE_V1_Y[4]->Fill(Rap,v1);
      }
      if(PID == 2212){
        PROFILE_V1_Y[5]->Fill(Rap,v1);
      }
      if(PID == -2212){
        PROFILE_V1_Y[6]->Fill(Rap,v1);
      }
      if(PID == 3122){
        PROFILE_V1_Y[7]->Fill(Rap,v1);
      }
      if(PID == -3122){
        PROFILE_V1_Y[8]->Fill(Rap,v1);
      }
      if(PID == 11111){
        PROFILE_V1_Y[9]->Fill(Rap,v1,W);
      }
      if(PID == -11111){
        PROFILE_V1_Y[10]->Fill(Rap,v1,W);
      }
      if(PID == 33333){
        PROFILE_V1_Y[11]->Fill(Rap,v1,W);
      }
      if(PID == 44444){
        PROFILE_V1_Y[12]->Fill(Rap,v1,W);
      }

     } // particle loop
    } // event loop

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V1_Y_] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "11111", "-11111", "33333", "44444" };
  int        hadron_index[_N_HISTOGRAMS_V1_Y_] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,    9   ,     10  ,   11   ,    12   };

 for(int ix =0; ix < _N_HISTOGRAMS_V1_Y_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v1_y-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
   }
  else{
     output_filename << "results/v1_eta-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_V1_Y[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V1_Y[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }


 // Slope calculation //
 for(int ix =0; ix < _N_HISTOGRAMS_V1_Y_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v1_y-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
   }
  else{
     output_filename << "results/v1_eta-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
  }
   output_filename << "_slope.dat";

   double xx_val[20] ;
   double yy_val[20] ;
   double yy_up_err[20] ;
   double yy_low_err[20] ;
   for(int id=0; id<20; id++){
    xx_val[id] = 0. ; 
    yy_val[id] = 0. ; 
    yy_up_err[id] = 0. ; 
    yy_low_err[id] = 0. ; 
   }
  int npoints = 0 ; 

   for(int i=1; i<= PROFILE_V1_Y[hadron_index[ix]]->GetNbinsX(); i++){
      if( fabs( PROFILE_V1_Y[hadron_index[ix]]->GetBinCenter(i) ) > 1.1 )
        continue ;
      xx_val[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinCenter(i)    ;
      yy_val[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i)   ;
      if( PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) > 0 ){
        yy_up_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) + 
                          PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ; 
        yy_low_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) - 
                          PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ; 
      }
      else{
        yy_up_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) - 
                          PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ; 
        yy_low_err[npoints] = PROFILE_V1_Y[hadron_index[ix]]->GetBinContent(i) + 
                          PROFILE_V1_Y[hadron_index[ix]]->GetBinError(i)   ; 
      }

      npoints += 1 ; 
   }

   double slope_val = fit_a_straight_line_and_get_slope(npoints,xx_val,yy_val) ;  
   double slope_up_err = fit_a_straight_line_and_get_slope(npoints,xx_val,yy_up_err) ;
   double slope_low_err = fit_a_straight_line_and_get_slope(npoints,xx_val,yy_low_err) ;
   mFile.open(output_filename.str().c_str(), std::ios::out );  
   mFile << npoints << "  " << slope_val << "  " << (slope_up_err - slope_val) << "  " << ( slope_val - slope_low_err ) << endl ;
           // Actual slope Error = Slope found from  points with ( value + error ) - (value)
   mFile.close();
 
 }
 // slope calculation end //

   for(int ixx=0; ixx < _N_HISTOGRAMS_V1_Y_ ; ixx++){
      PROFILE_V1_Y[ixx]->Clear(); 
   }


}





void observables::calculate_v2_vs_y_or_eta(int yflag, double psi2,  double pT_min, double pT_max ){
    if(yflag > 1){
      std::cout << "Calculating v2-y within pT range " 
                << pT_min << "-" << pT_max << " ... "
                << std::endl ;
    } 
    else{
      std::cout << "Calculating v2-eta within pT range " 
                << pT_min << "-" << pT_max << " ... "
                << std::endl ;
    }


   double Y_min   = -8. ; 
   double Y_max   =  8. ; 
   int    Y_bins  =  81 ; 
 
   const int   _N_HISTOGRAMS_V2_Y_ = 13 ; 
   TProfile*   PROFILE_V2_Y[_N_HISTOGRAMS_V2_Y_] ; 
   PROFILE_V2_Y[0] = new TProfile("PV0", "v2_y_charged_particle", Y_bins, Y_min, Y_max);
   PROFILE_V2_Y[1] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV1");  PROFILE_V2_Y[1]->SetTitle("pion_plus");
   PROFILE_V2_Y[2] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV2");  PROFILE_V2_Y[2]->SetTitle("pion_minus");
   PROFILE_V2_Y[3] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV3");  PROFILE_V2_Y[3]->SetTitle("kaon_plus");
   PROFILE_V2_Y[4] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV4");  PROFILE_V2_Y[4]->SetTitle("kaon_minus");
   PROFILE_V2_Y[5] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV5");  PROFILE_V2_Y[5]->SetTitle("proton");
   PROFILE_V2_Y[6] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV6");  PROFILE_V2_Y[6]->SetTitle("anti_proton");
   PROFILE_V2_Y[7] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV7");  PROFILE_V2_Y[7]->SetTitle("lambda");
   PROFILE_V2_Y[8] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV8");  PROFILE_V2_Y[8]->SetTitle("anti_lambda");
   PROFILE_V2_Y[9] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV9");  PROFILE_V2_Y[9]->SetTitle("deuteron");
   PROFILE_V2_Y[10] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV10");  PROFILE_V2_Y[10]->SetTitle("anti_deuteron");
   PROFILE_V2_Y[11] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV11");  PROFILE_V2_Y[11]->SetTitle("triton");
   PROFILE_V2_Y[12] = (TProfile*) PROFILE_V2_Y[0]->Clone("PV12");  PROFILE_V2_Y[12]->SetTitle("He3");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Pt = sqrt( Px * Px + Py * Py )  ;

      if(Pt > pT_max || Pt < pT_min)
        continue ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 

      double Rap ;


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

     double v2 = ( Px * Px - Py * Py ) / ( Pt * Pt ) ; 

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_V2_Y[0]->Fill(Rap,v2);
      }

      if(PID == 211){
        PROFILE_V2_Y[1]->Fill(Rap,v2);
      }
      if(PID == -211){
        PROFILE_V2_Y[2]->Fill(Rap,v2);
      }
      if(PID == 321){
        PROFILE_V2_Y[3]->Fill(Rap,v2);
      }
      if(PID == -321){
        PROFILE_V2_Y[4]->Fill(Rap,v2);
      }
      if(PID == 2212){
        PROFILE_V2_Y[5]->Fill(Rap,v2);
      }
      if(PID == -2212){
        PROFILE_V2_Y[6]->Fill(Rap,v2);
      }
      if(PID == 3122){
        PROFILE_V2_Y[7]->Fill(Rap,v2);
      }
      if(PID == -3122){
        PROFILE_V2_Y[8]->Fill(Rap,v2);
      }
      if(PID == 11111){
        PROFILE_V2_Y[9]->Fill(Rap,v2);
      }
      if(PID == -11111){
        PROFILE_V2_Y[10]->Fill(Rap,v2);
      }
      if(PID == 33333){
        PROFILE_V2_Y[11]->Fill(Rap,v2);
      }
      if(PID == 44444){
        PROFILE_V2_Y[12]->Fill(Rap,v2);
      }


     } // particle loop
    } // event loop

  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V2_Y_] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122" , "11111", "-11111", "33333", "44444" };
  int        hadron_index[_N_HISTOGRAMS_V2_Y_] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8   ,    9   ,     10   ,   11  ,    12   };

 for(int ix =0; ix < _N_HISTOGRAMS_V2_Y_; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v2_y-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
   }
  else{
     output_filename << "results/v2_eta-" << hadron_name[ix];
     output_filename << "_pT_" << pT_min << "_" << pT_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_V2_Y[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V2_Y[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V2_Y[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_V2_Y[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_V2_Y_ ; ixx++){
      PROFILE_V2_Y[ixx]->Clear(); 
   }


}



void observables::calculate_v2_pt( int yflag, double Rap_min, double Rap_max ){
    if(yflag > 1){
      std::cout << "Calculating v2-pT within y range " 
                << Rap_min << "-" << Rap_max << " ... "
                << std::endl ;
    } 
    else{
      std::cout << "Calculating v2-pT within eta range " 
                << Rap_min << "-" << Rap_max << " ... "
                << std::endl ;
    }

  const int _N_HISTOGRAMS_V2_PT = 13 ; 
  double v2_pt_bins[13] = {0.01, 0.1, 0.15, 0.3, 0.5, 0.75,1.0, 1.25, 1.5 , 1.75, 2.1, 2.5, 3.0} ; 

  TProfile*                 PROFILE_V2_PT[_N_HISTOGRAMS_V2_PT] ; 
  PROFILE_V2_PT[0] = new TProfile("V2PT00", "pion_pT_differential_v2", 12, v2_pt_bins);
  PROFILE_V2_PT[1] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT0");  PROFILE_V2_PT[1]->SetTitle("pion_plus");
  PROFILE_V2_PT[2] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT1");  PROFILE_V2_PT[2]->SetTitle("pion_minus");
  PROFILE_V2_PT[3] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT2");  PROFILE_V2_PT[3]->SetTitle("kaon_plus");
  PROFILE_V2_PT[4] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT3");  PROFILE_V2_PT[4]->SetTitle("kaon_minus");
  PROFILE_V2_PT[5] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT4");  PROFILE_V2_PT[5]->SetTitle("proton");
  PROFILE_V2_PT[6] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT5");  PROFILE_V2_PT[6]->SetTitle("anti_proton");
  PROFILE_V2_PT[7] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT6");  PROFILE_V2_PT[7]->SetTitle("Lambda");
  PROFILE_V2_PT[8] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT7");  PROFILE_V2_PT[8]->SetTitle("anti_Lambda");
  PROFILE_V2_PT[9] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT8");  PROFILE_V2_PT[9]->SetTitle("deuteron");
  PROFILE_V2_PT[10] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT9");  PROFILE_V2_PT[9]->SetTitle("anti_deuteron");
  PROFILE_V2_PT[11] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT10");  PROFILE_V2_PT[9]->SetTitle("triton");
  PROFILE_V2_PT[12] = (TProfile*) PROFILE_V2_PT[0]->Clone("V2PT11");  PROFILE_V2_PT[9]->SetTitle("He3");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 
      double v2 = ( Px * Px - Py * Py ) / ( Pt * Pt ) ;

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_V2_PT[0]->Fill(Pt,v2);
      }

      if(PID == 211){
        PROFILE_V2_PT[1]->Fill(Pt,v2);
      }
      if(PID == -211){
        PROFILE_V2_PT[2]->Fill(Pt,v2);
      }
      if(PID == 321){
        PROFILE_V2_PT[3]->Fill(Pt,v2);
      }
      if(PID == -321){
        PROFILE_V2_PT[4]->Fill(Pt,v2);
      }
      if(PID == 2212){
        PROFILE_V2_PT[5]->Fill(Pt,v2);
      }
      if(PID == -2212){
        PROFILE_V2_PT[6]->Fill(Pt,v2);
      }
      if(PID == 3122){
        PROFILE_V2_PT[7]->Fill(Pt,v2);
      }
      if(PID == -3122){
        PROFILE_V2_PT[8]->Fill(Pt,v2);
      }
      if(PID == 11111){
        PROFILE_V2_PT[9]->Fill(Pt,v2);
      }
      if(PID == -11111){
        PROFILE_V2_PT[10]->Fill(Pt,v2);
      }
      if(PID == 33333){
        PROFILE_V2_PT[11]->Fill(Pt,v2);
      }
      if(PID == 44444){
        PROFILE_V2_PT[12]->Fill(Pt,v2);
      }


     } // particle loop
    } // event loop


  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V2_PT] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "11111", "-11111", "33333", "44444" };
  int        hadron_index[_N_HISTOGRAMS_V2_PT] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8   ,   9   ,    10   ,    11   ,   12   };


 for(int ix =0; ix < _N_HISTOGRAMS_V2_PT; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v2_pt-" << hadron_name[ix];
     output_filename << "_y_" << Rap_min << "_" << Rap_max ;
   }
  else{
     output_filename << "results/v2_pt-" << hadron_name[ix];
     output_filename << "_eta_" << Rap_min << "_" << Rap_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_V2_PT[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V2_PT[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V2_PT[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_V2_PT[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_V2_PT ; ixx++){
     PROFILE_V2_PT[ixx]->Clear(); 
   }


}



void observables::calculate_v1_pt( int yflag, double Rap_min, double Rap_max ){
    if(yflag > 1){
      std::cout << "Calculating v1-pT within y range " 
                << Rap_min << "-" << Rap_max << " ... "
                << std::endl ;
    } 
    else{
      std::cout << "Calculating v1-pT within eta range " 
                << Rap_min << "-" << Rap_max << " ... "
                << std::endl ;
    }

  const int _N_HISTOGRAMS_V1_PT = 13 ; 
  double v1_pt_bins[13] = {0.01, 0.1, 0.15, 0.3, 0.5, 0.75,1.0, 1.25, 1.5 , 1.75, 2.1, 2.5, 3.0} ; 

  TProfile*                 PROFILE_V1_PT[_N_HISTOGRAMS_V1_PT] ; 
  PROFILE_V1_PT[0] = new TProfile("V2PT0", "pion_pT_differential_v2", 12, v1_pt_bins);
  PROFILE_V1_PT[1] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT0");  PROFILE_V1_PT[1]->SetTitle("pion_plus");
  PROFILE_V1_PT[2] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT1");  PROFILE_V1_PT[2]->SetTitle("pion_minus");
  PROFILE_V1_PT[3] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT2");  PROFILE_V1_PT[3]->SetTitle("kaon_plus");
  PROFILE_V1_PT[4] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT3");  PROFILE_V1_PT[4]->SetTitle("kaon_minus");
  PROFILE_V1_PT[5] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT4");  PROFILE_V1_PT[5]->SetTitle("proton");
  PROFILE_V1_PT[6] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT5");  PROFILE_V1_PT[6]->SetTitle("anti_proton");
  PROFILE_V1_PT[7] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT6");  PROFILE_V1_PT[7]->SetTitle("Lambda");
  PROFILE_V1_PT[8] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT7");  PROFILE_V1_PT[8]->SetTitle("anti_Lambda");
  PROFILE_V1_PT[9] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT8");  PROFILE_V1_PT[9]->SetTitle("deuteron");
  PROFILE_V1_PT[10] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT9");  PROFILE_V1_PT[10]->SetTitle("anti_deuteron");
  PROFILE_V1_PT[11] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT10");  PROFILE_V1_PT[11]->SetTitle("triton");
  PROFILE_V1_PT[12] = (TProfile*) PROFILE_V1_PT[0]->Clone("V2PT11");  PROFILE_V1_PT[12]->SetTitle("He3");

  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event = rif->get_event(ii) ; 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      double Px  = Event->get_particle(jj)->get_px()  ; 
      double Py  = Event->get_particle(jj)->get_py()  ; 
      double Pz  = Event->get_particle(jj)->get_pz()  ; 
      double E   = Event->get_particle(jj)->get_e()   ; 
      double P = sqrt( Px * Px + Py * Py + Pz * Pz )  ;
      double Rap ; 

      if( fabs(E-Pz) < 1E-10 )
        continue ; 

      if( fabs(P-Pz) < 1E-10 )
        continue ; 


      if (yflag > 0 )
        Rap = 0.5 * TMath::Log( ( E + Pz ) / (E - Pz) );
      else
        Rap = 0.5 * TMath::Log( ( P + Pz ) / (P - Pz) );

      double Pt = sqrt( Px * Px + Py * Py )  ;
      if(Rap > Rap_max || Rap < Rap_min)
        continue ; 
      double v1 = Px / Pt ; ;

     if(PID ==  211  || PID == -211  || PID ==  321 || PID == -321  || 
        PID == 2212  || PID == -2212 || PID == 3122 || PID == -3122 ||
        PID == 3222  || PID == -3222 || PID == 3112 || PID == -3112 ||
        PID == 3312  || PID == -3312 || PID == 3334 || PID == -3334  ){
        PROFILE_V1_PT[0]->Fill(Pt,v1);
      }

      if(PID == 211){
        PROFILE_V1_PT[1]->Fill(Pt,v1);
      }
      if(PID == -211){
        PROFILE_V1_PT[2]->Fill(Pt,v1);
      }
      if(PID == 321){
        PROFILE_V1_PT[3]->Fill(Pt,v1);
      }
      if(PID == -321){
        PROFILE_V1_PT[4]->Fill(Pt,v1);
      }
      if(PID == 2212){
        PROFILE_V1_PT[5]->Fill(Pt,v1);
      }
      if(PID == -2212){
        PROFILE_V1_PT[6]->Fill(Pt,v1);
      }
      if(PID == 3122){
        PROFILE_V1_PT[7]->Fill(Pt,v1);
      }
      if(PID == -3122){
        PROFILE_V1_PT[8]->Fill(Pt,v1);
      }
      if(PID == 11111){
        PROFILE_V1_PT[9]->Fill(Pt,v1);
      }
      if(PID == -11111){
        PROFILE_V1_PT[10]->Fill(Pt,v1);
      }
      if(PID == 33333){
        PROFILE_V1_PT[11]->Fill(Pt,v1);
      }
      if(PID == 44444){
        PROFILE_V1_PT[12]->Fill(Pt,v1);
      }

     } // particle loop
    } // event loop


  std::ofstream mFile;
  std::stringstream output_filename;

  std::string hadron_name[_N_HISTOGRAMS_V1_PT] = {"hpm", "211", "-211", "321", "-321", "2212", "-2212", "3122", "-3122", "11111", "-11111", "33333", "44444" };
  int        hadron_index[_N_HISTOGRAMS_V1_PT] = {  0,     1 ,     2 ,    3 ,     4 ,     5 ,      6 ,     7 ,      8  ,    9   ,     10  ,    11  ,    12   };


 for(int ix =0; ix < _N_HISTOGRAMS_V1_PT; ix++){
   output_filename.str("");

   if(yflag > 0 ){
     output_filename << "results/v1_pt-" << hadron_name[ix];
     output_filename << "_y_" << Rap_min << "_" << Rap_max ;
   }
  else{
     output_filename << "results/v1_pt-" << hadron_name[ix];
     output_filename << "_eta_" << Rap_min << "_" << Rap_max ;
  }
   output_filename << ".dat";

   mFile.open(output_filename.str().c_str(), std::ios::out );
   for(int i=1; i<= PROFILE_V1_PT[hadron_index[ix]]->GetNbinsX(); i++){
       mFile << PROFILE_V1_PT[hadron_index[ix]]->GetBinCenter(i) << "\t" << PROFILE_V1_PT[hadron_index[ix]]->GetBinContent(i) 
              << "\t" << PROFILE_V1_PT[hadron_index[ix]]->GetBinError(i) << std::endl;
   }
   mFile.close();
 }

   for(int ixx=0; ixx < _N_HISTOGRAMS_V1_PT ; ixx++){
     PROFILE_V1_PT[ixx]->Clear(); 
   }


}




double observables::fit_a_straight_line_and_get_slope(int n, double *x, double *y){
  double sum_x     = 0. ; 
  double sum_x_sqr = 0. ;
  double sum_y     = 0. ; 
  double sum_xy    = 0. ; 

  for(int i=0; i<n; i++){
    sum_x += x[i] ;
    sum_x_sqr += ( x[i] * x[i] ) ; 
    sum_y += y[i] ;
    sum_xy += ( x[i] * y[i] ) ;
  }

  return ( n * sum_xy - sum_x * sum_y ) / ( n * sum_x_sqr - sum_x * sum_x ) ;


}


void observables::perform_coalescence_of_p_n_to_form_deuteron(){
  int deutron_count = 0 ; 
  double proton_index_list[1000];
  double neutron_index_list[1000];
  double proton_coal_flag[1000];
  double neutron_coal_flag[1000];
  int num_of_protons = 0 ; 
  int num_of_neutrons = 0 ;
  for(int ii=0; ii<1000; ii++){
    proton_index_list[ii] = 0 ; 
    neutron_index_list[ii] = 0 ; 
    proton_coal_flag[ii] = 0 ; 
    neutron_coal_flag[ii] = 0 ; 
  } 
  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event   = rif->get_event(ii) ; 
    num_of_protons  = 0 ; 
    num_of_neutrons = 0 ;
    for(int ixx=0; ixx<1000; ixx++){
      proton_index_list[ixx] = 0 ; 
      neutron_index_list[ixx] = 0 ; 
      proton_coal_flag[ixx] = 0 ; 
      neutron_coal_flag[ixx] = 0 ; 
    } 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      if(PID == 2212 ){ 
        proton_index_list[num_of_protons] = jj ; 
        num_of_protons += 1 ; 
      }
      else if(PID == 2112 ){
        neutron_index_list[num_of_neutrons] = jj ;
        num_of_neutrons += 1 ; 
      } 
      else{
	continue ; 
      }
      
    } // particle loop
    
    // deutron production
    for(int ip=0; ip<num_of_protons; ip++){
      for(int in=0; in<num_of_neutrons; in++){
	
        int proton_index = proton_index_list[ip] ; 
        int neutron_index = neutron_index_list[in] ;

        if( proton_coal_flag[ip]>0 || neutron_coal_flag[in]>0 ){
          continue ; 
        } 
	
        double proton_px = Event->get_particle(proton_index)->get_px() ; 
        double proton_py = Event->get_particle(proton_index)->get_py() ; 
        double proton_pz = Event->get_particle(proton_index)->get_pz() ; 
        double proton_e  = Event->get_particle(proton_index)->get_e() ; 
        double proton_x  = Event->get_particle(proton_index)->get_x() ; 
        double proton_y  = Event->get_particle(proton_index)->get_y() ; 
        double proton_z  = Event->get_particle(proton_index)->get_z() ; 
        double proton_t  = Event->get_particle(proton_index)->get_t() ; 
	
        double neutron_px = Event->get_particle(neutron_index)->get_px() ; 
        double neutron_py = Event->get_particle(neutron_index)->get_py() ; 
        double neutron_pz = Event->get_particle(neutron_index)->get_pz() ; 
        double neutron_e  = Event->get_particle(neutron_index)->get_e() ; 
        double neutron_x  = Event->get_particle(neutron_index)->get_x() ; 
        double neutron_y  = Event->get_particle(neutron_index)->get_y() ; 
        double neutron_z  = Event->get_particle(neutron_index)->get_z() ; 
        double neutron_t  = Event->get_particle(neutron_index)->get_t() ;
	
        // determine deutron momentum
        double deutron_px = proton_px + neutron_px ; 
        double deutron_py = proton_py + neutron_py ; 
        double deutron_pz = proton_pz + neutron_pz ; 
        double deutron_e  = proton_e  + neutron_e  ; 
	
	
        // determine deutron velocity
        double deutron_vx =  deutron_px / deutron_e ; 
        double deutron_vy =  deutron_py / deutron_e ; 
        double deutron_vz =  deutron_pz / deutron_e ; 
	
        // proton and neutron momentum in the rest frame of deutron
        double proton_px_prime = -999 ; 
        double proton_py_prime = -1999 ; 
        double proton_pz_prime = -2999 ; 
        double proton_e_prime  = -3999 ; 
	
        double neutron_px_prime = -6999 ; 
        double neutron_py_prime = -7999 ; 
        double neutron_pz_prime = -8999 ; 
        double neutron_e_prime  = -9999 ; 
	
        lorentz_transformation_p(proton_e, proton_px, proton_py, proton_pz, deutron_vx,
                                 deutron_vy, deutron_vz, proton_e_prime, proton_px_prime,
                                 proton_py_prime, proton_pz_prime );
	
        lorentz_transformation_p(neutron_e, neutron_px, neutron_py, neutron_pz, deutron_vx,
                                 deutron_vy, deutron_vz, neutron_e_prime, neutron_px_prime,
                                 neutron_py_prime, neutron_pz_prime );    
	
        if (  ( proton_px_prime + neutron_px_prime ) > 1e-6  ) {
	  std::cout << "Something is wrong with cm frame: "
		    << "p1)x + p2)x  = " << ( proton_px_prime + neutron_px_prime ) << std::endl;
        }
        if (  ( proton_py_prime + neutron_py_prime ) > 1e-6  ) {
	  std::cout << "Something is wrong with cm frame: "
		    << "p1)y + p2)y  = " << ( proton_py_prime + neutron_py_prime ) << std::endl;
        }
        if (  ( proton_pz_prime + neutron_pz_prime ) > 1e-6  ) {
	  std::cout << "Something is wrong with cm frame: "
		    << "p1)z + p2)z  = " << ( proton_pz_prime + neutron_pz_prime ) << std::endl;
        }
	
        // proton and neutron position in the rest frame of deutron
        double proton_x_prime = -999 ; 
        double proton_y_prime = -1999 ; 
        double proton_z_prime = -2999 ; 
        double proton_t_prime = -3999 ; 
	
        double neutron_x_prime = -6999 ; 
        double neutron_y_prime = -7999 ; 
        double neutron_z_prime = -8999 ; 
        double neutron_t_prime = -9999 ; 
	
	
        lorentz_transformation_x(proton_t, proton_x, proton_y, proton_z, deutron_vx, deutron_vy, deutron_vz, 
				 proton_t_prime, proton_x_prime, proton_y_prime, proton_z_prime );
	
        lorentz_transformation_x(neutron_t, neutron_x, neutron_y, neutron_z, deutron_vx, deutron_vy, deutron_vz, 
				 neutron_t_prime, neutron_x_prime, neutron_y_prime, neutron_z_prime );      
	
	
        // determine proton velocity in CM frame(or in the deutron rest frame)
        double proton_vx_cm =  proton_px_prime / proton_e_prime ; 
        double proton_vy_cm =  proton_py_prime / proton_e_prime ; 
        double proton_vz_cm =  proton_pz_prime / proton_e_prime ; 
	
        // determine neutron velocity in CM frame(or in the deutron rest frame)
        double neutron_vx_cm = neutron_px_prime / neutron_e_prime ; 
        double neutron_vy_cm = neutron_py_prime / neutron_e_prime ; 
        double neutron_vz_cm = neutron_pz_prime / neutron_e_prime ; 
	
	
        // 3. Roll to the time, when the last hadron was born
        double tmax = std::max( proton_t_prime, neutron_t_prime );
        double proton_x_cm = proton_x_prime + ( tmax - proton_t_prime ) * proton_vx_cm ; 
        double proton_y_cm = proton_y_prime + ( tmax - proton_t_prime ) * proton_vy_cm ; 
        double proton_z_cm = proton_z_prime + ( tmax - proton_t_prime ) * proton_vz_cm ; 
	
        double neutron_x_cm = neutron_x_prime + ( tmax - neutron_t_prime ) * neutron_vx_cm ; 
        double neutron_y_cm = neutron_y_prime + ( tmax - neutron_t_prime ) * neutron_vy_cm ; 
        double neutron_z_cm = neutron_z_prime + ( tmax - neutron_t_prime ) * neutron_vz_cm ; 
	
        // generate a random test number and decide the possibility of sampling a deutron by comparing with wigner function.          
        double rho_1  = 1./ sqrt(2.) * ( proton_x_cm - neutron_x_cm ) ; 
        double rho_2  = 1./ sqrt(2.) * ( proton_y_cm - neutron_y_cm ) ; 
        double rho_3  = 1./ sqrt(2.) * ( proton_z_cm - neutron_z_cm ) ;
        //double rho_1  = ( proton_x_cm - neutron_x_cm ) ; 
        //double rho_2  = ( proton_y_cm - neutron_y_cm ) ; 
        //double rho_3  = ( proton_z_cm - neutron_z_cm ) ; 
        double rho_sq = rho_1 * rho_1 + rho_2 * rho_2 + rho_3 * rho_3 ;
	
        double prho_1  = sqrt(2.) * ( neutron_mass * proton_px_prime - proton_mass * neutron_px_prime ) / ( proton_mass + neutron_mass )   ; 
        double prho_2  = sqrt(2.) * ( neutron_mass * proton_py_prime - proton_mass * neutron_py_prime ) / ( proton_mass + neutron_mass )   ; 
        double prho_3  = sqrt(2.) * ( neutron_mass * proton_pz_prime - proton_mass * neutron_pz_prime ) / ( proton_mass + neutron_mass )   ; 
        //double prho_1  = 0.5 * ( proton_px_prime - neutron_px_prime )   ; 
        //double prho_2  = 0.5 * ( proton_py_prime - neutron_py_prime )   ; 
        //double prho_3  = 0.5 * ( proton_pz_prime - neutron_pz_prime )   ; 
        double prho_sq =  prho_1 * prho_1 + prho_2 * prho_2 + prho_3 * prho_3 ; 
	prho_sq *= (5.068*5.068) ; // converted GeV^2 to fm^{-2}
	
        // statistical factor g = 3/4
        double wigner_func = (3. / 4.) * 8. * exp( - ( rho_sq ) / ( deutron_sigma_rho * deutron_sigma_rho ) - prho_sq * deutron_sigma_rho * deutron_sigma_rho ) ;  
	
        
        /*
	  double test_random = 6. * rand->Rndm();
	  if( wigner_func > 6. ){
          std::cout << "[Warning] argument of exponential in wigner function is greater than 6 ..." << std::endl ; 
	  }
	  if( test_random > wigner_func ){
          continue ; 
	  }
	  else{      
          proton_coal_flag[ip] = 1 ;     
          neutron_coal_flag[in] = 1 ;     
          deutron_count += 1 ; 
          Event->add_particle(11111, 0.5*(proton_t+neutron_t), 0.5*(proton_x+neutron_x), 0.5*(proton_y+neutron_y), 0.5*(proton_z+neutron_z),
	  deutron_e, deutron_px, deutron_py, deutron_pz, 1.0) ;
	  }
        */
        
        if(wigner_func < 0.00001){
	  continue ; 
        }
        else{
          Event->add_particle(11111, 0.5*(proton_t+neutron_t), 0.5*(proton_x+neutron_x), 
                              0.5*(proton_y+neutron_y), 0.5*(proton_z+neutron_z),
                              deutron_e, deutron_px, deutron_py, deutron_pz, wigner_func );
        }
        
      } // loop on neutron
    } // loop on proton
    
  } // event loop
  std::cout << "deutron_count = " << deutron_count << std::endl ; 
  std::cout << "velocity is greater than 1 for " << velocity_failure_counter 
            << " times ..." << std::endl ; 
  
}



void observables::perform_coalescence_of_antip_antin_to_form_antideuteron(){
  // Throughoutt this function the variable named with proton 
  // should be replace by anti-proton and neutron should be 
  // replaced by anti-neutron.
  int deutron_count = 0 ; 
  double proton_index_list[1000];
  double neutron_index_list[1000];
  double proton_coal_flag[1000];
  double neutron_coal_flag[1000];
  int num_of_protons = 0 ; 
  int num_of_neutrons = 0 ;
  for(int ii=0; ii<1000; ii++){
    proton_index_list[ii] = 0 ; 
    neutron_index_list[ii] = 0 ; 
    proton_coal_flag[ii] = 0 ; 
    neutron_coal_flag[ii] = 0 ; 
  } 
  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event   = rif->get_event(ii) ; 
    num_of_protons  = 0 ; 
    num_of_neutrons = 0 ;
    for(int ixx=0; ixx<1000; ixx++){
      proton_index_list[ixx] = 0 ; 
      neutron_index_list[ixx] = 0 ; 
      proton_coal_flag[ixx] = 0 ; 
      neutron_coal_flag[ixx] = 0 ; 
    } 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      if(PID == -2212 ){ 
        proton_index_list[num_of_protons] = jj ; 
        num_of_protons += 1 ; 
      }
      else if(PID == -2112 ){
        neutron_index_list[num_of_neutrons] = jj ;
        num_of_neutrons += 1 ; 
      } 
      else{
	continue ; 
      }
      
    } // particle loop
    
    // deutron production
    for(int ip=0; ip<num_of_protons; ip++){
      for(int in=0; in<num_of_neutrons; in++){
	
        int proton_index = proton_index_list[ip] ; 
        int neutron_index = neutron_index_list[in] ;

        if( proton_coal_flag[ip]>0 || neutron_coal_flag[in]>0 ){
          continue ; 
        } 
	
        double proton_px = Event->get_particle(proton_index)->get_px() ; 
        double proton_py = Event->get_particle(proton_index)->get_py() ; 
        double proton_pz = Event->get_particle(proton_index)->get_pz() ; 
        double proton_e  = Event->get_particle(proton_index)->get_e() ; 
        double proton_x  = Event->get_particle(proton_index)->get_x() ; 
        double proton_y  = Event->get_particle(proton_index)->get_y() ; 
        double proton_z  = Event->get_particle(proton_index)->get_z() ; 
        double proton_t  = Event->get_particle(proton_index)->get_t() ; 
	
        double neutron_px = Event->get_particle(neutron_index)->get_px() ; 
        double neutron_py = Event->get_particle(neutron_index)->get_py() ; 
        double neutron_pz = Event->get_particle(neutron_index)->get_pz() ; 
        double neutron_e  = Event->get_particle(neutron_index)->get_e() ; 
        double neutron_x  = Event->get_particle(neutron_index)->get_x() ; 
        double neutron_y  = Event->get_particle(neutron_index)->get_y() ; 
        double neutron_z  = Event->get_particle(neutron_index)->get_z() ; 
        double neutron_t  = Event->get_particle(neutron_index)->get_t() ;
	
        // determine deutron momentum
        double deutron_px = proton_px + neutron_px ; 
        double deutron_py = proton_py + neutron_py ; 
        double deutron_pz = proton_pz + neutron_pz ; 
        double deutron_e  = proton_e  + neutron_e  ; 
	
	
        // determine deutron velocity
        double deutron_vx =  deutron_px / deutron_e ; 
        double deutron_vy =  deutron_py / deutron_e ; 
        double deutron_vz =  deutron_pz / deutron_e ; 
	
        // proton and neutron momentum in the rest frame of deutron
        double proton_px_prime = -999 ; 
        double proton_py_prime = -1999 ; 
        double proton_pz_prime = -2999 ; 
        double proton_e_prime  = -3999 ; 
	
        double neutron_px_prime = -6999 ; 
        double neutron_py_prime = -7999 ; 
        double neutron_pz_prime = -8999 ; 
        double neutron_e_prime  = -9999 ; 
	
        lorentz_transformation_p(proton_e, proton_px, proton_py, proton_pz, deutron_vx,
                                 deutron_vy, deutron_vz, proton_e_prime, proton_px_prime,
                                 proton_py_prime, proton_pz_prime );
	
        lorentz_transformation_p(neutron_e, neutron_px, neutron_py, neutron_pz, deutron_vx,
                                 deutron_vy, deutron_vz, neutron_e_prime, neutron_px_prime,
                                 neutron_py_prime, neutron_pz_prime );    
	
        if (  ( proton_px_prime + neutron_px_prime ) > 1e-6  ) {
	  std::cout << "Something is wrong with cm frame: "
		    << "p1)x + p2)x  = " << ( proton_px_prime + neutron_px_prime ) << std::endl;
        }
        if (  ( proton_py_prime + neutron_py_prime ) > 1e-6  ) {
	  std::cout << "Something is wrong with cm frame: "
		    << "p1)y + p2)y  = " << ( proton_py_prime + neutron_py_prime ) << std::endl;
        }
        if (  ( proton_pz_prime + neutron_pz_prime ) > 1e-6  ) {
	  std::cout << "Something is wrong with cm frame: "
		    << "p1)z + p2)z  = " << ( proton_pz_prime + neutron_pz_prime ) << std::endl;
        }
	
        // proton and neutron position in the rest frame of deutron
        double proton_x_prime = -999 ; 
        double proton_y_prime = -1999 ; 
        double proton_z_prime = -2999 ; 
        double proton_t_prime = -3999 ; 
	
        double neutron_x_prime = -6999 ; 
        double neutron_y_prime = -7999 ; 
        double neutron_z_prime = -8999 ; 
        double neutron_t_prime = -9999 ; 
	
	
        lorentz_transformation_x(proton_t, proton_x, proton_y, proton_z, deutron_vx, deutron_vy, deutron_vz, 
				 proton_t_prime, proton_x_prime, proton_y_prime, proton_z_prime );
	
        lorentz_transformation_x(neutron_t, neutron_x, neutron_y, neutron_z, deutron_vx, deutron_vy, deutron_vz, 
				 neutron_t_prime, neutron_x_prime, neutron_y_prime, neutron_z_prime );      
	
	
        // determine proton velocity in CM frame(or in the deutron rest frame)
        double proton_vx_cm =  proton_px_prime / proton_e_prime ; 
        double proton_vy_cm =  proton_py_prime / proton_e_prime ; 
        double proton_vz_cm =  proton_pz_prime / proton_e_prime ; 
	
        // determine neutron velocity in CM frame(or in the deutron rest frame)
        double neutron_vx_cm = neutron_px_prime / neutron_e_prime ; 
        double neutron_vy_cm = neutron_py_prime / neutron_e_prime ; 
        double neutron_vz_cm = neutron_pz_prime / neutron_e_prime ; 
	
	
        // 3. Roll to the time, when the last hadron was born
        double tmax = std::max( proton_t_prime, neutron_t_prime );
        double proton_x_cm = proton_x_prime + ( tmax - proton_t_prime ) * proton_vx_cm ; 
        double proton_y_cm = proton_y_prime + ( tmax - proton_t_prime ) * proton_vy_cm ; 
        double proton_z_cm = proton_z_prime + ( tmax - proton_t_prime ) * proton_vz_cm ; 
	
        double neutron_x_cm = neutron_x_prime + ( tmax - neutron_t_prime ) * neutron_vx_cm ; 
        double neutron_y_cm = neutron_y_prime + ( tmax - neutron_t_prime ) * neutron_vy_cm ; 
        double neutron_z_cm = neutron_z_prime + ( tmax - neutron_t_prime ) * neutron_vz_cm ; 
	
        // generate a random test number and decide the possibility of sampling a deutron by comparing with wigner function.          
        double rho_1  = 1./ sqrt(2.) * ( proton_x_cm - neutron_x_cm ) ; 
        double rho_2  = 1./ sqrt(2.) * ( proton_y_cm - neutron_y_cm ) ; 
        double rho_3  = 1./ sqrt(2.) * ( proton_z_cm - neutron_z_cm ) ;
        double rho_sq = rho_1 * rho_1 + rho_2 * rho_2 + rho_3 * rho_3 ;
	
        double prho_1  = sqrt(2.) * ( neutron_mass * proton_px_prime - proton_mass * neutron_px_prime ) / ( proton_mass + neutron_mass )   ; 
        double prho_2  = sqrt(2.) * ( neutron_mass * proton_py_prime - proton_mass * neutron_py_prime ) / ( proton_mass + neutron_mass )   ; 
        double prho_3  = sqrt(2.) * ( neutron_mass * proton_pz_prime - proton_mass * neutron_pz_prime ) / ( proton_mass + neutron_mass )   ; 
        double prho_sq =  prho_1 * prho_1 + prho_2 * prho_2 + prho_3 * prho_3 ; 
	prho_sq *= (5.068*5.068) ; // converted GeV^2 to fm^{-2}
	
        // statistical factor g = 3/4
        double wigner_func = (3. / 4.) * 8. * exp( - ( rho_sq ) / ( deutron_sigma_rho * deutron_sigma_rho ) - prho_sq * deutron_sigma_rho * deutron_sigma_rho ) ;  
        
        if(wigner_func < 0.00001){
	  continue ; 
        }
        else{
          Event->add_particle(-11111, 0.5*(proton_t+neutron_t), 0.5*(proton_x+neutron_x), 
                              0.5*(proton_y+neutron_y), 0.5*(proton_z+neutron_z),
                              deutron_e, deutron_px, deutron_py, deutron_pz, wigner_func );
        }
        
      } // loop on neutron
    } // loop on proton
    
  } // event loop
  
}


void observables::perform_coalescence_of_p_n_n_to_form_triton(){
  int deutron_count = 0 ; 
  double proton_index_list[1000];
  double neutron_index_list[1000];
  double proton_coal_flag[1000];
  double neutron_coal_flag[1000];
  int num_of_protons = 0 ; 
  int num_of_neutrons = 0 ;
  for(int ii=0; ii<1000; ii++){
    proton_index_list[ii] = 0 ; 
    neutron_index_list[ii] = 0 ; 
    proton_coal_flag[ii] = 0 ; 
    neutron_coal_flag[ii] = 0 ; 
  } 
  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event   = rif->get_event(ii) ; 
    num_of_protons  = 0 ; 
    num_of_neutrons = 0 ;
    for(int ixx=0; ixx<1000; ixx++){
      proton_index_list[ixx] = 0 ; 
      neutron_index_list[ixx] = 0 ; 
      proton_coal_flag[ixx] = 0 ; 
      neutron_coal_flag[ixx] = 0 ; 
    } 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      if(PID == 2212 ){ 
        proton_index_list[num_of_protons] = jj ; 
        num_of_protons += 1 ; 
      }
      else if(PID == 2112 ){
        neutron_index_list[num_of_neutrons] = jj ;
        num_of_neutrons += 1 ; 
      } 
      else{
	continue ; 
      }
      
    } // particle loop
    
    // deutron production
    for(int ip=0; ip<num_of_protons; ip++){
      for(int in=0; in<num_of_neutrons; in++){
        for(int it=0; it<num_of_neutrons; it++){
	  
          if(it == in){
            continue ;}
	  
	  int proton_index = proton_index_list[ip] ; 
	  int neutron1_index = neutron_index_list[in] ;
	  int neutron2_index = neutron_index_list[it] ;
	  
	  double proton_px = Event->get_particle(proton_index)->get_px() ; 
	  double proton_py = Event->get_particle(proton_index)->get_py() ; 
	  double proton_pz = Event->get_particle(proton_index)->get_pz() ; 
	  double proton_e  = Event->get_particle(proton_index)->get_e() ; 
	  double proton_x  = Event->get_particle(proton_index)->get_x() ; 
	  double proton_y  = Event->get_particle(proton_index)->get_y() ; 
	  double proton_z  = Event->get_particle(proton_index)->get_z() ; 
	  double proton_t  = Event->get_particle(proton_index)->get_t() ; 
	  
	  double neutron1_px = Event->get_particle(neutron1_index)->get_px() ; 
	  double neutron1_py = Event->get_particle(neutron1_index)->get_py() ; 
	  double neutron1_pz = Event->get_particle(neutron1_index)->get_pz() ; 
	  double neutron1_e  = Event->get_particle(neutron1_index)->get_e() ; 
	  double neutron1_x  = Event->get_particle(neutron1_index)->get_x() ; 
	  double neutron1_y  = Event->get_particle(neutron1_index)->get_y() ; 
	  double neutron1_z  = Event->get_particle(neutron1_index)->get_z() ; 
	  double neutron1_t  = Event->get_particle(neutron1_index)->get_t() ;
	  
	  double neutron2_px = Event->get_particle(neutron2_index)->get_px() ; 
	  double neutron2_py = Event->get_particle(neutron2_index)->get_py() ; 
	  double neutron2_pz = Event->get_particle(neutron2_index)->get_pz() ; 
	  double neutron2_e  = Event->get_particle(neutron2_index)->get_e() ; 
	  double neutron2_x  = Event->get_particle(neutron2_index)->get_x() ; 
	  double neutron2_y  = Event->get_particle(neutron2_index)->get_y() ; 
	  double neutron2_z  = Event->get_particle(neutron2_index)->get_z() ; 
	  double neutron2_t  = Event->get_particle(neutron2_index)->get_t() ;
	  
	  // determine triton momentum
	  double triton_px = proton_px + neutron1_px + neutron2_px ; 
	  double triton_py = proton_py + neutron1_py + neutron2_py ; 
	  double triton_pz = proton_pz + neutron1_pz + neutron2_pz ; 
	  double triton_e  = proton_e  + neutron1_e  + neutron2_e  ; 
	  
	  // determine triton velocity
	  double triton_vx =  triton_px / triton_e ; 
	  double triton_vy =  triton_py / triton_e ; 
	  double triton_vz =  triton_pz / triton_e ; 
	  
	  // proton and neutron momentum in the rest frame of deutron
	  double proton_px_prime = -999 ; 
	  double proton_py_prime = -1999 ; 
	  double proton_pz_prime = -2999 ; 
	  double proton_e_prime  = -3999 ; 
	  
	  double neutron1_px_prime = -6999 ; 
	  double neutron1_py_prime = -7999 ; 
	  double neutron1_pz_prime = -8999 ; 
	  double neutron1_e_prime  = -9999 ; 
	  
	  double neutron2_px_prime = -6999 ; 
	  double neutron2_py_prime = -7999 ; 
	  double neutron2_pz_prime = -8999 ; 
	  double neutron2_e_prime  = -9999 ; 
	  
	  
	  lorentz_transformation_p(proton_e, proton_px, proton_py, proton_pz, triton_vx,
				   triton_vy, triton_vz, proton_e_prime, proton_px_prime,
				   proton_py_prime, proton_pz_prime );
	  
	  lorentz_transformation_p(neutron1_e, neutron1_px, neutron1_py, neutron1_pz, triton_vx,
				   triton_vy, triton_vz, neutron1_e_prime, neutron1_px_prime,
				   neutron1_py_prime, neutron1_pz_prime );  
	  
	  lorentz_transformation_p(neutron2_e, neutron2_px, neutron2_py, neutron2_pz, triton_vx,
				   triton_vy, triton_vz, neutron2_e_prime, neutron2_px_prime,
				   neutron2_py_prime, neutron2_pz_prime );    
	  
	  
	  if (  ( proton_px_prime + neutron1_px_prime + neutron2_px_prime ) > 1e-6  ) {
            std::cout << "Something is wrong with cm frame: "
		      << "p1)x + p2)x  = " << ( proton_px_prime + neutron1_px_prime + neutron2_px_prime ) << std::endl;
	  }
	  if (  ( proton_py_prime + neutron1_py_prime + neutron2_py_prime ) > 1e-6  ) {
            std::cout << "Something is wrong with cm frame: "
		      << "p1)y + p2)y  = " <<  ( proton_py_prime + neutron1_py_prime + neutron2_py_prime ) << std::endl;
	  }
	  if (  ( proton_pz_prime + neutron1_pz_prime+ neutron2_pz_prime ) > 1e-6  ) {
            std::cout << "Something is wrong with cm frame: "
		      << "p1)z + p2)z  = " << ( proton_pz_prime + neutron1_pz_prime+ neutron2_pz_prime ) << std::endl;
	  }
	  
	  // proton and neutron position in the rest frame of deutron
	  double proton_x_prime = -999 ; 
	  double proton_y_prime = -1999 ; 
	  double proton_z_prime = -2999 ; 
	  double proton_t_prime = -3999 ; 
	  
	  double neutron1_x_prime = -6999 ; 
	  double neutron1_y_prime = -7999 ; 
	  double neutron1_z_prime = -8999 ; 
	  double neutron1_t_prime = -9999 ; 
	  
	  double neutron2_x_prime = -6999 ; 
	  double neutron2_y_prime = -7999 ; 
	  double neutron2_z_prime = -8999 ; 
	  double neutron2_t_prime = -9999 ; 
	  
	  lorentz_transformation_x(proton_t, proton_x, proton_y, proton_z, triton_vx, triton_vy, triton_vz, 
                                   proton_t_prime, proton_x_prime, proton_y_prime, proton_z_prime );
	  
	  lorentz_transformation_x(neutron1_t, neutron1_x, neutron1_y, neutron1_z, triton_vx, triton_vy, triton_vz, 
                                   neutron1_t_prime, neutron1_x_prime, neutron1_y_prime, neutron1_z_prime );      
	  
	  lorentz_transformation_x(neutron2_t, neutron2_x, neutron2_y, neutron2_z, triton_vx, triton_vy, triton_vz, 
                                   neutron2_t_prime, neutron2_x_prime, neutron2_y_prime, neutron2_z_prime );      
	  
	  // determine proton velocity in CM frame(or in the deutron rest frame)
	  double proton_vx_cm =  proton_px_prime / proton_e_prime ; 
	  double proton_vy_cm =  proton_py_prime / proton_e_prime ; 
	  double proton_vz_cm =  proton_pz_prime / proton_e_prime ; 
	  
	  // determine neutron velocity in CM frame(or in the deutron rest frame)
	  double neutron1_vx_cm = neutron1_px_prime / neutron1_e_prime ; 
	  double neutron1_vy_cm = neutron1_py_prime / neutron1_e_prime ; 
	  double neutron1_vz_cm = neutron1_pz_prime / neutron1_e_prime ; 
	  
	  double neutron2_vx_cm = neutron2_px_prime / neutron2_e_prime ; 
	  double neutron2_vy_cm = neutron2_py_prime / neutron2_e_prime ; 
	  double neutron2_vz_cm = neutron2_pz_prime / neutron2_e_prime ; 
	  
	  // 3. Roll to the time, when the last hadron was born
	  double tmax = std::max( std::max( proton_t_prime, neutron1_t_prime ) , neutron2_t_prime ) ;
	  double proton_x_cm = proton_x_prime + ( tmax - proton_t_prime ) * proton_vx_cm ; 
	  double proton_y_cm = proton_y_prime + ( tmax - proton_t_prime ) * proton_vy_cm ; 
	  double proton_z_cm = proton_z_prime + ( tmax - proton_t_prime ) * proton_vz_cm ; 
	  
	  double neutron1_x_cm = neutron1_x_prime + ( tmax - neutron1_t_prime ) * neutron1_vx_cm ; 
	  double neutron1_y_cm = neutron1_y_prime + ( tmax - neutron1_t_prime ) * neutron1_vy_cm ; 
	  double neutron1_z_cm = neutron1_z_prime + ( tmax - neutron1_t_prime ) * neutron1_vz_cm ; 
	  
	  double neutron2_x_cm = neutron2_x_prime + ( tmax - neutron2_t_prime ) * neutron2_vx_cm ; 
	  double neutron2_y_cm = neutron2_y_prime + ( tmax - neutron2_t_prime ) * neutron2_vy_cm ; 
	  double neutron2_z_cm = neutron2_z_prime + ( tmax - neutron2_t_prime ) * neutron2_vz_cm ; 
	  
	  
	  // generate a random test number and decide the possibility of sampling a deutron by comparing with wigner function.
	  double rho_1  = sqrt(1./2.) * ( neutron1_x_cm - neutron2_x_cm ) ; 
	  double rho_2  = sqrt(1./2.) * ( neutron1_y_cm - neutron2_y_cm ) ; 
	  double rho_3  = sqrt(1./2.) * ( neutron1_z_cm - neutron2_z_cm ) ; 
	  double rho_sq = rho_1 * rho_1 + rho_2 * rho_2 + rho_3 * rho_3 ;
	  
	  double prho_1 = sqrt(2.) * ( neutron_mass * neutron1_px_prime - neutron_mass * neutron2_px_prime ) / ( neutron_mass + neutron_mass ) ; 
	  double prho_2 = sqrt(2.) * ( neutron_mass * neutron1_py_prime - neutron_mass * neutron2_py_prime ) / ( neutron_mass + neutron_mass ) ; 
	  double prho_3 = sqrt(2.) * ( neutron_mass * neutron1_pz_prime - neutron_mass * neutron2_pz_prime ) / ( neutron_mass + neutron_mass ) ;
	  double prho_sq = prho_1 * prho_1 + prho_2 * prho_2 + prho_3 * prho_3 ; 
          
	  
	  double lambda_1  =  sqrt(2./3.) 
	    * ( (neutron_mass*neutron1_x_cm + neutron_mass*neutron2_x_cm)/(neutron_mass + neutron_mass) - proton_x_cm  ) ; 
	  double lambda_2  =  sqrt(2./3.) 
	    * ( (neutron_mass*neutron1_y_cm + neutron_mass*neutron2_y_cm)/(neutron_mass + neutron_mass) - proton_y_cm  ) ; 
	  double lambda_3  =  sqrt(2./3.) 
	    * ( (neutron_mass*neutron1_z_cm + neutron_mass*neutron2_z_cm)/(neutron_mass + neutron_mass) - proton_z_cm  ) ; 
	  double lambda_sq = lambda_1 * lambda_1 + lambda_2 * lambda_2 + lambda_3 * lambda_3 ;
	  
	  
	  double m1pm2pm3 = proton_mass + neutron_mass + neutron_mass ; 
	  double plambda_1= sqrt(3./2.) / ( m1pm2pm3 ) 
	    * (  proton_mass*( neutron1_px_prime+neutron2_px_prime ) - (neutron_mass+neutron_mass)*proton_px_prime ) ; 
	  double plambda_2= sqrt(3./2.) / ( m1pm2pm3 ) 
	    * (  proton_mass*( neutron1_py_prime+neutron2_py_prime ) - (neutron_mass+neutron_mass)*proton_py_prime ) ; 
	  double plambda_3= sqrt(3./2.) / ( m1pm2pm3 ) 
	    * (  proton_mass*( neutron1_pz_prime+neutron2_pz_prime ) - (neutron_mass+neutron_mass)*proton_pz_prime ) ; 
	  double plambda_sq =  plambda_1 * plambda_1 + plambda_2 * plambda_2 + plambda_3 * plambda_3 ; 
	  plambda_sq *= (5.068*5.068) ; // converted GeV^2 to fm^{-2}
	  
	  // statistical factor g = 1/4
	  double wigner_func = (1. / 4.) * 8. * 8. * exp( - ( rho_sq ) / ( triton_sigma_rho * triton_sigma_rho ) 
							  - ( lambda_sq ) / ( triton_sigma_lambda * triton_sigma_lambda ) 
							  - prho_sq * triton_sigma_rho * triton_sigma_rho
							  - plambda_sq * triton_sigma_lambda * triton_sigma_lambda  ) ;  
	  
	  if(wigner_func < 0.00001){
	    continue ; 
	  }
	  else{
	    Event->add_particle(33333, 0.33*(proton_t+neutron1_t+neutron2_t), 0.33*(proton_x+neutron1_x+neutron2_x), 
				0.33*(proton_y+neutron1_y+neutron2_y), 0.33*(proton_z+neutron1_z+neutron2_z),
				triton_e, triton_px, triton_py, triton_pz, wigner_func );
	  }
	  
	} // loop on neutron
      } // loop on neutron
    } // loop on proton
  } // event loop
  
}


void observables::perform_coalescence_of_p_p_n_to_form_He3(){
  // Throughout this function the variable named with proton 
  // should be replace by neutron and viceversa.
  // Triton should be replaced by Helium-3

  // coalescence sigma_rho and sigma_lambda 
  // for He3 is taken to be 0.139 according 
  // to  Phys. Rev. C 95, 054913 .

  // Proton and neutron masses value are interchanged 
  // This helps to modify the coalescence function of triton 
  // and to make it for He3 . 
  double ne_mass = proton_mass ; 
  double pr_mass = neutron_mass ; 

  int deutron_count = 0 ; 
  double proton_index_list[1000];
  double neutron_index_list[1000];
  double proton_coal_flag[1000];
  double neutron_coal_flag[1000];
  int num_of_protons = 0 ; 
  int num_of_neutrons = 0 ;
  for(int ii=0; ii<1000; ii++){
    proton_index_list[ii] = 0 ; 
    neutron_index_list[ii] = 0 ; 
    proton_coal_flag[ii] = 0 ; 
    neutron_coal_flag[ii] = 0 ; 
  } 
  
  int nEvents = rif->get_event_buffer_size() ; 
  for(int ii=0; ii<nEvents; ii++){
    events* Event   = rif->get_event(ii) ; 
    num_of_protons  = 0 ; 
    num_of_neutrons = 0 ;
    for(int ixx=0; ixx<1000; ixx++){
      proton_index_list[ixx] = 0 ; 
      neutron_index_list[ixx] = 0 ; 
      proton_coal_flag[ixx] = 0 ; 
      neutron_coal_flag[ixx] = 0 ; 
    } 
    int nParticles = Event->get_multiplicity_of_the_event();
    for(int jj=0; jj<nParticles; jj++){
      int    PID = Event->get_particle(jj)->get_pid() ; 
      if(PID == 2112 ){  // NOTE : neutrons are filled in proton array . 
        proton_index_list[num_of_protons] = jj ; 
        num_of_protons += 1 ; 
      }
      else if(PID == 2212 ){ // NOTE : protons are filled in neutron array .
        neutron_index_list[num_of_neutrons] = jj ;
        num_of_neutrons += 1 ; 
      } 
      else{
	continue ; 
      }
      
    } // particle loop
    
    // deutron production
    for(int ip=0; ip<num_of_protons; ip++){
      for(int in=0; in<num_of_neutrons; in++){
        for(int it=0; it<num_of_neutrons; it++){
	  
          if(it == in){
            continue ;}
	  
	  int proton_index = proton_index_list[ip] ; 
	  int neutron1_index = neutron_index_list[in] ;
	  int neutron2_index = neutron_index_list[it] ;
	  
	  double proton_px = Event->get_particle(proton_index)->get_px() ; 
	  double proton_py = Event->get_particle(proton_index)->get_py() ; 
	  double proton_pz = Event->get_particle(proton_index)->get_pz() ; 
	  double proton_e  = Event->get_particle(proton_index)->get_e() ; 
	  double proton_x  = Event->get_particle(proton_index)->get_x() ; 
	  double proton_y  = Event->get_particle(proton_index)->get_y() ; 
	  double proton_z  = Event->get_particle(proton_index)->get_z() ; 
	  double proton_t  = Event->get_particle(proton_index)->get_t() ; 
	  
	  double neutron1_px = Event->get_particle(neutron1_index)->get_px() ; 
	  double neutron1_py = Event->get_particle(neutron1_index)->get_py() ; 
	  double neutron1_pz = Event->get_particle(neutron1_index)->get_pz() ; 
	  double neutron1_e  = Event->get_particle(neutron1_index)->get_e() ; 
	  double neutron1_x  = Event->get_particle(neutron1_index)->get_x() ; 
	  double neutron1_y  = Event->get_particle(neutron1_index)->get_y() ; 
	  double neutron1_z  = Event->get_particle(neutron1_index)->get_z() ; 
	  double neutron1_t  = Event->get_particle(neutron1_index)->get_t() ;
	  
	  double neutron2_px = Event->get_particle(neutron2_index)->get_px() ; 
	  double neutron2_py = Event->get_particle(neutron2_index)->get_py() ; 
	  double neutron2_pz = Event->get_particle(neutron2_index)->get_pz() ; 
	  double neutron2_e  = Event->get_particle(neutron2_index)->get_e() ; 
	  double neutron2_x  = Event->get_particle(neutron2_index)->get_x() ; 
	  double neutron2_y  = Event->get_particle(neutron2_index)->get_y() ; 
	  double neutron2_z  = Event->get_particle(neutron2_index)->get_z() ; 
	  double neutron2_t  = Event->get_particle(neutron2_index)->get_t() ;
	  
	  // determine triton momentum
	  double triton_px = proton_px + neutron1_px + neutron2_px ; 
	  double triton_py = proton_py + neutron1_py + neutron2_py ; 
	  double triton_pz = proton_pz + neutron1_pz + neutron2_pz ; 
	  double triton_e  = proton_e  + neutron1_e  + neutron2_e  ; 
	  
	  // determine triton velocity
	  double triton_vx =  triton_px / triton_e ; 
	  double triton_vy =  triton_py / triton_e ; 
	  double triton_vz =  triton_pz / triton_e ; 
	  
	  // proton and neutron momentum in the rest frame of deutron
	  double proton_px_prime = -999 ; 
	  double proton_py_prime = -1999 ; 
	  double proton_pz_prime = -2999 ; 
	  double proton_e_prime  = -3999 ; 
	  
	  double neutron1_px_prime = -6999 ; 
	  double neutron1_py_prime = -7999 ; 
	  double neutron1_pz_prime = -8999 ; 
	  double neutron1_e_prime  = -9999 ; 
	  
	  double neutron2_px_prime = -6999 ; 
	  double neutron2_py_prime = -7999 ; 
	  double neutron2_pz_prime = -8999 ; 
	  double neutron2_e_prime  = -9999 ; 
	  
	  
	  lorentz_transformation_p(proton_e, proton_px, proton_py, proton_pz, triton_vx,
				   triton_vy, triton_vz, proton_e_prime, proton_px_prime,
				   proton_py_prime, proton_pz_prime );
	  
	  lorentz_transformation_p(neutron1_e, neutron1_px, neutron1_py, neutron1_pz, triton_vx,
				   triton_vy, triton_vz, neutron1_e_prime, neutron1_px_prime,
				   neutron1_py_prime, neutron1_pz_prime );  
	  
	  lorentz_transformation_p(neutron2_e, neutron2_px, neutron2_py, neutron2_pz, triton_vx,
				   triton_vy, triton_vz, neutron2_e_prime, neutron2_px_prime,
				   neutron2_py_prime, neutron2_pz_prime );    
	  
	  
	  if (  ( proton_px_prime + neutron1_px_prime + neutron2_px_prime ) > 1e-6  ) {
            std::cout << "Something is wrong with cm frame: "
		      << "p1)x + p2)x  = " << ( proton_px_prime + neutron1_px_prime + neutron2_px_prime ) << std::endl;
	  }
	  if (  ( proton_py_prime + neutron1_py_prime + neutron2_py_prime ) > 1e-6  ) {
            std::cout << "Something is wrong with cm frame: "
		      << "p1)y + p2)y  = " <<  ( proton_py_prime + neutron1_py_prime + neutron2_py_prime ) << std::endl;
	  }
	  if (  ( proton_pz_prime + neutron1_pz_prime+ neutron2_pz_prime ) > 1e-6  ) {
            std::cout << "Something is wrong with cm frame: "
		      << "p1)z + p2)z  = " << ( proton_pz_prime + neutron1_pz_prime+ neutron2_pz_prime ) << std::endl;
	  }
	  
	  // proton and neutron position in the rest frame of deutron
	  double proton_x_prime = -999 ; 
	  double proton_y_prime = -1999 ; 
	  double proton_z_prime = -2999 ; 
	  double proton_t_prime = -3999 ; 
	  
	  double neutron1_x_prime = -6999 ; 
	  double neutron1_y_prime = -7999 ; 
	  double neutron1_z_prime = -8999 ; 
	  double neutron1_t_prime = -9999 ; 
	  
	  double neutron2_x_prime = -6999 ; 
	  double neutron2_y_prime = -7999 ; 
	  double neutron2_z_prime = -8999 ; 
	  double neutron2_t_prime = -9999 ; 
	  
	  lorentz_transformation_x(proton_t, proton_x, proton_y, proton_z, triton_vx, triton_vy, triton_vz, 
                                   proton_t_prime, proton_x_prime, proton_y_prime, proton_z_prime );
	  
	  lorentz_transformation_x(neutron1_t, neutron1_x, neutron1_y, neutron1_z, triton_vx, triton_vy, triton_vz, 
                                   neutron1_t_prime, neutron1_x_prime, neutron1_y_prime, neutron1_z_prime );      
	  
	  lorentz_transformation_x(neutron2_t, neutron2_x, neutron2_y, neutron2_z, triton_vx, triton_vy, triton_vz, 
                                   neutron2_t_prime, neutron2_x_prime, neutron2_y_prime, neutron2_z_prime );      
	  
	  // determine proton velocity in CM frame(or in the deutron rest frame)
	  double proton_vx_cm =  proton_px_prime / proton_e_prime ; 
	  double proton_vy_cm =  proton_py_prime / proton_e_prime ; 
	  double proton_vz_cm =  proton_pz_prime / proton_e_prime ; 
	  
	  // determine neutron velocity in CM frame(or in the deutron rest frame)
	  double neutron1_vx_cm = neutron1_px_prime / neutron1_e_prime ; 
	  double neutron1_vy_cm = neutron1_py_prime / neutron1_e_prime ; 
	  double neutron1_vz_cm = neutron1_pz_prime / neutron1_e_prime ; 
	  
	  double neutron2_vx_cm = neutron2_px_prime / neutron2_e_prime ; 
	  double neutron2_vy_cm = neutron2_py_prime / neutron2_e_prime ; 
	  double neutron2_vz_cm = neutron2_pz_prime / neutron2_e_prime ; 
	  
	  // 3. Roll to the time, when the last hadron was born
	  double tmax = std::max( std::max( proton_t_prime, neutron1_t_prime ) , neutron2_t_prime ) ;
	  double proton_x_cm = proton_x_prime + ( tmax - proton_t_prime ) * proton_vx_cm ; 
	  double proton_y_cm = proton_y_prime + ( tmax - proton_t_prime ) * proton_vy_cm ; 
	  double proton_z_cm = proton_z_prime + ( tmax - proton_t_prime ) * proton_vz_cm ; 
	  
	  double neutron1_x_cm = neutron1_x_prime + ( tmax - neutron1_t_prime ) * neutron1_vx_cm ; 
	  double neutron1_y_cm = neutron1_y_prime + ( tmax - neutron1_t_prime ) * neutron1_vy_cm ; 
	  double neutron1_z_cm = neutron1_z_prime + ( tmax - neutron1_t_prime ) * neutron1_vz_cm ; 
	  
	  double neutron2_x_cm = neutron2_x_prime + ( tmax - neutron2_t_prime ) * neutron2_vx_cm ; 
	  double neutron2_y_cm = neutron2_y_prime + ( tmax - neutron2_t_prime ) * neutron2_vy_cm ; 
	  double neutron2_z_cm = neutron2_z_prime + ( tmax - neutron2_t_prime ) * neutron2_vz_cm ; 
	  
	  
	  // generate a random test number and decide the possibility of sampling a deutron by comparing with wigner function.
	  double rho_1  = sqrt(1./2.) * ( neutron1_x_cm - neutron2_x_cm ) ; 
	  double rho_2  = sqrt(1./2.) * ( neutron1_y_cm - neutron2_y_cm ) ; 
	  double rho_3  = sqrt(1./2.) * ( neutron1_z_cm - neutron2_z_cm ) ; 
	  double rho_sq = rho_1 * rho_1 + rho_2 * rho_2 + rho_3 * rho_3 ;
	  
	  double prho_1 = sqrt(2.) * ( ne_mass * neutron1_px_prime - ne_mass * neutron2_px_prime ) / ( ne_mass + ne_mass ) ; 
	  double prho_2 = sqrt(2.) * ( ne_mass * neutron1_py_prime - ne_mass * neutron2_py_prime ) / ( ne_mass + ne_mass ) ; 
	  double prho_3 = sqrt(2.) * ( ne_mass * neutron1_pz_prime - ne_mass * neutron2_pz_prime ) / ( ne_mass + ne_mass ) ;
	  double prho_sq = prho_1 * prho_1 + prho_2 * prho_2 + prho_3 * prho_3 ; 
          
	  
	  double lambda_1  =  sqrt(2./3.) 
	    * ( (ne_mass*neutron1_x_cm + ne_mass*neutron2_x_cm)/(ne_mass + ne_mass) - proton_x_cm  ) ; 
	  double lambda_2  =  sqrt(2./3.) 
	    * ( (ne_mass*neutron1_y_cm + ne_mass*neutron2_y_cm)/(ne_mass + ne_mass) - proton_y_cm  ) ; 
	  double lambda_3  =  sqrt(2./3.) 
	    * ( (ne_mass*neutron1_z_cm + ne_mass*neutron2_z_cm)/(ne_mass + ne_mass) - proton_z_cm  ) ; 
	  double lambda_sq = lambda_1 * lambda_1 + lambda_2 * lambda_2 + lambda_3 * lambda_3 ;
	  
	  
	  double m1pm2pm3 = pr_mass + ne_mass + ne_mass ; 
	  double plambda_1= sqrt(3./2.) / ( m1pm2pm3 ) 
	    * (  pr_mass*( neutron1_px_prime+neutron2_px_prime ) - (ne_mass+ne_mass)*proton_px_prime ) ; 
	  double plambda_2= sqrt(3./2.) / ( m1pm2pm3 ) 
	    * (  pr_mass*( neutron1_py_prime+neutron2_py_prime ) - (ne_mass+ne_mass)*proton_py_prime ) ; 
	  double plambda_3= sqrt(3./2.) / ( m1pm2pm3 ) 
	    * (  pr_mass*( neutron1_pz_prime+neutron2_pz_prime ) - (ne_mass+ne_mass)*proton_pz_prime ) ; 
	  double plambda_sq =  plambda_1 * plambda_1 + plambda_2 * plambda_2 + plambda_3 * plambda_3 ; 
	  plambda_sq *= (5.068*5.068) ; // converted GeV^2 to fm^{-2}
	  
	  // statistical factor g = 1/4
	  double wigner_func = (1. / 4.) * 8. * 8. * exp( - ( rho_sq ) / ( He3_sigma_rho * He3_sigma_rho ) 
							  - ( lambda_sq ) / ( He3_sigma_lambda * He3_sigma_lambda ) 
							  - prho_sq * He3_sigma_rho * He3_sigma_rho
							  - plambda_sq * He3_sigma_lambda * He3_sigma_lambda  ) ;  
	  
	  if(wigner_func < 0.00001){
	    continue ; 
	  }
	  else{
	    Event->add_particle(44444, 0.33*(proton_t+neutron1_t+neutron2_t), 0.33*(proton_x+neutron1_x+neutron2_x), 
				0.33*(proton_y+neutron1_y+neutron2_y), 0.33*(proton_z+neutron1_z+neutron2_z),
				triton_e, triton_px, triton_py, triton_pz, wigner_func );
	  }
	  
	} // loop on neutron
      } // loop on neutron
    } // loop on proton
  } // event loop
  
}







// general lorentz transformation when velocity direction is arbitrary. 
void observables::lorentz_transformation_x( double t, double x, double y, double z, 
                                         double vx, double vy, double vz,
                                         double &tprime, double &xprime,
                                         double &yprime, double &zprime ){

  double v2    = vx * vx + vy * vy + vz * vz ; 

  if( sqrt(v2) > 1.01 ){
    std::cout << "unphysical velocity (Not corrected) ... " << std::endl ;
    std::cout << "v square = " << v2 << std::endl ; 
    std::cout << "v = " << sqrt(v2) << std::endl ; 
    exit(1); 
  }

  if( v2 > 0.99999 ) {
    velocity_failure_counter += 1 ; 
    vx = vx * sqrt(0.99999) / sqrt(v2) ; 
    vy = vy * sqrt(0.99999) / sqrt(v2) ; 
    vz = vz * sqrt(0.99999) / sqrt(v2) ; 
    v2 = 0.99999 ;
  }

  double gamma = 1. / sqrt(1. - v2 ) ; 
  double vdotr = vx * x + vy * y + vz * z ; 

  tprime = gamma * t - gamma * vdotr ; 
  xprime = x - gamma * vx * t  + ( gamma - 1. ) * vx / v2 * vdotr ;  
  yprime = y - gamma * vy * t  + ( gamma - 1. ) * vy / v2 * vdotr ;  
  zprime = z - gamma * vz * t  + ( gamma - 1. ) * vz / v2 * vdotr ;  
 
  if(std::isnan(tprime) || std::isinf(tprime) ){
   std::cout << " tprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ;
    exit(1);  
  }
  if(std::isnan(xprime) || std::isinf(xprime) ){
   std::cout << " xprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ;
    exit(1);  
  }

  if(std::isnan(yprime) || std::isinf(yprime) ){
   std::cout << " yprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ; 
    exit(1); 
  }

  if(std::isnan(zprime) || std::isinf(zprime) ){
   std::cout << " zprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ; 
    exit(1); 
  }


}


void observables::lorentz_transformation_p( double E, double px, double py, double pz, 
                                         double vx, double vy, double vz,
                                         double &Eprime, double &pxprime,
                                         double &pyprime, double &pzprime ){

  double v2    = vx * vx + vy * vy + vz * vz ;

  if( sqrt(v2) > 1.01 ){
    std::cout << "unphysical velocity (Not corrected) ... " << std::endl ;
    std::cout << "v square = " << v2 << std::endl ; 
    std::cout << "v = " << sqrt(v2) << std::endl ; 
    exit(1); 
  }

  if( v2 > 0.99999 ) {
    vx = vx * sqrt(0.99999) / sqrt(v2) ; 
    vy = vy * sqrt(0.99999) / sqrt(v2) ; 
    vz = vz * sqrt(0.99999) / sqrt(v2) ; 
    v2 = 0.99999 ;
  }


  double gamma = 1. / sqrt(1. - v2 ) ; 
  double vdotp = vx * px + vy * py + vz * pz ; 

  Eprime = gamma * E - gamma * vdotp ; 
  pxprime = px - gamma * vx * E  + ( gamma - 1. ) * vx / v2 * vdotp ;  
  pyprime = py - gamma * vy * E  + ( gamma - 1. ) * vy / v2 * vdotp ;  
  pzprime = pz - gamma * vz * E  + ( gamma - 1. ) * vz / v2 * vdotp ;  


  if(std::isnan(Eprime) || std::isinf(Eprime) ){
   std::cout << " Eprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ;
    exit(1);  
  }
  if(std::isnan(pxprime) || std::isinf(pxprime) ){
   std::cout << " pxprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ; 
    exit(1); 
  }

  if(std::isnan(pyprime) || std::isinf(pyprime) ){
   std::cout << " pyprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ; 
    exit(1); 
  }

  if(std::isnan(pzprime) || std::isinf(pzprime) ){
   std::cout << " pzprime = inf/nan " << std::endl ;  
   std::cout << "vx = " << vx << "  vy = " << vy << "  vz = " 
             << vz << "  v2 = " << v2 << "  gamma = " << gamma 
             << std::endl ;
    exit(1);  
  }



}














