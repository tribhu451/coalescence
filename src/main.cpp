#include<iostream>
#include<fstream>
#include <vector>
#include "events.h"
#include "read_input_file.h"
#include "observables.h"
#include "inparams.h"
#include "read_pdg.h"
#include "inparams.h"


int main(int argc, char **argv){

  std::cout << ":::::::::::::::::::::::::::::::::::::::" << std::endl ; 
  std::cout << "::::::::::::   Coalescence   ::::::::::" << std::endl ; 
  std::cout << ":::::::::::::::::::::::::::::::::::::::" << std::endl ; 

  if(argc != 3){
    std::cout << "2 arguments required ..." << std::endl ;
    exit(1); 
  }

  input_paramters iparams ;
  read_parameters* r = new read_parameters();
  r->read_parameters_from_file(iparams, "input_parameters"); 

  read_pdg* RPDG = new read_pdg(iparams); 
  RPDG->read_and_store_particle_properties_with_decay_channels("PDG/pdg-urqmd_v3.3+_weak.dat");

  reso_decays* RD = new reso_decays(RPDG);
  read_input_file* RIF = new read_input_file(iparams, RD, argv[1], atof(argv[2]) );

  int nEvents = iparams.nEvents ; 
  

  // reading the input files //
  if(iparams.input_read_mode == 0 ){
    RIF->read_input_file_iSS_OSCAR(nEvents); }
  else if(iparams.input_read_mode == 1 ){
    RIF->read_particle_list_dat_from_urqmd(nEvents); }
  else if(iparams.input_read_mode == 2 ){
    RIF->read_particle_list_dat_from_urqmd_binary(nEvents); }
  else if(iparams.input_read_mode == 3 ){
    RIF->read_particle_list_dat_from_iSS_binary(nEvents); }
  else {
    std::cout << "reading mode not specified. Exiting ..." << std::endl ;  
    exit(1); 
  }

  observables* OBJ = new observables(iparams,RIF);

  std::cout << "Performing Coalescence ... " 
            << std::endl ; 
  OBJ->perform_coalescence();
  std::cout << "Coalescence done and deutrons " 
            << "are stored in particles list ... "
            << std::endl ; 

  OBJ->calculate_dnchdeta_eta(0.01,3);
  OBJ->calculate_dndy_y(0.01,3);
  OBJ->calculate_invariant_yield_vs_pt(1, -0.3, 0.3);
  OBJ->calculate_invariant_yield_vs_pt(1, -0.5, 0.5);
  OBJ->calculate_v1_vs_y_or_eta(1, 0, 0.15, 3 );
  OBJ->calculate_v1_vs_y_or_eta(1, 0, 0.2, 3 );
  OBJ->calculate_v1_vs_y_or_eta(1, 0, 0.4, 2 );
  OBJ->calculate_v2_pt( 0, -1.0, 1.0 );
  OBJ->calculate_v2_pt( 1, -0.5, 0.5 );
  
  return 0;
}





