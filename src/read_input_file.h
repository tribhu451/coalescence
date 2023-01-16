#pragma once
#include <iostream>
#include <fstream>
#include "events.h"
#include <vector>
#include <sstream>
#include <string>
#include "inparams.h"
#include "reso_decays.h"
#include <chrono>

using namespace std::chrono;


using std::istringstream ; 
class read_input_file{
  public :

    read_input_file(input_paramters &iparam_, reso_decays*, std::string path_for_input_file_to_be_read__ , int num_of_file_sets_to_be_read__ ) ;
    ~read_input_file();
    void read_input_file_iSS_OSCAR(int );
    void read_particle_list_dat_from_urqmd(int );
    void read_particle_list_dat_from_urqmd_binary(int TotalEvents) ; 
    void read_particle_list_dat_from_iSS_binary(int TotalEvents) ; 

    events* get_event(int xx){
      return &event_vector[xx] ; 
    }

    inline int get_event_buffer_size(){
       return event_vector.size() ; 
    }

  private :
    istringstream* iss;
    input_paramters &iparam ; 
    reso_decays* resonance_decays ; 
    char buff[400];
    std::vector<events> event_vector ;
    int get_PID_from_urqmd_MCID(int mcid, int iso) ;
    
    std::string path_for_input_file_to_be_read ;
    int num_of_file_sets_to_be_read ; 
 


};
