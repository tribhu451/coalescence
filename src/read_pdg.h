#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include "pdg_properties.h"
#include "inparams.h"

using std::cout ; 
using std::endl ; 
using std::fstream ; 
using std::string ; 
using std::istringstream;


class read_pdg{
  
 private :
  
  int get_total_number_of_lines_of_the_file(std::string );
  input_paramters &iparam ; 
  istringstream* iss;
  char   buff[400];

  std::vector<pdg_properties> particle_table;
  std::map<int, int>    particle_pid_index_map;
  void check_the_existance_of_all_mentioned_daughter_particles_in_decay_channels() ; 
  int does_same_pid_exist_previously(int ) ; 
  
 public :
  
  read_pdg(input_paramters &iparam);
  ~read_pdg();


  void read_and_store_particle_properties_with_decay_channels(std::string ) ;

 
  // ::: on particle type :::
  inline void add_particle(pdg_properties* part_type){
   particle_table.push_back(*part_type);                                      // put in vector
   particle_pid_index_map[part_type->get_pid()] = particle_table.size() - 1;  // set in map ( pid -to- index )
  }

  pdg_properties* get_particle_from_index(int a_index){
    return &(particle_table[a_index]);
  }

  pdg_properties* get_particle_from_pid(int a_pid){
    return &(particle_table[particle_pid_index_map[a_pid]]);
  }

  inline int get_particle_index_from_pid(int a_pid){
    return particle_pid_index_map[a_pid];
  }

  inline int get_total_particle_count(){
    return particle_table.size();
  }

  int does_the_particle_exist(int a_pid){
   return particle_pid_index_map.count(a_pid);
  }





  
};
