#pragma once
#include <iostream>
#include "particles.h"
#include <vector>


class events{
  public :
    events();
    ~events();
    void add_particle(int pid, double t, double x, double y, double z, double e, double px, double py, double pz);
    void add_particle(int pid, double t, double x, double y, double z, double e, double px, double py, double pz, double weight);
    void add_particle(particles* part) ; 

    particles* get_particle(int xx){
      return &particle_vector[xx] ; 
    }
  
   inline int get_multiplicity_of_the_event(){
     return particle_vector.size() ;
   }

   inline void clear_particle_vector(){
     particle_vector.clear() ; 
   }

  private :
    std::vector<particles> particle_vector ;

};
