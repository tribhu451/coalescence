#pragma once
#include "iostream"
#include "events.h"
#include "particles.h"
#include "read_pdg.h"
#include <vector>
#include <list>
#include <TRandom3.h>
#include <chrono>

using namespace std::chrono;

//#define M_PI 3.1415927

class reso_decays {


  public :
    reso_decays(read_pdg* )    ; 
    ~reso_decays()   ; 
    void perform_decays(events* Event) ;


  private :
   read_pdg* pdg ; 
   TRandom3* uniform_rand ; 
   void perform_two_body_decay(int mopid, double mom, double mow,
              double mot, double mox,  double moy,  double moz,
              double moe, double mopx, double mopy, double mopz, 
              double dt1m, double dt2m, int, int , particles* daughter1, particles* daughter2);
   void perform_three_body_decay(int mopid, double mom, double mow,
              double mot, double mox,  double moy,  double moz,
              double moe, double mopx, double mopy, double mopz, 
              double dt1m, double dt2m,  double dt3m , int , int , int , particles* daughter1, particles* daughter2, particles* daughter3 );


};
