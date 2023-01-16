#pragma once
#include <iostream>


class particles {

  public :
   particles();
   ~particles();
   inline void set_pid(int xx){
     PID = xx ; 
   }
   inline void set_e(double xx){
     E = xx ; 
   }
   inline void set_px(double xx){
     px = xx ; 
   }
   inline void set_py(double xx){
     py = xx ; 
   }
   inline void set_pz(double xx){
     pz = xx ; 
   }
   inline void set_t(double xx){
     t = xx ; 
   }
   inline void set_x(double xx){
     x = xx ; 
   }
   inline void set_y(double xx){
     y = xx ; 
   }
   inline void set_z(double xx){
     z = xx ; 
   }
   inline void set_weight(double xx){
     weight = xx ; 
   }


   inline int get_pid(){
     return PID ; 
   }
   inline double get_e(){
     return E ; 
   }
   inline double get_px(){
     return px ; 
   }
   inline double get_py(){
     return py ; 
   }
   inline double get_pz(){
     return pz ; 
   }
   inline double get_t(){
     return t ; 
   }
   inline double get_x(){
     return x ; 
   }
   inline double get_y(){
     return y ; 
   }
   inline double get_z(){
     return z ; 
   }
   inline double get_weight(){
     return weight ; 
   }


  private :
   int    PID ;
   double E   ;  
   double px  ;
   double py  ;
   double pz  ;
   double t   ; 
   double x   ; 
   double y   ; 
   double z   ; 
   double weight ; 

};
