#include "events.h"

events::events(){
}

events::~events(){
  particle_vector.clear() ; 
}

void events::add_particle(int pid, double t, double x, double y, double z, double e, double px, double py, double pz){
  particles* part = new particles();
  part->set_pid(pid);
  part->set_px(px);
  part->set_py(py);
  part->set_pz(pz);
  part->set_e(e);
  part->set_t(t);
  part->set_x(x);
  part->set_y(y);
  part->set_z(z);
  particle_vector.push_back(*part);
}


void events::add_particle(int pid, double t, double x, double y, double z, double e, double px, double py, double pz, double ww){
  particles* part = new particles();
  part->set_pid(pid);
  part->set_px(px);
  part->set_py(py);
  part->set_pz(pz);
  part->set_e(e);
  part->set_t(t);
  part->set_x(x);
  part->set_y(y);
  part->set_z(z);
  part->set_weight(ww);
  particle_vector.push_back(*part);
}



void events::add_particle(particles* part){
  particle_vector.push_back(*part);
}


