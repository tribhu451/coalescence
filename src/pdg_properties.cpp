#include "pdg_properties.h"

pdg_properties::pdg_properties() : pid(0),mass(0.), name ("Anonymous"),
width(0.), gspin(0.), net_baryon(0), net_strange(0), net_charm(0), 
net_bottom(0),gisospin(0), charge(0), number_of_decay_channels(0), stability(-1){
   mdecay_table = new decay_table() ; 
}

pdg_properties::~pdg_properties(){
}
