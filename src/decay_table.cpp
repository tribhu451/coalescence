#include "decay_table.h"

decay_table::decay_table(){
 mdecay_channel.clear() ; 
}

decay_table::~decay_table(){
 mdecay_channel.clear() ; 
}

void decay_table::add_decay_channel(decay_channel a_channel){
  mdecay_channel.push_back(a_channel);
 }

decay_channel* decay_table::get_decay_channel(int a_index ){
  return &(mdecay_channel[a_index]);
 }

int decay_table::get_channel_count(){
 return mdecay_channel.size() ; 
 }
