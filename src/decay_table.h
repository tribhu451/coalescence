#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<chrono>
#include<stdio.h>
#include<string.h>
#include <vector>
#include "decay_channel.h"

using std::cout;
using std::ofstream;
using std::endl;
using std::string;
using std::cin;
using std::fstream;

class decay_table{

private :
     std::vector<decay_channel>	mdecay_channel;


public :
 decay_table() ; 
 ~decay_table() ;
  void add_decay_channel(decay_channel ) ; 
  decay_channel* get_decay_channel(int ) ;
  int get_channel_count() ; 


};

