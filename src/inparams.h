#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<chrono>
#include<stdio.h>
#include<string.h>

using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::cin;
using std::fstream;
using std::ios;
using std::istringstream;

typedef struct{

int nEvents ; 
int input_read_mode ; 
int include_weak_decay  ; 


} input_paramters ; 



class read_parameters{

 public :
  
  read_parameters();
  ~read_parameters();
  void read_parameters_from_file(input_paramters &iparam, string input_file_name);
  
};

