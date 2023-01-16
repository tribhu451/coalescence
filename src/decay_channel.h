#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<chrono>
#include<stdio.h>
#include<string.h>
#include "decay_channel.h"

using std::cout;
using std::ofstream;
using std::endl;
using std::string;
using std::cin;
using std::fstream;

class decay_channel{

private :

double branching_ratio ;
int daughter_1 ; 
int daughter_2 ; 
int daughter_3 ; 
int daughter_4 ; 
int daughter_5 ; 


public :
 decay_channel(double , int, int, int , int, int ) ; 
 ~decay_channel() ; 

 double get_branching_ratio(){
  return branching_ratio ; 
 }

 double get_daughter_1(){
  return daughter_1 ; 
 }

 double get_daughter_2(){
  return daughter_2 ; 
 }

 double get_daughter_3(){
  return daughter_3 ; 
 }

 double get_daughter_4(){
  return daughter_4 ; 
 }

 double get_daughter_5(){
  return daughter_5 ; 
 }

 int is_two_body_decay(){
   if ( daughter_3 + daughter_4 + daughter_5  == 0 ) 
      return 1 ;
   else
      return 0 ; 
 }

};

