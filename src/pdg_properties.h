#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<chrono>
#include<stdio.h>
#include<string.h>
#include "decay_table.h"



using std::cout;
using std::ofstream;
using std::endl;
using namespace std::chrono;
using std::string;
using std::cin;
using std::fstream;
using std::ios;
using std::istringstream;

class pdg_properties{

private :

int      pid ; 
double   mass ; 
string   name ;
double   width ;
double   gspin ;
int      net_baryon ;
int      net_strange ;
int      net_charm ;
int      net_bottom ;
double   gisospin ;
int      charge ; 
int      number_of_decay_channels ;
int      stability ; 
decay_table* mdecay_table ;  



public :

 pdg_properties();
 ~pdg_properties();

 inline void set_pid(int xx){
  pid = xx ; 
 }

 inline void set_mass(double xx){
  mass = xx ; 
 }

 inline void set_name(string xx){
  name = xx ; 
 }

 inline void set_width(double xx){
  width = xx ; 
 }

 inline void set_gspin(double xx){
  gspin = xx ; 
 }

 inline void set_net_baryon(int xx){
  net_baryon = xx ; 
 }

 inline void set_net_strange(int xx){
  net_strange = xx ; 
 }

 inline void set_net_charm(int xx){
  net_charm = xx ; 
 }

 inline void set_net_bottom(int xx){
  net_bottom = xx ; 
 }

 inline void set_gisospin(double xx){
  gisospin = xx ; 
 }

 inline void set_charge(int xx){
  charge = xx ; 
 }

 inline void set_number_of_decay_channels(int xx){
  number_of_decay_channels = xx ; 
 }

 inline void set_stability(int xx){
  stability = xx ; 
 }

 inline int get_pid(){
  return pid ; 
 }

 inline double get_mass(){
  return mass ; 
 }

 inline string get_name(){
  return name ; 
 }

 inline double get_width(){
  return width ; 
 }

 inline double get_gspin(){
  return gspin ; 
 }

 inline int get_net_baryon(){
  return net_baryon ; 
 }

 inline int get_net_strange(){
  return net_strange ; 
 }

 inline int get_net_charm(){
  return net_charm ; 
 }

 inline int get_net_bottom(){
  return net_bottom ; 
 }

 inline double get_gisospin(){
  return gisospin ; 
 }

 inline int get_charge(){
  return charge ; 
 }

 inline int get_number_of_decay_channels(){
  return number_of_decay_channels ; 
 }
 
 inline int is_stable(){
  return stability ; 
 } 

 decay_table* get_decay_table(){
  return mdecay_table ; 
 }




};























