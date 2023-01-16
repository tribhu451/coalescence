#include "read_input_file.h"

read_input_file::read_input_file(input_paramters &iparam_, reso_decays* reso_decays_,
  std::string path_for_input_file_to_be_read__, int num_of_file_sets_to_be_read__) : 
iparam(iparam_), resonance_decays(reso_decays_){

path_for_input_file_to_be_read = path_for_input_file_to_be_read__ ;
num_of_file_sets_to_be_read = num_of_file_sets_to_be_read__ ;

}

read_input_file::~read_input_file(){
  event_vector.clear() ; 
}

void read_input_file::read_input_file_iSS_OSCAR(int TotalEvents) {

  std::cout << "reading OSCAR.DAT from iSS ..." 
              << std::endl ; 
  std::cout << "The file contains primordial particles info ..."
              << std::endl ; 
  std::string dummy ;
  int nch ;
  int index, pid ;
  double  x, y, e, px, py, pz, mass, t, z;

  std::ifstream file;
  file.open("OSCAR.DAT");
  if(!file){
    std::cout << "file not found." << std::endl;
    exit(1);
  }

  file.getline(buff,150) ;
  file.getline(buff,150) ;
  file.getline(buff,150) ;

  for(int ii=0; ii<TotalEvents; ii++){
    file.getline(buff,100);
    iss = new istringstream(buff);
    *iss >> dummy >> nch >> dummy >> dummy ;
    delete iss;
    std::cout << ii << "\t" << nch << std::endl ;
    events* Event = new events(); 
    for(int i=0; i<nch ; i++){
      file.getline(buff,300);
      iss = new istringstream(buff);
      *iss >> index >> pid >> px >> py >> pz >> e >> mass >> x >> y >> z >> t ;
      delete iss;
      // add the particle to the event
      Event->add_particle(pid, t, x, y, z,e, px, py, pz) ;
    } // particle loop

    event_vector.push_back(*Event);
  }

  if( iparam.include_weak_decay > 0 ){
    for(int ii=0; ii<TotalEvents; ii++){
       resonance_decays->perform_decays(get_event(ii)) ; 
     }
  }
  
  //std::cout << "Inside read_input_file, vector size : " << event_vector.size() << std::endl ; 
  //std::cout << "Inside read_input_file, vector size : " << &event_vector[0] << std::endl ; 
}


void read_input_file::read_particle_list_dat_from_urqmd(int TotalEvents){

  std::cout << "reading particle_list.dat from urqmd ..." 
              << std::endl ; 
 
  std::string dummy ;
  int nch ;
  double x, y, e, px, py, pz, mass, t, z;
  int temp_pid, temp_iso ; 

  std::ifstream file;
  file.open("particle_list.dat");
  if(!file){
    std::cout << "file not found." 
                << std::endl;
    exit(1);
  }

  for(int ii=0; ii<TotalEvents; ii++){

    for(int ij=0; ij<17; ij++){
      file.getline(buff,200) ;
    }

    file.getline(buff,100);
    iss = new istringstream(buff);
    *iss >> nch >> dummy ;
    delete iss;
    std::cout << ii << "\t" << nch << std::endl ;

    events* Event = new events();
 
    file.getline(buff,200) ; // read an useless line

    for(int i=0; i<nch ; i++){
      file.getline(buff,320);
      iss = new istringstream(buff);
      *iss       >> dummy >> dummy >> dummy >> dummy 
                 >> dummy >> dummy >> dummy >> dummy 
                 >> mass >> temp_pid >> temp_iso
                 >> dummy >> dummy >> dummy >> dummy
                 >> t >> x >> y >> z
                 >> e >> px >> py >> pz ;

      delete iss;
      // add the particle to the event
      Event->add_particle(get_PID_from_urqmd_MCID(temp_pid,temp_iso), t, x, y, z, e, px, py, pz) ;
    } // particle loop

    event_vector.push_back(*Event);
  }
 file.close() ; 

  if( iparam.include_weak_decay > 0 ){
    for(int ii=0; ii<TotalEvents; ii++){
       resonance_decays->perform_decays(get_event(ii)) ; 
     }
  }
  
}




void read_input_file::read_particle_list_dat_from_urqmd_binary(int TotalEvents){

  auto start = high_resolution_clock::now();

  std::cout << "reading particle_list.bin of urqmd ..." 
              << std::endl ; 
 
  std::string dummy ;
  int nch ;
  int pid ;
  double x, y, e, px, py, pz, mass, t, z ;
  int temp_iso ; 

  std::ifstream file;

  for(int input_file_index=0; input_file_index < num_of_file_sets_to_be_read ; input_file_index++ ){

   std::stringstream input_filename1;
   input_filename1.str(std::string());
   input_filename1 << path_for_input_file_to_be_read.c_str() ;
   input_filename1 << "/particle_list_set_";
   input_filename1 << input_file_index ;
   input_filename1 << ".bin";

   file.open(input_filename1.str().c_str(),std::ios::binary | std::ios::in);
   if(!file){
     std::cout << "file not found." 
                 << std::endl;
     exit(1);
   }


  std::cout << "reading file : " << input_filename1.str().c_str() << " ... " << std::endl ; 

  int ii=0 ; 
   while (ii < TotalEvents){
    file.read(reinterpret_cast<char *>(&nch), sizeof(int));
     if (file.eof()) break;

    if( (ii%500) == 0 ){
    std::cout << ii << "\t" << nch << std::endl ;
    }

    events* Event = new events();
 
    for(int i=0; i<nch ; i++){

        float particle_array[11];
          for (int i = 0; i < 11; i++) {
              float temp;
              file.read(reinterpret_cast<char *>(&temp), sizeof(float));
              particle_array[i] = temp;
          }
            mass = particle_array[0] ; 
            pid = particle_array[1] ; 
            temp_iso = particle_array[2] ; 
            t = particle_array[3] ; 
            x = particle_array[4] ; 
            y = particle_array[5] ; 
            z = particle_array[6] ; 
            e = particle_array[7] ; 
            px = particle_array[8] ; 
            py = particle_array[9] ; 
            pz = particle_array[10] ; 

            // only reading protons, neutrons and their anti-particles
            /*
            int score = 0 ; 
            if(get_PID_from_urqmd_MCID(pid,temp_iso) == 2212 )  { score += 1 ; }
            if(get_PID_from_urqmd_MCID(pid,temp_iso) == -2212 ) { score += 1 ; }
            if(get_PID_from_urqmd_MCID(pid,temp_iso) == 2112 )  { score += 1 ; }
            if(get_PID_from_urqmd_MCID(pid,temp_iso) == -2112 ) { score += 1 ; }
            if(score == 0 ) { continue ; }
            */

      // add the particle to the event
      Event->add_particle(get_PID_from_urqmd_MCID(pid,temp_iso), t, x, y, z, e, px, py, pz) ;
    } // particle loop

    event_vector.push_back(*Event);
    ii++ ; 
  }

  file.close() ;
 } // input_file_index loop


 
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);
  std::cout << "UrQMD binary file reading finished in " << duration.count() 
  	    << " sec.  ... "  << std::endl;


  if( iparam.include_weak_decay > 0 ){

    start = high_resolution_clock::now();
    std::cout << "Decaying resonances ... " << std::endl ; 

    for(int ievents=0; ievents<get_event_buffer_size(); ievents++){
       resonance_decays->perform_decays(get_event(ievents)) ; 
     }

    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "Resonances decay finished in " << duration.count() 
  	    << " sec.  ... "  << std::endl;

  }
  
}



void read_input_file::read_particle_list_dat_from_iSS_binary(int TotalEvents){

  auto start = high_resolution_clock::now();

  std::cout << "reading particle_list.bin of iSS ..." 
              << std::endl ; 
 
  std::string dummy ;
  int nch ;
  int pid ;
  double x, y, e, px, py, pz, mass, t, z ;


  std::ifstream file;

  for(int input_file_index=0; input_file_index < num_of_file_sets_to_be_read ; input_file_index++ ){

   std::stringstream input_filename1;
   input_filename1.str(std::string());
   input_filename1 << path_for_input_file_to_be_read.c_str() ;
   input_filename1 << "/primordial_particle_list_set_";
   input_filename1 << input_file_index ;
   input_filename1 << ".bin";

   file.open(input_filename1.str().c_str(),std::ios::binary | std::ios::in);
   if(!file){
     std::cout << "file not found." 
                 << std::endl;
     exit(1);
   }


  std::cout << "reading file : " << input_filename1.str().c_str() << " ... " << std::endl ; 

  int ii=0 ; 
   while (ii < TotalEvents){
    file.read(reinterpret_cast<char *>(&nch), sizeof(int));
     if (file.eof()) break;

    if( (ii%500) == 0 ){
    std::cout << ii << "\t" << nch << std::endl ;
    }

    events* Event = new events();
 
    for(int i=0; i<nch ; i++){

        float particle_array[10];
          for (int i = 0; i < 10; i++) {
              float temp;
              file.read(reinterpret_cast<char *>(&temp), sizeof(float));
              particle_array[i] = temp;
          }
            mass = particle_array[0] ; 
            pid = particle_array[1] ; 
            t = particle_array[2] ; 
            x = particle_array[3] ; 
            y = particle_array[4] ; 
            z = particle_array[5] ; 
            e = particle_array[6] ; 
            px = particle_array[7] ; 
            py = particle_array[8] ; 
            pz = particle_array[9] ; 

      // add the particle to the event
      Event->add_particle(pid, t, x, y, z, e, px, py, pz) ;
    } // particle loop

    event_vector.push_back(*Event);
    ii++ ; 
  }

  file.close() ;
 } // input_file_index loop


 
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);
  std::cout << "iSS binary file reading finished in " << duration.count() 
  	    << " sec.  ... "  << std::endl;


  if( iparam.include_weak_decay > 0 ){

    start = high_resolution_clock::now();
    std::cout << "Decaying resonances ... " << std::endl ; 

    for(int ievents=0; ievents<get_event_buffer_size(); ievents++){
       resonance_decays->perform_decays(get_event(ievents)) ; 
     }

    stop = high_resolution_clock::now();
    duration = duration_cast<seconds>(stop - start);
    std::cout << "Resonances decay finished in " << duration.count() 
  	    << " sec.  ... "  << std::endl;

  }
  
}



int read_input_file::get_PID_from_urqmd_MCID(int mcid, int iso){

  if (mcid == 101 && iso == 2 )
    return 211 ;                         // pi^{+}
  else if (mcid == 101 && iso == 0 )
    return 111 ;                         // pi^{0}
  else if (mcid == 101 && iso == -2 )
    return -211 ;                        // pi^{-}
  else if (mcid == 106 && iso == 1 )
    return 321 ;                         // k^{+}
  else if (mcid == 106 && iso == -1 )
    return 311 ;                         // k^{0}
  else if (mcid == -106 && iso == 1 )
    return -311 ;                        // anti k^{0}
  else if (mcid == -106 && iso == -1 )
    return -321 ;                        // k^{-}
  else if (mcid == 109 && iso == 0 )
    return 333 ;                         // phi(1020)
  else if (mcid == 102 && iso == 0 )
    return 221 ;                         // eta
  else if (mcid == 100 && iso == 0 )
    return 22 ;                          // photon
  else if (mcid == 1 && iso == 1 )
    return 2212 ;                        // p 
  else if (mcid == 1 && iso == -1 )
    return 2112 ;                        // n 
  else if (mcid == -1 && iso == -1 )
    return -2212 ;                       // anti p
  else if (mcid == -1 && iso == 1 )
    return -2112 ;                       // anti n
  else if (mcid == 40 && iso == 2 )
    return 3222 ;                        // Sigma^+
  else if (mcid == -40 && iso == -2 )
    return -3222 ;                       // anti Sigma^-
  else if (mcid == 40 && iso == 0 )
    return 3212 ;                        // Sigma^0
  else if (mcid == -40 && iso == 0 )
    return -3212 ;                       // anti Sigma^0
  else if (mcid == 40 && iso == -2 )
    return 3112 ;                        // Sigma^-
  else if (mcid == -40 && iso == 2 )
    return -3112 ;                       // anti Sigma^+
  else if (mcid == 49 && iso == -1 )
    return 3312 ;                        // Xi^-
  else if (mcid == 49 && iso == 1 )
    return 3322 ;                        // Xi^0
  else if (mcid == -49 && iso == -1 )
    return -3322 ;                       // anti Xi^0
  else if (mcid == -49 && iso == 1 )
    return -3312 ;                       // anti Xi^+
  else if (mcid == 27 && iso == 0 )
    return 3122 ;                        // Labmda
  else if (mcid == -27 && iso == 0 )
    return -3122 ;                       // anti Labmda
  else if (mcid == 55 && iso == 0 )
    return 3334 ;                        // Omega
  else if (mcid == -55 && iso == 0 )
    return -3334 ;                       // anti Omega
  else                  
    return 0 ; 

}



