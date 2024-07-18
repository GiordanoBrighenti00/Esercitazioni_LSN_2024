/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
  
  int equilibrium_time = 100; // number of steps necesssary to reach thermal equilibrium
  
  if (-1 > 0){ // Data Blocking: serve per il 7.1
    

    for(int i=0; i < SYS.get_nbl()+1; i++){ //loop over blocks
      if(i == 0){
        for(int j=0; j < equilibrium_time; j++){ //loop over steps in the first block: the transient-to-equilibrium one
          SYS.step();
        }
      }else{
        for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
          SYS.step();
          SYS.measure();
        }
        SYS.averages(i);
        SYS.block_reset(i);

        if(-1 < 0){ // prevent the creation of too many files if a high number of blocks is chosen
          //SYS.write_XYZ(i);  // writing the actual XYZ configuration at the end of each block
        }
      }  
    }
    SYS.finalize();
  }

  if(+1>0){
    for(int j=0; j < equilibrium_time; j++){ //loop over steps in the first block: the transient-to-equilibrium one
      SYS.step();
    }
    for(int i=0; i < SYS.get_nbl(); i++){
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
      }
      
      SYS.averages(i);
      SYS.block_reset(i);
      if(i%5000 == 0){
        cout << "Completati " << i << endl;
      }
    }
    SYS.finalize();
  }

  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
