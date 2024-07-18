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

  // VERSION 1: STANDARD CODE: DO NOT CHANGE

  if(-1 > 0){
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
        if(j%10 == 0){
          //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          nconf++;
        }
      }
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
    SYS.finalize();
  }

  // VERSION 2.A: MODIFIED CODE: SIMULATION BEFORE EQUILIBRATION TIME
  if (-1 > 0){
    int equilibrium_time = 250; // Number of Verlet integration steps necessary to reach thermal equilibrium at T = 1.1

    for(int i=0; i < equilibrium_time; i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
        if(j%10 == 0){
          //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          nconf++;
        }
      }
      //SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
    //SYS.finalize();
  }

  // VERSION 2.B: MODIFIED CODE: SIMULATION AFTER EQUILIBRATION TIME
  if (-1 > 0){

    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
        if(j%10 == 0){
          //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          nconf++;
        }
      }
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
    SYS.finalize();
  }

  // VERSION 3: UNIFIED CODE: an initial block performs the transient to thermal equilibrium, followed by standard blocks that accumulate statistics
  if (+1 > 0){
    int equilibrium_time = 1; // number of steps necesssary to reach thermal equilibrium at T = 1.1

    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
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

        if(i < 100){ // prevent the creation of too many files if a high number of blocks is chosen
          SYS.write_XYZ(i);  // writing the actual XYZ configuration at the end of each block
        }
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
