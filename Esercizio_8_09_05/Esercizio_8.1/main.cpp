#include <stdio.h>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include "random.h"
#include "simulator.h"

using namespace std;

// ---------- INIZIO DELLE FUNZIONI -------------

double min(double x, double y){
  if(x < y){
    return x;
  }else{
    return y;
  }
}

// --------- FINE DELLE FUNZIONI ----------------

int main(int argc, char *argv[]) {
  Random rnd;
  int seed[5];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   for(int i=0; i<20; i++){
      //cout << rnd.Rannyu() << endl;
   } 
  rnd.SaveSeed();

  // Esercizio 8.1: Campionare la distribuzione di probabilità Psi_{sigma, mu}(x) per poi calcolare il valore di aspettazione <H>:
  // Ovvero: bisogna integrare la funzione H psi/psi = V(x) - 0.5 * psi''(x)/psi(x) sulla misura p(x)*dx = |psi(x)|^2*dx
  // L'importance sampling è direttamente campionato con Metropolis

  double sigma = 0.12;
  double mu = 0.45;

  // Parametri del Metropolis
  int N1 = 2000;
  int N = 15000;
  int M = 20;
  double L = 1.25; // valore del passo ottimale che rende l'acceptance circa il 50%
  double x_c = 0.25;

  Simulator SYS(mu, sigma, N1, N, M , L);
  SYS.SetRandom(rnd);
  SYS.SetOutputFiles("output_positions.txt", "output_stats.txt");
  SYS.SetStartingPoint(x_c);
  SYS.Metro_Walk();

  return 0;
}