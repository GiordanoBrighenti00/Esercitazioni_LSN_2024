#include <stdio.h>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include "random.h"
#include "simulator.h"

using namespace std;


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

  // Esercizio 8.2: Effettuare un simulated annealing per identificare i parametri (mu, sigma) che minimizzano il valore di aspettazione dell'energia. 
  //Partendo da un set di parametri imposti in ingresso, si calcola il valore <E> con Metropolis campionato su |psi(x)|^2 (ciclo interno).
  //Si propone una mossa di traslazione del punto (mu, beta), la cui probabilità di accettazione dipende dal confronto tra le probabilità di 
//Boltzmann associate alle energie <E> prima e dopo il passo. Si ripete iterativamente tale campionamento Metropolis (ciclo intermedio).
   // Si ripetono iterativamente i blocchi così definiti, variando per ogni nuovo batch la temperatura, accumulando statistica attraverso le medie globali (ciclo esterno)

// Parametri di partenza: argomenti della classe Simulator
  double sigma = -0.94;
  double mu = 1.59;

double T = 0.91; // Temperatura iniziale (normalizzata secondo Boltzmann)

  //double L = 1.25; // valore del passo ottimale che permette di campionare la psi
  double x_c = 0.25;

  Simulator SYS(mu, sigma);
  SYS.SetRandom(rnd);
  SYS.SetStartingPoint(x_c);
  SYS.SetStep(0.025, 0.025, 2.05); // Passo ideale x è 1.95
  SYS.Simulated_Annealing(T, 0.8, "output_instant.txt", "output_stats.txt", 8000, 200, 700, 2000);
  //SYS.Metro_Walk(20, 10000, 1500, "output_instant.txt", "output_stats.txt");
   // alpha = 0.8 è ideale
  return 0;
}