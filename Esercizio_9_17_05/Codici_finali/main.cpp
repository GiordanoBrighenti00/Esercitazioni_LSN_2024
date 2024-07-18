#include <stdio.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <armadillo>
#include "random.h"
#include "Genetic.h"

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

   // Esercizio 9.1: Risolvere il problema del commesso viaggiatore per 34 città disposte su una circonferenza, 
   //utilizzando gli algoritmi genetici

   /// --- ORA I PUNTI SONO GENERATI ALL'INTERNO DEL QUADRATO ---

   arma::rowvec x_cities_sq(34, arma::fill::zeros);
   arma::rowvec y_cities_sq(34, arma::fill::zeros);

   cout << "# City" << "\t" << "x" << "\t" << "y" << endl;

   for(int i=0; i<34; i++){
       x_cities_sq[i] = rnd.Rannyu(-0.5,0.5);
       y_cities_sq[i] = rnd.Rannyu(-0.5,0.5);

       cout << i << "\t" <<  x_cities_sq[i] << "\t " << y_cities_sq[i] << endl;
   }

   // Generazione di 34 posizioni sulla circonferenza di raggio unitario con distribuzione uniforme
   int N1 = 100;

   Genetic G34SQ(N1, x_cities_sq, y_cities_sq);
   G34SQ.Set_Random(rnd);
   G34SQ.Set_Mutation_Probability(0.116);
   G34SQ.Set_Crossover_Probability(0.792);
   G34SQ.Initialize();
    cout << endl <<  "population 0";
   G34SQ.Print_Population();
   //G34.Get_Specimen(4);

   ofstream out3;
   out3.open("100_Stats_square.txt");
   out3 << "# Generation" << "\t" << "Best_cost" << "\t" << "Average_cost_best" << "\t" <<  "Stdev_cost_best" << endl;

   for(int i=0; i<1000; i++){
      cout << endl <<  "population # " << i+1;
      G34SQ.Evolution_Step();
      G34SQ.Print_Population();
      cout << endl;
      //G34.Get_Specimen(4);

      // --- Calcolo la statistica per stamparla in output
      int N_stat = floor(N1/2);

      double best_cost = G34SQ.Get_Cost(0);
      double average_cost_best = 0;
      double average_cost_best_2 = 0;
      double stdev_cost_best = 0;
      for(int j=0; j < N_stat; j++){
         average_cost_best = average_cost_best + G34SQ.Get_Cost(j);
      }
      average_cost_best = average_cost_best / (N_stat);

      for(int j=0; j < N_stat; j++){
         average_cost_best_2 = average_cost_best_2 + (G34SQ.Get_Cost(j) - average_cost_best)*(G34SQ.Get_Cost(j) - average_cost_best);
      }
      stdev_cost_best = sqrt(average_cost_best_2 / (N_stat-1) );

      if(i%20 == 0){
         out3 << i  << "\t" << best_cost << "\t" << average_cost_best << "\t" << stdev_cost_best << endl; 
      }

   }

   out3.close();

   ofstream out4;
   out4.open("100_Final_sequence_square.txt");
   for(int j = 0; j < 34; j++ ){
      int num = G34SQ.Get_City(0,j);
      out4 <<  num << "\t" << x_cities_sq[num-1] << "\t" << y_cities_sq[num-1] << endl;
   }
   out4.close();

   return 0;
   
   // --- ORA LO FACCIO SUL BORDO DELLA CIRCONFERENZA ---
   
   // Creazione di due vettori che contengono, rispettivamente, le coordinate x e y delle città

   arma::rowvec x_cities(34, arma::fill::zeros);
   arma::rowvec y_cities(34, arma::fill::zeros);

   cout << "# City" << "\t" << "x" << "\t" << "y" << endl;

   for(int i=0; i<34; i++){
       double theta = rnd.Rannyu(0,2*M_PI);
       x_cities[i] = cos(theta);
       y_cities[i] = sin(theta);

       cout << i << "\t" <<  x_cities[i] << "\t " << y_cities[i] << endl;
   }
  
   // Generazione di 34 posizioni sulla circonferenza di raggio unitario con distribuzione uniforme
   int N = 45;
   
   Genetic G34(N, x_cities, y_cities);
   G34.Set_Random(rnd);
   G34.Set_Mutation_Probability(0.17);
   G34.Set_Crossover_Probability(0.80);
   G34.Initialize();
    cout << endl <<  "population 0";
   G34.Print_Population();
   //G34.Get_Specimen(4);

   ofstream out;
   out.open("45_Stats.txt");
   out << "# Generation" << "\t" << "Best_cost" << "\t" << "Average_cost_best" << "\t" <<  "Stdev_cost_best" << endl;
   
   for(int i=0; i<1000; i++){
      cout << endl <<  "population # " << i+1;
      G34.Evolution_Step();
      G34.Print_Population();
      cout << endl;
      //G34.Get_Specimen(4);

      // --- Calcolo la statistica per stamparla in output
      int N_stat = floor(N/2);
      
      double best_cost = G34.Get_Cost(0);
      double average_cost_best = 0;
      double average_cost_best_2 = 0;
      double stdev_cost_best = 0;
      for(int j=0; j < N_stat; j++){
         average_cost_best = average_cost_best + G34.Get_Cost(j);
      }
      average_cost_best = average_cost_best / (N_stat);
      
      for(int j=0; j < N_stat; j++){
         average_cost_best_2 = average_cost_best_2 + (G34.Get_Cost(j) - average_cost_best)*(G34.Get_Cost(j) - average_cost_best);
      }
      stdev_cost_best = sqrt(average_cost_best_2 / (N_stat-1) );

      if(i%20 == 0){
         out << i  << "\t" << best_cost << "\t" << average_cost_best << "\t" << stdev_cost_best << endl; 
      }
      
   }
   
   out.close();

   ofstream out2;
   out2.open("45_Final_sequence.txt");
   for(int j = 0; j < 34; j++ ){
      int num = G34.Get_City(0,j);
      out2 <<  num << "\t" << x_cities[num-1] << "\t" << y_cities[num-1] << endl;
   }
   out2.close();

   // Generazione di 34 posizioni sulla circonferenza di raggio unitario con distribuzione uniforme
   N = 100;

   Genetic G34a(N, x_cities, y_cities);
   G34a.Set_Random(rnd);
   G34a.Set_Mutation_Probability(0.05);
   G34a.Set_Crossover_Probability(0.80);
   G34a.Initialize();
    cout << endl <<  "population 0";
   G34a.Print_Population();
   //G34.Get_Specimen(4);

   ofstream outa;
   outa.open("100_Stats.txt");
   outa << "# Generation" << "\t" << "Best_cost" << "\t" << "Average_cost_best" << "\t" <<  "Stdev_cost_best" << endl;

   for(int i=0; i<1000; i++){
      cout << endl <<  "population # " << i+1;
      G34a.Evolution_Step();
      G34a.Print_Population();
      cout << endl;
      //G34.Get_Specimen(4);

      // --- Calcolo la statistica per stamparla in output
      int N_stat = floor(N/2);

      double best_cost = G34a.Get_Cost(0);
      double average_cost_best = 0;
      double average_cost_best_2 = 0;
      double stdev_cost_best = 0;
      for(int j=0; j < N_stat; j++){
         average_cost_best = average_cost_best + G34a.Get_Cost(j);
      }
      average_cost_best = average_cost_best / (N_stat);

      for(int j=0; j < N_stat; j++){
         average_cost_best_2 = average_cost_best_2 + (G34a.Get_Cost(j) - average_cost_best)*(G34a.Get_Cost(j) - average_cost_best);
      }
      stdev_cost_best = sqrt(average_cost_best_2 / (N_stat-1) );

      if(i%20 == 0){
         outa << i  << "\t" << best_cost << "\t" << average_cost_best << "\t" << stdev_cost_best << endl; 
      }

   }

   outa.close();

   ofstream out2a;
   out2a.open("100_Final_sequence.txt");
   for(int j = 0; j < 34; j++ ){
      int num = G34a.Get_City(0,j);
      out2a <<  num << "\t" << x_cities[num-1] << "\t" << y_cities[num-1] << endl;
   }
   out2a.close();

   // Generazione di 34 posizioni sulla circonferenza di raggio unitario con distribuzione uniforme
   N = 100;

   Genetic G34b(N, x_cities, y_cities);
   G34b.Set_Random(rnd);
   G34b.Set_Mutation_Probability(0.1);
   G34b.Set_Crossover_Probability(0.80);
   G34b.Initialize();
    cout << endl <<  "population 0";
   G34a.Print_Population();
   //G34.Get_Specimen(4);

   ofstream outb;
   outb.open("100_alt_Stats.txt");
   outb << "# Generation" << "\t" << "Best_cost" << "\t" << "Average_cost_best" << "\t" <<  "Stdev_cost_best" << endl;

   for(int i=0; i<1000; i++){
      cout << endl <<  "population # " << i+1;
      G34b.Evolution_Step();
      G34b.Print_Population();
      cout << endl;
      //G34.Get_Specimen(4);

      // --- Calcolo la statistica per stamparla in output
      int N_stat = floor(N/2);

      double best_cost = G34b.Get_Cost(0);
      double average_cost_best = 0;
      double average_cost_best_2 = 0;
      double stdev_cost_best = 0;
      for(int j=0; j < N_stat; j++){
         average_cost_best = average_cost_best + G34b.Get_Cost(j);
      }
      average_cost_best = average_cost_best / (N_stat);

      for(int j=0; j < N_stat; j++){
         average_cost_best_2 = average_cost_best_2 + (G34b.Get_Cost(j) - average_cost_best)*(G34b.Get_Cost(j) - average_cost_best);
      }
      stdev_cost_best = sqrt(average_cost_best_2 / (N_stat-1) );

      if(i%20 == 0){
         outb << i  << "\t" << best_cost << "\t" << average_cost_best << "\t" << stdev_cost_best << endl; 
      }

   }

   outb.close();

   ofstream out2b;
   out2b.open("100_alt_Final_sequence.txt");
   for(int j = 0; j < 34; j++ ){
      int num = G34b.Get_City(0,j);
      out2b <<  num << "\t" << x_cities[num-1] << "\t" << y_cities[num-1] << endl;
   }
   out2b.close();
   

   
  return 0;
}