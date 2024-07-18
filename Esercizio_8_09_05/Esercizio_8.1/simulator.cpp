#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "random.h"
#include "simulator.h"

using namespace std;

// Costruttore di default
Simulator :: Simulator(){}

// Costruttore con parametri in input
Simulator :: Simulator(double mu, double sigma, int N1, int N, int M, double L){
  _mu = mu;
  _sigma = sigma;
  _N1 = N1;
  _N = N;
  _M = M;
  _L = L;
}

// Distruttore di default
Simulator :: ~Simulator(){}

// --- INIZIO DEI METODI ---

void Simulator :: SetOutputFiles(string filename1, string filename2){
  _filename1 = filename1;
  _filename2 = filename2;
}

void Simulator :: SetRandom(Random& rnd){
  _rnd = rnd;
}

void Simulator :: SetStartingPoint(double x_c){
  _x = x_c;
}

void Simulator :: SetPosition(double x){
  _x = x;
}

double Simulator :: GetPosition(){
  return _x;
}

double Simulator :: V_Potential(double x){
  return x*x*x*x - (5.0)*x*x/(2.0);
}

double Simulator :: Psi_T(double x){  // Funzione d'onda Psi_{mu,sigma}(x)
  double exp1 = exp(-(x - _mu)*(x - _mu)/(2.0 * _sigma * _sigma));
  double exp2 = exp(-(x + _mu)*(x + _mu)/(2.0 * _sigma * _sigma));
  return exp1 + exp2;
}

double Simulator :: Lapl_T(double x){ // Applicazione dell'energia cinetica a Psi_{mu,sigma}(x) (Schrodinger)
  double h_cut = 1.0;
  double mass = 1.0;

  double coeff1 = -(1.0)/(_sigma * _sigma) + (x - _mu)*(x - _mu)/(_sigma * _sigma * _sigma * _sigma);
  double coeff2 = -(1.0)/(_sigma * _sigma) + (x + _mu)*(x + _mu)/(_sigma * _sigma * _sigma * _sigma);
  double exp1 = exp(-(x - _mu)*(x - _mu)/(2.0 * _sigma * _sigma));
  double exp2 = exp(-(x + _mu)*(x + _mu)/(2.0 * _sigma * _sigma));

  return - (h_cut * h_cut)*(coeff1*exp1 + coeff2*exp2)/(2.0*mass);
}

double Simulator :: E_integrand(double x){
  return this->V_Potential(x) + this->Lapl_T(x)/(this->Psi_T(x) );
}

double Simulator :: Metro_Step(){
  // Genero il candidato punto successivo del rando walk con distribuzione di probabilità:
  // p(x) = 1/(2L); p(y) = 1/(2L); p(z) = 1/(2L)
  double d_x = _rnd.Rannyu(-_L,_L);

  // Propongo il candidato della nuova posizione del random walk
  double x_new = _x + d_x;

  // Computo della probabilità di accettazione
  double A = 0;
  A = min(1.0, Psi_T(x_new)*Psi_T(x_new) / (Psi_T(_x)*Psi_T(_x)) );

  double test = _rnd.Rannyu(0,1);

  //Accettazione
  if(test < A){ 
    this->SetPosition(x_new);
  }
  return A;
}

void Simulator :: Metro_Walk(){

  ofstream out1(_filename1);
  ofstream out2(_filename2);

  out1 << "#" << "\t" << "X" << "\t " <<  endl;
  out2 << "#" << "\t" << "Current_x" << "\t" << "Ave_x" << "\t" << "Stdev_x" << "\t" << "Ave_accept" << "\t" << "Current_E" << "\t" << "Ave_E" << "\t" << "Stdev_E" <<   endl;

  double ave_x = 0;    // Inizializzo la media globale della posizione sui blocchi 
  double ave_x2 = 0;
  
  double ave_A = 0;  // Inizializzo la media globale sui blocchi della probabilità di accettazione

  double ave_E = 0;    // Inizializzo la media globale dell'energia sui blocchi 
  double ave_E2 = 0;

  // Inizio del blocco di equilibrazione, in cui non devo salvare medie globali

  for(int i=0; i < _N1; i++){ // Si lavora sul passo i-esimo del batch di equilibrazione
    double A = this->Metro_Step();
  }
  // Fine del blocco di equilibrazione

  // Inizio dei blocchi di evoluzione post-equilibraizone, dove accumulo statistica attraverso medie globali

  for(int j = 0; j < _M; j++){// Si lavora nel batch j-esimo
    double mean = 0;
    double r  = 0;
    double E = 0; // Valore di aspettazione dell'energia

    for(int i=0; i < _N; i++){ // Si lavora sul passo i-esimo del batch j-esimo
      double A = Metro_Step();

      out1 << j * _N + i << "\t"<< _x << "\t " << endl;

      r = r + _x;
      mean = mean + A;
      E = E + this->E_integrand(_x);

    }
    ave_x = ave_x + (r / _N);
    ave_x2 = ave_x2 + (r / _N)*(r / _N);
    
    ave_A = ave_A + mean / _N;

    ave_E = ave_E + (E / _N);
    ave_E2 = ave_E2 + (E / _N)*(E / _N);

    // Su una sola riga, vengono stampati...
    out2 <<  (j+1) << "\t"                                                        // il numero del blocco
      <<  (r/_N)  << "\t"                                                         // la posizione media sul blocco corrente          
      <<  (ave_x/(j+1))  << "\t"                                                  // la media globale della posizione fino al blocco corrente
      << sqrt(ave_x2 / (j+1) - (ave_x/(j+1))*(ave_x/(j+1))) /sqrt(j+1) << "\t"     // la rispettiva deviazione standard della media
      <<  (ave_A/(j+1)) << "\t"                                                   // la media globale della probabilità di accettazione fino al blocco corrente
      <<  (E / _N)  << "\t"                                                     // il valore di aspettazione dell'energia sul blocco corrente
      <<  (ave_E / (j+1))  << "\t"                                                // la media globale dell'energia fino al blocco corrente 
      <<  sqrt(ave_E2 / (j+1) - (ave_E/(j+1))*(ave_E/(j+1))) /sqrt(j+1) << "\t"    // la rispettiva deviazione standard della media
      << endl;
  }
  // Fine del data-blocking
  out1.close();
  out2.close();
}



