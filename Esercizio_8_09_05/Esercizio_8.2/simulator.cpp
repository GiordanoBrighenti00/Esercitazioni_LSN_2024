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
Simulator :: Simulator(double mu, double sigma){
  _mu = mu;
  _sigma = sigma;

  // Default values per i passi di campionamento L
  _L1 = 0.6;  // Passo di campionamento di mu per la distribuzione exp^(-beta*<E>) secondo SA
  _L2 = 0.6;  // Passo di campionamento di sigma per la distribuzione exp^(-beta*<E>) secondo SA
  _L3 = 0.25;  // Passo di campionamento di x per la distribuzione |Psi(x)|^2 secondo Metropolis

  // Default value per il valore di aspettazione dell'energia <E>
  _E = 0.0;
}

// Distruttore di default
Simulator :: ~Simulator(){}

// --- INIZIO DEI METODI ---

void Simulator :: SetRandom(Random& rnd){
  _rnd = rnd;
}

void Simulator :: SetStartingPoint(double x_c){
  _x = x_c;
}

void Simulator :: SetPosition(double x){
  _x = x;
}

void Simulator :: SetStep(double L1, double L2, double L3){
  _L1 = L1;
  _L2 = L2;
  _L3 = L3;
}

double Simulator :: GetPosition(){
  return _x;
}

double Simulator :: min(double x, double y){
  if(x < y){
    return x;
  }else{
    return y;
  }
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
  double d_x = _rnd.Rannyu(-_L3,_L3);

  // Propongo il candidato della nuova posizione del random walk
  double x_new = _x + d_x;

  // Computo della probabilità di accettazione
  double A = 0;
  A = this->min(1.0, Psi_T(x_new)*Psi_T(x_new) / (Psi_T(_x)*Psi_T(_x)) );
  
  double test = _rnd.Rannyu(0,1);

  //Accettazione
  if(test < A){ 
    this->SetPosition(x_new);
  }
  return A;
}

double Simulator :: Step_Annealing(double beta, double K){
  double d_mu = _rnd.Rannyu(- _L1, _L1);
  double d_sigma = _rnd.Rannyu(- _L2, _L2);

  // --- OLD: Calcolo dell'energia nella posizione vecchia sullo spazio (mu,sigma)

  double E = 0;

  for (int k = 0; k < K; k++){  // Calcolo il valore di aspettazione dell'energia per (mu, beta) fissati
    double A = Metro_Step();  // Provvede di per sè a proporre ed approvare la nuova mossa nello spazio x, eseguendo Metropolis secondo la distr. |Psi(x)|^2
    E = E + this->E_integrand(_x);
  }
  E = E / K;
  double exp1 = exp(-beta*E);

  // Effettuo forzatamente le mosse per creare la nuova posizione (mu, sigma)
  _mu = _mu + d_mu;
  _sigma = _sigma + d_sigma;

   // --- NEW: Calcolo dell'energia nella potenziale nuova posizione sullo spazio (mu,sigma)

  double E_new = 0;

  for (int k = 0; k < K; k++){  // Calcolo il valore di aspettazione dell'energia per (mu, beta) fissati
    double A = Metro_Step();  // Provvede di per sè a proporre ed approvare la nuova mossa nello spazio x, eseguendo Metropolis secondo la distr. |Psi(x)|^2
    E_new = E_new + this->E_integrand(_x); 
  }
  E_new = E_new / K;
  double exp2 = exp(-beta*E_new);

  if(beta  > 75){
    cout << "WARNING: CRITICAL COOLING --- UNDERFLOW RISK!!!" << "\t" << "beta = " << "\t" << beta  << endl;
  }

  double A = this->min(1.0, exp2 / exp1);

  double test = _rnd.Rannyu(0,1);
  if(test > A){ // Se l'accettazione è rifiutata, si ritorna nella posizione di partenza (mu, sigma), ovvero lo step è statico
    _mu = _mu - d_mu;
    _sigma = _sigma - d_sigma;
  }else{
    E = E_new;
  }

  _E = E;
  return A;
}

void Simulator :: Simulated_Annealing(double T, double alpha, string filename1, string filename2, double N1, double M, double N, double K){
  // Variabili usate solamente all'interno del metodo Simulated_Annealing:
  // T: temperatura iniziale
  // alpha: tasso di raffreddamento del Simulated Annealing
  // N1: numero di passi di equilibrazione, senza accumulo di statistica, con lieve raffreddamento ad ogni step
  // M: numero di blocchi
  // N: numero di passi di Simulated Annealing sullo spazio dei parametri per ciascun blocco, a temperatura fissata: si campionano con Metropolis secondo la distr. epx(-beta*<E>)
  // K: numero di punti su cui si integra il valore di aspettazione dell'energia <E>: si campionano con Metropolis secondo la distr. |Psi(x)|^2

  ofstream out1(filename1);
  ofstream out2(filename2);

  out1 << "#SA" << "\t" << "SA_Energy" << endl;
  out2 << "#Blk" << "\t" << "Blk_Energy" << "\t" << "Global_energy" << "\t" << "Stdev_Energy" << "\t" << "Current_mu" << "\t" << "Current sigma" << "\t" << "Blk_Acceptance" <<   endl;
  
  double beta = 1/T; // beta iniziale

  // --- Inizio del blocco di equlibrazione ---

  for(int i = 0; i < N1; i++){ // passo i-esimo dle blocco di equilibrazione
    double A = this->Step_Annealing(beta, K);
    if(i%50 == 0){
      beta = beta / 0.99;
    }
  }

  // Slot delle medie globali:
  double ave_E = 0;
  double ave_E2 = 0;

  for(int j = 0; j < M; j++){  // Mi muovo sul blocco j-esimo: ogni batch ha una temperatura fissata
    double SA_E = 0; // Energia <E> mediata sugni N passi del blocco j-esimo del simulated annealing
    double SA_A = 0;
    
    for(int i = 0; i < N; i++){  // Mi muovo sul passo i-esimo del blocco j-esimo del Simulated Annealing: ho uno step Montecarlo in cui mi sposto sullo spazio dei parametri (mu,beta)
      double A = this->Step_Annealing(beta,K);

      SA_E = SA_E + _E;
      SA_A = SA_A + A;
      out1 << j * N + i << "\t" << _E << "\t" << endl;
       
    }
    SA_E = SA_E / N;
    SA_A = SA_A / N;
    
    ave_E = ave_E + SA_E;
    ave_E2 = ave_E2 + SA_E*SA_E;
    
    // Sono stampati sulla stessa riga...
    out2 << j+1 << "\t"                   // il numero del blocco 
      << SA_E << "\t"                   // il valore di aspettazione dell'energia <E> mediato sui passi del simulated Annealing del blocco j-.esimo
      << ave_E / (j+1) << "\t"          // la media globale di <E> fino al blocco j-esimo
      << sqrt(ave_E2 / (j+1) - (ave_E/(j+1))*(ave_E/(j+1)) )/(sqrt(j+1)) << "\t"  // la rispettiva incertezza
      << _mu << "\t"                    // il parametro mu al termine del blocco j-esimo del Simulated Annealing
      << _sigma << "\t"                 // il parametro sigma al termine del blocco j-esimo del Simulated Annealing
      << SA_A << "\t"                   // la probabilità di accettazione del blocco j-esimo
      << endl;

    
    beta = beta + (alpha * (1.0 - j/M) ); // Al termine del batch j-esimo, raffreddo il sistema
    _L1 = _L1 * (1 - (j/M)*(j/M) * alpha);
    _L2 = _L2 * (1 - (j/M)*(j/M) * alpha);
  }

  
}



