#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

#include "random.h"
#include "Genetic.h"

using namespace std;

Genetic::Genetic(int N, arma::rowvec x, arma::rowvec y){
  // Costruttore: inserire il numero N di cromosomi che costituiscono una popolazione
  // Inserire le coordinate delle 34 città, da caricare in un vettore di ascisse (x) e un vettore di ordinate (y)
  _N = N;
  _M = 1; // Contatore del numero di newborn generati
  _x_cities = x;
  _y_cities = y;

  // Creazione di una popolazione di deafault, con tutti i croimosomi identici [1, 2, .. , 34]
  // Serve per partire da una popolazione di cromosomi nel caso ci si dimentichi di richiamare
  // Initialize()
  arma::mat pop(_N, 34, arma::fill::zeros);

  // Popolazione nuova: qui vengono salvati i nuovi esemplari prima che sostituiscano la popolazione precedente
  arma::mat newpop(_N, 34, arma::fill::zeros);
  
  for(int i=0; i<_N; i++){
    for(int j=0; j<34; j++){
      pop(i,j) = j+1;
    }
  }
  _pop = pop;
  _newpop = pop;

  // Inizializzo gli indici dei cromosomi selezionati per la riproduzione ad ogni step di evoluzione
  _mate_1 = 0;
  _mate_2 = 1;

  // Inserisco il 100% come probabilità di mutazione dei nuovi individui di default
  _p_m = 1.0;
  // Inserisco il 100% come probabilità di crossover dei nuovi individui di default
  _p_c = 1.0;
}

Genetic :: ~Genetic(){}

void Genetic:: Set_Random(Random& rnd){
  _rnd = rnd;
}

void Genetic :: Set_Mutation_Probability(double prob){
  _p_m = prob;
}

void Genetic :: Set_Crossover_Probability(double prob){
  _p_c = prob;
}

void Genetic :: Initialize(){
  // Inizializzazione della popolazione:generare N cromosomi, ciascuno costituito da una 
  // permutazione dei numeri [1, ... , 34]: attenzione: la città 1 deve rimanere sempre nella posizione 0

  for(int k = 0; k < _N; k++){
    arma::rowvec rcopy = _pop.row(double(k));
    rcopy = Perm_FY(rcopy);
    _pop.row(k) = rcopy;
  }
  
}

arma::rowvec Genetic :: Perm_FY(arma::rowvec vec){ // Questa funzione vuole in input un vettore riga
  // Permutazione pseudo-casuale di un cromosoma di 34 geni, con l'algoritmo di Fischer-Yates:
  // Tutte le possibili permutazioni devono avere la stessa probabilità di essere estratte
  for(int i = 1; i < 34; i++){
    // Estrazione della posizione da scambiare
    double test = _rnd.Rannyu(i,34);
    int pos = double(floor(test));
    std::swap(vec[i],vec[pos]); 
  }
  return vec;
}

//void Genetic :: Swap(double &x, double &y){
 // double z = y;
  //y = x;
 // x = z;
//}

void Genetic :: Get_Specimen(int l){
  // Stampare un cromosoma 
  
  arma::rowvec rcopy = _pop.row(l);
  for(int i = 0; i < 34; i++){
    //cout << i << "\t";
    cout << rcopy(i) << "\t";
  }
  cout << endl;
}

double Genetic :: Euclid_Norm(double x, double y){
  // Norma euclidea di un vettore in 2D, che può consistere nel vettore distanza tra due punti
  double norm = sqrt(x*x + y*y);
  if(norm == 0){
    cout << "Fatal Error: Null Euclidean Norm!" << endl;
  }
  return norm;
}

double Genetic :: Cost_Path(arma::rowvec vec){
  // Calcolare il costo del percorso di un cromosoma: serve per valutare la fitness
  double path = 0;
  for(int i = 0; i < 33; i++){
    // Somma di distanze euclidee tra città adiacenti, percorse nell'ordine riportato dal cromosoma
    int start_city = vec[i] - 1 ; // Salva l'indice vettoriale associato al numero della città di partenza in posizione i-esima nel cromosoma
    int end_city = vec[i+1] - 1; // Salva l'indice vettoriale associato al numero della città di arrivo in posizione i-esima nel cromosoma
    // Spiegazione: vec[i] = [1, 2, ..., 33, 34] , cioè è il numero identificativo della città, riportato nella posizione i-esima del cromosoma
    path = path + this->Euclid_Norm( _x_cities[end_city] - _x_cities[start_city] , _y_cities[end_city] - _y_cities[start_city] );
  }
  
  // Aggiungo al cammino la distanza tra la città di arrivo e la città di partenza, siccome il percorso è chiuso
  int end_city = int(vec[33]) - 1 ;
  int start_city = int(vec[0]) - 1; 
  if(start_city != 0){
    // start_city deve essere sempre la città 1, le cui coordinate sono indicate nell'indice 0 dei vettori cities
    cout << "FATAL ERROR! Check it out!!!" << endl; 
  }
  path = path + this->Euclid_Norm( _x_cities[end_city] - _x_cities[start_city] , _y_cities[end_city] - _y_cities[start_city] );
  return path;
}

double Genetic :: Get_Cost(int l){
  // restituisce il costo del cromosoma l-esimo
  return this->Cost_Path(_pop.row(l));
}

void Genetic :: Sort(){
  // Con un algoritmo di tipo Select_Sort, si riordinano i cromosomi in
  for(int i = 0; i < _N - 1; i++){
    for(int j = i + 1; j < _N; j++){
      arma::rowvec chrom_prev= _pop.row(i);
      arma::rowvec chrom_next = _pop.row(j);

      if(this->Cost_Path(chrom_prev) > this->Cost_Path(chrom_next) ){
        std::swap(chrom_prev, chrom_next);
      }
      _pop.row(i) = chrom_prev;
      _pop.row(j) = chrom_next;
    }

  }
}

void  Genetic :: Selection(){
  // Con un algoritmo di tipo Select_Sort, si riordinano i cromosomi in ordine crescente di Cost_Path,
  // stampando nel fratttempo un vettore di probabilità che indica la probabilità di selezionare i vari cromosomi
  arma::rowvec cost_vec(_N, arma::fill::zeros);
  arma::rowvec prob_vec(_N, arma::fill::zeros);
  double tot_prob = 0;

  for(int i = 0; i < _N - 1; i++){
    for(int j = i + 1; j < _N; j++){
      arma::rowvec chrom_prev= _pop.row(i);
      arma::rowvec chrom_next = _pop.row(j);

      if(this->Cost_Path(chrom_prev) > this->Cost_Path(chrom_next) ){
        std::swap(chrom_prev, chrom_next);
        cost_vec[i] = this->Cost_Path(chrom_prev);
        cost_vec[i+1] = this->Cost_Path(chrom_next);
        //tot_prob = tot_prob + 1 / cost_vec[i];
        //prob_vec[i] = tot_prob;
      }
      _pop.row(i) = chrom_prev;
      _pop.row(j) = chrom_next;
    }
    //tot_prob = tot_prob + 1 / cost_vec[i];
    //prob_vec[i] = tot_prob;

  }

  for(int i = 0; i < _N; i++){
    //cost_vec[i] = this->Cost_Path(_pop.row(i));
    tot_prob = tot_prob + 1 / cost_vec[i];
    prob_vec[i] = tot_prob;
  }
  
  tot_prob = tot_prob + 1 / cost_vec[_N  - 1];
  prob_vec[_N - 1] = tot_prob;

  // Ho creato un vettore di probabilità in cui ogni elemento i-esimo è la somma delle probabilità (non normalizzate)
  // dei cromosomi, ordinati con costo crescente, fino al cromosoma i-esimo
  double test_1 = _rnd.Rannyu(0,tot_prob);
  // Trovare l'elemento più piccolo di prob_vec che sia maggiore di test_1, al fine di identificare il primo cromosoma selezionato
  int specimen_1 = Find_Index(test_1 , prob_vec);

  double test_2 = _rnd.Rannyu(0,tot_prob);
  int specimen_2 = Find_Index(test_2 , prob_vec);

  _mate_1 = specimen_1;
  _mate_2 = specimen_2;

  // Il metodo Selection() ordina tutti gli N cromosomi in ordine crescente di costo e, inoltre, 
  // riporta gli indici dei due cromosomi selezionati per la riproduzione
}

int Genetic :: Find_Index(double threshold, arma::rowvec vec){
  // Identifica il più piccolo elemento di un vettore (ordinato in modo crescente) che supera una certa soglia, riportandone l'indice
  double index = 0;
  for(int i = 0; i < _N; i++){
    if(vec[i] > threshold){
      index = i;
      break;
    }
  }
  return index;
  
}

void Genetic :: Crossover_Reproduction(){
  // Riproduzione con Crossover: vengono generati due nuovi individui che vanno a sostituire 
  // i due cromosomi nella popolazione precedente con minor fitness (cioè l'N-1-esimo e l'N esimo)
  // Innanzitutto, estrarre l'indice in cui avviene il Crossover, che deve essere da 1 a 33
  int crossover_point = floor(_rnd.Rannyu(1,34)); //la città 1 deve essere sempre la prima
  // Nuovi individui
  arma::rowvec newborn_1(34, arma::fill::zeros);
  arma::rowvec newborn_2(34, arma::fill::zeros);
  
  for(int i = 0; i < crossover_point; i++){
    newborn_1[i] = _pop.row(_mate_1)[i];
  }

  for(int i = 0; i < crossover_point; i++){
    newborn_2[i] = _pop.row(_mate_2)[i];
  }

  // L'obiettivo è ricostruire la seconda parte del cromosoma mate_1 (che chiamo slab), 
  // riportando le città mancanti nell'ordine in cui sono riportate nel cromosoma mate_2
  int size = 34 - crossover_point;
  // Vettore di salvataggio dello slab del newborn_1:
  arma::rowvec slab_1(size, arma::fill::zeros);
  //Vettore di salvataggio dello slab del newborn_2:
  arma::rowvec slab_2(size, arma::fill::zeros);
  
  for(int i = 0; i < size; i++){
    slab_1[i] = _pop.row(_mate_1)[crossover_point + i];
  }
  int k = 0;
  for(int i = 0; i < 34; i++){
    
    for(int j = 0; j < size; j++){
      if(_pop.row(_mate_2)[i] == _pop.row(_mate_1)[j + crossover_point] ){
        slab_1[k] = _pop.row(_mate_2)[i];
        k ++;
      }
    }
  }

  // L'obiettivo è ricostruire la seconda parte del cromosoma mate_2 (che chiamo slab), 
  // riportando le città mancanti nell'ordine in cui sono riportate nel cromosoma mate_1
  for(int i = 0; i < size; i++){
    slab_2[i] = _pop.row(_mate_2)[crossover_point + i];
  }
  k = 0;
  for(int i = 0; i < 34; i++){
    
    for(int j = 0; j < size; j++){
      if(_pop.row(_mate_1)[i] == _pop.row(_mate_2)[j + crossover_point] ){
        slab_2[k] = _pop.row(_mate_1)[i];
        k ++;
      }
    }
  }

  for(int i = 0; i < size; i++){
    newborn_1[i + crossover_point] = slab_1[i];
  }

  // Completato il newborn_1

  for(int i = 0; i < size; i++){
    newborn_2[i + crossover_point] = slab_2[i];
  }

  // Completato il newborn 2

  //check
  for(int i = 0; i < size - 1; i++){
    for(int j = i + 1; j < size; j++){
      if(slab_1[i] == slab_1[j] ){
        cout << "FATAL ERROR! Incorrect slab!!!" << endl;
      }
      if(slab_2[i] == slab_2[j]){
        cout << "FATAL ERROR! Incorrect slab!!!" << endl;
      }
    }
  }

  double test_crossover = _rnd.Rannyu(0,1);
 
  if(test_crossover < _p_c){
    _newpop.row(_M) = newborn_1;
    if(_N > _M + 1 ){// Impedisco di riempire uno slot che non c'è
      _newpop.row(_M + 1) = newborn_2;
    }
  }
  
  
}

void Genetic :: Mutation(){
  double test_1 = _rnd.Rannyu(0,1);
  if(test_1 < _p_m){
    this->Mutation_1(_M);
  }

  double test_2 = _rnd.Rannyu(0,1);
  if(test_2 < _p_m){// Impedisco di riempire uno slot che non c'è
    if(_N > _M + 1){
      this->Mutation_1(_M + 1);
    }
  }

  double test_3 = _rnd.Rannyu(0,1);
  if(test_3 < _p_m){
    this->Mutation_2(_M);
    
  }

  double test_4 = _rnd.Rannyu(0,1);
  if(test_4 < _p_m){
    if(_N > _M + 1 ){// Impedisco di riempire uno slot che non c'è
      this->Mutation_2(_M + 1);
    }
  }

  double test_5 = _rnd.Rannyu(0,1);
  if(test_5 < _p_m){
    this->Mutation_3(_M);
  }

  double test_6 = _rnd.Rannyu(0,1);
  if(test_6 < _p_m){
    if(_N > _M + 1 ){// Impedisco di riempire uno slot che non c'è
      this->Mutation_3(_M + 1);
    }
  }
}

void Genetic :: Mutation_1(int index){
  // Mutation: per ogni cromosoma, estrarre due numeri casuali tra 0 e 34 per determinare la posizione
  int mutation_point_1 = floor(_rnd.Rannyu(1,34));
  int mutation_point_2 = floor(_rnd.Rannyu(1,34));
  while(mutation_point_1 == mutation_point_2){ 
    // Ripeto l'estrazione del secondo punto di mutazione ogni volta che capita che sia uguale al primo
    mutation_point_2 = floor(_rnd.Rannyu(1,34));
  }

  std::swap( _newpop.row(index)[mutation_point_1] , _newpop.row(index)[mutation_point_2] );

  int mutation_point_3 = floor(_rnd.Rannyu(1,34));
  int mutation_point_4 = floor(_rnd.Rannyu(1,34));
  while(mutation_point_3 == mutation_point_4){ 
    // Ripeto l'estrazione del secondo punto di mutazione ogni volta che capita che sia uguale al primo
    mutation_point_4 = floor(_rnd.Rannyu(1,34));
  }
}

void Genetic :: Mutation_2(int index){
  // Definisco la sequenza da traslare
  int start_point = floor(_rnd.Rannyu(1, 33)); // Pto iniziale in pos. da 1 a 31
  int end_point = floor(_rnd.Rannyu(start_point, 33)); // Pto finale in pos. da start a 32
  int slab_size = 1 + end_point - start_point; // va da 1 a 32
  // Definisco la lunghezza di traslazione
  int shift_size = floor(_rnd.Rannyu(1, 33 - end_point)); // Va da 1 a 32

  // Ricopio la sequenza del vettore di partenza da riportare a SX

  // Nuovo cromosoma
  arma::rowvec new_chrom(34, arma::fill::zeros);
  new_chrom = _newpop.row(index);

  // Inserisco la sequenza traslata
  for(int i = 0; i < slab_size; i++){
    new_chrom[start_point + shift_size + i] = _newpop.row(index)[start_point + i];
  }

  // Inserisco la sequenza riportata a SX
  for(int i=0; i < shift_size; i++){
    new_chrom[start_point + i] = _newpop.row(index)[end_point + 1 + i];
  }
  
  _newpop.row(index) = new_chrom;
}

void Genetic :: Mutation_3(int index){
  // Inversione speculare di una sequenza di città all'interno del cromosoma
  int start_point = floor(_rnd.Rannyu(1, 34));
  int end_point = floor(_rnd.Rannyu(start_point, 34));
  int length = end_point - start_point + 1;
  // Definisco il nuovo cromosoma
  arma::rowvec new_chrom(34, arma::fill::zeros);
  new_chrom = _newpop.row(index);
  // Inverto la sequenza di città
  for(int i = 0; i < length; i++){
    new_chrom[i + start_point] = _newpop.row(index)[end_point - i];
  }
  _newpop.row(index) = new_chrom;
}

void Genetic :: Evolution_Step(){
  // Step evolutivo: a partire da una popolazione di cromosomi, vengono effettuati i seguenti processi:
  _M = 0;
  while(_M < _N){
    // Selezione: ordinamento per fitness e selezione dei cromosomi da far riprodurre
    this->Selection();
    // Crossover: costruisco i due nuovi individui con il crossover e li metto al posto dei due cromosomi
    // della popolazione precedente meno performanti
    this->Crossover_Reproduction();
    // Mutation: effettuo la mutazione dei due nuovi individui, i quali sono ancora collocati negli 
    // ultimi due indici delle popolazione
    this->Mutation();
    _M = _M + 2;
    
  }
  _pop = _newpop;

    
}

void Genetic :: Print_Population(){
  // Sort: ordinare tutta la nuova popolazione, inclusi i nuovi membri, nell'ordine crescente per costo
  this->Sort();
  cout << endl;
  cout << "#Chromosome" << "\t" << "Cost" << endl;
  for(int i=0; i<_N; i++){
    cout << i << "\t" << this->Get_Cost(i) << endl;
  }
}

double Genetic :: Get_City(int i, int j){
  return _pop.row(i)[j];
}

