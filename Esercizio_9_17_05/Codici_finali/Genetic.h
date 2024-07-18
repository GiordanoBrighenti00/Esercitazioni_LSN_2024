#ifndef __Simulator__
#define __Simulator__

#include <string>
#include <armadillo>

using namespace std;

class Genetic{

  private:
    // Posizioni delle città: non vengono modificate durante l'esecuzione
    arma::rowvec _x_cities;
    arma::rowvec _y_cities;
    // Popolazione di individui: ciascuno è un vettore di dimensione N
    int _N;
    arma::mat _pop;
    arma::mat _newpop;
    // Indici dei cromosomi da far accoppiare per ogni step evolutivo
    int _mate_1;
    int _mate_2;
    // Numero di generazioni M: quante volte devo creare coppie di nuovi individui
    int _M;
    // Generatore di numeri pseudo-casuali
    Random _rnd;
    // Probabilità di mutazione di un nuovo individuo
    double _p_m;
    // Probabilità di crossover durante la riproduzione
    double _p_c;

  protected:

  public:
    // Costruttore di default: riceve le dimensioni della popolazione,
    // le coordinate fissate delle città e imposta di default tutti i cromosomi con 
    // la permutazione banale [1, 2, ... , 34], mette M = 1
    Genetic(int N, arma::rowvec x , arma::rowvec y); 
    // Distruttore di default
    ~Genetic();
    // Inserire il generatore di numeri pseudocasuali
    void Set_Random(Random& rnd);
    // Inserire la probabilità di mutazione per ciascun nuovo individuo
    void Set_Mutation_Probability(double prob);
    // Inserire la probabilità di crossover per ogni riproduzione
    void Set_Crossover_Probability(double prob);
    // Method to save the seed to a file
    void Initialize();
    // Permutatore Fischer-Yates
    arma::rowvec Perm_FY(arma::rowvec vec);
    // Operatore di scambio di due variabili (è utilizzato per scambiare due elementi in un vec)
    //void Swap(double& x, double& y);
    // Accesso al neurone l-esimo della popolazione
    void Get_Specimen(int l);
    // Funzione per calcolare la norma Euclidea
    double Euclid_Norm(double x, double y);
    // Funzione che calcola la fitness di un cromosoma
    double Cost_Path(arma::rowvec vec);
    // Funione che riporta il costo-path del cromosoma l-esimo
    double Get_Cost(int l);
    // Funzione di ordinamento della popolazione in base al costo-path
    void Sort();
    // Funzione di selezione di due cromosomi 
    void Selection();
    // Funzione che riporta l'indice del cromosoma selezionato secondo la probabilità data dal fitness
    int Find_Index(double threshold, arma::rowvec vec);
    // Operatore di Riproduzione con Crossover: crea due nuovi individudi
    void Crossover_Reproduction();
    // Operatore di Mutazione sui nuovi individui
    void Mutation();
    void Mutation_1(int index);
    void Mutation_2(int index);
    void Mutation_3(int index);
    // Step di evoluzione: effettua le operazioni di Selezione, Crossover e Mutazione
    void Evolution_Step();
    // Stampare a video la popolazione: riportare una tabella con numero del cromosoma e il suo cost_path
    void Print_Population();
    // Funzione del tutto generale: restituisce l'indice identificativo della città in
    // posizione j-esima all'interno del cromosoma i-esimo
    double Get_City(int i, int j);


};


#endif