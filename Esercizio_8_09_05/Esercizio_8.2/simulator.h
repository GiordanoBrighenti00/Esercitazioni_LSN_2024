#ifndef __Simulator__
#define __Simulator__

#include <string>

using namespace std;

class Simulator{

private:
  double _mu, _sigma, _x;  // Parametri _mu e _sigma
  double _L1, _L2, _L3;  // Passi di integrazione per il Metropolis
  double _E;  // Valore di aspettazione dell'energia <E>
  Random _rnd;
public:
  //Costruttori
  Simulator();
  Simulator(double mu, double sigma);
  //Distruttori
  ~Simulator();
  //Metodi
  void SetRandom(Random& rnd);
  void SetStartingPoint(double x_c);
  void SetPosition(double x);
  void SetStep(double L1, double L2, double L3);
  double GetPosition();
  double Psi_T(double x);
  //double GetExpectEValue();
  double Metro_Step();
  void Metro_Walk(double M, double N, double N1, string filename1, string filename2);
  double Step_Annealing(double beta, double K);
  void Simulated_Annealing(double T, double alpha, string filename1, string filename2, double N1, double M, double N, double K);
  double V_Potential(double x);
  double min(double x, double y);
  double Lapl_T(double x);
  double E_integrand(double x);


};

#endif // __Simulator__