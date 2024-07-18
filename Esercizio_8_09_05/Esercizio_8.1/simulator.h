#ifndef __Simulator__
#define __Simulator__

#include <string>

using namespace std;

class Simulator{

private:
  double _mu, _sigma, _x;  // Parametri _mu e _sigma
  double _L;  // Passo di campionamento del Metropolis
  int _N1, _N, _M;  // N1: passi del batch di equilibrazione, N: passi del batch; M: # di batch
  string _filename1, _filename2; //filename1: file di salvataggio delle posizioni istantanee , filename2: file di salvatggio della statistica con il data-blocking
  Random _rnd;
public:
  //Costruttori
  Simulator();
  Simulator(double mu, double sigma, int N1, int N, int M, double L);
  //Distruttori
  ~Simulator();
  //Metodi
  void SetOutputFiles(string filename1, string filename2);
  void SetRandom(Random& rnd);
  void SetStartingPoint(double x_c);
  void SetPosition(double x);
  double GetPosition();
  double Psi_T(double x);
  //double GetExpectEValue();
  double Metro_Step();
  void Metro_Walk();
  double V_Potential(double x);
  double Lapl_T(double x);
  double E_integrand(double x);


};

#endif // __Simulator__