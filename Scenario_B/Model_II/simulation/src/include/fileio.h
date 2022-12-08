#include <string>
using namespace std;

typedef struct  {
  double beta;
  double phi0;
  double gamma;
  double sM;
  double c;
  double sigma;
  double alpha;
  double u;
  double N0;
  double s0;
} rescueinput;

rescueinput filetoinput(string filepath) ;
