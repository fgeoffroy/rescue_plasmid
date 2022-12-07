#include <string>
using namespace std;

typedef struct  {
  double beta;
  double phi;
  double sM;
  double c;
  double sigma;
  double alpha;
  double N0;
  double s0;
  double tEst;
} rescueinput;

rescueinput filetoinput(string filepath) ;
