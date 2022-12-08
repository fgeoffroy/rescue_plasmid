#include <string>
using namespace std;

typedef struct  {
  double beta;
  double phi;
  double sM;
  double c;
  double alpha;
  double u;
  double N0;
  double s0;
} rescueinput;

rescueinput filetoinput(string filepath) ;
