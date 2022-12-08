#include <string>
using namespace std;

typedef struct  {
  double beta;
  double gamma;
  double sP;
  double sM;
  double c;
  double sigma;
  double alpha;
  double u;
  double s0;
} rescueinput;

rescueinput filetoinput(string filepath) ;
