#include "fileio.h"
#include <fstream> // To use ifstream
#include <sstream>
#include <stdio.h>
#include <iostream>

rescueinput filetoinput(string filepath) {
  string var, line, item;
  ifstream f(filepath);

  rescueinput input;

  if (f.good()!=true) std::cout << "File error" << '\n';
  while (getline(f,line)) {
      var=line.substr(line.find('#')+1,line.find('#')+3);
      line=line.substr(0,line.find('#'));
      if (var.find("beta")!=string::npos) input.beta=stod(line);
      if (var.find("phi")!=string::npos) input.phi=stod(line);
      if (var.find("sM")!=string::npos) input.sM=stod(line);
      if (var.find("c")!=string::npos) input.c=stod(line);
      if (var.find("sigma")!=string::npos) input.sigma=stod(line);
      if (var.find("alpha")!=string::npos) input.alpha=stod(line);
      if (var.find("N0")!=string::npos) input.N0=stod(line);
      if (var.find("s0")!=string::npos) input.s0=stod(line);
      if (var.find("tEst")!=string::npos) input.tEst=stod(line);
  }

  return input;
}
