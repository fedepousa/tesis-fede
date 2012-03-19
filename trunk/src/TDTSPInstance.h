#ifndef TDTSPINSTANCE_HH
#define TDTSPINSTANCE_HH

#include <fstream>
#include <vector>
#include <string>

#include <math.h>
#include <stdio.h>
#include <string.h>
using namespace std;

#include "defs.h"

#define cube vector<vector<vector <double>>>

class TDTSPInstance {

  public:
    TDTSPInstance(char* source);
    int getn();
    cube &getCosts();

  private:
    void _createCube();

    ifstream _source;
    int _n;
    cube _costs;

};
#endif
