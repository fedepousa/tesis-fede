#include <stdlib.h>
//#include <iostream.h>
#include <iostream>
#include <string>
#include "lib/TDTSPInstance.h"


TDTSPInstance::TDTSPInstance(char* source)
:_source(source)
{ 
  if(_source.fail()){
    cout << "Error en nombre de archivo" << endl;
     exit(0);
  }
  string aux;
  double value;
  int length = strlen(source);
  size_t p=string::npos;
  
  switch(source[length-3]){
     case 't': 
        if(source[length-4] == '.') {
           //Es una instancia .tsp 
           cout << "Entro tsp" << endl;

           while(p==string::npos && !_source.eof()){
              getline(_source,aux);
              p = aux.find("DIMENSION");       
           }
           if(p==string::npos) cout << "Algo funcion� mal...";
           _n = atoi(aux.substr(p+11).c_str());
           p=string::npos;
           while(p==string::npos && !_source.eof()){
              getline(_source,aux);
              p = aux.find("EDGE_WEIGHT_TYPE");       
           }
           if(p==string::npos) cout << "Algo funcion� mal...";
           if(aux.find("EUC_2D")!=string::npos){
              p=string::npos;
              while(p==string::npos && !_source.eof()){
                 getline(_source,aux);
                 p = aux.find("NODE_COORD_SECTION");       
              }
              vector<double> X(_n,0.0);
              vector<double> Y(_n,0.0);

              for (int i = 0; i < _n; i++) _source >> aux >> X[i] >> Y[i];
              _createCube();
              for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                value = floor(sqrt((X[i]-X[j])*(X[i]-X[j])+(Y[i]-Y[j])*(Y[i]-Y[j]))+0.5); 
                for (int k = 0; k < _n; k++) {
                   _costs[i][j][k] = value;
                   _costs[j][i][k] = value;
                }
              }
           } else if (aux.find("GEO")!=string::npos){
              p=string::npos;
              while(p==string::npos && !_source.eof()){
                 getline(_source,aux);
                 p = aux.find("NODE_COORD_SECTION");       
              }
              vector<double> X(_n,0.0);
              vector<double> Y(_n,0.0);
              double Pi=3.141592;
              for (int i = 0; i < _n; i++) {
                 _source >> aux >> X[i] >> Y[i];
                 int deg = (int) X[i];
                 double min = X[i] - deg;
                 X[i] = Pi * (deg + 5.0 * min / 3.0) / 180.0;

                 deg = (int) Y[i];
                 min = Y[i] - deg;
                 Y[i] = Pi * (deg + 5.0 * min / 3.0) / 180.0;
              }
              _createCube();
              double value1, value2, value3, value;
              for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                 value1 = cos(Y[i]-Y[j]);
                 value2 = cos(X[i]-X[j]);
                 value3 = cos(X[i]+X[j]);
                 value = (int)(6378.388 * acos(0.5*((1.0+value1)*value2 - (1.0-value1)*value3))+1.0);
                 for (int k = 0; k < _n; k++) {
                    _costs[i][j][k] = value;
                    _costs[j][i][k] = value;
                 }
              }
           } else if(aux.find("ATT")!=string::npos){
              p=string::npos;
              while(p==string::npos && !_source.eof()){
                 getline(_source,aux);
                 p = aux.find("NODE_COORD_SECTION");       
              }
              vector<double> X(_n,0.0);
              vector<double> Y(_n,0.0);
              for (int i = 0; i < _n; i++) _source >> aux >> X[i] >> Y[i];
              _createCube();
              for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                 value = ceil(sqrt(((X[i]-X[j])*(X[i]-X[j])+(Y[i]-Y[j])*(Y[i]-Y[j]))/10.0)); 
                 for (int k = 0; k < _n; k++) {
                    _costs[i][j][k] = value;
                    _costs[j][i][k] = value;
                 }
              }
           } else if(aux.find("EXPLICIT")!=string::npos){
              //cout<<"EXPLICIT"<<endl;
              p=string::npos;
              while(p==string::npos && !_source.eof()){
                 getline(_source,aux);
                 p = aux.find("EDGE_WEIGHT_FORMAT");       
              }
              if(aux.find("LOWER_DIAG_ROW")!=string::npos){
                 p=string::npos;
                 while(p==string::npos && !_source.eof()){
                    getline(_source,aux);
                    p = aux.find("EDGE_WEIGHT_SECTION");       
                 }
                 _createCube();
                 for (int i = 0; i < _n; i++) for (int j = 0; j <= i; j++) {
                    _source >> value;
                    for (int k = 0; k < _n; k++) {
                       _costs[i][j][k] = value;
                       _costs[j][i][k] = value;
                    }
                 }
              } else if(aux.find("UPPER_ROW")!=string::npos){
                 p=string::npos;
                 while(p==string::npos && !_source.eof()){
                    getline(_source,aux);
                    p = aux.find("EDGE_WEIGHT_SECTION");       
                 }
                 _createCube();
                 for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                    _source >> value;
                    for (int k = 0; k < _n; k++) {
                       _costs[i][j][k] = value;
                       _costs[j][i][k] = value;
                     }
                 }           
              } else if(aux.find("FULL_MATRIX")!=string::npos){
                 p=string::npos;
                 while(p==string::npos && !_source.eof()){
                    getline(_source,aux);
                    p = aux.find("EDGE_WEIGHT_SECTION");       
                 }
                 _createCube();
                 for (int i = 0; i < _n; i++) for (int j = 0; j < _n; j++) {
                    _source >> value;
                    for (int k = 0; k < _n; k++) _costs[i][j][k] = value;
                 }
              }
           } else {
              cout << "�Ehh qu�? " << endl;
           }
        } else {
           //Es una instancia .atsp (TSP asym) 
           cout << "Entro atsp" << endl; 
           getline(_source,aux);
           getline(_source,aux);
           getline(_source,aux);
           _source >> aux;
           _source >> _n;
           getline(_source,aux);
           getline(_source,aux);
           getline(_source,aux);
           getline(_source,aux);
           _createCube();

           for (int i = 0; i < _n; i++) for (int j = 0; j < _n; j++) {
              _source >> value;
              for (int k = 0; k < _n; k++) _costs[i][j][k] = value;
           }
        }
        break;

     case 'd':
        if (source[length - 2] != 'm') {
           //Es una instancia .dat (random)
           cout << "Entro dat" << endl; 
           _source >> _n;
           _createCube();
           for (int i = 0; i < _n; i++) for (int j = 0; j < _n; j++) for (int k = 0; k < _n; k++) _source >> _costs[i][j][k];
        } else {
           if (source[length - 4] == 'a') {
              //Es una instancia .atsp transformada a DMP
              cout << "Entro admp" << endl; 
              string aux;
              getline(_source,aux);
              getline(_source,aux);
              getline(_source,aux);
              _source >> aux;
              _source >> _n;
              _n++;
              getline(_source,aux);
              getline(_source,aux);
              getline(_source,aux);
              getline(_source,aux);
              _createCube();

              for (int i = 1; i < _n; i++) for (int j = 1; j < _n; j++) {
                 _source >> value;
                 for (int k = 1; k < _n - 1; k++) _costs[i][j][k] = (_n - 1 - k)*value;
              }
           } else {
              //Es una instancia .tsp transformada a DMP

              while(p==string::npos && !_source.eof()){
                 getline(_source,aux);
                 p = aux.find("DIMENSION");       
              }
              if(p==string::npos) cout << "Algo funcion� mal...";
              _n = atoi(aux.substr(p+11).c_str());

              p=string::npos;
              while(p==string::npos && !_source.eof()){
                 getline(_source,aux);
                 p = aux.find("EDGE_WEIGHT_TYPE");       
              }
              if(p==string::npos) cout << "Algo funcion� mal...";
              if(aux.find("EUC_2D")!=string::npos){
                 p=string::npos;
                 while(p==string::npos && !_source.eof()){
                    getline(_source,aux);
                    p = aux.find("NODE_COORD_SECTION");       
                 } 
                 vector<double> X(_n,0.0);
                 vector<double> Y(_n,0.0);

                 for (int i = 0; i < _n; i++) _source >> aux >> X[i] >> Y[i];
                 _createCube();
                 for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                    value = floor(sqrt((X[i]-X[j])*(X[i]-X[j])+(Y[i]-Y[j])*(Y[i]-Y[j]))+0.5); 
                    for (int k = 0; k < _n; k++) {
                       _costs[i][j][k] = (_n - k)*value;
                       _costs[j][i][k] = (_n - k)*value;
                    }
                 }
              } else if (aux.find("GEO")!=string::npos){
                 p=string::npos;
                 while(p==string::npos && !_source.eof()){
                    getline(_source,aux);
                    p = aux.find("NODE_COORD_SECTION");       
                 } 
                 vector<double> X(_n,0.0);
                 vector<double> Y(_n,0.0);
                 double Pi=3.141592;
                 for (int i = 0; i < _n; i++) {
                    _source >> aux >> X[i] >> Y[i];
                    cout << i << ": " << X[i] << " - " << Y[i] << endl; 
                    int deg = (int) X[i];
                    double min = X[i] - deg;
                    X[i] = Pi * (deg + 5.0 * min / 3.0) / 180.0;
                    deg = (int) Y[i];
                    min = Y[i] - deg;
                    Y[i] = Pi * (deg + 5.0 * min / 3.0) / 180.0;
                 }
                 _createCube();
                 double value1, value2, value3, value;
                 for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                    value1 = cos(Y[i]-Y[j]);
                    value2 = cos(X[i]-X[j]);
                    value3 = cos(X[i]+X[j]);
                    value = (int)(6378.388 * acos(0.5*((1.0+value1)*value2 - (1.0-value1)*value3))+1.0);
                    for (int k = 0; k < _n; k++) {
                       _costs[i][j][k] = (_n - k)*value;
                       _costs[j][i][k] = (_n - k)*value;
                    }
                 }        
              } else if(aux.find("ATT")!=string::npos){
                 p=string::npos;
                 while(p==string::npos && !_source.eof()){
                    getline(_source,aux);
                    p = aux.find("NODE_COORD_SECTION");       
                 }
                 vector<double> X(_n,0.0);
                 vector<double> Y(_n,0.0);
                 for (int i = 0; i < _n; i++) _source >> aux >> X[i] >> Y[i];
                 _createCube();
                 for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                    value = ceil(sqrt(((X[i]-X[j])*(X[i]-X[j])+(Y[i]-Y[j])*(Y[i]-Y[j]))/10.0)); 
                    for (int k = 0; k < _n; k++) {
                       _costs[i][j][k] = (_n - k)*value;
                       _costs[j][i][k] = (_n - k)*value;
                    }
                 }
              } else if(aux.find("EXPLICIT")!=string::npos){
                 //cout<<"EXPLICIT"<<endl;
                 p=string::npos;
                 while(p==string::npos && !_source.eof()){
                    getline(_source,aux);
                    p = aux.find("EDGE_WEIGHT_FORMAT");       
                 }
                 if(aux.find("LOWER_DIAG_ROW")!=string::npos){
                    p=string::npos;
                    while(p==string::npos && !_source.eof()){
                       getline(_source,aux);
                       p = aux.find("EDGE_WEIGHT_SECTION");       
                    }
                    _createCube();
                    for (int i = 0; i < _n; i++) for (int j = 0; j <= i; j++) {
                       _source >> value;
                       for (int k = 0; k < _n; k++) {
                          _costs[i][j][k] = (_n - k)*value;
                          _costs[j][i][k] = (_n - k)*value;
                       }
                    }
                 } else if(aux.find("UPPER_ROW")!=string::npos){
                    p=string::npos;
                    while(p==string::npos && !_source.eof()){
                       getline(_source,aux);
                       p = aux.find("EDGE_WEIGHT_SECTION");       
                    }
                    _createCube();
                    for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
                       _source >> value;
                       //Modificado: Tiene que llegar hasta _n - 2. (en _n - 1 vale cero) 
                       for (int k = 0; k < _n; k++) {
                          _costs[i][j][k] = (_n - k)*value;
                          _costs[j][i][k] = (_n - k)*value;
                       }
                    }               
                 } else if(aux.find("FULL_MATRIX")!=string::npos){
                    p=string::npos;
                    while(p==string::npos && !_source.eof()){
                       getline(_source,aux);
                       p = aux.find("EDGE_WEIGHT_SECTION");       
                    }
                    _createCube();
                    for (int i = 0; i < _n; i++) for (int j = 0; j < _n; j++) {
                       _source >> value;
                       for (int k = 0; k < _n; k++) _costs[i][j][k] = (_n - k)*value;
                    }
                 }
              } else {
                 cout << "�Ehh qu�? " << endl;
              }
           }
        }      
        break;

     case 'v':
        cout << "Entro VRP" << endl; 
        while(p==string::npos && !_source.eof()){
           getline(_source,aux);
           p = aux.find("DIMENSION");
        }
        if(p==string::npos) cout << "Algo funcion� mal...";
        _n = atoi(aux.substr(p+11).c_str());
        p=string::npos;
        while(p==string::npos && !_source.eof()){
           getline(_source,aux);
           p = aux.find("EDGE_WEIGHT_TYPE");
        }
        if(p==string::npos) cout << "Algo funcion� mal...";

        _createCube();

        if(aux.find("RANA")!=string::npos){
           getline(_source,aux);
           getline(_source,aux);
           for (int i = 0; i < _n; i++) for (int j = 1; j < _n; j++) if(i!=j){   
   	      _source >> value;
              for (int k = 0; k < _n; k++) _costs[i][j][k] = (_n - 1 - k)*value;
           }
           for (int j = 1; j < _n; j++) for (int k = 0; k < _n ; k++) _costs[j][0][k] = 0;
        }

        if(aux.find("RANS")!=string::npos){
           getline(_source,aux);
           getline(_source,aux);
           for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {   
   	      _source >> value;
              for (int k = 0; k < _n; k++) { 
                 _costs[i][j][k] = (_n - 1 - k)*value;
                 _costs[j][i][k] = (_n - 1 - k)*value;
              }
           }
           for (int j = 1; j < _n; j++) for (int k = 0; k < _n ; k++) _costs[j][0][k] = 0;
        }
		
        if(aux.find("EUC_2D")!=string::npos){
           getline(_source,aux);
           getline(_source,aux);
           vector<double> X(_n,0.0);
           vector<double> Y(_n,0.0);
           for (int i = 0; i < _n; i++) _source >> aux >> X[i] >> Y[i];
           for (int i = 0; i < _n; i++) for (int j = i+1; j < _n; j++) {
              value = floor(sqrt( (X[i]-X[j]) * (X[i]-X[j]) + (Y[i]-Y[j]) * (Y[i]-Y[j]) )+0.5 ); 
              for (int k = 0; k < _n; k++) {
                 _costs[i][j][k] = (_n - 1 - k)*value;
                 _costs[j][i][k] = (_n - 1 - k)*value;
              }
           }
           for (int j = 1; j < _n; j++) for (int k = 0; k < _n ; k++) _costs[j][0][k] = 0;
        }
        break;

     case 'T':
        //Es una instancia .TXT (sched)
        cout << "Entro TXT" << endl;
        _source >> _n;
        _n++;
        _createCube();
        getline(_source,aux);
        vector<double> pj(_n);
        for (int i = 1; i < _n; i++) _source >> pj[i];

        double trash;
        for (int i = 1; i < _n; i++) _source >> trash;

        for (int j = 1; j < _n; j++) {
           _source >> value;
           _costs[0][j][0] = _n*value;
        }

        for (int i = 1; i < _n; i++) for (int j = 1; j < _n; j++) {
            _source >> value;
           for (int k = 1; k < _n - 1; k++) _costs[i][j][k] = (_n - k)*(value + pj[j]);
        }
        break;

//     default:
//        break;
  }
}

int TDTSPInstance::getn() {
  return _n;
}

cube& TDTSPInstance::getCosts() {
  return _costs;
}

void TDTSPInstance::_createCube() {
  vector<double> row(_n,0.0);
  vector<vector<double> > matriz(_n);
  cube cubo(_n);
  for (int i = 0; i < _n; i++) {
    matriz[i] = row;
  }

  for (int i = 0; i < _n; i++) {
    cubo[i] = matriz;
  }

  _costs = cubo;
}
