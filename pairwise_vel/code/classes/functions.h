#ifndef FUNCTIONS_H
#define FUNCTIONS_H


#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>


#ifndef PI
#define PI 3.14159265
#endif

using namespace std;

typedef vector< vector<double> > matrix2d;


template <class T>
T vec_min (const vector<T> &vec){
  T ris = 1.e+100;
  for (int i_n=0;i_n<(int)(vec.size());i_n++){
    if (vec[i_n]<=ris)
      ris = vec[i_n];
  }
  return(ris);
}


template <class T>
T vec_max (const vector<T> &vec){
  T ris = -1.e+100;
  for (int i_n=0;i_n<(int)(vec.size());i_n++){
    if (vec[i_n]>=ris)
      ris = vec[i_n];
  }
  return(ris);
}


/////////////////////////////////////////////////////////////////////////////////////////
/////////          MATRICES          ////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
template <class T>
class matrix{
private:
  int nr,nc;
  vector<vector<T> > vecvec;
public:
  matrix(const int &nr, const int &nc){
    vecvec =  vector<vector<T> >(nr, vector<T>(nc,0));
  };
  T &operator ()(const int &i, const int &j){return vecvec[i][j];};
  T operator ()(const int &i, const int &j) const {return vecvec[i][j];};
};


template <class T>
class matrix3d{
private:
  int nr,nc,nd;
  vector<vector<vector<T> > > vecvecvec;
public:
  matrix3d(const int &nr, const int &nc, const int &nd){
    vecvecvec =  vector<vector<vector<T> > >(nr, vector<vector<T> >(nc, vector<T>(nd,0)));
  };
  T &operator ()(const int &i, const int &j, const int &k){return vecvecvec[i][j][k];};
  T operator ()(const int &i, const int &j, const int &k) const {return vecvecvec[i][j][k];};
};


//////////////////////////////////////////////////////////////////////////////////////////
//////////          LINKED LIST          /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void llist (const int &, const int &, const double &, const vector<double> &, const vector<double> &, const vector<double> &, vector<int> &, matrix3d<int> &);


double sep(const double &, const double &, const double &, const double &, const double &, const double &);
double vel(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &);
void print_time();
void get_time(int &, int &, int &, int &, int &, int &);

#endif
