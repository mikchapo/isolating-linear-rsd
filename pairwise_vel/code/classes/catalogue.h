#ifndef CATALOGUE_H
#define CATALOGUE_H


#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <cmath>
#include <iomanip>

#include <algorithm>




using namespace std;


class catalogue{
 protected:
  string fname;
 public:
  void set_fname(const string fnamein){fname = fnamein;}
};


class halocat: public catalogue{
 public:
  vector<long int> id;
  vector<int> descid, np, pid;
  vector<double> m200b, vmax, vrms, r200b, rs, x, y, z, vx, vy, vz;
  void fill_cat();
  void fill_cat_v_lin();
  void fill_cat_v_nl();
  // af refers to Aemulus format, the catalogue format made to be compatible with the emular code
  void fill_cat_af();
  void fill_cat_af_v_lin();
  void fill_cat_af_v_nl();
};


class partcat: public catalogue{
 public:
  vector<double> x, y, z, vx, vy, vz;
  void fill_cat();
};


#endif
