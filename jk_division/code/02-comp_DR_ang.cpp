#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <bitset>

#include "classes/functions.h"
#include "classes/catalogue.h"









int main (void){


  ostringstream sfname;
  string fname;
  ifstream pfr;
  ofstream pfw;


  ////////////////////////////////////////////////////////////////
  //////////          SELECT FIELD: NGC or SGC          //////////
  ////////////////////////////////////////////////////////////////
  string data_input, rand_input, bw_input, output_root, par_output, fib_output;
  int nReg;
  cin>>data_input;
  cin>>rand_input;
  cin>>bw_input;
  cin>>output_root;
  cin>>par_output;
  cin>>fib_output;
  cin>>nReg;

  cout<<"Parameter Loading Complete"<<endl;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  lrg gal;
  sfname<<data_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  gal.set_fname (fname);
  gal.fill_cat_pip();

  int cntg = (int)(gal.id.size());

  int nx = 31;
  int ny = 60;

  vector<long int> bw_ids(cntg);
  matrix<int> bw_weights(cntg, ny);

  load_bw(bw_input, ny, bw_ids, bw_weights);

  cout<<"Catalogue Loading Complete"<<endl;

  /////////////////////////////////////////////////////////////
  //////////          LOAD RANDOM CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rnd;
  sfname<<rand_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  rnd.set_fname (fname);
  rnd.fill_cat_jk();

  int cntr = (int)(rnd.ra.size());

  cout<<"Random Catalogue Loaded"<<endl;


  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double tmin = 0.01;
  double tmax = 3.;
  double st = 0.075;
  tmin = tmin * pow(10., -st/2.);
  int nBint = (int)(1./st*log10(tmax/tmin))+1;
  tmax = tmin * pow(10., nBint*st);
  cout<<"Binning Complete"<<endl;


  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  // Include pairwise weight from Mohammad and Percival 2021 to correct mismatch between auto and cross-pairs
  double jk_weight = 1. - nReg / (2. + sqrt(2.) * (nReg - 1.));

  matrix<double> DR_par(nReg+2, nBint);
  matrix<double> DR_fib(nReg+2, nBint);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++) {
    double iip_weight;
    if (gal.fiber[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
#pragma omp atomic
        sum_fbr+=(int)(__builtin_popcount(bw_weights(i,is)));
      }
      iip_weight = (double)(sum_cov)/(double)(sum_fbr);
    }
    for (int j=0;j<cntr;j++) {
      double pair_sep = ang_sep(gal.ra[i], gal.dec[i], rnd.ra[j], rnd.dec[j]);
      if (pair_sep>tmin && pair_sep<tmax) {
        int bint = (int)(log10(pair_sep/tmin)/st);
        for (int k=0;k<nReg+2;k++) {
          if ((gal.jk[i]!=k && rnd.jk[j]!=k && gal.jk[i]!=-1 && rnd.jk[j]!=-1) || k==(nReg+1)) {
#pragma omp atomic
            DR_par(k, bint) += 1.;
            if (gal.fiber[i]==1) {
#pragma omp atomic
              DR_fib(k, bint) += iip_weight;
            }
          } else if ((gal.jk[i]!=k || rnd.jk[j]!=k) && (gal.jk[i]!=-1 && rnd.jk[j]!=-1)) {
#pragma omp atomic
            DR_par(k, bint) += jk_weight;
            if (gal.fiber[i]==1) {
#pragma omp atomic
              DR_fib(k, bint) += iip_weight * jk_weight;
            }
          }
        }
      }
    }
    if (((i+1)%5000)==0){
      cout<<"I = "<<setfill('0')<<setw(6)<<(i+1)<<"\tNg = "<<cntg;
      print_time();
    }
  }

  cout<<"Pair Counting Complete"<<endl;

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  for (int k=0;k<nReg+2;k++) {
    ostringstream dirname;
    dirname << k;
    sfname<<output_root + dirname.str() + par_output;
    fname = sfname.str();
    cout<<"Par "<<k<<" fname: "<<fname<<endl;
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    for (int i=0;i<nBint;i++) {
      pfw<<tmin*pow(10., ((double)(i) + 0.5)*st)<<"\t"<<setprecision(6)<<fixed<<DR_par(k ,i)<<endl;
    }
    pfw.close();

    sfname<<output_root + dirname.str() + fib_output;
    fname = sfname.str();
    cout<<"Fib "<<k<<" fname: "<<fname<<endl;
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    for (int i=0;i<nBint;i++) {
      pfw<<tmin*pow(10., ((double)(i) + 0.5)*st)<<"\t"<<setprecision(6)<<fixed<<DR_fib(k, i)<<endl;
    }
    pfw.close();
  }
  
  cout<<"Saving Complete"<<endl;

  return (0);
}
