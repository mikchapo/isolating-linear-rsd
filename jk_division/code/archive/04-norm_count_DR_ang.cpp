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
  string data_input, rand_input, bw_input, output_root, norm_output;
  int nReg;
  cin>>data_input;
  cin>>rand_input;
  cin>>bw_input;
  cin>>output_root;
  cin>>norm_output;
  cin>>nReg;


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
  cout<<"bw_file: "<<bw_input<<endl;

  load_bw(bw_input, ny, bw_ids, bw_weights);

  /////////////////////////////////////////////////////////////
  //////////          LOAD RANDOM CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rnd;
  sfname<<rand_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  rnd.set_fname(fname);
  rnd.fill_cat_jk();

  int cntr = (int)(rnd.ra.size());


  ///////////////////////////////////////////////////
  //////////         WEIGHT COUNTS         //////////
  ///////////////////////////////////////////////////
  // Include pairwise weight from Mohammad and Percival 2021 to correct mismatch between auto and cross-pairs
  double jk_weight = 1. - nReg / (2. + sqrt(2.) * (nReg - 1.));

  vector<double> DR_norm_par(nReg+2);
  vector<double> DR_norm_fib(nReg+2);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
    double iip_weight;
    if (gal.fiber[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
#pragma omp atomic
        sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is)));
      }
      iip_weight = (double)(sum_cov)/(double)(sum_fbr);
    }
    for (int j=0;j<cntr;j++) {
      for (int k=0;k<nReg+2;k++) {
        if ((gal.jk[i]!=k && rnd.jk[j]!=k && gal.jk[i]!=-1 && rnd.jk[j]!=-1) || k==(nReg+1)) {
#pragma omp atomic
          DR_norm_par[k] += 1.;
          if (gal.fiber[i]==1) {
#pragma omp atomic
            DR_norm_fib[k] += iip_weight;
          }
        } else if ((gal.jk[i]!=k || rnd.jk[j]!=k) && (gal.jk[i]!=-1 && rnd.jk[j]!=-1)) {
#pragma omp atomic
          DR_norm_par[k] += jk_weight;
          if (gal.fiber[i]==1) {
#pragma omp atomic
            DR_norm_fib[k] += iip_weight * jk_weight;
          }
        }
      }
    }
  }

  cout<<"Finished DR ";
  print_time();

  for (int k=0;k<nReg+2;k++) {
    cout<<"DR_norm_par "<<k<<": "<<DR_norm_par[k]<<endl<<"DR_norm_fib "<<k<<": "<<DR_norm_fib[k]<<endl;
  }



  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  for (int k=0;k<nReg+2;k++) {
    ostringstream dirname;
    dirname << k;
    sfname<<output_root + dirname.str() + norm_output;
    fname = sfname.str();
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    pfw<<"DR_norm_par\tDR_norm_fib"<<endl<<setprecision(6)<<fixed<<DR_norm_par[k]<<"\t"<<DR_norm_fib[k];
    pfw.close();
  }

  return (0);
}
