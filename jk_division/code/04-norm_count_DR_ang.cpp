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
  double alpha = nReg / (2. + sqrt(2.) * (nReg - 1.));

  vector<double> D_eff_par(nReg+2);
  vector<double> D_eff_fib(nReg+2);
  vector<double> R_eff(nReg+2);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
    for (int k=0;k<nReg+2;k++) {
      if ((gal.jk[i]!=k && gal.jk[i]!=-1) || k==(nReg+1)) {
#pragma omp atomic
        D_eff_par[k] += 1.;
      }
    }
    if (gal.fiber[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
#pragma omp atomic
        sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is)));
      }
//    if (sum_fbr==0) {
//      cout<<"Object "<<gal.id[i]<<" sum_fbr=0"<<endl;
//    }
      double iip_weight = (double)(sum_cov)/(double)(sum_fbr);
      for (int k=0;k<nReg+2;k++) {
        if ((gal.jk[i]!=k && gal.jk[i]!=-1) || k==(nReg+1)) {
#pragma omp atomic
          D_eff_fib[k] += iip_weight;
        }
      }
    }
  }

  for (i=0;i<cntr;i++){
    for (int k=0;k<nReg+2;k++) {
      if ((rnd.jk[i]!=k && rnd.jk[i]!=-1) || k==(nReg+1)) {
        R_eff[k] += 1.;
      }
    }
  }

  cout<<"Finished DR ";
  print_time();

  vector<double> DR_norm_par(nReg+2);
  vector<double> DR_norm_fib(nReg+2);

  double D_eff_par_reg, D_eff_fib_reg, R_eff_reg;

  for (int k=0;k<nReg+2;k++) {
    // Calculate the effective number of objects in the removed region
    if (k==(nReg+1)) {
      D_eff_par_reg = 0.;
      D_eff_fib_reg = 0.;
      R_eff_reg = 0.;
    } else {
      D_eff_par_reg = D_eff_par[nReg] - D_eff_par[k];
      D_eff_fib_reg = D_eff_fib[nReg] - D_eff_fib[k];
      R_eff_reg = R_eff[nReg] - R_eff[k];
   }

    // Remove the correct amount of cross pairs from the norm count
    DR_norm_par[k] = D_eff_par[nReg]*R_eff[nReg] - alpha*(D_eff_par_reg*R_eff[nReg] + D_eff_par[nReg]*R_eff_reg - 2*D_eff_par_reg*R_eff_reg) - D_eff_par_reg*R_eff_reg;
    DR_norm_fib[k] = D_eff_fib[nReg]*R_eff[nReg] - alpha*(D_eff_fib_reg*R_eff[nReg] + D_eff_fib[nReg]*R_eff_reg - 2*D_eff_fib_reg*R_eff_reg) - D_eff_fib_reg*R_eff_reg;
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
