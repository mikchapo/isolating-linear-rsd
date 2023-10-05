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
  cout<<"Reading Input"<<endl;

  string data_input, bw_input, output_root, dd_norm_output;
  int nReg;
  cin>>data_input;
  cin>>bw_input;
  cin>>output_root;
  cin>>dd_norm_output;
  cin>>nReg;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////

  cout<<"Loading data cat"<<endl;

  lrg gal;
  sfname<<data_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  gal.set_fname (fname);
  gal.fill_cat_pip();

  cout<<"Loading bw weights"<<endl;

  int cntg = (int)(gal.id.size());

  int nx = 31;
  int ny = 60;


  vector<long int> bw_ids(cntg);
  matrix<int> bw_weights(cntg, ny);
  cout<<"bw_input: "<<bw_input<<endl;

  load_bw(bw_input, ny, bw_ids, bw_weights);



  ///////////////////////////////////////////////////
  //////////         WEIGHT COUNTS         //////////
  ///////////////////////////////////////////////////
  // Include pairwise weight from Mohammad and Percival 2021 to correct mismatch between auto and cross-pairs
  double jk_weight = 1. - nReg / (2. + sqrt(2.) * (nReg - 1.));

  vector<double> DD_norm(nReg+2);
  int i;
  cout<<"Calculating DD norm"<<endl;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg;i++){
    if (gal.clustering[i]==1) {
      for (int j=i+1;j<cntg;j++){
        if (gal.clustering[j]==1) {
          int sum_fbr = 0;
          int sum_cov = nx*ny;
          for (int is=0;is<ny;is++) {
#pragma omp atomic
            sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is) & bw_weights(j, is)));
          }
          double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
          for (int k=0;k<nReg+2;k++) {
            if ((gal.jk[i]!=k && gal.jk[i]!=-1 && gal.jk[j]!=k && gal.jk[j]!=-1) || k==(nReg+1)) {
#pragma omp atomic
              DD_norm[k] += gal.syst[i]*gal.noz[i]*gal.fkp[i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight;
            } else if ((gal.jk[i]!=k || gal.jk[j]!=k) && (gal.jk[i]!=-1 && gal.jk[j]!=-1)) {
#pragma omp atomic
              DD_norm[k] += gal.syst[i]*gal.noz[i]*gal.fkp[i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight*jk_weight;
            }
          }
        }
      }
    }
    if ((i+1)%10000==0) {
      cout<<"Completed "<<i<<" galaxies";
      print_time();
    }
  }
  cout<<"Finished DD ";
  print_time();

  for (int k=0;k<nReg+2;k++) {
    cout<<k<<" DD_norm: "<<DD_norm[k]<<endl;
  }


  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  cout<<"Saving output"<<endl;

  for (int k=0;k<nReg+2;k++) {
    ostringstream dirname;
    dirname << k;
    sfname<<output_root + dirname.str() + dd_norm_output;
    fname = sfname.str();
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    pfw<<"DD_norm"<<endl<<setprecision(6)<<fixed<<DD_norm[k];
    pfw.close();
  }

  return (0);
}
