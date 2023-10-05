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

  string data_input_ngc, data_input_sgc, bw_input_ngc, bw_input_sgc, output_root, dd_norm_output;
  int nReg;
  cin>>data_input_ngc;
  cin>>data_input_sgc;
  cin>>bw_input_ngc;
  cin>>bw_input_sgc;
  cin>>output_root;
  cin>>dd_norm_output;
  cin>>nReg;

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////

  cout<<"Loading data cat"<<endl;

  lrg gal_ngc;
  sfname<<data_input_ngc;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal_ngc.set_fname (fname);
  gal_ngc.fill_cat_pip();

  lrg gal_sgc;
  sfname<<data_input_sgc;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  gal_sgc.set_fname (fname);
  gal_sgc.fill_cat_pip();
  
  cout<<"Loading bw weights"<<endl;

  int nx = 31;
  int ny = 60;
 
  int cntg_ngc = (int)(gal_ngc.id.size());


  vector<long int> bw_ids_ngc(cntg_ngc);
  matrix<int> bw_weights_ngc(cntg_ngc, ny);
  cout<<"bw_input_ngc: "<<bw_input_ngc<<endl;

  load_bw(bw_input_ngc, ny, bw_ids_ngc, bw_weights_ngc);


  int cntg_sgc = (int)(gal_sgc.id.size());


  vector<long int> bw_ids_sgc(cntg_sgc);
  matrix<int> bw_weights_sgc(cntg_sgc, ny);
  cout<<"bw_input_sgc: "<<bw_input_sgc<<endl;

  load_bw(bw_input_sgc, ny, bw_ids_sgc, bw_weights_sgc);



  ///////////////////////////////////////////////////
  //////////         WEIGHT COUNTS         //////////
  ///////////////////////////////////////////////////
  // Include pairwise weight from Mohammad and Percival 2021 to correct mismatch between auto and cross-pairs
  double jk_weight = 1. - nReg / (2. + sqrt(2.) * (nReg - 1.));

  vector<double> DD_norm(nReg+2);
  int i;

  cout<<"Calculating DD norm"<<endl;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntg_ngc;i++){
    if (gal_ngc.clustering[i]==1) {
      for (int j=0;j<cntg_sgc;j++){
        if (gal_sgc.clustering[j]==1) {
          int sum_fbr = 0;
          int sum_cov = nx*ny;
          for (int is=0;is<ny;is++) {
#pragma omp atomic
            sum_fbr+=(int)(__builtin_popcount(bw_weights_ngc(i, is) & bw_weights_sgc(j, is)));
          }
          double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
          for (int k=0;k<nReg+2;k++) {
            if ((gal_ngc.jk[i]!=k && gal_ngc.jk[i]!=-1 && gal_sgc.jk[j]!=k && gal_sgc.jk[j]!=-1) || k==(nReg+1)) {
#pragma omp atomic
              DD_norm[k] += gal_ngc.syst[i]*gal_ngc.noz[i]*gal_ngc.fkp[i]*gal_sgc.syst[j]*gal_sgc.noz[j]*gal_sgc.fkp[j]*pip_weight;
            } else if ((gal_ngc.jk[i]!=k || gal_sgc.jk[j]!=k) && (gal_ngc.jk[i]!=-1 && gal_sgc.jk[j]!=-1)) {
#pragma omp atomic
              DD_norm[k] += gal_ngc.syst[i]*gal_ngc.noz[i]*gal_ngc.fkp[i]*gal_sgc.syst[j]*gal_sgc.noz[j]*gal_sgc.fkp[j]*pip_weight*jk_weight;
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
    pfw<<"DD_norm"<<endl<<setprecision(6)<<fixed<<DD_norm[k]<<endl;
    pfw.close();
  }
  
  return (0);
}
