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
  string data_input, rand_input, bw_input, output_root, dr_norm_output;
  int nReg;
  cin>>data_input;
  cin>>rand_input;
  cin>>bw_input;
  cin>>output_root;
  cin>>dr_norm_output;
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
  cout<<"bw_input: "<<bw_input<<endl;

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

  vector<double> D_eff(nReg+2);
  vector<double> R_eff(nReg+2);
  vector<double> DR_norm(nReg+2);
  int i;   
  for (i=0;i<cntg;i++){
    if (gal.clustering[i]==1) {
      int sum_fbr = 0;
      int sum_cov = nx*ny;
      for (int is=0;is<ny;is++) {
        sum_fbr+=(int)(__builtin_popcount(bw_weights(i, is)));
      }
      double iip_weight = (double)(sum_cov)/(double)(sum_fbr);
      for (int k=0;k<nReg+2;k++) {
        if ((gal.jk[i]!=k && gal.jk[i]!=-1) || k==(nReg+1)) {
          D_eff[k] += gal.syst[i]*gal.noz[i]*gal.fkp[i]*iip_weight;
        }
      }
    }
  }
  for (i=0;i<cntr;i++){
    for (int k=0;k<nReg+2;k++) {
      if ((rnd.jk[i]!=k && rnd.jk[i]!=-1) || k==(nReg+1)) {
        R_eff[k] += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
      }
    }
  }
  cout<<"Finished DR ";
  print_time();

  double D_eff_reg, R_eff_reg;

  for (int k=0;k<nReg+2;k++) {
    // Calculate the effective number of objects in the remvoed region
    if (k==(nReg+1)) {
      D_eff_reg = 0.;
      R_eff_reg = 0.;
    } else {
      D_eff_reg = D_eff[nReg] - D_eff[k];
      R_eff_reg = R_eff[nReg] - R_eff[k];
    }

    // Remove the correct amount of cross pairs from the norm count
    DR_norm[k] = D_eff[nReg]*R_eff[nReg] - alpha*(D_eff_reg*R_eff[nReg] + D_eff[nReg]*R_eff_reg - 2*D_eff_reg*R_eff_reg) - D_eff_reg*R_eff_reg;
    cout<<k<<" D_eff: "<<D_eff[k]<<" R_eff: "<<R_eff[k]<<" DR_norm: "<<DR_norm[k]<<endl;
  }
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  for (int k=0;k<nReg+2;k++) {
    ostringstream dirname;
    dirname << k;
    sfname<<output_root + dirname.str() + dr_norm_output;
    fname = sfname.str();
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    pfw<<"D_eff\tR_eff\tDR_norm"<<endl<<setprecision(6)<<fixed<<D_eff[k]<<"\t"<<R_eff[k]<<"\t"<<DR_norm[k];
    pfw.close();
  }
  
  return (0);
}
