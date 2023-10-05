#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <algorithm>

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

  string rand_input, output_root, rr_norm_output;
  int nReg;
  cin>>rand_input;
  cin>>output_root;
  cin>>rr_norm_output;
  cin>>nReg;


  /////////////////////////////////////////////////////////////
  //////////          LOAD RANDOM CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  cout<<"Loading rand cat"<<endl;

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

  vector<double> R_eff_1(nReg+2);
  vector<double> R_eff_2(nReg+2);
  vector<double> RR_norm(nReg+2);
  int i;
  
  cout<<"Calculating RR norm"<<endl;
  for (i=0;i<cntr;i++){
    for (int k=0;k<nReg+2;k++) {
      if ((rnd.jk[i]!=k && rnd.jk[i]!=-1) || k==(nReg+1)) {
        R_eff_1[k] += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
        R_eff_2[k] += rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i]*rnd.syst[i]*rnd.cp[i]*rnd.noz[i]*rnd.fkp[i];
      }
    }
  }

  double R_eff_1_reg, R_eff_2_reg;

  for (int k=0;k<nReg+2;k++) {
    // Calculate the effective number of objects in the remvoed region
    if (k==(nReg+1)) {
      R_eff_1_reg = 0.;
      R_eff_2_reg = 0.;
    } else {
      R_eff_1_reg = R_eff_1[nReg] - R_eff_1[k];
      R_eff_2_reg = R_eff_2[nReg] - R_eff_2[k];
    }
    RR_norm[k] = (R_eff_1[nReg]*R_eff_1[nReg] - R_eff_2[nReg]) / 2. - alpha*(R_eff_1_reg*R_eff_1[k]) - (R_eff_1_reg*R_eff_1_reg - R_eff_2_reg) / 2.;
    cout<<k<<" R_eff_1: "<<R_eff_1[k]<<" R_eff_2: "<<R_eff_2[k]<<" RR_norm: "<<RR_norm[k]<<endl;
  }
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  cout<<"Saving output"<<endl;

  for (int k=0;k<nReg+2;k++) {
    ostringstream dirname;
    dirname << k;
    sfname<<output_root + dirname.str() + rr_norm_output;
    fname = sfname.str();
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    pfw<<"R_eff_1\tR_eff_2\tRR_norm"<<endl<<setprecision(6)<<fixed<<R_eff_1[k]<<"\t"<<R_eff_2[k]<<"\t"<<RR_norm[k];
    pfw.close();
  }
  
  return (0);
}
