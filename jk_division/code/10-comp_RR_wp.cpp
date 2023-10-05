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
  string rand_input, output_root, rr_output;
  int nReg;
  cin>>rand_input;
  cin>>output_root;
  cin>>rr_output;
  cin>>nReg;
  

  /////////////////////////////////////////////////////////////
  //////////          LOAD GALAXY CATALOGUE          //////////
  /////////////////////////////////////////////////////////////
  randcat rand;
  sfname<<rand_input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
    
  rand.set_fname (fname);
  rand.fill_cat_jk();
    
  int cntr = (int)(rand.ra.size());


  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double r_perp_min = 0.1;
  double r_perp_max = 60.;
  double log_r_perp_min = log10(r_perp_min);
  double log_r_perp_max = log10(r_perp_max);
  int n_bin_perp = 180;
  double dr_perp = ((r_perp_max-r_perp_min)/(double)n_bin_perp);
  double log_dr_perp = ((log_r_perp_max - log_r_perp_min)/(double)n_bin_perp);
 
  double r_pi_min = 0.;
  double r_pi_max = 150.;
  int n_bin_pi = 150;
  double dr_pi = (r_pi_max - r_pi_min)/ (double)n_bin_pi;



 

  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////
  double xmin = vec_min (rand.x);
  double ymin = vec_min (rand.y);
  double zmin = vec_min (rand.z);

  double xmax = vec_max (rand.x);
  double ymax = vec_max (rand.y);
  double zmax = vec_max (rand.z);

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  const int M = 300;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cntr);
  matrix3d<int> label(M, M, M);
  llist (cntr, M, l, rand.x, rand.y, rand.z, lst, label);
    
  double rmax = sqrt(r_perp_max*r_perp_max + r_pi_max*r_pi_max)+10.;
   

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  // Include pairwise weight from Mohammad and Percival 2021 to correct mismatch between auto and cross-pairs
  double jk_weight = 1. - nReg / (2. + sqrt(2.) * (nReg - 1.));

  matrix3d<double> RR(nReg+2,n_bin_perp, n_bin_pi);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntr;i++){
    int i_1=floor(((rand.x[i]-xmin)-rmax)/l);
    int i_2=floor(((rand.x[i]-xmin)+rmax)/l);
    i_1 = max (0,i_1);
    i_2 = min (M-1,i_2);
    int j_1=floor(((rand.y[i]-ymin)-rmax)/l);
    int j_2=floor(((rand.y[i]-ymin)+rmax)/l);
    j_1 = max (0,j_1);
    j_2 = min (M-1,j_2);
    int k_1=floor(((rand.z[i]-zmin)-rmax)/l);
    int k_2=floor(((rand.z[i]-zmin)+rmax)/l);
    k_1 = max (0,k_1);
    k_2 = min (M-1,k_2);
    for (int il=i_1;il<=i_2;il++){
      for (int jl=j_1;jl<=j_2;jl++){
        for (int kl=k_1;kl<=k_2;kl++){
          int j=label(il,jl,kl);
          while (j!=0){
            if (i<j){
              double r_pi = pi(rand.x[i], rand.y[i], rand.z[i], rand.x[j], rand.y[j], rand.z[j]);
              if (r_pi>r_pi_min && r_pi<r_pi_max){
      	        double r_perp = rp(rand.x[i], rand.y[i], rand.z[i], rand.x[j], rand.y[j], rand.z[j], r_pi);
      		  
      	        if (r_perp>r_perp_min && r_perp<r_perp_max){
      	          int bin_perp = (int)((log10(r_perp)-log_r_perp_min)/log_dr_perp);
      	          int bin_pi = (int)((r_pi - r_pi_min)/dr_pi);
                  for (int k=0;k<nReg+2;k++) {
                    if ((rand.jk[i]!=k && rand.jk[i]!=-1 && rand.jk[j]!=k && rand.jk[j]!=-1) || k==(nReg+1)) {
      #pragma omp atomic
      	              RR(k, bin_perp, bin_pi) = RR(k, bin_perp, bin_pi) + rand.syst[i]*rand.cp[i]*rand.noz[i]*rand.fkp[i]*rand.syst[j]*rand.cp[j]*rand.noz[j]*rand.fkp[j];
                    } else if ((rand.jk[i]!=k || rand.jk[j]!=k) && (rand.jk[i]!=-1 && rand.jk[j]!=-1)) {
      #pragma omp atomic
      	              RR(k, bin_perp, bin_pi) = RR(k, bin_perp, bin_pi) + rand.syst[i]*rand.cp[i]*rand.noz[i]*rand.fkp[i]*rand.syst[j]*rand.cp[j]*rand.noz[j]*rand.fkp[j]*jk_weight;
                    }
                  }
      	        }
              }
      	    }
      	    j = lst[j];
      	  }
        }
      }
    }
    if ((i%10000)==0){
      cout<<"i = "<<setfill('0')<<setw(6)<<i<<"\tNr = "<<cntr;
      print_time();
    }
  }
    
    

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  for (int k=0;k<nReg+2;k++) {
    ostringstream dirname;
    dirname << k;
    sfname<<output_root + dirname.str() + rr_output;
    fname = sfname.str();
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    for (int i=0;i<n_bin_perp;i++)
      for (int j=0;j<n_bin_pi;j++)
        pfw<<pow(10., i*log_dr_perp + log_r_perp_min)<<"\t"<<j*dr_pi<<"\t"<<setprecision(6)<<fixed<<RR(k,i,j)<<endl;
    pfw.close();
  }
  
  return (0);
}
