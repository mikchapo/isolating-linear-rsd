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
  string data_input, bw_input, input_root, auw_input, output_root, dd_output;
  int nReg;
  cin>>data_input;
  cin>>bw_input;
  cin>>input_root;
  cin>>auw_input;
  cin>>output_root;
  cin>>dd_output;
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

  int nx = 31;   // Constant value for PIP weighting
  int ny = 60;   // Number of integers used to store the bitwise weights

  vector<long int> bw_ids(cntg);
  matrix<int> bw_weights(cntg, ny);

  load_bw(bw_input, ny, bw_ids, bw_weights);


  cout<<"Catalogue Loading Complete"<<endl;

  ///////////////////////////////////////////////////////////////////////
  //////////          LINEAR BINNING ON pair_sep AND pair_mu          //////////
  ///////////////////////////////////////////////////////////////////////
  double sepmin = 0.1;
  double sepmax = 60.;
  double logSepmin = log10(sepmin);
  double logSepmax = log10(sepmax);
  int nBinsep = 180;
  double dsep = ((sepmax-sepmin)/(double)nBinsep);
  double logDsep = ((logSepmax - logSepmin)/(double)nBinsep);

  double mumin = 0.;
  double mumax = 1.;
  double dmu = 0.01;
  int nBinmu = (int)((mumax-mumin)/dmu);

  double tmin = 0.01;
  double tmax = 3.;
  double st = 0.075;
  tmin = tmin * pow(10., -st/2.);
  int nBint = (int)(1./st*log10(tmax/tmin))+1;
  tmax = tmin * pow(10., nBint*st);

  cout<<"Binning Complete"<<endl;

  matrix<double> auw_matrix(nReg+2, nBint);
  load_auw_matrix(input_root, auw_input, nReg, auw_matrix);

  for (int p=0;p<nReg+2;p++) {
    for (int q=0;q<nBint;q++) {
      cout<<"JK: "<<p<<" theta: "<<q<<" auw: "<<auw_matrix(p, q)<<endl;
    }
  }



  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////
  double xmin = vec_min (gal.x);
  double ymin = vec_min (gal.y);
  double zmin = vec_min (gal.z);

  cout<<"Found coord mins"<<endl;

  double xmax = vec_max (gal.x);
  double ymax = vec_max (gal.y);
  double zmax = vec_max (gal.z);

  cout<<"Found coord maxs"<<endl;

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  cout<<"Found coord ranges"<<endl;

  const int M = 300;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cntg);
  matrix3d<int> label(M, M, M);
  llist (cntg, M, l, gal.x, gal.y, gal.z, lst, label);

  cout<<"Linked List Setup Complete"<<endl;



  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  // Include pairwise weight from Mohammad and Percival 2021 to correct mismatch between auto and cross-pairs
  double jk_weight = 1. - nReg / (2. + sqrt(2.) * (nReg - 1.));

  matrix3d<double> DD(nReg+2,nBinsep,nBinmu);
  int _i;
  // int output_counter = 0;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (_i=0;_i<cntg;_i++){
    if (gal.clustering[_i]==1) {
      int i_1=floor(((gal.x[_i]-xmin)-sepmax)/l);
      int i_2=floor(((gal.x[_i]-xmin)+sepmax)/l);
      i_1 = max (0,i_1);
      i_2 = min (M-1,i_2);
      int j_1=floor(((gal.y[_i]-ymin)-sepmax)/l);
      int j_2=floor(((gal.y[_i]-ymin)+sepmax)/l);
      j_1 = max (0,j_1);
      j_2 = min (M-1,j_2);
      int k_1=floor(((gal.z[_i]-zmin)-sepmax)/l);
      int k_2=floor(((gal.z[_i]-zmin)+sepmax)/l);
      k_1 = max (0,k_1);
      k_2 = min (M-1,k_2);
      for (int il=i_1;il<=i_2;il++){
        for (int jl=j_1;jl<=j_2;jl++){
          for (int kl=k_1;kl<=k_2;kl++){
            int j=label(il,jl,kl);
            while (j!=0){
              if (_i<j && gal.clustering[j]==1) {
                double pair_sep = sep(gal.x[_i], gal.y[_i], gal.z[_i], gal.x[j], gal.y[j], gal.z[j]);
                if (pair_sep>sepmin && pair_sep<sepmax){
                  double r_pi = pi(gal.x[_i], gal.y[_i], gal.z[_i], gal.x[j], gal.y[j], gal.z[j]);
                  double pair_mu = r_pi / pair_sep;

                  if (pair_mu>mumin && pair_mu<mumax){
      	            int binsep = (int)((log10(pair_sep)-logSepmin)/logDsep);
      	            int binmu = (int)((pair_mu-mumin)/dmu);

                    double pair_sep_ang = ang_sep(gal.ra[_i], gal.dec[_i], gal.ra[j], gal.dec[j]);
                    int bint = (int)(log10(pair_sep_ang/tmin)/st);

                    int sum_fbr = 0;
                    int sum_cov = nx*ny;
                    for (int is=0;is<ny;is++) {
#pragma omp atomic
                      sum_fbr+=(int)(__builtin_popcount(bw_weights(_i,is) & bw_weights(j,is)));
                    }
                    double pip_weight = (double)(sum_cov)/(double)(sum_fbr);
                    // if (output_counter<10) {
                    //   cout<<"i: "<<_i<<" j: "<<j<<" ID_i: "<<gal.id[_i]<<" ID_j: "<<gal.id[j]<<" binsep: "<<binsep<<" binmu: "<<binmu<<" pip_weight: "<<pip_weight<<endl;
                    // }
                    for (int k=0;k<nReg+2;k++) {
                      // if (output_counter<10) {
                      //   cout<<"k: "<<k<<" gal.jk[i]: "<<gal.jk[_i]<<" gal.jk[j]: "<<gal.jk[j]<<endl;
                      // }
                      if ((gal.jk[_i]!=k && gal.jk[_i]!=-1 && gal.jk[j]!=k && gal.jk[j]!=-1) || k==(nReg+1)) {
                        double ang_weight = 1.;
                        if (bint < 0) {
                          ang_weight = auw_matrix(k, 0);
                        } else if (bint < nBint) {
                          ang_weight = auw_matrix(k, bint);
                        }
                        // if (output_counter<10) {
                        //  cout<<"k: "<<k<<" ang_weight: "<<ang_weight<<" tot_weight: "<<gal.syst[_i]*gal.noz[_i]*gal.fkp[_i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight*ang_weight<<endl;
                        // }
#pragma omp atomic
                        DD(k,binsep,binmu) = DD(k,binsep,binmu) + gal.syst[_i]*gal.noz[_i]*gal.fkp[_i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight*ang_weight;
                      } else if ((gal.jk[_i]!=k || gal.jk[j]!=k) && (gal.jk[_i]!=-1 && gal.jk[j]!=-1)) {
                        double ang_weight = 1.;
                        if (bint < 0) {
                          ang_weight = auw_matrix(k, 0);
                        } else if (bint < nBint) {
                          ang_weight = auw_matrix(k, bint);
                        }
#pragma omp atomic
                        DD(k,binsep,binmu) = DD(k,binsep,binmu) + gal.syst[_i]*gal.noz[_i]*gal.fkp[_i]*gal.syst[j]*gal.noz[j]*gal.fkp[j]*pip_weight*ang_weight*jk_weight;
                      }
                    }
// pragma omp atomic
                    // output_counter += 1;
                  }
                }
              }
              j = lst[j];
            }
          }
        }
      }
    }
    if ((_i%10000)==0){
      cout<<"I = "<<setfill('0')<<setw(6)<<_i<<"\tNg = "<<cntg;
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
    sfname<<output_root + dirname.str() + dd_output;
    fname = sfname.str();
    sfname.str("");
    sfname.clear();
    pfw.open (fname.c_str(), ios::out);
    for (int i=0;i<nBinsep;i++) {
      for (int j=0;j<nBinmu;j++) {
        pfw<<pow(10., i*logDsep + logSepmin)<<"\t"<<j*dmu<<"\t"<<setprecision(6)<<fixed<<DD(k,i,j)<<endl;
      }
    }
    pfw.close();
  }
  
  cout<<"Saving Complete"<<endl;

  return (0);
}
