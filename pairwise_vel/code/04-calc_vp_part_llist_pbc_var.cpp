// Calculate the mean pairwise velocity of particles in separation bins
// v0.1.0, 2021-09-14 - Code started

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


  string input, output;
//  int nBinsep;
  cin>>input;
  cin>>output;
//  cin>>nBinsep;

  partcat part;
  sfname<<input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  part.set_fname (fname);
  part.fill_cat();

  int cntp = (int)(part.x.size());
  cout<<"N_part: "<<cntp<<endl;

  ///////////////////////////////////////////////////////////////////////
  //////////          LOGARITHMIC BINNING ON pair_sep               //////////
  ///////////////////////////////////////////////////////////////////////
  double sepmin = 0.01;
  double sepmax = 100.;
  double logSepmin = log10(sepmin);
  double logSepmax = log10(sepmax);
  int nBinsep = 80;
//  double dsep = ((sepmax-sepmin)/(double)nBinsep);
  double logDsep = ((logSepmax - logSepmin)/(double)nBinsep);


  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////
  double xmin = vec_min (part.x);
  double ymin = vec_min (part.y);
  double zmin = vec_min (part.z);

  double xmax = vec_max (part.x);
  double ymax = vec_max (part.y);
  double zmax = vec_max (part.z);

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  const double Lbox = 1050.;
  const int M = 100;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cntp);
  matrix3d<int> label(M, M, M);
  llist (cntp, M, l, part.x, part.y, part.z, lst, label);


  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  vector<double> vp(nBinsep);
  vector<double> vp2(nBinsep);
  vector<double> Nvp(nBinsep);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cntp;i++){
    int i_1=floor(((part.x[i]-xmin)-sepmax)/l);
    int i_2=floor(((part.x[i]-xmin)+sepmax)/l);

    int j_1=floor(((part.y[i]-ymin)-sepmax)/l);
    int j_2=floor(((part.y[i]-ymin)+sepmax)/l);

    int k_1=floor(((part.z[i]-zmin)-sepmax)/l);
    int k_2=floor(((part.z[i]-zmin)+sepmax)/l);
    // Loop over all boxes within sepmax of the target
    int im, jm, km;
    double  xoff, yoff, zoff;
    for (int il=i_1;il<=i_2;il++){
      if (il < 0) {
        im = M + il;
        xoff = -Lbox;
      } else if (il > M-1) {
        im = il - M;
        xoff = Lbox;
      } else {
        im = il;
        xoff = 0.;
      }
      for (int jl=j_1;jl<=j_2;jl++){
        if (jl < 0) {
          jm = M + jl;
          yoff = -Lbox;
        } else if (jl > M-1) {
          jm = jl - M;
          yoff = Lbox;
        } else {
          jm = jl;
          yoff = 0.;
        }
        for (int kl=k_1;kl<=k_2;kl++){
          if (kl < 0) {
            km = M + kl;
            zoff = -Lbox;
          } else if (kl > M-1) {
            km = kl - M;
            zoff = Lbox;
          } else {
            km = kl;
            zoff = 0.;
          }
          // Start j as the index of the most recent object to be placed in box im, jm, km
          int j=label(im,jm,km);
          // If j==0 then all objects in the box have been looped over and the loop moves to the next box
          while (j!=0){
            if (i<j){
              double pair_sep = sep(part.x[i], part.y[i], part.z[i], part.x[j] + xoff, part.y[j] + yoff, part.z[j] + zoff);
              if (pair_sep>sepmin && pair_sep<sepmax){
                double pair_vel = vel(pair_sep, part.x[i], part.y[i], part.z[i], part.x[j] + xoff, part.y[j] + yoff, part.z[j] + zoff, part.vx[i], part.vy[i], part.vz[i], part.vx[j], part.vy[j], part.vz[j]);
                int binsep = (int)((log10(pair_sep)-logSepmin)/logDsep);
                if (binsep == nBinsep) {
                  binsep = nBinsep - 1;
                }
#pragma omp atomic
                Nvp[binsep] += 1.;
#pragma omp atomic
	        vp[binsep] += pair_vel;
#pragma omp atomic
                vp2[binsep] += pair_vel*pair_vel;
              }
            }
            // Change j to the index of the the previous object placed in box im, jm, km
	    j = lst[j];
	  } // end while loop
        }
      }
    }
    if ((i%10000)==0){
      cout<<"i = "<<setfill('0')<<setw(8)<<i<<"\tNp = "<<cntp;
      print_time();
    }
  }

  cout<<"Finished loops"<<endl;

  for (int i=0; i<nBinsep; i++) {
    vp2[i] = (vp2[i] - vp[i] * vp[i] / Nvp[i]) / (Nvp[i] - 1);
    vp[i] = vp[i] / Nvp[i];
  }

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  pfw.open (fname.c_str(), ios::out);
  for (int i=0; i<nBinsep; i++)
    pfw<<pow(10., i*logDsep + logSepmin)<<"\t"<<setprecision(6)<<fixed<<vp[i]<<"\t"<<vp2[i]<<"\t"<<Nvp[i]<<endl;
  pfw.close();

  cout<<"Finished Code"<<endl;
  return (0);
}


// Change Log
