// Calculate the mean pairwise velocity of halos in separation and mass bins
// v0.1.3, 2021-06-16 - Changed to logarithmic separation binning, changed mass and sep bins limits, added periodic boundary conditions

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
//  int M;
  cin>>input;
//  cin>>M;
  cin>>output;

  bool all_halos = true;

  halocat halo;
  sfname<<input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  halo.set_fname (fname);
  halo.fill_cat();

  int cnth = (int)(halo.id.size());
  cout<<"N_halo: "<<cnth<<endl;

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

  double logMmin = 10.50;
  double logMmax = 15.;
  double logDM = 0.50;
  int nBinM = round((logMmax - logMmin) / logDM);

  /////////////////////////////////////////////////////
  //////////         LINKED LIST          /////////////
  /////////////////////////////////////////////////////
  double xmin = vec_min (halo.x);
  double ymin = vec_min (halo.y);
  double zmin = vec_min (halo.z);

  double xmax = vec_max (halo.x);
  double ymax = vec_max (halo.y);
  double zmax = vec_max (halo.z);

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  const int M = 100;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cnth);
  matrix3d<int> label(M, M, M);
  llist (cnth, M, l, halo.x, halo.y, halo.z, lst, label);




  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix<double> vp(nBinsep, nBinM);
  matrix<double> Nvp(nBinsep, nBinM);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cnth;i++){
    if (halo.pid[i]==-1 || all_halos) {
      int binM = (int)((log10(halo.m200b[i])-logMmin)/logDM);
      // i_1 and i_2 give the indices of the furthest boxes to either side along the x-axis that can be within sepmax of the current object
      int i_1=floor(((halo.x[i]-xmin)-sepmax)/l);
      int i_2=floor(((halo.x[i]-xmin)+sepmax)/l);
      // Make sure i_1 and i_2 are within the indices of the list
      i_1 = max (0,i_1);
      i_2 = min (M-1,i_2);
      // Repeat for y and z axes
      int j_1=floor(((halo.y[i]-ymin)-sepmax)/l);
      int j_2=floor(((halo.y[i]-ymin)+sepmax)/l);
      j_1 = max (0,j_1);
      j_2 = min (M-1,j_2);
      int k_1=floor(((halo.z[i]-zmin)-sepmax)/l);
      int k_2=floor(((halo.z[i]-zmin)+sepmax)/l);
      k_1 = max (0,k_1);
      k_2 = min (M-1,k_2);
      // Loop over all boxes within sepmax of the target
      for (int il=i_1;il<=i_2;il++){
        for (int jl=j_1;jl<=j_2;jl++){
          for (int kl=k_1;kl<=k_2;kl++){
            // Start j as the index of the most recent object to be placed in box il, jl, kl
            int j=label(il,jl,kl);
            // If j==0 then all objects in the box have been looped over and the loop moves to the next box
            while (j!=0){
              if (i<j && (halo.pid[j]==-1 || all_halos)){
                int binMj = (int)((log10(halo.m200b[j])-logMmin)/logDM);
                if (binMj == binM) {
                  double pair_sep = sep(halo.x[i], halo.y[i], halo.z[i], halo.x[j], halo.y[j], halo.z[j]);
                  if (pair_sep>sepmin && pair_sep<sepmax){
                    double pair_vel = vel(pair_sep, halo.x[i], halo.y[i], halo.z[i], halo.x[j], halo.y[j], halo.z[j], halo.vx[i], halo.vy[i], halo.vz[i], halo.vx[j], halo.vy[j], halo.vz[j]);
                    int binsep = (int)((log10(pair_sep)-logSepmin)/logDsep);
//	            int binsep = (int)((pair_sep-sepmin) / dsep);
#pragma omp atomic
	            vp(binsep, binM) = vp(binsep, binM) + pair_vel;
#pragma omp atomic
                    Nvp(binsep, binM) += 1.;
                  }
                }
	      }
              // Change j to the index of the the previous object placed in box il, jl, kl
	      j = lst[j];
	    } // end while loop
          }
        }
      }
    }
    if ((i%100000)==0){
      cout<<"i = "<<setfill('0')<<setw(8)<<i<<"\tNh = "<<cnth;
      print_time();
    }
  }

  for (int i=0; i<nBinsep; i++) {
    for (int j=0; j<nBinM; j++) {
      vp(i, j) = vp(i, j) / Nvp(i, j);
    }
  }


  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();
  int mass_index = fname.find("MASS");
  for (int i=0; i<nBinM; i++) {
    ostringstream mass_min, mass_max;
    mass_min << fixed << setprecision(2) << (i*logDM + logMmin);
    mass_max << fixed << setprecision(2) << ((i+1)*logDM + logMmin);
    string mass_label = mass_min.str() + "-" + mass_max.str();
    if (i==0) {
      fname.replace(mass_index, 4, mass_label);
    } else {
      fname.replace(mass_index, 11, mass_label);
    }
    pfw.open (fname.c_str(), ios::out);
    for (int j=0; j<nBinsep; j++)
      pfw<<pow(10., j*logDsep + logSepmin)<<"\t"<<setprecision(6)<<fixed<<vp(j, i)<<"\t"<<Nvp(j, i)<<endl;
//       pfw<<(j*dsep + sepmin)<<"\t"<<setprecision(6)<<fixed<<vp(j, i)<<"\t"<<Nvp(j, i)<<endl;
    pfw.close();
  }

  return (0);
}


// Change Log
// v0.1.2, 2021-06-16 - Added mass binning
// v0.1.1, 2021-06-15 - Added check to remove subhalos
// v0.1.0, 2021-06-15 - Copied from 03-comp_eBOSS_RR_llist.cpp and modified
