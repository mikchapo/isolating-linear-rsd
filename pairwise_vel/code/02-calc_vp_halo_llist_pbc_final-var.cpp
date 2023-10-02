// Calculate the mean pairwise velocity of halos in separation and mass bins
// v0.1.7, 2022-08-10 - Increased number of mass bins

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


  string input, output, fill_type;
//  int nBinsep;
  cin>>input;
  cin>>output;
  cin>>fill_type;
//  cin>>nBinsep;

  bool all_halos = true;

  halocat halo;
  sfname<<input;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  halo.set_fname (fname);
  if (fill_type=="v-lin") {
    halo.fill_cat_af_v_lin();
  } else if (fill_type=="v-nl") {
    halo.fill_cat_af_v_nl();
  } else {
    halo.fill_cat_af();
  }

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

//  double logMmin = 10.50;
  double logMmax = 15.;
  double logMmin = 12.;
//  double logMmax = 14.;
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

  const double Lbox = 1050.;
  const int M = 100;
  const double l = max(max(dx, dy), dz)/((double)(M));
  vector<int> lst(cnth);
  matrix3d<int> label(M, M, M);
  llist (cnth, M, l, halo.x, halo.y, halo.z, lst, label);


//  cout<<"X: "<<halo.x[0]<<"\t"<<"Y: "<<halo.y[0]<<"\t"<<"Z: "<<halo.z[0]<<"\t"<<"VX: "<<halo.vx[0]<<"\t"<<"VY: "<<halo.vy[0]<<"\t"<<"VZ: "<<halo.vz[0]<<endl;

  ///////////////////////////////////////////////////
  //////////          PAIR COUNTS          //////////
  ///////////////////////////////////////////////////
  matrix<double> vp(nBinsep, nBinM);
  matrix<double> vp2(nBinsep, nBinM);
  matrix<double> Nvp(nBinsep, nBinM);
  int i;
#pragma omp parallel for schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0;i<cnth;i++){
    if ((halo.pid[i]==-1 || all_halos) && log10(halo.m200b[i]) < logMmax && log10(halo.m200b[i]) > logMmin) {
      int binM = (int)((log10(halo.m200b[i])-logMmin)/logDM);
      // i_1 and i_2 give the indices of the furthest boxes to either side along the x-axis that can be within sepmax of the current object
      int i_1=floor(((halo.x[i]-xmin)-sepmax)/l);
      int i_2=floor(((halo.x[i]-xmin)+sepmax)/l);
      // Make sure i_1 and i_2 are within the indices of the list
//      i_1 = max (0,i_1);
//      i_2 = min (M-1,i_2);
      // Repeat for y and z axes
      int j_1=floor(((halo.y[i]-ymin)-sepmax)/l);
      int j_2=floor(((halo.y[i]-ymin)+sepmax)/l);
//      j_1 = max (0,j_1);
//      j_2 = min (M-1,j_2);
      int k_1=floor(((halo.z[i]-zmin)-sepmax)/l);
      int k_2=floor(((halo.z[i]-zmin)+sepmax)/l);
//      k_1 = max (0,k_1);
//      k_2 = min (M-1,k_2);
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
              if (i<j && (halo.pid[j]==-1 || all_halos)){
                int binMj = (int)((log10(halo.m200b[j])-logMmin)/logDM);
                if (binMj == binM) {
                  double pair_sep = sep(halo.x[i], halo.y[i], halo.z[i], halo.x[j] + xoff, halo.y[j] + yoff, halo.z[j] + zoff);
                  if (pair_sep>sepmin && pair_sep<sepmax){
                    double pair_vel = vel(pair_sep, halo.x[i], halo.y[i], halo.z[i], halo.x[j] + xoff, halo.y[j] + yoff, halo.z[j] + zoff, halo.vx[i], halo.vy[i], halo.vz[i], halo.vx[j], halo.vy[j], halo.vz[j]);
                    int binsep = (int)((log10(pair_sep)-logSepmin)/logDsep);
                    if (binsep==nBinsep) {
                      binsep = nBinsep - 1;
                    }
//	            int binsep = (int)((pair_sep-sepmin) / dsep);
//                    cout<<"i: "<<i<<"\t"<<"j: "<<j<<"\t"<<"binsep: "<<binsep<<"\t"<<"binM: "<<binM<<"\t"<<"pair_vel: "<<pair_vel<<endl;
#pragma omp atomic
                    Nvp(binsep, binM) += 1.;
#pragma omp atomic
	            vp(binsep, binM) += pair_vel;
#pragma omp atomic
	            vp2(binsep, binM) += pair_vel*pair_vel;
                  }
                }
	      }
              // Change j to the index of the the previous object placed in box im, jm, km
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

  cout<<"Finished loops"<<endl;

  for (int i=0; i<nBinsep; i++) {
    for (int j=0; j<nBinM; j++) {
//      cout<<"Initial vals for "<<i<<", "<<j<<" - vp: "<<vp(i, j)<<"\tvp2: "<<vp2(i, j)<<"\tNvp: "<<Nvp(i, j)<<endl;
      vp2(i, j) = (vp2(i, j) - vp(i, j) * vp(i, j) / Nvp(i, j)) / (Nvp(i, j) - 1);
      vp(i, j) = vp(i, j) / Nvp(i, j);
//      cout<<"Final vals for "<<i<<", "<<j<<" - vp: "<<vp(i, j)<<"\tvp2: "<<vp2(i, j)<<"\tNvp: "<<Nvp(i, j)<<endl;
    }
  }

  cout<<"Finished calculation"<<endl;

  //////////////////////////////////////////////////////////////
  //////////          STORE COUNTS IN A FILE          //////////
  //////////////////////////////////////////////////////////////
  sfname<<output;
  fname = sfname.str();
  sfname.str("");
  sfname.clear();

  int fill_type_index = fname.find("FILL_TYPE");
  fname.replace(fill_type_index, 9, fill_type);

  cout<<"Added fill_type"<<endl;

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
      pfw<<pow(10., j*logDsep + logSepmin)<<"\t"<<setprecision(6)<<fixed<<vp(j, i)<<"\t"<<vp2(j, i)<<"\t"<<Nvp(j, i)<<endl;
//       pfw<<(j*dsep + sepmin)<<"\t"<<setprecision(6)<<fixed<<vp(j, i)<<"\t"<<Nvp(j, i)<<endl;
    pfw.close();
  }

  cout<<"Finished Code"<<endl;
  return (0);
}


// Change Log
// v0.1.6, 2021-10-05 - Update fill_cat functinos to Aemulus format
// v0.1.5, 2021-08-19 - Debugged mass binning and nan output
// v0.1.4, 2021-08-03 - Added linear and non-linear velocity loading
// v0.1.3, 2021-06-16 - Changed to logarithmic separation binning, changed mass and sep bins limits, added periodic boundary conditions
// v0.1.2, 2021-06-16 - Added mass binning
// v0.1.1, 2021-06-15 - Added check to remove subhalos
// v0.1.0, 2021-06-15 - Copied from 03-comp_eBOSS_RR_llist.cpp and modified
