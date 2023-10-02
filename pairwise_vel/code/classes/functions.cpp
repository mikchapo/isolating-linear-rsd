//#ifndef FUNCTIONS_C__
//#define FUNCTIONS_C__


#include "functions.h"




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////          COMPUTE sep          //////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double sep (const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2){
  return (sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)));
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////          COMPUTE pairwise vel          //////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double vel (const double &sep, const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2, const double &vx1, const double &vy1, const double &vz1, const double &vx2, const double &vy2, const double &vz2){
  return (((vx1-vx2)*(x1-x2)+(vy1-vy2)*(y1-y2)+(vz1-vz2)*(z1-z2))/sep);
//  double pair_vel = (((vx1-vx2)*(x1-x2)+(vy1-vy2)*(y1-y2)+(vz1-vz2)*(z1-z2))/sep);
//  if (isnan(pair_vel)) {
//    cout<<"NaN Vel!!"<<endl;
//    cout<<"vx1: "<<vx1<<"\tvx2: "<<vx2<<"\tx1: "<<x1<<"\tx2: "<<x2<<"\tvy1: "<<vy1<<"\tvy2: "<<vy2<<"\ty1: "<<y1<<"\ty2: "<<y2<<"\tvz1: "<<vz1<<"\tvz2: "<<vz2<<"\tz1: "<<z1<<"\tz2: "<<z2<<"\tsep: "<<sep<<endl;
//  }
//  return pair_vel;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////          LINKED LIST          /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void llist (const int &N, const int &M, const double &l, const vector<double> &x, const vector<double> &y, const vector<double> &z, vector<int> &lst, matrix3d<int> &label){
  int g,h,p,ind;
  // Parameters
  // N - number of objects in the catalogue
  // M - number of boxes along each axis
  // l - length of a single box along the longest axis
  // x - vector of the x-coordinates
  // y - vector of the y-coordinates
  // z - vector of the z-coordinates
  // lst - vector of length N
  // label - 3d matrix of size (M, M, M)


  // Minimum values for each coordinate
  double xmin = vec_min(x);
  double ymin = vec_min(y);
  double zmin = vec_min(z);

  // Initialize lst to 0 for all entries
  for (g=0;g<N;g++)
    lst[g]=0;

  // Initialize label to 0 for all entries
  for (g=0;g<M;g++)
    for (h=0;h<M;h++)
      for (p=0;p<M;p++)
	label(g,h,p)=0;

  for (ind=0;ind<N;ind++){
    // Find the box indices g, h, p for each axis for object ind
    g=floor((x[ind]-xmin)/l);
    h=floor((y[ind]-ymin)/l);
    p=floor((z[ind]-zmin)/l);
    // Make sure each index is within the proper limits (which it should be by definition)
    g=min(max(g,0),M-1);
    h=min(max(h,0),M-1);
    p=min(max(p,0),M-1);
    // label stores the index of the most recent object to be placed in it, and lst stores the index of the object to go in the box before this one
    lst[ind]=label(g,h,p);
    label(g,h,p)=ind;
  }
  //  cout<<label[4][1][2]<<endl;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////          PRINT CURRENT TIME          //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void print_time(){

  //  cout<<time()<<endl;
  time_t currentTime;
  struct tm *localTime;

  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time

  int Day    = localTime->tm_mday;
  int Month  = localTime->tm_mon + 1;
  int Year   = localTime->tm_year + 1900;
  int Hour   = localTime->tm_hour;
  int Min    = localTime->tm_min;
  int Sec    = localTime->tm_sec;
  cout << "\t\t\tTime: " << Hour << ":" << Min << ":" << Sec << "\t"<< Day << "/" << Month << "/" << Year << std::endl;

}



void get_time(int &Day, int &Month, int &Year, int &Hour, int &Min, int &Sec){

  //  cout<<time()<<endl;
  time_t currentTime;
  struct tm *localTime;

  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time

  Day    = localTime->tm_mday;
  Month  = localTime->tm_mon + 1;
  Year   = localTime->tm_year + 1900;
  Hour   = localTime->tm_hour;
  Min    = localTime->tm_min;
  Sec    = localTime->tm_sec;

  }

//#endif
