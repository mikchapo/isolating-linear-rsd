//#ifndef FUNCTIONS_C__
//#define FUNCTIONS_C__


#include "functions.h"









//////////////////////////////////////////////////////////////////////////////////////////
//////////          COUNT LINES NOT STARTING WITH #          /////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int cnt_lines (const string &fname){

  int ris = 0;

  string line;

  ifstream pf;

  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while(getline (pf, line)){
	if (line[0]!='#')
	  ris++;
      }
    }
    pf.close();
  }
  else
    efopen(fname);

  return (ris);
}






//////////////////////////////////////////////////////////////////////////////////////////
//////////          MESSAGE ERROR OPENING FILE          //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void efopen (const string &fname){
  cout<<"Error opening '"<<fname<<"'"<<endl;
}







//////////////////////////////////////////////////////////////////////////////////////////
//////////          INTERPOLATION          ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void comp_pow_law_interp_pars (const double &x0, const double &x1, const double &y0, const double &y1, double &parrzin, double &pargin){
  pargin = log(y0/y1)/log(x0/x1);
  parrzin = exp(log(x0)-log(y0)/(log(y0/y1)/log(x0/x1)));
}

void comp_lin_interp_pars (const double &x0, const double &x1, const double &y0, const double &y1, double &parmin, double &parqin){
  parmin = (y1-y0)/(x1-x0);
  parqin = y0-(y1-y0)/(x1-x0)*x0;
}

double interp_pl (const double &val, const vector<double> &vx, const vector<double> &vy){

  double risw,gw,r0w,mw,qw;
  
  if (val>(vx[0]-0.0000001) && val<(vx[(int)(vx.size())-2]+0.0000001)){
    for (int iw=0;iw<(int)(vx.size()-1);iw++){
      if (val>(vx[iw]-0.0000001) && val<(vx[iw+1]+0.0000001)){
	comp_pow_law_interp_pars (vx[iw], vx[iw+1], vy[iw], vy[iw+1], r0w, gw);
	//	comp_lin_interp_pars (vx[iw], vx[iw+1], vy[iw], vy[iw+1], mw, qw);
      }
      //    int binw = (int)((val-vec_min(vx))/(vx[1]-vx[0]));
      //    comp_pow_law_interp_pars (vx[binw], vx[binw+1], vy[binw], vy[binw+1], r0w, gw);
    }
  }
  if (val<vx[0])
    comp_pow_law_interp_pars (vx[0], vx[1], vy[0], vy[1], r0w, gw);
   
  if (val>vx[(int)(vx.size()-2)]){
    comp_pow_law_interp_pars (vx[(int)(vx.size()-2)], vx[(int)(vx.size()-1)], vy[(int)(vx.size()-2)], vy[(int)(vx.size()-1)], r0w, gw);
    //    comp_lin_interp_pars (vx[(int)(vx.size()-2)], vx[(int)(vx.size()-1)], vy[(int)(vx.size()-2)], vy[(int)(vx.size()-1)], mw, qw);
  }

  risw = pow ((val/r0w),gw);

  if (risw!=risw || isinf(risw)==1){
    if (val>(vx[0]-0.0000001) && val<(vx[(int)(vx.size())-2]+0.0000001)){
      for (int iw=0;iw<(int)(vx.size()-1);iw++){
	if (val>(vx[iw]-0.0000001) && val<(vx[iw+1]+0.0000001))
	  comp_lin_interp_pars (vx[iw], vx[iw+1], vy[iw], vy[iw+1], mw, qw);	
      }
    }
    if (val<vx[0])
      comp_lin_interp_pars (vx[0], vx[1], vy[0], vy[1], mw, qw);
    if (val>vx[(int)(vx.size()-2)])
      comp_lin_interp_pars (vx[(int)(vx.size()-2)], vx[(int)(vx.size()-1)], vy[(int)(vx.size()-2)], vy[(int)(vx.size()-1)], mw, qw);
    risw = mw*val+qw;
  }
  return (risw);
}




double interp_lin (const double &val, const vector<double> &vx, const vector<double> &vy){

  double risw,mw,qw;
  
  
  if (val>(vx[0]-0.0000001) && val<(vx[(int)(vx.size())-2]+0.0000001)){
    for (int iw=0;iw<(int)(vx.size()-1);iw++){
      if (val>(vx[iw]-0.0000001) && val<(vx[iw+1]+0.0000001))
	comp_lin_interp_pars (vx[iw], vx[iw+1], vy[iw], vy[iw+1], mw, qw);	
    }
  }
  if (val<vx[0])
    comp_lin_interp_pars (vx[0], vx[1], vy[0], vy[1], mw, qw);
  if (val>vx[(int)(vx.size()-2)])
    comp_lin_interp_pars (vx[(int)(vx.size()-2)], vx[(int)(vx.size()-1)], vy[(int)(vx.size()-2)], vy[(int)(vx.size()-1)], mw, qw);
  risw = mw*val+qw;
  
  return (risw);
}







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////          COMPUTE pi          //////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double pi(const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2){
  //  return (fabs(z2-z1));
  return (fabs((x1*x1-x2*x2)+(y1*y1-y2*y2)+(z1*z1-z2*z2))/sqrt((x1+x2)*(x1+x2)+(y1+y2)*(y1+y2)+(z1+z2)*(z1+z2)));
}







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////          COMPUTE rp          //////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double rp (const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2, const double &pi){
  //  return (sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
  return(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)-pi*pi));
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////          COMPUTE sep          //////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double sep (const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2){
  return (sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)));
}


double ang_sep_approx (const double &ra1, const double &dec1, const double &ra2, const double &dec2){
  return (sqrt((ra2-ra1)*(ra2-ra1)+(dec2-dec1)*(dec2-dec1)));
}


double ang_sep (const double &ra1, const double &dec1, const double &ra2, const double &dec2){
  return (acos(sin(dec1*PI/180)*sin(dec2*PI/180)+cos(dec1*PI/180)*cos(dec2*PI/180)*cos((ra1-ra2)*PI/180))*180/PI);
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////          LINKED LIST          /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void llist (const int &N, const int &M, const double &l, const vector<double> &x, const vector<double> &y, const vector<double> &z, vector<int> &lst, matrix3d<int> &label){
  int g,h,p,ind;

  double xmin = vec_min(x);
  double ymin = vec_min(y);
  double zmin = vec_min(z);


  for (g=0;g<N;g++)
    lst[g]=0;
  
  for (g=0;g<M;g++)
    for (h=0;h<M;h++)
      for (p=0;p<M;p++)
	label(g,h,p)=0;
  
  for (ind=0;ind<N;ind++){
    g=floor((x[ind]-xmin)/l);
    h=floor((y[ind]-ymin)/l);
    p=floor((z[ind]-zmin)/l);
    g=min(max(g,0),M-1);
    h=min(max(h,0),M-1);
    p=min(max(p,0),M-1);
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

//////////////////////////////////////////////////////////////////////////////////////////
//////////          PIP Weights          //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void load_bw (const string &fname, const int &nb, vector<long int> &ids_in, matrix<int> &weights){

  ifstream pf;
  string temp;
  pf.open (fname.c_str(), ios::in);
  int _i = 0;
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *temp1 = strdup(temp.c_str());
          char *p = strtok (temp1," ,\t");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," ,\t");
          }
          int _j = 0;
          while (_j<nb){
            weights(_i,_j) = (int)(vec[_j+4]);
            _j++;
          }
          ids_in[_i] = (long int)(vec[0]);
          vec.clear();
          free (p);
          free (temp1);
          _i++;
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening "<<fname<<endl;
}


void load_auw_matrix (const string &input_root, const string &auw_input, const int &nReg, matrix<double> &auw_matrix){

  for (int _k=0;_k<nReg+2;_k++) {
    ifstream pf;
    string temp;
    ostringstream dirname;
    dirname << _k;
    ostringstream sfname;
    sfname<<input_root + dirname.str() + auw_input;
    string fname;
    fname = sfname.str();
    pf.open (fname.c_str(), ios::in);
    int _i = 0;
    if (pf){
      while (!pf.eof()){
        while (getline (pf, temp)){
          if (temp[0]!='#' && temp[0]!='r'){
            vector<double> vec;
            char *temp1 = strdup(temp.c_str());
            char *p = strtok (temp1," ,\t");
            while (p!= NULL){
              vec.push_back(atof(p));
              p = strtok (NULL," ,\t");
            }
            auw_matrix(_k,_i) = vec[1];
            vec.clear();
            free (p);
            free (temp1);
            _i+=1;
          }
        }
      }
      pf.close();
    }
    else
      cout<<"Error opening "<<fname<<endl;
  }
}



//#endif
