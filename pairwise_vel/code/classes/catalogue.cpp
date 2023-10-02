#include "catalogue.h"




void halocat::fill_cat(){

  ifstream pf;
  string temp;


  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
	  id.push_back(vec[0]);
	  descid.push_back(vec[1]);
          m200b.push_back(vec[2]);
          vmax.push_back(vec[3]);
          vrms.push_back(vec[4]);
          r200b.push_back(vec[5]);
	  rs.push_back(vec[6]);
	  np.push_back(vec[7]);
          x.push_back(vec[8]);
	  y.push_back(vec[9]);
          z.push_back(vec[10]);
          vx.push_back(vec[11]);
	  vy.push_back(vec[12]);
          vz.push_back(vec[13]);
          pid.push_back(vec[14]);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}

void halocat::fill_cat_v_lin(){

  ifstream pf;
  string temp;


  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
	  id.push_back(vec[0]);
	  descid.push_back(vec[1]);
          m200b.push_back(vec[2]);
          vmax.push_back(vec[3]);
          vrms.push_back(vec[4]);
          r200b.push_back(vec[5]);
	  rs.push_back(vec[6]);
	  np.push_back(vec[7]);
          x.push_back(vec[8]);
	  y.push_back(vec[9]);
          z.push_back(vec[10]);
          vx.push_back(vec[15]);
	  vy.push_back(vec[16]);
          vz.push_back(vec[17]);
          pid.push_back(vec[14]);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}

void halocat::fill_cat_v_nl(){

  ifstream pf;
  string temp;


  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
	  id.push_back(vec[0]);
	  descid.push_back(vec[1]);
          m200b.push_back(vec[2]);
          vmax.push_back(vec[3]);
          vrms.push_back(vec[4]);
          r200b.push_back(vec[5]);
	  rs.push_back(vec[6]);
	  np.push_back(vec[7]);
          x.push_back(vec[8]);
	  y.push_back(vec[9]);
          z.push_back(vec[10]);
          vx.push_back(vec[18]);
	  vy.push_back(vec[19]);
          vz.push_back(vec[20]);
          pid.push_back(vec[14]);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}

void halocat::fill_cat_af(){

  ifstream pf;
  string temp;


  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
	  id.push_back(vec[0]);
	  descid.push_back(-1);
          m200b.push_back(vec[1]);
          r200b.push_back(vec[2]);
          vmax.push_back(vec[3]);
          vrms.push_back(0);
	  rs.push_back(0);
	  np.push_back(0);
          x.push_back(vec[4]);
	  y.push_back(vec[5]);
          z.push_back(vec[6]);
          vx.push_back(vec[7]);
	  vy.push_back(vec[8]);
          vz.push_back(vec[9]);
          pid.push_back(-1);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}

void halocat::fill_cat_af_v_lin(){

  ifstream pf;
  string temp;


  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
	  id.push_back(vec[0]);
	  descid.push_back(-1);
          m200b.push_back(vec[1]);
          r200b.push_back(vec[2]);
          vmax.push_back(vec[3]);
          vrms.push_back(0);
	  rs.push_back(0);
	  np.push_back(0);
          x.push_back(vec[4]);
	  y.push_back(vec[5]);
          z.push_back(vec[6]);
          vx.push_back(vec[10]);
	  vy.push_back(vec[11]);
          vz.push_back(vec[12]);
          pid.push_back(-1);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}

void halocat::fill_cat_af_v_nl(){

  ifstream pf;
  string temp;


  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
	  id.push_back(vec[0]);
	  descid.push_back(-1);
          m200b.push_back(vec[1]);
          r200b.push_back(vec[2]);
          vmax.push_back(vec[3]);
          vrms.push_back(0);
	  rs.push_back(0);
	  np.push_back(0);
          x.push_back(vec[4]);
	  y.push_back(vec[5]);
          z.push_back(vec[6]);
          vx.push_back(vec[7] - vec[10]);
	  vy.push_back(vec[8] - vec[11]);
          vz.push_back(vec[9] - vec[12]);
          pid.push_back(-1);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}


void partcat::fill_cat(){

  ifstream pf;
  string temp;


  pf.open (fname.c_str(), ios::in);
  if (pf){
    while (!pf.eof()){
      while (getline (pf, temp)){
        if (temp[0]!='#' && temp[0]!='r'){
          vector<double> vec;
          char *p;
          char* temp1 = strdup(temp.c_str());
          p = strtok (temp1," W P ,\t\n");
          while (p!= NULL){
            vec.push_back(atof(p));
            p = strtok (NULL," W P ,\t");
          }
          x.push_back(vec[4]);
	  y.push_back(vec[5]);
          z.push_back(vec[6]);
          vx.push_back(vec[7]);
	  vy.push_back(vec[8]);
          vz.push_back(vec[9]);
          vec.clear();
          free (p);
          free (temp1);
        }
      }
    }
    pf.close();
  }
  else
    cout<<"Error opening '"<<fname<<"'!!!"<<endl;
}
