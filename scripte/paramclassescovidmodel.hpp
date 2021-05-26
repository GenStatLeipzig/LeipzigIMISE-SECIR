//#include <iostream>
  using namespace Rcpp;
class compl_param
{
  public:
    double lam;
  double lam0;
  double lam1;
  double blood_vol;
  double pdcyclophosphamide;
  double pdrelprocarbazine;
  double pdreladriamycine;
  double pdreletoposide;
  double fp;
  double transitt;
  double ststpervol;
  double degrxyz;
  double kpkcprocarbazine1;
  
  //double dt;
};
class tech_param
{
  public:
    double t00;
  double tend;
  double dt;
};
class compl_compartments
{
  public:
    CharacterVector names;
  int  Sc_ind;
  int  Ec_ind;
  int  Dc_ind;
  int  Rc_ind;
  int  StCare_ind;
  int  baddefined;
  
  IntegerVector  IAc_ind;
  IntegerVector  ISc_ind;
  IntegerVector  Cc_ind;
  //mutant:
  int  EMuc_ind;
  IntegerVector  IAMuc_ind;
  IntegerVector  ISMuc_ind;
  
  int IAcnum;
  int  IScnum;
  //mutant:
  int IAMucnum;
  int ISMucnum;
  int Ccnum;
  
  
  int NumOfCompartments;
  int NumofVariables;
  
  
  
    
};

class Compl_variables
{ 
  public:
  int num_var;
  int num_comp;
  
  DoubleVector IAc;//(num_var+10);//12
  DoubleVector ISc;//(num_var+10);//12
  //mutant:
  DoubleVector IAMuc;//(num_var+10);//12
  DoubleVector ISMuc;//(num_var+10);//12
  DoubleVector Cc;//(num_var+10);//12
  
  DoubleVector dIAcdt;//(num_var+10);//12
  DoubleVector dIScdt;//(num_var+10);//12
  DoubleVector dCcdt;//(num_var+10);//12
  
  DoubleVector dAdt;// (num_var);//c++ values are -1 shifted
  //DoubleVector death_exp_ind;
  //DoubleVector crit_exp_ind;
   
  Compl_variables(DoubleVector u,compl_compartments compartments, CharacterVector compnames);//constructor
  ~Compl_variables();//deconstructor
};


class integral_par
{
  public:
  int ksimem_lab;

  //basic parameters:
  double   r1;
  double   r2;
  double   r3;
  double   r4;
  double   rb4;
  double   r5;
  double   r6;
  double   r7;
  double   r8;
  double   r9;
   double   psymp;
  double   pcrit;
  double   pdeath;
  double   influx;
  double   Iainit;
  double   Iamuinit;//mutant
  double mur;//mutant
  int date_mu;//mutant
  double   procentmeas;
  double   influx_start;
  double   blockeff1;
  double   blockeff2;
  double   blockeff3;
  double   blockeff4;
  double   blockeff5;
  double   blockeff6;
  double   blockeff7;
  double   blockeff8;
  double   blockeff9;
  double   blockeff10;
  
  double   unlockeff1;
  double   unlockeff2;
  double   unlockeff3;
  double   unlockeff4;
  double   unlockeff5;
  double   unlockeff6;
  double   unlockeff7;
  double   unlockeff8;
  double   unlockeff9;
  double   blockeff_r2;
  double   powcrit;
  double   death_exp1;
  double   death_exp2;
  double   death_exp3;
  double   death_exp4;
  double   death_exp5;
  double   death_exp6;
  double   death_exp7;
  double   death_exp8;
  double   death_exp9;
  double   death_exp10;
  double   death_exp11;
  double   crit_exp1;
  double   crit_exp2;
  double   crit_exp3;
  double   crit_exp4;
  double   crit_exp5;
  double   crit_exp6;
  double   crit_exp7;
  double   crit_exp8;
  double   crit_exp9;
  double   crit_exp10;
  double   crit_exp11;
  double   crit_exp12;
  double   spdeath;
  double   death_int1;
  double   death_int2;
  double   death_int3;
  double   death_int4;
  double   death_int5;
  double   death_int6;
  double   death_int7;
  double   death_int8;
  double   death_int9;
  double   death_int10;
  double   death_int11;
  double   crit_int1;
  double   crit_int2;
  double   crit_int3;
  double   crit_int4;
  double   crit_int5;
  double   crit_int6;
  double   crit_int7;
  double   crit_int8;
  double   crit_int9;
  double   crit_int10;
  double   crit_int11;
  double   crit_int12;
//special parameters:
  int parsim_label;
  int   delt;
  int   r81;
  int   delaylock;
  int   nsympt;
  int  delayunlock;
  int  unlockinflux;
  int   DelayData;
  double wecumul;
  int   deldistrgam;
  int   parsim_treat;
  int   crit_prior_num;
  int   death_prior_num;
  int   deltapcritpdeath;
  //other parameters:
  int Unlock1;
  int Unlock2;
  int Unlock3;
  int Unlock4;
  int Unlock5;
  int Unlock6;
  int Unlock7;
  int Unlock8;
  int Unlock9;
  int Lock1;
  int Lock2;
  int Lock3;
  int Lock4;
  int Lock5;
  int Lock6;
  int Lock7;
  int Lock8;
  int Lock9;
  int Lock10;
  
  //int DateStart;
  double Sc0;
  int speclength;
  
  //unique to c++-vesion parameters;
  int num_steps;
  int num_locks; 
  int num_unlocks;
  // ??? DateStart
//DoubleVector CMEB_nor;
 

};

//NumericMatrix ind;
//NumericVector num;
//int nummax;//maximal number of subcompartment