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

  double   blockeff_r2;
  double   powcrit;
 
 
  double   spdeath;


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
  //poly parameters:
  double   mat_contact11;
  double   mat_contact12;
  double   mat_contact13;
  double   mat_contact14;
  double   mat_contact22;
  double   mat_contact23;
  double   mat_contact24;
  double   mat_contact33;
  double   mat_contact34;
  double   mat_contact44;
  double   rinflux1;
  double   rinflux2;
  double   rinflux3;
  double   rinflux4;
  double   rtr_1;
  double   rtr_2;
  double   rtr_3;
  double   rtr_4;
  double   prop1_0;
  double   prop2_0;
  double   prop3_0;
  double   prop4_0;
  double   rdeathage_1;
  double   rdeathage_2;
  double   rdeathage_3;
  double   rdeathage_4;
    
  
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