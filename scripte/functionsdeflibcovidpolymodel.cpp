#include <Rcpp.h>
#include <paramclassescovidpolymodel.hpp>
  #include <vector>
  #include <math.h>   /* tanh, log */ 
  
  //!!!!DoubleVector initiation is very important for the output. Without it the output is numeric(0)
//#include <rk4.cpp>
  //#include <rk4.hpp>  
  using namespace std;//for cout<<!!! 
  using namespace Rcpp;
double sum_vect( DoubleVector v);
integral_par assignparams(DoubleVector parval,CharacterVector parnames, int numpar, integral_par par00);
compl_compartments assigncompartments(CharacterVector compnames,NumericMatrix compind, DoubleVector compnum, DoubleVector comppar);
 
//integral_par assignparams_erythro(DoubleVector parval,CharacterVector parnames, int numpar);
DoubleVector dvectcppzerr( int num_var);
DoubleVector TransformParametersVeccpp(DoubleVector parval, CharacterVector parnames, NumericMatrix transform_label_LB_UB, int numpar);
//[[Rcpp::export]]
//////////////////////////constructparamscomp:
  integral_par assignparams(DoubleVector parval,CharacterVector parnames, int numpar, integral_par par00)
{
  // cout<<parnames<<endl;
//integral_par par00;
  for(int ind_par =0; ind_par<numpar; ind_par++)
    {
    //  cout<<"line 20, functionsdeflib; ind_par, numpar = "<<ind_par <<", "<<numpar <<endl;
    if(strcmp(parnames(ind_par), "r1")==0){
      par00.r1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r2")==0){
      par00.r2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r3")==0){
      par00.r3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r4")==0){
      par00.r4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rb4")==0){
      par00.rb4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r5")==0){
      par00.r5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r6")==0){
      par00.r6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r7")==0){
      par00.r7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r8")==0){
      par00.r8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r9")==0){
      par00.r9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "psymp")==0){
      par00.psymp = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "pcrit")==0){
      par00.pcrit = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "pdeath")==0){
      par00.pdeath = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "influx")==0){
      par00.influx = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Iainit")==0){
      par00.Iainit = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Iamuinit")==0){ //Mutant!!!
      par00.Iamuinit = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mur")==0){ //Mutant!!!
      par00.mur = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "date_mu")==0){ //Mutant!!!
      par00.date_mu = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "procentmeas")==0){
      par00.procentmeas = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "influx_start")==0){
      par00.influx_start = parval(ind_par);
    }
    /*if(strcmp(parnames(ind_par), "blockeff1")==0){
      par00.blockeff1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff2")==0){
      par00.blockeff2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff3")==0){
      par00.blockeff3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff4")==0){
      par00.blockeff4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff5")==0){
      par00.blockeff5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff6")==0){
      par00.blockeff6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff7")==0){
      par00.blockeff7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff8")==0){
      par00.blockeff8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff9")==0){
      par00.blockeff9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff10")==0){
      par00.blockeff10 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff1")==0){
      par00.unlockeff1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff2")==0){
      par00.unlockeff2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff3")==0){
      par00.unlockeff3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff4")==0){
      par00.unlockeff4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff5")==0){
      par00.unlockeff5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff6")==0){
      par00.unlockeff6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff7")==0){
      par00.unlockeff7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff8")==0){
      par00.unlockeff8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockeff9")==0){
      par00.unlockeff9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "blockeff_r2")==0){
      par00.blockeff_r2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "powcrit")==0){
      par00.powcrit = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp1")==0){
      par00.death_exp1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp2")==0){
      par00.death_exp2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp3")==0){
      par00.death_exp3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp4")==0){
      par00.death_exp4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp5")==0){
      par00.death_exp5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp6")==0){
      par00.death_exp6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp7")==0){
      par00.death_exp7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp8")==0){
      par00.death_exp8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp9")==0){
      par00.death_exp9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp10")==0){
      par00.death_exp10 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_exp11")==0){
      par00.death_exp11 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp1")==0){
      par00.crit_exp1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp2")==0){
      par00.crit_exp2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp3")==0){
      par00.crit_exp3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp4")==0){
      par00.crit_exp4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp5")==0){
      par00.crit_exp5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp6")==0){
      par00.crit_exp6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp7")==0){
      par00.crit_exp7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp8")==0){
      par00.crit_exp8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp9")==0){
      par00.crit_exp9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp10")==0){
      par00.crit_exp10 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp11")==0){
      par00.crit_exp11 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_exp12")==0){
      par00.crit_exp12 = parval(ind_par);
    }
     */
    if(strcmp(parnames(ind_par), "spdeath")==0){
      par00.spdeath = parval(ind_par);
    }
    /*if(strcmp(parnames(ind_par), "death_int1")==0){
      par00.death_int1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int2")==0){
      par00.death_int2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int3")==0){
      par00.death_int3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int4")==0){
      par00.death_int4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int5")==0){
      par00.death_int5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int6")==0){
      par00.death_int6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int7")==0){
      par00.death_int7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int8")==0){
      par00.death_int8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int9")==0){
      par00.death_int9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int10")==0){
      par00.death_int10 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_int11")==0){
      par00.death_int11 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int1")==0){
      par00.crit_int1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int2")==0){
      par00.crit_int2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int3")==0){
      par00.crit_int3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int4")==0){
      par00.crit_int4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int5")==0){
      par00.crit_int5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int6")==0){
      par00.crit_int6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int7")==0){
      par00.crit_int7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int8")==0){
      par00.crit_int8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int9")==0){
      par00.crit_int9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int10")==0){
      par00.crit_int10 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int11")==0){
      par00.crit_int11 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_int12")==0){
      par00.crit_int12 = parval(ind_par);
    }
     */
    if(strcmp(parnames(ind_par), "parsim_label")==0){
      par00.parsim_label = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "delt")==0){
      par00.delt = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "r81")==0){
      par00.r81 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "delaylock")==0){
      par00.delaylock = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "nsympt")==0){
      par00.nsympt = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "delayunlock")==0){
      par00.delayunlock = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "unlockinflux")==0){
      par00.unlockinflux = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "DelayData")==0){
      par00.DelayData = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "wecumul")==0){
      par00.wecumul = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "deldistrgam")==0){
      par00.deldistrgam = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "parsim_treat")==0){
      par00.parsim_treat = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "crit_prior_num")==0){
      par00.crit_prior_num = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "death_prior_num")==0){
      par00.death_prior_num = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "deltapcritpdeath")==0){
      par00.deltapcritpdeath = parval(ind_par);
    }
    /*if(strcmp(parnames(ind_par), "DateStart")==0){
      par00.DateStart = parval(ind_par);
    }*/
    /*if(strcmp(parnames(ind_par), "Unlock1")==0){
      par00.Unlock1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock2")==0){
      par00.Unlock2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock3")==0){
      par00.Unlock3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock4")==0){
      par00.Unlock4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock5")==0){
      par00.Unlock5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock6")==0){
      par00.Unlock6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock7")==0){
      par00.Unlock7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock8")==0){
      par00.Unlock8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Unlock9")==0){
      par00.Unlock9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock1")==0){
      par00.Lock1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock2")==0){
      par00.Lock2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock3")==0){
      par00.Lock3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock4")==0){
      par00.Lock4 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock5")==0){
      par00.Lock5 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock6")==0){
      par00.Lock6 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock7")==0){
      par00.Lock7 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock8")==0){
      par00.Lock8 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock9")==0){
      par00.Lock9 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "Lock10")==0){
      par00.Lock10 = parval(ind_par);
    }*/
    
    if(strcmp(parnames(ind_par), "mat_contact11")==0){
      par00.mat_contact11 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact12")==0){
      par00.mat_contact12 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact13")==0){
      par00.mat_contact13 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact14")==0){
      par00.mat_contact14 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact22")==0){
      par00.mat_contact22 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact23")==0){
      par00.mat_contact23 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact24")==0){
      par00.mat_contact24 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact33")==0){
      par00.mat_contact33 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact34")==0){
      par00.mat_contact34 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "mat_contact44")==0){
      par00.mat_contact44 = parval(ind_par);
    }
    
    if(strcmp(parnames(ind_par), "rinflux1")==0){
      par00.rinflux1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rinflux2")==0){
      par00.rinflux2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rinflux3")==0){
      par00.rinflux3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rinflux4")==0){
      par00.rinflux4 = parval(ind_par);
    }
    
    if(strcmp(parnames(ind_par), "rtr_1")==0){
      par00.rtr_1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rtr_2")==0){
      par00.rtr_2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rtr_3")==0){
      par00.rtr_3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rtr_4")==0){
      par00.rtr_4 = parval(ind_par);
    }
    
    
    if(strcmp(parnames(ind_par), "prop1_0")==0){
      par00.prop1_0 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "prop2_0")==0){
      par00.prop2_0 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "prop3_0")==0){
      par00.prop3_0 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "prop4_0")==0){
      par00.prop4_0 = parval(ind_par);
    }
    
    if(strcmp(parnames(ind_par), "rdeathage_1")==0){
      par00.rdeathage_1 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rdeathage_2")==0){
      par00.rdeathage_2 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rdeathage_3")==0){
      par00.rdeathage_3 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "rdeathage_4")==0){
      par00.rdeathage_4 = parval(ind_par);
    }
   
    
    
    if(strcmp(parnames(ind_par), "Sc0")==0){
      par00.Sc0 = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "speclength")==0){
      par00.speclength = parval(ind_par);
    }
    //unique to c++ version parameters:
    
    if(strcmp(parnames(ind_par), "num_steps")==0){
      par00.num_steps = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "num_locks")==0){
      par00.num_locks = parval(ind_par);
    }
    if(strcmp(parnames(ind_par), "num_unlocks")==0){
      par00.num_unlocks = parval(ind_par);
    }
  }
  return par00;
  }
compl_compartments assigncompartments(CharacterVector compnames,NumericMatrix compind, DoubleVector compnum, DoubleVector comppar)
{
  compl_compartments compartments00;
  compartments00.NumOfCompartments = comppar(0);//!!! 2 par
  compartments00.NumofVariables = comppar(1);
  int ind_subc;
  for(int ind_comp =0; ind_comp<compartments00.NumOfCompartments; ind_comp++)
  { 
    if(strcmp(compnames(ind_comp), "Sc"  )==0){
      //compartments00.Scnum = compnum(ind_comp);
      compartments00.Sc_ind = compind(ind_comp,0);
    //  cout<<"def comp! l294 "<< compartments00.Sc_ind <<endl;
    }
    if(strcmp(compnames(ind_comp), "Ec"  )==0){
      //compartments00.Ecnum = compnum(ind_comp);
      compartments00.Ec_ind = compind(ind_comp,0);
    } 
    //mutant:
    if(strcmp(compnames(ind_comp), "EMuc"  )==0){
      //compartments00.Ecnum = compnum(ind_comp);
      compartments00.EMuc_ind = compind(ind_comp,0);
    }
    if(strcmp(compnames(ind_comp), "Dc"  )==0){
      //compartments00.Dcnum = compnum(ind_comp);
      compartments00.Dc_ind = compind(ind_comp,0);
    } 
    if(strcmp(compnames(ind_comp), "Rc"  )==0){
      //compartments00.Rcnum = compnum(ind_comp);
      compartments00.Rc_ind = compind(ind_comp,0);
    } 
    if(strcmp(compnames(ind_comp), "StCare"  )==0){
      //compartments00.StCare_indnum = compnum(ind_comp);
      compartments00.StCare_ind = compind(ind_comp,0);
    } 
    
    if(strcmp(compnames(ind_comp), "IAc"  )==0){
      compartments00.IAcnum = compnum(ind_comp);
      DoubleVector temp_v1(compnum(ind_comp));
      for(ind_subc = 0;ind_subc<compnum(ind_comp);ind_subc++){
        temp_v1(ind_subc)=compind(ind_comp,ind_subc);
      }
      compartments00.IAc_ind = temp_v1;
    } 
    if(strcmp(compnames(ind_comp), "ISc"  )==0){
      compartments00.IScnum = compnum(ind_comp);
      DoubleVector temp_v2(compnum(ind_comp));
      for(ind_subc = 0;ind_subc<compnum(ind_comp);ind_subc++){
        temp_v2(ind_subc)=compind(ind_comp,ind_subc);
      }
      compartments00.ISc_ind = temp_v2;
    } 
    if(strcmp(compnames(ind_comp), "Cc"  )==0){
      compartments00.Ccnum = compnum(ind_comp);
      DoubleVector temp_v3(compnum(ind_comp));
      for(ind_subc = 0;ind_subc<compnum(ind_comp);ind_subc++){
        temp_v3(ind_subc)=compind(ind_comp,ind_subc);
      }
      compartments00.Cc_ind = temp_v3;
    } 
    //mutant:
    if(strcmp(compnames(ind_comp), "IAMuc"  )==0){
      compartments00.IAMucnum = compnum(ind_comp);
      DoubleVector temp_v1(compnum(ind_comp));
      for(ind_subc = 0;ind_subc<compnum(ind_comp);ind_subc++){
        temp_v1(ind_subc)=compind(ind_comp,ind_subc);
      }
      compartments00.IAMuc_ind = temp_v1;
    } 
    //mutant:
    if(strcmp(compnames(ind_comp), "ISMuc"  )==0){
      compartments00.ISMucnum = compnum(ind_comp);
      DoubleVector temp_v2(compnum(ind_comp));
      for(ind_subc = 0;ind_subc<compnum(ind_comp);ind_subc++){
        temp_v2(ind_subc)=compind(ind_comp,ind_subc);
      }
      compartments00.ISMuc_ind = temp_v2;
    } 
    
    //cout<<"compartments00.etoposidenum = "<<compartments00.etoposidenum<<endl;
  }
  /////////////////////
    return compartments00;
}                          

Compl_variables::Compl_variables(DoubleVector u,compl_compartments compartments, CharacterVector compnames)
{
  num_comp = compnames.length();
  num_var = u.length();
  //cout<<"functionsdeflibcovidmodel, l 501, num_var =  "<<num_var<<endl;
  dAdt = dvectcppzerr(num_var);
  IAc  = dvectcppzerr(compartments.IAcnum);//(num_var+10);//12
  ISc = dvectcppzerr(compartments.IScnum);//(num_var+10);//12
  for(int ind_comp =0;ind_comp< num_comp; ind_comp++)
  {
    if(compnames(ind_comp)=="EMuc")
    {
      //cout<<"functionsdeflibcovidmodel, l 505, compartments.IAMucnum =  "<<compartments.IAMucnum<<endl;
      IAMuc  = dvectcppzerr(compartments.IAMucnum);//(num_var+10);//12
      //cout<<"functionsdeflibcovidmodel, l 507, compartments.ISMucnum =  "<<compartments.ISMucnum<<endl;
      ISMuc = dvectcppzerr(compartments.ISMucnum);//(num_var+10);//12
      //cout<<"functionsdeflibcovidmodel, l 509, compartments.Ccnum =  "<<compartments.Ccnum<<endl;
    }
  }
  
  Cc = dvectcppzerr(compartments.Ccnum);//(num_var+10);//12
  //cout<<"functionsdeflibcovidmodel, l 511, num_var =  "<<num_var<<endl;
  dIAcdt = dvectcppzerr(compartments.IAcnum);;//(num_var+10);//12
  dIScdt = dvectcppzerr(compartments.IScnum);;//(num_var+10);//12
  dCcdt = dvectcppzerr(compartments.Ccnum);;//(num_var+10);//12
 // cout<<"functionsdeflibcovidmodel, l 515 "<<endl;
  
  
}
Compl_variables::~Compl_variables()//destructor, takes no action
{
  
}
DoubleVector dvectcppzerr( int num_v)
{
  DoubleVector v_out(num_v);
  
  for(int ind_subc = 0;ind_subc<num_v;ind_subc++)
  {
    v_out(ind_subc) = 0.0;
  }
  return v_out;
}

  DoubleVector TransformParametersVeccpp(DoubleVector parval, CharacterVector parnames, NumericMatrix transform_label_LB_UB, int numpar)
  {
    DoubleVector outputval(numpar); 
    for(int ind_par =0; ind_par<numpar; ind_par++)
    {
      if(transform_label_LB_UB(ind_par,0)==0)
      {
        outputval(ind_par) = parval(ind_par);
      }else if(transform_label_LB_UB(ind_par,0)==1)
      {
        outputval(ind_par) = exp(parval(ind_par));
      }else if(transform_label_LB_UB(ind_par,0)==2)
      {
        outputval(ind_par) = transform_label_LB_UB(ind_par,1) + exp(parval(ind_par))*(transform_label_LB_UB(ind_par,2)-transform_label_LB_UB(ind_par,1))/(1+exp(parval(ind_par)));
      }
    }
    return outputval;
  }
 