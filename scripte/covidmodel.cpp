#include <Rcpp.h>
//#include <functionslib.cpp> 
#include <functionsdeflibcovidmodel.cpp>// After revision!!!//double include!?  functionsdeflib.cpp functionsdeflibnewrevised1.cpp

#include <vector>
#include <math.h>   /* tanh, log */ 

//!!!!DoubleVector initiation is very important for the output. Without it the output is numeric(0)
//#include <rk4.cpp>
//#include <rk4.hpp> 
using namespace std;//for cout<<!!! 
using namespace Rcpp;

//DoubleVector in_out_from_amplif(double a);
double X_p(double r_1,double r_2,double p_1);
//[[Rcpp::export]]
//////////////////////////constructparamscomp:
  DoubleVector covidmodel(DoubleVector u, integral_par param_, compl_compartments compartments, CharacterVector compnames, Compl_variables compl_varstart)  
{
  //cout<<"integralode, line 342, param_.Cmkcnorsum"<<endl;
  int num_var = u.length();
  int num_comp = compnames.length();
  
    
    
  DoubleVector dAdt (num_var);//(num_var+1)c++ values are -1 shifted
  //cout<<"num_var+1 = "<<num_var+1<<", dAdt.length()"<<dAdt.length()<<endl;
  //DoubleVector KsiPred(2);
  for(int ind_subc = 0;ind_subc<(num_var);ind_subc++)
  {
    dAdt(ind_subc) = 0; 
    //dAdt(ind_subc)=-10000000000000000;
  }
  
  
  
    
 
  // cout<<"integralode, line 484"<<endl;

  for(int ind_subc = 0;ind_subc<compartments.IAcnum;ind_subc++)
  {
    compl_varstart.IAc(ind_subc)=std::max(0.0,u(compartments.IAc_ind(ind_subc)));
  }
  for(int ind_subc = 0;ind_subc<compartments.IScnum;ind_subc++)
  {
    compl_varstart.ISc(ind_subc)=std::max(0.0,u(compartments.ISc_ind(ind_subc)));
  } 
  for(int ind_subc = 0;ind_subc<compartments.Ccnum;ind_subc++)
  {
    compl_varstart.Cc(ind_subc)=std::max(0.0,u(compartments.Cc_ind(ind_subc)));
  }
  double Sc=std::max(0.0,u(compartments.Sc_ind));///demargfactor;//circulated band
  double Ec=std::max(0.0,u(compartments.Ec_ind));///demargfactor;//circulated segmented
  double Dc=std::max(0.0,u(compartments.Dc_ind));
  double Rc=std::max(0.0,u(compartments.Rc_ind));
  double StCare = std::max(0.0,u(compartments.StCare_ind));
  double pIA1IS1  = X_p(param_.rb4, param_.r4, param_.psymp );
  double pIS1Cc1 = param_.pcrit;
  double pCc1Dc = param_.pdeath;
  
  double dScdt = -param_.influx -param_.r1*(Sc/param_.Sc0)*(compl_varstart.IAc(0)+compl_varstart.IAc(1)+compl_varstart.IAc(2)) - param_.r2*(Sc/param_.Sc0)*(compl_varstart.ISc(0)+compl_varstart.ISc(1)+compl_varstart.ISc(2));
  double dEcdt = param_.r1*(Sc/param_.Sc0)*(compl_varstart.IAc(0)+compl_varstart.IAc(1)+compl_varstart.IAc(2)) + param_.r2*(Sc/param_.Sc0)*(compl_varstart.ISc(0)+compl_varstart.ISc(1)+compl_varstart.ISc(2))-param_.r3*Ec;
  compl_varstart.dIAcdt(0) = param_.influx + param_.r3*Ec - (pIA1IS1*param_.rb4 + (1-pIA1IS1)*param_.r4)*compl_varstart.IAc(0);
  compl_varstart.dIAcdt(1) = (1-pIA1IS1)*param_.r4*compl_varstart.IAc(0)-param_.r4*compl_varstart.IAc(1);
  compl_varstart.dIAcdt(2) = param_.r4*compl_varstart.IAc(1)-param_.r4*compl_varstart.IAc(2); 
  compl_varstart.dIScdt(0) = pIA1IS1*param_.rb4*compl_varstart.IAc(0)-(pIS1Cc1*param_.r6 + (1-pIS1Cc1)*param_.r5)*compl_varstart.ISc(0);
  compl_varstart.dIScdt(1) = (1-pIS1Cc1)*param_.r5*compl_varstart.ISc(0)-(1-param_.spdeath)*param_.r5*compl_varstart.ISc(1)-param_.r8*param_.spdeath*compl_varstart.ISc(1);
  compl_varstart.dIScdt(2) = (1-param_.spdeath)*param_.r5*compl_varstart.ISc(1)-param_.r5*compl_varstart.ISc(2);
  compl_varstart.dCcdt(0) = pIS1Cc1*param_.r6*compl_varstart.ISc(0)-(pCc1Dc*param_.r8 + (1-pCc1Dc)*param_.r7)*compl_varstart.Cc(0);
  compl_varstart.dCcdt(1) = (1-pCc1Dc)*param_.r7*compl_varstart.Cc(0)-param_.r7*compl_varstart.Cc(1);
  compl_varstart.dCcdt(2) = param_.r7*compl_varstart.Cc(1)-param_.r7*compl_varstart.Cc(2);
  double dDcdt = pCc1Dc*param_.r8*compl_varstart.Cc(0)+param_.r8*param_.spdeath*compl_varstart.ISc(1);
  double dStCaredt = param_.r7*compl_varstart.Cc(2) - param_.r9*StCare;
  double dRcdt = param_.r5*compl_varstart.ISc(2) + param_.r9*StCare + param_.r4*compl_varstart.IAc(2);
  
        //cout<<"covidmodel, line 77, dEcdt = "<< dEcdt<<", dScdt = "<<dScdt<<"yyy = "<<param_.r1*(Sc/param_.Sc0)*(compl_varstart.IAc(0)+compl_varstart.IAc(1)+compl_varstart.IAc(2))<<", param_.r1 ="<<param_.r1<<", (Sc/param_.Sc0)= "<<(Sc/param_.Sc0)<<", param_.Sc0= "<<param_.Sc0<<", (compl_varstart.IAc(0)+compl_varstart.IAc(1)+compl_varstart.IAc(2)) = "<<(compl_varstart.IAc(0)+compl_varstart.IAc(1)+compl_varstart.IAc(2))<<endl;
        dAdt(compartments.Ec_ind)=dEcdt;
        dAdt(compartments.Sc_ind)=dScdt;
        dAdt(compartments.Dc_ind)=dDcdt;
        dAdt(compartments.Rc_ind)= dRcdt;
        dAdt(compartments.StCare_ind)=dStCaredt;
        for (int ind_comp=0;ind_comp < compartments.IAcnum; ind_comp++){
          dAdt(compartments.IAc_ind(ind_comp))= compl_varstart.dIAcdt(ind_comp);
        }
        for (int ind_comp=0;ind_comp < compartments.IScnum; ind_comp++){
        dAdt(compartments.ISc_ind(ind_comp))= compl_varstart.dIScdt(ind_comp);
        }
        for (int ind_comp=0;ind_comp < compartments.Ccnum; ind_comp++){
        dAdt(compartments.Cc_ind(ind_comp))= compl_varstart.dCcdt(ind_comp);
        }
  
  
        
        
        return dAdt;//   !!!!dAdt2;//u2;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
  
  double X_p(double r_1,double r_2,double p_1)
  {
    return r_2*p_1/(r_1*(1-p_1)+r_2*p_1);//out_=x+y+c*x*y/(abs(x)+abs(y));
  }