#include <Rcpp.h>        
//#include <paramclassesintegral.hpp>
#include <vector>               
#include <math.h>   /* tanh, log */     
#include <covidmodel.cpp>                       
using namespace std;//for cout<<!!!     
using namespace Rcpp;   
double max_of_vect(DoubleVector v); 
double DelayEff(int t0,double par1,double par2,int t1,int delaylock);
NumericMatrix GenerateUpdateTreatStruct(int num_steps, int num_locks, int num_unlocks, integral_par param_ );

//[[Rcpp::export]]
NumericMatrix solvecovidmodel( DoubleVector parval,CharacterVector parnames, int numpar, DoubleVector u0,
                               int num_steps, int numvar, int num_locks, int num_unlocks,
                                CharacterVector compnames,NumericMatrix compind, 
                                DoubleVector compnum, DoubleVector comppar)//NumericMatrix subintervals, NumericMatrix subintervals
{
  //u0 is the initial value!
  // int num_steps, int numvar
  
  integral_par par00; //compl_param
  compl_compartments compartments00; 
  compartments00 = assigncompartments(compnames,compind, compnum, comppar);
  //cout<<"l 70, start"<<endl; 
  par00 = assignparams(parval,parnames, numpar); 
  Compl_variables compl_varstart(u0,compartments00,compnames);
  
  
  if(par00.parsim_label==4)
  {
    par00.r2 = par00.r1;
  }
  if(par00.parsim_label==3)
  {
    par00.r4 = par00.r5;
  }
  if(par00.parsim_label==2)
  {
    par00.r2 = par00.r1;
    par00.r4 = par00.r5;
  }
  DoubleVector u1(u0.length());
  NumericMatrix sim_arr(num_steps,numvar);
  for(int indl =0;  indl < u0.length(); indl++)
  {
    sim_arr(1, indl) = u0(indl);
    u1(indl) = u0(indl);
  }
  return sim_arr;
}



///////////////////////
///////////////////////
double max_of_vect(DoubleVector v)
{
  double max_val=v(0);
  for(int indi = 0; indi< v.length(); indi++)
  {
    max_val=fmax(max_val, v(indi));
  }
  return max_val;
}
NumericMatrix GenerateUpdateTreatStruct(int num_steps, int num_locks, int num_unlocks, integral_par param_ )
{
  NumericMatrix treat_arr(num_steps,2);
  for(int indi = 0; indi< num_steps; indi++)
  {
    treat_arr(indi,1)  = param_.r1;
    treat_arr(indi,1) = param_.r2;
    if(indi>0)
    {
      treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff1, param_.Lock1,param_.delaylock);
      if(num_locks>1)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff2, param_.Lock2,param_.delaylock);
      }
      if(num_locks>2)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff3, param_.Lock3,param_.delaylock);
      }
      if(num_locks>3)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff4, param_.Lock4,param_.delaylock);
      }
      if(num_locks>4)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff5, param_.Lock5,param_.delaylock);
      }
      if(num_locks>5)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff6, param_.Lock6,param_.delaylock);
      }
      if(num_locks>6)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff7, param_.Lock7,param_.delaylock);
      }
      treat_arr(indi,1) = treat_arr(indi,1)*DelayEff(indi,1,param_.unlockeff1,param_.Unlock1,param_.delayunlock);
      
      if(num_unlocks>1)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff2, param_.Unlock2,param_.delayunlock);
      }
      if(num_unlocks>2)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff3, param_.Unlock3,param_.delayunlock);
      }
      if(num_unlocks>3)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff4, param_.Unlock4,param_.delayunlock);
      }
      if(num_unlocks>4)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff5, param_.Unlock5,param_.delayunlock);
      }
      if(num_unlocks>5)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff6, param_.Unlock6,param_.delayunlock);
      }
      if(num_unlocks>6)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff7, param_.Unlock7,param_.delayunlock);
      }
      if(num_unlocks>7)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff8, param_.Unlock8,param_.delayunlock);
      }
      if(num_unlocks>8)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff9, param_.Unlock9,param_.delayunlock);
      }
      /*for(int ind_locks= 1; ind_locks <num_locks; ind_locks++)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1,par2=param_[[paste('blockeff',ind_locks,sep='')]],
                                                             t1=param_[[paste('Lock',ind_locks,sep='')]],
                                                                      param_,delaylock);
      }
      for(int ind_unlocks= 1; ind_unlocks <num_locks; ind_unlocks++)
      {
        treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,par1=1,par2=param_[[paste('unlockeff',ind_unlocks,sep='')]],
                                                             t1= param_[[paste('Unlock',ind_unlocks,sep='')]],
                                                                       delaylock=param_$delayunlock);
      }*/
      
      
      if(param_.parsim_treat==0){
        treat_arr(indi,2) = treat_arr(indi,2)*DelayEff(indi,1,param_.blockeff_r2,param_.Lock2,
                                                              param_.delaylock);
      }
    }
    
    
  }
  
  
 
  
  return treat_arr;//colnames(treat_arr)<-c('r1','r2')
}


double DelayEff(int t0,double par1,double par2,int t1,int delaylock)
{
  double out_=0;
  if(t0<=t1)
  {
    out_ = par1;
  }else{
    if(t0>(t1+delaylock))
    {
      out_ = par2;
    }else{
      out_ = (par1+(par2-par1)*(double)((t0-t1))/((double)delaylock));
    }
  }  
  return out_;
}
