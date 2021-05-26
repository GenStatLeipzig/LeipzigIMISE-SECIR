#include <Rcpp.h>        
//#include <paramclassesintegral.hpp>
#include <vector>                
#include <math.h>   /* tanh, log */     
#include <covidmodel.cpp>                       
using namespace std;//for cout<<!!!      
using namespace Rcpp;   
double max_of_vect(DoubleVector v); 
double DelayEff(int t0,double par1,double par2,int t1,int delaylock);
NumericMatrix GenerateUpdateTreatStruct( integral_par param_, DoubleVector blockeff_vec,DoubleVector unlockeff_vec,
                                         IntegerVector Lock_vec, IntegerVector Unlock_vec, int num_steps); 
NumericMatrix GenerateUpdateTreatStructDir( integral_par param_, DoubleVector blockeff_vec,DoubleVector unlockeff_vec,
                                            IntegerVector Lock_vec, IntegerVector Unlock_vec, int num_steps);

 
//[[Rcpp::export]]
NumericMatrix solvecovidmodel( DoubleVector parval,CharacterVector parnames, int numpar, DoubleVector u0,
                               int num_steps, int numvar, int num_locks, int num_unlocks,
                                CharacterVector compnames,NumericMatrix compind, 
                                DoubleVector compnum, DoubleVector comppar,
                                DoubleVector crit_exp, DoubleVector death_exp,
                                IntegerVector death_int_round,IntegerVector crit_int_round,
                                DoubleVector blockeff_vec,DoubleVector unlockeff_vec,
                                IntegerVector Lock_vec, IntegerVector Unlock_vec, int treatpardir)//NumericMatrix subintervals, NumericMatrix subintervals
{
  //u0 is the initial value!
  // int num_steps, int numvar
   
  integral_par par00; //compl_param
  compl_compartments compartments00;  
  compartments00 = assigncompartments(compnames,compind, compnum, comppar);
  /*cout<<"l 29, start, compnames= "<<compnames<<endl;
  cout<<"l 30, start, compind= "<<compind<<endl;  
  cout<<"l 31, start, compnum= "<<compnum<<endl; 
  cout<<"l 32, start, comppar= "<<comppar<<endl; */
  par00 = assignparams(parval,parnames, numpar, par00); 
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
  DoubleVector nextav_(u0.length());
  NumericMatrix sim_arr(num_steps,numvar+5);
  for(int indl =0;  indl < u0.length(); indl++)
  {
    sim_arr(0, indl) = u0(indl);
    u1(indl) = u0(indl);
  }
  NumericMatrix treat_arr;
  if(treatpardir==0)
  {
    treat_arr=GenerateUpdateTreatStruct(par00, blockeff_vec, unlockeff_vec, Lock_vec, Unlock_vec, num_steps );
  }else if(treatpardir==1)
  {
    treat_arr=GenerateUpdateTreatStructDir(par00, blockeff_vec, unlockeff_vec, Lock_vec, Unlock_vec, num_steps );
  }
  
  //cout<<"l56, sim_arr(0,0) = "<<sim_arr(0,0)<<", param_.r1,2 = "<<  par00.r1<<", "<<par00.r2<<", treat_arr= "<<treat_arr(0,0)<<","<<treat_arr(0,1)<<", 1:"<<treat_arr(1,0)<<","<<treat_arr(1,1)<<endl;
  int death_exp_n = death_exp.length()+2;
  int crit_exp_n = crit_exp.length()+2;
  double influx0 = par00.influx;
  double pcrit0 = par00.pcrit;
  double pdeath0 = par00.pdeath;
  double spdeath0 = par00.spdeath;
  IntegerVector intervals_death(death_exp_n);
  IntegerVector intervals_crit(crit_exp_n);
  DoubleVector pdeath_v(num_steps);//norm_death_exp_d
  DoubleVector pcrit_v(num_steps);//norm_crit_exp_d
  DoubleVector DICplus(num_steps);
  intervals_death(0) = 1;
  intervals_crit(0) =1;  
  //cout<<"l69"<<endl;///
  for(int ind_e = 0; ind_e <(death_exp_n-2); ind_e++)
  {  
    intervals_death(ind_e+1) = std::min(num_steps,intervals_death(ind_e)+death_int_round(ind_e));  
  } 
  for(int ind_e = 0; ind_e <(crit_exp_n-2); ind_e++)
  {  
    intervals_crit(ind_e+1) = std::min(num_steps,intervals_crit(ind_e)+crit_int_round(ind_e));
  }
  //delta of the last interval:
  //if('deltapcritpdeath'%in%names(param_))
  //{ 
  intervals_crit(crit_exp_n-2) = std::min(num_steps,par00.speclength- par00.deltapcritpdeath+1);//#instead xtot$num_steps-...
  intervals_death(death_exp_n-2) = std::min(num_steps,par00.speclength- par00.deltapcritpdeath+1);//#instead xtot$num_steps-...
  //}
  intervals_death(death_exp_n-1) = num_steps;
  intervals_crit(crit_exp_n-1) = num_steps;
  /*cout<<"l 81, intervals_death(death_exp_n-2) = "<<intervals_death(death_exp_n-2)<<", death_exp_n = "<<death_exp_n<<", num_steps= "<<num_steps<<endl; 
  cout<<"l 82, intervals_crit(crit_exp_n-2) = "<<intervals_crit(crit_exp_n-2)<<", "<<crit_exp_n<<endl;
  cout<<"l 83, intervals_death(death_exp_n-1) = "<<intervals_death(death_exp_n-1)<<endl; 
  cout<<"l 84, intervals_crit(crit_exp_n-1) = "<<intervals_crit(crit_exp_n-1)<<", par00.deltapcritpdeath = "<<par00.deltapcritpdeath<<", par00.speclength- par00.deltapcritpdeath+1 = "<<par00.speclength- par00.deltapcritpdeath+1<<endl;
  cout<<"l 85, par00.Sc0 = "<< par00.Sc0<<" par00.Lock7 = "<< par00.Lock7<<" par00.blockeff_r2 = "<< par00.blockeff_r2<<" par00.unlockeff2 = "<< par00.unlockeff2<<" par00.death_exp8 = "<< par00.death_exp8<<endl;
  cout<<"l 86, par00.crit_int11 = "<< par00.crit_int11<<" par00.crit_int12 = "<< par00.crit_int12<<" par00.parsim_label = "<< par00.parsim_label<<" par00.deltapcritpdeath = "<< par00.deltapcritpdeath<<" par00.Unlock1 = "<< par00.Unlock1<<endl;
  */
  //intervals_crit<-pmin(intervals_crit,par00.num_steps)
  //intervals_death <- pmin(intervals_death,par00.num_steps )
  //intervals_crit[length(intervals_crit)]<-par00.num_steps
  //intervals_death[length(intervals_death)]<-par00.num_steps
  //cout<<"l97"<<endl;
  ///
  for(int ind_e=0; ind_e <(death_exp_n-1);ind_e++)
  {
    if(ind_e==0)
    {
      for(int ind_int=(intervals_death(ind_e)-1); ind_int <(intervals_death(ind_e+1));ind_int++)
      {
        pdeath_v(ind_int) = pdeath0;//1;
      }
      
//norm_death_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)] <- 1
    }else{
      for(int ind_int=(intervals_death(ind_e)-1); ind_int <(intervals_death(ind_e+1));ind_int++)
      {
        pdeath_v(ind_int) = pdeath0*death_exp(ind_e-1);
      }
    }
    
  }
  ///
  
  //cout<<"l119"<<endl;
  for(int ind_e=0; ind_e <(crit_exp_n-1);ind_e++)
  {
    if(ind_e==0)
    {
      for(int ind_int=(intervals_crit(ind_e)-1); ind_int <(intervals_crit(ind_e+1));ind_int++)
      {
        pcrit_v(ind_int) = pcrit0;//1;
      }
      
      //norm_crit_exp_d[((ind_e-1)*num_rep_crit_exp+1):(ind_e*num_rep_crit_exp)] <- 1
    }else{
      for(int ind_int=(intervals_crit(ind_e)-1); ind_int <(intervals_crit(ind_e+1));ind_int++)
      {
        pcrit_v(ind_int) = pcrit0*crit_exp(ind_e-1);
      }
    }
  }
  //cout<<"l137"<<endl;
  double pIA1IS1  = X_p(par00.rb4, par00.r4, par00.psymp );
  sim_arr(0,numvar) = std::max(0.0,u1(compartments00.ISc_ind(0)));
  DICplus(0) = 0;
  sim_arr(0,numvar+1) = treat_arr(0,0);
  sim_arr(0,numvar+2) = treat_arr(0,1);
  sim_arr(0,numvar+3) = pcrit_v(0);
  sim_arr(0,numvar+4) = pdeath_v(0);
  //cout<<"l142, sim_arr(0,0) = "<<sim_arr(0,0)<<endl; 
  for(int indi =1; indi<num_steps; indi++)
  {
    par00.pcrit = pcrit_v(indi);
    par00.pdeath = pdeath_v(indi);
    par00.spdeath = spdeath0*pdeath_v(indi)/pcrit0;
    //cout<<"l148"<<endl;
    sim_arr(indi,numvar+3) = pcrit_v(indi);
    sim_arr(indi,numvar+4) = pdeath_v(indi);
    sim_arr(indi,numvar+1) = treat_arr(indi,0);
    sim_arr(indi,numvar+2) = treat_arr(indi,1);
    //cout<<"l153"<<endl;
    if(indi<par00.influx_start)
    {
      par00.influx = 0;
    }else{
      if(par00.Unlock1<=indi)
      {
        par00.influx = DelayEff(indi, 0, par00.unlockinflux, par00.Unlock1, par00.delayunlock);
      }else
      { 
        par00.influx = DelayEff(indi,influx0,0,Lock_vec(0),par00.delaylock);//Lock_vec instead par00.Lock1!!!
      }
    }
    par00.r1 = treat_arr(indi,0); 
    par00.r2 = treat_arr(indi,1);    
    //nextav_= ModelScholzManyCompNew(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
    //cout<<"l169, par00.r1 = "<<par00.r1<<", par00.r2 = "<<par00.r2<<", par00.influx = "<<par00.influx<<", <par00.Lock1 = "<<par00.Lock1<<endl;//cout<<"l169, compartments00.Sc_ind = "<<compartments00.Sc_ind<<", compartments00.Ec_ind = "<<compartments00.Ec_ind<<endl;//                                delt= param_.delt);  
    nextav_= covidmodel(u1, par00, compartments00,  compnames, compl_varstart);
    double help_v =  pIA1IS1*par00.rb4*std::max(0.0,u1(compartments00.IAc_ind(0)));//previous step
    for(int indl =0;  indl < u0.length(); indl++)
    {
      sim_arr(indi, indl) = u1(indl)+par00.delt*nextav_(indl);//par00.delt*nextav_(indl)  Not necessary yet!!!!
      u1(indl) = sim_arr(indi, indl);
    } 
    DICplus(indi)=help_v+(1-pIA1IS1)*par00.spdeath*std::max(0.0,u1(compartments00.ISc_ind(1)));
    sim_arr(indi,numvar) = sim_arr(indi-1,numvar)+DICplus(indi-1);
    //cout<<"l183,  indi= "<<indi<<", acc eff = "<< sim_arr(indi,numvar)<<", help_v "<<help_v<<", "<<(1-pIA1IS1)*par00.spdeath*std::max(0.0,u1(compartments00.ISc_ind(1)))<<", param_$spdeath= "<<par00.spdeath<<endl;
    /*if(indi<9)
    {  
       cout<<"l176, nextav_ = "<<nextav_<<endl;// cout<<"l176, compartments00.Sc_ind = "<<compartments00.Sc_ind<<", compartments00.Ec_ind"<<endl;//cout<<"l176, nextav_ = "<<nextav_<<endl;
    }*/
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
NumericMatrix GenerateUpdateTreatStruct(integral_par param_, DoubleVector blockeff_vec, DoubleVector unlockeff_vec, IntegerVector Lock_vec, IntegerVector Unlock_vec, int num_steps )
{
  //int num_steps =  param_.num_steps;
  int num_locks =  param_.num_locks;//blockeff_vec.length();//param_.num_locks; 
  int num_unlocks =  param_.num_unlocks;//unlockeff_vec.length();
  
  NumericMatrix treat_arr(num_steps,2);
  for(int indi = 0; indi< num_steps; indi++)
  {
    treat_arr(indi,0)  = param_.r1;
    treat_arr(indi,1) = param_.r2;
    if(param_.parsim_treat==0){
      treat_arr(0,1) = treat_arr(0,1)*DelayEff(0,1,param_.blockeff_r2,param_.Lock2,
                param_.delaylock);
    }else{
      treat_arr(0,1) = param_.blockeff_r2*treat_arr(0,0);
    }
    
    //cout<<"line 205, indi = "<<indi<<", treat_arr = "<<treat_arr(indi,0)<<", "<<treat_arr(indi,1)<<endl; 
    if(indi>0)
    {
      for(int indil = 0; indil< num_locks; indil++)
      {
        treat_arr(indi,0)= treat_arr(indi,0)*DelayEff(indi,1.0, blockeff_vec(indil), Lock_vec(indil),param_.delaylock);
      }
     // cout<<"line 212, indi = "<<indi<<", treat_arr = "<<treat_arr(indi,0)<<", "<<treat_arr(indi,1)<<endl;
      for(int indiu = 0; indiu< num_unlocks; indiu++)
      {
        treat_arr(indi,0)= treat_arr(indi,0)*DelayEff(indi,1.0, unlockeff_vec(indiu), Unlock_vec(indiu),param_.delaylock);
      }
      //cout<<"line 217, indi = "<<indi<<", treat_arr = "<<treat_arr(indi,0)<<", "<<treat_arr(indi,1)<<endl;
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
        treat_arr(indi,1) = treat_arr(indi,1)*DelayEff(indi,1,param_.blockeff_r2,param_.Lock2,
                                                              param_.delaylock);
      }else{
        treat_arr(indi,1) = param_.blockeff_r2*treat_arr(indi,0);
      }
      //cout<<"line 238, indi = "<<indi<<", treat_arr = "<<treat_arr(indi,0)<<", "<<treat_arr(indi,1)<<endl;
    }
    //cout<<"line 240, indi = "<<indi<<", treat_arr = "<<treat_arr(indi,0)<<", "<<treat_arr(indi,1)<<endl;
    
  }
  //cout<<"line 243, param_.r1,2 = "<<  param_.r1<<", "<<param_.r2<<endl;
  
 
  
  return treat_arr;//colnames(treat_arr)<-c('r1','r2')
}

NumericMatrix GenerateUpdateTreatStructDir(integral_par param_, DoubleVector blockeff_vec, DoubleVector unlockeff_vec, IntegerVector Lock_vec, IntegerVector Unlock_vec, int num_steps )
{
  //int num_steps =  param_.num_steps;
  int num_locks =  param_.num_locks;//blockeff_vec.length();//param_.num_locks; 
  int num_unlocks =  param_.num_unlocks;//unlockeff_vec.length();
  
  NumericMatrix treat_arr(num_steps,2);
  treat_arr(0,0)  = param_.r1;
  if(param_.parsim_treat==0){
    treat_arr(0,1) = treat_arr(0,1)*DelayEff(0,1,param_.blockeff_r2,param_.Lock2,
              param_.delaylock);
  }else{
    treat_arr(0,1) = param_.blockeff_r2*treat_arr(0,0);
  }
  int ind_lock = 0;
  int ind_unlock = 0;
  double prev_val =1;
  double help_val;
  for(int indi = 0; indi< num_steps; indi++)
  {
    if((indi>=Lock_vec(ind_lock))&&(indi <= (Lock_vec(ind_lock)+param_.delaylock)))
    {
      help_val = DelayEff(indi,prev_val, blockeff_vec(ind_lock), Lock_vec(ind_lock),param_.delaylock);
      treat_arr(indi,0)= treat_arr(0,0)*help_val;
      treat_arr(indi,1)= treat_arr(0,1)*help_val;
      if(indi == (Lock_vec(ind_lock)+param_.delaylock))
      {
        prev_val = blockeff_vec(ind_lock);
        if(ind_lock<(num_locks-1))
        {
          ind_lock = ind_lock+1;
        }
      }
    }else if((indi>=Unlock_vec(ind_unlock))&&(indi <= (Unlock_vec(ind_unlock)+param_.delayunlock)))
    {
      help_val = DelayEff(indi,prev_val,unlockeff_vec(ind_unlock),
                           Unlock_vec(ind_unlock),
                           param_.delayunlock);
      treat_arr(indi,0)= treat_arr(0,0)*help_val;
      treat_arr(indi,1)= treat_arr(0,1)*help_val;
      if(indi == (Unlock_vec(ind_unlock)+param_.delayunlock))
      {
        prev_val = unlockeff_vec(ind_unlock);
        if(ind_unlock<(num_unlocks-1))
        {
          ind_unlock = ind_unlock+1;
        }
      }
    }else if(indi > 0)//c++!!!
    {
      treat_arr(indi,0) = treat_arr(indi-1,0);
      treat_arr(indi,1) = treat_arr(indi-1,1);
    }
    
    //cout<<"l336, indi treat_arr(indi,0,1) = "<<indi<<", "<<treat_arr(indi,0)<<", "<<treat_arr(indi,1)<<", param_.r1,2 = "<<  param_.r1<<", "<<param_.r2<<", param_.blockeff_r2 = "<<param_.blockeff_r2<<", prev_val = "<<prev_val<<endl;
    
  }
  //cout<<"line 243, param_.r1,2 = "<<  param_.r1<<", "<<param_.r2<<endl;
  
  
  
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
