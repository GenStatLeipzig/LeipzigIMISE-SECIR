CreateSimAllels<-function(ind_indiv, xtot,main_struct,allel_eff,allel_fr,predt0)
{ #xtot a list of all relevant parameters
  #allel_eff is an effect of the mutant allel on infection rate in % - 30 or 50%
  #allel_fr - the frequency of mutant allel at the last week of observations in %
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct) 
  num_treats <- param_$num_locks+param_$num_unlocks
  ind_max_treat <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("treat",num_treats,sep ='')
  ind_max_treat_num <- (1:length(ind_max_treat))[ind_max_treat]
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct) 
  completedates0 <- main_struct$completedates[[ind_indiv]]
  xtot_long<-xtot
  xtot_long$num_steps<-xtot$num_steps+predt0
  stoch_obj0$num_steps <- xtot_long$num_steps#Only one time!
  xtot_long<-InitParam(x=xtot_long$phiopt$indiv[[ind_indiv]],ind_indiv=ind_indiv,xtot_long,stoch_obj0,main_struct)
  xtot_long$phiopt$indiv[[ind_indiv]] <- UpdateRandomParameters(xtot_long$phiopt$indiv[[ind_indiv]],param_,main_struct)
  
  normal_val <- xtot_long$phiopt$indiv[[ind_indiv]][ind_max_treat_num] 
  
  param_<-ConstructParam(x=xtot_long$phiopt,ind_indiv,main_struct) 
  main_struct$completedates[[ind_indiv]]<-(main_struct$completedates[[ind_indiv]][1]+(0:(xtot_long$num_steps-1)))
  
  sim_res_long<-SimCOVIDStochasticKIAdvTreat(xtot_long,main_struct,ind_indiv)#SimCOVIDStochasticAdvTreat
  sim_res_long<-NameOutput(sim_res_long,main_struct)
  
  xtot_long$phiopt$indiv[[ind_indiv]][ind_max_treat_num]  <- normal_val+log(1+allel_eff/100)
  
  xtot_long<-InitParam(x=xtot_long$phiopt$indiv[[ind_indiv]],ind_indiv,xtot_long,stoch_obj0,main_struct)
  #xtot$num_steps<-xtot$num_steps+predt
  if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))#(main_struct$label_treat_info==1)
  {
    #sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
    #Rcpp:
    #sim_res<-SimulCppCovidModelIndiv(param_,xtot,ind_indiv,main_struct0)
    param_<- IndivUpdateParam(x=xtot_long$phiopt,ind_indiv,main_struct,param_=main_struct$param_)
    sim_res_longallel<-SimulCppCovidModel(param_,xtot_long,ind_indiv,main_struct0)
  }else{
    sim_res_longallel<-SimCOVIDStochastic(xtot_long,main_struct,ind_indiv)#ind_indiv=1
  }
  
  #xtot_long$phiopt$indiv[[ind_indiv]] <- UpdateRandomParameters(xtot_long$phiopt$indiv[[ind_indiv]],param_,main_struct)
  #param_<-ConstructParam(x=xtot_long$phiopt,ind_indiv,main_struct) 
  
  #sim_res_longallel<-SimCOVIDStochasticKIAdvTreat(xtot_long,main_struct,ind_indiv)#SimCOVIDStochasticAdvTreat
  sim_res_longallel<-NameOutput( sim_res_longallel,main_struct)
  
  return(((100-allel_fr)/100)*sim_res_long+(allel_fr/100)*sim_res_longallel)
  
}