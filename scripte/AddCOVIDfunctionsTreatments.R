
SimCOVIDStochasticKIAdvTreat<- function(xtot,main_struct,ind_indiv)
{
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct) 
  # param_$speclength<- main_struct$speclength
  if(!'spdeath'%in%names(param_))
  {
    param_$spdeath<-0
  }
  param_$spdeath0 <- param_$spdeath
  if(param_$parsim_label==4)
  {
    param_$r2 <- param_$r1
  }
  if(param_$parsim_label==3)
  {
    param_$r4 <- param_$r5
  }
  if(param_$parsim_label==2)
  {
    param_$r2 <- param_$r1
    param_$r4 <- param_$r5
  }
  sim_arr<-array(dim=c(xtot$num_steps,main_struct$numvar))
  sim_arr[1,]<- xtot$init
  if(main_struct$LogSc==1)
  {
    sim_arr[1,main_struct$compartments[[ind_indiv]]$Sc_ind]<- log(xtot$init[main_struct$compartments[[ind_indiv]]$Sc_ind])
  }
  sim_darray<-array(dim=c(xtot$num_steps,main_struct$numvar))
  sim_darray[1, ]<- rep.int(x= 0, times= main_struct$numvar)
  if(main_struct$treatpardir==2)
  {
    treat_arr <- GenerateTreatStructKIDir(xtot,main_struct,ind_indiv,param_)
  }
  else if(main_struct$treatpardir==1)
  {
    treat_arr <- GenerateTreatStructDir(xtot,main_struct,ind_indiv,param_)
  }else if(main_struct$treatpardir==0)
  {  
    treat_arr <- GenerateTreatStruct(xtot,main_struct,ind_indiv,param_)
  }
  ###
  
  if("DateStepsPDeath"%in%names(main_struct$covariates))
  {
    death_exp_ind<-grep(pattern='death_exp',x=names(param_))
    crit_exp_ind<-grep(pattern='crit_exp',x=names(param_))
    death_exp_n <- main_struct$num_pdeath+1
    crit_exp_n <- main_struct$num_pcrit+1
    
    norm_death_exp_d <- vector(length=length(main_struct$completedates[[ind_indiv]]))
    norm_crit_exp_d <- vector(length=length(main_struct$completedates[[ind_indiv]]))
    for(ind_e in 1:death_exp_n)
    {
      if(ind_e==1)
      {
        good_ind<-(main_struct$covariates$DateStepsPDeath[[ind_indiv]][1]>=main_struct$completedates[[ind_indiv]])
        norm_death_exp_d[good_ind] <- 1
      }else{
        good_ind<-(main_struct$covariates$DateStepsPDeath[[ind_indiv]][ind_e-1]<main_struct$completedates[[ind_indiv]])
        if(ind_e <death_exp_n)
        {  
          good_ind <- good_ind&(main_struct$covariates$DateStepsPDeath[[ind_indiv]][ind_e]>=main_struct$completedates[[ind_indiv]])
        }
        norm_death_exp_d[good_ind]<-param_[[death_exp_ind[ind_e-1]]]
        #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
      }
      
    }
    
    for(ind_e in 1:crit_exp_n)
    {
      if(ind_e==1)
      {
        good_ind<-(main_struct$covariates$DateStepsPCrit[[ind_indiv]][1]>=main_struct$completedates[[ind_indiv]])
        norm_crit_exp_d[good_ind] <- 1
      }else{
        good_ind<-(main_struct$covariates$DateStepsPCrit[[ind_indiv]][ind_e-1]<main_struct$completedates[[ind_indiv]])
        if(ind_e <crit_exp_n)
        {  
          good_ind <- good_ind&((main_struct$covariates$DateStepsPCrit[[ind_indiv]][ind_e]>=main_struct$completedates[[ind_indiv]]))
        }
        norm_crit_exp_d[good_ind]<-param_[[crit_exp_ind[ind_e-1]]]
        #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
      }
      
    }
    
  }else{
    
    death_exp_ind<- main_struct$death_exp_ind# grep(pattern='death_exp',x=names(param_))
    crit_exp_ind<-main_struct$crit_exp_ind# grep(pattern='crit_exp',x=names(param_))
    death_exp_n <- main_struct$death_exp_n# length(death_exp_ind)+2
    crit_exp_n <- main_struct$crit_exp_n# length(crit_exp_ind)+2
    #norm_death_exp_d<-main_struct$age_struct$norm_death_exp_d
    #norm_crit_exp_d<-main_struct$age_struct$norm_death_exp_d
    norm_death_exp_d <- vector(length = xtot$num_steps)
    norm_crit_exp_d <- vector(length = xtot$num_steps)
    intervals_death <- vector(length = death_exp_n)
    intervals_crit <- vector(length = crit_exp_n)
    intervals_death[1] <- 1
    intervals_crit[1] <- 1
    if("death_int1"%in%names(param_)||('dates_treatments'%in%names(main_struct)))
    {
      help_str <- KICreateIntervals(param_,xtot,main_struct)
      intervals_crit <- help_str$intervals_crit
      intervals_death <- help_str$intervals_death
    }else{
      label_KI <- 0
      num_rep_death_exp<-floor(main_struct$speclength/death_exp_n)
      num_rep_crit_exp<-floor(main_struct$speclength/crit_exp_n)
      for(ind_e in 1:(death_exp_n-2))
      {  
        intervals_death[ind_e+1]<- ind_e*num_rep_death_exp
      } 
      for(ind_e in 1:(crit_exp_n-2))
      {  
        intervals_crit[ind_e+1]<- ind_e*num_rep_crit_exp
      }
      intervals_death[death_exp_n]<- main_struct$speclength
      intervals_crit[death_exp_n]<- main_struct$speclength
    }
    
    if(death_exp_n>0)
    {
      #num_rep_death_exp<-floor(main_struct$speclength/death_exp_n)#floor(length(norm_death_exp_d)/death_exp_n)
      for(ind_e in 1:(death_exp_n-1))
      {
        if(ind_e==1)
        {
          norm_death_exp_d[intervals_death[ind_e]:intervals_death[ind_e+1]] <- 1
          #norm_death_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)] <- 1
        }else{
          norm_death_exp_d[intervals_death[ind_e]:intervals_death[ind_e+1]] <- param_[[death_exp_ind[ind_e-1]]]
          #norm_death_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[death_exp_ind[ind_e-1]]]
          #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
        }
        
      }
      #norm_death_exp_d[(death_exp_n*num_rep_death_exp):length(main_struct$age_struct$norm_death_exp_d)]<-param_[[death_exp_ind[death_exp_n-1]]]
      #norm_crit_exp_d[(death_exp_n*num_rep_death_exp):length(main_struct$age_struct$norm_death_exp_d)]<-param_[[norm_crit_exp_d[death_exp_n]]]
    }
    #
    if(crit_exp_n>0)
    {
      #num_rep_crit_exp<-floor(main_struct$speclength/crit_exp_n)#floor(length(norm_crit_exp_d)/crit_exp_n)
      
      for(ind_e in 1:(crit_exp_n-1))
      {
        if(ind_e==1)
        {
          norm_crit_exp_d[intervals_crit[ind_e]:intervals_crit[ind_e+1]] <- 1
          #norm_crit_exp_d[((ind_e-1)*num_rep_crit_exp+1):(ind_e*num_rep_crit_exp)] <- 1
        }else{
          norm_crit_exp_d[intervals_crit[ind_e]:intervals_crit[ind_e+1]] <- param_[[crit_exp_ind[ind_e-1]]]
          #norm_crit_exp_d[((ind_e-1)*num_rep_crit_exp+1):(ind_e*num_rep_crit_exp)]<-param_[[crit_exp_ind[ind_e-1]]]
        }
        
        #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
      }
      #norm_death_exp_d[(death_exp_n*num_rep_death_exp):length(main_struct$age_struct$norm_death_exp_d)]<-param_[[death_exp_ind[death_exp_n]]]
      #norm_crit_exp_d[(crit_exp_n*num_rep_crit_exp):length(norm_crit_exp_d)]<-param_[[crit_exp_ind[crit_exp_n-1]]]
    }
  }
  if('date_mu'%in%colnames(main_struct$dates_treatments))
  {  
    mut_step <-as.numeric(as.Date(main_struct$dates_treatments[ind_indiv,'date_mu'],format = "%d.%m.%Y")-param_$DateStart)
  }else{
    mut_step <- (-1)
  }
  param_$influx0<-param_$influx
  param_$pcrit0 <- param_$pcrit
  param_$pdeath0 <- param_$pdeath
  #pcrit_v<-pmin(1,param_$pcrit0*main_struct$age_struct$critical_exp_d)
  pcrit_v <- pmin(0.95,param_$pcrit0*norm_crit_exp_d)#^param_$powcrit0.5
  pdeath_v <- pmin(0.66,param_$pdeath0*norm_death_exp_d)#^(1-param_$powcrit) 0.5
  if(length(pcrit_v)<xtot$num_steps)
  {
    pcrit_v[(length(pcrit_v)+1):dim(sim_arr)[1]] <- pcrit_v[(length(pcrit_v))]
    pdeath_v[(length(pdeath_v)+1):dim(sim_arr)[1]] <- pdeath_v[(length(pdeath_v))]
  }
  #if(max(param_$pcrit0*main_struct$age_struct$critical_exp_d)>1)
  #{
  #   pcrit_v<-param_$pcrit0*main_struct$age_struct$critical_exp_d/max(pcrit_v)
  # } 
  #pdeath_v<- pmin(1,param_$pdeath0*main_struct$age_struct$death_exp_d)#pmin(1,param_$pdeath0*main_struct$age_struct$death_exp_d/pcrit_v)#param_$pdeath0*main_struct$age_struct$death_exp_d/pcrit_v
  
  # if(max(pdeath_v)>1)
  # {
  #   pdeath_v <- pdeath_v/max(pdeath_v)
  # }
  for(indi in 2:xtot$num_steps){
    if(indi == (mut_step+1))#if(indi >= (mut_step+1))#if(indi == (mut_step+1))
    {
      #Variant1:
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[1]] <- param_$Iamuinit
      #Variant2 accumulating:
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[1]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[1]] + param_$Iamuinit/7
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[2]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[2]] + param_$Iamuinit/7
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[3]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[3]] + param_$Iamuinit/7
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[1]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[1]] + param_$Iamuinit/7
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[2]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[2]] + param_$Iamuinit/7
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[3]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[3]] + param_$Iamuinit/7
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$EMuc_ind] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$EMuc_ind]+param_$Iamuinit/7
      #Variant3 (proportional):
      sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[1]] <- param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAc_ind[1]]
      sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[2]] <- param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAc_ind[2]]
      sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[3]] <- param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAc_ind[3]]
      sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[1]] <- param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISc_ind[1]]
      sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[2]] <- param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISc_ind[2]]
      sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[3]] <- param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISc_ind[3]]
      sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$EMuc_ind] <- param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$Ec_ind]
      
      #Variant3 (influx, proportional):
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[1]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[1]] + param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAc_ind[1]]
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[2]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[2]] + param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAc_ind[2]]
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[3]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAMuc_ind[3]] + param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$IAc_ind[3]]
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[1]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[1]] + param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISc_ind[1]]
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[2]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[2]] + param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISc_ind[2]]
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[3]] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISMuc_ind[3]] + param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$ISc_ind[3]]
      #sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$EMuc_ind] <- sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$EMuc_ind] + param_$Iamuinit*sim_arr[indi-1,main_struct$compartments[[ind_indiv]]$Ec_ind]
      
    }
    if('age_struct'%in%names(main_struct))
    {
      main_struct$age_struct$norm_death_exp_d
      param_$pcrit<-  pcrit_v[indi]#min(1,param_$pcrit0*main_struct$age_struct$critical_exp_d[indi])
      param_$spdeath <- param_$spdeath0*pdeath_v[indi]/param_$pcrit0#param_$death0
      ##good variant:
      #param_$pcrit<-min(1,param_$pcrit0*main_struct$age_struct$norm_death_exp_d[indi]^0.5)
      #param_$pdeath<-min(1,param_$pdeath0*main_struct$age_struct$norm_death_exp_d[indi]^0.5)
      # if(param_$spdeath0==0)
      #{
      param_$pdeath<-pdeath_v[indi]#min(1,param_$pdeath0*main_struct$age_struct$death_exp_d[indi]/param_$pcrit)#/main_struct$age_struct$norm_critical_exp_d[indi])
      # }else{
      #   param_$pdeath <- param_$pdeath0
      #}
      
      param_$age_lab<-1
    }else{
      param_$age_lab<-0
    }
    
    
    if(indi<param_$influx_start)
    {
      param_$influx<-0
    }else{
      #if(param_$Unlock1<=indi)
      #{
      #  param_$influx<-DelayEff(t0=indi,par1=0,par2=param_$unlockinflux,t1=param_$Unlock1,delaylock=param_$delayunlock)
      #}else
     # { 
        param_$influx<-DelayEff(t0=indi,par1=param_$influx0,par2=0,t1=param_$date_treat1,delaylock=param_$delaylock)
      #}
      
    }#if(param_$Lock1<=indi)
    
    param_$r1 <- treat_arr[indi,'r1']
    param_$r2 <- treat_arr[indi,'r2']
    if(main_struct$model_opt==1)
    {
      nextav_<-ModelScholzOrig(param_,states=sim_arr[indi-1,],delt= param_$delt)
    }
    if(main_struct$model_opt==2)
    {
      nextav_<-ModelScholzOrigManyComp(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                                       delt= param_$delt)#ModelScholzIA2(param_,states=sim_arr[indi-1,],delt= param_$delt)
    }  
    if(main_struct$model_opt==3)
    {
      nextav_<-ModelScholzManyCompNew(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                                      delt= param_$delt)
    } 
    if((main_struct$model_opt==31)||(main_struct$model_opt==32))
    {
      nextav_<-ModelScholzManyCompMut(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                                      delt= param_$delt)
    } 
    if(main_struct$model_opt==4)
    {
      if(main_struct$LogSc==0)
      { 
        nextav_<-ModelSIR(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                          delt= param_$delt)
      }else{
        nextav_<-ModelSIRLog(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                             delt= param_$delt)
      }
    } 
    
    
    # nextav_<-ModelScholzOrig(param_,states=sim_arr[indi-1,],delt= param_$delt)#,u=stoch_obj$umatr[indi-1,],mod_type=main_struct$mod_type
    sim_arr[indi,]<-sim_arr[indi-1,]+param_$delt*nextav_+xtot$stox[indi-1,]#rnorm(main_struct$numvar)*
    sim_darray[indi, ]<-nextav_
  }
  ###output:
  if(main_struct$LogSc==0)
  { 
    
    return(cbind(sim_arr,SumIncreasIS(param_,sim_arr,param_$delt,main_struct$compartments[[ind_indiv]],model_opt=main_struct$model_opt,pdeath_v),treat_arr,pcrit_v,pdeath_v))
  
    }else{
    sim_arr0<-sim_arr
    sim_arr0[,main_struct$compartments[[ind_indiv]]$Sc_ind]<-exp(sim_arr0[,main_struct$compartments[[ind_indiv]]$Sc_ind])
    cbind(sim_arr0,SumIncreasIS(param_,sim_arr0,param_$delt,main_struct$compartments[[ind_indiv]],model_opt=main_struct$model_opt,pdeath_v),treat_arr)
  }
}
##################
SimCOVIDStochasticAdvTreat<- function(xtot,main_struct,ind_indiv)
{
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct) 
 # param_$speclength<- main_struct$speclength
  if(!'spdeath'%in%names(param_))
  {
    param_$spdeath<-0
  }
  param_$spdeath0 <- param_$spdeath
  if(param_$parsim_label==4)
  {
    param_$r2 <- param_$r1
  }
  if(param_$parsim_label==3)
  {
    param_$r4 <- param_$r5
  }
  if(param_$parsim_label==2)
  {
    param_$r2 <- param_$r1
    param_$r4 <- param_$r5
  }
  sim_arr<-array(dim=c(xtot$num_steps,main_struct$numvar))
  sim_arr[1,]<- xtot$init
  if(main_struct$LogSc==1)
  {
    sim_arr[1,main_struct$compartments[[ind_indiv]]$Sc_ind]<- log(xtot$init[main_struct$compartments[[ind_indiv]]$Sc_ind])
  }
  sim_darray<-array(dim=c(xtot$num_steps,main_struct$numvar))
  sim_darray[1, ]<- rep.int(x= 0, times= main_struct$numvar)
  if(main_struct$treatpardir==1)
  {
    treat_arr<-GenerateTreatStructDir(xtot,main_struct,ind_indiv,param_)
  }else if(main_struct$treatpardir==0)
  {  
    treat_arr<-GenerateTreatStruct(xtot,main_struct,ind_indiv,param_)
  }
  ###
  
  if("DateStepsPDeath"%in%names(main_struct$covariates))
  {
    death_exp_ind<-grep(pattern='death_exp',x=names(param_))
    crit_exp_ind<-grep(pattern='crit_exp',x=names(param_))
    death_exp_n <- main_struct$num_pdeath+1
    crit_exp_n <- main_struct$num_pcrit+1
    
    norm_death_exp_d <- vector(length=length(main_struct$completedates[[ind_indiv]]))
    norm_crit_exp_d <- vector(length=length(main_struct$completedates[[ind_indiv]]))
    for(ind_e in 1:death_exp_n)
    {
      if(ind_e==1)
      {
        good_ind<-(main_struct$covariates$DateStepsPDeath[[ind_indiv]][1]>=main_struct$completedates[[ind_indiv]])
        norm_death_exp_d[good_ind] <- 1
      }else{
        good_ind<-(main_struct$covariates$DateStepsPDeath[[ind_indiv]][ind_e-1]<main_struct$completedates[[ind_indiv]])
        if(ind_e <death_exp_n)
        {  
          good_ind <- good_ind&(main_struct$covariates$DateStepsPDeath[[ind_indiv]][ind_e]>=main_struct$completedates[[ind_indiv]])
        }
        norm_death_exp_d[good_ind]<-param_[[death_exp_ind[ind_e-1]]]
        #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
      }
      
    }
     
    for(ind_e in 1:crit_exp_n)
    {
      if(ind_e==1)
      {
        good_ind<-(main_struct$covariates$DateStepsPCrit[[ind_indiv]][1]>=main_struct$completedates[[ind_indiv]])
        norm_crit_exp_d[good_ind] <- 1
      }else{
        good_ind<-(main_struct$covariates$DateStepsPCrit[[ind_indiv]][ind_e-1]<main_struct$completedates[[ind_indiv]])
        if(ind_e <crit_exp_n)
        {  
          good_ind <- good_ind&((main_struct$covariates$DateStepsPCrit[[ind_indiv]][ind_e]>=main_struct$completedates[[ind_indiv]]))
        }
        norm_crit_exp_d[good_ind]<-param_[[crit_exp_ind[ind_e-1]]]
        #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
      }
      
    }
    
  }else{
    
    death_exp_ind<- main_struct$death_exp_ind# grep(pattern='death_exp',x=names(param_))
    crit_exp_ind<-main_struct$crit_exp_ind# grep(pattern='crit_exp',x=names(param_))
    death_exp_n <- main_struct$death_exp_n# length(death_exp_ind)+2
    crit_exp_n <- main_struct$crit_exp_n# length(crit_exp_ind)+2
    #norm_death_exp_d<-main_struct$age_struct$norm_death_exp_d
    #norm_crit_exp_d<-main_struct$age_struct$norm_death_exp_d
    norm_death_exp_d <- vector(length = xtot$num_steps)
    norm_crit_exp_d <- vector(length = xtot$num_steps)
    intervals_death <- vector(length = death_exp_n)
    intervals_crit <- vector(length = crit_exp_n)
    intervals_death[1] <- 1
    intervals_crit[1] <- 1
    if("death_int1"%in%names(param_))
    {
      help_str <- KICreateIntervals(param_,xtot,main_struct)
      intervals_crit <- help_str$intervals_crit
      intervals_death <- help_str$intervals_death
    }else{
      label_KI <- 0
      num_rep_death_exp<-floor(main_struct$speclength/death_exp_n)
      num_rep_crit_exp<-floor(main_struct$speclength/crit_exp_n)
      for(ind_e in 1:(death_exp_n-2))
      {  
        intervals_death[ind_e+1]<- ind_e*num_rep_death_exp
      } 
      for(ind_e in 1:(crit_exp_n-2))
      {  
        intervals_crit[ind_e+1]<- ind_e*num_rep_crit_exp
      }
      intervals_death[death_exp_n]<- main_struct$speclength
      intervals_crit[death_exp_n]<- main_struct$speclength
    }
    
    if(death_exp_n>0)
    {
      #num_rep_death_exp<-floor(main_struct$speclength/death_exp_n)#floor(length(norm_death_exp_d)/death_exp_n)
      for(ind_e in 1:(death_exp_n-1))
      {
        if(ind_e==1)
        {
          norm_death_exp_d[intervals_death[ind_e]:intervals_death[ind_e+1]] <- 1
          #norm_death_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)] <- 1
        }else{
          norm_death_exp_d[intervals_death[ind_e]:intervals_death[ind_e+1]] <- param_[[death_exp_ind[ind_e-1]]]
          #norm_death_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[death_exp_ind[ind_e-1]]]
          #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
        }
        
      }
      #norm_death_exp_d[(death_exp_n*num_rep_death_exp):length(main_struct$age_struct$norm_death_exp_d)]<-param_[[death_exp_ind[death_exp_n-1]]]
      #norm_crit_exp_d[(death_exp_n*num_rep_death_exp):length(main_struct$age_struct$norm_death_exp_d)]<-param_[[norm_crit_exp_d[death_exp_n]]]
    }
    #
    if(crit_exp_n>0)
    {
      #num_rep_crit_exp<-floor(main_struct$speclength/crit_exp_n)#floor(length(norm_crit_exp_d)/crit_exp_n)
      
      for(ind_e in 1:(crit_exp_n-1))
      {
        if(ind_e==1)
        {
          norm_crit_exp_d[intervals_crit[ind_e]:intervals_crit[ind_e+1]] <- 1
          #norm_crit_exp_d[((ind_e-1)*num_rep_crit_exp+1):(ind_e*num_rep_crit_exp)] <- 1
        }else{
          norm_crit_exp_d[intervals_crit[ind_e]:intervals_crit[ind_e+1]] <- param_[[crit_exp_ind[ind_e-1]]]
          #norm_crit_exp_d[((ind_e-1)*num_rep_crit_exp+1):(ind_e*num_rep_crit_exp)]<-param_[[crit_exp_ind[ind_e-1]]]
        }
        
        #norm_crit_exp_d[((ind_e-1)*num_rep_death_exp+1):(ind_e*num_rep_death_exp)]<-param_[[crit_exp_ind[ind_e]]]
      }
      #norm_death_exp_d[(death_exp_n*num_rep_death_exp):length(main_struct$age_struct$norm_death_exp_d)]<-param_[[death_exp_ind[death_exp_n]]]
      #norm_crit_exp_d[(crit_exp_n*num_rep_crit_exp):length(norm_crit_exp_d)]<-param_[[crit_exp_ind[crit_exp_n-1]]]
    }
  }
  
  
  param_$influx0<-param_$influx
  param_$pcrit0 <- param_$pcrit
  param_$pdeath0 <- param_$pdeath
  #pcrit_v<-pmin(1,param_$pcrit0*main_struct$age_struct$critical_exp_d)
  pcrit_v <- pmin(0.95,param_$pcrit0*norm_crit_exp_d)#^param_$powcrit0.5
  pdeath_v <- pmin(0.66,param_$pdeath0*norm_death_exp_d)#^(1-param_$powcrit) 0.5
  if(length(pcrit_v)<xtot$num_steps)
  {
    pcrit_v[(length(pcrit_v)+1):dim(sim_arr)[1]] <- pcrit_v[(length(pcrit_v))]
    pdeath_v[(length(pdeath_v)+1):dim(sim_arr)[1]] <- pdeath_v[(length(pdeath_v))]
  }
  #if(max(param_$pcrit0*main_struct$age_struct$critical_exp_d)>1)
  #{
 #   pcrit_v<-param_$pcrit0*main_struct$age_struct$critical_exp_d/max(pcrit_v)
 # } 
  #pdeath_v<- pmin(1,param_$pdeath0*main_struct$age_struct$death_exp_d)#pmin(1,param_$pdeath0*main_struct$age_struct$death_exp_d/pcrit_v)#param_$pdeath0*main_struct$age_struct$death_exp_d/pcrit_v
 
  # if(max(pdeath_v)>1)
 # {
 #   pdeath_v <- pdeath_v/max(pdeath_v)
 # }
  for(indi in 2:xtot$num_steps){
    if('age_struct'%in%names(main_struct))
    {
      main_struct$age_struct$norm_death_exp_d
      param_$pcrit<-  pcrit_v[indi]#min(1,param_$pcrit0*main_struct$age_struct$critical_exp_d[indi])
      param_$spdeath <- param_$spdeath0*pdeath_v[indi]/param_$pcrit0#param_$death0
      ##good variant:
      #param_$pcrit<-min(1,param_$pcrit0*main_struct$age_struct$norm_death_exp_d[indi]^0.5)
      #param_$pdeath<-min(1,param_$pdeath0*main_struct$age_struct$norm_death_exp_d[indi]^0.5)
     # if(param_$spdeath0==0)
      #{
        param_$pdeath<-pdeath_v[indi]#min(1,param_$pdeath0*main_struct$age_struct$death_exp_d[indi]/param_$pcrit)#/main_struct$age_struct$norm_critical_exp_d[indi])
     # }else{
     #   param_$pdeath <- param_$pdeath0
      #}
      
     param_$age_lab<-1
    }else{
      param_$age_lab<-0
    }
   
    
    if(indi<param_$influx_start)
    {
      param_$influx<-0
    }else{
      if(param_$Unlock1<=indi)
      {
        param_$influx<-DelayEff(t0=indi,par1=0,par2=param_$unlockinflux,t1=param_$Unlock1,delaylock=param_$delayunlock)
      }else
      { 
        param_$influx<-DelayEff(t0=indi,par1=param_$influx0,par2=0,t1=param_$Lock1,delaylock=param_$delaylock)
      }
      
    }#if(param_$Lock1<=indi)
    
    param_$r1 <- treat_arr[indi,'r1']
    param_$r2 <- treat_arr[indi,'r2']
    if(main_struct$model_opt==1)
    {
      nextav_<-ModelScholzOrig(param_,states=sim_arr[indi-1,],delt= param_$delt)
    }
    if(main_struct$model_opt==2)
    {
      nextav_<-ModelScholzOrigManyComp(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                                       delt= param_$delt)#ModelScholzIA2(param_,states=sim_arr[indi-1,],delt= param_$delt)
    }  
    if(main_struct$model_opt==3)
    {
      nextav_<-ModelScholzManyCompNew(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                                      delt= param_$delt)
    } 
    if((main_struct$model_opt==31)||(main_struct$model_opt==32))
    {
      nextav_<-ModelScholzManyCompMut(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                                      delt= param_$delt)
    } 
    if(main_struct$model_opt==4)
    {
      if(main_struct$LogSc==0)
      { 
        nextav_<-ModelSIR(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                          delt= param_$delt)
      }else{
        nextav_<-ModelSIRLog(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                             delt= param_$delt)
      }
    } 
    
    
    # nextav_<-ModelScholzOrig(param_,states=sim_arr[indi-1,],delt= param_$delt)#,u=stoch_obj$umatr[indi-1,],mod_type=main_struct$mod_type
    sim_arr[indi,]<-sim_arr[indi-1,]+param_$delt*nextav_+xtot$stox[indi-1,]#rnorm(main_struct$numvar)*
    sim_darray[indi, ]<-nextav_
  }
  ###output:
  if(main_struct$LogSc==0)
  { 
    
    cbind(sim_arr,SumIncreasIS(param_,sim_arr,param_$delt,main_struct$compartments[[ind_indiv]],model_opt=main_struct$model_opt,pdeath_v),treat_arr,pcrit_v,pdeath_v)
  }else{
    sim_arr0<-sim_arr
    sim_arr0[,main_struct$compartments[[ind_indiv]]$Sc_ind]<-exp(sim_arr0[,main_struct$compartments[[ind_indiv]]$Sc_ind])
    cbind(sim_arr0,SumIncreasIS(param_,sim_arr0,param_$delt,main_struct$compartments[[ind_indiv]],model_opt=main_struct$model_opt,pdeath_v),treat_arr)
  }
}
#################
#################
GenerateTreatStruct<- function(xtot,main_struct,ind_indiv,param_)
{
  #param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)
  treat_arr <- array(dim=c(xtot$num_steps,2))
  colnames(treat_arr)<-c('r1','r2')
  
  treat_arr[,'r1']<-rep.int(x=param_$r1, times = xtot$num_step)
  treat_arr[,'r2']<-rep.int(x=param_$r2, times = xtot$num_step)
  treat_arr <- UpdateTreatStruct(treat_arr,num_steps=xtot$num_steps,param_,main_struct)
  
  return(treat_arr)
}
UpdateTreatStruct<- function(treat_arr,num_steps,param_,main_struct)
{
  all_unlocks <- main_struct$all_unlocks#names(param_)[grep(pattern = 'Unlock', x = names(param_))]
  all_locks <- main_struct$all_locks#names(param_)[grep(pattern = 'Lock', x = names(param_))]
  #all_locksr2 <- names(param_)[grep(pattern = 'Lock', x = names(param_))]
  
  num_locks <- main_struct$num_locks#length(all_locks)
  num_unlocks <- main_struct$num_unlocks#length(all_unlocks)
  
  for(indi in 2:num_steps){
    for(ind_locks in 1:num_locks)
    {
      treat_arr[indi,'r1']<- treat_arr[indi,'r1']*DelayEff(t0=indi,par1=1,par2=param_[[paste('blockeff',ind_locks,sep='')]],
                            t1=param_[[paste('Lock',ind_locks,sep='')]],
                            delaylock=param_$delaylock)
    }
    for(ind_unlocks in 1:num_unlocks)
    {
      treat_arr[indi,'r1']<- treat_arr[indi,'r1']*DelayEff(t0=indi,par1=1,par2=param_[[paste('unlockeff',ind_unlocks,sep='')]],
                            t1= param_[[paste('Unlock',ind_unlocks,sep='')]],
                            delaylock=param_$delayunlock)
    }
    if(param_$parsim_treat==0){
      treat_arr[indi,'r2'] <- treat_arr[indi,'r2']*DelayEff(t0=indi,par1=1,par2=param_$blockeff_r2,t1=param_$Lock2,
                                                            delaylock=param_$delaylock)
    }
    
  }
  if(param_$parsim_treat==1)
  {
    treat_arr[,'r2'] <- param_$blockeff_r2*treat_arr[,'r1']
  }
  return(treat_arr)
}
###
GenerateTreatStructDir<-function(xtot,main_struct,ind_indiv,param_)#SimulIntegralModelShort
{
  num_steps <- xtot$num_steps
  treat_arr <- array(dim=c(xtot$num_steps,2))
  colnames(treat_arr)<-c('r1','r2')
  
  #treat_arr[,'r1']<-rep.int(x=1, times = xtot$num_step)#param_$r1
  #treat_arr[,'r2']<-rep.int(x=1, times = xtot$num_step)#param_$r2
  treat_arr[1,'r1']<-1
  treat_arr[1,'r2']<-1
  compl_parn <- c(names(param_),'num_locks','num_unlocks')#c(names(param_)[names(param_)!='DateStart'],'num_locks','num_unlocks')
  compl_parl <- length(compl_parn)
  compl_parv <- vector(length = compl_parl)
  indc<-0
  for(indi in 1: length(param_))
  {
    if(names(param_)[indi]!='DateStart')
    {
      #indc <- indc+1
      compl_parv[indi] <- param_[[indi]]
    }else{
      compl_parv[indi] <- 0#
    }
    
  }
  # blockeff_vec <- compl_parv[grep(pattern='blockeff',x=compl_parn)]
  
  compl_parv[compl_parl-1]<- main_struct$num_locks
  compl_parv[compl_parl]<- main_struct$num_unlocks
  blockeff_vec <- compl_parv[main_struct$blockeff_ind]
  unlockeff_vec <- compl_parv[main_struct$unlockeff_ind]
  Lock_vec <- compl_parv[main_struct$all_locks_ind]#-1 for c++!!!
  Unlock_vec <- compl_parv[main_struct$all_unlocks_ind]#-1 for c++!!!
  prev_val<-1
  ind_lock <-1
  ind_unlock <- 1
  num_locks <- main_struct$num_locks
  num_unlocks <- main_struct$num_unlocks
  for(indi in 2:num_steps){
    if((indi>=Lock_vec[ind_lock])&&(indi <= (Lock_vec[ind_lock]+param_$delaylock)))
    {
      
      treat_arr[indi,'r1'] <- DelayEff(t0=indi,par1=prev_val,par2=blockeff_vec[ind_lock],
                                       t1=Lock_vec[ind_lock],
                                       delaylock=param_$delaylock)
      if(indi == (Lock_vec[ind_lock]+param_$delaylock))
      {
        prev_val<-blockeff_vec[ind_lock]
        if(ind_lock<num_locks)
        {
          ind_lock <- ind_lock+1
        }
      }
    }
    else if((indi>=Unlock_vec[ind_unlock])&&(indi <= (Unlock_vec[ind_unlock]+param_$delayunlock)))
    {
      
      treat_arr[indi,'r1'] <-DelayEff(t0=indi,par1=prev_val,par2=unlockeff_vec[ind_unlock],
                                      t1=Unlock_vec[ind_unlock],
                                      delaylock=param_$delayunlock)
      if(indi == (Unlock_vec[ind_unlock]+param_$delayunlock))
      {
        prev_val<-unlockeff_vec[ind_unlock]
        if(ind_unlock<num_unlocks)
        {
          ind_unlock <- ind_unlock+1
        }
      }
    }else if(indi > 1)
    {
      treat_arr[indi,'r1'] <- treat_arr[indi-1,'r1']
    }
    
  }
  treat_arr[,'r1']<- param_$r1*treat_arr[,'r1']
  treat_arr[,'r2'] <- param_$blockeff_r2*treat_arr[,'r1']
  return(treat_arr)
}
###
GenerateTreatStructKIDir<-function(xtot,main_struct,ind_indiv,param_)#SimulIntegralModelShort
{
  num_steps <- xtot$num_steps
  treat_arr <- array(dim=c(xtot$num_steps,2))
  colnames(treat_arr)<-c('r1','r2')
  
  #treat_arr[,'r1']<-rep.int(x=1, times = xtot$num_step)#param_$r1
  #treat_arr[,'r2']<-rep.int(x=1, times = xtot$num_step)#param_$r2
  treat_arr[1,'r1']<-1
  treat_arr[1,'r2']<-1
  compl_parn <- c(names(param_),'num_treats')#c(names(param_)[names(param_)!='DateStart'],'num_locks','num_unlocks')
  compl_parl <- length(compl_parn)
  compl_parv <- vector(length = compl_parl)
  indc<-0
  non_num_names <- c('DateStart','age_spec_param_names','age_spec_phi','age_spec_param','age_spec_dim','age_spec_id')
  for(indi in 1: length(param_))
  {
    
    if(!(names(param_)[indi]%in%non_num_names))#!='DateStart')
    {
      #indc <- indc+1
      compl_parv[indi] <- param_[[indi]]
      if(length(param_[[indi]])>1)
      {  
        print(c(indi,names(param_)[indi],length(param_[[indi]])))
      }
    }else{
      compl_parv[indi] <- 0#
    }
    
  }
  # blockeff_vec <- compl_parv[grep(pattern='blockeff',x=compl_parn)]
  
  #compl_parv[compl_parl-1]<- main_struct$num_locks
  compl_parv[compl_parl]<- main_struct$treat_par_num
  treat_vec <- compl_parv[main_struct$treat_par_ind]
  
  LockUnlock_vec <- round(compl_parv[main_struct$treat_date_ind])#-1 for c++!!!
  prev_val<-1
  ind_treat <-1
  
  num_treats <- main_struct$treat_par_num
  #num_unlocks <- main_struct$treat_par_num
  for(indi in 2:num_steps){
    if((indi>=LockUnlock_vec[ind_treat])&&(indi <= (LockUnlock_vec[ind_treat]+param_$delaylock)))
    {
      
      treat_arr[indi,'r1'] <- DelayEff(t0=indi,par1=prev_val,par2=treat_vec[ind_treat],
                                       t1=LockUnlock_vec[ind_treat],
                                       delaylock=param_$delaylock)
      if(indi == (LockUnlock_vec[ind_treat]+param_$delaylock))
      {
        prev_val<-treat_vec[ind_treat]
        if(ind_treat<num_treats)
        {
          ind_treat <- ind_treat+1
        }
      }
    }
    else if(indi > 1)
    {
      treat_arr[indi,'r1'] <- treat_arr[indi-1,'r1']
    }
    
  }
  treat_arr[,'r1']<- param_$r1*treat_arr[,'r1']
  treat_arr[,'r2'] <- param_$blockeff_r2*treat_arr[,'r1']
  return(treat_arr)
}
###
KICreateIntervals <- function(param_,xtot,main_struct)
{ 
  death_exp_ind<- main_struct$death_exp_ind#grep(pattern='death_exp',x=names(param_))
  crit_exp_ind<- main_struct$crit_exp_ind#grep(pattern='crit_exp',x=names(param_))
  death_exp_n <- main_struct$death_exp_n# length(death_exp_ind)+2
  crit_exp_n <- main_struct$crit_exp_n# length(crit_exp_ind)+2
  death_int_ind<- main_struct$death_int_ind#grep(pattern='death_int',x=names(param_))
  crit_int_ind<- main_struct$crit_int_ind#grep(pattern='crit_int',x=names(param_))
  intervals_death <- vector(length = death_exp_n)
  intervals_crit <- vector(length = crit_exp_n)
  intervals_death[1] <- 1
  intervals_crit[1] <- 1
  if('dates_treatments'%in%names(main_struct))
  {
    for(ind_e in 1:(death_exp_n-2))
    {  
      intervals_death[ind_e+1]<- round(param_[[paste('death_date',ind_e,sep='')]])
    } 
    for(ind_e in 1:(crit_exp_n-2))
    {  
      intervals_crit[ind_e+1]<- round(param_[[paste('crit_date',ind_e,sep='')]])
    }
  }
  else
  {
    for(ind_e in 1:(death_exp_n-2))
    {  
      intervals_death[ind_e+1]<- intervals_death[ind_e]+round(param_[[paste('death_int',ind_e,sep='')]])
    } 
    for(ind_e in 1:(crit_exp_n-2))
    {  
      intervals_crit[ind_e+1]<- intervals_crit[ind_e]+round(param_[[paste('crit_int',ind_e,sep='')]])
    }
  }
  #delta of the last interval:
  if('deltapcritpdeath'%in%names(param_))
  { 
    intervals_crit[crit_exp_n-1]<- param_$speclength- param_$deltapcritpdeath+1#instead xtot$num_steps-...
    intervals_death[death_exp_n-1]<- param_$speclength- param_$deltapcritpdeath+1#instead xtot$num_steps-...
  }
  intervals_death[death_exp_n]<- xtot$num_steps
  intervals_crit[crit_exp_n]<- xtot$num_steps
  intervals_crit<-pmin(intervals_crit,xtot$num_steps)
  intervals_death <- pmin(intervals_death,xtot$num_steps )
  intervals_crit[length(intervals_crit)]<-xtot$num_steps
  intervals_death[length(intervals_death)]<-xtot$num_steps
  
  out_ <- list()
  out_$intervals_crit <- intervals_crit
  out_$intervals_death <- intervals_death
  return(out_)
}
DirFromIndirTreatParams <- function(x,treat_arr,param_,main_struct)
{  
  num_locks <- main_struct$num_locks
  num_unlocks <- main_struct$num_unlocks
  num_steps <- dim(treat_arr)[1]
  dir_lock_vect<-vector(length = num_locks)
  dir_unlock_vect<-vector(length = num_unlocks)
  for(indi in 1:num_steps)
  {  
    for(ind_locks in 1:num_locks)
    {
      if(indi==(param_[[paste('Lock',ind_locks,sep='')]]+param_$delaylock+1))
      {
        dir_lock_vect[ind_locks] <- treat_arr[indi,'vr1']/treat_arr[1,'vr1']
      }
      
    }
    for(ind_unlocks in 1:num_unlocks)
    {
      if(indi==(param_[[paste('Unlock',ind_unlocks,sep='')]]+param_$delayunlock+1))
      {
        dir_unlock_vect[ind_unlocks] <- treat_arr[indi,'vr1']/treat_arr[1,'vr1']
      }
    }
  }
  out_ <- list()
  out_$dir_lock_vect <- dir_lock_vect
  out_$dir_unlock_vect <- dir_unlock_vect
  out_$x <-x
  
  return(out_)
}
AddControle<- function(path_d,cont_)
{
  temp_indiv<-read.csv(paste(path_d,'/indivpar.csv',sep = ''), sep = ";",dec = ',')
  temp_date<-read.csv(paste(path1,'/DatesUncertainty.csv',sep = ''), sep = ";",dec = ',')
  
}