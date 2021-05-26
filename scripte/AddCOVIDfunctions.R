InfluxesCalc <- function(param_,states,compartments,delt)#u, ,mod_type DifModelFribergLog
{
  #states are logaithm of states for the cells and chemotherapy effects:
  #u is a vector of treatments: cyclophosphamide, doxorub, etoposide, platelets, filgrastim
  #num_transit
  #if(mod_type==2){
  #   num_transit<-1
  # }
  
  #num_sys_comp<- num_transit+6#
  #dScdt <- vector(length = compartments$Scnum)
  dIAcdt <- vector(length = compartments$IAcnum)
  dIScdt <- vector(length = compartments$IScnum)
  dCcdt <- vector(length = compartments$Ccnum)
  #dDcdt <- vector(length = compartments$Dcnum)
  # dRc1dt <- vector(length = compartments$Rc1num)
  # dRc2dt <- vector(length = compartments$Rc2num)
  
  Sc<- max(states[compartments$Sc_ind],0)#exp(states[num_sys_comp+1])+
  IAc <- states[compartments$IAc_ind]#exp(states[num_sys_comp+2])+
  ISc <- states[compartments$ISc_ind]#exp(states[num_sys_comp+3])+
  
  Cc <- states[compartments$Cc_ind]#exp(states[num_sys_comp+5])+ basic level is neglidgible!!! for the thrombo model!!!
  Dc <- states[compartments$Dc_ind]
  Rc1 <- states[compartments$Rc1_ind]
  Rc2 <- states[compartments$Rc2_ind]#exp(states[num_sys_comp+4])+
  # if('r81'%in%names(param_))
  # {
  #   add_death<- param_$r81*Cc/(param_$capacity+Cc)
  # }else{
  add_death<-0
  #}
  
  dScdt <- -param_$influx -param_$r1*(Sc/param_$Sc0)*(IAc[2]+IAc[3]) - param_$r2*(Sc/param_$Sc0)*(ISc[1]+ISc[2]+ISc[3])
  dIAcdt[1] <- param_$r1*param_$r1*(Sc/param_$Sc0)*(IAc[2]+IAc[3]) + param_$r2*param_$r2*(Sc/param_$Sc0)*(ISc[1]+ISc[2]+ISc[3]) -param_$r4*IAc[1]
  dIAcdt[2] <- param_$influx + param_$r4*IAc[1]-param_$r3*IAc[2]#*param_$psymp*
  dIAcdt[3] <- param_$r3*(1-param_$psymp)*IAc[2]-param_$rd4*IAc[3]#param_$influx + 
  
  dIScdt[1] <- param_$r3*param_$psymp*IAc[2]-param_$r5*ISc[1]
  dIScdt[2] <- param_$r5*ISc[1]-(param_$r5+param_$r6)*ISc[2]
  dIScdt[3] <- param_$r5*ISc[2]-(param_$r5+param_$r6)*ISc[3]
  dCcdt[1] <- param_$r6*(ISc[2]+ISc[3])-(param_$r7 + param_$r8+add_death)*Cc[1]
  dCcdt[2] <- param_$r7*Cc[1]-(param_$r7 + param_$r8+add_death)*Cc[2]
  dCcdt[3] <- param_$r7*Cc[2]-(param_$r7 + param_$r8+add_death)*Cc[3]
  dDcdt <- param_$r8*(Cc[1]+Cc[2]+Cc[3])+add_death*sum(Cc)
  dRc1dt <- param_$r5*ISc[3] + param_$r7*Cc[3]
  dRc2dt <-   param_$rd4*IAc[3]
  out_<-vector(length=compartments$NumofVariables)
  out_[compartments$Sc_ind]<- dScdt
  out_[compartments$IAc_ind] <- dIAcdt
  out_[compartments$ISc_ind] <- dIScdt
  
  out_[compartments$Cc_ind] <- dCcdt
  out_[compartments$Dc_ind] <- dDcdt
  out_[compartments$Rc1_ind] <- dRc1dt
  out_[compartments$Rc2_ind] <- dRc2dt
  out_<-delt*out_#c(dScdt,dIAcdt,dIScdt, dCcdt,dDcdt,dRc1dt,dRc2dt)###,dedrug)#dycyclo,dydoxo,dyetop,dpl
  out_
 
  
  
  out_<-vector(length=compartments$NumofVariables)
  out_[compartments$Sc_ind]<- dScdt
  out_[compartments$IAc_ind] <- dIAcdt
  out_[compartments$ISc_ind] <- dIScdt
  
  out_[compartments$Cc_ind] <- dCcdt
  out_[compartments$Dc_ind] <- dDcdt
  out_[compartments$Rc1_ind] <- dRc1dt
  out_[compartments$Rc2_ind] <- dRc2dt
  out_<-delt*out_#c(dScdt,dIAcdt,dIScdt, dCcdt,dDcdt,dRc1dt,dRc2dt)###,dedrug)#dycyclo,dydoxo,dyetop,dpl
  #if(param_$add_out)
  #{
    InflIAc1r1 <- param_$r1*(Sc/param_$Sc0)*(IAc[2]+IAc[3])
    InflIAcd1r2 <- param_$r2*(Sc/param_$Sc0)*(ISc[1]+ISc[2]+ISc[3])
    OutlIAcd1<-param_$r4*IAc[1]
    add_out<-delt*c(InflIAc1r1,InflIAcd1r2,param_$influx,OutlIAcd1,param_$r3*(1-param_$psymp)*IAc[2])
    out_list<-list()
    out_list$diff <- out_
    out_list$add_out <-add_out
  #}
    out_list
}


AddSimCOVIDStochastic<- function(xtot,main_struct,ind_indiv)
{
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct) 
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
  sim_arr <- array(dim=c(xtot$num_steps,main_struct$numvar))
  add_inf <- array(dim=c(xtot$num_steps,5))
  add_inf[1,] <- rep.int(x=0, times=5)
  sim_arr[1,] <- xtot$init
  sim_darray<-array(dim=c(xtot$num_steps,main_struct$numvar))
  sim_darray[1, ]<- rep.int(x= 0, times= main_struct$numvar)
  param_$r10<-param_$r1*param_$blockeff
  param_$r20<-param_$r2*param_$blockeff2# param_$r20<-param_$r2*param_$blockeff
  
  param_$r11<-param_$r10*param_$blockeff2
  param_$r21<-param_$r20*param_$blockeff2
  param_$influx0<-param_$influx
  for(indi in 2:xtot$num_steps){
    #if(indi>=48)
    #{
    #  print(indi)
    # } 
    if(indi <=max(unique(main_struct$measurements[[ind_indiv]]$data[1,])))
    {  
      param_$blockeff30<-param_$blockeff3*indi
    }else
    {
      param_$blockeff30<-max(unique(main_struct$measurements[[ind_indiv]]$data[1,]))
    }
    if(indi<param_$influx_start)
    {
      param_$influx<-0
    }else{
      param_$influx<-DelayEff(t0=indi,par1=param_$influx0,par2=0,t1=param_$Lock1,delaylock=param_$delaylock)
    }#if(param_$Lock1<=indi)
    #{
    #  param_$influx<-0
    #}
    
    #if(param_$Lock2<=indi)
    #{
    #  param_$r1 <- param_$r10
    #  param_$r2 <- param_$r20
    # }
    #!!!!One shutdown!
    # if(param_$Lock3<=indi){#!!!!!
    #   param_$r1<-DelayEff(t0=indi,par1=param_$r10,par2=param_$r11,t1=param_$Lock3,delaylock=param_$delaylock)
    #   param_$r2<-DelayEff(t0=indi,par1=param_$r20,par2=param_$r21,t1=param_$Lock3,delaylock=param_$delaylock)
    # }else{
    param_$r1<-DelayEff(t0=indi,par1=param_$r1,par2=param_$r10,t1=param_$Lock2,delaylock=param_$delaylock)
    param_$r2<-DelayEff(t0=indi,par1=param_$r2,par2=param_$r20,t1=param_$Lock2,delaylock=param_$delaylock)
    
    #}
    
    # if(param_$Lock3<=indi)
    #{
    #  param_$r1 <- param_$r11
    #  param_$r2 <- param_$r21
    #}
    #param_$r1<- (1-param_$block)*param_$r1+param_$r1*param_$block*param_$blockeff
    #param_$r2<- (1-param_$block)*param_$r2+param_$r2*param_$block*param_$blockeff
    #if(main_struct$model_opt==1)
    #{
    #  nextav_<-ModelScholzOrig(param_,states=sim_arr[indi-1,],delt= param_$delt)
    #}
    if(main_struct$model_opt==2)
    {
      #nextav_<-ModelScholzOrigManyComp(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
       #                                delt= param_$delt)#ModelScholzIA2(param_,states=sim_arr[indi-1,],delt= param_$delt)
      out_list<-InfluxesCalc(param_,states=sim_arr[indi-1,],compartments=main_struct$compartments[[ind_indiv]],
                   delt= param_$delt)
      nextav_ <- out_list$diff
      add_out <- out_list$add_out 
      add_inf[indi,] <- add_out
      
      }  
    #if(main_struct$model_opt==3)
    #{
    #  nextav_<-ModelScholzHolger(param_,states=sim_arr[indi-1,],delt= param_$delt)
    #} 
    
    # nextav_<-ModelScholzOrig(param_,states=sim_arr[indi-1,],delt= param_$delt)#,u=stoch_obj$umatr[indi-1,],mod_type=main_struct$mod_type
    sim_arr[indi,]<-sim_arr[indi-1,]+nextav_+xtot$stox[indi-1,]#rnorm(main_struct$numvar)*
    sim_darray[indi, ]<-nextav_
  }
  ###output:
  list_out <-list()
  list_out$sim_res <- cbind(sim_arr,SumIncreasIS(param_,sim_arr,param_$delt,main_struct$compartments[[ind_indiv]]))
  list_out$sim_res <- NameOutput(list_out$sim_res,main_struct)
  colnames(add_inf)<-c('InfluxIA1fromr1','InfluxIA1fromr2','InfluxIA2OtherLands','InfluxIA2fromIA1','InfluxIA3fromIA2')
    #InflIAc1r1,InflIAcd1r2,OutlIAcd1,param_$r4*IAc[1],param_$r3*(1-param_$psymp)*IAc[2]
  list_out$add_inf <- add_inf
  list_out$sim_darray <- sim_darray
  colnames(list_out$sim_darray)<-colnames(list_out$sim_res)[1:main_struct$numvar]
  ##output:
  list_out
}