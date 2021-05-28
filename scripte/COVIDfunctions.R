ModelSIRLog<- function(param_,states,compartments,delt)
{
  LSc<- states[compartments$Sc_ind]#max(states[compartments$Sc_ind],0)
  Sc<- exp(LSc)
  Ic <- states[compartments$Ic_ind]
  Rc <- states[compartments$Ic_ind]
  Nc <- Ic+Sc+Rc
  #dScdt <- param_$v*Nc-param_$beta*Sc*Ic/max(Nc,1)-param_$mu*Sc
  dLScdt <- -param_$beta*Ic/max(Nc,1)-param_$mu  #param_$v*Nc-no
  dIcdt <- param_$beta*Sc*Ic/max(Nc,1)- param_$gam*Ic -param_$mu*Ic
  dRcdt<- param_$gam*Ic-param_$mu*Rc
  out_<-vector(length=compartments$NumofVariables)
  out_[compartments$Sc_ind]<- dLScdt
  out_[compartments$Ic_ind] <- dIcdt
  out_[compartments$Rc_ind] <- dRcdt
  
  out_#delt*
}
ModelSIR<- function(param_,states,compartments,delt)
{
  Sc<- states[compartments$Sc_ind]#max(states[compartments$Sc_ind],0)
  Ic <- states[compartments$Ic_ind]
  Rc <- states[compartments$Rc_ind]
  Nc <- Ic+Sc+Rc
  dScdt <- param_$v*Nc-param_$beta*Sc*Ic/Nc-param_$mu*Sc#max(Nc,1
  dIcdt <- param_$beta*Sc*Ic/Nc- param_$gam*Ic -param_$mu*Ic
  dRcdt<- param_$gam*Ic-param_$mu*Rc
  out_<-vector(length=compartments$NumofVariables)
  out_[compartments$Sc_ind]<- dScdt
  out_[compartments$Ic_ind] <- dIcdt
  out_[compartments$Rc_ind] <- dRcdt
  
  out_#delt*
}
ModelScholzOrigManyComp <- function(param_,states,compartments,delt)#u, ,mod_type DifModelFribergLog
{
  
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
  dIAcdt[1] <- p_infect(param_$r1,param_)*param_$r1*(Sc/param_$Sc0)*(IAc[2]+IAc[3]) + p_infect(param_$r2,param_)*param_$r2*(Sc/param_$Sc0)*(ISc[1]+ISc[2]+ISc[3]) -param_$r4*IAc[1]
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
  out_<-out_#delt*  #c(dScdt,dIAcdt,dIScdt, dCcdt,dDcdt,dRc1dt,dRc2dt)###,dedrug)#dycyclo,dydoxo,dyetop,dpl
  out_
}####
ModelScholzManyCompNew<- function(param_,states,compartments,delt)#u, ,mod_type DifModelFribergLog
{
 
  Sc<- max(states[compartments$Sc_ind],0)#exp(states[num_sys_comp+1])+
  Ec<- max(states[compartments$Ec_ind],0)#
  IAc <- states[compartments$IAc_ind]#exp(states[num_sys_comp+2])+
  ISc <- states[compartments$ISc_ind]#exp(states[num_sys_comp+3])+
  
  Cc <- states[compartments$Cc_ind]#exp(states[num_sys_comp+5])+ basic level is neglidgible!!! for the thrombo model!!!
  Dc <- states[compartments$Dc_ind]
  Rc <- states[compartments$Rc_ind]#now Rc1,RC2
  StCare <- states[compartments$StCare_ind]#exp(states[num_sys_comp+4])+
  
  dIAcdt <- vector(length = compartments$IAcnum)
  dIScdt <- vector(length = compartments$IScnum) 
  dCcdt <- vector(length = compartments$Ccnum)
  #X_p(r_1 = ,r_2 = ,p_1 = )
  pIA1IS1 <-X_p(r_1 = param_$rb4,r_2 = param_$r4,p_1 = param_$psymp )
  if(param_$age_lab==0)
  {  
    pIS1Cc1 <-X_p(r_1 = param_$r6,r_2 = param_$r5,p_1 = param_$pcrit )
    pCc1Dc <-X_p(r_1 = param_$r8,r_2 = param_$r7,p_1 = param_$pdeath)
  }else{
    pIS1Cc1 <- param_$pcrit
    pCc1Dc <- param_$pdeath
  }
  
  dScdt <- -param_$influx -param_$r1*(Sc/param_$Sc0)*(IAc[1]+IAc[2]+IAc[3]) - param_$r2*(Sc/param_$Sc0)*(ISc[1]+ISc[2]+ISc[3])
  dEcdt <- param_$r1*(Sc/param_$Sc0)*(IAc[1]+IAc[2]+IAc[3]) + param_$r2*(Sc/param_$Sc0)*(ISc[1]+ISc[2]+ISc[3])-param_$r3*Ec
  dIAcdt[1] <- param_$influx + param_$r3*Ec - (pIA1IS1*param_$rb4 + (1-pIA1IS1)*param_$r4)*IAc[1]
  dIAcdt[2] <- (1-pIA1IS1)*param_$r4*IAc[1]-param_$r4*IAc[2]#*param_$psymp*
  dIAcdt[3] <- param_$r4*IAc[2]-param_$r4*IAc[3]#param_$influx + 
  dIScdt[1] <- pIA1IS1*param_$rb4*IAc[1]-(pIS1Cc1*param_$r6 + (1-pIS1Cc1)*param_$r5)*ISc[1]
  dIScdt[2] <- (1-pIS1Cc1)*param_$r5*ISc[1]-(1-param_$spdeath)*param_$r5*ISc[2]-param_$r8*param_$spdeath*ISc[2]
  dIScdt[3] <- (1-param_$spdeath)*param_$r5*ISc[2]-param_$r5*ISc[3]
  dCcdt[1] <- pIS1Cc1*param_$r6*ISc[1]-(pCc1Dc*param_$r8 + (1-pCc1Dc)*param_$r7)*Cc[1]
  dCcdt[2] <- (1-pCc1Dc)*param_$r7*Cc[1]-param_$r7*Cc[2]
  dCcdt[3] <- param_$r7*Cc[2]-param_$r7*Cc[3]
  dDcdt <- pCc1Dc*param_$r8*Cc[1]+param_$r8*param_$spdeath*ISc[2]
  dStCaredt <- param_$r7*Cc[3] - param_$r9*StCare
  dRcdt <- param_$r5*ISc[3] + param_$r9*StCare + param_$r4*IAc[3]
  
  out_<-vector(length=compartments$NumofVariables)
  out_[compartments$Sc_ind]<- dScdt
  out_[compartments$Ec_ind]<- dEcdt
  out_[compartments$IAc_ind] <- dIAcdt
  out_[compartments$ISc_ind] <- dIScdt
  
  out_[compartments$Cc_ind] <- dCcdt
  out_[compartments$Dc_ind] <- dDcdt
  out_[compartments$Rc_ind] <- dRcdt
  out_[compartments$StCare_ind] <- dStCaredt
  out_<-out_#delt*  ##c(dScdt,dIAcdt,dIScdt, dCcdt,dDcdt,dRc1dt,dRc2dt)###,dedrug)#dycyclo,dydoxo,dyetop,dpl
  out_
}
SumIncreasIS <- function(param_,states,delt,compartments,model_opt,pdeath_v)
{
  spdeath_v <- param_$spdeath0*pdeath_v/param_$pcrit0
  if(model_opt<=2)
  {
    if('psymp'%in%names(param_))
    {
      IAc2 <- states[,compartments$IAc_ind[2]]
      DICplus <- c(0,delt*param_$r3*param_$psymp*IAc2[1:(length(IAc2)-1)])#delt*param_$r3*param_$psymp*IAc2
    }else{
      IAc <- states[,2]
      
      DICplus<-c(0,delt*param_$r3*IAc[1:(length(IAc)-1)])#delt*param_$r3*IAc
      
    }
  }else if(model_opt ==3){
    IAc1 <- states[,compartments$IAc_ind[1]]
    pIA1IS1 <-X_p(r_1 = param_$rb4,r_2 = param_$r4,p_1 = param_$psymp )
    DICplus<-c(0,delt*pIA1IS1*param_$rb4*IAc1[1:(length(IAc1)-1)])+(1-pIA1IS1)*spdeath_v*states[,compartments$ISc_ind[2]]#unregistred Is after the death!!!!
  }else if((model_opt ==31)||(model_opt ==32)){
    IAc1 <- states[,compartments$IAc_ind[1]]
    IAMuc1 <- states[,compartments$IAMuc_ind[1]]
    pIA1IS1 <-X_p(r_1 = param_$rb4,r_2 = param_$r4,p_1 = param_$psymp )
    DICplus<-c(0,delt*pIA1IS1*param_$rb4*IAc1[1:(length(IAc1)-1)])+(1-pIA1IS1)*spdeath_v*states[,compartments$ISc_ind[2]]#unregistred Is after the death!!!!
    DICplus<-DICplus+c(0,delt*pIA1IS1*param_$rb4*IAMuc1[1:(length(IAMuc1)-1)])+(1-pIA1IS1)*spdeath_v*states[,compartments$ISMuc_ind[2]]#unregistred Is after the death!!!!
    
  }else if(model_opt ==4){
    Sc<- pmax(states[,compartments$Sc_ind],0)
    Ic <- states[,compartments$Ic_ind]
    Rc <- states[,compartments$Rc_ind]
    Nc <- Ic+Sc+Rc
    DICplus<-c(0,(delt*Sc*param_$beta*Ic/pmax(Nc,1))[1:(length(Ic)-1)])
  }
  out_<-vector(length=length(DICplus))
  if(model_opt ==31)
  {
    out_[1]<-states[1,compartments$ISc_ind[1]]+states[1,compartments$ISMuc_ind[1]]
  }
  else if(model_opt !=4)
  {
    out_[1]<-states[1,compartments$ISc_ind[1]]
  }else{
    out_[1]<-states[1,compartments$Ic_ind[1]]
  }
  
  for(indi in 2:length(DICplus) )
  {
    
    out_[indi]<-out_[indi-1]+DICplus[indi-1]
  }
  ####
  out_
}
###
SumIncreasISMut <- function(param_,states,delt,compartments,model_opt,pdeath_v)
{
  spdeath_v <- param_$spdeath0*pdeath_v/param_$pcrit0
  if(model_opt<=2)
  {
    if('psymp'%in%names(param_))
    {
      IAc2 <- states[,compartments$IAc_ind[2]]
      DICplus <- c(0,delt*param_$r3*param_$psymp*IAc2[1:(length(IAc2)-1)])#delt*param_$r3*param_$psymp*IAc2
    }else{
      IAc <- states[,2]
      
      DICplus<-c(0,delt*param_$r3*IAc[1:(length(IAc)-1)])#delt*param_$r3*IAc
      
    }
  }else if(model_opt ==3){
    IAc1 <- states[,compartments$IAc_ind[1]]
    pIA1IS1 <-X_p(r_1 = param_$rb4,r_2 = param_$r4,p_1 = param_$psymp )
    DICplus<-c(0,delt*pIA1IS1*param_$rb4*IAc1[1:(length(IAc1)-1)])+(1-pIA1IS1)*spdeath_v*states[,compartments$ISc_ind[2]]#unregistred Is after the death!!!!
  }else if(model_opt ==31){
    IAc1 <- states[,compartments$IAc_ind[1]]
    IAMuc1 <- states[,compartments$IAMuc_ind[1]]
    pIA1IS1 <-X_p(r_1 = param_$rb4,r_2 = param_$r4,p_1 = param_$psymp )
    DICplus<-c(0,delt*pIA1IS1*param_$rb4*IAc1[1:(length(IAc1)-1)])+(1-pIA1IS1)*spdeath_v*states[,compartments$ISc_ind[2]]#unregistred Is after the death!!!!
    DICplus<-cbind(DICplus,c(0,delt*pIA1IS1*param_$rb4*IAMuc1[1:(length(IAMuc1)-1)])+(1-pIA1IS1)*spdeath_v*states[,compartments$ISMuc_ind[2]])#unregistred Is after the death!!!!
    
  }else if(model_opt ==4){
    Sc<- pmax(states[,compartments$Sc_ind],0)
    Ic <- states[,compartments$Ic_ind]
    Rc <- states[,compartments$Rc_ind]
    Nc <- Ic+Sc+Rc
    DICplus<-c(0,(delt*Sc*param_$beta*Ic/pmax(Nc,1))[1:(length(Ic)-1)])
  }
  out_<-vector(length=length(DICplus))
  if(model_opt ==31)
  {
    out_<-array(dim=dim(DICplus))
    out_[1,1]<-states[1,compartments$ISc_ind[1]]+states[1,compartments$ISMuc_ind[1]]
    out_[1,2]<-states[1,compartments$ISMuc_ind[1]]
  }
  else if(model_opt !=4)
  {
    out_[1]<-states[1,compartments$ISc_ind[1]]
  }else{
    out_[1]<-states[1,compartments$Ic_ind[1]]
  }
  
  for(indi in 2:length(IAc1))#length(DICplus) )
  {
    if(model_opt ==31)
    {
      out_[indi,]<-out_[indi-1,]+DICplus[indi-1,]
    }else{
      out_[indi]<-out_[indi-1]+DICplus[indi-1]
    }
    
  }
  ####
  out_
}
###

TransPhitoKsi<-function(phi,parsim_label)
{
  ksi<-vector(length=9)
  if(parsim_label==0)#no parsimony
  {  
    ksi[1:9]<-exp(phi[1:9])
  }
  if(parsim_label==1)#r1=r2
  {  
    ksi[1:2]<-exp(phi[1])
    ksi[3:9]<-exp(phi[2:8])
  }
  if(parsim_label==2) 
  {  
    ksi[1:2]<-exp(phi[1])#r1=r2,
    #r5=r4= 1/14
    #r7=1/(7*4)
    #r1,r3,r6,r8 to estimate
    #r6 depends on the age!!!!
    #ksi[]
    ksi[3:9]<-exp(phi[2:8])
    
  }
  #ksi[16]<-log(2)#
  ###
  ksi
}

#############
#############

SimCOVIDStochastic<- function(xtot,main_struct,ind_indiv)
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
  sim_arr<-array(dim=c(xtot$num_steps,main_struct$numvar))
  sim_arr[1,]<- xtot$init
  if(main_struct$LogSc==1)
  {
    sim_arr[1,main_struct$compartments[[ind_indiv]]$Sc_ind]<- log(xtot$init[main_struct$compartments[[ind_indiv]]$Sc_ind])
  }
  sim_darray<-array(dim=c(xtot$num_steps,main_struct$numvar))
  sim_darray[1, ]<- rep.int(x= 0, times= main_struct$numvar)
  param_$r10<-param_$r1
  param_$r20<-param_$r2
  param_$r11<-param_$r1*param_$blockeff
  param_$r21<-param_$r2*param_$blockeff2# param_$r20<-param_$r2*param_$blockeff
  param_$r12<-param_$r11*param_$unlockeff
  param_$r13<-param_$r11*param_$unlockeff*param_$unlockeff1
  param_$r14<-param_$r11*param_$unlockeff*param_$unlockeff1*param_$unlockeff2
  #param_$r22<-param_$r20*param_$blockeff2
  
  #earlier:
  #param_$r10<-param_$r1*param_$blockeff
  #param_$r20<-param_$r2*param_$blockeff2# param_$r20<-param_$r2*param_$blockeff
  #param_$r11<-param_$r10*param_$blockeff2
  #param_$r21<-param_$r20*param_$blockeff2
  param_$influx0<-param_$influx
  for(indi in 2:xtot$num_steps){
    
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
      if(param_$Unlock<=indi)
      {
        param_$influx<-DelayEff(t0=indi,par1=0,par2=param_$unlockinflux,t1=param_$Unlock,delaylock=param_$delayunlock)
      }else
       { 
          param_$influx<-DelayEff(t0=indi,par1=param_$influx0,par2=0,t1=param_$Lock1,delaylock=param_$delaylock)
       }
    
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
      #Earlier: param_$r1<-DelayEff(t0=indi,par1=param_$r1,par2=param_$r10,t1=param_$Lock2,delaylock=param_$delaylock)
      # earlier: param_$r2<-DelayEff(t0=indi,par1=param_$r2,par2=param_$r20,t1=param_$Lock2,delaylock=param_$delaylock)
      
    if(param_$Unlock<=indi)
    {
      if(param_$Unlock2<=indi)
      {
        param_$r1<-DelayEff(t0=indi,par1=param_$r13,par2=param_$r14,t1=param_$Unlock,delaylock=param_$delayunlock)
      }else if(param_$Unlock1<=indi)
      {
        param_$r1<-DelayEff(t0=indi,par1=param_$r12,par2=param_$r13,t1=param_$Unlock,delaylock=param_$delayunlock)
      }else{
        param_$r1<-DelayEff(t0=indi,par1=param_$r11,par2=param_$r12,t1=param_$Unlock,delaylock=param_$delayunlock)
      }
      
    }else
    {
      param_$r1<-DelayEff(t0=indi,par1=param_$r10,par2=param_$r11,t1=param_$Lock2,delaylock=param_$delaylock)
    }
    param_$r2<-DelayEff(t0=indi,par1=param_$r20,par2=param_$r21,t1=param_$Lock2,delaylock=param_$delaylock)
   # print(c(indi,param_$Lock2, param_$Unlock,param_$r1,param_$r2))
    #}
    
   # if(param_$Lock3<=indi)
    #{
    #  param_$r1 <- param_$r11
    #  param_$r2 <- param_$r21
    #}
    #param_$r1<- (1-param_$block)*param_$r1+param_$r1*param_$block*param_$blockeff
    #param_$r2<- (1-param_$block)*param_$r2+param_$r2*param_$block*param_$blockeff
    #if(indi>=49)
    #{
    #  print(indi)
   # } 
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
    cbind(sim_arr,SumIncreasIS(param_,sim_arr,param_$delt,main_struct$compartments[[ind_indiv]],model_opt=main_struct$model_opt))
  }else{
    sim_arr0<-sim_arr
    sim_arr0[,main_struct$compartments[[ind_indiv]]$Sc_ind]<-exp(sim_arr0[,main_struct$compartments[[ind_indiv]]$Sc_ind])
    cbind(sim_arr0,SumIncreasIS(param_,sim_arr0,param_$delt,main_struct$compartments[[ind_indiv]],model_opt=main_struct$model_opt))
  }
}
#################
#################

CreateStochProtocolDataPrimitive<-function(main_struct,ind_indiv)#,delta)#delta in hours; no umatr!!!!
{
  measurements<- main_struct$measurements[[ind_indiv]]
  out_<-list()
  relev_times <- measurements$date#round(ProtocolNew$parttimes)
  out_$alldates <- sort(unique(measurements$date))##????
  all_times<-sort(unique(as.numeric(measurements$date-min(relev_times))))#min(relev_times):max(relev_times)
  num_steps<- length(main_struct$completedates[[ind_indiv]])#max(all_times)+1#length(all_times)
  out_$num_steps<-num_steps
  #umatr<-array(dim=c(num_steps,5))
  statesm<-vector(length=num_steps)
  #print(all_times)
  
  out_$all_times<- all_times
  ###measurements:
  relev_times<- measurements$date#!!!round(measurements$data[1,])
  statesm<-array(dim=c(num_steps,main_struct$YTYPES$num))#vector(length=num_steps)
  #print(all_times)
  ind_ytype0<-0
  for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
  {
    ind_ytype0 <- ind_ytype0+1
    for(indi0 in 1:length(all_times))#num_steps
    {
     # print('start')
     # print(c(0,0,indi0,ind_ytype0,dim(statesm),length(out_$alldates)))
      indi<-all_times[indi0]+1
      statesm[indi,ind_ytype0] <- NA
      cand_time<-out_$alldates[indi0]#all_times[indi]
     # print(c(indi,indi0,ind_ytype0,dim(statesm),length(out_$alldates)))
     # print(cand_time)
      if(cand_time%in%relev_times){
        relev_indi<-(1:length(relev_times))[cand_time==relev_times]
        #print(c(ind_ytype0,length(relev_times)))
        for(ind_meas in 1:length(relev_indi) )
        {  
          if((!is.na(measurements$data[2,relev_indi][ind_meas]))&(measurements$datatypes[relev_indi][ind_meas]==ind_ytype))
          {
            statesm[indi,ind_ytype0] <- measurements$data[2,relev_indi][ind_meas]
           # print(c(-10000,length(relev_times),ind_ytype0))
          }
        }
      }
    }
  }
  out_$statesm<-statesm
  ###
  out_
}

InitParam<-function(x,ind_indiv,phiopt,stoch_obj,main_struct)
{
  phiopt<-generate_listvCOVID(x,ind_indiv,xtot=phiopt,main_struct)#generate_listvect(x,phiopt,main_struct)
  param_<-ConstructParam(x=phiopt,ind_indiv,main_struct)###???ind_indiv=1
  
  #InitParam
  xtot<-list()
  xtot$phiopt <- phiopt
  
  start_ind <- min((1:length(stoch_obj$statesm))[!is.na(stoch_obj$statesm)])#1
  #xtot$x$stst<-c(rep.int(x=stoch_obj$statesm[start_ind]*xtot$x$kcirc/xtot$x$k,times= (num_transit+1)),stoch_obj$statesm[start_ind])
  
  
  #rep.int(x=0,times=5),
  if(main_struct$numindiv==1)
  {
    xtot<-InitUpdate(xtot,main_struct)
  }else 
  {
    if(main_struct$num_estim_indiv[ind_indiv]==0)
    {  
      xtot<-InitUpdate(xtot,main_struct)
    }else{
      xtot<-InitUpdateMultiple(xtot,ind_indiv,main_struct)
    }
  }
  xtot$sig <-0.1#residual error
  xtot$stochsd<-rep.int(x=0,times=main_struct$numvar)#main_struct$numvar
  xtot$stox <- vect2arr(x = rep.int(x=0,times = main_struct$numvar*stoch_obj$num_steps),nrows = stoch_obj$num_steps,ncols=main_struct$numvar)
  #xtot$x$kydoxo <- log(2)/0.5
  #xtot$x$ketop <- log(2)/0.5
  xtot$num_steps <- stoch_obj$num_steps
  xtot$prior <- RevTransformParametersAll(main_struct$priorval,main_struct)[main_struct$subs_indiv_num]#prior can be only those, that individually estimated!
  xtot$om <- main_struct$priorom[main_struct$subs_indiv_num]#
  
  ###
  xtot                      
}
UpdateParam<-function(x,ind_indiv,xtot,stoch_obj,main_struct)
{
  xtot$phiopt<-generate_listvCOVID(x,ind_indiv,xtot=xtot$phiopt,main_struct)#generate_listvect(x,phiopt,main_struct)
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)###???ind_indiv=1
   start_ind <- min((1:length(stoch_obj$statesm))[!is.na(stoch_obj$statesm)])#1
  #xtot$x$stst<-c(rep.int(x=stoch_obj$statesm[start_ind]*xtot$x$kcirc/xtot$x$k,times= (num_transit+1)),stoch_obj$statesm[start_ind])
  
  
  #rep.int(x=0,times=5),0
  if(main_struct$numindiv==1)
  {
    xtot<-InitUpdate(xtot,main_struct)
  }else
  {
    if(main_struct$num_estim_indiv[ind_indiv]==0)
    {  
      xtot<-InitUpdate(xtot,main_struct)
    }else{
      xtot<-InitUpdateMultiple(xtot,ind_indiv,main_struct)
    }
  }
  
   xtot$sig <-0.1#residual error
   xtot$stochsd<-rep.int(x=0,times=main_struct$numvar)#main_struct$numvar
   xtot$stox <- vect2arr(x = rep.int(x=0,times = main_struct$numvar*stoch_obj$num_steps),nrows = stoch_obj$num_steps,ncols=main_struct$numvar)
   #xtot$x$kydoxo <- log(2)/0.5
   #xtot$x$ketop <- log(2)/0.5
   xtot$num_steps <- stoch_obj$num_steps
  ###
  xtot                      
}
UpdateTreatParam<-function(ind_indiv,xtot,stoch_obj,main_struct)
{
   param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)###???ind_indiv=1
   if(main_struct$numindiv==1)
   {
     xtot<-InitUpdate(xtot,main_struct)
   }else
   {  
    if(main_struct$num_estim_indiv[ind_indiv]==0)
    {  
      xtot<-InitUpdate(xtot,main_struct)
    }else{
      xtot<-InitUpdateMultiple(xtot,ind_indiv,main_struct)
    }
  }
  xtot$stox <- vect2arr(x = rep.int(x=0,times = main_struct$numvar*stoch_obj$num_steps),nrows = stoch_obj$num_steps,ncols=main_struct$numvar)
  #xtot$x$kydoxo <- log(2)/0.5
  #xtot$x$ketop <- log(2)/0.5
  xtot$num_steps <- stoch_obj$num_steps
  ###
  xtot                      
}
Deriveylabels<-function(ind_ytype,label_cummul)#!!!instead complex_par!!!!
{
  ##no demargination1!!!! Otherwise it is complicated
  if(label_cummul==1)
  {
    add_name<-"Cumulative"
  }
  else{
    add_name <- "Daily"
  }
  out_ <-switch (ind_ytype, 
                 "Critical" = {'Daily critical states'},
                 "Recovery1" = {'Registered recoveries'},
                 "Death" = {paste(add_name, ' deaths',sep="")},
                 "Total" = {paste(add_name, ' registered cases',sep="")} 
                 
  );
  out_
}
vect2arr<-function(x,nrows,ncols)
{
  yout<-array(dim=c(nrows,ncols))
  for(indr in 1:nrows)
  {
    yout[indr,] <- x[((indr-1)*ncols+1):(indr*ncols)]
  }
  yout
}
arr2vect<-function(x,nrows,ncols)
{
  yout<-vector(length =nrows*ncols)
  for(indr in 1:nrows)
  {
    yout[((indr-1)*ncols+1):(indr*ncols)] <- x[indr,]
  }
  yout
}
#######
######
CalculateOutput<-function(S,ind_ytype,param_)#!!!instead complex_par!!!! complex_par,compartments,
{
  ##no demargination1!!!! Otherwise it is complicated
   
  #S[,compartments$]
  out_ <-switch (ind_ytype, 
                 "Total" = {param_$procentmeas*S[,'SumIsc'] },#{apply(S[,c('Cc','Dc', 'Rc1')],1,sum)+param_$procentmeas*S[,'ISc'] },#  ,'Rc0'
                 "Death" = {
                   #S[,'Dc']
                   if("Dc"%in%dimnames(S)[[2]]){#(is.vector(S[,'Cc']))
                     S[,'Dc']
                   }else{
                     1:dim(S)[1]#1+0*S[,1]#compartments$Ccnum
                     
                   }
                   },
                 "Recovery1" = {S[,'Rc1']},
                 "Critical" = {
                   if("Cc"%in%dimnames(S)[[2]]){#(is.vector(S[,'Cc']))
                     S[,'Cc']
                   }else{
                     if("Cc1"%in%dimnames(S)[[2]]){
                        apply(S[,paste('Cc',1:3,sep='')],1,sum)#compartments$Ccnum
                     }
                     else{
                       1:dim(S)[1]#1+0*S[,1]###!!!!!29.04.2020
                     }
                   }
                   },#S[,compartments$CEPO_cent_ind
                 "dTotal" = {param_$procentmeas*S[,'SumIsc'] },
                 "pcrit"= {S[,'pcrit']},
                 "pdeath"= {S[,'pdeath']},
                 "vr1"= {S[,'vr1']},
                 "vr2"= {S[,'vr2']},
                 "vrs"= {param_$spdeath*S[,'vr1']/param_$pcrit}
                 #  ,'Rc0'
                 #complex_par$paramnum$FERR_P_nor must be read from the initial data!!!!!
  );
  out_
}
GetPredictFromSimul<-function(tall,tobs,S,param_,ytype,dataytypes)#measurements$dataytypes,
{  #main_struct$logres(ind_ytypenum) option is abolished!!!!
  out_<-list()
  res_daily<-list()
  res_cumul<-list()
  timesi<-list()
  predicted_daily<- list()
  predicted_cumul<- list()
  for (ind_ytype0 in 1:ytype$num)
  {#main_struct$YTYPES$num#ytype$num
    ind_ytype <- ytype$list[[ind_ytype0]];
    ind_ytypenum <- ytype$listnum[ind_ytype0];
    if(sum(dataytypes==ind_ytype)>0)
    {  
      tobsi<-tobs[dataytypes==ind_ytype]
      timesi[[ind_ytypenum]] <- tobsi
      if(ind_ytype=="Critical")
      {  
        res_daily[[ind_ytypenum]] <- CalculateOutput(S,ind_ytype,param_)#CalculateOutput(S,i,ind_ytype1,main_struct,param_,opt_,compartments);
        res_cumul[[ind_ytypenum]] <- SumCalcl(res_daily[[ind_ytypenum]])
      }
      if(ind_ytype=="Death")
      {
        res_cumul[[ind_ytypenum]] <- Deathreportscumul(S,param_)
        res_daily[[ind_ytypenum]] <- DiffCalclAdv(res_cumul[[ind_ytypenum]],av_num= param_$DelayData)#c(0,
        #res_[[ind_ytypenum]]<-DeathreportsTransform(S,param_)
      }
      if(ind_ytype=="Total")
      {
        res_cumul[[ind_ytypenum]]<-NewCasesreportscumul(S,param_)###!!!04.01.2021
        res_daily[[ind_ytypenum]] <- DiffCalclAdv(res_cumul[[ind_ytypenum]],av_num= param_$DelayData)#c(0,
        #res_[[ind_ytypenum]]<-DeathreportsTransform(S,param_)
      }
      predicted_daily[[ind_ytypenum]] <- extractsolrelevindapprox(tall=tall,t=tobsi,ytall=res_daily[[ind_ytypenum]])#interp1(param_$datatime,res_ ,tobsi);#t);
      #if(ind_ytype!="Critical")
      #{
        predicted_cumul[[ind_ytypenum]] <- extractsolrelevindapprox(tall=tall,t=tobsi,ytall=res_cumul[[ind_ytypenum]])
      #}
    }
    
    
  }
  out_$res_cumul <- res_cumul
  out_$res_daily <- res_daily
  out_$predicted_cumul<-predicted_cumul
  out_$predicted_daily<-predicted_daily
  out_$timesi<-timesi
  #output!
  out_
}
 
GetResidErrFromPredicted<-function(simul_struct,ytype,measurements,main_struct,param_)#measurements$dataytypes,
{  #main_struct$logres[ind_ytypenum] option is abolished!!!!
  out_<-list()
  #out_$indivLL1 <- 0;
  #out_$indivLL10 <- 0;
  out_$quadrresidvar <- list();
  out_$quadrresidvarvect <- list();
  out_$residvarvect <- list()
  out_$dates<-list()
  out_$quadrresid <- vector(length = ytype$num)#!!!! 25.03.2021      ###vector(length = main_struct$YTYPES$num)#list();
  all_times<- sort(unique(measurements$date))#!!!!unique(measurements$data[1,])
  for (ind_ytype0 in 1:ytype$num)
  {#main_struct$YTYPES$num#ytype$num
    ind_ytype <- ytype$list[[ind_ytype0]];
    ind_ytypenum <- ytype$listnum[ind_ytype0];
    #acurr <- main_struct$acurr[ind_ytypenum];
    observedi<-measurements$data[measurements$namesmeasdata=="Y",measurements$datatypes==ind_ytype];#main_struct$measurements[[ind_indiv]]$datatypes!!!
    if(main_struct$expdist[ind_ytypenum]==1)
    {
      label_transf <- 'log'
    }else if (main_struct$grade_resid_err[ind_ytypenum]==1)
    {
      label_transf <- 'grade'
    }else{
      label_transf <- 'normal'
    }
    if((main_struct$label_cummul==0)&(ind_ytype!="Critical"))
    {
      #if(ind_ytype=="Total")
      #{  
      ##observedi <- DiffCalcl(zzdense)
        predi <- simul_struct$predicted_daily[[ind_ytypenum]]#!!!! 02.02.2021 DiffCalclAdv(simul_struct$predicted_daily[[ind_ytypenum]],av_num= param_$DelayData)#c(0,
      #}#else{
       # predi <- DiffCalcl1(simul_struct$predicted[[ind_ytypenum]])#Deaths
     # }
      observedi <- AverSim(observedi,av_num= param_$DelayData)#AverAdv
      if(ind_ytype=="Death")
      {
        predi[length(predi)]<-predi[length(predi)-1]#!!!21.0.2021!!!
      }
    }else{
      if(main_struct$label_cummul==0)
      {  
        predi <- simul_struct$predicted_daily[[ind_ytypenum]]
      
      }else {#else if(ind_ytype!="Critical")
        predi <- simul_struct$predicted_cumul[[ind_ytypenum]]
      }
    }
    #observedi[[ind_ytypenum]] <- measurements$data[measurements$namesmeasdata=="Y",measurements$datatypes==ind_ytype];#main_struct$measurements[[ind_indiv]]$datatypes!!!
    out_$dates[[ind_ytypenum]] <-measurements$date[measurements$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
    #xx <- measurements$date[measurements$datatypes==ytype$list[ind_ytype0]]#measurements$data[1,measurements$datatypes==ytype$list[ind_ytype0]]
    #relev_indi<-(all_times%in%xx)
    if(main_struct$week_label==1)
    {
      predi <- TransformtoWeek(x=predi,num_weeks = main_struct$num_weeks, erstsonntagind = main_struct$erstsonntagindi,
                               label_cummul=main_struct$label_cummul)
      observedi <- TransformtoWeek(x=observedi,num_weeks = main_struct$num_weeks, erstsonntagind = main_struct$erstsonntagindi,
                                   label_cummul=main_struct$label_cummul)
    }
    out_$quadrresidvarvect[[ind_ytypenum]] <- (TransformSol(y=predi,label_transf ,main_struct,ind_ytype,ind_ytypenum,main_struct$label_cummul)-
                                                 TransformSol(y = observedi,label_transf ,main_struct,ind_ytype,ind_ytypenum,main_struct$label_cummul))^2
    out_$residvarvect[[ind_ytypenum]] <- (TransformSol(y=predi,label_transf ,main_struct,ind_ytype,ind_ytypenum,main_struct$label_cummul)-
                                                 TransformSol(y = observedi,label_transf ,main_struct,ind_ytype,ind_ytypenum,main_struct$label_cummul))
    if((ind_ytype=="Critical")||(ind_ytype=="Death"))
    {
      if(ind_ytype=="Critical")
      {  
        if(main_struct$week_label==1)
        {
          relevi <-pmin(1,TransformtoWeek(x=measurements$RelevClinData*1,num_weeks = main_struct$num_weeks, erstsonntagind = main_struct$erstsonntagindi,
                                          label_cummul=main_struct$label_cummul))
        }else{
          relevi <- (measurements$RelevClinData*1)
        }
      }else if(ind_ytype=="Death")
      {
        relevi <- rep.int(x=1, times= length(out_$quadrresidvarvect[[ind_ytypenum]]))
        relevi[length(out_$quadrresidvarvect[[ind_ytypenum]])]<-0
      }
      
      
      out_$quadrresidvarvect[[ind_ytypenum]]<-out_$quadrresidvarvect[[ind_ytypenum]]*relevi
      out_$residvarvect[[ind_ytypenum]]<-out_$residvarvect[[ind_ytypenum]]*relevi#(measurements$RelevClinData*1)
      out_$quadrresidvar[[ind_ytypenum]]<- sum(out_$quadrresidvarvect[[ind_ytypenum]],na.rm=T)
      out_$quadrresid[ind_ytypenum] <- (out_$quadrresidvar[[ind_ytypenum]]/sum(relevi,na.rm = T))^0.5#sum(measurements$RelevClinData*1)
    }
    else{
      out_$quadrresidvar[[ind_ytypenum]]<- sum(out_$quadrresidvarvect[[ind_ytypenum]])
      out_$quadrresid[ind_ytypenum] <- (out_$quadrresidvar[[ind_ytypenum]]/length(observedi))^0.5;
    }
    if(ind_ytype=="Critical")
    {
      out_$numobs[[ind_ytypenum]] <- sum(measurements$RelevClinData)
    }else
    {  
      out_$numobs[[ind_ytypenum]] <- length(observedi)
    }
  }
  ###Output: 
  out_
}
TransformSol<-function(y,label_transf,main_struct,ind_ytype,ind_ytypenum,label_cummul)
{
  switch(label_transf,
         'log' = {
           if(label_cummul==1)
           {
             sensi<-main_struct$cumulexpsens[ind_ytypenum]
           }else{
             sensi<-main_struct$expsens[ind_ytypenum]
           }
           out_<-log(y+sensi)
         },
         'normal' = {
           out_ <- y^main_struct$transformres[ind_ytypenum]
         },
         'grade' = {
           #if(ind_ytype == 'PLC')
           #{
           out_<-ThromboGradeFuncVect(y,ind_ytype)
           #}
         }
         
  )
  out_
}
RevTransformSol<-function(y,label_transf,main_struct,ind_ytype,ind_ytypenum,label_cummul)
{
  switch(label_transf,
         'log' = {
           if(label_cummul==1)
           {
             sensi<-main_struct$cumulexpsens[ind_ytypenum]
           }else{
             sensi<-main_struct$expsens[ind_ytypenum]
           }
           out_<-exp(y)-sensi
         },
         'normal' = {
           out_ <- y^(1/main_struct$transformres[ind_ytypenum])
         },
         'grade' = {
           #if(ind_ytype == 'PLC')
           #{
           out_<-ThromboGradeFuncVect(y,ind_ytype)
           #}
         }
         
  )
  out_
} 
nnlfromx<-function(x,xtot,main_struct, ind_indiv)
{
  xtot$phiopt <- generate_listvect(x,xtot$phiopt,main_struct)
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)
  datatypesall<-as.character(main_struct$measurements[[ind_indiv]]$datatypes)
  tt<-main_struct$measurements[[ind_indiv]]$data[1,indtot]
  if(main_struct$label_treat_info==1)
  {
    sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
  }else{
    sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
  }
  #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
  #sim_res<-NameOutput(sim_res,main_struct)
  colnames(sim_res) <- main_struct$colnames_sim_res
  
  
  temp_predict<-GetPredictFromSimul(tall = tt,tobs= main_struct$measurements[[ind_indiv]]$data[1,] ,S = sim_res,param_,
                                    ytype = main_struct$YTYPE[[ind_indiv]],dataytypes = datatypesall)
  
  if((main_struct$label_cummul==0)&(ind_ytype!="Critical"))
  {
    measurementsi <- main_struct$Daily$measurements[[ind_indiv]]
  }else{
    measurementsi <- main_struct$measurements[[ind_indiv]]
  }
  temp_resid<-GetResidErrFromPredicted(simul_struct=temp_predict,ytype = main_struct$YTYPE[[ind_indiv]],measurements=measurementsi,main_struct,param_)
  
  out_<-0
  for(indi in 1:length(temp_resid$quadrresid))
  {
   out_ <- out_+main_struct$weights[indi]*temp_resid$quadrresid[indi]/(main_struct$acurr[indi]^ main_struct$transformres[indi])^2
  }
  ###
   out_
}
nnlfromxtot<-function(xtot,main_struct, ind_indiv)
{
  out_<- list()
  param_ <- main_struct$param_
  ###The old slow variant:
  #param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)
  #new variant, 20.11.2020!!!! for SimulCppCovidModel c++:
    param_<- IndivUpdateParam(x=xtot$phiopt,ind_indiv,main_struct,param_=main_struct$param_)
  #datatypesall<-as.character(main_struct$measurements[[ind_indiv]]$datatypes)
  #indtot<-main_struct$measurements[[ind_indiv]]$datatypes=='Total'
  if(main_struct$model_opt==32)
  {
    main_struct$muonedose <- 1
    main_struct$muvariant <- 'proportional'
    #main_struct$model_opt <- 31#mutant compartments
    sim_res<-SimCOVIDStochasticPolyMuKIAdvTreat(xtot,main_struct,ind_indiv)
    sim_res<-NamePolyOutput(sim_res_poly,main_struct)
  }else
  {  
    if(main_struct$label_treat_info==2)
    {
      sim_res<-SimulCppCovidModel(param_,xtot,ind_indiv,main_struct)
    }
    else if(main_struct$label_treat_info==1)
    {
      #sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
      #Rcpp!!!!:
      sim_res<-SimulCppCovidModel(param_,xtot,ind_indiv,main_struct)
      #sim_res<-SimulCppCovidModelIndiv(param_,xtot,ind_indiv,main_struct)
    }else{
      sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
    }
    sim_res<-NameOutput(sim_res,main_struct)
  }  #???colnames(sim_res) <- main_struct$colnames_sim_res
  #main_struct$measurements[[ind_indiv]]$data[1,]
  label_cummul0 <- main_struct$label_cummul#
  main_struct$label_cummul<-0
  temp_resid0 <- ResidCalc(sim_res,param_,main_struct, ind_indiv)
  main_struct$label_cummul<-1
  temp_resid1 <- ResidCalc(sim_res,param_,main_struct, ind_indiv)
  out_$residvarvect<- temp_resid0$residvarvect
  out_$residvarvect_cumul<- temp_resid1$residvarvect
  
  out_$dates <- temp_resid0$dates
  out_$dates_cumul <- temp_resid1$dates
  main_struct$label_cummul<- label_cummul0#Original value #
  out_$quadrresid <- temp_resid0$quadrresid#vector(length=length(temp_resid$quadrresid))
  out_$quadrresid_cumul <- temp_resid1$quadrresid#
  out_$nLL_tot<-0
  out_$resid_term<-0
  out_$resid_term_cumul<-0
  out_$qresid_term <- 0
  out_$qresid_term_cumul <- 0
  for(indi in 1:length(temp_resid0$quadrresid))
  {
    out_$qresid_term <-out_$qresid_term + temp_resid0$numobs[[indi]]*((temp_resid0$quadrresid[indi]/(main_struct$acurr[indi]^ main_struct$transformres[indi]))^2)
    #out_$quadrresid[indi]<-temp_resid$quadrresid[indi]
    out_$resid_term <- out_$resid_term + as.numeric(main_struct$weights[indi])*temp_resid0$numobs[[indi]]*(log(main_struct$acurr[indi]))
    out_$nLL_tot <- out_$nLL_tot+as.numeric(main_struct$weights[indi])*temp_resid0$numobs[[indi]]*(log(main_struct$acurr[indi])+(temp_resid0$quadrresid[indi]/(main_struct$acurr[indi]^ main_struct$transformres[indi]))^2)
  }
  for(indi in 1:length(temp_resid1$quadrresid))
  {
    #out_$quadrresid[indi]<-temp_resid$quadrresid[indi]
    if(main_struct$YTYPE[[ind_indiv]]$list[indi]!="Critical")
    {  
      out_$qresid_term_cumul <-out_$qresid_term_cumul+ as.numeric(main_struct$weights[indi])*temp_resid1$numobs[[indi]]*(temp_resid1$quadrresid[indi]/(main_struct$acurr_cumul[indi]^ main_struct$transformres[indi]))^2
      out_$resid_term_cumul <- out_$resid_term_cumul + as.numeric(main_struct$weights[indi])*temp_resid1$numobs[[indi]]*(log(main_struct$acurr_cumul[indi]))
      out_$nLL_tot <- out_$nLL_tot+param_$wecumul*as.numeric(main_struct$weights[indi])*temp_resid1$numobs[[indi]]*(log(main_struct$acurr_cumul[indi])+(temp_resid1$quadrresid[indi]/(main_struct$acurr_cumul[indi]^ main_struct$transformres[indi]))^2)
    }
  }
  ###
  if('prior'%in%names(xtot))
  {
    if('crit_prior_num'%in%names(param_))
    {
      ind_prev_crit <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("crit_exp",param_$crit_prior_num-1:4,sep ='')
      ind_prev_criti <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("crit_exp",param_$crit_prior_num-1,sep ='')
      xtot$prior[main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])]==paste("crit_exp",param_$crit_prior_num,sep ='')] <- xtot$phiopt$indiv[[ind_indiv]][ind_prev_criti]#mean(xtot$phiopt$indiv[[ind_indiv]][ind_prev_crit])
      xtot$om[main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])]==paste("crit_exp",param_$crit_prior_num,sep ='')] <-    sd(xtot$phiopt$indiv[[ind_indiv]][ind_prev_crit])
      
      
      ind_prev_death <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("death_exp",param_$death_prior_num-1:4,sep ='')
      ind_prev_deathi <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("death_exp",param_$death_prior_num-1,sep ='')
      xtot$prior[main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])]==paste("death_exp",param_$death_prior_num,sep ='')] <- xtot$phiopt$indiv[[ind_indiv]][ind_prev_deathi]#mean(xtot$phiopt$indiv[[ind_indiv]][ind_prev_death])
      xtot$om[main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])]==paste("death_exp",param_$death_prior_num,sep ='')] <-    sd(xtot$phiopt$indiv[[ind_indiv]][ind_prev_death])
      
      
      ind_prev_crit0 <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("crit_exp",1:(param_$crit_prior_num-1),sep ='')
      ind_prev_death0 <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("death_exp",1:(param_$death_prior_num-1),sep ='')
      
     # out_$nLL_tot <- out_$nLL_tot+param_$dynam_we*(sd(xtot$phiopt$indiv[[ind_indiv]][ind_prev_death0])+sd(xtot$phiopt$indiv[[ind_indiv]][ind_prev_crit0]))
      out_$nLL_tot <- out_$nLL_tot+param_$dynam_we*(SpecialSd(DiffCalcl(xtot$phiopt$indiv[[ind_indiv]][ind_prev_death0]))+SpecialSd(DiffCalcl(xtot$phiopt$indiv[[ind_indiv]][ind_prev_crit0])))
    }
    
    indiv_names <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])]
    start_date <- min(unique(main_struct$measurements[[ind_indiv]]$date))
    num_indiv_par <- length(xtot$phiopt$indiv[[ind_indiv]])
    num_treats <- param_$num_locks+param_$num_unlocks
    ind_date_prev_treat <- indiv_names %in% paste("date_treat",1:(num_treats-1),sep ='')
    ind_date_prev_treat_num <- (1:num_indiv_par)[ind_date_prev_treat]
    bad_region1<-as.numeric(c(as.Date("2020-12-19")-start_date,as.Date("2021-01-19")-start_date))
    bad_region2<-as.numeric(as.Date("2020-08-01")-start_date)
    dati1 <- exp(xtot$phiopt$indiv[[ind_indiv]][ind_date_prev_treat_num])
    selected_treat_dates_ind <- ((dati1<bad_region1[1])|(dati1>bad_region1[2]))&(dati1>bad_region2)
    ind_prev_treat <- indiv_names %in% paste("treat",1:(num_treats-1),sep ='')
    ind_prev_treat_num <- (1:num_indiv_par)[ind_prev_treat]
    
    relev_treat_indi <- as.numeric(sub(pattern = 'date_treat',replacement = '', x=indiv_names[ ind_date_prev_treat_num[selected_treat_dates_ind]]))
    relev_treat_indi1 <- relev_treat_indi[2:length(relev_treat_indi)][DiffCalcl(relev_treat_indi)==1]
    relev_treat_coeff_num <- ind_prev_treat_num[relev_treat_indi1]
    sd_tr <- SpecialSd(xtot$phiopt$indiv[[ind_indiv]][relev_treat_coeff_num]-xtot$phiopt$indiv[[ind_indiv]][relev_treat_coeff_num-1])
    out_$nLL_tot <- out_$nLL_tot+param_$dynam_we*sd_tr
    ind_prior <- main_struct$prior[main_struct$subs_indiv_num]==1
    out_$prior<- sum(((xtot$phiopt$indiv[[ind_indiv]]-xtot$prior)/xtot$om)[ind_prior]^2)
    out_$nLL_tot <- out_$nLL_tot+out_$prior
  }
  ###
  spdeath_v <- param_$spdeath*sim_res[,'pdeath']/param_$pcrit
  simdailydeathIs2<- c(0,sim_res[1:(xtot$num_steps-1),'ISc2']*param_$r8*spdeath_v[2:(xtot$num_steps)])#c(0,(sim_res_long1$ISc2[1:(stoch_obj0$num_steps-1)]*param_$r8*spdeath_v[2:(stoch_obj0$num_steps)]))
  simdailydeathCrit<- c(0,(sim_res[2:(xtot$num_steps),'pdeath']*param_$r8*sim_res[1:(xtot$num_steps-1),'Cc1']))
  fr_Is2 <- sum( simdailydeathIs2)/(sum( simdailydeathIs2)+sum(simdailydeathCrit))
  ##??? do not remember!!! 27.05.2021
  #out_$add_penalty <- ((fr_Is2-0.5)/0.05)^2
  #out_$nLL_tot<- out_$nLL_tot+ out_$add_penalty
  ##consistency:
  out_$co_term<-0
  if("death_int1"%in%names(param_))
  {
    help_str <- KICreateIntervals(param_,xtot,main_struct)
    intervals_crit <- help_str$intervals_crit
    intervals_death <- help_str$intervals_death
    if(!all(diff(intervals_death)>0))
    {
      out_$co_term<-10^8
    }
    if(!all(diff(intervals_crit)>0))
    {
      out_$co_term<-10^8
    }
    if(!all(diff(intervals_death)>2))
    {
      out_$co_term<-10^5
    }
    if(!all(diff(intervals_crit)>2))
    {
      out_$co_term<-10^5
    }
    if(!all(diff(intervals_death)>4))
    {
      out_$co_term<-10^4
    }
    if(!all(diff(intervals_crit)>4))
    {
      out_$co_term<-10^4
    }
  }
  if("treat_date_names"%in%names(main_struct))
  {
    for(indi_p in 2:length(main_struct$treat_date_names))
    {  
      if((param_[[main_struct$treat_date_names[indi_p]]]<= (param_[[main_struct$treat_date_names[indi_p-1]]]+2)))
      {
        out_$co_term<-10^8
      }else if((param_[[main_struct$treat_date_names[indi_p]]]<= (param_[[main_struct$treat_date_names[indi_p-1]]]+3)))
      {
        out_$co_term<-  out_$co_term + 10^6
      } 
      else if((param_[[main_struct$treat_date_names[indi_p]]]<= (param_[[main_struct$treat_date_names[indi_p-1]]]+4)))
      {
        out_$co_term<-  out_$co_term + 10^5
      } 
    }
    if(param_[[main_struct$treat_date_names[indi_p]]]>xtot$num_steps)
    {
      out_$co_term<-10^8
    }
    out_$unrealpdeath<-sum(sim_res[,'pdeath']>0.6)*10 #>0.8
  }
  out_$nLL_tot<- out_$nLL_tot+ out_$co_term+out_$unrealpdeath
  return(out_)
}
func_pop_thrombo<-function(x,xtot,main_struct, ind_indiv)
{  
  -nnlfromx(x,xtot,main_struct, ind_indiv)
}
GoalFuncComp<-function(xtot,main_struct, opt_mod, label_opt,nLL_indiv,ind_indiv1)
{
  #param_ <- main_struct$param_
  ygoalf<-list()
  main_struct$erstsonntagindi <- main_struct$erstsonntagind[ind_indiv]
  if(main_struct$model_opt==32)
  {
    sol_ <- nnlfromxtotpoly(xtot,main_struct, ind_indiv1)
  }else{
    sol_ <- nnlfromxtot(xtot,main_struct, ind_indiv1)
  }
  
  ygoalf$nLL_tot <-   sol_$nLL_tot
  ygoalf$nLL_indiv <-  ygoalf$nLL_tot
  ygoalf$quadrresid <- sol_$quadrresid
  ygoalf$resid_term <- sol_$resid_term 
  ygoalf$qresid_term <- sol_$qresid_term
  ygoalf$qresid_term_cumul <- sol_$qresid_term_cumul
  
  ygoalf$quadrresid_cumul <- sol_$quadrresid_cumul
  ygoalf$resid_term_cumul <- sol_$resid_term_cumul 
  ygoalf$residvarvect <- sol_$residvarvect
  ygoalf$residvarvect_cumul  <- sol_$residvarvect_cumul
  
  ygoalf$dates <- sol_$dates
  ygoalf$dates_cumul  <- sol_$dates_cumul
  
  ygoalf$prior <- sol_$prior
  #ygoalf$add_penalty <- sol_$add_penalty
  ygoalf$co_term <- sol_$co_term
  ygoalf$unrealpdeath <- sol_$unrealpdeath
  ygoalf
}
InitUpdate <- function(xtot,main_struct)
{
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv=1,main_struct)
  xtot$init<-c(param_$Sc0,rep.int(x=0,times= (main_struct$numvar-1)))#rep.int(x=0,times=5),0
  if("Iainit"%in%main_struct$names)
  {
    if(main_struct$model_opt==1)
    {
      xtot$init[2]<-param_$Iainit
    }
    if(main_struct$model_opt==2)
    {
      xtot$init[main_struct$compartments[[1]]$IAc_ind[1]] <-param_$Iainit# For COVID every Id have the same compartments number!!!!
      #xtot$init[2]<- (1-param_$pres)*param_$Iainit
      #xtot$init[3]<-param_$pres*param_$Iainit
    }
    if(main_struct$model_opt==3)
    {
      xtot$init[main_struct$compartments[[1]]$Ec_ind[1]] <-param_$Iainit# For COVID every Id have the same compartments number!!!!
      #xtot$init[2]<- (1-param_$pres)*param_$Iainit
      #xtot$init[3]<-param_$pres*param_$Iainit
    }
    if(main_struct$model_opt==4)
    {
      xtot$init[main_struct$compartments[[1]]$Ic_ind[1]] <-param_$Iainit#
    }
  }
  
  xtot
}
InitUpdateMultiple <- function(xtot,ind_indiv,main_struct)
{
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)
  xtot$init<-c(param_$Sc0,rep.int(x=0,times= (main_struct$numvar-1)))#rep.int(x=0,times=5),0
  if("Iainit"%in%main_struct$names)
  {
    if(main_struct$model_opt==1)
    {
      xtot$init[2]<-param_$Iainit
    }
    if(main_struct$model_opt==2)
    {
      xtot$init[main_struct$compartments[[ind_indiv]]$IAc_ind[1]] <-param_$Iainit# For COVID every Id have the same compartments number!!!!
      #xtot$init[2]<- (1-param_$pres)*param_$Iainit
      #xtot$init[3]<-param_$pres*param_$Iainit
    }
    if(main_struct$model_opt==3)
    {
      xtot$init[main_struct$compartments[[ind_indiv]]$Ec_ind[1]] <-param_$Iainit# For COVID every Id have the same compartments number!!!!
      #xtot$init[2]<- (1-param_$pres)*param_$Iainit
      #xtot$init[3]<-param_$pres*param_$Iainit
    }
    if(main_struct$model_opt==4)
    {
      xtot$init[main_struct$compartments[[ind_indiv]]$Ic_ind[1]] <-param_$Iainit#
    }
  }
  
  xtot
}
CollectPredictions <- function(predt,xtot,outpop,main_struct,stoch_obj0,ind_indiv,param_,schutd_sc,alpha_)#alpha_=0.05
{
  xtot$num_steps<-xtot$num_steps+predt
  stoch_obj0$num_steps <- xtot$num_steps
  xtot$stox <- vect2arr(x = rep.int(x=0,times = main_struct$numvar*stoch_obj0$num_steps),nrows = stoch_obj0$num_steps,ncols=main_struct$numvar)
  completedates0 <- main_struct$completedates[[ind_indiv]]
  main_struct$completedates[[ind_indiv]]<-(main_struct$completedates[[ind_indiv]][1]+(0:(xtot$num_steps-1)))
  
  main_struct$num_weeks<- (xtot$num_steps- main_struct$erstsonntagind[ind_indiv]+1)%/%7
  out_str <- list()
  out_str$altern_sc<-list()
  num_sol<-sum(outpop$accepted)
  out_str$accepted<-outpop$accepted
  if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))#(main_struct$label_treat_info==1)
  {  
    out_str$arr<-array(dim=c(xtot$num_steps,(main_struct$numvar+5),num_sol))
  }else{
    {  
      out_str$arr<-array(dim=c(xtot$num_steps,(main_struct$numvar+1),num_sol))
    }
  }
  ind_ <- 0
  if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))#(main_struct$label_treat_info==1)
  {
    mn <- vect2arr(x = rep.int(x=0,times = (main_struct$numvar+1)*xtot$num_steps),nrows = xtot$num_steps,ncols=(main_struct$numvar+3))#0*xtot$stox#rep.int(x=0,times=xtot$num_steps)
  }else{
    mn <- vect2arr(x = rep.int(x=0,times = (main_struct$numvar+1)*xtot$num_steps),nrows = xtot$num_steps,ncols=(main_struct$numvar+1))#0*xtot$stox#rep.int(x=0,times=xtot$num_steps)
  }
  sn <- mn#0*xtot$stox#rep.int(x=0,times=xtot$num_steps)
  for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
  {
    out_str[[ind_ytype]]<-list()
    out_str[[ind_ytype]]$mn <- rep.int(x=0, times = xtot$num_steps)
    out_str[[ind_ytype]]$sn <- rep.int(x=0, times = xtot$num_steps) 
    out_str[[ind_ytype]]$PerDay<-list()
  }
  sim_resy<-list()
  sim_resyp<-list()
  if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))
  {
    #sim_res_best<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
    #Rcpp:
    #sim_res<-SimulCppCovidModelIndiv(param_,xtot,ind_indiv,main_struct)
    param_<- IndivUpdateParam(x=xtot$phiopt,ind_indiv,main_struct,param_=main_struct$param_)
    sim_res_best<-SimulCppCovidModel(param_,xtot,ind_indiv,main_struct)
  }else{
    sim_res_best<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
  }
  #sim_res_best<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
  sim_res_best<-NameOutput(sim_res_best,main_struct)
  #colnames(sim_res_best) <- main_struct$colnames_sim_res
  xtot_alt_sc<-xtot
  if(!is.null(schutd_sc))
  {
    num_alt_scenarios<-length(schutd_sc)
  }else{
    num_alt_scenarios<-0
  }
  #if(main_struct$names[main_struct$subs_popi_num[[ind_indiv]]])
  ###alternative scenarios are skipped!!!
  #sim_res_full_shutd<-SimCOVIDStochastic(xtot,main_struct,ind_indiv=1)
  #sim_res<-NameOutput(sim_res,main_struct)
  #fits_res$xtot
  
  for(indi in 1:dim(outpop$candapp)[1] )
  {
    if(outpop$accepted[indi]==1)
    {
      ind_<-ind_+1
      xtot<-InitParam(x=outpop$candapp[indi,],ind_indiv,xtot,stoch_obj0,main_struct)
      #xtot$num_steps<-xtot$num_steps+predt
      if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))#(main_struct$label_treat_info==1)
      {
        #sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
        #Rcpp:
        #sim_res<-SimulCppCovidModelIndiv(param_,xtot,ind_indiv,main_struct)
        param_<- IndivUpdateParam(x=xtot$phiopt,ind_indiv,main_struct,param_=main_struct$param_)
        sim_res<-SimulCppCovidModel(param_,xtot,ind_indiv,main_struct)
      }else{
        sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
      }
      #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
      sim_res<-NameOutput(sim_res,main_struct)
      #colnames(sim_res) <- main_struct$colnames_sim_res
      #sim_res <- sim_res[,1:(main_struct$numvar+3)] #10.11.2020 ????
      sim_resp<-sim_res
      out_str$arr[,,ind_]<-sim_res
      #skip mn, sn!!! 10.11.2020
      #sn<-sn+sim_res+((indi*sim_res-mn)^2)/(indi*(indi+1))
     # mn<-mn+sim_res
      
      for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
      {
        if(ind_ytype=="Death")
        {
          sim_resy[[ind_ytype]]<-Deathreportscumul(S=sim_res[,],param_)
        }else if(ind_ytype=="Total")
          {
            sim_resy[[ind_ytype]]<-NewCasesreportscumul(S=sim_res,param_)###!!!04.01.2021
            #res_[[ind_ytypenum]]<-DeathreportsTransform(S,param_)
          }
        else  
        {
          sim_resy[[ind_ytype]]<- CalculateOutput(sim_res,ind_ytype,param_)
        }
        #skip mn, sn!!! 10.11.2020
        #out_str[[ind_ytype]]$sn <- out_str[[ind_ytype]]$sn + sim_resy[[ind_ytype]] + ((indi*sim_resy[[ind_ytype]]-out_str[[ind_ytype]]$mn)^2)/(indi*(indi+1))
        #out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mn+sim_resy[[ind_ytype]]
      }
      sim_resyp <- sim_resy
    }else{
      #skip mn, sn!!! 10.11.2020
      #sn<-sn+sim_resp+((indi*sim_resp-mn)^2)/(indi*(indi+1))
      #mn<-mn+sim_resp
      #for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
      #{
      #  out_str[[ind_ytype]]$sn <- out_str[[ind_ytype]]$sn  + sim_resyp[[ind_ytype]] + ((indi*sim_resyp[[ind_ytype]]-out_str[[ind_ytype]]$mn)^2)/(indi*(indi+1))
      #  out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mn+sim_resyp[[ind_ytype]]
     # }
      #out_str$arr[,,ind_]<-out_str$arr[,,ind_-1]#!!!!
    }
    print(c(ind_,indi))
  }
  out_str$mean <- mn/dim(outpop$candapp)[1]
  out_str$sd <- (sn/(dim(outpop$candapp)[1]+1))^0.5
  for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
  {
    #skip mn, sn!!! 10.11.2020
    #out_str[[ind_ytype]]$sn <- (out_str[[ind_ytype]]$sn/(dim(outpop$candapp)[1]+1))^0.5
    #out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mn/dim(outpop$candapp)[1]
    ###Ovewrite!!!
    
    arr3d0<-out_str$arr
    print('times of CalcCondifenceFun, ind_ytype = ')
    print(ind_ytype)
    if(ind_ytype!="Critical")
    {
      print(', cumul:')
      ptm <- proc.time()
      aaa<-AdvCalcConfidenceFun(outpop=outpop,arr3d=arr3d0,fun_=sum,arr_names=colnames(sim_res),
                             ind_ytype=ind_ytype,alpha_,param_,label_cummul=1,ind_indiv,main_struct)
      print((proc.time() - ptm))
      arr3d1 <- arr3d0#DiffCalclArr(out_str$arr) !!!10.04.2020, complicated, redo!!!!
    }
    print('times of CalcCondifenceFun, ind_ytype = ')
    print(ind_ytype)
    print(', daily:')
      
    ptm <- proc.time()
   
    bbb<-AdvCalcConfidenceFun(outpop=outpop,arr3d=arr3d1,fun_=sum,arr_names=colnames(sim_res),
                           ind_ytype=ind_ytype,alpha_,param_,label_cummul=0,ind_indiv,main_struct)
    print((proc.time() - ptm))
    if(ind_ytype!="Critical")
    {  
      out_str[[ind_ytype]]$mean <- aaa$mn#NameOutput(out_str$mean,main_struct)
      out_str[[ind_ytype]]$sd <-  aaa$sd#NameOutput(out_str$sd,main_struct)
      out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mean
    }
    out_str[[ind_ytype]]$PerDay$mean <- bbb$mn#NameOutput(out_str$mean,main_struct)
    out_str[[ind_ytype]]$PerDay$sd <-  bbb$sd#NameOutput(out_str$sd,main_struct)
    if(length(alpha_)==1)
    {  
      out_str[[ind_ytype]]$PerDay$lb <-  bbb$lb
      out_str[[ind_ytype]]$PerDay$ub <-  bbb$ub
      if(ind_ytype!="Critical")
      {  
        out_str[[ind_ytype]]$lb <-  aaa$lb
        out_str[[ind_ytype]]$ub <-  aaa$ub
      }
    }else{
      out_str[[ind_ytype]]$PerDay$quantile <-  bbb$quantile
      out_str[[ind_ytype]]$PerDay$week <- bbb$week
      if(ind_ytype!="Critical")
      {  
        out_str[[ind_ytype]]$quantile <-  aaa$quantile
        out_str[[ind_ytype]]$week <- aaa$week
      }
     }
  }
  out_str$week_dates <- aaa$week_dates
  ##
  arr3d0<-out_str$arr
  ptm <- proc.time()
  print('pcrit:')
  ind_ytype <-  'pcrit'
  aaa<-AdvCalcConfidenceFun(outpop=outpop,arr3d=arr3d0,fun_=sum,arr_names=colnames(sim_res),
                            ind_ytype= ind_ytype,alpha_,param_,label_cummul=0,ind_indiv,main_struct)
  out_str[[ind_ytype]]$mean <- aaa$mn#NameOutput(out_str$mean,main_struct)
  out_str[[ind_ytype]]$sd <-  aaa$sd#NameOutput(out_str$sd,main_struct)
  
  out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mean 
  
  #out_str[[ind_ytype]]$PerDay$mean <- bbb$mn#NameOutput(out_str$mean,main_struct)
 # out_str[[ind_ytype]]$PerDay$sd <-  bbb$sd#NameOutput(out_str$sd,main_struct)
  if(length(alpha_)==1)
  {  
    out_str[[ind_ytype]]$lb <-  aaa$lb
    out_str[[ind_ytype]]$ub <-  aaa$ub
  }else{
    out_str[[ind_ytype]]$quantile <-  aaa$quantile
    out_str[[ind_ytype]]$week <- aaa$week
  }
  print((proc.time() - ptm))
  ptm <- proc.time()
  print('pdeath:')
  ind_ytype <- 'pdeath'
  aaa<-AdvCalcConfidenceFun(outpop=outpop,arr3d=arr3d0,fun_=sum,arr_names=colnames(sim_res),
                            ind_ytype=ind_ytype,alpha_,param_,label_cummul=0,ind_indiv,main_struct)
  if(length(alpha_)==1)
  {  
    out_str[[ind_ytype]]$lb <-  aaa$lb
    out_str[[ind_ytype]]$ub <-  aaa$ub
  }else{
    out_str[[ind_ytype]]$quantile <-  aaa$quantile
    out_str[[ind_ytype]]$week <- aaa$week
  }
  print((proc.time() - ptm))
  ###
  print('vr1:')
  ind_ytype <- 'vr1'
  aaa<-AdvCalcConfidenceFun(outpop=outpop,arr3d=arr3d0,fun_=sum,arr_names=colnames(sim_res),
                            ind_ytype=ind_ytype,alpha_,param_,label_cummul=0,ind_indiv,main_struct)
  if(length(alpha_)==1)
  {  
    out_str[[ind_ytype]]$lb <-  aaa$lb
    out_str[[ind_ytype]]$ub <-  aaa$ub
  }else{
    out_str[[ind_ytype]]$quantile <-  aaa$quantile
    out_str[[ind_ytype]]$week <- aaa$week
  }
  print((proc.time() - ptm))
  ###
  print('vr2:')
  ind_ytype <- 'vr2'
  aaa<-AdvCalcConfidenceFun(outpop=outpop,arr3d=arr3d0,fun_=sum,arr_names=colnames(sim_res),
                            ind_ytype=ind_ytype,alpha_,param_,label_cummul=0,ind_indiv,main_struct)
  if(length(alpha_)==1)
  {  
    out_str[[ind_ytype]]$lb <-  aaa$lb
    out_str[[ind_ytype]]$ub <-  aaa$ub
  }else{
    out_str[[ind_ytype]]$quantile <-  aaa$quantile
    out_str[[ind_ytype]]$week <- aaa$week
  }
  
  print('vrs:')
  ind_ytype <- 'vrs'
  aaa<-AdvCalcConfidenceFun(outpop=outpop,arr3d=arr3d0,fun_=sum,arr_names=colnames(sim_res),
                            ind_ytype=ind_ytype,alpha_,param_,label_cummul=0,ind_indiv,main_struct)
  if(length(alpha_)==1)
  {  
    out_str[[ind_ytype]]$lb <-  aaa$lb
    out_str[[ind_ytype]]$ub <-  aaa$ub
  }else{
    out_str[[ind_ytype]]$quantile <-  aaa$quantile
    out_str[[ind_ytype]]$week <- aaa$week
  }
  
  print((proc.time() - ptm))
  
  out_str
}
CollectEstimations <- function(outpop,main_struct,ind_indiv,parsim_label)
{
  if(sum(outpop$results==max(outpop$results))>1)
  {  
    xres<-outpop$candapp[outpop$results==max(outpop$results),][1,]
  }else
  {
    xres<-outpop$candapp[outpop$results==max(outpop$results),]
  }
  accepted_good <-outpop$candapp[outpop$accepted==1,]
  corr_matr<-cor(accepted_good)
  colnames(corr_matr)<-main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]
  rownames(corr_matr)<-main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]
  all_orig_param<-array(dim=dim(accepted_good))
  best_res_orig<-vector(length=length(xres)) 
  for(ind_par0 in 1:length(xres))
  {
    if(main_struct$numpoppar>0)
    {
      ind_par <- main_struct$names[main_struct$subs_popi_num[[ind_indiv]]][ind_par0]
    }
    else
    {  
      ind_par<-main_struct$names[main_struct$subs_indivi[ind_indiv,]==1][ind_par0]
    }
    for(ind_row in 1:dim(accepted_good)[1])
    {  
      all_orig_param[ind_row,ind_par0]<-TransformParameters(accepted_good[ind_row,ind_par0],main_struct, ind_par)
    }
    best_res_orig[ind_par0]<-TransformParameters(xres[ind_par0],main_struct, ind_par)
  }
  
  mean_st_err_best<-list()
  mean_st_all_param<-list()
  if(main_struct$numpoppar>0)
  {
    mean_st_err_best$Parameter <- main_struct$names[main_struct$subs_popi_num[[ind_indiv]]]
  }else
  {
    mean_st_err_best$Parameter <-  main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]
  }
  mean_st_err_best$best_mu <- xres
  mean_st_err_best$mean_mu <- apply(accepted_good,2,mean)
  mean_st_err_best$std_mu <- apply(accepted_good,2,sd)
  mean_st_err_best$best <- best_res_orig
  mean_st_err_best$mean <- apply(all_orig_param,2,mean)
  mean_st_err_best$sterr <- apply(all_orig_param,2,sd)
  
  mean_st_all_param$Parameter <- main_struct$names
  mean_st_all_param$parsimony<-main_struct$names
  if(param_$parsim_label==4)
  {
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r2'] <- 'r1'
  }
  if(param_$parsim_label==3)
  {
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r4'] <- 'r5'
  }
  if(param_$parsim_label==2)
  {
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r2'] <- 'r1'
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r4'] <- 'r5'
  }  
  mean_st_all_param$estimated<-rep.int(x=0,times= main_struct$numallpar)
  for(ind_par0 in 1:main_struct$numallpar)
  {
    ind_par<-main_struct$names[ind_par0]
    if(ind_par %in%mean_st_err_best$Parameter)
    {
      mean_st_all_param$estimated[ind_par0]<-1
    }
  }   
      
      
  mean_st_all_param$LB <- main_struct$LB
  mean_st_all_param$UB <- main_struct$UB
  mean_st_all_param$Transform <- vector(length=length(main_struct$transform_label))
  mean_st_all_param$Transform[main_struct$transform_label==0]<-'Normal'
  mean_st_all_param$Transform[main_struct$transform_label==1]<-'Lognormal'
  mean_st_all_param$Transform[main_struct$transform_label==2]<-'Logit'
  mean_st_all_param$LB[main_struct$transform_label==0]<- - Inf
  mean_st_all_param$UB[main_struct$transform_label==0]<-  Inf
  mean_st_all_param$LB[main_struct$transform_label==1] <- 0
  for(ind_names in c('best_mu','mean_mu','std_mu','best','mean','sterr'))
  {
    
    mean_st_all_param[[ind_names]]<- rep.int(x=NA,times= main_struct$numallpar)#vector(length=main_struct$numallpar)
    if(ind_names=='mean')
    {
      mean_st_all_param[[ind_names]]<-main_struct$init_mean
    }
    for(ind_par0 in 1:main_struct$numallpar)
    {
      ind_par<-main_struct$names[ind_par0]
      if(ind_par %in%mean_st_err_best$Parameter)
      {
        mean_st_all_param[[ind_names]][ind_par0]<-mean_st_err_best[[ind_names]][ind_par ==mean_st_err_best$Parameter]
      }
      

    }
  }
  
  if(param_$parsim_label==4)
  {
    mean_st_all_param$mean[mean_st_all_param$Parameter=='r2'] <- mean_st_all_param$mean[mean_st_all_param$Parameter=='r1']
  }
  if(param_$parsim_label==3)
  {
    mean_st_all_param$mean[mean_st_all_param$Parameter=='r4'] <- mean_st_all_param$mean[mean_st_all_param$Parameter=='r5']
  }
  if(param_$parsim_label==2)
  {
    mean_st_all_param$mean[mean_st_all_param$Parameter=='r2'] <- mean_st_all_param$mean[mean_st_all_param$Parameter=='r1']
    mean_st_all_param$mean[mean_st_all_param$Parameter=='r4'] <- mean_st_all_param$mean[mean_st_all_param$Parameter=='r5']
    
  }  
  
  
  write.table(x=mean_st_err_best,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"EstimationResults",main_struct$covariates$Land[ind_indiv],".csv",sep=''))
  write.table(x=mean_st_all_param,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"EstimationResultsAllParameters",main_struct$covariates$Land[ind_indiv],".csv",sep=''))
  
  
  ###Output:
  mean_st_err_best
}

CollectEstimationsTogether <- function(fits_resi,main_struct,parsim_label)
{
  
  mean_st_all_param<-list()
  mean_st_all_param$Parameter <- main_struct$names
  mean_st_all_param$parsimony<-main_struct$names
  if(param_$parsim_label==4)
  {
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r2'] <- 'r1'
  }
  if(param_$parsim_label==3)
  {
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r4'] <- 'r5'
  }
  if(param_$parsim_label==2)
  {
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r2'] <- 'r1'
    mean_st_all_param$parsimony[mean_st_all_param$Parameter=='r4'] <- 'r5'
  }  
  mean_st_all_param$estimated<-rep.int(x=0,times= main_struct$numallpar)
  mean_st_all_param$LB <- main_struct$LB
  mean_st_all_param$UB <- main_struct$UB
  mean_st_all_param$Transform <- vector(length=length(main_struct$transform_label))
  mean_st_all_param$Transform[main_struct$transform_label==0]<-'Normal'
  mean_st_all_param$Transform[main_struct$transform_label==1]<-'Lognormal'
  mean_st_all_param$Transform[main_struct$transform_label==2]<-'Logit'
  mean_st_all_param$LB[main_struct$transform_label==0]<- - Inf
  mean_st_all_param$UB[main_struct$transform_label==0]<-  Inf
  mean_st_all_param$LB[main_struct$transform_label==1] <- 0
  estimated_params <- main_struct$names[main_struct$subs_indivi[ind_indiv=1,]==1]#not important which ind_indiv
  for(ind_par0 in 1:main_struct$numallpar)
  {
    ind_par<-main_struct$names[ind_par0]
    if(ind_par %in%estimated_params)
    {
      mean_st_all_param$estimated[ind_par0]<-1
    }
  }  
  
  mean_st_all_param1 <- mean_st_all_param#other order of fields
  for(ind_indiv in 1:main_struct$numindiv)
  {
    outpop<-fits_resi[[ind_indiv]]$outpop
    if(sum(outpop$results==max(outpop$results))>1)
    {  
      xres<-outpop$candapp[outpop$results==max(outpop$results),][1,]
    }else
    {
      xres<-outpop$candapp[outpop$results==max(outpop$results),]
    }
    all_orig_param<-array(dim=dim(outpop$candapp))
    best_res_orig<-vector(length=length(xres))
    for(ind_par0 in 1:length(xres))
    {
      ind_par<-main_struct$names[main_struct$subs_indivi[ind_indiv,]==1][ind_par0]
      
      for(ind_row in 1:dim(outpop$candapp)[1])
      {  
        all_orig_param[ind_row,ind_par0]<-TransformParameters(outpop$candapp[ind_row,ind_par0],main_struct, ind_par)
      }
      best_res_orig[ind_par0]<-TransformParameters(xres[ind_par0],main_struct, ind_par)
    }
     
    #####
    landi <- main_struct$covariates$Land[[ind_indiv]]
    for(ind_names in paste(c('best_mu','mean_mu','std_mu','best','mean','sterr'),landi,sep=''))
    {
      mean_st_all_param[[ind_names]]<- rep.int(x=NA,times= main_struct$numallpar)#vector(length=main_struct$numallpar)
      if(ind_names==paste('mean',landi,sep=''))
      {
        mean_st_all_param[[ind_names]]<-main_struct$init_mean
      }
    }
    mean_st_all_param[[paste('best_mu',landi,sep='')]][main_struct$subs_indiv_num] <- xres
    mean_st_all_param[[paste('mean_mu',landi,sep='')]][main_struct$subs_indiv_num] <- apply(outpop$candapp,2,mean)
    mean_st_all_param[[paste('std_mu',landi,sep='')]][main_struct$subs_indiv_num] <- apply(outpop$candapp,2,sd)
    mean_st_all_param[[paste('best',landi,sep='')]][main_struct$subs_indiv_num] <- best_res_orig
    mean_st_all_param[[paste('mean',landi,sep='')]][main_struct$subs_indiv_num] <- apply(all_orig_param,2,mean)
    mean_st_all_param[[paste('sterr',landi,sep='')]][main_struct$subs_indiv_num] <- apply(all_orig_param,2,sd)
    
    if(param_$parsim_label==4)
    {
      mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r2'] <- mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r1']
    }
    if(param_$parsim_label==3)
    {
      mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r4'] <- mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r5']
    }
    if(param_$parsim_label==2)
    {
      mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r2'] <- mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r1']
      mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r4'] <- mean_st_all_param[[paste('mean',landi,sep='')]][mean_st_all_param$Parameter=='r5']
      
    }
    
    #####
  }
  
  
  for(ind_names in c('best_mu','mean_mu','std_mu','best','mean','sterr') )
  {
    for(ind_indiv in 1:main_struct$numindiv)
    {
      landi <- main_struct$covariates$Land[[ind_indiv]]
      mean_st_all_param1[[paste(ind_names,landi,sep='')]]<-mean_st_all_param[[paste(ind_names,landi,sep='')]]
      
    }
    
  }
    
  
  write.table(x=mean_st_all_param1,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"EstimationResultsAllParametersAllLands.csv",sep=''))
  
  ###Output:
  mean_st_all_param
}

DelayEff<-function(t0,par1,par2,t1,delaylock)
{
  if(t0<=t1)
  {
    out_ <- par1
  }else{
    if(t0>(t1+delaylock))
    {
      out_ <- par2
    }else{
      out_<- as.numeric(par1+(par2-par1)*(t0-t1)/delaylock)
    }
  }  
  out_
}
NameOutput<-function(sim_res,main_struct)
{
  if(main_struct$model_opt==1)
  {
    colnames(sim_res)<-c('Sc','IAc','ISc', 'Cc','Dc', 'Rc1','Rc2','SumIsc')
  }
  else if(main_struct$model_opt==2)
  {
    #colnames(sim_res)<-c('Sc','IAc1','IAc2','ISc', 'Cc','Dc', 'Rc1','Rc2')
    compartments<-main_struct$compartments[[1]]
    colnames_sim_res <-vector(length=compartments$NumofVariables+1)
    colnames_sim_res[compartments$Sc_ind]<- 'Sc'
    colnames_sim_res[compartments$IAc_ind] <- c('IAc1','IAc2','IAc3')
    colnames_sim_res[compartments$ISc_ind] <- c('ISc1','ISc2','ISc3')
    colnames_sim_res[compartments$Cc_ind] <- paste('Cc',1:compartments$Ccnum,sep='')
    colnames_sim_res[compartments$Dc_ind] <- 'Dc'
    colnames_sim_res[compartments$Rc1_ind] <- 'Rc1'
    colnames_sim_res[compartments$Rc2_ind] <- 'Rc2'
    colnames_sim_res[compartments$NumofVariables+1]<-'SumIsc'
    colnames(sim_res) <- colnames_sim_res
  }
  else if(main_struct$model_opt==3)
  {
    #colnames(sim_res)<-c('Sc','IAc','ISc','ISc2', 'Cc','Dc', 'Rc1','Rc2','Rc0') 
    compartments<-main_struct$compartments[[1]]
    colnames_sim_res <-vector(length=compartments$NumofVariables+1)
    colnames_sim_res[compartments$Sc_ind]<- 'Sc'
    colnames_sim_res[compartments$Ec_ind]<- 'Ec'
    colnames_sim_res[compartments$IAc_ind] <- c('IAc1','IAc2','IAc3')
    colnames_sim_res[compartments$ISc_ind] <- c('ISc1','ISc2','ISc3')
    colnames_sim_res[compartments$Cc_ind] <- paste('Cc',1:compartments$Ccnum,sep='')
    colnames_sim_res[compartments$Dc_ind] <- 'Dc'
    colnames_sim_res[compartments$Rc_ind] <- 'Rc'
    colnames_sim_res[compartments$StCare_ind] <- 'StCare'
    colnames_sim_res[compartments$NumofVariables+1]<-'SumIsc'
    if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))
    {
      colnames_sim_res<-c(colnames_sim_res,c('vr1','vr2','pcrit','pdeath'))
    }  
    colnames(sim_res) <- colnames_sim_res
  }
  else if((main_struct$model_opt==31)||(main_struct$model_opt==32))
  {
    #colnames(sim_res)<-c('Sc','IAc','ISc','ISc2', 'Cc','Dc', 'Rc1','Rc2','Rc0') 
    compartments<-main_struct$compartments[[1]]
    colnames_sim_res <-vector(length=compartments$NumofVariables+1)
    colnames_sim_res[compartments$Sc_ind]<- 'Sc'
    colnames_sim_res[compartments$Ec_ind]<- 'Ec'
    colnames_sim_res[compartments$IAc_ind] <- c('IAc1','IAc2','IAc3')
    colnames_sim_res[compartments$ISc_ind] <- c('ISc1','ISc2','ISc3')
    colnames_sim_res[compartments$EMuc_ind]<- 'EMuc'
    colnames_sim_res[compartments$IAMuc_ind] <- c('IAMuc1','IAMuc2','IAMuc3')
    colnames_sim_res[compartments$ISMuc_ind] <- c('ISMuc1','ISMuc2','ISMuc3')
    colnames_sim_res[compartments$Cc_ind] <- paste('Cc',1:compartments$Ccnum,sep='')
    colnames_sim_res[compartments$Dc_ind] <- 'Dc'
    colnames_sim_res[compartments$Rc_ind] <- 'Rc'
    colnames_sim_res[compartments$StCare_ind] <- 'StCare'
    colnames_sim_res[compartments$NumofVariables+1]<-'SumIsc'
    if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))
    {
      colnames_sim_res<-c(colnames_sim_res,c('vr1','vr2','pcrit','pdeath'))
    }  
    colnames(sim_res) <- colnames_sim_res
  }
  else if(main_struct$model_opt==4)
  {
    #colnames(sim_res)<-c('Sc','IAc','ISc','ISc2', 'Cc','Dc', 'Rc1','Rc2','Rc0') 
    compartments<-main_struct$compartments[[1]]
    colnames_sim_res <-vector(length=compartments$NumofVariables+1)
    colnames_sim_res[compartments$Sc_ind]<- 'Sc'
    colnames_sim_res[compartments$Ic_ind]<- 'Ic'
    colnames_sim_res[compartments$Rc_ind]<- 'Rc'
    
    colnames_sim_res[compartments$NumofVariables+1]<-'SumIsc'
    colnames(sim_res) <- colnames_sim_res
  }
  ###
  sim_res
}
TransformDates<-function(times_date,t_start)
{
  times_num<-unique(as.numeric(times_date-t_start))#min(relev_times):max(relev_times)
  times_num
}
CompleteDates <- function(times_date)
{
  min_date<-min(times_date)
  times_num<-TransformDates(times_date,t_start=min_date)
  
  min_date+(min(times_num):max(times_num))
}
RunFit <- function(xtot,main_struct,nonmem_options,stoch_obj0)
{
  resid_err_count<-array(dim=c(2,length(main_struct$acurr)))
  for(ind_sim in 1:2)
  {
    outpop <-MCMCDomesticHaarGeneral(ind_indiv =1,xtot,main_struct,nonmem_options)
    if(sum(outpop$results==max(outpop$results))>1)
    {  
      xres<-outpop$candapp[outpop$results==max(outpop$results),][1,]
    }else
    {
      xres<-outpop$candapp[outpop$results==max(outpop$results),]
    }
    xtot<-InitParam(x=xres,ind_indiv=1,xtot$phiopt,stoch_obj0,main_struct)#outpop3$candapp[2000,]
    param_<-ConstructParam(x=xtot$phiopt,ind_indiv=1,main_struct)
    
    
    ggg<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                      nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=1)#
    
    #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv=1)
    if(main_struct$label_treat_info==1)
    {
      sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv=1)
    }else{
      sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv=1)#ind_indiv=1
    }
    #sim_res<-NameOutput(sim_res,main_struct)
    colnames(sim_res) <- main_struct$colnames_sim_res
    main_struct$acurr<-ggg$quadrresid
    resid_err_count[ind_sim,] <- main_struct$acurr
  }
  corr_matr<-cor(outpop$candapp)
  colnames(corr_matr)<-main_struct$names[main_struct$subs_popi_num[[1]]]
  rownames(corr_matr)<-main_struct$names[main_struct$subs_popi_num[[1]]]
  serr<-apply(outpop$candapp,2,sd)
  names(serr)<-colnames(corr_matr)
  out_<-list()
  
  #out_$main_struct <- main_struct
  out_$sim_res <- sim_res
  out_$xtot <- xtot
  out_$xres <- xres
  out_$outpop <- outpop
  out_$resid_err_count <- resid_err_count
  out_$corr_matr <- corr_matr
  out_$serr <- serr
  ######
  out_
}
RunFitIndiv <- function(xtot,ind_indiv,main_struct,nonmem_options,stoch_obj0)
{
  resid_err_count<-array(dim=c(2,length(main_struct$acurr)))
  for(ind_sim in 1:2)
  {
    outpop <-MCMCDomesticHaarGeneral(ind_indiv,xtot,main_struct,nonmem_options)
    if(sum(outpop$results==max(outpop$results))>1)
    {  
      xres<-outpop$candapp[outpop$results==max(outpop$results),][1,]
    }else
    {
      xres<-outpop$candapp[outpop$results==max(outpop$results),]
    }
    #xtot<- InitParam(x=xres,ind_indiv ,phiopt,stoch_obj0,main_struct)
    xtot<- UpdateParam(x=xres,ind_indiv,xtot,stoch_obj0,main_struct)#outpop3$candapp[2000,]
    param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)#ind_indiv=1
    
     
    ggg<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                      nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)#ind_indiv1=1
    
    if(main_struct$label_treat_info==1)
    {
      sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
    }else{
      sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
    }
    #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
    #sim_res<-NameOutput(sim_res,main_struct)
    colnames(sim_res) <- main_struct$colnames_sim_res
    main_struct$acurr<-ggg$quadrresid
    resid_err_count[ind_sim,] <- main_struct$acurr
  }
  corr_matr<-cor(outpop$candapp)
  colnames(corr_matr) <- main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]
  rownames(corr_matr) <- main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]
  serr<-apply(outpop$candapp,2,sd)
  names(serr)<-colnames(corr_matr)
  xtot$phiopt$resid[ind_indiv,] <- ggg$quadrresid#fits_res$resid_err_count[dim(fits_res$resid_err_count)[1],]
  out_<-list()
  
  #out_$main_struct <- main_struct
  out_$sim_res <- sim_res
  out_$xtot <- xtot
  out_$xres <- xres
  out_$outpop <- outpop
  out_$resid_err_count <- resid_err_count
  out_$corr_matr <- corr_matr
  out_$serr <- serr
  ######
  out_
}
CalcConfidenceFun<-function(outpop,arr3d,fun_,arr_names,ind_ytype,alpha_,param_,label_cummul,ind_indiv,main_struct)
{
  sub_arr <-switch (ind_ytype, 
                    "Critical" = {c('Cc1','Cc2','Cc3')},
                    "Recovery1" = {'Rc1'},
                    "Death" = {'Dc'},
                    "Total" = {'SumIsc'} 
  )
  out_<-list()
  out_$mn<-vector(length=dim(arr3d)[1])
  out_$sd<-vector(length=dim(arr3d)[1])
  if(length(alpha_)==1)
  {  
    out_$lb<-vector(length=dim(arr3d)[1])
    out_$ub<-vector(length=dim(arr3d)[1])
  }else{
    out_$quantile<- array(dim=c(length=dim(arr3d)[1],length(alpha_)))
  }
  help_arr<-array(dim=c(length(outpop$accepted),length(sub_arr)))
  help_arr <- array(dim=c(dim(arr3d)[1],length(outpop$accepted)))
  #help_arr[1]<-vecti[1]
  #ind_<-1
  ind_t <- 0
  for(indi in 1:dim(outpop$candapp)[1] )
  {
    
    ####
    
    if(outpop$accepted[indi]==1)
    {
      ind_t<-ind_t+1
      #    help_vect[indi] <- vecti[ind_]
    }#else{
    #  help_vect[indi]<- help_vect[indi-1]
    #}
    arr_t<-(arr3d[,,ind_t])
    colnames(arr_t)<- arr_names
    if(length(sub_arr)>1){
      vecti<-apply(arr_t[,sub_arr],1,fun_)
    }else
    {
      if(sub_arr=='SumIsc')
      {  
        vecti<-param_$procentmeas*arr_t[,sub_arr]#fun_( is sum!!!!
      }else if(sub_arr=='Dc')
      {
        #vecti<-Deathreportscumul(S=arr_t,param_)
        vecti<-Deathreportsonlycumul(S=arr_t,param_)
      }else{
        vecti<-arr_t[,sub_arr]#fun_( is sum!!!!
      }
    }
    vecti <- SimpDailyFromCumulative(x = vecti,label_cummul,ind_ytype,av_num = param_$DelayData)
    help_arr[,indi] <- vecti#[ind_]
     
    
  }
  
  for(ind_time in 1:dim(arr3d)[1])
  {
    help_vect <-help_arr[ind_time,]
    out_$mn[ind_time] <-mean(help_vect)
    out_$sd[ind_time] <-sd(help_vect)
    if(length(alpha_)==1)
    {  
      out_$lb[ind_time] <- sort(help_vect)[round((alpha_/2)*dim(outpop$candapp)[1],digits = 0)]
      out_$ub[ind_time] <- sort(help_vect)[round((1-alpha_/2)*dim(outpop$candapp)[1],digits = 0)]
    }else{
      for(ind_q in 1: length(alpha_))
      {
        out_$quantile[ind_time,ind_q ] <- sort(help_vect)[round((alpha_[ind_q])*dim(outpop$candapp)[1],digits = 0)]
      }
    }
  }
  
  
  
  out_
}
AdvCalcConfidenceFun<-function(outpop,arr3d,fun_,arr_names,ind_ytype,alpha_,param_,label_cummul,ind_indiv,main_struct)
{
  if(ind_ytype %in% main_struct$YTYPE[[ind_indiv]]$list)
  {  
    ind_ytypenum<- (1:length(main_struct$YTYPE[[ind_indiv]]$list))[main_struct$YTYPE[[ind_indiv]]$list==ind_ytype]
  
    if((label_cummul==1))
    {  
      sd_rand <- main_struct$acurr_cumul[ind_ytypenum]
    }else{
      sd_rand <- main_struct$acurr[ind_ytypenum]
    }
    if(is.na(sd_rand))
    {
      sd_rand <- 1
    }
  
    if(main_struct$expdist[ind_ytypenum]==1)
    {
      label_transf <- 'log'
    }else if (main_struct$grade_resid_err[ind_ytypenum]==1)
    {
      label_transf <- 'grade'
    }else{
      label_transf <- 'normal'
    }
  }else{
    sd_rand<- 0
  }
  out_<-list()
  out_$mn<-vector(length=dim(arr3d)[1])
  out_$sd<-vector(length=dim(arr3d)[1])
  out_$week <- list()
  if(length(alpha_)==1)
  {  
    out_$lb<-vector(length=dim(arr3d)[1])
    out_$ub<-vector(length=dim(arr3d)[1])
  }else{
    out_$quantile<- array(dim=c(length=dim(arr3d)[1],length(alpha_)))
    out_$week$quantile<- array(dim=c(length=main_struct$num_weeks,length(alpha_)))
  }
  #help_arr<-array(dim=c(length(outpop$accepted),length(sub_arr)))
  help_arr <- array(dim=dim(arr3d)[c(1,3)])#array(dim=c(dim(arr3d)[1],length(outpop$accepted)))
  help_arr_week <- array(dim=c(main_struct$num_weeks,dim(arr3d)[3]))
  #help_arr[1]<-vecti[1]
  #ind_<-1
  ind_t <- 0
  ###
  num_samp<-dim(outpop$candapp)[1]
  ind_accepted<-(1:num_samp)[outpop$accepted==1]
  num_accepted<-DiffCalcl(ind_accepted)
  N_accepted <- sum(outpop$accepted)
  if(max(ind_accepted)==num_samp)
  {
    num_accepted<-c(num_accepted,1)
  }else{
    num_accepted<-c(num_accepted,num_samp-max(ind_accepted)+1)
  }
  
  prob_accepted<-num_accepted/sum(num_accepted)
  #out_str[[paste('sol',ind_ytype,sep='')]] <- array(dim=dim(out_str$arr)[c(1,3)])
  bad_sol <- rep.int(x=0, times=N_accepted)
  for(ind_s in 1:N_accepted)
  {
    s_i<-arr3d[,,ind_s]
    colnames(s_i) <- arr_names
    
    
    if(ind_ytype=="Death")
    {
      help_arr[,ind_s]<-Deathreportscumul(S=s_i,param_)
    }else if(ind_ytype=="Total")
    {
      help_arr[,ind_s]<-NewCasesreportscumul(S=s_i,param_)###!!!04.01.2021
      #res_[[ind_ytypenum]]<-DeathreportsTransform(S,param_)
    }
    else{
      help_arr[,ind_s]<- CalculateOutput(S=s_i,ind_ytype,param_)#pmax(CalculateOutput(S=s_i,ind_ytype,param_),0)
    }
    if(sum(help_arr[,ind_s]<0))
    {
      bad_sol[ind_s] <- 1
    }
    if(sd_rand>0)
    {  
      help_arr[,ind_s] <- SimpDailyFromCumulative(x = help_arr[,ind_s],label_cummul,ind_ytype,av_num = param_$DelayData)
    }
    if(sum(is.na(help_arr[,ind_s]))>0)
    {
      ind_s
    }
  }
  
  #test_ <- sample(x=out_str[[paste('sol',ind_ytype,sep='')]][298,], size=100000, replace = TRUE, prob = num_accepted/sum(num_accepted))
  
  ###
  ind_0 <- main_struct$erstsonntagind[ind_indiv]
  week_resid<-0
  week_ind <- 0
  num_samp_all <- 100000
  help_vect_week<-rep.int(x=0,times=num_samp_all)
  for(ind_time in 1:dim(arr3d)[1])
  {
    
    help_vect0 <- sample(x=help_arr[ind_time,bad_sol==0], size=num_samp_all, replace = TRUE, prob = num_accepted[bad_sol==0]/sum(num_accepted[bad_sol==0]))
    if(sd_rand>0)
    {
      pert_ <- rnorm(num_samp_all,mean=0,  sd= sd_rand)
      help_vectp<- pmax(TransformSol(y=help_vect0,label_transf ,main_struct,ind_ytype,ind_ytypenum,label_cummul)+pert_,0)
      help_vect<- RevTransformSol(y=help_vectp,label_transf ,main_struct,ind_ytype,ind_ytypenum,label_cummul)
    }else{
      help_vect<- help_vect0
    }
    if(ind_time>=ind_0)
    {
      week_resid <- week_resid+1
      
      if(label_cummul==0)
      {
        help_vect_week <- help_vect_week+help_vect
      }else{
        help_vect_week <- help_vect
      }
      if(week_resid==7)
      {
        week_ind<- week_ind+1
        if(length(alpha_)>1){
          for(ind_q in 1: length(alpha_))
          {
            out_$week$quantile[week_ind,ind_q ] <- sort(help_vect_week)[round((alpha_[ind_q])*num_samp_all,digits = 0)]
          }
        }
        week_resid <- 0
        help_vect_week<-help_vect_week*0
      }
      
    }
    out_$mn[ind_time] <-mean(help_vect)
    out_$sd[ind_time] <-sd(help_vect)
    if(length(alpha_)==1)
    {  
      out_$lb[ind_time] <- sort(help_vect)[round((alpha_/2)*num_samp_all,digits = 0)]
      out_$ub[ind_time] <- sort(help_vect)[round((1-alpha_/2)*num_samp_all,digits = 0)]
    }else{
      for(ind_q in 1: length(alpha_))
      {
        out_$quantile[ind_time,ind_q ] <- sort(help_vect)[round((alpha_[ind_q])*num_samp_all,digits = 0)]
        if(is.na(out_$quantile[ind_time,ind_q ])>0)
        {
          ind_q
        }
      }
    }
  }
  out_$week_dates<-main_struct$completedate[[ind_indiv]][main_struct$erstsonntagind[ind_indiv]+7*(0:(main_struct$num_weeks-1))]
  
  
  out_
}
p_infect<-function(x,param_)
{
  
  #if(x==0)
  #{
  #  out_<-0
  #}
  #else{
  #out_<-x/(x+param_$sat_trans)#
  #}
x
}
X_p<-function(r_1,r_2,p_1)
{
  
  r_2*p_1/(r_1*(1-p_1)+r_2*p_1)
}
DailyFromCumulative<- function(x,label_cummul,ind_ytype,av_num)
{
  if((label_cummul==0)&(ind_ytype!="Critical"))
  {
    #zzdense<-DiffCalcl(zzdense)
    if(ind_ytype=="Total")
    {  
      out_ <- DiffCalclAdv(x,av_num)
    }else{
      out_ <- DiffCalcl1(x)
    }
  }else{
    out_ <- x
  }
  ####
  out_
}
SimpDailyFromCumulative<- function(x,label_cummul,ind_ytype,av_num)
{
  if((label_cummul==0)&(ind_ytype!="Critical"))
  {
    #zzdense<-DiffCalcl(zzdense)
    out_ <- DiffCalcl1(x)
  }else{
    out_ <- x
  }
  ####
  out_
}
CreateSimStruct<- function(sim_res,predt0,main_struct,ind_indiv)
{
  out_<-list()
  all_times<-unique(main_struct$measurements[[ind_indiv]]$date)
  xx0<-c(all_times,(max(all_times)+1):(max(all_times)+predt0))
  out_$Date<-xx0
  for(indi in colnames(sim_res))
  {
    out_[[indi]]<- sim_res[,indi]
  }
  
  out_
  
}
#generate_listvect
generate_listvCOVID<-function(x,ind_indiv,xtot,main_struct)#activation function of the first compartment
{
  if(main_struct$numpoppar>0)
  {    
    xtot$pop<-x[1:main_struct$numpoppar]
  }
  #sum_i<-main_struct$numpoppar
  #for(indi in 1:main_struct$numindiv)
  #{
  if(main_struct$numindiv>0)#??? it was > 1
  {  
    if(main_struct$num_estim_indiv[ind_indiv]>0)
    {  
      xtot$indiv[[ind_indiv]]<-x#[(1+sum_i):(sum_i+main_struct$num_estim_indiv[indi])] 
      #sum_i <- sum_i+main_struct$num_estim_indiv[indi]
    }  
  #}
  }
  xtot
}
PrepareSimul<- function(ind_indiv,xtot,main_struct,param_)#xtot_long,main_struct,ind_indiv
{ 
  if(main_struct$label_treat_info==1)
  {
    sim_res_long<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
    #sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
  }else{
    sim_res_long<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
  }
  #sim_res_long<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
  #sim_res_long<-NameOutput(sim_res_long,main_struct)
  colnames(sim_res_long) <- main_struct$colnames_sim_res
  sim_res_long1 <- CreateSimStruct(sim_res=sim_res_long,predt0,main_struct,ind_indiv)
  
  ind_ytype <- "Total"
  zzdense<-CalculateOutput(S=sim_res_long,ind_ytype ,param_)#sim_res[relev_indi,]
  zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul=0,ind_ytype,av_num = param_$DelayData)
  sim_res_long1$DailySumIsc <- zzdense
  
  
  ind_ytype <- "Death"
  zzdense<-CalculateOutput(S=sim_res_long,ind_ytype ,param_)#sim_res[relev_indi,]
  zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul=0,ind_ytype,av_num = param_$DelayData)
  sim_res_long1$DailyDeath <- zzdense
  return(sim_res_long1)
}
RunSaem <- function(xtot,ind_indiv,main_struct,nonmem_options,stoch_obj0)
{
  resid_err_count<-array(dim=c(2,length(main_struct$acurr)))
  for(ind_sim in 1:2)
  {
    outpop <-MCMCDomesticHaarGeneral(ind_indiv,xtot,main_struct,nonmem_options)
    if(sum(outpop$results==max(outpop$results))>1)
    {  
      xres<-outpop$candapp[outpop$results==max(outpop$results),][1,]
    }else
    {
      xres<-outpop$candapp[outpop$results==max(outpop$results),]
    }
    #xtot<- InitParam(x=xres,ind_indiv ,phiopt,stoch_obj0,main_struct)
    xtot<- UpdateParam(x=xres,ind_indiv,xtot,stoch_obj0,main_struct)#outpop3$candapp[2000,]
    param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct)#ind_indiv=1
    
    
    ggg<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                      nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)#ind_indiv1=1
    if(main_struct$label_treat_info==1)
    {
      sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
    }else{
      sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
    }
    #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
    #sim_res<-NameOutput(sim_res,main_struct)
    colnames(sim_res) <- main_struct$colnames_sim_res
    main_struct$acurr<-ggg$quadrresid
    resid_err_count[ind_sim,] <- main_struct$acurr
  }
  corr_matr<-cor(outpop$candapp)
  colnames(corr_matr) <- main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]
  rownames(corr_matr) <- main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]
  serr<-apply(outpop$candapp,2,sd)
  names(serr)<-colnames(corr_matr)
  xtot$phiopt$resid[ind_indiv,] <- ggg$quadrresid#fits_res$resid_err_count[dim(fits_res$resid_err_count)[1],]
  out_<-list()
  
  #out_$main_struct <- main_struct
  out_$sim_res <- sim_res
  out_$xtot <- xtot
  out_$xres <- xres
  out_$outpop <- outpop
  out_$resid_err_count <- resid_err_count
  out_$corr_matr <- corr_matr
  out_$serr <- serr
  ######
  out_
}
ResidCalc<-function(sim_res,param_,main_struct, ind_indiv)
{
  datatypesall<-as.character(main_struct$measurements[[ind_indiv]]$datatypes)
  if(main_struct$label_cummul==0)#&(ytype!="Critical")
  {
    measurementsi <- main_struct$Daily$measurements[[ind_indiv]]
    tobsii <- main_struct$Daily$measurements[[ind_indiv]]$date
    tt<- sort(unique(main_struct$Daily$measurements[[ind_indiv]]$date))#not [indtot]!!! other measurements can be later!!!#main_struct$measurements[[ind_indiv]]$data[1,indtot]
  }else{
    measurementsi <- main_struct$measurements[[ind_indiv]]
    tobsii <-main_struct$measurements[[ind_indiv]]$date
    tt<- sort(unique(main_struct$measurements[[ind_indiv]]$date))#not [indtot]!!! other measurements can be later!!!#main_struct$measurements[[ind_indiv]]$data[1,indtot]
  }
  tt<-tt[1]+0:as.numeric(max(tt)-min(tt))##Every day must be included!!!!
  temp_predict<-GetPredictFromSimul(tall = tt,tobs= tobsii ,S = sim_res,param_,
                                    ytype = main_struct$YTYPE[[ind_indiv]],dataytypes = datatypesall)
  temp_resid<-GetResidErrFromPredicted(simul_struct=temp_predict,ytype = main_struct$YTYPE[[ind_indiv]],measurements=measurementsi,main_struct,param_)
  return(temp_resid)
}
find_dist<-function(x)
{
  x1<-sort(unique(x))
  num_points <- length(x1)
  sum_x<-sum(x1)
  out_ <- array(dim = c(2,num_points))
  for(ind_x in 1:num_points)
  {
    out_[,ind_x] <- c(x1[ind_x],sum(x==x1[ind_x])/length(x))
  }
  return(out_)
}
KullbackLeiblerDistanceLognorm <-function(x,emp_d,cens_max)
{
  if(sum(x<0)==0)
  {  
  spec2 <-dlnorm(emp_d[1,], meanlog =x[1],sdlog = x[2], log = FALSE)
  spec0<-dlnorm(1:cens_max, meanlog =x[1],sdlog = x[2], log = FALSE)#minimal delay is one day!!!!
  #if(!is.finite(dgamma(0, shape=x[1],scale = x[2], log = FALSE)))
  #{
  #  spec0[1]<-
  #}
  spec2<-spec2/sum(spec0)
  spec1 <- emp_d[2,]
  nonz_ind<-(spec1>0)&(spec2>0)
  #print(c(x,(sum(spec1*log(spec1/spec2), na.rm = FALSE))))
  out_v<-spec1[nonz_ind]*log(spec1[nonz_ind]/spec2[nonz_ind])
  out_t<-sum(out_v[!is.infinite(out_v)], na.rm = FALSE)
  }else
  {
    out_t <- 100000
  }
  if(is.na(out_t))
  {  
    print(c(x,out_t))
  }
  return(out_t)
}
ComposeKullbackLeiblerGamma <-function(x,ti)#,emp_d,cens_max
{
  scale1 <- x[5]#as.numeric((20^2)/30)
  shape1<-  x[6]#as.numeric((30^2)/(sc*20^2))
  
  scalei <- x[3]#as.numeric((7.5^2)/12)
  shapei<-  x[4]#as.numeric((12^2)/(7.5^2))
  
  scale0 <- x[1]#as.numeric((0.8^2)/2)
  shape0<-  x[2]#as.numeric((2^2)/(sc*0.8^2))
  thetf <- x[7]#0.09
  thetf1 <- x[8]#0.65
  out_<-thetf*dgamma(ti, shape=shape0,scale = scale0, log = FALSE)+thetf1*dgamma(ti, shape=shapei,scale = scalei, log = FALSE)
  out_<- out_+(1-thetf-thetf1)*dgamma(ti, shape=shape1,scale = scale1, log = FALSE)
  
  return(out_)
}
ComposeKullbackLeiblerDistanceGamma <-function(x,emp_d,cens_max)
{
  if((sum(x<0)==0)&((x[7]+x[8])<1))
  {  
  spec2 <- ComposeKullbackLeiblerGamma(x,ti=emp_d[1,])
  spec0 <- ComposeKullbackLeiblerGamma(x,ti=1:cens_max)
  spec2<-spec2/sum(spec0)
  spec1 <- emp_d[2,]
  out_ <- sum(spec1*log(spec1/spec2), na.rm = FALSE)
  }else{
    out_ <- 10000
  }
  #print(c(x,(sum(spec1*log(spec1/spec2), na.rm = FALSE))))
  return(out_)
}
JointComposeKullbackLeiblerDistanceGamma <-function(x,distr_par,distr_m_l,scope_)
{
  out_ <- 0
  num_subj <- 0
  for(ind_m in scope_[1]:scope_[2])
  {
     emp_disti <- find_dist(distr_m_l[[ind_m]])
    num_subji <- as.numeric(distr_par[ind_m,'num'])
    num_subj<-num_subj+num_subji
    out_ <- out_ + num_subji*ComposeKullbackLeiblerDistanceGamma(x,emp_d=emp_disti,cens_max=as.numeric(distr_par[ind_m,'untilmax']))
    
  }
  out_<-out_/num_subj
  print(out_)
  #if(sum(x<0)>0)
  #{
 #   out_ <- 100000
 # }else
 # {  
 #   out_<-out_/num_subj
  #}
  return(out_)
}

KullbackLeiblerDistanceGamma <-function(x,emp_d,cens_max)
{
  
  spec2 <-dgamma(emp_d[1,], shape=x[1],scale = x[2], log = FALSE)
  spec0<-dgamma(1:cens_max, shape=x[1],scale = x[2], log = FALSE)#minimal delay is one day!!!!
  #if(!is.finite(dgamma(0, shape=x[1],scale = x[2], log = FALSE)))
  #{
  #  spec0[1]<-
  #}
  spec2<-spec2/sum(spec0)
  spec1 <- emp_d[2,]
  #print(c(x,(sum(spec1*log(spec1/spec2), na.rm = FALSE))))
  return(sum(spec1*log(spec1/spec2), na.rm = FALSE))
}
JointKullbackLeiblerDistanceGamma <-function(x,distr_par,distr_m_l,scope_)
{
  out_ <- 0
  num_subj <- 0
  for(ind_m in scope_[1]:scope_[2])
  {
    #Not necessary:
    #scalei <- as.numeric((sc*distr_par[ind_m,'sd']^2)/distr_par[ind_m,'mean'])
    #shapei<-  as.numeric((distr_par[ind_m,'mean']^2)/(sc*distr_par[ind_m,'sd']^2))
    emp_disti <- find_dist(distr_m_l[[ind_m]])
    num_subji <- as.numeric(distr_par[ind_m,'num'])
    num_subj<-num_subj+num_subji
    out_ <- out_ + num_subji*KullbackLeiblerDistanceGamma(x,emp_d=emp_disti,cens_max=as.numeric(distr_par[ind_m,'untilmax']))
    
  }
  if(sum(x<0)>0)
  {
    out_ <- 100000
  }else
  {  
  out_<-out_/num_subj
  }
  return(out_)
}
JointKullbackLeiblerDistanceLognorm <-function(x,distr_par,distr_m_l,scope_)
{
  out_ <- 0
  num_subj <- 0
  for(ind_m in scope_[1]:scope_[2])
  {
    
    emp_disti <- find_dist(distr_m_l[[ind_m]])
    num_subji <- as.numeric(distr_par[ind_m,'num'])
    num_subj<-num_subj+num_subji
    out_ <- out_ + num_subji*KullbackLeiblerDistanceLognorm(x,emp_d=emp_disti,cens_max=as.numeric(distr_par[ind_m,'untilmax']))
    
  }
  if(sum(x<0)>0)
  {
    out_ <- 100000
  }else
  {  
    out_<-out_/num_subj
  }
  print(c(x,out_))
  #out_<-out_/num_subj
  return(out_)
}
####
NewCasesreportsTransformDaily<-function(S,param_)
{
  
  
  param_$meanii <- 1.7745677
  param_$stdii <- 0.5305291
  zzdense<-S[,'SumIsc']#CalculateOutput(S,ind_ytype="Death",param_)#sim_res[relev_indi,]
  zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul=0,ind_ytype="Total",av_num = param_$DelayData)
  
  newcases_measured<-vector(length=length(zzdense))
  newcases_measured[1]<-0
  newcases_measured0<-rep.int(x=0, times=length(zzdense))
  sumil0<-sum(dlnorm(1:21, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE))#censored 21 days
  gam_v<-vector(length=(length(zzdense)-1))
  for(indi2 in 1:length(zzdense))
  {
    gam_v[indi2] <- dlnorm(indi2, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE)/sumil0
  }
  gam_v[22:length(gam_v)]<-0
  for(indi in 2:length(zzdense))
  {
    
    for(indi1 in 1:(indi-1))
    {
      newcases_measured[indi]<-newcases_measured[indi]+zzdense[indi1]*gam_v[indi-indi1]#dgamma(indi-indi1, shape=param_$shapei,scale = param_$scalei, log = FALSE)/sumi0
    }
    
  }
  return(param_$procentmeas*newcases_measured)
}
####
NewCasesreportscumul<-function(S,param_)
{
  
  newcases_measured<-NewCasesreportsTransformDaily(S,param_)
  out_<-rep.int(x=0,times=length(newcases_measured))
  for(indi in 2:length(newcases_measured))
  {
    out_[indi]<-out_[indi-1]+newcases_measured[indi]
  }
  return(out_)
}
####
DeathreportsTransform1<-function(S,param_)
{
  #param_$shapei <- 1.403895
  #param_$scalei <- 10.289708#!!!
  param_$shapei <- 1.345948
  param_$scalei <- 10.854349
  
  param_$meanii <- 2.407824
  param_$stdii <- 1.054611
  zzdense<-CalculateOutput(S,ind_ytype="Death",param_)#sim_res[relev_indi,]
  zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul=0,ind_ytype="Death",av_num = param_$DelayData)
  
  death_measured<-vector(length=length(zzdense))
  death_measured[1]<-0
  death_measured0<-rep.int(x=0, times=length(zzdense))
  sumi0<-sum(dgamma(1:1000, shape=param_$shapei,scale = param_$scalei, log = FALSE))
  sumil0<-sum(dlnorm(1:1000, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE))
  gam_v<-vector(length=(length(zzdense)-1))
  if(param_$deldistrgam==1)
  {  
    for(indi2 in 1:length(zzdense))
    {
      gam_v[indi2] <- dgamma(indi2, shape=param_$shapei,scale = param_$scalei, log = FALSE)/sumi0
    }
  }else{
    for(indi2 in 1:length(zzdense))
    {
      gam_v[indi2] <- dlnorm(indi2, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE)/sumil0
    }
  } 
  for(indi in 2:length(zzdense))
  {
    
    for(indi1 in 1:(indi-1))
    {
      death_measured[indi]<-death_measured[indi]+zzdense[indi1]*gam_v[indi-indi1]#dgamma(indi-indi1, shape=param_$shapei,scale = param_$scalei, log = FALSE)/sumi0
    }
    
  }
  return(death_measured)
}
Deathreportscumul1<-function(S,param_)
{
  
  death_measured<-DeathreportsTransform1(S,param_)
  out_<-rep.int(x=0,times=length(death_measured))
  for(indi in 2:length(death_measured))
  {
    out_[indi]<-out_[indi-1]+death_measured[indi]
  }
  return(out_)
}
CreateDistribDelayfunction <- function(S,param_)
{
  #param_$shapei <- 1.403895
  #param_$scalei <- 10.289708#!!!
  
  param_$shapei <- 1.345948
  param_$scalei <- 10.854349
  param_$meanii <- 2.407824
  param_$stdii <- 1.054611
  ##compose:
  #scope is 7:94
  param_$compdelay <- c(0.27185318, 6.85339157, 4.55802254, 2.56815863, 13.78361428,  2.08600438,  0.08077174, 0.60164892)
  #scope is 22:86
  param_$compdelay <- c(0.24994480,  7.28648036,  4.59225565,  2.50958035, 14.98464565,  1.75996584,  0.08300311,  0.58052233)
  
  t2<-dim(S)[1]#length(zzdense)
  
  #death_measured0<-rep.int(x=0, times=t2)
  sumi0<-sum(dgamma(1:1000, shape=param_$shapei,scale = param_$scalei, log = FALSE))
  sumil0<-sum(dlnorm(1:1000, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE))
  sumic0 <- sum(ComposeKullbackLeiblerGamma(x=param_$compdelay, ti=1:1000))
  gam_v<-vector(length=(t2-1))
  sgam_v<-rep.int(x=0, times=t2)
  for(indi2 in 1:t2)
  {
    #sumi0<-sum(dgamma(1:(t2-indi2), shape=param_$shapei,scale = param_$scalei, log = FALSE))
    if(param_$deldistrgam==1)
    {
      gam_v[indi2] <- dgamma(indi2, shape=param_$shapei,scale = param_$scalei, log = FALSE)/sumi0
    }else if (param_$deldistrgam==2){
      gam_v[indi2] <- dlnorm(indi2, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE)/sumil0
    }else
    {
      gam_v[indi2] <- ComposeKullbackLeiblerGamma(x=param_$compdelay, ti=indi2)/sumic0
    }
    #if(indi2>1)
    #{
    #   sgam_v[indi2]<-sgam_v[indi2-1]+gam_v[indi2]
    # }else{
    #   sgam_v[indi2]<-gam_v[indi2]
    # }
    
  }
  return(gam_v)
}


CreateNormDistribDelayfunction <- function(S,param_)
{
  #param_$shapei <- 1.403895
  #param_$scalei <- 10.289708#!!!
  
  me<-param_$mean_delay#(7+9+8+7+6+6+7)/7
  sd_=4
  sumi0 <-sum(dnorm(1:100,mean=me,sd=4))#4
  
  t2<-dim(S)[1]#length(zzdense)
  
  #death_measured0<-rep.int(x=0, times=t2)
  #sumi0<-sum(dgamma(1:1000, shape=param_$shapei,scale = param_$scalei, log = FALSE))
  #sumil0<-sum(dlnorm(1:1000, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE))
  #sumic0 <- sum(ComposeKullbackLeiblerGamma(x=param_$compdelay, ti=1:1000))
  gam_v<-vector(length=(t2-1))
  sgam_v<-rep.int(x=0, times=t2)
  for(indi2 in 1:t2)
  {
    gam_v[indi2] <- dnorm(indi2,mean=me,sd=sd_)/sumi0
    
  }
  return(gam_v)
}


CreateDistribDelaySimple <- function(t2)
{
  #param_$shapei <- 1.403895
  #param_$scalei <- 10.289708#!!!
  param_ <- list()
  param_$shapei <- 1.345948
  param_$scalei <- 10.854349
  param_$meanii <- 2.407824
  param_$stdii <- 1.054611
  ##compose:
  #scope is 7:94
  param_$compdelay <- c(0.27185318, 6.85339157, 4.55802254, 2.56815863, 13.78361428,  2.08600438,  0.08077174, 0.60164892)
  #scope is 22:86
  param_$compdelay <- c(0.24994480,  7.28648036,  4.59225565,  2.50958035, 14.98464565,  1.75996584,  0.08300311,  0.58052233)
  
  #t2<-dim(S)[1]#length(zzdense)
  
  #death_measured0<-rep.int(x=0, times=t2)
  sumi0<-sum(dgamma(1:1000, shape=param_$shapei,scale = param_$scalei, log = FALSE))
  sumil0<-sum(dlnorm(1:1000, meanlog = param_$meanii,sdlog = param_$stdii, log = FALSE))
  sumic0 <- sum(ComposeKullbackLeiblerGamma(x=param_$compdelay, ti=1:1000))
  gam_v<-vector(length=(t2-1))
  sgam_v<-rep.int(x=0, times=t2)
  for(indi2 in 1:t2)
  {
    gam_v[indi2] <- ComposeKullbackLeiblerGamma(x=param_$compdelay, ti=indi2)/sumic0
  }
  return(gam_v)
}
DeathreportsTransform<-function(S,param_)
{
  #gam_v<-CreateDistribDelayfunction(S,param_)
  gam_v <- CreateNormDistribDelayfunction(S,param_)
  zzdense<-CalculateOutput(S,ind_ytype="Death",param_)#sim_res[relev_indi,]
  zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul=0,ind_ytype="Death",av_num = param_$DelayData)
  t2<-dim(S)[1]#length(zzdense)
  death_measured<-vector(length=t2)
  death_measured[1]<-0
  for(indi in 2:(t2-1))
  {
    
    for(indi1 in 1:(indi-1))#for(indi1 in 1:(t2-indi))
    {
      #death_measured[indi]<-death_measured[indi]+zzdense[indi]*gam_v[indi1]#/sgam_v[(t2-indi)]#dgamma(indi-indi1, shape=param_$shapei,scale = param_$scalei, log = FALSE)/sumi0
    ##ECDC data:
      death_measured[indi]<-death_measured[indi]+zzdense[indi1]*gam_v[indi-indi1]
    }
    #print(c(indi,death_measured[indi]))
  }
  death_measured[t2]<-0#t1<t2
  ###alternative:
  #death_measured1<-vector(length=t2)
  #death_measured1[1]<-0
  #for(indi in 2:(t2))
  #{
    
    #for(indi1 in 1:(indi-1))
    #{
    #  death_measured1[indi]<-death_measured1[indi]+zzdense[indi1]*gam_v[indi-indi1]#/sgam_v[(t2-indi)]#dgamma(indi-indi1, shape=param_$shapei,scale = param_$scalei, log = FALSE)/sumi0
    #}
    #print(c(indi,death_measured[indi]))
  #}
  #death_measured1[t2]<-0#t1<t2
  
  
  return(death_measured)
}
Deathreportscumul<-function(S,param_)
{
  
  death_measured<-DeathreportsTransform(S,param_)
  out_<-rep.int(x=0,times=length(death_measured))
  for(indi in 2:length(death_measured))
  {
    out_[indi]<-out_[indi-1]+death_measured[indi]
  }
  return(out_)
}
###
Deathreportsonly<-function(S,param_)
{
  #gam_v<-CreateDistribDelayfunction(S,param_)
  gam_v <- CreateNormDistribDelayfunction(S,param_)
  zzdense<-CalculateOutput(S,ind_ytype="Death",param_)#sim_res[relev_indi,]
  zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul=0,ind_ytype="Death",av_num = param_$DelayData)
  t2<-dim(S)[1]#length(zzdense)
  
  death_reported<-vector(length=t2)
  death_reported[1]<-0
  for(indi in 2:(t2-1))
  {
    
    for(indi1 in 1:(indi-1))
    {
      death_reported[indi]<-death_reported[indi]+zzdense[indi1]*gam_v[indi-indi1]#/sgam_v[(t2-indi)]#dgamma(indi-indi1, shape=param_$shapei,scale = param_$scalei, log = FALSE)/sumi0
    }
    #print(c(indi,death_measured[indi]))
  }
  #death_measured[t2]<-0#t1<t2
  return(death_reported)
}
Deathreportsonlycumul<-function(S,param_)
{
  
  death_measured<-Deathreportsonly(S,param_)
  death_measured1<-DeathreportsTransform(S,param_)
  out_<-rep.int(x=0,times=length(death_measured))
  #xxhelp<-rep.int(x=0,times=length(death_measured))
  for(indi in 2:length(death_measured))
  {
    out_[indi]<-out_[indi-1]+death_measured[indi]
 #   xxhelp[indi]<- xxhelp[indi-1]+death_measured1[indi]
  }
  return(out_) 
}
###############
###############
resid_errf<-function(main_struct,ind_indiv)
{  
  measurements <- main_struct$Daily$measurements[[ind_indiv]]
  resiv<-list()
  observedi <- list()
  observedi1 <- list()
  ytype<-main_struct$YTYPE[[ind_indiv]]
  for (ind_ytype0 in 1:main_struct$YTYPE[[ind_indiv]]$num)
  { 
    ind_ytype <- ytype$list[[ind_ytype0]];
    ind_ytypenum <- ytype$listnum[ind_ytype0];
    observedi[[ind_ytype0]]<-measurements$data[measurements$namesmeasdata=="Y",measurements$datatypes==ind_ytype];#main_struct$measurements[[ind_indiv]]$datatypes!!!
    observedi1[[ind_ytype0]] <- AverSim(observedi[[ind_ytype0]],av_num= 7)
    resiv[[ind_ytype0]]<-abs(observedi1[[ind_ytype0]]-observedi[[ind_ytype0]])
  }
  out_<-list()
  out_$observedi<-observedi
  out_$observedi1<-observedi1
  out_$resiv<-resiv
  return(out_)
}  
resid_errfg<-function(x,main_struct,ind_indiv,ind_ytype0)
{
  measurements <- main_struct$Daily$measurements[[ind_indiv]]
  ytype<-main_struct$YTYPE[[ind_indiv]]
  ind_ytype <- ytype$list[[ind_ytype0]];
  ind_ytypenum <- ytype$listnum[ind_ytype0];
  observedi<-measurements$data[measurements$namesmeasdata=="Y",measurements$datatypes==ind_ytype];#main_struct$measurements[[ind_indiv]]$datatypes!!!
  observedi1 <- AverSim(observedi,av_num= 7)
  resi<-observedi^x-observedi1^x
  observedii<-observedi^x
  
  if(ind_ytype=='Critical')
  {
    out_<-(as.numeric(lm(resi[measurements$RelevClinData]~observedii[measurements$RelevClinData])$coefficients[2]))
  }else{
    out_<-(as.numeric(lm(resi~observedii)$coefficients[2]))
  }
  
  return(out_)
}
ComplexCumulDeathreportsTransform<-function(S,param_)
{
  t20<-dim(S)[1]
  #death_measured<-vector(length=t20)
  death_measuredcu<-rep.int(x=0, times =t20)#death_measuredcu[1]<-0
  for(t2 in 3:t20)
  {
    death_measuri <- DeathreportsTransform(S[1:t2,],param_)
    death_measuredcu[t2]<-sum(death_measuri)
  }
  
  return(death_measuredcu)
}
#plot(out_$observedi[[ind_ytype0]],out_$resiv[[ind_ytype0]],log='y')
#plot(observedi[[ind_ytype0]])
#lines(abs(observedi1[[ind_ytype0]]-out_$observedi[[1]]))
#lines(observedi1[[ind_ytype0]])
#lines(resiv[[ind_ytype0]])
#plot(log(out_$observedi[[1]]),log(out_$resiv[[1]]))
#lm(log(out_$resiv[[ind_ytype0]])~log(out_$observedi[[ind_ytype0]]))#-1

#lines(log(out_$observedi[[1]]),log(out_$observedi[[1]])*0.7952)
#lines(log(out_$observedi[[1]]),log(out_$observedi[[1]])*0.7952)
SimFromIndivMCMC<-function(ind_indiv,outpopi,xtot,stoch_obj0,main_struct)
{
  #ind_indiv<-1
  if(sum(outpopi[[ind_indiv]]$results==max(outpopi[[ind_indiv]]$results))>1)
  {  
    xres<-outpopi[[ind_indiv]]$candapp[outpopi[[ind_indiv]]$results==max(outpopi[[ind_indiv]]$results),][1,]
  }else
  {
    xres<-outpopi[[ind_indiv]]$candapp[outpopi[[ind_indiv]]$results==max(outpopi[[ind_indiv]]$results),]
  }
  stoch_obj0 <-CreateStochProtocolDataPrimitive(main_struct,ind_indiv)
  xtot<- UpdateParam(x=xres,ind_indiv,xtot,stoch_obj0,main_struct)
  #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
  if(main_struct$label_treat_info==1)
  {
    sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
  }else{
    sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
  }
  #sim_res<-NameOutput(sim_res,main_struct)
  colnames(sim_res) <- main_struct$colnames_sim_res
  return(sim_res)
}
SpecialSd<- function(x)
{
  (sum(x^2)/length(x))^0.5
}
InteravalPartition<- function(x)
{
  nx<- length(x)
  out_<-array(dim = c(nx,4))
  colnames(out_)<-c('sign','start','end','length')
  ni <- 1
  indi <- 2
  ind_cl<- 0
  while(indi<=nx)
  {  
    while( (x[indi]==x[ni])&(indi<=nx))
    {
      indi<- indi +1
    }
    ind_cl <- ind_cl+1
    out_[ind_cl,]<-c(x[ni],ni,indi-1,indi-ni)
    ni<-indi
    indi<- indi+1
  }
  return(out_[1:ind_cl,])
  
}

InteravalPartitionAdv<- function(ygoal,ind_ytype,cumul_label)
{
  
  if(cumul_label==1)
  {
    arr_i<-InteravalPartition(sign(ygoal$residvarvect_cumul[[ind_ytype]]))
    relev_dates <- ygoal$dates_cumul[[ind_ytype]]
  }else{
    arr_i<-InteravalPartition(sign(ygoal$residvarvect[[ind_ytype]]))
    relev_dates <- ygoal$dates[[ind_ytype]]
  }
  
  out_<- list()
  for(ind_n in colnames(arr_i))
  {
    out_[[ind_n]]<- arr_i[,ind_n]
  }
  
  out_$start_date<-relev_dates[out_$start]
  out_$end_date<-relev_dates[out_$end]
  
  return(out_)
  
}

UpdateMainstuctVparams <- function(param_,main_struct,label_data)
{
  
  if(label_data=='LandKI')
  {  
    main_struct$unlockeff_ind <- grep(pattern='unlockeff',x=names(param_))
    main_struct$num_unlocks <- length(main_struct$unlockeff_ind)
    
    main_struct$all_unlocks_ind <- grep(pattern = 'Unlock', x = names(param_))
    main_struct$all_locks_ind <- grep(pattern = 'Lock', x = names(param_))
    
    main_struct$all_unlocks <- names(param_)[main_struct$all_unlocks_ind]
    main_struct$all_locks <- names(param_)[main_struct$all_locks_ind]
    #main_struct$num_unlocks <- length(main_struct$all_unlocks_ind)
    main_struct$num_locks <- length(main_struct$all_locks_ind)
    main_struct$blockeff_ind <- grep(pattern='blockeff',x=names(param_))[1:main_struct$num_locks]#exclude blockeff_r2
    #main_struct$num_locks <- length(main_struct$blockeff_ind)
    main_struct$death_exp_ind<-grep(pattern='death_exp',x=names(param_))
    main_struct$crit_exp_ind<-grep(pattern='crit_exp',x=names(param_))
    main_struct$death_exp_n <- length(main_struct$death_exp_ind)+2
    main_struct$crit_exp_n <- length(main_struct$crit_exp_ind)+2
    main_struct$death_int_ind<-grep(pattern='death_int',x=names(param_))
    main_struct$crit_int_ind<-grep(pattern='crit_int',x=names(param_))
  }else if((label_data=='LandKIAdv')||(label_data=='LandMuKIAdv')||(label_data=='LandMuPolyKIAdv')||(label_data=='LandMuPolyKIAdvDiag'))
  {
    main_struct$death_exp_ind<-grep(pattern='death_exp',x=names(param_))
    main_struct$crit_exp_ind<-grep(pattern='crit_exp',x=names(param_))
    main_struct$death_exp_n <- length(main_struct$death_exp_ind)+2
    main_struct$crit_exp_n <- length(main_struct$crit_exp_ind)+2
    main_struct$death_date_ind<-grep(pattern='death_date',x=names(param_))
    main_struct$crit_date_ind<-grep(pattern='crit_date',x=names(param_))
   } 
  return(main_struct)
}
RecursiveCovMatrixAv <- function(arr_)
{
  dims<-dim(arr_)
  out_<-list()
  #out_$mu <- arr_[1,]
  #out_$cov_ <-  array(data = 0,dim =c(dims[2],dims[2]))
  out_$mu <- apply(arr_[1:2,],2,mean)
  out_$cov_ <- cov(arr_[1:2,])
  for(indi in 3:dims[1] )
  {
    out_ <- RecursiveCovMatrixAvUpd(xn=arr_[indi,],mup = out_$mu,covp = out_$cov_,nn = indi)
  }
  return(out_)
}
RecursiveCovMatrixAvUpd <- function(xn,mup,covp,nn)
{
  #previous mu, cov for n-1
  out_<-list()
  out_$mu <- ((nn-1)/nn)*mup+xn/nn
  delt_ <- xn-mup#out_$mu
  out_$cov_ <- ((nn-2)/(nn-1))*covp+(1/nn)*delt_%*%t(delt_)
  return(out_)
}
TransformtoWeek <- function(x,num_weeks,erstsonntagind,label_cummul)
{  
   week_resid<-0
  week_ind <- 0
  #num_samp_all <- 100000
  #help_vect_week<-rep.int(x=0,times=num_samp_all)
  
  y_out <- vector(length = num_weeks)
  for(ind_w in 1:num_weeks )
  {
    if(label_cummul==0)
    {
      y_out[ind_w] <- mean(x[((ind_w-1)*7+erstsonntagind):(7*ind_w+erstsonntagind)],na.rm=T)
    }else{
      y_out[ind_w] <- x[7*ind_w]
    }
  }
  return(y_out)
}
UpdateRandomParameters<-function(x,param_,main_struct)
{
  ind_prev_criti <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("crit_exp",param_$crit_prior_num-1,sep ='')
  x[main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])]==paste("crit_exp",param_$crit_prior_num,sep ='')] <- x[ind_prev_criti]#prev_crit_val#mean(xtot$phiopt$indiv[[ind_indiv]][ind_prev_crit])
  ind_prev_deathi <- main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])] %in% paste("death_exp",param_$death_prior_num-1,sep ='')
  x[main_struct$names[(main_struct$subs_indivi_num[[ind_indiv]])]==paste("death_exp",param_$death_prior_num,sep ='')] <- x[ind_prev_deathi]#mean(xtot$phiopt$indiv[[ind_indiv]][ind_prev_death])
  
  return(x)
}
ShiftLogicalVector <- function(x,del)
{
  first_t<-1
  while(x[first_t]==FALSE)
  {
    first_t <- first_t+1
  }
  x[(first_t-del):(first_t-1)]<-TRUE
  x[(length(x)-(del-1)):length(x)]<- FALSE
  return(x)
}
CompareLimits<-function(x,ind_indiv,main_struct)
{  
  
  if(main_struct$numindivpar1>0)
  {  
    limits_d<- array(dim=c(main_struct$num_estim_indiv[ind_indiv],5))
    for(ind_par in 1:main_struct$num_estim_indiv[ind_indiv])
    {
      # print(ind_par)
      namei <- main_struct$names_indiv_opt_parameters[[ind_indiv]][ind_par]#main_struct$names[main_struct$subs_indivi_num[[ind_indiv]]][ind_par]
      transformed_par <- TransformParameters(x[ind_par],main_struct,namei)
      param_[[namei]]<-transformed_par
      ind1 <- which(main_struct$names==namei)
      transform_label <-toString(main_struct$transform_label[ind1])
      if(transform_label=='0')
      {
        limits_d[ind_par,]<-c(param_[[namei]],-Inf,Inf)
      }else if(transform_label=='1')
      {
        limits_d[ind_par,]<-c(param_[[namei]],param_[[namei]],Inf,1,Inf)
      }else if(transform_label=='2')
      {
        limits_d[ind_par,]<-c(param_[[namei]],param_[[namei]]-main_struct$LB[ind1],main_struct$UB[ind1]-param_[[namei]],(param_[[namei]]-main_struct$LB[ind1])/param_[[namei]],(main_struct$UB[ind1]-param_[[namei]])/param_[[namei]])
      }
    }
  }
  colnames(limits_d)<-c('param','LB','UB','RLB','RUB')
  rownames(limits_d) <- main_struct$names_indiv_opt_parameters[[ind_indiv]]
  return(limits_d)
}
PerturbPar<-function(x,delt,ind_indiv,par_name,main_struct,transform_opt,rel_opt)
{
  out_ <-x
  #ind1<- which(main_struct$names==par_name)
  ind_par <-(1:main_struct$num_estim_indiv[ind_indiv])[main_struct$names_indiv_opt_parameters[[ind_indiv]]==par_name]
  if(transform_opt==1)
  {
    if(rel_opt==0)
    {  
      out_[ind_par]<-x[ind_par]+delt
    }else if(rel_opt==1)
    {
      out_[ind_par]<-x[ind_par]*delt
    }
  }else{
    transformed_par <- TransformParameters(x[ind_par],main_struct,par_name)
    if(rel_opt==0)
    {
      transformed_par1<- transformed_par+delt
    }else if(rel_opt==1)
    {
      transformed_par1<- transformed_par*delt
    }
    out_[ind_par] <- RevTransformParameters(transformed_par1,main_struct,par_name)
    
  }
  return(out_)
}
Sensitivanalysis<-function(xtot,delt,ind_indiv,main_struct,transform_opt,rel_opt)
{
  goalf0<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                       nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
  xtot0<-xtot
  sensitivity_r<- array(dim=c(main_struct$num_estim_indiv[ind_indiv],6))
  for(ind_par in 1:main_struct$num_estim_indiv[ind_indiv])
  {
    # print(ind_par)
    namei <- main_struct$names_indiv_opt_parameters[[ind_indiv]][ind_par]#main_struct$names[main_struct$subs_indivi_num[[ind_indiv]]][ind_par]
    print(namei)
    xpert <-PerturbPar(x=xtot$phiopt$indiv[[ind_indiv]],delt=delt,ind_indiv,par_name=namei,main_struct,transform_opt=0,rel_opt=1)
    xtot<- UpdateParam(x=xpert,ind_indiv,xtot,stoch_obj0,main_struct)
    goalf1<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                         nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
    sensitivity_r[ind_par,1] <- goalf1$nLL_tot-goalf0$nLL_tot
    sensitivity_r[ind_par,2] <- goalf1$prior-goalf0$prior
    sensitivity_r[ind_par,3] <- goalf1$co_term-goalf0$co_term
    sensitivity_r[ind_par,4] <- goalf1$unrealpdeath-goalf0$unrealpdeath
    sensitivity_r[ind_par,5] <- goalf1$qresid_term-goalf0$qresid_term
    sensitivity_r[ind_par,6] <- goalf1$qresid_term_cumul-goalf0$qresid_term_cumul
    xtot <- xtot0
  }
  colnames(sensitivity_r)<-c('nLL_tot','prior','bad_intervals','unrealpdeath','qresid','qresid_cumul')
  rownames(sensitivity_r) <- main_struct$names_indiv_opt_parameters[[ind_indiv]]
  return(sensitivity_r)
}