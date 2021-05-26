ModelScholzManyCompMut<- function(param_,states,compartments,delt)#u, ,mod_type DifModelFribergLog
{
  
  Sc<- max(states[compartments$Sc_ind],0)#exp(states[num_sys_comp+1])+
  Ec<- max(states[compartments$Ec_ind],0)#
  IAc <- states[compartments$IAc_ind]#exp(states[num_sys_comp+2])+
  ISc <- states[compartments$ISc_ind]#exp(states[num_sys_comp+3])+
  EMuc<- max(states[compartments$EMuc_ind],0)#
  IAMuc <- states[compartments$IAMuc_ind]#exp(states[num_sys_comp+2])+
  ISMuc <- states[compartments$ISMuc_ind]#exp(states[num_sys_comp+3])+
  
  Cc <- states[compartments$Cc_ind]#exp(states[num_sys_comp+5])+ basic level is neglidgible!!! for the thrombo model!!!
  Dc <- states[compartments$Dc_ind]
  Rc <- states[compartments$Rc_ind]#now Rc1,RC2
  StCare <- states[compartments$StCare_ind]#exp(states[num_sys_comp+4])+
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
  dScdt <- dScdt+param_$mur*(-param_$r1*(Sc/param_$Sc0)*(IAMuc[1]+IAMuc[2]+IAMuc[3]) - param_$r2*(Sc/param_$Sc0)*(ISMuc[1]+ISMuc[2]+ISMuc[3]))
  dCcdt[1] <- pIS1Cc1*param_$r6*(ISc[1]+ISMuc[1])-(pCc1Dc*param_$r8 + (1-pCc1Dc)*param_$r7)*Cc[1]
  dCcdt[2] <- (1-pCc1Dc)*param_$r7*Cc[1]-param_$r7*Cc[2]
  dCcdt[3] <- param_$r7*Cc[2]-param_$r7*Cc[3]
  dDcdt <- pCc1Dc*param_$r8*Cc[1]+param_$r8*param_$spdeath*(ISc[2]+ISMuc[2])
  dStCaredt <- param_$r7*Cc[3] - param_$r9*StCare
  dRcdt <- param_$r5*(ISc[3]+ISMuc[3]) + param_$r9*StCare + param_$r4*(IAc[3]+IAMuc[3])
  
  out_<-vector(length=compartments$NumofVariables)
  out_ <- InfectCompDerDeriv(out_,compartments,param_,pIA1IS1,pIS1Cc1,Sc,Ec,IAc,ISc,mur=1)
  out_ <- InfectCompDerDeriv(out_,compartments,param_,pIA1IS1,pIS1Cc1,Sc,Ec=EMuc,IAc=IAMuc,ISc=ISMuc,mur = param_$mur)
  out_[compartments$Sc_ind]<- dScdt
  #out_[compartments$Ec_ind]<- dEcdt
  #out_[compartments$IAc_ind] <- dIAcdt
  #out_[compartments$ISc_ind] <- dIScdt
  
  out_[compartments$Cc_ind] <- dCcdt
  out_[compartments$Dc_ind] <- dDcdt
  out_[compartments$Rc_ind] <- dRcdt
  out_[compartments$StCare_ind] <- dStCaredt
  #out_<-out_#delt*  ##c(dScdt,dIAcdt,dIScdt, dCcdt,dDcdt,dRc1dt,dRc2dt)###,dedrug)#dycyclo,dydoxo,dyetop,dpl
  return(out_)
}
InfectCompDerDeriv<-  function(out_,compartments,param_,pIA1IS1,pIS1Cc1,Sc,Ec,IAc,ISc,mur)
{
  ###fuck
  if(mur==1)
  {  
    dIAcdt <- vector(length = compartments$IAcnum)
    dIScdt <- vector(length = compartments$IScnum)  
   
  }else{
    dIAcdt <- vector(length = compartments$IAMucnum)
    dIScdt <- vector(length = compartments$ISMucnum) 
    
  }
  if(mur==1)
  {
    dEcdt <- param_$influx
  }else{
    dEcdt <- 0
  }  
  dEcdt <- dEcdt+ mur*param_$r1*(Sc/param_$Sc0)*(IAc[1]+IAc[2]+IAc[3]) + mur*param_$r2*(Sc/param_$Sc0)*(ISc[1]+ISc[2]+ISc[3])-param_$r3*Ec
  dIAcdt <- vector(length = compartments$IAcnum)
  dIScdt <- vector(length = compartments$IScnum)
 # if(mur==1)
 # {
 #   dIAcdt[1] <- param_$influx
 # }else{
 #   dIAcdt[1] <- 0
 # }  
  dIAcdt[1] <-  param_$r3*Ec - (pIA1IS1*param_$rb4 + (1-pIA1IS1)*param_$r4)*IAc[1]#dIAcdt[1] +
  dIAcdt[2] <- (1-pIA1IS1)*param_$r4*IAc[1]-param_$r4*IAc[2]#*param_$psymp*
  dIAcdt[3] <- param_$r4*IAc[2]-param_$r4*IAc[3]#param_$influx + 
  dIScdt[1] <- pIA1IS1*param_$rb4*IAc[1]-(pIS1Cc1*param_$r6 + (1-pIS1Cc1)*param_$r5)*ISc[1]
  dIScdt[2] <- (1-pIS1Cc1)*param_$r5*ISc[1]-(1-param_$spdeath)*param_$r5*ISc[2]-param_$r8*param_$spdeath*ISc[2]
  dIScdt[3] <- (1-param_$spdeath)*param_$r5*ISc[2]-param_$r5*ISc[3]
  if(mur==1)
  {  
    out_[compartments$Ec_ind]<- dEcdt
    out_[compartments$IAc_ind] <- dIAcdt
    out_[compartments$ISc_ind] <- dIScdt
  }else{
    out_[compartments$EMuc_ind]<- dEcdt
    out_[compartments$IAMuc_ind] <- dIAcdt
    out_[compartments$ISMuc_ind] <- dIScdt 
   } 
  return(out_)
}
