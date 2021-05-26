CreateComplexParameters<-function(x,ind_indiv_relev,param_,main_struct)#SimulIntegralModelShort  ind_biol,ind_indiv
{#par_indiv,i,main_struct
  #Actual!!!!!
  compartments <-main_struct$compartments[[ind_indiv_relev]]
  param_$blood_vol<-main_struct$blood_vol[ind_indiv_relev]#blood volune isa covariate of a data (in virtual it is different from the real one)
  #param_<-ConstructParam(x,ind_indiv,main_struct) #best_par<-ConstructParam(par_indiv=par_indiv,main_struct=main_struct)
  param_<-ParsimParam(param_,compartments)
  param_$pegtype<-1#default pegtype!!!!
  #Absolete parameter cpgbkdormrev from an earlier PGB model:
  #param_$cpgbkdormrev <- param_$asnor/((2^(param_$TPGBnorf/param_$cyclepgbmax)-2^(param_$TPGBnorf/param_$cyclepgbnor))/(2^(param_$TPGBnorf/param_$cyclepgbnor)-1))
  param_$nncmdorm <- compartments$Ccmdormnum
  param_$nncm <- compartments$Ccmnum
  if(sum(names(param_)=='pdcytost')>0)
  {main_struct$fun_name <- 'KomarBoneRemodelingranuloLymphoIntegralFastBandCeCycVPKPD3'}
  else
  {main_struct.fun_name <- 'KomarBoneRemodelingranuloLymphoIntegralFastBandCeCycVPKPD2'}
  
  # param_<-ParsimParam(param_,main_struct$specialpar)
  param_ <- SetToSteadyState(ind_indiv_relev,main_struct,param_)#must depend on relevant measurements. earlier: #ind_indiv
  
  
  #main_struct$now_biol<-0#!!!!Currently!!!1
  
  
  param_$nncmdorm <- compartments$Ccmdormnum
  param_$Acmnor <- 2^param_$nncmdt
  param_$Ccmdormlab <- 1
  param_<-InitiGranuloThromboKomarovaUpd(ind_indiv_relev,param_,compartments,main_struct)#ind_indiv_relev=ind_indiv
  
  param_$NG4=compartments$CG4num
  param_$NG5=compartments$CG5num
  param_$NG6=compartments$CG6num
  #  AA=zeros(compartments$NumofVariables,1)
  # param_$COBp_nor =COBp_nor 
  
  param_$listactivsyscomp=c(compartments$Cs_ind,compartments$Ccm_ind,compartments$CPGB_ind,compartments$Cmkc_ind[c(1,2,3,5,7,9,11)])
  param_$Cplcnorsum<-sum(param_$Cplcnor)
  param_$Cplsnorsum<-sum(param_$Cplsnor)
  param_$Cmkcnorsumtot<-sum(param_$Cmkcnor)
  param_$Cmkcnorsum<-sum(param_$Cmkcnor[c(1:(compartments$Cmkcnum-2),compartments$Cmkcnum)])
  param_$CG5norsum<-sum(param_$CG5nor)
  param_$CG4norsum<-sum(param_$CG4nor)
  param_$CG6norsum<-sum(param_$CG6nor)
  
  param_$CCmnorsum<-sum(param_$CCmnor)
  param_$CCmdormnorsum<-sum(param_$Ccmdorm)
  param_$CG6outnorlast<-param_$CG6outnor[compartments$CG6num]
  
  complex_par<-StratifyParam(param_)
  ###Outout:
  complex_par
}
InitiateComplexParameters<-function(x,ind_indiv,main_struct)#SimulIntegralModelShort
{#par_indiv,i,main_struct
  #Actual!!!!!
  #compartments <-main_struct$compartments[[ind_indiv]] + We do not need compartments hier!!!!
  param_<-ConstructParam(x,ind_indiv,main_struct) 
  
  ###Outout:
  param_
}