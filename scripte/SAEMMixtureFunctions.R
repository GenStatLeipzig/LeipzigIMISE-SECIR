SAEM_SMStepPecial<-function(gam_,outpopi,prev_pop,nonmem_options,relev_ind,plausible_ind,num_par,xtot,main_struct,stoch_objp)#notreat_ind,
{
  curr_pop <- prev_pop
  #burn<-nonmem_options$burn
  #num_sim<- dim(outpopi[[relev_ind[1]]]$candapp)[1]#nonmem_options$m
  num_sim <- sum(plausible_ind)
  plausible_num <- (1:length(plausible_ind))[plausible_ind]
  if(gam_==1)
  {  
    curr_pop$s1<-array(dim=c(max(relev_ind),num_par))
    curr_pop$s2 <- 0
    curr_pop$s5 <- 0
    curr_pop$s51 <- 0
  }else{
    curr_pop$s1 <- prev_pop$s1
    curr_pop$s2 <- prev_pop$s2
    curr_pop$s5 <- prev_pop$s5
    curr_pop$s51 <- prev_pop$s51
  }
  
  help_s5 <- array(dim=c(max(relev_ind),main_struct$YTYPES$num))#array(dim=c(length(relev_ind),main_struct$YTYPES$num))
  help_s51 <- array(dim=c(max(relev_ind),main_struct$YTYPES$num))#array(dim=c(length(relev_ind),main_struct$YTYPES$num))                 
  
  res_simstochp<-list()
  num_steps_tot<-0
  num_steps_no_treat<-0
  num_steps_meas<-0
  num_indiv <- length(relev_ind)
  num_measi <- array(dim=c(max(relev_ind),main_struct$YTYPES$num))#array(dim=c(length(relev_ind),main_struct$YTYPES$num))
  for(ind_indiv in relev_ind)
  {
    help_s5[ind_indiv,] <- rep.int(x=0, times = main_struct$YTYPES$num)
    help_s51[ind_indiv,] <- rep.int(x=0, times = main_struct$YTYPES$num)
    curr_pop$s1[ind_indiv,] <- StochUpd(x=prev_pop$s1[ind_indiv,],m=outpopi[[ind_indiv]]$candapp[plausible_num,1:num_par],gam_)
    #curr_pop$s1[ind_indiv,] <- prev_pop$s1[ind_indiv,]+gam_*(apply(outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],2,mean)-prev_pop$s1[ind_indiv,])
    #ind_measi <- (!is.na(stoch_objp[[ind_indiv]]$statesm))
    #num_measi <- sum(ind_measi)
    #stoch_obj0 <-CreateStochProtocolDataPrimitive(main_struct,ind_indiv)
    
    for(ind_pat in plausible_num)#(burn+1):num_sim
    {
      #xtot[[ind_indiv]]<-generate_listvect2(xx=outpopi[[ind_indiv]]$candapp[ind_pat,],xtot[[ind_indiv]],stoch_objp[[ind_indiv]],main_struct)
      #res_simstochp[[ind_indiv]]<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
      #help_s5 <- help_s5 + sum(((log(stoch_objp[[ind_indiv]]$statesm[ind_measi])-res_simstochp[[ind_indiv]][ind_measi,main_struct$num_transit+2]))^2)
      
      #xtot<-InitParam(x=outpopi[[ind_indiv]]$candapp[ind_pat,],ind_indiv,
      #                xtot,stoch_obj0,main_struct)#
      #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
      #sim_res<-NameOutput(sim_res,main_struct0)
      #xtot<-UpdateParam(x=fits_resi[[ind_indiv]]$xres,ind_indiv ,xtot,stoch_obj0,main_struct0)
      #xtot <- UpdateTreatParam(ind_indiv,xtot,stoch_obj0,main_struct)
      #temp_resid <- ResidCalc(sim_res,param_,main_struct, ind_indiv)
      help_s5[ind_indiv,] <- help_s5[ind_indiv,] + (outpopi[[ind_indiv]]$residapp[ind_pat,]^2)*as.numeric(main_struct$measurements[[ind_indiv]]$length_data)#temp_resid$numobs#temp_resid$quadrresid*temp_resid$numobs
      help_s51[ind_indiv,] <- help_s51[ind_indiv,] + (outpopi[[ind_indiv]]$residapp_cumul[ind_pat,]^2)*as.numeric(main_struct$measurements[[ind_indiv]]$length_data)#temp_resid$numobs#temp_resid$quadrresid*temp_resid$numobs
      #goals<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
      #                    nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)#  
    }
    num_measi[ind_indiv,] <- as.numeric(main_struct$measurements[[ind_indiv]]$length_data)#temp_resid$numobs
    #num_steps_meas <- num_steps_meas+num_measi#population residual error is calculated only for untreated subjects!!! 
    
  }
  #average on chains:
  help_s5 <- help_s5/num_sim#(num_sim-burn)
  help_s51 <- help_s51/num_sim#(num_sim-burn)
  #approximation step:
  
  curr_pop$s5<- curr_pop$s5+gam_*(help_s5-curr_pop$s5)
  curr_pop$s51<- curr_pop$s51+gam_*(help_s51-curr_pop$s51)
  if(num_indiv>1)
  { 
    curr_pop$mu<-apply(curr_pop$s1[relev_ind,],2,mean)
    ind_indiv <- relev_ind[1]
    help1<-AddMatrVectRows(m=outpopi[[ind_indiv]]$candapp[plausible_num,1:num_par],x= - curr_pop$mu)^2#(burn+1):num_sim
    for(ind_indiv in relev_ind[2:length( relev_ind)])
    {
      help1 <- help1 + AddMatrVectRows(m=outpopi[[ind_indiv]]$candapp[plausible_num,1:num_par],x= - curr_pop$mu)^2#(burn+1):num_sim
    }
    curr_pop$s2 <- StochUpd(x=curr_pop$s2,m=help1,gam_)#mean on chains
    curr_pop$om<- (curr_pop$s2/num_indiv)^0.5
  }
  #if(sum(relev_ind%in%notreat_ind)>0)
  #{  
  curr_pop$sig_resid <- (curr_pop$s5/num_measi)^0.5#num_steps_meas
  curr_pop$sig_resid_cumul <- (curr_pop$s51/num_measi)^0.5#num_steps_meas
  #}   
  
  #s3,s31,s4,s41
  #curr_pop$stochsd<- (curr_pop$s3/(num_steps_tot))^0.5
  #curr_pop$stochsd[main_struct$num_transit+2]<- (curr_pop$s31/(num_steps_no_treat))^0.5
  #curr_pop$sig2<-vector(length = dim(xtot[[relev_ind[1]]]$stox)[2])
  ##if(sum(relev_ind%in%notreat_ind)>0)
  ##{  
  # curr_pop$sig2[1:3]<- ((curr_pop$s4/(num_steps_no_treat))^0.5)[1:3]#only untreated patients are scored
  ## }
  #curr_pop$sig2[4:dim(xtot[[relev_ind[1]]]$stox)[2]]<- ((help_s41/num_indiv)^0.5)[4:dim(xtot[[relev_ind[1]]]$stox)[2]]
  return(curr_pop)
}