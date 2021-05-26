SimulCppCovidModel<-function(param_,xtot,ind_indiv,main_struct)#SimulIntegralModelShort
{
  
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
  
  
  
  if(main_struct$treatpardir==1)
  {  
    compl_parv[compl_parl-1]<- main_struct$num_locks
    compl_parv[compl_parl]<- main_struct$num_unlocks
    Sloc <- solvecovidmodel(  parval = compl_parv, parnames = compl_parn, numpar = compl_parl, u0=xtot$init,
                   num_steps = xtot$num_steps, numvar= length(xtot$init), num_locks=main_struct$num_locks, num_unlocks=main_struct$num_unlocks,
                   compnames =  main_struct$compartments[[ind_indiv]]$names,
                   compind = main_struct$compartments[[ind_indiv]]$compind-1,#!!!not hear, bolus is updated in R comp_ind must be -1 shifted, c++, -1
                   compnum = main_struct$compartments[[ind_indiv]]$compnum, 
                   comppar = c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),#=c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),
                   crit_exp = compl_parv[main_struct$crit_exp_ind], death_exp = compl_parv[main_struct$death_exp_ind],
                   death_int_round = round(compl_parv[main_struct$death_int_ind]),
                   crit_int_round =round(compl_parv[main_struct$crit_int_ind]),
                   blockeff_vec = compl_parv[main_struct$blockeff_ind], unlockeff_vec = compl_parv[main_struct$unlockeff_ind],
                   Lock_vec = compl_parv[main_struct$all_locks_ind]-1, Unlock_vec = compl_parv[main_struct$all_unlocks_ind]-1,
                   treatpardir=main_struct$treatpardir)#c++ starts with
  
  }else if(main_struct$model_opt==31)
  {
    compl_parv[compl_parl-1]<- param_$num_locks
    compl_parv[compl_parl]<- param_$num_unlocks
    Sloc <- solvecovidmodelMuKI(  parval = compl_parv, parnames = compl_parn, numpar = compl_parl, u0=xtot$init,
                                num_steps = xtot$num_steps, numvar= length(xtot$init), num_locks=param_$num_locks, num_unlocks=param_$num_unlocks,
                                compnames =  main_struct$compartments[[ind_indiv]]$names,
                                compind = main_struct$compartments[[ind_indiv]]$compind-1,#!!!not hear, bolus is updated in R comp_ind must be -1 shifted, c++, -1
                                compnum = main_struct$compartments[[ind_indiv]]$compnum, 
                                comppar = c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),#=c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),
                                crit_exp = compl_parv[main_struct$crit_exp_ind], death_exp = compl_parv[main_struct$death_exp_ind],
                                death_int_round = round(compl_parv[main_struct$death_date_ind]),
                                crit_int_round =round(compl_parv[main_struct$crit_date_ind]),
                                treat_vec = compl_parv[main_struct$treat_par_ind],
                                LockUnlock_vec = round(compl_parv[main_struct$treat_date_ind]-1),
                                treatpardir=main_struct$treatpardir)#c++ starts with
  }else
  {
    compl_parv[compl_parl-1]<- param_$num_locks
    compl_parv[compl_parl]<- param_$num_unlocks
    Sloc <- solvecovidmodelKI(  parval = compl_parv, parnames = compl_parn, numpar = compl_parl, u0=xtot$init,
                              num_steps = xtot$num_steps, numvar= length(xtot$init), num_locks=param_$num_locks, num_unlocks=param_$num_unlocks,
                              compnames =  main_struct$compartments[[ind_indiv]]$names,
                              compind = main_struct$compartments[[ind_indiv]]$compind-1,#!!!not hear, bolus is updated in R comp_ind must be -1 shifted, c++, -1
                              compnum = main_struct$compartments[[ind_indiv]]$compnum, 
                              comppar = c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),#=c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),
                              crit_exp = compl_parv[main_struct$crit_exp_ind], death_exp = compl_parv[main_struct$death_exp_ind],
                              death_int_round = round(compl_parv[main_struct$death_date_ind]),
                              crit_int_round =round(compl_parv[main_struct$crit_date_ind]),
                              treat_vec = compl_parv[main_struct$treat_par_ind],
                              LockUnlock_vec = round(compl_parv[main_struct$treat_date_ind]-1),
                              treatpardir=main_struct$treatpardir)#c++ starts with
  }
  return(Sloc)
}
SimulCppCovidModelIndiv<-function(param_,xtot,ind_indiv,main_struct)#SimulIntegralModelShort
{
  
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
  
  numindivpar <- main_struct$numindivpar1
  numpoppar <- main_struct$numpoppar
  if(numpoppar>0)
  {
    popparamnames <- main_struct$names_pop_opt_parameters
    popparameterspsi <- xtot$phiopt$pop
    poptransform_label_LB_UB <- main_struct$poptransform_label_LB_UB
  }else{
    popparamnames <- c('a','b')
    popparameterspsi <- rep.int(x=0,times =2)
    poptransform_label_LB_UB <- array(data = 0, dim =c(2,2))
  }
  if(numindivpar>0)
  {
    indivparamnames <- main_struct$names_indiv_opt_parameters[[ind_indiv]]
    subs_indivi_num <- main_struct0$subs_indivi_num[[ind_indiv]]
    indivparameterspsi <- xtot$phiopt$indiv[[ind_indiv]]
    indivtransform_label_LB_UB <- main_struct$indivtransform_label_LB_UB[[ind_indiv]]
  }else{
    indivparamnames <- c('a','b')
    subs_indivi_num <- c(0,0)
    indivparameterspsi <- rep.int(x=0,times =2)
    indivtransform_label_LB_UB <- array(data = 0, dim =c(2,2))
  }
  Sloc <- solvecovidmodelindiv(  parval = compl_parv, parnames = compl_parn, numpar = compl_parl, u0=xtot$init,
        num_steps = xtot$num_steps, numvar= length(xtot$init), num_locks=main_struct$num_locks, num_unlocks=main_struct$num_unlocks,
        compnames =  main_struct$compartments[[ind_indiv]]$names,
        compind = main_struct$compartments[[ind_indiv]]$compind-1,# -1!!!not hear, bolus is updated in R comp_ind must be -1 shifted, c++, -1
        compnum = main_struct$compartments[[ind_indiv]]$compnum, 
        comppar = c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),#=c(main_struct$compartments[[ind_indiv]]$NumOfCompartments,main_struct$compartments[[ind_indiv]]$NumofVariables),
        crit_exp_ind = main_struct$crit_exp_ind-1, death_exp_ind = main_struct$death_exp_ind-1,
        death_int_ind = main_struct$death_int_ind-1, crit_int_ind = main_struct$crit_int_ind-1,
        blockeff_ind = main_struct$blockeff_ind-1, unlockeff_ind = main_struct$unlockeff_ind-1,
        all_locks_ind = main_struct$all_locks_ind-1, all_unlocks_ind = main_struct$all_unlocks_ind-1,
        numpoppar = numpoppar, numindivpar = numindivpar,  
        popparamnames=popparamnames, ind_opt_parameters = main_struct$ind_opt_parameters , popparameterspsi = popparameterspsi,
        indivparamnames = indivparamnames, subs_indivi_num = subs_indivi_num, indivparameterspsi = indivparameterspsi,
        poptransform_label_LB_UB= poptransform_label_LB_UB,
        indivtransform_label_LB_UB =indivtransform_label_LB_UB)#c++ starts with
  
  return(Sloc)
}