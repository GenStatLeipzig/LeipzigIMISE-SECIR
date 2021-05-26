#Integral Model
hook_jeeves_run_together_indiv_pop<-function(xtot,main_struct, opt_mod)
{
  writing_options <- list()
  if(main_struct$pop_modelling==1)#must be defined!!!
  {
    output_struct<-step_big_hook_jeeves_nonmem(xtot,main_struct, opt_mod)
    xtot <- output_struct$xtot # not params as in the matlab version!!!!!!!!
    writing_options$fval <- output_struct$fval 
    #output_struct$optim_options <- optim_options
    writing_options$intermediate=0;
    writing_options$pop_modelling=1;
    WriteHierarchicalResults(output_struct,xtot,results_file = main_struct1$path_results ,main_struct,optim_options,writing_options)
    
  }else{
    if(main_struct$subjects_subgroup==0){
      seqrelev=1:main_struct$numindiv;#later redefine
    }else{
      seqrelev=main_struct$relevant_subjects;
    }
    
    for (ind_indiv in 1:main_struct$numindiv)
    {
      if(ind_indiv%in%seqrelev)
      {
        main_struct$curr_ind_indiv=ind_indiv;
        writing_options$curr_ind_indiv=ind_indiv;
        output_struct<-step_big_hook_jeeves_nonmem(xtot,main_struct, opt_mod)
      }
    }
    
  }
}

step_big_hook_jeeves_nonmem<-function(xtot,main_struct, opt_mod)# the matlab name is: hook_jeeves_TPO_together_indiv_pop
{
  main_struct<-UpdateMainstuctforParametersEstimations(main_struct)
  
  if(main_struct$subjects_subgroup==0)
  {
    seqrelev <- 1:main_struct$numindiv;#later redefine
  }else{
    seqrelev <- main_struct$relevant_subjects;
  }
    ind_loc=0;
    out1<-InitOutputArrayss (xtot=phiopt,main_struct=main_struct,optim_options  = main_struct$optim_options)
    optim_options<-main_struct$optim_options
    err_loc_step<-vector(length=optim_options$num_it_u)
    while(ind_loc<optim_options$num_it)
    {#&(err_min<400))#for ind_loc=1:num_it
      ind_loc=ind_loc+1;
      hjloop<-step_me_hook_jeeves_nonmem(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,seqrelev=seqrelev)
      
      xtot <- hjloop$xtot
      optim_options <- hjloop$optim_options
      main_struct$optim_options<-optim_options 
      out1<- UpdateOutputArrays(ind_loc,output_struct=out1,xtot=xtot,main_struct=main_struct,optim_options  = optim_options,writing_options=list()) 
      WriteHierarchicalResults(output_struct=out1,xtot=xtot,results_file = main_struct$resultspath ,main_struct=main_struct,optim_options  = optim_options,writing_options=list())
      err_loc_step[ind_loc]=optim_options$nLL_tot_min
    }#!!! for  while(ind_loc<optim_options.num_it)
     
     #steps with stopping criteria:
     #1.  goal function is not too small;
     #2. the step is less than maximal one
     #3. The difference between GF of this step and that of 4 steps before is substantial
     num_steps_to_compare_back=4;
     if (num_steps_to_compare_back<=ind_loc)
     {
        num_steps_to_compare_back = ind_loc-1;
     }
     
     if (ind_loc>num_steps_to_compare_back) # This line added by OA, May 30, 2011
     {#while((ind_loc<param0.num_it_u)&((abs(err_loc_step(ind_loc-4)-err_loc_step[ind_loc])>param_.adj0)|(abs(err_loc_step(ind_loc-4)-err_loc_step[ind_loc])>param0.adj*abs(err_loc_step[ind_loc]))))
     #(err_loc_step[ind_loc]>optim_options$ref_err_min)&
       while(((ind_loc<optim_options$num_it_u)&(abs(err_loc_step[ind_loc-num_steps_to_compare_back]-err_loc_step[ind_loc])>optim_options$adj*abs(err_loc_step[ind_loc]))))
       {
         ind_loc=ind_loc+1;
        
         hjloop<-step_me_hook_jeeves_nonmem(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,seqrelev=seqrelev)
         
         xtot <- hjloop$xtot
         optim_options <- hjloop$optim_options
         main_struct$optim_options<-optim_options 
         out1<- UpdateOutputArrays(ind_loc,output_struct=out1,xtot=xtot,main_struct=main_struct,optim_options  = optim_options,writing_options=list()) 
         WriteHierarchicalResults(output_struct=out1,xtot=xtot,results_file = main_struct$resultspath ,main_struct=main_struct,optim_options  = optim_options,writing_options=list())
         err_loc_step[ind_loc]=optim_options$nLL_tot_min
       }
     } # This line added by OA, May 30, 2011
   output_struct<-list()
   output_struct$xtot <- xtot# not params as in the matlab version!!!!!!!!
   output_struct$fval <-  optim_options$nLL_tot_min
   output_struct$optim_options <- optim_options
   }


step_sm_hook_jeeves_nonmem<-function(optim_options,ind_par_loop,xtot,main_struct,label_opt, opt_mod,ind_indiv)
{
  #optim_options is a field of main_struct!!!!
  output_struct<-list()
  #optim_options<-main_struct$optim_options
  nLL_tot_min <- optim_options$nLL_tot_min;
  nLL_indiv_min <- optim_options$nLL_indiv_min;
  x0 <- xtot;
  if(label_opt=='pop')
    {
    eps_vect <- optim_options$pop$eps_vect;
    par_a <- optim_options$pop$LB;
    par_b <- optim_options$pop$UB;
    ind_indiv <- -10;#virtual input...
  }else if(label_opt=='indiv')
  {
    eps_vect <- optim_options$indiv$eps_vect[[ind_indiv]];
    par_a <- optim_options$indiv$LB[[ind_indiv]];
    par_b <- optim_options$indiv$UB[[ind_indiv]];
  }
  else if(label_opt=='resid')
  {  
    eps_vect <- optim_options$resid$eps_vect[[ind_indiv]];
    par_a <- optim_options$resid$LB[[ind_indiv]];
    par_b <- optim_options$resid$UB[[ind_indiv]];
  } 
  
  par_unpert <- UpdateField(x0,label_opt,ind_par_loop,'vect',x0,ind_indiv);
  par_pos <- par_unpert+eps_vect[ind_par_loop];
  par_neg <- par_unpert-eps_vect[ind_par_loop];
  
  #whether positive and negative perturbations falls into allowed boundaries
  #[par_a,par_b]
  if((par_pos<par_b[ind_par_loop]) & (par_neg>par_a[ind_par_loop]))
  {    #if positive and negative perturbations lay in allowed interval
    #positive perturbation in the ind_par_loop-th component:
    #x0[ind_par_loop] <- par_pos;
    x0 <- UpdateField(par_pos,label_opt,ind_par_loop,'str',x0,ind_indiv);
    #cand_sim <- evalindivLLadv1(x0,ind_indiv,main_struct,[],opt_mod);
    #cand_sim <- funevalTogetherLLadvThromboResest(x0,main_struct,opt_mod);
    #if(label_opt == 'pop'){#!!!!27.11.2018, no need to separate these cases!!!!!
      cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
                  #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min);
    #(xtot,main_struct,label_opt, opt_mod,nLL_indiv,ind_indiv)
   # }else{
   #   cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
                  #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min,ind_indiv);
   # }
    nLL_tot_pos <-  cand_sim$nLL_tot;#feval(fitfun,x0);
    nLL_indiv_pos <-  cand_sim$nLL_indiv;
    #     nLL_pos  <-  feval(fitfun,x0);
    
    if(optim_options$err_orient<0)
    {
      nLL_tot_pos <- 1/(nLL_tot_pos+1e-10);
    }
    #negative perturbation in the ind_par_loop-th component:
    x0 <- UpdateField(par_neg,label_opt,ind_par_loop,'str',x0,ind_indiv);
    #if(label_opt == 'pop')
    #{#!!!!27.11.2018, no need to separate these cases!!!!!
      cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
                #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min);
    #}else
    #{
    #  cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
                  #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min,ind_indiv);
    #}
    nLL_tot_neg <-  cand_sim$nLL_tot;#feval(fitfun,x0);
    nLL_indiv_neg <-  cand_sim$nLL_indiv;
    # nLL_neg <-  feval(fitfun,x0);
    if(optim_options$err_orient<0)
    {
      nLL_tot_neg <- 1/(nLL_tot_neg+1e-10);
    }
    #par_pert is the minimal of positive and negative perturbations:
    if(nLL_tot_pos<nLL_tot_neg)
    {
      par_pert <- par_pos;
      nLL_tot_pert <- nLL_tot_pos;
      nLL_indiv_pert <- nLL_indiv_pos;
    }else
    {
      par_pert <- par_neg;
      nLL_tot_pert <- nLL_tot_neg;
      nLL_indiv_pert <- nLL_indiv_neg;
    }
  }else
  {
    if(par_pos<par_b[ind_par_loop])#if positive perturbation is in allowed boundary:
    {  
      par_pert <- par_pos;
      x0 <- UpdateField(par_pos,label_opt,ind_par_loop,'str',x0,ind_indiv);
      #if(label_opt == 'pop')#!!!!27.11.2018, no need to separate these cases!!!!!
     # {
        cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
        #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min);
     # }else{
    #    cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
        #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min,ind_indiv);
    #  }
      nLL_tot_pos <-  cand_sim$nLL_tot;#feval(fitfun,x0);
      nLL_indiv_pos <-  cand_sim$nLL_indiv;
      if(optim_options$err_orient<0)
      {
        nLL_tot_pos <- 1/(nLL_tot_pos+1e-10);
      }
      nLL_tot_pert <- nLL_tot_pos;
      nLL_indiv_pert <- nLL_indiv_pos;
    }else
    {
      if(par_neg>par_a[ind_par_loop]) #if negative perturbation in allowed boundary:
      {
        par_pert <- par_neg;
        x0 <- UpdateField(par_neg,label_opt,ind_par_loop,'str',x0,ind_indiv);
        #if(label_opt == 'pop')#!!!!27.11.2018, no need to separate these cases!!!!!
        #{
          cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
          #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min);
        #}else{
        #  cand_sim <- GoalFuncComp(xtot=x0,main_struct,opt_mod,label_opt,nLL_indiv = nLL_indiv_min,ind_indiv1 = ind_indiv)
          #feval(fitfun,x0,main_struct,label_opt,opt_mod,nLL_indiv_min,ind_indiv);
        #}
        nLL_tot_neg <-  cand_sim$nLL_tot;#feval(fitfun,x0);
        nLL_indiv_neg <-  cand_sim$nLL_indiv;
        if(optim_options$err_orient<0)
        {
          nLL_tot_neg <- 1/(nLL_tot_neg+1e-10);
        }
        nLL_tot_pert <- nLL_tot_neg;
        nLL_indiv_pert <- nLL_indiv_neg;
      }else
      {   #both perturbations are out of the allowed interval
        if(eps_vect[ind_par_loop]>optim_options$eps_min)
        {
          nLL_tot_pert <- nLL_tot_min+1.0;#the result is bad
          #x0[ind_par_loop] <-   par_unpert;
          x0 <- UpdateField(par_unpert,label_opt,ind_par_loop,'str',x0,ind_indiv);
        }
      }
    }
  }
  
  if(nLL_tot_min>nLL_tot_pert)
  {
    #if smallest perturbed goal function  is lower than the unperturbed
    x0 <- UpdateField(par_pert,label_opt,ind_par_loop,'str',x0,ind_indiv);
    eps_vect[ind_par_loop] <-  eps_vect[ind_par_loop]*optim_options$update_1;
    nLL_tot_min <- nLL_tot_pert;
    nLL_indiv_min <- nLL_indiv_pert;
  }else{
    x0 <- UpdateField(par_unpert,label_opt,ind_par_loop,'str',x0,ind_indiv);
    if((nLL_tot_min<=nLL_tot_pert) & (eps_vect[ind_par_loop]>optim_options$eps_min))
    {
      eps_vect[ind_par_loop] <-  eps_vect[ind_par_loop]/optim_options$update_2;
    }
  }
  xtot <- x0;
  optim_options$nLL_indiv_min <- nLL_indiv_min;
  optim_options$nLL_tot_min <- nLL_tot_min;
  'small steps:'
  label_opt
  print(c(ind_par_loop,nLL_tot_min/1000))#print(c(ind_par_loop, nLL_tot_min/1000,nLL_indiv_min/1000)) !!! not necessary to print!!!
  
  if(label_opt== 'pop')
  {
    optim_options$pop$eps_vect <- eps_vect;
  }else if(label_opt == 'indiv'){
    optim_options$indiv$eps_vect[[ind_indiv]] <- eps_vect; 
  }else if(label_opt== 'resid'){
    optim_options$resid$eps_vect[[ind_indiv]] <- eps_vect;
  }
  
  output_struct$xtot <- x0
  output_struct$optim_options <- optim_options
  ###output:
  output_struct
  
}

step_me_hook_jeeves_nonmem<-function(xtot,main_struct,opt_mod,seqrelev)
{
  optim_options<-main_struct$optim_options
  if(main_struct$pop_modelling==0)#Instead of external options!!!!
  {
    ind_indiv <- main_struct.curr_ind_indiv
    seqrelev<-c(ind_indiv)
  }
  if(main_struct$numpoppar>0)
  {
    label_opt='pop'
    cand_sim<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                           nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=1)#für pop Option kann ind_indiv1 beliebeig sein
    optim_options$nLL_tot_min <- cand_sim$nLL_tot
    optim_options$nLL_indiv_min <- cand_sim$nLL_indiv#!!!!Corrected/added 23.11.2018
    #optim_options$nLL_output_tot_min <- cand_sim$nLL_tot#!!!!Corrected/added 28.11.2018
    for (ind_par_loop in 1:main_struct$numpoppar)
    {#length(xtot)
      main_struct$relev_input_indiv <- main_struct$helplistid[main_struct$subs_popisubs[,ind_par_loop]==1]
      if(ind_par_loop>14)
      {
        ind_par_loop
      }
      temp_struct <- step_sm_hook_jeeves_nonmem(optim_options,ind_par_loop,xtot,main_struct,label_opt,opt_mod,ind_indiv = -1);#in diesem Fall ind_indiv kann beliebig sein
      xtot <- temp_struct$xtot
      optim_options <- temp_struct$optim_options
      #step_sm_hook_jeeves_nonmem(optim_options,ind_par_loop,xtot,main_struct,label_opt, model_opt,ind_indiv)
      #step_sm_hook_jeeves_nonmem(optim_options,ind_par_loop,xtot,ind_indiv,main_struct,label_opt,model_opt)
      
      #if (optim_options$display==1)!!!! not yet,may be problematic for big sets
      #{
        ##print(c(ind_loc, ind_par_loop,optim_options$nLL_tot_min))###in the big loop
        #print(c(optim_options$nLL_indiv_min))!!!! not yet,may be problematic for big sets
      #}
    }
    #optim_options$nLL_output_tot_min <- optim_options$nLL_tot_min#!!!!Corrected/added 28.11.2018
  }
  if(main_struct$numpoppar==0)
  {
    optim_options$nLL_indiv_min=rep.int(x=0,times=main_struct$numindiv)#initialiation, when no population parameters!!!!
    optim_options$nLL_tot_min <- 0
    #optim_options$nLL_output_tot_min <- 0#!!!!Corrected/added 28.11.2018
  }
  label_opt <- 'indiv'
  for (ind_indiv in 1:main_struct$numindiv)
  {
    if(ind_indiv%in%seqrelev)
    {
      if(main_struct$numpoppar==0)
      {
        cand_sim<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                               nLL_indiv=optim_options$nLL_indiv_min,ind_indiv1=ind_indiv)#für pop Option kann ind_indiv1 beliebeig sein, not rep.int(x=0,times=main_struct$numindiv)!!!
        optim_options$nLL_tot_min<- cand_sim$nLL_tot
        optim_options$nLL_indiv_min[ind_indiv] <- cand_sim$nLL_indiv[ind_indiv]#!!!! very important for the individual cases when numpoppar=0
      }
      for (ind_par_loop in 1:main_struct$num_estim_indiv[ind_indiv])
      {#length(xtot)
        temp_struct <- step_sm_hook_jeeves_nonmem(optim_options,ind_par_loop,xtot,main_struct,label_opt,opt_mod,ind_indiv);
        xtot <- temp_struct$xtot
        optim_options <- temp_struct$optim_options
        
        #step_sm_hook_jeeves_nonmem(optim_options,ind_par_loop,xtot,ind_indiv,main_struct,label_opt,model_opt)
        if (optim_options$display==1)
        {
          #print(c(ind_loc, ind_par_loop,optim_options$nLL_tot_min))###in the big loop
          #print(c(optim_options$nLL_indiv_min))#!!!! not yet
        }
      }
    }
  }
  label_opt='resid'
  for (ind_indiv in 1:main_struct$numindiv)
  {
    if(ind_indiv %in% seqrelev)
    {
      #if(sum(strcmp(fields(main_struct),'act_err'))>0)
      if(main_struct$numpoppar==0)
      {
        cand_sim<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                               nLL_indiv=optim_options$nLL_indiv_min ,ind_indiv1=ind_indiv)#, not rep.int(x=0,times=main_struct$numindiv)!!!!für pop Option kann ind_indiv1 beliebeig sein
        optim_options$nLL_tot_min<- cand_sim$nLL_tot
        optim_options$nLL_indiv_min[ind_indiv] <- cand_sim$nLL_indiv[ind_indiv]#!!!! very important for the individual cases when numpoppar=0
      }
      relevant_output_names <- main_struct$outcomesnames[main_struct$act_err==1];
      ind_all_resid <- 1:main_struct$numtypesresidinit;
      relevant_output_ind <- ind_all_resid [main_struct$act_err==1]
      optimized_output_ind <- ind_all_resid[main_struct$opt_err==1]
      for (ind_par_loop in 1:dim(xtot$resid)[2])#size(main_struct$residinit,2)#length(xtot)
      {
        if((sum(relevant_output_ind[ind_par_loop]==main_struct$YTYPE[[ind_indiv]]$listnum)>0)&(sum(relevant_output_ind[ind_par_loop]==optimized_output_ind)>0))
        {
          print(c(ind_par_loop,relevant_output_ind[ind_par_loop]))
          relevant_output_names[ind_par_loop]
          temp_struct <- step_sm_hook_jeeves_nonmem(optim_options,ind_par_loop,xtot,main_struct,label_opt,opt_mod,ind_indiv);
          xtot <- temp_struct$xtot
          optim_options <- temp_struct$optim_options
          
          #step_sm_hook_jeeves_nonmem(optim_options,ind_par_loop,xtot,ind_indiv,main_struct,label_opt,opt_mod)
          #if (optim_options$display==1)#strange, gives errors!!!!!!!!
          #{
          # print(c(ind_loc, ind_par_loop,optim_options$nLL_tot_min))####in the big loop
          #  print(c(optim_options$nLL_indiv_min))
         # }
        }
      }
      
    }
  }
  #optim_options$nLL_output_tot_min <- sum(optim_options$nLL_indiv_min)#!!!!Corrected/added 28.11.2018
  xtot<-EstimDistribIndivPar(xtot,main_struct)
  output_struct<-list()
  output_struct$xtot <- xtot 
  output_struct$optim_options <- optim_options
  ###output:
  output_struct
}
