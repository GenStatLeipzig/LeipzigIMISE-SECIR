IndivUpdateParam<-function(x,ind_indiv,main_struct,param_)#xtot$indiv  par_indiv 02.04.2020
{  
  
  if(main_struct$numpoppar>0)#if there exist population parameters!!!!
  {
    ind_opt_parameters <- allpar[main_struct$subs_pop_opt>0]
    transformed_par<-vector(length = main_struct$numpoppar)
    for(ind_par in 1:main_struct$numpoppar)#main_struct$numpoppar
    {
      namei <- main_struct$names_pop_opt_parameters[ind_par]#main_struct$names[ind_opt_parameters][ind_par]
      transformed_par[ind_par] <- TransformParameters(x$pop[ind_par],main_struct,namei)
      param_[[namei]]<-transformed_par[ind_par] 
    }
  }#####individual parameters:
  transformed_par<-vector(length = main_struct$num_estim_indiv[ind_indiv])
  if(main_struct$numindivpar1>0)
  {  
    for(ind_par in 1:main_struct$num_estim_indiv[ind_indiv])
    {
      # print(ind_par)
      namei <- main_struct$names_indiv_opt_parameters[[ind_indiv]][ind_par]#main_struct$names[main_struct$subs_indivi_num[[ind_indiv]]][ind_par]
      transformed_par[ind_par] <- TransformParameters(x$indiv[[ind_indiv]][ind_par],main_struct,namei)
      param_[[namei]]<-transformed_par[ind_par] 
    }
  }
  
  
  return(param_)
}
ConstructParam<-function(x,ind_indiv,main_struct)#xtot$indiv  par_indiv 02.04.2020
{  
  param_<-list()
  for(ind_par in 1:main_struct$numallpar)
  {
    param_[[main_struct$names[ind_par]]]<-main_struct$init_mean[ind_par]
  }
  param_$label_poly_age <- main_struct$label_poly_age
  if(main_struct$label_poly_age>0)
  {
    param_$age_spec_id <- main_struct$age_spec_id
    param_$age_spec_param_num <- main_struct$age_spec_param_num
    
    
    param_$age_spec_param_names  <- main_struct$age_spec_param_names
    param_$age_spec_dim <- main_struct$age_spec_dim
    
    if(ind_indiv%in%main_struct$age_spec_id)
    {
      dim(main_struct$age_spec_phi[[ind_indiv]])
      param_$age_spec_phi <-array(dim=main_struct$age_spec_dim)
      param_$age_spec_param <-array(dim=main_struct$age_spec_dim)
      param_$age_spec_phi[1,]  <- x$indiv[[ind_indiv]][1:main_struct$numindivpar1][main_struct$names_indiv_opt_parameters[[ind_indiv]]%in%main_struct$age_spec_param_names]#main_struct$age_spec_phi
      for(ind_agei in 2:main_struct$age_spec_dim[1])
      {
        param_$age_spec_phi[ind_agei,] <- x$indiv[[ind_indiv]][main_struct$numindivpar1+main_struct$age_spec_dim[2]*(ind_agei-2)+(1:main_struct$age_spec_dim[2])]
      }
      for(ind_agei in 1:main_struct$age_spec_dim[1])
      {
        for(ind_poly in 1:main_struct$age_spec_dim[2])
        {
          param_$age_spec_param[ind_agei,ind_poly]<- TransformParameters(x=param_$age_spec_phi[ind_agei,ind_poly],main_struct,
                                                                         name = main_struct$age_spec_param_names[ind_poly])
        }
      }
        
    }
    
    
  }
  allpar<-1:main_struct$numallpar
  #####population parameters:
  ind_opt_parameters <- allpar[main_struct$subs_pop_opt>0]#!!!! popi???
  
  if(main_struct$numpoppar>0)#if there exist population parameters!!!!
  {
    transformed_par<-vector(length = main_struct$numpoppar)
    for(ind_par in 1:main_struct$numpoppar)#main_struct$numpoppar
    {
      namei<-main_struct$names[ind_opt_parameters][ind_par]
      transformed_par[ind_par] <- TransformParameters(x$pop[ind_par],main_struct,namei)
      param_[[namei]]<-transformed_par[ind_par] 
    }
  }#####individual parameters:
  transformed_par<-vector(length = main_struct$num_estim_indiv[ind_indiv])
  if(main_struct$numindivpar1>0)
  {  
    for(ind_par in 1:main_struct$num_estim_indiv[ind_indiv])
    {
     # print(ind_par)
      namei<-main_struct$names[main_struct$subs_indivi_num[[ind_indiv]]][ind_par]
      transformed_par[ind_par] <- TransformParameters(x$indiv[[ind_indiv]][ind_par],main_struct,namei)
      param_[[namei]]<-transformed_par[ind_par] 
    }
  }
  ####no assignments for individual fixed values, they mast be equal to init_mean!!!!!
  #####fixed parameters:
  for(ind_par in 1:main_struct$numspecialpar)#main_struct$numpoppar
  {
    namei<-as.character(main_struct$specialpar$names[ind_par])
    param_[[namei]]<-main_struct$specialpar$val[ind_par] 
  }
  param_$DateStart <- min(main_struct$measurements[[ind_indiv]]$date)
  if(main_struct$label_treat_info==0)
  {
    for(ind_cov in 1:3)
    {
      if(paste('DateLock',ind_cov,sep = '')%in%names(main_struct$covariates))
      {
        param_[[paste('Lock',ind_cov,sep = '')]] <- as.numeric(main_struct$covariates[[paste('DateLock',ind_cov,sep = '')]][main_struct$covariates$ID==ind_indiv] - param_$DateStart)
      }
    }  
    #for(ind_cov in 1:3)
    #{
    if('Unlock'%in%names(main_struct$covariates))
    {
      param_[['Unlock']] <- as.numeric(main_struct$covariates[['Unlock']][main_struct$covariates$ID==ind_indiv] - param_$DateStart)
    }
    if('Unlock1'%in%names(main_struct$covariates))
    {
      param_[['Unlock1']] <- as.numeric(main_struct$covariates[['Unlock1']][main_struct$covariates$ID==ind_indiv] - param_$DateStart)
    }
    if('Unlock2'%in%names(main_struct$covariates))
    {
      param_[['Unlock2']] <- as.numeric(main_struct$covariates[['Unlock2']][main_struct$covariates$ID==ind_indiv] - param_$DateStart)
    }
  #}  
  }else if(main_struct$label_treat_info==1){
    for(ind_cov in 1:main_struct$num_unlocks)
    {
      param_[[paste('Unlock',ind_cov,sep='')]] <- as.numeric(main_struct$covariates[[paste('Unlock',ind_cov,sep='')]][main_struct$covariates$ID==ind_indiv] - param_$DateStart)
    }
    for(ind_cov in 1:main_struct$num_locks)
    {
      param_[[paste('Lock',ind_cov,sep = '')]] <- as.numeric(main_struct$covariates[[paste('DateLock',ind_cov,sep = '')]][main_struct$covariates$ID==ind_indiv] - param_$DateStart)
    }
  }
  param_$crit_prior_num <- length(grep(pattern='crit_exp',x=names(param_)))
  param_$death_prior_num <- length(grep(pattern='death_exp',x=names(param_)))
  
  param_$Sc0 <-  main_struct$covariates$Sc0[ind_indiv] 
  param_$speclength<- main_struct$speclength
  
  if('dates_treatments'%in%names(main_struct))
  {
    num_locks <- 0
    num_unlocks <- 0
    prev_val <- 1
    for(ind_tr in 1:main_struct$treat_par_num)
    {
      if(param_[[paste('treat',ind_tr,sep = '')]]>=prev_val)
      {
        num_unlocks<- num_unlocks+1
        param_[[paste('Unlock',num_unlocks,sep='')]]<-param_[[paste('date_treat',ind_tr,sep = '')]]
      }else{
        num_locks <- num_locks+1
        param_[[paste('Lock',num_locks,sep='')]]<-param_[[paste('date_treat',ind_tr,sep = '')]]
      }
      prev_val <- param_[[paste('treat',ind_tr,sep = '')]]
    }
    param_$num_locks <- num_locks
    param_$num_unlocks <- num_unlocks
  }
  if('const'%in%names(main_struct$indiv))
  {
    for(ind_names in colnames(main_struct$indiv$const))
    {
      param_[[ind_names]]<- main_struct$indiv$const[ind_indiv,ind_names]
    }
  }
  
  return(param_)
}
ConstructParamVect<-function(par_indiv,main_struct)#xtot$indiv
{  
  param_<-vector(length=main_struct$numallpar)
  for(ind_par in 1:main_struct$numallpar)
  {
    if(main_struct$subs_indiv[ind_par]==0)
    {
      param_[ind_par]<-main_struct$init_mean[ind_par]
    }else{
      good_ind<-main_struct$subs_indiv_num==ind_par
      param_[ind_par]<-TransformParameters(par_indiv[good_ind],main_struct,main_struct$names[main_struct$subs_indiv_num][good_ind])
    }
    
  }
  param_
}
##################
ConstructParamVectDir<-function(par_indiv)#xtot$indiv
{  
  param_<-vector(length=length(names(par_indiv)))
  for(ind_par in 1:length(names(par_indiv)))
  {
    #print(ind_par)
    param_[ind_par]<-par_indiv[[ind_par]]
  }
  param_
}
####################
RevTransformParametersVect<-function(x,transform_label,LB,UB)
{
  y=vector(length=length(x))
  for (ind1  in 1:length(x))
  {
    y[ind1]<- switch (toString(transform_label), 
                      "0"  = {x[ind1]},
                      "1"  = {log(x[ind1])},
                      "2"  = {lt=(x[ind1]-LB)/(UB-LB);
                      log(lt/(1-lt))}
    )
    #print(c(ind1,y[ind1],toString(main_struct$transform_label[ind1]),toString(main_struct$transform_label[ind1])==0,toString(main_struct$transform_label[ind1])==1,toString(main_struct$transform_label[ind1])==2))
    #if(toString(main_struct$transform_label[ind1])==2){
    #  print(ind1,lt,log(lt/(1-lt)))
    #}
  }
  y
}
RevTransformParametersAll<-function(x,main_struct)
{
  y=vector(length=length(x))
  for (ind1  in 1:length(x))
  {
    y[ind1]<- switch (toString(main_struct$transform_label[ind1]), 
                      "0"  = {x[ind1]},
                      "1"  = {log(x[ind1])},
                      "2"  = {lt=(x[ind1]-main_struct$LB[ind1])/(main_struct$UB[ind1]-main_struct$LB[ind1]);
                      log(lt/(1-lt))}
                      )
    #print(c(ind1,y[ind1],toString(main_struct$transform_label[ind1]),toString(main_struct$transform_label[ind1])==0,toString(main_struct$transform_label[ind1])==1,toString(main_struct$transform_label[ind1])==2))
    #if(toString(main_struct$transform_label[ind1])==2){
    #  print(ind1,lt,log(lt/(1-lt)))
    #}
  }
  y
}
RevTransformParameters<-function(x,main_struct,name)
{
  ind1=which(main_struct$names==name)
  y<-switch (toString(main_struct$transform_label[ind1]), 
          "0"  = {x},
          "1"  = {log(x)},
          "2"  = {lt=(x-main_struct$LB[ind1])/(main_struct$UB[ind1]-main_struct$LB[ind1]);
          log(lt/(1-lt))}
  )
  y
}

############################33
TransformParametersVect<-function(x,transform_label,LB,UB)
{
  y=vector(length=length(x))
  for (ind1  in 1:length(x))
  {
    y[ind1]<- switch (toString(transform_label[ind1]), 
                      "0"  = {x[ind1]},
                      "1"  = {exp(x[ind1])},
                      "2"  = {LB[ind1]+exp(x[ind1])*(UB[ind1]-LB[ind1])/(1+exp(x[ind1]))}
                      
    )
    
  }
  y
}
TransformParametersAll<-function(x,main_struct)
{
  y=vector(length=length(x))
  for (ind1  in 1:length(x))
  {
    y[ind1]<- switch (toString(main_struct$transform_label[ind1]), 
                      "0"  = {x[ind1]},
                      "1"  = {exp(x[ind1])},
                      "2"  = {main_struct$LB[ind1]+exp(x[ind1])*(main_struct$UB[ind1]-main_struct$LB[ind1])/(1+exp(x[ind1]))}
                      
    )
    
  }
  y
}
TransformParameters<-function(x,main_struct,name)
{
  ind1=which(main_struct$names==name)
  yy<-switch (toString(main_struct$transform_label[ind1]), 
          "0"  = {x},
          "1"  = {exp(x)},
          "2"  = {main_struct$LB[ind1]+exp(x)*(main_struct$UB[ind1]-main_struct$LB[ind1])/(1+exp(x))}
  )
  yy
}

################################3
extractsol<-function(tall,t,y)
{
  nt<-length(t)
  mt<-max(t)
  ntall<-length(tall)
  yout<-matrix(nrow = nt, ncol = ncol(y) )
  t_ind<-1
  for(tall_ind in 1:ntall)
  {
    #print(tall[tall_ind])
    #print(t_ind)
    if((tall[tall_ind]>t[t_ind])&&(t_ind<nt)&&(tall[tall_ind]<=mt))
    {
      t_ind<-t_ind+1
    }else{
      if(tall[tall_ind]==t[t_ind]){
        yout[t_ind,]<-y[tall_ind,]
      }
    }
  }
  
  yout
}

################
extractrelevind<-function(tall,t)#only when equal 
{
  nt<-length(t)
  mt<-max(t)
  ntall<-length(tall)
  yout<-vector(length = nt)
  t_ind<-1
  for(tall_ind in 1:ntall)
  {
    #print(tall[tall_ind])
    #print(t_ind)
    if((tall[tall_ind]>t[t_ind])&&(t_ind<nt)&&(tall[tall_ind]<=mt))
    {
      t_ind<-t_ind+1
    }else{
      if(tall[tall_ind]==t[t_ind]){
          yout[t_ind]<-tall_ind
      }
    }
  }
  yout
}
#############################
extractrelevind<-function(tall,t)#only when equal 
{
  nt<-length(t)
  mt<-max(t)
  ntall<-length(tall)
  yout<-vector(length = nt)
  t_ind<-1
  for(tall_ind in 1:ntall)
  {
    #print(tall[tall_ind])
    #print(t_ind)
    if((tall[tall_ind]>t[t_ind])&&(t_ind<nt)&&(tall[tall_ind]<=mt))
    {
      t_ind<-t_ind+1
    }else{
      if(tall[tall_ind]==t[t_ind]){
        yout[t_ind]<-tall_ind
      }
    }
  }
  yout
}
#######
extractrelevindapprox<-function(tall,t,tol_)#only when equal 
{
  nt<-length(t)
  mt<-max(t)
  ntall<-length(tall)
  yout<-vector(length = nt)
  t_ind<-1
  for(tall_ind in 1:ntall)
  {
    #print(tall[tall_ind])
    #print(t_ind)
    if((abs(tall[tall_ind]-t[t_ind])<tol_)&&(t_ind<=nt)&&(tall[tall_ind]<=mt))#if((tall[tall_ind]>t[t_ind])&&(t_ind<nt)&&(tall[tall_ind]<=mt))
    {
      yout[t_ind]<-tall_ind
      t_ind<-t_ind+1
    }
  }
  yout
}
############
extractsolrelevindapprox<-function(tall,t_,ytall)#,tol_)#only when equal 
{
  if(!is.numeric(tall))#!!!Date format 02.04.2020
  { 
    yout<-vector(length = length(t_))
    tmin<-min(tall)
    tall<-TransformDates(tall,tmin)
    t_<-TransformDates(t_,tmin)
    for(t_ind in 1:length(t_))
    {
      tall_ind<- (t_[t_ind]==tall)
      yout[t_ind]<-  ytall[tall_ind]# ytall[tall[tall_ind]+1]????!!!09.04.2020
      
    }
  }else{  
  nt<-length(t_)
  mt<-max(t_)
  ntall<-length(tall)
  yout<-vector(length = nt)
  t_ind<-1
     
    if((max(tall)<max(t_))&(max(tall)>max(t_)-0.0000001))#a little numerical uncertainty
    {
      tall[tall==max(tall)] <- max(t_)
    }
    
    for(tall_ind in 1:(ntall-1))
    {
      #print(tall[tall_ind])
      #print(t_ind)
      #t[(t_ind:length(t))]==t[t_ind]
      if(((tall[tall_ind]<=t_[t_ind])&&(tall[tall_ind+1]>=t_[t_ind]))&&(t_ind<=nt)&&(tall[tall_ind]<=mt))#if((abs(tall[tall_ind]-t[t_ind])<tol_)&&(t_ind<=nt)&&(tall[tall_ind]<=mt))
      {
        deltx<-as.numeric(tall[tall_ind+1]-tall[tall_ind])
        delty<-(ytall[tall_ind+1]-ytall[tall_ind])
        yout[t_ind]<- ytall[tall_ind]+(t_[t_ind]-tall[tall_ind])*delty/deltx#tall_ind
        t_ind<-t_ind+1
      }
    }
    yout<-yout[1:(t_ind-1)]#????
  }  
  
  ###output:
  yout
}
extractsolrelevindapproxmatr1<-function(tall,t,ytall)#,tol_)#only when equal 
{
  if((max(tall)<max(t))&(max(tall)>max(t)-0.0000001))#a little numerical uncertainty
  {
    tall[tall==max(tall)] <- max(t)
  }
  nt<-length(t)
  mt<-max(t)
  ntall<-length(tall)
  yout<-array(dim = c(nt,dim(ytall)[2]))
  t_ind<-1
  t_help<- vector(length = nt)
  for(tall_ind in 1:(ntall-1))
  {
    #print(tall[tall_ind])
    #print(t_ind)
    if(((tall[tall_ind]<=t[t_ind])&&(tall[tall_ind+1]>=t[t_ind]))&&(t_ind<=nt)&&(tall[tall_ind]<=mt))#if((abs(tall[tall_ind]-t[t_ind])<tol_)&&(t_ind<=nt)&&(tall[tall_ind]<=mt))
    {
      deltx<-(tall[tall_ind+1]-tall[tall_ind])
      delty<-(ytall[tall_ind+1,]-ytall[tall_ind,])
      yout[t_ind,]<- ytall[tall_ind,]+(t[t_ind]-tall[tall_ind])*delty/deltx#tall_ind
      t_help[t_ind] <- tall[tall_ind]
      t_ind<-t_ind+1
    }
  }
  yout<-yout[1:(t_ind-1),]
  ###output:
  yout
}




extractsolrelevindapproxmatr<-function(tall,t,ytall)#,tol_)#only when equal 
{
  if((max(tall)<max(t))&(max(tall)>max(t)-0.0000001))#a little numerical uncertainty
  {
    tall[tall==max(tall)] <- max(t)
  }
  nt<-length(t)
  mt<-max(t)
  ntall<-length(tall)
  yout<-array(dim = c(nt,dim(ytall)[2]))
  #t_ind<-1
  tall_ind <- 1
  t_help<- vector(length = nt)
  for(t_ind in 1:nt)
  {
    #print(tall[tall_ind])
    #print(t_ind)
    good_ind<-(tall<=t[t_ind])
    tall_ind<-max((1:ntall)[good_ind])
    
    if((tall_ind+1)<=ntall)#if((abs(tall[tall_ind]-t[t_ind])<tol_)&&(t_ind<=nt)&&(tall[tall_ind]<=mt))
    {
      deltx<-(tall[tall_ind+1]-tall[tall_ind])
      delty<-(ytall[tall_ind+1,]-ytall[tall_ind,])
      yout[t_ind,]<- ytall[tall_ind,]+(t[t_ind]-tall[tall_ind])*delty/deltx#tall_ind
      
      
      t_help[t_ind] <- tall[tall_ind]
      tall_ind<-tall_ind+1
    }else{
      yout[t_ind,]<- ytall[tall_ind,]
    }
    if(sum(is.na(yout[t_ind,]))>0)
    {
      print(t_ind)
    }
  }
  #yout<-yout[1:(t_ind-1),]
  ###output:
  yout
}



############ Parsimony function:
ParsimParam<-function(param,compartments)#only when equal function(param,specialpar)
{
  
  if(sum(names(param)=="kMGB4")==0)
  {
    param$kMGB4 <-1
  }
  if(sum(compartments$names=="GCSFEndomem")>0)
  {
    param$labelegscsfendomem <- 1
  }else{
    param$labelegscsfendomem <- 0
  }
  if(sum(compartments$names=="Ksimem")>0)
  {
    param$ksimem_lab <- 1
  }else{
    param$ksimem_lab <- 0
  }
  
  
  if(sum(names(param)=="relcarbopacl")>0)
  {
    param$pdcarbo <- param$pdpacl*param$relcarbopacl 
  }
  ################Integral:
  rel_kdormrev_asnor=1/param$kdormfactor
  param$rel_kdormrev_asnor=rel_kdormrev_asnor
  if(param$kdormrevsep==0)
  {
    param$kdormrev=param$asnor*rel_kdormrev_asnor
  }
  #####Granulo:
  if(sum(names(param)=='DGCSFFP')>0)
  {
    if(param$DGCSFFP==1)
    {
      param$DGCSFP <- param$DGCSFF
    }
  }
  if(sum(names(param)=='relparam')>0)
  {
    if(param$relparam==1)
    {
      #param$TPGBminf <- param$TPGBminfrel*param$TPGBnorf
      #param$TPGBmaxf <- param$TPGBmaxfrel*param$TPGBnorf
  
      param$TG4minf <- param$TG4minfrel*param$TG4norf
      param$TG4maxf <- param$TG4maxfrel*param$TG4norf
  
      param$TG5minf <- param$TG5minfrel*param$TG5norf
      param$TG5maxf <- param$TG5maxfrel*param$TG5norf
  
      param$TG6minf <- param$TG6minfrel*param$TG6norf
      param$TG6maxf <- param$TG6maxfrel*param$TG6norf
  
      param$AG6minf <- param$AG6minfrel*param$AG6norf
  
    }
  }
  ###Omitted:param$alphN=real(param$alphN);!!!????
  param <- PKParamPharmODe(param,compartments)##!!!!!????
  if(param$Tproplsingle==1)
  {
    param$Tpropl2 <- param$Tpropl
  }
  rel_kdormrev_asnor <- 1/param$kdormfactor
  param$rel_kdormrev_asnor <- rel_kdormrev_asnor
  if(param$kdormrevsep==0)
  {
    param$kdormrev <- param$asnor*rel_kdormrev_asnor
  }
  
  #print(tall[tall_ind])
  #print(t_ind)
  pars_option<-param$p_1_pars#specialpar$val[specialpar$names=="p_1_pars"]
  if(pars_option==1)
  {
    param$p_16_1_norm <- param$p_8_1_norm
    param$p_32_1_norm <- param$p_8_1_norm
    param$p_16_2_norm <- param$p_8_2_norm
  }else{
    if(pars_option==2){
      param$p_16_1_norm <- param$p_8_1_norm
      param$p_32_1_norm <- param$p_8_1_norm
      param$p_16_2_norm <- param$p_8_2_norm
      param$p_32_2_norm <- param$p_8_2_norm
    }else{
      if(pars_option==3){
        param$p_16_1_norm <- param$p_8_1_norm
        param$p_32_1_norm <- param$p_8_1_norm
        param$Trev_dorm32_norm <- param$Trev_dorm16_norm
        param$p_16_2_norm <- param$p_8_2_norm
        param$p_32_2_norm <- param$p_8_2_norm
      }
    }
    
  }
  
  if(sum(names(param)=='rgammap')>0)
  {
    param$gammap <- param$rgammap*param$bacm
  }
  if(sum(names(param)=='ksegcircnor.f')==0)
  {
    param$ksegcircnor.f <- param$kbandcircnor.f
  }
  if(sum(names(param)=='ksegcircmax.f')==0)
  {
    param$ksegcircmax.f <- param$kbandcircmax.f
  }
  if(sum(names(param)=='rwplc')>0)
  {
    param$wplc <- param$rwplc*param$km
  }
  if(sum(names(param)=='cssuppsame')>0)
  {
    if(param$cssuppsame==1)
    {
      param$cssuppmkc <- param$cssupp
    }
  
  }        
  ######
  if(sum(names(param)=='gra_t_pars')>0)
  {
    if(param$gra_t_pars==1)
    {#switch(param$gra_t_pars, #if
      #1= {param$btg5f = param$btg4f;
      param$btg6f <- param$btg4f
      param$TG5minfrel <- param$TG4minfrel
      param$TG6minfrel <- param$TG4minfrel
      param$TG5maxfrel <- param$TG4maxfrel
      param$TG6maxfrel <- param$TG4maxfrel
      param$ksegcircb.f <- param$kbandcircb.f
    }
    else if (param$gra_t_pars==11)
    {#mostly used
      #                 param$btg5f = param$btg4f;
      #                 param$btg6f = param$btg4f;
      ##11= {param$TG5minfrel = param$TG4minfrel;
      param$TG6minfrel <- param$TG4minfrel
      param$TG5maxfrel <- param$TG4maxfrel
      param$TG6maxfrel <- param$TG4maxfrel
      param$ksegcircb.f <- param$kbandcircb.f
    }
    else if(param$gra_t_pars==12)
    {#Better!!
      ##12= {param$TG5minfrel = param$TG4minfrel;
      param$TG6minfrel <- param$TG4minfrel
      param$TG5maxfrel <- param$TG4maxfrel
      param$TG6maxfrel <- param$TG4maxfrel
    }
    else if(param$gra_t_pars==13)
    {#Better!!
      ##12= {param$TG5minfrel = param$TG4minfrel;
      param$TG6minfrel <- param$TG4minfrel
      param$TG5maxfrel <- param$TG4maxfrel
      param$TG6maxfrel <- param$TG4maxfrel
      param$TG5maxfrel <- param$TG4maxfrel
      param$btg6f <- param$btg5f
    }
      #param$ksegcircb.f = param$kbandcircb.f;
    else if(param$gra_t_pars==2)
    {
      ##2= {param$btg5f = param$btg4f;
      param$btg6f <- param$btg4f
      param$kbandcircb.f <- param$btg4f;
      param$TG5minfrel <- param$TG4minfrel
      param$TG6minfrel <- param$TG4minfrel
      param$TG5maxfrel <- param$TG4maxfrel
      param$TG6maxfrel <- param$TG4maxfrel
      #param$TPGBminfrel <- param$TCGminfrel
      #param$TPGBmaxfrel <- param$TCGmaxfrel
      param$ksegcircb.f <- param$kbandcircb.f
    }
    else if(param$gra_t_pars==3)
    {
      ##3={param$btg5f = param$btg4f;
      param$btg6f <- param$btg4f
      #param$kbandcircb.f = param$btg4f;
      param$TG5minfrel <- param$TG4minfrel
      param$TG6minfrel <- param$TG4minfrel
      param$TG5maxfrel <- param$TG4maxfrel
      param$TG6maxfrel <- param$TG4maxfrel
     # param$TPGBminfrel <- param$TCGminfrel
      #param$TPGBmaxfrel <- param$TCGmaxfrel;
      #param$ksegcircb.f <- param$kbandcircb.f; 
    } 
    else if(param$gra_t_pars==4)
    {
      ##4={param$btg5f = param$btg4f;
      param$btg6f <- param$btg4f;
      param$TG5minfrel <- param$TG4minfrel;
      param$TG6minfrel <- param$TG4minfrel;
      param$TG5maxfrel <- param$TG4maxfrel;
      param$TG6maxfrel <- param$TG4maxfrel;
      #param$ksegcircb.f = param$kbandcircb.f;
    }
    param$factmech = param$blood_vol*(3/2)*230*(10^9)*30
    
  }       
  #####
  ###???Necessary???:
  param = CalcParamPloidyMechanisticSpecParamShort(param);
  if(sum(names(param)=='parssenslympho')>0)
  {
    if(param$parssenslympho==1)
    {
      param$cyclelymphb <- param$blymphrec
    }
  }
  if(sum(names(compartments)=='lymphmem_ind')>0)
  {
    param$complexlymphmem <- 1
  }else
  {
    param$complexlymphmem <- 0
  }
  
  ##output:
  param
}

  
CalcParamPloidyMechanisticSpecParamShort<-function(param)#function(param,specialpar)
{#param_ -<param!!!!02.08.2018!!!!
  param$errorlabel <- 0
  rel_4_8 <- param$rel_4_8 #specialpar$val[specialpar$names=="rel_4_8"]# param$rel_4_8
  rel_8_16 <- param$rel_8_16#specialpar$val[specialpar$names=="rel_8_16"]# param$rel_8_16
  rel_16_32<- param$rel_16_32#specialpar$val[specialpar$names=="rel_16_32"]# <- param$rel_16_32
  rel_32_64<- param$rel_32_64#specialpar$val[specialpar$names=="rel_32_64"]# <- param$rel_32_64
  rel_8_4 <- 1/rel_4_8
  rel_16_8 <- 1/rel_8_16
  rel_32_16 <- 1/rel_16_32
  rel_64_32 <- 1/rel_32_64
  rel_16_4 <- rel_16_8*rel_8_4
  rel_32_4 <- rel_16_4*rel_32_16
  rel_64_4 <- rel_32_4*rel_64_32
  
  param$rel_8_4 <- rel_8_4
  param$rel_16_8 <- rel_16_8
  param$rel_32_16 <- rel_32_16
  param$rel_64_32 <- rel_64_32
  param$rel_16_4 <- rel_16_4
  param$rel_32_4 <- rel_32_4
  param$rel_64_4 <- rel_64_4
  ###output:
  param
}
###############
###############
PKParamPharmODe<- function(param_,compartments)
{
  if(sum(names(compartments)=='paclitaxel_ind')>0)
  {
     #     param_$kpacl = param_$Clpacl/ param_$Vpacl1;
     #     param_$kpacl12 = param_$Qpacl2/param_$Vpacl1;
     #     param_$kpacl13 = param_$Qpacl3/param_$Vpacl1;
     #     param_$kpacl21 = param_$Qpacl2/param_$Vpacl2;
     #     param_$kpacl31 = param_$Qpacl3/param_$Vpacl3;
     if(sum(names(param_)=='paclitaxeltotal')>0)
     {
       if(param_$paclitaxeltotal==1)
       {
         #VMELpacl micromol/h
     
         param_$kpacl21 <- param_$kpacl210# is given. No V2!!!!!
         param_$kpacl13 <- param_$Qpacl3/param_$Vpacl1
         param_$kpacl31 <- param_$Qpacl3/param_$Vpacl3
       }else
       {
         param_$kpacl <- param_$Clpacl/param_$Vpacl1
         param_$kpacl12 <- param_$Qpacl2/param_$Vpacl1
         param_$kpacl13 <- param_$Qpacl3/param_$Vpacl1
         param_$kpacl21 <- param_$Qpacl2/param_$Vpacl2
         param_$kpacl31 <- param_$Qpacl3/param_$Vpacl3   
       }
    }else
    {
      param_$kpacl <- param_$Clpacl/param_$Vpacl1
      param_$kpacl12 <- param_$Qpacl2/param_$Vpacl1
      param_$kpacl13 <- param_$Qpacl3/param_$Vpacl1
      param_$kpacl21 <- param_$Qpacl2/param_$Vpacl2
      param_$kpacl31 <- param_$Qpacl3/param_$Vpacl3   
    }
  }
     if(sum(names(compartments)=='carboplatin_ind')>0)
     {
       param_$kcarbo <- param_$Clcarbo/param_$Vcarbo1
       param_$kcarbo12 <- param_$Qcarbo2/param_$Vcarbo1
       param_$kcarbo21 <- param_$Qcarbo2/param_$Vcarbo2  
     }
     ###output:
   param_
 }
###############
###############
thets <-function(Crel)
{  
  if (Crel<=1)
    {out_ <- 2/Crel^0.6}
  else
  {out_ <- 2}
  out_
}

in_out_from_amplif<-function(a)
{
  if(a==1)
  {
    in_ <- log(2)
  }
  else
  {
    if (a==0)
    {
      in_ <- 0
    }
    else
    {
      in_ <- (a-1)/log2(a)
    }
  }
  if (a!=0)
  {
    out_ <- a/in_}
  else
    out_ <- 0;
  end;
  c(in_, out_)
}


SetToSteadyState<- function (ind_indiv ,main_struct,param_)
{
  if (sum('PLC'==main_struct$YTYPE[[ind_indiv]]$list)==0&(sum('TPO'==main_struct$YTYPE[[ind_indiv]]$list)==0))
  {
    param_$x0 <- 1;
    param_$xtpo0 <- 1;
    param_$xtpo00 <- 1;
    param_$xcm0 <- 1;
    param_$xmkc0 <- 1;            
  }
  if(sum(names(main_struct)=='set_to_steady_state')>0){
    if(main_struct.set_to_steady_state==1)
    {
      param_$x0 <- 1
      param_$xtpo0 <- 1
      param_$xtpo00 <- 1
      param_$xcm0 <- 1
      param_$xmkc0 <- 1  
      param_$xtpo0 <- 1
      param_$xanc0 <- 1
      param_$xplc0 <- 1
      param_$xleuko0 <- 1
      param_$xloc0 <- 1
    }
  }
  ####output:
  param_
}

StratifyParam<-function(param_)
{
  output_<-list()
  paramlist<-list()
  paramnum<-list()
  paramnames<-list()
  for(ind_name in names(param_))
  {
    if(length(param_[[ind_name]])>1)
    {
      paramlist[[ind_name]]<-param_[[ind_name]]
    }else
    {
      paramnum[[ind_name]]<-param_[[ind_name]]
      paramnames[[ind_name]]<-names(param_)[ind_name]
    }
  }
  paramv<-vector(length=length(names(paramnum)))
  paramn<-names(paramnum)
  for(ind_name in 1:length(names(paramnum)))
  {
    #print(c(ind_name,names(param_)[ind_name]))
    paramv[ind_name]<-paramnum[[ind_name]]
    #paramn[ind_name]<-paramnames[[ind_name]]
  }
  output_$paramlist <- paramlist
  output_$paramnum <- paramnum
  output_$paramv <-paramv #output_$paramv <- c(paramv,0,10,0.5)#initialization of t00, tend and dt 
  output_$paramn <- paramn#output_$paramn <- c(paramn,'t00','tend','dt')
  output_$length_paramv <- length(paramv)#output_$length_paramv <- length(paramv)+3
  ####output:
  output_
}
EstimDistribIndivPar<-function(xtot,main_struct)#
{
  if(main_struct$numindiv>0)
  {
    if(main_struct$subjects_subgroup==0){#05.04.2016  only selected subjects contribute to the goal function!!!!
      seqrelev=1:main_struct$numindiv;#later redefine
    }else{
      seqrelev=main_struct$relevant_subjects;
    }
    main_struct$currav<-rep.int(x=0, times = main_struct$numindiv)#[1:main_struct$numindiv]=0;
    
    xtot$av<-rep.int(x=NA, times = main_struct$numindivpar1)#(1:main_struct$numindivpar1)=NA;    
    xtot$om<-rep.int(x=NA, times = main_struct$numindivpar1)#(1:main_struct$numindivpar1)=NA;
    for (ind_par in 1:main_struct$numindivpar1){
      temp_v=NULL;#[]
      temp_sum_ind=0;
      temp_sum=0;
      temp_sum_q=0;
      for (ind_temp in 1:length(seqrelev)){#main_struct$numindiv
        ind_relev=(main_struct$subs_indiv_num[ind_par]==main_struct$subs_indivi_num[[seqrelev[ind_temp]]]);
        if(sum(ind_relev)>0){
          temp_sum_ind=temp_sum_ind+1;
          temp_sum=temp_sum+xtot$indiv[[seqrelev[ind_temp]]][ind_relev==1];#main_struct$subs_indivi_num[[seqrelev[ind_temp]]](ind_relev);
          temp_sum_q=temp_sum_q+(xtot$indiv[[seqrelev[ind_temp]]][ind_relev==1])^2;
          temp_v=c(temp_v,xtot$indiv[[seqrelev[ind_temp]]][ind_relev==1]);
        }
      }
      if(temp_sum_ind>0){
        xtot$av[ind_par]=temp_sum/temp_sum_ind;
      }
      if(temp_sum_ind>1){
        xtot$om[ind_par]=(temp_sum_q/(temp_sum_ind-1)-(xtot$av[ind_par]^2)*temp_sum_ind/(temp_sum_ind-1))^0.5;
      }
      xtot$numindpar[ind_par]=temp_sum_ind;
    }
    
  }
  ####current output!!!! Must be updated
  xtot
}



UpdateField <-function(x,label_opt,ind,label_class,y,ind_indiv)
{
  if(label_class=='str')
  {
    #updates the ind -th element of either pop, indiv or resid field with x-value
    if(label_opt=='pop'){
      y$pop[ind]=x;
    }else if(label_opt=='indiv'){
      y$indiv[[ind_indiv]][ind]=x;
    }else if(label_opt=='resid')
    {
      y$resid[ind_indiv,ind]=x;
    }
  }else if(label_class=='vect')
  {
    #extracts the ind -th element of either pop, indiv or resid field. x is a structure.
    #clear y; !!!! Not necessary
    if(label_opt=='pop')
    {
      y=x$pop[ind];
    }
    else if(label_opt=='indiv')
    {
      y=x$indiv[[ind_indiv]][ind];
    }
    else if(label_opt=='resid')
    {
      y=x$resid[ind_indiv,ind]
    }
  }
  #output:
  y
}

#in the contrast to Matlab version, paths and file names were excluded from this function!!! optim option is always field of main_struct!!!!!
UpdateMainstuctforParametersEstimations<-function(main_struct)
{
  
  main_struct$optim_options$pop<-list()
  main_struct$optim_options$indiv<-list()
  main_struct$optim_options$resid<-list()
  main_struct$optim_options$indiv$LB<-list()
  main_struct$optim_options$indiv$UB<-list()
  main_struct$optim_options$indiv$eps_vect<-list()
  
  main_struct$optim_options$resid$LB<-list()
  main_struct$optim_options$resid$UB<-list()
  main_struct$optim_options$resid$eps_vect<-list()
  
  LB<-main_struct$optim_options$LB#Val[main_struct$optim_options$Name=="LB"]
  UB<-main_struct$optim_options$UB#Val[main_struct$optim_options$Name=="UB"]
  eps_init<-main_struct$optim_options$eps_init#Val[main_struct$optim_options$Name=="eps_init"]
  LB_resid<-main_struct$optim_options$LB_resid#Val[main_struct$optim_options$Name=="LB_resid"]
  UB_resid<-main_struct$optim_options$UB_resid#Val[main_struct$optim_options$Name=="UB_resid"]
  
  main_struct$optim_options$pop$LB<-rep.int(x=LB, times = main_struct$numpoppar)
  main_struct$optim_options$pop$UB<-rep.int(x=UB, times = main_struct$numpoppar)
  main_struct$optim_options$pop$eps_vect<-rep.int(x= eps_init, times = main_struct$numpoppar)
  for(ind_indiv in 1: main_struct$numindiv)
  {
    main_struct$optim_options$indiv$LB[[ind_indiv]] <-rep.int(x=LB, times = main_struct$num_estim_indiv[[ind_indiv]])
    main_struct$optim_options$indiv$UB[[ind_indiv]] <-rep.int(x=UB, times = main_struct$num_estim_indiv[[ind_indiv]])
    main_struct$optim_options$indiv$eps_vect[[ind_indiv]] <-rep.int(x=eps_init, times = main_struct$num_estim_indiv[[ind_indiv]])
    
    main_struct$optim_options$resid$LB[[ind_indiv]] <-rep.int(x=LB_resid, times = dim(main_struct$residinit)[2])
    main_struct$optim_options$resid$UB[[ind_indiv]] <-rep.int(x=UB_resid, times = dim(main_struct$residinit)[2])
    main_struct$optim_options$resid$eps_vect[[ind_indiv]] <- 0.5*main_struct$residinit[ind_indiv,]
  }
  #output:
  main_struct
}

new_dir<-"R:/Blutmodelle/Yuri/modelle/R files/IntegralModel1/DataRevisedNext/"

UpdateDefinitionsFromResults<-function(main_struct, relev_row,new_dir,old_dir)
{
  allfiles <-list.files(old_dir)#main_struct$resultspath
  for(ind_f in allfiles)
  {
    file.copy(paste(old_dir,ind_f,sep =""), paste(new_dir,ind_f,sep =""))
  }
  
  indiv_data<-read.csv(file=paste(new_dir,'indivpar.csv',sep=""), sep = ";",dec = main_struct$dec)
  resid_data<-read.csv(file=paste(new_dir,'resid.csv',sep=""), sep = ";",dec = main_struct$dec)
###############Read results until prescribed row:
  for(ind_indiv in 1: main_struct$numindiv)
  {
    tempori <- read.csv(file=paste(main_struct$resultspath,'maini',ind_indiv,'.csv',sep=""), sep = ";",dec = main_struct$dec)
    indiv_data[ind_indiv,2:dim(indiv_data)[2]] <- tempori[relev_row,]
    ###
    temporresidi <- read.csv(file=paste(main_struct$resultspath,'residi',ind_indiv,'.csv',sep=""), sep = ";",dec = main_struct$dec)
    resid_data[ind_indiv,2:dim(resid_data)[2]] <- temporresidi[relev_row,]
  }
  write.table(x = indiv_data, file = paste(new_dir,'indivpar.csv',sep=""), sep = ";" ,  dec = main_struct$dec, row.names=F)
  write.table(x = resid_data, file = paste(new_dir,'resid.csv',sep=""), sep = ";" ,  dec = main_struct$dec, row.names=F)
  main_data<-read.csv(file=paste(old_dir,'main.csv',sep=""), sep = ";")#!!!!!!new_dir
  #main_data$init_mean[2:(dim(indiv_data)[2]-1)] <- as.numeric(indiv_data[1,2:dim(indiv_data)[2]])#Any is good
  main_data$init_mean <- as.numeric(indiv_data[1,2:dim(indiv_data)[2]])#Any is good
  main_data$init_mean[main_data$names%in%c("xleuko0","xloc0","xanc0","xplc0","xtpo0")]<-1
  write.table(x =main_data, file = paste(new_dir,'main.csv',sep=""), sep = ";" ,  dec = main_struct$dec, row.names=F)
  ###output:
  new_dir
} 
Find_rep_array<-function(indiv_meas){
  #input_array$data is a two-dimentional data array, first row is a time, second 
  out_put<-array(dim=dim(input_array$data))
  all_outputs<-levels(indiv_meas$datatypes)
  for(ind_out in all_outputs)
  {
    relev_ind<- (indiv_meas$datatypes==ind_out)
    repeated_ind<-Find_rep_vect(indiv_meas$data[1,relev_ind])
  }
}
Find_rep_vect<-function(indiv_vect){
  l<-length(indiv_vect)
  out_<-rep.int(x=0,times = l)
  for(indi in 1:(l-1)){
    if(indiv_vect[indi]==indiv_vect[indi+1])
    {
      out_[indi+1]<-out_[indi]+1
    }
  }
  ###
  out_
}