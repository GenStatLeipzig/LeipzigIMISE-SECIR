MCMCHaararoundAvFaststep<-function(x,ind_indiv,xtot,main_struct,nonmem_options,num_var,ind_prev_res,nlocal,prev_cluster_str)#,compl_structindiv
{
  label_complete <- (-1)#0#1
  
  max_ellips_bound<-qchisq(p=0.975, df=num_var, ncp = 0, lower.tail = TRUE, log.p = FALSE)#0.975
  min_ellips_bound<-qchisq(p=0.025, df=num_var, ncp = 0, lower.tail = TRUE, log.p = FALSE)#0.025
  accept_min <- nonmem_options$accept_min
  we_min<- 0.0001#0.00000001
  stab_level <- nonmem_options$stab_level
  alpha_clust <- nonmem_options$alpha_clust#10#1/0.05#0.0125#40#exp(4)
  alpha_nll_clust <- nonmem_options$alpha_nll_clust
  min_steps_cov <- 2#4# 2#2#4
  step_stab <-5
  #No need in nonmem_options$metric_vect !!!!
  metric_array1<- vector(length = 4)
  max_av_cluster <- 100#100
  burnin_cluster <- nonmem_options$burnin_cluster#10#5010#50#10
  gam_mix <- nonmem_options$gam_mix
  min_cluster_length <- 10#5#10#50#10
  metric_array<- array(dim =  c(nlocal,13) )#-1 #13.12.2020!!!!
  #colnames(metric_array) <- c('j','PrevPiVal','CurrPiVal','RatioAccepted','mahalanobis_d','RatioCnd','DeltaRatio','DeltaNll','ind_cluster','num_steps','instab_val','instab_distr','alpha')
  colnames(metric_array)  <- c('j','PrevPiVal','CurrPiVal','RatioAccepted','mahalanobis_d','RatioCnd','DeltaRatio','DeltaNll','ind_cluster','num_steps','ind_prev_cluster','new_clust_source','alpha')
  cluster_str <- list()
  if(is.null(prev_cluster_str))
  {
    length_prev_cluster_str <- 0
    
  }else{
    length_prev_cluster_str <- dim(prev_cluster_str$cluster_values)[1]
  }
  max_cluster_num <- 5*(length_prev_cluster_str+ceiling( nlocal/(min_cluster_length)))#+burnin_cluster !!!!???
  cluster_str$cluster_values <- array(dim=c(max_cluster_num,8))
  colnames(cluster_str$cluster_values) <- c('cluster_val','nll_val','ind_clast','num_steps','we','av_val','hist','adjust_clust_val')#c('max_goal','max_nll','num','val')
  
  cluster_str$cluster_av <- array(dim=c(ceiling( length_prev_cluster_str+nlocal/(min_cluster_length)),num_var))#+burnin_cluster
  cluster_str$cluster_cov <- array(dim=c(num_var,num_var,ceiling( length_prev_cluster_str+nlocal/(min_cluster_length))))#+burnin_cluster
  cluster_str$adapt_acc_cov <- array(dim=c(num_var,num_var,ceiling(length_prev_cluster_str+ nlocal/(min_cluster_length))))#+burnin_cluster
  cluster_str$virt_acc_av <- array(dim=c(ceiling( length_prev_cluster_str+nlocal/(min_cluster_length)),num_var))#+burnin_cluster
  #ind_cluster <- length_prev_cluster_str
  num_clusters <- length_prev_cluster_str#ind_cluster
  #x <- x_start
  if(!is.null(prev_cluster_str))
  {
    start_clust_label <- 0
    cluster_str$cluster_values[1:length_prev_cluster_str,] <- prev_cluster_str$cluster_values
    cluster_str$cluster_av[1:length_prev_cluster_str,] <- prev_cluster_str$cluster_av[1:length_prev_cluster_str,]
    cluster_str$cluster_cov[,,1:length_prev_cluster_str] <- prev_cluster_str$cluster_cov[,,1:length_prev_cluster_str]
    cluster_str$adapt_acc_cov[,,1:length_prev_cluster_str] <- prev_cluster_str$adapt_acc_cov[,,1:length_prev_cluster_str]
    cluster_str$virt_acc_av[1:length_prev_cluster_str,] <- prev_cluster_str$virt_acc_av[1:length_prev_cluster_str,]
    if(num_clusters>1)
    {
      help_arr_cl<- vector(length = num_clusters)
      for(ind_c1 in 1:num_clusters)
      {
        help_arr_cl[ind_c1] <- MetricsMultivar(ve1=cluster_str$cluster_av[ind_c1,],
                                       ve2=x,
                                       mat=cluster_str$adapt_acc_cov[,,ind_c1])
      }
      ind_cluster <- (1:num_clusters)[help_arr_cl==min(help_arr_cl)][1]#AssigngClustNum(rel_we=cluster_str$cluster_values[1:num_clusters,'we'])
    }else{
      ind_cluster <- 1
    }
    
  }else{
    start_clust_label <- 1
    num_clusters <- 1
    ind_cluster <- 1
    cluster_str$virt_acc_av[ind_cluster,] <- nonmem_options$prevmean
    cluster_str$cluster_av[ind_cluster,] <- nonmem_options$prevmean
    cluster_str$cluster_cov[,,ind_cluster] <- nonmem_options$prevcov#acc_cov
    cluster_str$adapt_acc_cov[,,ind_cluster] <- cluster_str$cluster_cov[,,ind_cluster]+diag(nonmem_options$eps_chol,c(num_var,num_var))
    cluster_str$cluster_values[ind_cluster,'we'] <- 1
    
  }
  #proposed_var<- t(chol(cluster_str$adapt_acc_cov[,,ind_cluster]))
  if(nonmem_options$advance==0 )
  {  
    const_term<-num_var/2*log(2*pi)
    limi <- -log(0.05)/0.5
    nonmem_options$sd <- (nonmem_options$sd0/num_var)^0.5#Chol is equivalent to square root!!!
    sdsc<- (nonmem_options$sd0/num_var)
    nonmem_options$update_tann<-1
    nonmem_options$Tanninit<-nonmem_options$sd
    Tann <- rep.int(x=nonmem_options$Tanninit,times=nlocal+1)#1
    Tann[1] <- nonmem_options$Tanninit#1
  }else if(nonmem_options$advance==1 )
  {
    nonmem_options$sd <- 1#(nonmem_options$sd0/num_var)^0.5#Chol is equivalent to square root!!!
    sdsc<- 1#(nonmem_options$sd0/num_var)
    nonmem_options$update_tann<-1
    #nonmem_options$Tanninit<-nonmem_options$sd
    Tann <- rep.int(x=1,times=nlocal+1)#1
    Tann[10*(1:(nlocal/10))] <- 1.1
    Tann[20*(1:(nlocal/20))] <- 1.2
    Tann[50*(1:(nlocal/50))] <- 1.5
    Tann[100*(1:ceiling(nlocal/100))] <- 2
    Tann[200*(1:ceiling(nlocal/100))] <- 5
  }
  if(is.null(x))
  {
    proposed_var<- t(chol(cluster_str$adapt_acc_cov[,,ind_cluster]))
    x <- as.vector(ProposalDistInit(cluster_str$cluster_av[ind_cluster,],
                                    proposed_var*Tann[1]))
  }
  #print_ind <-  1000*(1:10)#1000*(1:(nonmem_options$m/1000)) 
  xtot<- UpdateParam(x=x,ind_indiv,xtot,stoch_obj0,main_struct)
  #num_var <- length(x)
  
  if(ind_prev_res==0)###????
  {
    start_label<-1
    #nonmem_options$prevres<-array(dim = c(ceiling(nonmem_options$m/nonmem_options$jump)+1,num_var))
    #ind_prev_res<-0
  }else{
    start_label<-0
    
  }
  #nonmem_options$adapt_acc_cov <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
  
  
  candappadv <- array(dim=c(nlocal,num_var))
  results <- vector(length = nlocal)
  accepted <- vector(length = nlocal)
  naccepted <- vector(length = nlocal)
  
  #print(length(Tann))
  #cand[1,]<-x
  candappadv[1,]<-x
  
  accepted[1]<- 1
  naccepted[1] <- 1
  #y<-func_pop_thrombo(x,xtot,main_struct)
  #ind_indiv <- 1
  ygoalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                       nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1= ind_indiv)
  y<- (- ygoalf$nLL_tot)
  y20<-  ygoalf$quadrresid
  if(main_struct$model_opt!=32)
  {  
    residapp <- array(dim=c(nlocal,length(ygoalf$quadrresid)))
    residapp_cumul <- array(dim=c(nlocal,length(ygoalf$quadrresid)))
    residapp[1,]<-y20
    residapp_cumul[1,] <- ygoalf$quadrresid_cumul
  }else{
    residapp <- list()
    residapp_cumul <- list()
    poly_s <- c( "0_14", "15_34", "35_59", "60_79", "80","all")#"15_59"
    for(ind_poly in poly_s)
    {
      residapp[[ind_poly]] <- array(dim=c(nlocal,length(ygoalf$quadrresid[[ind_poly]])))
      residapp_cumul[[ind_poly]] <- array(dim=c(nlocal,length(ygoalf$quadrresid[[ind_poly]])))
      residapp[[ind_poly]][1,]<-y20[[ind_poly]]
      residapp_cumul[[ind_poly]][1,] <- ygoalf$quadrresid_cumul[[ind_poly]]
    }
    
  }
  nonmem_options$nLL_indiv<-ygoalf$nLL_indiv
  
  j_accepted <- 1
  j <- 1
  results[1] <- y
  metric_array[j_accepted,1] <- 1
  metric_array[1,3] <- CalculatePiClever(x_v=x,cluster_str,num_clusters,ind_cluster,label_complete)
  metric_array[1,2] <- metric_array[1,3]#previous value
  metric_array[1,'mahalanobis_d'] <- MetricsMultivar(ve1=cluster_str$cluster_av[ind_cluster,],
                                                        ve2=x,
                                                        mat=cluster_str$adapt_acc_cov[,,ind_cluster])
  metric_array[1,'ind_cluster'] <- ind_cluster
  if(metric_array[1,3]==(-10^6))
  {
    cond_new_cluster <- TRUE
    #print(c(cluster_str$cluster_values[ind_cluster,1],alpha_unbound,alpha_unbound0))
    ind_prev_res <- min(burnin_cluster,ind_prev_res)
    
    #acc_av <- candappadv[j_accepted,]
    #acc_cov <- 0.2*acc_cov0+0.8*acc_cov
    ind_cluster_prev <- ind_cluster
    num_clusters <- num_clusters+1 #!!!!!
    ind_cluster <- num_clusters
    
    
    #metric_array[j_accepted,10] <- ind_prev_res-(burnin_cluster-1)
    #metric_array[j_accepted,11] <- (results[j_accepted])#Simple!!! without maximum 
    cluster_str$cluster_values[ind_cluster,1] <- metric_array[j_accepted,4]
    cluster_str$cluster_values[ind_cluster,2] <- results[j_accepted]#max(results[j_accepted],cluster_str$cluster_values[ind_cluster,1])
    cluster_str$cluster_values[ind_cluster,3] <- ind_cluster#j_accepted
    
    cluster_str$cluster_values[1:num_clusters,'we'] <- AssigngClustWeights(x = cluster_str$cluster_values[1:ind_cluster,2],gam_mix)
    if(cluster_str$cluster_values[num_clusters,'we']==0)
    {
      cluster_str$cluster_values[num_clusters,'we']<-0.0001
      cluster_str$cluster_values[1:num_clusters,'we'] <- cluster_str$cluster_values[1:num_clusters,'we']/sum(cluster_str$cluster_values[1:num_clusters,'we'])
    }
    cluster_str$cluster_values[ind_cluster,'hist']<- (-9)
    cluster_str$cluster_values[ind_cluster,'adjust_clust_val']<-cluster_str$cluster_values[ind_cluster,1] 
    #cluster_str$cluster_values[ind_cluster,1] <- max(metric_array[1:j_accepted,4])
    cluster_str$cluster_av[ind_cluster,] <- candappadv[j_accepted,]
    ##updates from previous cluster:
    cluster_str$virt_acc_av[ind_cluster,] <- cluster_str$virt_acc_av[ind_cluster_prev,]
    cluster_str$cluster_cov[,,ind_cluster] <- cluster_str$cluster_cov[,,ind_cluster_prev]
    cluster_str$cluster_values[ind_cluster,'av_val']<- metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,'adjust_clust_val']
    cluster_str$adapt_acc_cov[,,ind_cluster]  <- cluster_str$cluster_cov[,,ind_cluster]+diag(nonmem_options$eps_chol,c(num_var,num_var))
    #cluster_str$cluster_values[ind_cluster,'av_val']<- log(alpha_unbound0)#cluster_str$cluster_values[ind_cluster,1]#1
    metric_array[1,3] <- CalculatePiClever(x_v=x,cluster_str,num_clusters,ind_cluster,label_complete)
    metric_array[1,2] <- metric_array[1,3]#previous value
    metric_array[1,'mahalanobis_d'] <- MetricsMultivar(ve1=cluster_str$cluster_av[ind_cluster,],
                                                       ve2=x,
                                                       mat=cluster_str$adapt_acc_cov[,,ind_cluster])
    
  }
  
  #cluster_str$cluster_values[ind_cluster,4] <- 1
  if(start_clust_label==1)
  {  
    cluster_str$cluster_values[ind_cluster,4] <- 1###!!!!
    cluster_str$cluster_values[ind_cluster,2] <- results[j_accepted]
    cluster_str$cluster_values[ind_cluster,3] <- ind_cluster#j_accepted
    cluster_str$cluster_values[ind_cluster,1] <- metric_array[1,3]+results[1]
    cluster_str$cluster_values[ind_cluster,'adjust_clust_val']<- cluster_str$cluster_values[ind_cluster,1]
  }
  
  metric_array[1,4] <- cluster_str$cluster_values[ind_cluster,1]#results[1]+cluster_str$cluster_values[ind_cluster,1]
  
  ind_iter <-1 
  #xcenter <- candapp[1:j-1,]
  
  stop_cond <- F
  
  #num_clusters <-1
  while(!stop_cond) #for(j in 2:nlocal)
  {
    j<- j+1
    low_accept<- 0
    if((j_accepted==1)&&(j>accept_min))
    {
      low_accept <- 1
    }
 #   print(j)
   # if(j_accepted>1) ##12.01.2021 !!!!
   # {
      #print(metric_array[1:j_accepted,])
    
      if((j-metric_array[j_accepted,1])>accept_min)
      {
        low_accept <- 1
      }
   # print(c(j,low_accept,y,y1))  
    #}
    if(low_accept==1)##Exclude now---
    {
      print('Low aceptance!!!!!')
      proposed_var<- t(chol(cluster_str$adapt_acc_cov[,,ind_cluster]))
      acc_av <- cluster_str$cluster_av[ind_cluster,] 
      x3 <- as.vector(ProposalDistInit(cluster_str$cluster_av[ind_cluster,],
                                       proposed_var*Tann[j]))
      xtot3<- UpdateParam(x=x3,ind_indiv,xtot,stoch_obj0,main_struct)
      
      
      #xtot1<-InitUpdate(xtot1,main_struct)
      ygoalf3<-GoalFuncComp(xtot=xtot3,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                            nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
      y3<- - ygoalf3$nLL_tot
      metric_array3 <- CalculatePiClever(x_v=x3,cluster_str,num_clusters,ind_cluster,label_complete)
      z3 <- metric_array3+y3
      if(y3>y)
      {
        x1 <- x3
        y1 <- y3
        low_accept <- 2
      } else if((!is.na(z3))&&(z3>metric_array[j_accepted,4]-10000)){#100  10000
        xxx0<-x
        yyy0 <- y
        x<-x3
        y<-y3
      }
     
    }
    #####
    # nonmem_options abschaffen!!!!
    #label_opt1="pop"
    ############
   # if(num_clusters > (length_prev_cluster_str+0.5*(nlocal/(min_cluster_length))))#13.12.2020!!!!
   # {
   #   alpha_clust<-40
      #nonmem_options$eps_chol <- 2*(1e-08)
      #burnin_cluster <-50
   # }
   # if(num_clusters > (length_prev_cluster_str+0.75*(nlocal/(min_cluster_length))))#13.12.2020!!!!
   # {
   #   alpha_clust<-80
      #nonmem_options$eps_chol <- 2*(1e-08)
      #burnin_cluster <-50
  #  }
   # if(num_clusters > (length_prev_cluster_str+0.9*(nlocal/(min_cluster_length))))#13.12.2020!!!!
   # {
    #  alpha_clust<-1000
      #stop_cond <- TRUE
    #}
    if(num_clusters >= (length_prev_cluster_str+0.9*(nlocal/(min_cluster_length))))#13.12.2020!!!!
    {
      #alpha_clust<-1000
      stop_cond <- TRUE
    }
    if(num_clusters>1)
    {
      ind_cluster <- AssigngClustNum(rel_we=cluster_str$cluster_values[1:num_clusters,'we'])
    }else{
      ind_cluster <- 1
    }
    proposed_var<- t(chol(cluster_str$adapt_acc_cov[,,ind_cluster]))
    acc_av <- cluster_str$cluster_av[ind_cluster,] 
    if(low_accept < 2)
   {   
    x1 <- as.vector(ProposalDistInit(cluster_str$cluster_av[ind_cluster,],
                                    proposed_var*Tann[j]))
    
    xtot1<-xtot 
    xtot1<- UpdateParam(x=x1,ind_indiv,xtot,stoch_obj0,main_struct)
    #if(ind_indiv==0)
    #{  
   #   xtot1$phiopt<-generate_listvect(x1,xtot$phiopt,main_struct)
   # }
    #else{  
    #  xtot1$phiopt<-generate_listvCOVID(x1,ind_indiv,xtot$phiopt,main_struct)
    #}
    #if(main_struct$num_estim_indiv[ind_indiv]==0)
    #{  
   #   xtot1<-InitUpdate(xtot1,main_struct)
    #}else{
    #  xtot1<-InitUpdateMultiple(xtot1,ind_indiv,main_struct)
    #}
    #xtot1<-InitUpdate(xtot1,main_struct)
    
    ygoalf1<-GoalFuncComp(xtot=xtot1,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                          nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
    y1<- - ygoalf1$nLL_tot
    y2<- ygoalf1$quadrresid
    }
    #ptm <- proc.time()
     
    metric_array1[3] <- CalculatePiClever(x_v=x1,cluster_str,num_clusters,ind_cluster,label_complete)
    #if(label_complete==1)
    #{
      New_Pi <- CalculatePiClever(x_v=x,cluster_str,num_clusters,ind_cluster=metric_array[j_accepted,'ind_cluster']
                                    ,label_complete)#13.12.2020!!!
      #metric_array[j_accepted,2] <- metric_array[j_accepted,3]
      metric_array[j_accepted,3] <- New_Pi#CalculatePiClever(x_v=x,cluster_str,num_clusters)##13.12.2020!!! !!!update the older solution 09.12.2020
      metric_array[j_accepted,4] <- y+metric_array[j_accepted,3]
      
    #}
     metric_array[j_accepted,6] <- metric_array1[3] 
     instab_val <- metric_array[j_accepted,3]-metric_array[j_accepted,2]#abs(metric_array[j_accepted,3]-metric_array[j_accepted,2])#!!!!14.12.2020#abs(metric_array[j_accepted,3]-New_Pi)#13.12.2020!!!
    instab_distr <- FALSE#(instab_val>stab_level)||(metric_array1[3]==(-10^6))||(New_Pi==(-10^6))#5 2!!!! 4 10 #13.12.2020!!!
    #metric_array[j_accepted,11] <- instab_val#13.12.2020!!!!
    #metric_array[j_accepted,12] <- instab_distr#13.12.2020!!!!
    
    #if(instab_val>0)
    #{
      # print("instable!!!")
      # print(c(j,instab_val,instab_distr,NA,NA,NA,cluster_str$cluster_values[ind_cluster,4] ))
    #}
    #else if(){
    #  print("stable!!!")
    #  print(c(j,instab_val,instab_distr,NA,NA,NA,Inf,Inf,cluster_str$cluster_values[ind_cluster,4] ))
    #  }
    if(is.infinite(metric_array1[3])||is.infinite(metric_array[j_accepted,3]))
    {
      j_accepted
    }
    #metric_array[j_accepted,2] <-  metric_array[j_accepted,3]
    #metric_array1[1:3] <- MetricArraySimpUpd(ve_prev=acc_av,ve_curr=x1,#!!!x1, not x 08.12.2020
    #                                         cov_matr=nonmem_options$adapt_acc_cov )
    #metric_array0_1_3 <- MetricArraySimpUpd(ve_prev=acc_av,ve_curr=x,#!!!update the older solution 08.12.2020
    #                                         cov_matr=nonmem_options$adapt_acc_cov )
    metric_array1[4] <- y1+metric_array1[3]
    #metric_array[j_accepted,4] <- y+metric_array[j_accepted,3]  #!!!update the older solution 08,9.12.2020 We must compare Pi based on the same distribution!
    
    ###We found alternative to this:!!!!18.12.2020
    ###new_chain<-StartNewChaincond(metr_val=metric_array[j_accepted,4],num_clusters,cluster_str)
    new_chain<- 0
    if(start_clust_label==1)
    {
      
      cluster_str$cluster_values[ind_cluster,'av_val']<- metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,'adjust_clust_val'] #log(alpha_unbound0)#cluster_str$cluster_values[ind_cluster,1]#1
      cluster_str$cluster_values[ind_cluster,'hist']<- (6)
      alpha_unbound0 <- (exp(metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,'adjust_clust_val'] ))#(exp(metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,1] ))
      alpha_unbound_pi0 <- (exp(results[j_accepted]-cluster_str$cluster_values[ind_cluster,'av_val'] ))#(exp(metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,1] ))
      
      
    }
     #<- rep.int(x=0,times=4)
   # print(c(j_accepted,10*(proc.time() - ptm),ind_prev_res))
    metric_array[j_accepted,7] <- metric_array1[4]-metric_array[j_accepted,4]
    metric_array[j_accepted,8] <- y1-y
    if(is.na(metric_array1[4]))
    {
      alpha_<-0
      alpha_unbound<-0
    }else{
      alpha_ <- min(1,(exp(metric_array1[4]-metric_array[j_accepted,4])))#min(1,(exp(-y1+y))^(1/T_ann(j)))
      alpha_unbound <- (exp(metric_array1[4]-metric_array[j_accepted,4]))
    }
    #if(Tann[j]>1)
   # {
   #   print(c(Tann[j],alpha_,metric_array1[4]-metric_array[j_accepted,4],y1-y))
   # }
    if(is.na(y1))#((is.na(y1))||((low_accept==1)&&(alpha_==0)))#if((is.na(y1))||((low_accept==1)&&(alpha_==0)))
    {
      alpha_unbound_pi <-0
    }else{
      alpha_unbound_pi <-exp(y1-y)
    }
    if(low_accept>=1)
    {
      print('Low aceptance!!!!!')
    }
    #metric_array[j_accepted,5]<-new_chain
   # if((new_chain==1)||((y1-y)<(-15)))
    if((new_chain==1)&&(alpha_>=1))
    {
      alpha_unbound <-exp(y1-y)
      alpha_ <- min(1,alpha_unbound)
      
    } 
    if(instab_distr)#13.12.2020!!!!:
    {
      alpha_unbound <-exp(y1-y)
      alpha_ <- min(1,alpha_unbound)
    }
    if(is.na(alpha_))
    {
      alpha_
    }
    metric_array[j_accepted,'alpha'] <- alpha_#13.12.2020!!!!
    ##may be good:
    #print('y,y1,alpha:')
    #print(c(y,y1))
    #print(c(y1-y,alpha_))
    # if(is.na(alpha_))
    #{
    #  print(alpha_)
    #}
    label_accept<-0
    if(alpha_==1)
    {
      label_accept <-1
      
      label_accept <-1
    }else{
      add_guess <- runif(1)
      #generate uniformly distributed random number add_guess
      #("additional guess") and accept the candidate if it exeeds
      #alpha_; otherwise reject it and use the preceeding value.
      #Update the fields.
      if(add_guess<= alpha_)
      {
        
        label_accept <- 1
      }else{
       
        label_accept <- 0
      }
    }
    if(label_accept==1)
    {
      
      if((y1-y)<(-10))
      {  
        paradox <- 1
      }
      ind_prev_res <- ind_prev_res+1
      j_accepted <- j_accepted+1
      metric_array[j_accepted,1] <- j
      candappadv[j_accepted,]<-x1
      results[j_accepted]<-y1
      metric_array[j_accepted,3:4] <- metric_array1[3:4]  ###Not here!!! 14.12.2020
      metric_array[j_accepted,2] <-  metric_array[j_accepted-1,3]#first pi value  ###Not here!!! 14.12.2020
      if(is.na(metric_array[j_accepted,4] ))
      {
        print("Fuck")
      }
      accepted[j]<-1
      naccepted[j_accepted]<-1
      nonmem_options$nLL_indiv<- ygoalf1$nLL_indiv
      if(main_struct$model_opt!=32)
      {  
        residapp[j_accepted,]<-y2
        residapp_cumul[j_accepted,]<- ygoalf1$quadrresid_cumul
      }else{
        for(ind_poly in poly_s)
        {
          residapp[[ind_poly]][j_accepted,]<-y2[[ind_poly]]
          residapp_cumul[[ind_poly]][j_accepted,] <- ygoalf1$quadrresid_cumul[[ind_poly]]
        }
      }
      ###
      xtot<- xtot1
      x <- x1 ##19-12-2020
      y <- results[j_accepted]##19-12-2020
      
      metric_array[j_accepted,'mahalanobis_d'] <- MetricsMultivar(ve1=cluster_str$cluster_av[ind_cluster,],
                                                         ve2=x,
                                                         mat=cluster_str$adapt_acc_cov[,,ind_cluster])
      
      ###
    }
    else
    {
      #candapp[j,]<-candapp[j-1,]
      naccepted[j_accepted]<-naccepted[j_accepted]+1
      #results[j]<-results[j-1]
      accepted[j]<-0
      #residapp[j,]<-residapp[j-1,]
      #residapp_cumul[j,]<- residapp_cumul[j-1,]
    }
    #x<-candappadv[j_accepted,]#starting center for the next iteration ##19-12-2020
    #y <- results[j_accepted]##19-12-2020
    
    #xtot$phiopt<-generate_listvect(candappadv[j_accepted,],xtot$phiopt,main_struct) ##19-12-2020
    
    #print("results:")
    #if(j%in%print_ind)
    # {  
    #  print(c(j,j_accepted,results[j_accepted],naccepted[j_accepted],alpha_))#too long:,mean(accepted[1:j]))),Tann[j+1]
    # }
    
    if(label_accept==1)#(burnin_cluster-2)min_steps_cov
    {
      if((y1-y)<(-15))
      {
        print('!!!!')
        print(y)
      }

      #13.12.2020!!!!, much shorter:
      alpha_unbound0 <- (exp(metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,'adjust_clust_val'] ))#(exp(metric_array[j_accept
      alpha_unbound_pi0 <- (exp(results[j_accepted]-cluster_str$cluster_values[ind_cluster,'av_val'] ))#(exp(metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,1] ))
      #alpha_unbound00 <- ((metric_array[j_accepted,4]-cluster_str$cluster_values[ind_cluster,'av_val'] )/cluster_str$cluster_values[ind_cluster,'hist'])
      if(instab_distr)#13.12.2020!!!!:
      {
        alpha_unbound0 <- alpha_unbound#alpha_unbound_pi0
        cond_new_cluster <- (alpha_unbound0>alpha_nll_clust)
      }else{
        cond_new_cluster <- alpha_unbound0>alpha_clust#(alpha_unbound0>(exp(log(alpha_clust)+instab_val)))#&&(cluster_str$cluster_values[ind_cluster,4]>1)#!!!13.12.2020
      }
      cond_new_cluster <- ((metric_array[j_accepted,'mahalanobis_d']/sdsc)>max_ellips_bound)||((metric_array[j_accepted,'mahalanobis_d']/sdsc)<min_ellips_bound)
     #   x
      #Classical algorithm!!!!22.12.2020
      metric_array[j_accepted,"new_clust_source"]<-cond_new_cluster# 1 id from mahalanobis_d
      cond_better_clust <- ((results[j_accepted] - cluster_str$cluster_values[ind_cluster,2] )>=log(alpha_nll_clust))
      cond_new_cluster <- cond_new_cluster||((results[j_accepted] - cluster_str$cluster_values[ind_cluster,2] )>=log(alpha_nll_clust))#When nLL is much better then the cluster's one
      
      ####
      if(cond_new_cluster)
      {
        help_arr_cl<- vector(length = num_clusters)
        for(ind_c1 in 1:num_clusters)
        {
          help_arr_cl[ind_c1] <- MetricsMultivar(ve1=cluster_str$cluster_av[ind_c1,],
                                                 ve2=x,
                                                 mat=cluster_str$adapt_acc_cov[,,ind_c1])
        }
        ind_cluster0 <- (1:num_clusters)[help_arr_cl==min(help_arr_cl)][1]#AssigngClustNum(rel_we=cluster_str$cluster_values[1:num_clusters,'we'])
        help_mahalanobis_d <- MetricsMultivar(ve1=cluster_str$cluster_av[ind_cluster0,],
                                                                    ve2=x,
                                                                    mat=cluster_str$adapt_acc_cov[,,ind_cluster0])
        cond_new_cluster0 <- (help_mahalanobis_d>max_ellips_bound)||(help_mahalanobis_d<min_ellips_bound)
        cond_new_cluster0 <- cond_new_cluster||((results[j_accepted] - cluster_str$cluster_values[ind_cluster0,2] )>=log(alpha_nll_clust))#When nLL is much better then the cluster's one
        if(!cond_new_cluster0)
        {
          cond_new_cluster0
        }
      }
      ####
      
      if(is.na(cond_new_cluster)||(is.infinite(cond_new_cluster)))#13.12.2020!!!!:
      {
        print("Fuckkkkkkkkkkkkkkkk")
        cond_new_cluster <- FALSE
      }
      #cond_new_cluster <- (!instab_distr)&&(cond_new_cluster) 13.12.2020!!!!:???
      prev_cluster_length<- cluster_str$cluster_values[ind_cluster,4]
      #cond_new_cluster <- cond_new_cluster&&(cluster_str$cluster_values[ind_cluster,4] >= min_cluster_length)
      metric_array[j_accepted,'ind_prev_cluster']<- ind_cluster
      if(cond_new_cluster)#1.644854  1.959964
      {
        #print(c(cluster_str$cluster_values[ind_cluster,1],alpha_unbound,alpha_unbound0))
        ind_prev_res <- min(burnin_cluster,ind_prev_res)
       print(c(min_ellips_bound,metric_array[j_accepted,'mahalanobis_d'],max_ellips_bound,results[j_accepted] - cluster_str$cluster_values[ind_cluster,2]))
        #acc_av <- candappadv[j_accepted,]
        #acc_cov <- 0.2*acc_cov0+0.8*acc_cov
        ind_cluster_prev <- ind_cluster
        num_clusters <- num_clusters+1 #!!!!!
        ind_cluster <- num_clusters
        cluster_str$cluster_values[ind_cluster,'hist']<-cluster_str$cluster_values[ind_cluster_prev,'hist']
        if(metric_array[j_accepted,"new_clust_source"])
        {
          if(Tann[j]==1)
          {
            cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+1
          }else if(Tann[j]==1.2)
          {
            cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+2
          }else if(Tann[j]==1.5)
          {
            cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+3
          }else if(Tann[j]==2)
          {
            cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+4
          }
          else if(Tann[j]==5)
          {
            cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+5
          }
          if(cond_better_clust)
          {
            cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+7
          }
      }else
      {
        if(Tann[j]==1)
        {
          cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+9
        }else{
          cluster_str$cluster_values[ind_cluster,'hist'] <- 10*cluster_str$cluster_values[ind_cluster,'hist']+8
        }
      }
        #metric_array[j_accepted,10] <- ind_prev_res-(burnin_cluster-1)
        #metric_array[j_accepted,11] <- (results[j_accepted])#Simple!!! without maximum
        
        cluster_str$cluster_values[ind_cluster,1] <- metric_array[j_accepted,4]
        cluster_str$cluster_values[ind_cluster,2] <- results[j_accepted]#max(results[j_accepted],cluster_str$cluster_values[ind_cluster,1])
        cluster_str$cluster_values[ind_cluster,3] <- ind_cluster#j_accepted
        
        cluster_str$cluster_values[ind_cluster,4] <- 1###!!!!
        if(is.na(cluster_str$cluster_values[ind_cluster,1]))
        {
          metric_array[j_accepted,4]
        }
        cluster_str$cluster_values[1:num_clusters,'we'] <- AssigngClustWeights(x = cluster_str$cluster_values[1:ind_cluster,2],gam_mix)
        cluster_str$cluster_values[ind_cluster,'av_val']<- log(alpha_unbound0)#cluster_str$cluster_values[ind_cluster,1]#1
        #cluster_str$cluster_values[ind_cluster,'hist']<- 0
        cluster_str$cluster_values[ind_cluster,'adjust_clust_val']<-cluster_str$cluster_values[ind_cluster,1] 
        #cluster_str$cluster_values[ind_cluster,1] <- max(metric_array[1:j_accepted,4])
        cluster_str$cluster_av[ind_cluster,] <- candappadv[j_accepted,]
        ##updates from previous cluster:
        cluster_str$virt_acc_av[ind_cluster,] <- cluster_str$virt_acc_av[ind_cluster_prev,]
        cluster_str$cluster_cov[,,ind_cluster] <- cluster_str$cluster_cov[,,ind_cluster_prev]
      }else
      {
        
        ###
        
        
        #cluster_str$adapt_acc_cov[,,ind_cluster] <- cluster_str$cluster_cov[,,ind_cluster]+diag(nonmem_options$eps_chol,c(num_var,num_var))
        #cluster_str$cluster_av[ind_cluster,] <- candappadv[j_accepted,]
        #proposed_var<- t(chol(nonmem_options$adapt_acc_cov))#acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
        ###
        
        #metric_array[j_accepted,11] <- (results[j_accepted])#Simple!!! without maximum 
        #cluster_str$cluster_values[ind_cluster,1] <- max(cluster_str$cluster_values[ind_cluster,1],metric_array[j_accepted,4])
        #cluster_str$cluster_values[ind_cluster,2] <- max(results[j_accepted],cluster_str$cluster_values[ind_cluster,1])
        #cluster_str$cluster_values[ind_cluster,3] <- ind_prev_res-(min(ind_prev_res,burnin_cluster)-1)
        cluster_str$cluster_values[ind_cluster,4] <- cluster_str$cluster_values[ind_cluster,4]+1###!!!!
      
         # help_cluster_cov_av <-  RecursiveCovMatrixAvUpd(xn= log(alpha_unbound0),#metric_array[j_accepted,4],# results[j_accepted] metric_array[j_accepted,1]#nonmem_options$prevres[ind_prev_res+1,],
         #                                      mup = cluster_str$cluster_values[ind_cluster,'av_val'],
         #                                      covp = cluster_str$cluster_values[ind_cluster,'hist'],
         #                                      nn = cluster_str$cluster_values[ind_cluster,4])#ind_prev_res+1!!!???
        
        #cluster_str$cluster_values[ind_cluster,'av_val'] <- help_cluster_cov_av$mu
        #cluster_str$cluster_values[ind_cluster,'hist'] <- help_cluster_cov_av$cov^0.5
        #cluster_str$cluster_cov[,,ind_cluster] <- nonmem_options$adapt_acc_cov#This is always updated!
        metric_arrayhelp1_3 <- CalculatePiClever(x_v=cluster_str$cluster_av[ind_cluster,],cluster_str,num_clusters,ind_cluster,label_complete)
        cluster_str$cluster_values[ind_cluster,'adjust_clust_val']<- cluster_str$cluster_values[ind_cluster,2]+metric_arrayhelp1_3
      }
      metric_array[j_accepted,10]<- cluster_str$cluster_values[ind_cluster,4]
      ###
      ##update current clusters covariance matrix and corresponding proposal distribution's element:
      if((low_accept==0)||(cond_new_cluster))
     {   
        help_cov_av<-  RecursiveCovMatrixAvUpd(xn= candappadv[j_accepted,],#nonmem_options$prevres[ind_prev_res+1,],
                                               mup = cluster_str$virt_acc_av[ind_cluster,],
                                               covp = cluster_str$cluster_cov[,,ind_cluster],
                                               nn = max(burnin_cluster,prev_cluster_length))#max(burnin_cluster,cluster_str$cluster_values[ind_cluster,4]))
        max(burnin_cluster,prev_cluster_length)
      
      #ind_prev_res+1!!!???  num_steps
        cluster_str$virt_acc_av[ind_cluster,] <- help_cov_av$mu#  changes by changing cluster!!!
        cluster_str$cluster_cov[,,ind_cluster] <-  help_cov_av$cov_
        
        cluster_str$adapt_acc_cov[,,ind_cluster] <- cluster_str$cluster_cov[,,ind_cluster]+diag(nonmem_options$eps_chol,c(num_var,num_var))
        #cluster_str$cluster_av[ind_cluster,] <- candappadv[j_accepted,] No!!!! 20.12.2020
        proposed_var<- t(chol(cluster_str$adapt_acc_cov[,,ind_cluster]))#acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
        #cluster_str$cluster_cov[,,ind_cluster] <- cluster_str$adapt_acc_cov[,,ind_cluster]#acc_cov
      ###
      }else{
        print("bad values")
      }
      metric_array[j_accepted,'ind_cluster']<- ind_cluster
      if(main_struct$label_branch >1)
      {
        print(c(j,j_accepted,cluster_str$cluster_values[ind_cluster,'num_steps'],ind_cluster,metric_array[j_accepted,4], results[j_accepted],alpha_unbound,alpha_unbound0,cluster_str$cluster_values[ind_cluster,'adjust_clust_val']))# alpha_unbound00 ,alpha_unbound0))#c(4:8) metric_array[j_accepted,c(4:6)]
        print(c(j_accepted,metric_array[j_accepted,4]-metric_array[j_accepted-1,4],results[j_accepted]-results[j_accepted-1],new_chain))
      }
        
      j_accepted
      
       
    }
    #}
    
    #-must be updated at the end!!!
    #nonmem_options$adapt_acc_cov <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
    #nonmem_options$proposed_var<- t(chol(nonmem_options$adapt_acc_cov))#acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
    
    stop_cond <- stop_cond||(j>=nlocal)#((j>=nlocal)&&metric_array[j_accepted,8])
  }
  
  ###output:
  out_ <- list()
  # out_$cand <- cand
  out_$metric_array <- metric_array[1:j_accepted,]
  out_$results <-results[1:j_accepted]
  out_$total_accepted <- j_accepted
  out_$accepted <- accepted
  out_$naccepted <- naccepted[1:j_accepted]
  out_$Tann <- Tann
  out_$candapp <- candappadv[1:j_accepted,]
  if(main_struct$model_opt!=32)
  {  
    out_$residapp <- residapp[1:j_accepted,]
    out_$residapp_cumul <- residapp_cumul[1:j_accepted,]
  }else{
    for(ind_poly in poly_s)
    {
      out_$residapp[[ind_poly]] <- residapp[[ind_poly]][1:j_accepted,]
      out_$residapp_cumul[[ind_poly]] <- residapp_cumul[[ind_poly]][1:j_accepted,]
    }
  }
  #out_$proposed_var <- nonmem_options$proposed_var
  #out_$prevmean <- virt_acc_av#acc_av#help_cov_av$mu
  #out_$prevcov <- acc_cov#help_cov_av$cov_
  iclusters <- 1:num_clusters
  iclusters_r <- iclusters[cluster_str$cluster_values[1:num_clusters,'we']>we_min]
   
  out_$num_clusters <- length(iclusters_r)#num_clusters
  if(out_$num_clusters>1)#
  {  
  out_$cluster_str$cluster_values <- cluster_str$cluster_values[1:num_clusters,][iclusters_r,]
  out_$cluster_str$cluster_av <- cluster_str$cluster_av[1:num_clusters,][iclusters_r,]
  out_$cluster_str$cluster_cov <- cluster_str$cluster_cov[,,1:num_clusters][,,iclusters_r]
  out_$cluster_str$adapt_acc_cov <- cluster_str$adapt_acc_cov[,,1:num_clusters][,,iclusters_r]
  out_$cluster_str$virt_acc_av <- cluster_str$virt_acc_av[1:num_clusters,][iclusters_r,]
  }else{
    out_$cluster_str$cluster_values <- array( data = cluster_str$cluster_values[iclusters_r,],dim =c(1,dim(cluster_str$cluster_values)[2]),
                                              dimnames = dimnames(cluster_str$cluster_values))
    out_$cluster_str$cluster_av <- array( data = cluster_str$cluster_av[iclusters_r,],dim =c(1,dim(cluster_str$cluster_av)[2]))
    out_$cluster_str$cluster_cov <- array( data = cluster_str$cluster_cov[,,iclusters_r],dim =c(dim(cluster_str$cluster_cov)[1:2],1))
    out_$cluster_str$adapt_acc_cov <- array( data = cluster_str$adapt_acc_cov[,,iclusters_r],dim =c(dim(cluster_str$adapt_acc_cov)[1:2],1))
    out_$cluster_str$virt_acc_av <- array( data = cluster_str$virt_acc_av[iclusters_r,],dim =c(1,dim(cluster_str$virt_acc_av)[2]))
  }
  out_$cluster_str$cluster_values[,'we']<-out_$cluster_str$cluster_values[,'we']/sum(out_$cluster_str$cluster_values[,'we'])
  out_$ind_prev_res <- ind_prev_res
  out_$alpha_unbound0 <- alpha_unbound0
  out_$lastclustval <- cluster_str$cluster_values[num_clusters,1]
  #out_$ <-
  ###
  out_
  
} 

MCMCHaararoundAvGlobal<-function(x,ind_indiv,xtot,main_struct,nonmem_options,num_var,nlocal,nglobal,ind_prev_res,prev_cluster_str)#,compl_structindiv
{
  #nonmem_options$prevcov <- precov
  ntot<- nlocal*nglobal
  print('global steps =')
  print(nglobal)
  print('local steps =')
  print(nlocal)
  ptm <- proc.time()
  if(is.null(prev_cluster_str))
  { 
  localMCMCi <-MCMCHaararoundAvFaststep(x=x,ind_indiv,xtot,main_struct,nonmem_options,num_var=num_var,
                                         ind_prev_res=ind_prev_res,nlocal=nlocal,prev_cluster_str=NULL)#,compl_structindiv
  }else{
    localMCMCi <-MCMCHaararoundAvFaststep(x=x,ind_indiv,xtot,main_struct,nonmem_options,num_var=num_var,
                                          ind_prev_res=ind_prev_res,nlocal=nlocal,prev_cluster_str=prev_cluster_str)#,compl_structindiv
  }
  print((proc.time() - ptm))
 # write.table(x=localMCMCi$candapp, sep = ";" ,  dec = main_struct$dec, 
 #             row.names=F,col.names = F,file=paste(main_struct$resultspath,
 #                                                  "MCMCStore.csv",sep=''), append = F)
 # write.table(x=cbind(localMCMCi$naccepted,localMCMCi$results),sep = ";" ,  dec = main_struct$dec, 
#              row.names=F,col.names = F,file=paste(main_struct$resultspath,
#                                                   "ResultsStore.csv",sep=''), append = F)
  
 # nonmem_options$prevcov<-localMCMCi$prevcov
 # nonmem_options$prevmean <- localMCMCi$prevmean
 # nonmem_options$proposed_var <- localMCMCi$proposed_var
  nonmem_options$metric_vect <- localMCMCi$metric_array[dim(localMCMCi$metric_array)[1],]
  nonmem_options$alpha_unbound0 <- localMCMCi$alpha_unbound0 
  nonmem_options$lastclustval <- localMCMCi$lastclustval 
  #diagi <- diag(localMCMCi$prevcov)
  #print(c(min(diagi),mean(diagi),max(diagi))^0.5)
  #diagi <- diag(localMCMCi$proposed_var)
  #print(c(min(diagi),mean(diagi),max(diagi)))#^0.5 
  print(c(min(apply(localMCMCi$candapp,2,sd)),(mean(apply(localMCMCi$candapp,2,var))^0.5),max(apply(localMCMCi$candapp,2,sd))))
  glob_results <- localMCMCi$results
  glob_naccepted  <- localMCMCi$naccepted
  glob_candapp  <- localMCMCi$candapp
  glob_residapp  <- localMCMCi$residapp
  glob_residapp_cumul  <- localMCMCi$residapp_cumul
  glob_metric_array <- localMCMCi$metric_array 
  glob_cluster_str <- localMCMCi$cluster_str 
  cluster_cov0<- localMCMCi$cluster_str$cluster_cov
  ind_prev_resi<-0  
  #numi_local_steps<-vector(length=ntot)
  #numi_local_steps[1]<-length(localMCMCi$results)
  ind_local <- 1
  for(ind_local in 1:nglobal)#251while(ind_local <=nglobal)#for(ind_local in 1:nglobal)#251
  {
   # if((nlocal/10)<length(localMCMCi$results))
   # {
      
   #   nonmem_options$burnin_cluster <-  100
   #   nonmem_options$gam_mix <- 0.75
   #   nonmem_options$alpha_clust <- 1/0.025#10
   # }
    #if((nlocal/20)>length(localMCMCi$results))
    #{
   #   nonmem_options$burnin_cluster <- 10
   #   nonmem_options$gam_mix <- 0.5
   #   nonmem_options$alpha_clust <- 20#10
   # }
    localMCMCi0<-localMCMCi
    #ind_local <- ind_local+1
    ind_prev_resi <- ind_prev_resi+length(localMCMCi$results)
    ind_prev_res <- ind_prev_res+length(localMCMCi$results)
    if(length(localMCMCi$results)>1)
    {
      localMCMCi <- localMCMCi0
     # nonmem_options$sd0 <- num_var#length(xtot$phiopt$indiv[[ind_indiv]])#13.12.2020!!!!
    }
    
      ptm <- proc.time()
      localMCMCi <- MCMCHaararoundAvFaststep(x=localMCMCi$candapp[length(localMCMCi$results),],#x=NULL,
                                             ind_indiv,xtot,main_struct,nonmem_options,num_var=num_var,
                                             ind_prev_res=localMCMCi$ind_prev_res,nlocal=nlocal,prev_cluster_str=glob_cluster_str)#,compl_structindiv
      
  
      
      print(c(ind_local,(proc.time() - ptm),length(localMCMCi$results),ind_prev_resi))
    if(length(localMCMCi$results)>1)
    {
      if(localMCMCi$total_accepted>1)
      {  
        glob_results <- c(glob_results,localMCMCi$results)
        glob_naccepted  <- c(glob_naccepted,localMCMCi$naccepted)
        glob_candapp  <- rbind(glob_candapp,localMCMCi$candapp)
        if(main_struct$model_opt!=32)
        {  
          glob_residapp  <- rbind(glob_residapp, localMCMCi$residapp)
          glob_residapp_cumul  <- rbind(glob_residapp_cumul,localMCMCi$residapp_cumul)
        }else{ 
          for(ind_poly in poly_s)
          {
            glob_residapp[[ind_poly]]  <- rbind(glob_residapp[[ind_poly]], localMCMCi$residapp[[ind_poly]])
            glob_residapp_cumul[[ind_poly]]  <- rbind(glob_residapp_cumul[[ind_poly]],localMCMCi$residapp_cumul[[ind_poly]])
          }
        }
        glob_metric_array <- rbind(glob_metric_array,localMCMCi$metric_array) 
        glob_cluster_str <- localMCMCi$cluster_str
        
        #glob_cluster_str$cluster_values <- rbind(glob_cluster_str$cluster_values,localMCMCi$cluster_str$cluster_values)
        #glob_cluster_str$cluster_av <- rbind(glob_cluster_str$cluster_av,localMCMCi$cluster_str$cluster_av)
        #glob_cluster_str$cluster_cov<-array(dim=dim(localMCMCi$cluster_str$cluster_cov)+c(0,0,dim(cluster_cov0)[3]))
        #for(indi in 1:dim(cluster_cov0)[3])
        #{
        #  glob_cluster_str$cluster_cov[,,indi] <- cluster_cov0[,,indi]
        #}
       # for(indi in 1:dim(localMCMCi$cluster_str$cluster_cov)[3])
       # {
       #   glob_cluster_str$cluster_cov[,,dim(cluster_cov0)[3]+indi] <- localMCMCi$cluster_str$cluster_cov[,,indi]
       # }
        cluster_cov0<- glob_cluster_str$cluster_cov
        
        
        nonmem_options$proposed_var <- localMCMCi$proposed_var
        #numi_local_steps[ind_local+1]<-length(localMCMCi$results)
       # nonmem_options$prevcov<-localMCMCi$prevcov
       # nonmem_options$prevmean <- localMCMCi$prevmean
       # nonmem_options$metric_vect <- localMCMCi$metric_array[dim(localMCMCi$metric_array)[1],]
        nonmem_options$alpha_unbound0 <- localMCMCi$alpha_unbound0 
        nonmem_options$lastclustval <- localMCMCi$lastclustval 
     #   write.table(x=localMCMCi$candapp, sep = ";" ,  dec = main_struct$dec, 
     #               row.names=F,col.names = F,file=paste(main_struct$resultspath,
     #                                                    "MCMCStore.csv",sep=''), append = T)
      #  write.table(x=cbind(localMCMCi$naccepted,localMCMCi$results),sep = ";" ,  dec = main_struct$dec, 
      #              row.names=F,col.names = F,file=paste(main_struct$resultspath,
      #                                                   "ResultsStore.csv",sep=''), append = T)
      }
      #diagi <- diag(localMCMCi$prevcov)
      #print(c(min(diagi),mean(diagi),max(diagi))^0.5)
     # diagi <- diag(localMCMCi$proposed_var)
     # print(c(min(diagi),mean(diagi),max(diagi)))#^0.5
      print(c(min(apply(localMCMCi$candapp,2,sd)),(mean(apply(localMCMCi$candapp,2,var))^0.5),max(apply(localMCMCi$candapp,2,sd))))
      print(c(min(apply(glob_candapp,2,sd)),(mean(apply(glob_candapp,2,var))^0.5),max(apply(glob_candapp,2,sd))))
      glob_candapp
    }else{
      print('no accepted values!')
    }
  }
  output_ <- list()
  output_$results <- glob_results
  output_$candapp <- glob_candapp
  output_$naccepted <- glob_naccepted
  output_$residapp <- glob_residapp
  output_$residapp_cumul <- glob_residapp_cumul
  output_$prevcov <- nonmem_options$prevcov
  output_$prevmean <- nonmem_options$prevmean
  output_$proposed_var <- nonmem_options$proposed_var
  output_$metric_array <- glob_metric_array
  output_$metric_array <- glob_metric_array
  output_$cluster_str <- glob_cluster_str#$cluster_values
  output_$ind_prev_res <- localMCMCi$ind_prev_res
  output_$alpha_unbound0 <- localMCMCi$alpha_unbound0 
  output_$lastclustval <- localMCMCi$lastclustval 
  return(output_)
}



MetricArrayUpd<-function(metric_vect_prev,ind_prev_res,ve_prev,ve_curr,cov_matr,results)
{
  metric_vect <- vector(length=8)
  
  metric_vect[1]<- MetricsMultivar(ve1=ve_curr,
                                   ve2=ve_prev,
                                   mat=cov_matr)
  metric_vect[2]<- determinant(cov_matr)$modulus
  metric_vect[3]<- 0.5*(metric_vect[1]+metric_vect[2])
  metric_vect[4] <- results[j_accepted]+metric_vect[3]
  help_metric_cov_av<-  RecursiveCovMatrixAvUpd(xn= metric_vect[4],#nonmem_options$prevres[ind_prev_res+1,],
                                                mup = metric_vect_prev[5],
                                                covp = metric_vect_prev[6],
                                                nn = ind_prev_res)#ind_prev_res+1!!!???
  metric_vect[5] <- help_metric_cov_av$mu
  metric_vect[6] <-  help_metric_cov_av$cov_^0.5
  metric_vect[7]<- metric_vect[4]-metric_vect[5]
  metric_vect[8]<-  metric_vect[7]/metric_vect[6]
  
  return(metric_vect)
}
MetricArrayInitiate<-function(metric_vect_prev,ind_prev_res,matr_prev,ve_curr,cov_matr,results)
{
  j_accepted <- 4#!!!
  metric_array <- array(dim = c(j_accepted,8))
  
  metric_array[1:j_accepted,2]<- determinant(cov_matr)$modulus
  for(ind_j in 1:  j_accepted)
  {
    metric_array[ind_j,1]<- MetricsMultivar(ve1=matr_prev[ind_j,],ve2=ve_curr,
                                            mat=nonmem_options$adapt_acc_cov)
    
    metric_array[ind_j,3]<- 0.5*(metric_array[ind_j,1]+metric_array[ind_j,2])
    metric_array[ind_j,4] <- results[ind_j]+metric_array[ind_j,3]
  }
  metric_array[1:j_accepted,5]<- mean(metric_array[1:j_accepted,4])
  metric_array[1:j_accepted,7]<- metric_array[1:j_accepted,4]-metric_array[1:j_accepted,5]
  metric_array[1:j_accepted,6]<- sd(metric_array[1:j_accepted,4])
  metric_array[1:j_accepted,8]<-  metric_array[1:j_accepted,7]/metric_array[1:j_accepted,6]
  
  
  return(metric_array)
}

OriginalMCMCHaararoundAvFaststep<-function(x,ind_indiv,xtot,main_struct,nonmem_options,num_var,ind_prev_res,nlocal,prev_cluster_str)#,compl_structindiv
{
  alpha_clust <- 1/0.025#0.0125#40#exp(4)
  min_steps_cov <- 5#4
  #No need in nonmem_options$metric_vect !!!!
  metric_array1<- vector(length = 4)
  max_av_cluster <- 100#100
  burnin_cluster <-5#5
  min_cluster_length <- 0#10#10
  cluster_str <- list()
  max_cluster_num <- length(prev_cluster_str)+ceiling( nlocal/(min_cluster_length+burnin_cluster))
  cluster_str$cluster_values <- array(dim=c(max_cluster_num,4))
  colnames(cluster_str$cluster_values)<-c('max_goal','max_nll','num','val')
  cluster_str$cluster_av <- array(dim=c(ceiling( nlocal/(min_cluster_length+burnin_cluster)),num_var))
  cluster_str$cluster_cov <- array(dim=c(num_var,num_var,ceiling( nlocal/(min_cluster_length+burnin_cluster))))
  const_term<-num_var/2*log(2*pi)
  limi <- -log(0.05)/0.5
  metric_array<- array(dim =  c(nlocal,10) )#-1
  #print_ind <-  1000*(1:10)#1000*(1:(nonmem_options$m/1000)) 
  xtot<- UpdateParam(x=x,ind_indiv,xtot,stoch_obj0,main_struct)
  #num_var <- length(x)
  if(ind_prev_res==0)
  {
    if('prevcov'%in%names(nonmem_options))
    {
      acc_cov <- nonmem_options$prevcov
    }else{
      acc_cov <- array(data = 0, dim=c(num_var,num_var))
    }
    acc_av <- x
    
    
    #nonmem_options$prevres<-array(dim = c(ceiling(nonmem_options$m/nonmem_options$jump)+1,num_var))
    #ind_prev_res<-0
  }else{
    
    #ind_prev_res<- dim(nonmem_options$prevres)[1]
    #nonmem_options$prevres<- rbind(nonmem_options$prevres,array(dim = c(ceiling(nonmem_options$m/nonmem_options$jump)+1,num_var)))
    ##no prev_cov!
    acc_av <- nonmem_options$prevmean
    acc_cov <- nonmem_options$prevcov
  }
  nonmem_options$sd <- (nonmem_options$sd0/num_var)^0.5#Chol is equivalent to square root!!!
  nonmem_options$update_tann<-1
  nonmem_options$Tanninit<-nonmem_options$sd
  
  ##num_chains is currently 1
  #xx<-generate_pop_vector(xtot,main_struct1)
  #cand <- array(dim=c(nlocal,num_var))
  candappadv <- array(dim=c(nlocal,num_var))
  results <- vector(length = nlocal)
  accepted <- vector(length = nlocal)
  naccepted <- vector(length = nlocal)
  Tann <- rep.int(x=nonmem_options$Tanninit,times=nlocal+1)#1
  #print(length(Tann))
  #cand[1,]<-x
  candappadv[1,]<-x
  Tann[1] <- nonmem_options$Tanninit#1
  accepted[1]<- 1
  naccepted[1] <- 1
  #y<-func_pop_thrombo(x,xtot,main_struct)
  #ind_indiv <- 1
  ygoalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                       nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1= ind_indiv)
  y<- (- ygoalf$nLL_tot)
  y20<-  ygoalf$quadrresid
  residapp <- array(dim=c(nlocal,length(ygoalf$quadrresid)))
  residapp_cumul <- array(dim=c(nlocal,length(ygoalf$quadrresid)))
  residapp[1,]<-y20
  residapp_cumul[1,] <- ygoalf$quadrresid_cumul
  nonmem_options$nLL_indiv<-ygoalf$nLL_indiv
  results[1] <- y
  
  metric_array[1,1:3] <- MetricArraySimpUpd(ve_prev=acc_av,ve_curr=x,
                                            cov_matr=nonmem_options$adapt_acc_cov <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var)))
  metric_array[1,4] <- results[1]+metric_array[1,3]
  ind_iter <-1 
  #xcenter <- candapp[1:j-1,]
  j_accepted <- 1
  stop_cond <- F
  j <- 1
  num_clusters <-1
  while(!stop_cond) #for(j in 2:nlocal)
  {
    j<- j+1
    
    #####
    # if(nonmem_options$type[j]==0)#population
    # {
    #    label_opt1="pop"
    #    x0<-x
    #  }else{
    #   label_opt1="indiv"
    #    ind_indiv<-nonmem_options$type[j]
    #    x0<- generate_indiv_vector (xtot,main_struct,ind_indiv)
    #  }
    #label_opt1="pop"
    ############
    
    x1<-ProposalDistInit(acc_av,nonmem_options$proposed_var*Tann[j])#proposal dist acc_av
    #cand[j,]<-x1
    #y1<-func_pop_thrombo(x1,xtot,main_struct)
    xtot1<-xtot 
    if(ind_indiv==0)
    {  
      xtot1$phiopt<-generate_listvect(x1,xtot$phiopt,main_struct)
    }
    else{  
      xtot1$phiopt<-generate_listvCOVID(x1,ind_indiv,xtot$phiopt,main_struct)
    }
    if(main_struct$num_estim_indiv[ind_indiv]==0)
    {  
      xtot1<-InitUpdate(xtot1,main_struct)
    }else{
      xtot1<-InitUpdateMultiple(xtot1,ind_indiv,main_struct)
    }
    #xtot1<-InitUpdate(xtot1,main_struct)
    ygoalf1<-GoalFuncComp(xtot=xtot1,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                          nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
    y1<- - ygoalf1$nLL_tot
    y2<- ygoalf1$quadrresid
    
    #ptm <- proc.time()
    metric_array1[1:3] <- MetricArraySimpUpd(ve_prev=acc_av,ve_curr=x,
                                             cov_matr=nonmem_options$adapt_acc_cov <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var)))
    metric_array1[4] <- y1+metric_array1[3]
    #<- rep.int(x=0,times=4)
    # print(c(j_accepted,10*(proc.time() - ptm),ind_prev_res))
    if(ind_prev_res<5)
    {
      if(is.na(y1))
      {
        alpha_<-0
        alpha_unbound<-0
      }else{
        alpha_ <- min(1,(exp(y1-y)))#min(1,(exp(-y1+y))^(1/T_ann(j)))
        alpha_unbound <- exp(y1-y)
      }
    }else
    {
      if(is.na(metric_array1[4]))
      {
        alpha_<-0
        alpha_unbound<-0
      }else{
        alpha_ <- min(1,(exp(metric_array1[4]-metric_array[j_accepted,4])))#min(1,(exp(-y1+y))^(1/T_ann(j)))
        alpha_unbound <- (exp(metric_array1[4]-metric_array[j_accepted,4]))
      }
    }
    if(is.na(alpha_))
    {
      alpha_
    }
    ##may be good:
    #print('y,y1,alpha:')
    #print(c(y,y1))
    #print(c(y1-y,alpha_))
    # if(is.na(alpha_))
    #{
    #  print(alpha_)
    #}
    label_accept<-0
    if(alpha_==1)
    {
      j_accepted<- j_accepted+1
      candappadv[j_accepted,]<-x1
      results[j_accepted]<-y1
      metric_array[j_accepted,1:4] <- metric_array1[1:4]
      accepted[j]<-1
      naccepted[j_accepted]<-1
      nonmem_options$nLL_indiv<- ygoalf1$nLL_indiv
      residapp[j_accepted,]<-y2
      residapp_cumul[j_accepted,]<- ygoalf1$quadrresid_cumul
      label_accept <-1
    }else{
      add_guess <- runif(1)
      #generate uniformly distributed random number add_guess
      #("additional guess") and accept the candidate if it exeeds
      #alpha_; otherwise reject it and use the preceeding value.
      #Update the fields.
      if(add_guess<= alpha_)
      {
        j_accepted <- j_accepted+1
        candappadv[j_accepted,]<-x1
        results[j_accepted]<-y1
        metric_array[j_accepted,1:4] <- metric_array1[1:4]
        accepted[j]<-1
        naccepted[j_accepted]<-1
        nonmem_options$nLL_indiv<- ygoalf1$nLL_indiv
        residapp[j_accepted,]<-y2
        residapp_cumul[j_accepted,]<- ygoalf1$quadrresid_cumul
        label_accept <- 1
      }else{
        #candapp[j,]<-candapp[j-1,]
        naccepted[j_accepted]<-naccepted[j_accepted]+1
        #results[j]<-results[j-1]
        accepted[j]<-0
        #residapp[j,]<-residapp[j-1,]
        #residapp_cumul[j,]<- residapp_cumul[j-1,]
        label_accept <- 0
      }
    }
    x<-candappadv[j_accepted,]#starting center for the next iteration
    y <- results[j_accepted]
    
    xtot$phiopt<-generate_listvect(candappadv[j_accepted,],xtot$phiopt,main_struct)
    
    #print("results:")
    #if(j%in%print_ind)
    # {  
    #  print(c(j,j_accepted,results[j_accepted],naccepted[j_accepted],alpha_))#too long:,mean(accepted[1:j]))),Tann[j+1]
    # }
    if(label_accept==1)
    {  
      ind_prev_res <- ind_prev_res+1
    }
    if((j_accepted>min_steps_cov)&&(label_accept==1))#(burnin_cluster-2)
    {
      
      if((j_accepted==min_steps_cov)&(ind_prev_res==0))#initialization (burnin_cluster-1)
      {
        # ind_prev_res <- min_steps_cov#1+ind_prev_res#!!!#burnin_cluster-2
        acc_av <- apply(candappadv[1:min_steps_cov,],2,mean)#(burnin_cluster-1)
        acc_cov <- cov(candappadv[1:min_steps_cov,])#nonmem_options$prevcov#cov(candappadv[1:2,])
        nonmem_options$adapt_acc_cov <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
        #metric_array[1:j_accepted,2]<- determinant(nonmem_options$adapt_acc_cov)$modulus
        
        #metric_array[j_accepted,5]<- mean(metric_array[1:j_accepted,4])
        #metric_array[1:j_accepted,7]<- 0#metric_array[1:j_accepted,4]-metric_array[1:j_accepted,5]
        #metric_array[1:j_accepted,6]<- sd(metric_array[1:j_accepted,4])
        #metric_array[1:j_accepted,8]<-  metric_array[1:j_accepted,7]/metric_array[1:j_accepted,6]
        #metric_array[1:j_accepted,9]<- num_clusters
        #metric_array[1:j_accepted,10] <- 1:j_accepted
        cluster_str$cluster_values[num_clusters,1] <- max(metric_array[1:j_accepted,4])
        cluster_str$cluster_values[num_clusters,2] <- max(results[1:j_accepted])
        cluster_str$cluster_values[num_clusters,3] <- j_accepted
        cluster_str$cluster_values[num_clusters,4] <- 1###!!!!
        cluster_str$cluster_av[num_clusters,] <- acc_av
        cluster_str$cluster_cov[,,num_clusters] <- acc_cov
      }else
      {  
        # nonmem_options$prevres<-rbind(nonmem_options$prevres,candapp[subseti,])
        
        #nonmem_options$prevres[ind_prev_res+1,] <- candapp[j,]
        help_cov_av<-  RecursiveCovMatrixAvUpd(xn= candappadv[j_accepted,],#nonmem_options$prevres[ind_prev_res+1,],
                                               mup = acc_av,covp = acc_cov,
                                               nn = ind_prev_res)#ind_prev_res+1!!!???
        acc_av <- help_cov_av$mu
        acc_cov <-  help_cov_av$cov_
        #newres<-nonmem_options$prevres#rbind(nonmem_options$prevres,candapp[subseti,])
        #upd_cov<-FastCovAdv(arr=candapp[subseti,],prevcov=nonmem_options$proposed_var,steps_=nonmem_options$step_init)
        #nonmem_options$step_init <- nonmem_options$step_init+1
        
        
        #nonmem_options$proposed_var<- t(chol(cov(newres[1:ind_prev_res,])+diag(nonmem_options$eps_chol,dim(cov(newres)))))
        nonmem_options$adapt_acc_cov <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
        nonmem_options$proposed_var<- t(chol(nonmem_options$adapt_acc_cov))#acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
        
        #metric_array[j_accepted,1]<- MetricsMultivar(ve1=candappadv[j_accepted,],
        #                                             ve2=acc_av,#candappadv[j_accepted-1,],
        #                                             mat=nonmem_options$adapt_acc_cov)
        #metric_array[j_accepted,2]<- determinant(nonmem_options$adapt_acc_cov)$modulus
        #metric_array[j_accepted,3]<- 0.5*(metric_array[j_accepted,1]+metric_array[j_accepted,2])
        #metric_array[j_accepted,4] <- results[j_accepted]+metric_array[j_accepted,3]
        
        ##Not necessary now:
        #help_metric_cov_av<-  RecursiveCovMatrixAvUpd(xn= metric_array[j_accepted,4],#nonmem_options$prevres[ind_prev_res+1,],
        #                                              mup = metric_array[j_accepted-1,5],
        #                                              covp = metric_array[j_accepted-1,6],
        #                                              nn = min(max_av_cluster,ind_prev_res))#ind_prev_res+1!!!???
        #metric_array[j_accepted,5] <- help_metric_cov_av$mu
        #metric_array[j_accepted,6] <-  help_metric_cov_av$cov_^0.5
        #metric_array[j_accepted,7]<- metric_array[j_accepted,4]-metric_array[j_accepted,5]
        #metric_array[j_accepted,8]<-  metric_array[j_accepted,7]/metric_array[j_accepted,6]
        cond_new_cluster <- (alpha_unbound>alpha_clust)&&(ind_prev_res>(min_cluster_length+burnin_cluster) )
        if(ind_prev_res<5)
        {
          alpha_unbound0 <- alpha_unbound
        }else{
          alpha_unbound0 <- (exp(metric_array1[4]-cluster_str$cluster_values[num_clusters,1]))
        }
        
        #cond_new_cluster <- (alpha_unbound0>alpha_clust)&&(ind_prev_res>(min_cluster_length+burnin_cluster) )
        
        #cond_new_cluster <- ((metric_array[j_accepted,8]>1.644854)&&(ind_prev_res>(min_cluster_length+burnin_cluster) ))
        if(cond_new_cluster)#1.644854  1.959964
        {
          print(c(cluster_str$cluster_values[num_clusters,1],alpha_unbound0,alpha_unbound))
          ind_prev_res <- burnin_cluster
          #help_cov_av<-  RecursiveCovMatrixAvUpd(xn= candappadv[j_accepted,],#nonmem_options$prevres[ind_prev_res+1,],
          #                                        mup = acc_av,covp = acc_cov,
          #                                       nn = ind_prev_res)#ind_prev_res+1!!!???
          #acc_cov0 <-  help_cov_av$cov_
          #print(c(min(diag(acc_cov)^0.5),mean(diag(acc_cov)^0.5),max(diag(acc_cov)^0.5)))
          # print(c(min(diag(acc_cov0)^0.5),mean(diag(acc_cov0)^0.5),max(diag(acc_cov0)^0.5)))
          #acc_cov <-acc_cov0
          #acc_cov1 <- cov(candappadv[(j_accepted-1):j_accepted,])
          acc_av <- candappadv[j_accepted,]
          #acc_cov <- 0.2*acc_cov0+0.8*acc_cov
          num_clusters <- num_clusters+1
          metric_array[j_accepted,10] <- ind_prev_res-(burnin_cluster-1)
          #metric_array[j_accepted,11] <- (results[j_accepted])#Simple!!! without maximum 
          cluster_str$cluster_values[num_clusters,1] <- metric_array[j_accepted,4]
          cluster_str$cluster_values[num_clusters,2] <- results[j_accepted]#max(results[j_accepted],cluster_str$cluster_values[num_clusters,1])
          cluster_str$cluster_values[num_clusters,3] <- ind_prev_res-(burnin_cluster-1)
          cluster_str$cluster_values[num_clusters,4] <- 1###!!!!
        }else{
          metric_array[j_accepted,10] <- j_accepted
          #metric_array[j_accepted,11] <- (results[j_accepted])#Simple!!! without maximum 
          cluster_str$cluster_values[num_clusters,1] <- max(cluster_str$cluster_values[num_clusters,1],metric_array[j_accepted,4])
          cluster_str$cluster_values[num_clusters,2] <- max(results[j_accepted],cluster_str$cluster_values[num_clusters,1])
          cluster_str$cluster_values[num_clusters,3] <- ind_prev_res-(burnin_cluster-1)
          cluster_str$cluster_values[num_clusters,4] <- 1###!!!!
        }
        cluster_str$cluster_av[num_clusters,] <- acc_av
        cluster_str$cluster_cov[,,num_clusters] <- acc_cov
        metric_array[j_accepted,9]<- num_clusters
        if(main_struct$label_branch >1)
        {
          print(c(j,j_accepted,ind_prev_res,num_clusters,metric_array[j_accepted,4], results[j_accepted],alpha_unbound))#,alpha_unbound0))#c(4:8) metric_array[j_accepted,c(4:6)]
        }
        
      }
      
      
    }
    #}
    stop_cond <- (j>=nlocal)#((j>=nlocal)&&metric_array[j_accepted,8])
  }
  
  ###output:
  out_ <- list()
  # out_$cand <- cand
  out_$metric_array <- metric_array[1:j_accepted,]
  out_$results <-results[1:j_accepted]
  out_$total_accepted <- j_accepted
  out_$accepted <- accepted
  out_$naccepted <- naccepted[1:j_accepted]
  out_$Tann <- Tann
  out_$candapp <- candappadv[1:j_accepted,]
  out_$residapp <- residapp[1:j_accepted,]
  out_$residapp_cumul <- residapp_cumul[1:j_accepted,]
  out_$proposed_var <- nonmem_options$proposed_var
  out_$prevmean <- acc_av#help_cov_av$mu
  out_$prevcov <- acc_cov#help_cov_av$cov_
  out_$num_clusters <- num_clusters
  out_$cluster_str$cluster_values <- cluster_str$cluster_values[1:out_$num_clusters,]
  out_$cluster_str$cluster_av <- cluster_str$cluster_av[1:out_$num_clusters,]
  out_$cluster_str$cluster_cov <- cluster_str$cluster_cov[,,1:out_$num_clusters]
  out_$ind_prev_res <- ind_prev_res
  #out_$ <-
  ###
  out_
  
} 
AssigngClustWeights<-function(x,gam_mix)
{
  max_x <- max(x)
  we1 <- exp(x-max_x)^gam_mix#0.3#^0.3#0.75
  wen <- we1/sum(we1)
  return(wen)
}
AssigngClustNum<-function(rel_we)
{
  n_clust<- length(rel_we)
  func_guess <- runif(1)
  
  rel_wec <- 0#rel_we[1]
  ind_cl <- 0
  while((func_guess > rel_wec)&&(ind_cl <=n_clust))#for(ind_cl in 1: n_clust)
  {
    ind_cl <- ind_cl+1
    rel_wec <- rel_wec+rel_we[ind_cl]
    if (func_guess <= rel_wec)
    {
      i_clust <- ind_cl
    }
  }
  
  return(i_clust)
}
CalculatePi<-function(x_v,cluster_str,num_clusters)
{
  out_ <- 0
  for(ind_c in 1: num_clusters)
  {
    metric_array0_1_3 <- MetricArraySimpUpd(ve_prev=cluster_str$cluster_av[ind_c,],
                                            ve_curr=x_v,#!!!update the older solution 08.12.2020
                                            cov_matr=cluster_str$adapt_acc_cov[,,ind_c])
    out_ <- as.bigz(out_ + as.bigz(cluster_str$cluster_values[ind_c,'we']*exp(-metric_array0_1_3[3])))
  }
  if((out_==0)||(is.na(out_)||(is.infinite(out_))))#
  {
    out_1<- -(10^6)#very big: 1/0
  }
  else{
    out_1<- (-log(out_))
  }
  return(out_1)
}
ShowPi<-function(x_v,cluster_str,num_clusters)
{
  out_ <- array(dim=c(num_clusters,4))
  colnames(out_) = c('we','val','det','resid')
  for(ind_c in 1: num_clusters)
  {
    help_metr <- MetricArraySimpUpd(ve_prev=cluster_str$cluster_av[ind_c,],
                       ve_curr=x_v,#!!!update the older solution 08.12.2020
                       cov_matr=cluster_str$adapt_acc_cov[,,ind_c])
    out_[ind_c,'we'] <-cluster_str$cluster_values[ind_c,'we'] 
    out_[ind_c,'resid'] <- help_metr[1]#MetricsMultivar(ve1=x_v,
                      #ve2=cluster_str$cluster_av[ind_c,],
                     # mat=cluster_str$adapt_acc_cov[,,ind_c])
    out_[ind_c,'det'] <-  help_metr[2]#determinant(cluster_str$adapt_acc_cov[,,ind_c])$modulus
    out_[ind_c,'val'] <-  help_metr[3]#
    
  }
  return(out_)
}
CalculatePiClever<-function(x_v,cluster_str,num_clusters,ind_cluster,label_complete)
{
  if(label_complete==1)
   { 
    out_ <- 0
    log_clust_dens<- vector(length = num_clusters)
    for(ind_c in 1: num_clusters)
    {
      metric_array0_1_3 <- MetricArraySimpUpd(ve_prev=cluster_str$cluster_av[ind_c,],
                                              ve_curr=x_v,#!!!update the older solution 08.12.2020
                                              cov_matr=cluster_str$adapt_acc_cov[,,ind_c])
      if(cluster_str$cluster_values[ind_c,'we']==0)
      {
        log_clust_dens[ind_c] <- NA
      }else{
        log_clust_dens[ind_c]<-log(cluster_str$cluster_values[ind_c,'we'])-metric_array0_1_3[3]
      }
      
      #out_ <- as.bigz(out_ + as.bigz(cluster_str$cluster_values[ind_c,'we']*exp(-metric_array0_1_3[3])))
    }
    max_log_val<-max(log_clust_dens,na.rm=T)
    relev_ind1<- !is.na(log_clust_dens)
    if(sum(relev_ind1)>0)
    {  
      relev_ind<- log_clust_dens[relev_ind1]>max_log_val+(-10*log(10))
      log_clust_dens_relev<-log_clust_dens[relev_ind1][relev_ind]
      min_log_val<- min(log_clust_dens_relev)
      term_ <- sum(exp(log_clust_dens_relev-min_log_val))
      out_1 <- -(min_log_val+log(term_))
    }
    else{
      out_1 <-0
    }
  }else if(label_complete==0){
    out_1 <- 0.5*MetricsMultivar(ve1=cluster_str$cluster_av[ind_cluster,],
                             ve2=x_v,
                             mat=cluster_str$adapt_acc_cov[,,ind_cluster])
  }else if(label_complete==-1)
  {
    out_1 <- MetricArraySimpUpd(ve_prev=cluster_str$cluster_av[ind_cluster,],
                               ve_curr=x_v,#!!!update the older solution 08.12.2020
                               cov_matr=cluster_str$adapt_acc_cov[,,ind_cluster])[3]
  }
  return(out_1)
}
StartNewChaincond<-function(metr_val,num_clusters,cluster_str)
{
  out_ <- 0
  if(num_clusters>1)
  {
    plausib_vals <-cluster_str$cluster_values[1:num_clusters,'adjust_clust_val']>=(metr_val-log(0.1))
    if((sum(cluster_str$cluster_values[1:num_clusters,'we'][plausib_vals])<0.05)&&sum(plausib_vals)>0)
    {
      out_ <- 1
    }
  }
  
  return(out_)
}
MetricArraySimpUpd<-function(ve_prev,ve_curr,cov_matr)
{
  metric_vect <- vector(length=3)
  
  metric_vect[1]<- MetricsMultivar(ve1=ve_curr,
                                   ve2=ve_prev,
                                   mat=cov_matr)
  metric_vect[2]<- determinant(cov_matr)$modulus
  metric_vect[3]<- 0.5*(metric_vect[1]+metric_vect[2])
  
  
  
  return(metric_vect)
}
UpdateClusters<-function(cluster_str, ind_indiv, xtot,  main_struct, gam_mix)
{
  label_complete <- 0
  num_clusters <- dim(cluster_str$cluster_values)[1]
  for(ind_c in 1: num_clusters)
  {
    
    if(ind_indiv==0)
    {  
      xtot$phiopt<-generate_listvect(cluster_str$cluster_av[ind_c,],xtot$phiopt,main_struct)
    }
    else{  
      xtot$phiopt<-generate_listvCOVID(cluster_str$cluster_av[ind_c,],ind_indiv,xtot$phiopt,main_struct)
    }
    if(main_struct$num_estim_indiv[ind_indiv]==0)
    {  
      xtot<-InitUpdate(xtot,main_struct)
    }else{
      xtot<-InitUpdateMultiple(xtot,ind_indiv,main_struct)
    }
    #xtot<-InitUpdate(xtot,main_struct)
    ygoalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                          nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
    y1<- - ygoalf$nLL_tot
    print(c(ind_c,num_clusters,y1))
    metric_array1<- vector(length = 4)
    metric_array1[3] <- CalculatePiClever(x_v=cluster_str$cluster_av[ind_c,],cluster_str,num_clusters,ind_c,label_complete)
    metric_array1[4] <- y1+metric_array1[3]
    cluster_str$cluster_values[ind_c,1] <- metric_array1[4]
    cluster_str$cluster_values[ind_c,2] <- y1 
    cluster_str$cluster_values[ind_c,'adjust_clust_val']<-cluster_str$cluster_values[ind_c,1] 
  }
  if(num_clusters>1)
  {  
  cluster_str$cluster_values[1:num_clusters,'we'] <- AssigngClustWeights(x = cluster_str$cluster_values[1:ind_c,2],gam_mix)
  }
  return(cluster_str)
}
SelectPlausibleSamples<-function(resultsv,delt_plausible,min_sample_fr,min_sample_size)
{
  max_res <- max(resultsv)
  sortv<- sort(resultsv,decreasing = T)
  res_1 <- sortv[min_sample_size]
  res_2 <- sortv[ceiling(min_sample_fr*length(resultsv))]
  plausible_ind <- (resultsv>=(max(resultsv)-delt_plausible))|(resultsv>=res_1)|(resultsv>=res_2)
  return(plausible_ind)
}
ClustHistShort<-function(clust_obj,digi)
{
  num_clust<-length(clust_obj$cluster_values[,'hist'])
  for(ind_clust in 1:num_clust)
  {
    helpi <- clust_obj$cluster_values[ind_clust,'hist']/(10^digi)
    if(helpi>1)
    { 
    clust_obj$cluster_values[ind_clust,'hist']<- (helpi-floor(helpi))*(10^digi)
    }
  }
  return(clust_obj)
}