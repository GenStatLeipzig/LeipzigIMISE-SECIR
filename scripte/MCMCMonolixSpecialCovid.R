MCMCDomesticHaarFaststep<-function(x,ind_indiv,xtot,main_struct,nonmem_options,num_var,ind_prev_res,nlocal)#,compl_structindiv
{
  alpha_clust <- exp(10)#3
  const_term<-num_var/2*log(2*pi)
  limi <- -log(0.05)/0.5
  metric_array<- array(dim =  c(nlocal,6) )#-1
  #print_ind <-  1000*(1:10)#1000*(1:(nonmem_options$m/1000)) 
  xtot<- UpdateParam(x=x,ind_indiv,xtot,stoch_obj0,main_struct0)
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
  ind_iter <-1 
  #xcenter <- candapp[1:j-1,]
  j_accepted <- 1
  for(j in 2:nlocal)
  {
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
    label_opt1="pop"
    ############
      if(nonmem_options$from_mu==1)
      {
        x1<-ProposalDistInit(acc_av,nonmem_options$proposed_var*Tann[j])
      }else{
        x1<-ProposalDistInit(x,nonmem_options$proposed_var*Tann[j])
      }
      
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
      
      if(is.na(y1))
      {
        alpha_<-0
        alpha_unbound<-0
      }else{
        alpha_ <- min(1,(exp(y1-y)))#min(1,(exp(-y1+y))^(1/T_ann(j)))
        alpha_unbound <- exp(y1-y)
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
      if((j_accepted==2))
      {
        ind_prev_res <- ind_prev_res+1
        if(ind_prev_res==1)
        {  
          
          acc_av <- apply(candappadv[1:j_accepted,],2,mean)
          acc_cov <- nonmem_options$prevcov#cov(candappadv[1:j_accepted,])
        }
        if((ind_prev_res>nonmem_options$ignore)&(main_struct$label_branch >0))
        {  
          #nonmem_options$proposed_var<- t(chol(cov(newres[1:ind_prev_res,])+diag(nonmem_options$eps_chol,dim(cov(newres)))))
          adapt_acc_cov0 <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
          proposed_var0 <- t(chol(adapt_acc_cov0))#acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
          
          metric_array[j_accepted,1]<- MetricsMultivar(ve1=candappadv[j_accepted,],
                                                       ve2=candappadv[j_accepted-1,],
                                                       mat= adapt_acc_cov0)
          metric_array[j_accepted,2]<- determinant(adapt_acc_cov0)$modulus
          metric_array[j_accepted,3]<- 0.5*(metric_array[j_accepted,1]+metric_array[j_accepted,2])
          metric_array[j_accepted,4] <- results[j_accepted]+metric_array[j_accepted,3]
          j_rel <- j_accepted-1
          if(j_accepted==2)
          {
            metric_array[j_accepted,5]<-metric_array[j_accepted,4]
          }else{
            metric_array[j_accepted,5] <- metric_array[j_accepted-1,5]*(j_rel-1)/j_rel+metric_array[j_accepted,4]/j_rel
          }
          metric_array[j_accepted,6]<-metric_array[j_accepted,4]-metric_array[j_accepted,5]
          if(main_struct$label_branch >1)
          {
            print(c(ind_prev_res,j_accepted,metric_array[j_accepted,]))
          }
        }
      }  
      else if(j_accepted>2)#&&(label_accept==1))
      {#((j_accepted>3)&&(label_accept==1))
        #if((j_accepted==4)&(ind_prev_res==0))#initialization
        #{
        #  ind_prev_res <- 3+ind_prev_res#!!!
        #  acc_av <- apply(candappadv[1:4,],2,mean)
        #  acc_cov <- cov(candappadv[1:4,])
          
        #}else
        #{  
          # nonmem_options$prevres<-rbind(nonmem_options$prevres,candapp[subseti,])
        
          ind_prev_res <- ind_prev_res+1
         # acc_av <- candappadv[j_accepted,]
         # help_cov_av<-  RecursiveCovMatrixAvUpd(xn= candappadv[j_accepted,],#nonmem_options$prevres[ind_prev_res+1,],
          #                                       mup = acc_av,covp = acc_cov,
        #                                         nn = ind_prev_res)
       #   acc_cov <-  help_cov_av$cov_
          #nonmem_options$prevres[ind_prev_res+1,] <- candapp[j,]
         # if(alpha_clust < alpha_unbound)
         # {
            #ind_prev_res <-5 #Fuck!
         #   acc_av <- candappadv[j_accepted,]
         #   help_cov_av<-  RecursiveCovMatrixAvUpd(xn= candappadv[j_accepted,],#nonmem_options$prevres[ind_prev_res+1,],
          #                                         mup = acc_av,covp = acc_cov,
          #                                         nn = ind_prev_res+1)
         #   acc_cov <-  help_cov_av$cov_
          #}else
      #    {  
          if('mem'%in%names(nonmem_options))
          {
            help_cov_av<-  RecursiveCovMatrixAvUpd(xn= candappadv[j_accepted,],#nonmem_options$prevres[ind_prev_res+1,],
                                                   mup = acc_av,covp = acc_cov,
                                                   nn = min(nonmem_options$mem,ind_prev_res+1))#ind_prev_res+1!!!???
          }
          else{
            help_cov_av<-  RecursiveCovMatrixAvUpd(xn= candappadv[j_accepted,],#nonmem_options$prevres[ind_prev_res+1,],
                                                   mup = acc_av,covp = acc_cov,
                                                   nn = ind_prev_res+1)#ind_prev_res+1!!!???
          }
          #test_num <- max(abs(help_cov_av$mu-apply(candapp[1:j,],2,mean)))
          #if(test_num>0.0000000001)
          #{
          #  print(c(j,test_num))
          #}
          
          acc_av <- help_cov_av$mu
          acc_cov <-  help_cov_av$cov_
          nonmem_options$adapt_acc_cov <- acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
          nonmem_options$proposed_var<- t(chol(nonmem_options$adapt_acc_cov))#acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))
          
     #   }
          
        #newres<-nonmem_options$prevres#rbind(nonmem_options$prevres,candapp[subseti,])
        #upd_cov<-FastCovAdv(arr=candapp[subseti,],prevcov=nonmem_options$proposed_var,steps_=nonmem_options$step_init)
        #nonmem_options$step_init <- nonmem_options$step_init+1
        if((ind_prev_res>nonmem_options$ignore)&(main_struct$label_branch >0))
        {  
          #nonmem_options$proposed_var<- t(chol(cov(newres[1:ind_prev_res,])+diag(nonmem_options$eps_chol,dim(cov(newres)))))
          
          metric_array[j_accepted,1]<- MetricsMultivar(ve1=candappadv[j_accepted,],
                                                      ve2=candappadv[j_accepted-1,],
                                                      mat=nonmem_options$adapt_acc_cov)
          metric_array[j_accepted,2]<- determinant(nonmem_options$adapt_acc_cov)$modulus
          metric_array[j_accepted,3]<- 0.5*(metric_array[j_accepted,1]+metric_array[j_accepted,2])
          metric_array[j_accepted,4] <- results[j_accepted]+metric_array[j_accepted,3]
          j_rel <- j_accepted-1#j_accepted-3
          #if(j_accepted==4)
          #{
         #   metric_array[j_accepted,5]<-metric_array[j_accepted,4]
          #}else{
            metric_array[j_accepted,5] <- metric_array[j_accepted-1,5]*(j_rel-1)/j_rel+metric_array[j_accepted,4]/j_rel
          #}
          metric_array[j_accepted,6]<-metric_array[j_accepted,4]-metric_array[j_accepted,5]
          if(main_struct$label_branch >1)
          {
            print(c(j,ind_prev_res,j_accepted,metric_array[j_accepted,], results[j_accepted]))
          }
        }
      }
      }
    #}
    
  }
  
  ###output:
  out_ <- list()
 # out_$cand <- cand
  out_$metric_array <- metric_array[1:j_accepted,]
  out_$results <-results[1:j_accepted]
  out_$accepted <- accepted
  out_$naccepted <- naccepted[1:j_accepted]
  out_$Tann <- Tann
  out_$candapp <- candappadv[1:j_accepted,]
  out_$residapp <- residapp[1:j_accepted,]
  out_$residapp_cumul <- residapp_cumul[1:j_accepted,]
  out_$proposed_var <- nonmem_options$proposed_var
  out_$prevmean <- acc_av#help_cov_av$mu
  out_$prevcov <- acc_cov#help_cov_av$cov_
  out_$total_accepted <- j_accepted
  #out_$ <-
  ###
  out_
  
} 
####
MCMCDomesticHaarGlobal<-function(x,ind_indiv,xtot,main_struct,nonmem_options,num_var,nlocal,nglobal,ind_prev_res)#,compl_structindiv
{
  #nonmem_options$prevcov <- precov
  ntot<- nlocal*nglobal
  print('global steps =')
  print(nglobal)
  print('local steps =')
  print(nlocal)
  ptm <- proc.time()
  localMCMCi <- MCMCDomesticHaarFaststep(x=x,ind_indiv,xtot,main_struct,nonmem_options,num_var=num_var,
                                         ind_prev_res=ind_prev_res,nlocal=nlocal)#,compl_structindiv
  print((proc.time() - ptm))
  write.table(x=localMCMCi$candapp, sep = ";" ,  dec = main_struct0$dec, 
              row.names=F,col.names = F,file=paste(main_struct0$resultspath,
                                                   "MCMCStore.csv",sep=''), append = F)
  write.table(x=cbind(localMCMCi$naccepted,localMCMCi$results),sep = ";" ,  dec = main_struct0$dec, 
              row.names=F,col.names = F,file=paste(main_struct0$resultspath,
                                                   "ResultsStore.csv",sep=''), append = F)
  
  nonmem_options$prevcov<-localMCMCi$prevcov
  nonmem_options$prevmean <- localMCMCi$prevmean
  nonmem_options$proposed_var <- localMCMCi$proposed_var
  diagi <- diag(localMCMCi$prevcov)
  print(c(min(diagi),mean(diagi),max(diagi))^0.5)
  diagi <- diag(localMCMCi$proposed_var)
  print(c(min(diagi),mean(diagi),max(diagi))^0.5)
  glob_results <- localMCMCi$results
  glob_naccepted  <- localMCMCi$naccepted
  glob_candapp  <- localMCMCi$candapp
  glob_residapp  <- localMCMCi$residapp
  glob_residapp_cumul  <- localMCMCi$residapp_cumul
  glob_metric_array <- localMCMCi$metric_array 
  ind_prev_resi<-0  
  #numi_local_steps<-vector(length=ntot)
  #numi_local_steps[1]<-length(localMCMCi$results)
  for(ind_local in 1:nglobal)#251
  {
    ind_prev_resi <- ind_prev_resi+length(localMCMCi$results)
    ind_prev_res <- ind_prev_res+length(localMCMCi$results)
    if(localMCMCi$total_accepted>1)#(length(localMCMCi$results)>1)
    {
      ptm <- proc.time()
      localMCMCi <- MCMCDomesticHaarFaststep(x=localMCMCi$candapp[length(localMCMCi$results),],
                                             ind_indiv,xtot,main_struct0,nonmem_options,num_var=num_var,
                                             ind_prev_res=ind_prev_res,nlocal=nlocal)#,compl_structindiv
      print(c(ind_local,(proc.time() - ptm),length(localMCMCi$results),ind_prev_resi))
      glob_results <- c(glob_results,localMCMCi$results)
      glob_naccepted  <- c(glob_naccepted,localMCMCi$naccepted)
      glob_candapp  <- rbind(glob_candapp,localMCMCi$candapp)
      glob_residapp  <- rbind(glob_residapp, localMCMCi$residapp)
      glob_residapp_cumul  <- rbind(glob_residapp_cumul,localMCMCi$residapp_cumul)
      glob_metric_array <- rbind(glob_metric_array,localMCMCi$metric_array) 
      nonmem_options$proposed_var <- localMCMCi$proposed_var
      #numi_local_steps[ind_local+1]<-length(localMCMCi$results)
      nonmem_options$prevcov<-localMCMCi$prevcov
      nonmem_options$prevmean <- localMCMCi$prevmean
      write.table(x=localMCMCi$candapp, sep = ";" ,  dec = main_struct0$dec, 
                  row.names=F,col.names = F,file=paste(main_struct0$resultspath,
                                                       "MCMCStore.csv",sep=''), append = T)
      write.table(x=cbind(localMCMCi$naccepted,localMCMCi$results),sep = ";" ,  dec = main_struct0$dec, 
                  row.names=F,col.names = F,file=paste(main_struct0$resultspath,
                                                       "ResultsStore.csv",sep=''), append = T)
      diagi <- diag(localMCMCi$prevcov)
      print(c(min(diagi),mean(diagi),max(diagi))^0.5)
      diagi <- diag(localMCMCi$proposed_var)
      print(c(min(diagi),mean(diagi),max(diagi))^0.5)
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
  return(output_)
}
  ####



MCMCDomesticHaarGeneral<-function(ind_indiv,xtot,main_struct,nonmem_options)#,compl_structindiv
{
  print_ind <-  1000*(1:10)#1000*(1:(nonmem_options$m/1000)) 
  if(main_struct$numpoppar>0)#(ind_indiv==0)#population only
  {  
    x<-generate_pop_vector(xtot$phiopt,main_struct)
  }
  else{
    x<-xtot$phiopt$indiv[[ind_indiv]]#generate_listvCOVID(x,ind_indiv,xtot,main_struct)
  }
  num_var <- length(x)
  if(is.null(nonmem_options$prevres))
  {
    acc_av <- x
    acc_cov <- array(data = 0, dim=c(num_var,num_var))
    nonmem_options$prevres<-array(dim = c(ceiling(nonmem_options$m/nonmem_options$jump)+1,num_var))
    ind_prev_res<-0
  }else{
    ind_prev_res<- dim(nonmem_options$prevres)[1]
    nonmem_options$prevres<- rbind(nonmem_options$prevres,array(dim = c(ceiling(nonmem_options$m/nonmem_options$jump)+1,num_var)))
    acc_av <- nonmem_options$prevmean
    acc_cov <- nonmem_options$prevcov
  }
  if('sd0'%in%names(nonmem_options))
  {
    nonmem_options$sd <- (nonmem_options$sd0/num_var)^0.5#Chol is equivalent to square root!!!
    nonmem_options$update_tann<-1
    nonmem_options$Tanninit<-nonmem_options$sd
  }else{
    nonmem_options$update_tann<-1.1
    nonmem_options$Tanninit <- 1
  }
    
  ##num_chains is currently 1
  #xx<-generate_pop_vector(xtot,main_struct1)
  cand <- array(dim=c(nonmem_options$m,num_var))
  candapp <- array(dim=c(nonmem_options$m,num_var))
  results <- vector(length = nonmem_options$m)
  accepted <- vector(length = nonmem_options$m)
  Tann <- rep.int(x=nonmem_options$Tanninit,times=nonmem_options$m+1)#1
  #compl_struct <- GenerateMCMCSpecStruct(xtot,main_struct,nonmem_options,candapp,results,accepted)
  print(length(Tann))
  #cand[1,]<-x
  candapp[1,]<-x
  Tann[1] <- nonmem_options$Tanninit#1
  accepted[1]<- 1
  #y<-func_pop_thrombo(x,xtot,main_struct)
  #ind_indiv <- 1
  ygoalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                          nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1= ind_indiv)
  y<- (- ygoalf$nLL_tot)
  y20<-  ygoalf$quadrresid
  residapp <- array(dim=c(nonmem_options$m,length(ygoalf$quadrresid)))
  residapp_cumul <- array(dim=c(nonmem_options$m,length(ygoalf$quadrresid)))
  residapp[1,]<-y20
  residapp_cumul[1,] <- ygoalf$quadrresid_cumul
  nonmem_options$nLL_indiv<-ygoalf$nLL_indiv
  results[1] <- y
  ind_iter <-1 
  #xcenter <- candapp[1:j-1,]
  for(j in 2:nonmem_options$m)
  {
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
    label_opt1="pop"
    ############
    #if(j==2)
    #{
    #  x1<-ProposalDistInit(x,nonmem_options$proposed_var*Tann[j])
   # }else{
   #   x1<-ProposalDistrMultiv(prev_chain=candapp[1:j-1,],omega_=nonmem_options$proposed_var*Tann[j],nonmem_options$ignore,opt_=j-10*floor(j/10))
   # } 
    if(((j-2)==round((j-2)/20)*20)&(j>nonmem_options$ignore)&(nonmem_options$num_indiv_runs>0))
    {
      compl_structindiv <- GenerateMCMCSpecIndivStruct(main_struct,nonmem_options)
      #print("individual separately:")
      #WriteMCMCResults(output_mcmc=outpop,main_struct,dir_name=nonmem_options$dir_name)
      mcmc_indiv_output<- list()
      #if(ind_iter<nonmem_options$num_indiv_runs)
      #{
        for(ind_indiv in 1:main_struct$numindiv)
        {
          mcmc_indiv_output[[ind_indiv]] <- MCMCchain(xtot,main_struct,type_mcmc=ind_indiv,compl_struct=compl_structindiv,nonmem_options,ind_iter=1)
          xtot<-update_listindivvect(mcmc_indiv_output[[ind_indiv]]$candapp[nonmem_options$num_steps$indiv,],xtot,main_struct,ind_indiv)#generate_listvect(x1,xtot,main_struct)
          nonmem_options$nLL_indiv <- mcmc_indiv_output[[ind_indiv]]$nLL_indiv
          
        }
        x<- generate_pop_vector(xtot,main_struct)
        x1<- generate_pop_vector(xtot,main_struct)
        ygoalf1<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                              nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
        y1 <- - ygoalf1$nLL_tot
        y <- y1
        y2<- ygoalf1$quadrresid
        residapp[j,]<-y2
        residapp_cumul[j,]<- ygoalf1$quadrresid_cumul
        candapp[j,]<-x1
        results[j]<-y1
        accepted[j]<-1
        nonmem_options$nLL_indiv<- ygoalf1$nLL_indiv
        #nonmem_options$prevres<-rbind(nonmem_options$prevres,x2)
        ind_iter <- ind_iter +1
        alpha_ <- 1#not determined in this case
      #}
      
    }else{
      x1<-ProposalDistInit(x,nonmem_options$proposed_var*Tann[j])
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
      
      if(is.na(y1))
      {
        alpha_<-0
      }else{
         alpha_ <- min(1,(exp(y1-y)))#min(1,(exp(-y1+y))^(1/T_ann(j)))
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
      if(alpha_==1)
      {
        candapp[j,]<-x1
        results[j]<-y1
        accepted[j]<-1
        nonmem_options$nLL_indiv<- ygoalf1$nLL_indiv
        residapp[j,]<-y2
        residapp_cumul[j,]<- ygoalf1$quadrresid_cumul
      }else{
        add_guess <- runif(1)
        #generate uniformly distributed random number add_guess
        #("additional guess") and accept the candidate if it exeeds
        #alpha_; otherwise reject it and use the preceeding value.
        #Update the fields.
        if(add_guess<= alpha_)
        {
          candapp[j,]<-x1
          results[j]<-y1
          accepted[j]<-1
          nonmem_options$nLL_indiv<- ygoalf1$nLL_indiv
          residapp[j,]<-y2
          residapp_cumul[j,]<- ygoalf1$quadrresid_cumul
        }else{
          candapp[j,]<-candapp[j-1,]
          results[j]<-results[j-1]
          accepted[j]<-0
          residapp[j,]<-residapp[j-1,]
          residapp_cumul[j,]<- residapp_cumul[j-1,]
        }
      }
      x<-candapp[j,]#starting center for the next iteration
      y <- results[j]
    }
    
    
    xtot$phiopt<-generate_listvect(candapp[j,],xtot$phiopt,main_struct)
    ###time consuming, currently not used:
   # if(j>nonmem_options$temp_check)
   # {
   #  # print("annealing criteria:")
   #  # print(mean(accepted[(j-nonmem_options$temp_check+1):j]))
   #   if((mean(accepted[(j-nonmem_options$temp_check+1):j])<nonmem_options$acceptl)&(Tann[j]>nonmem_options$minvar))
   #   {
   #     Tann[(j+1):nonmem_options$m]<-Tann[j]/nonmem_options$update_tann#1.1#1.2
   #   }
   #   if(mean(accepted[(j-nonmem_options$temp_check+1):j])>nonmem_options$acceptu)
  #    {
   #     Tann[(j+1):nonmem_options$m]<-Tann[j]*nonmem_options$update_tann#1.1#1.2
   #   }
      
   # }
    #print("results:")
    if(j%in%print_ind)
    {  
    print(c(j,results[j],accepted[j],alpha_))#too long:,mean(accepted[1:j]))),Tann[j+1]
    }
    if(nonmem_options$write == 1)
    {
      if(j==round(j/20)*20)
      {
        print("writing")
        #WriteMCMCResults(output_mcmc=outpop,main_struct,dir_name=nonmem_options$dir_name)
        WriteMCMCIntermResults(j,results,accepted,Tann,candapp,main_struct,dir_name=nonmem_options$dir_name)
        mcmc_indiv_output<- list()
        
      }
    } 
    if((j==round(j/nonmem_options$jump)*nonmem_options$jump)&(j>3*nonmem_options$jump))
    {
      if(nonmem_options$jump>1)
      {  
        subseti<-SubsetVect(1:j,step_=nonmem_options$jump)
      }else{
        subseti<-1:j
      }
      if((j==4*nonmem_options$jump)&(ind_prev_res==0))#initialization
      {
        nonmem_options$prevres[1:4,]<- candapp[subseti,]
        
        ind_prev_res <- 3+ind_prev_res#!!!
        acc_av <- apply(nonmem_options$prevres[1:4,],2,mean)
        acc_cov <- cov(nonmem_options$prevres[1:4,])
        
      }else
      {  
     # nonmem_options$prevres<-rbind(nonmem_options$prevres,candapp[subseti,])
        ind_prev_res <- ind_prev_res+1
        nonmem_options$prevres[ind_prev_res+1,] <- candapp[j,]
        help_cov_av<-  RecursiveCovMatrixAvUpd(xn=nonmem_options$prevres[ind_prev_res+1,],
                                               mup = acc_av,covp = acc_cov,
                                               nn = ind_prev_res+1)
        #test_num <- max(abs(help_cov_av$mu-apply(candapp[1:j,],2,mean)))
        #if(test_num>0.0000000001)
        #{
        #  print(c(j,test_num))
        #}
        acc_av <- help_cov_av$mu
        acc_cov <-  help_cov_av$cov_
      }
      if('mem'%in%names(nonmem_options))
      {
        mem_step <- ceiling(nonmem_options$mem/nonmem_options$jump)
        if((dim(nonmem_options$prevres)[1])-mem_step+1>0)
        {  
         nonmem_options$prevres<-nonmem_options$prevres[((dim(nonmem_options$prevres)[1])-mem_step+1):dim(nonmem_options$prevres)[1],]
        }
      }
      newres<-nonmem_options$prevres#rbind(nonmem_options$prevres,candapp[subseti,])
      #upd_cov<-FastCovAdv(arr=candapp[subseti,],prevcov=nonmem_options$proposed_var,steps_=nonmem_options$step_init)
      #nonmem_options$step_init <- nonmem_options$step_init+1
      if(ind_prev_res>nonmem_options$ignore)
      {  
        #nonmem_options$proposed_var<- t(chol(cov(newres[1:ind_prev_res,])+diag(nonmem_options$eps_chol,dim(cov(newres)))))
        nonmem_options$proposed_var<- t(chol(acc_cov+diag(nonmem_options$eps_chol,c(num_var,num_var))))
        
      }
    }
    
    
  }
  
  ###output:
  out_ <- list()
  out_$cand <- cand
  out_$results <-results
  out_$accepted <- accepted
  out_$Tann <- Tann
  out_$candapp <- candapp
  out_$residapp <- residapp
  out_$residapp_cumul <- residapp_cumul
  out_$proposed_var <- nonmem_options$proposed_var
  out_$prevmean <- help_cov_av$mu
  out_$prevcov <- help_cov_av$cov_
  #out_$ <-
  ###
  out_
   
}

MCMCchain<-function(xtot,main_struct,type_mcmc,compl_struct,nonmem_options,ind_iter)
{
  out_ <- list()
  if(type_mcmc==0)#population
  {
    label_opt1="pop"
    x0 <- generate_pop_vector(xtot,main_struct)
    Tann<-compl_struct$Tann
    num_steps<-nonmem_options$num_steps
    proposed_var<- compl_struct$proposed_var
    candapp <- compl_struct$candapp
    results <- compl_struct$results
    accepted <- compl_struct$accepted
    prevres <- compl_struct$prevres
    ind_indiv<-1
  }else{
    ind_runi<-((ind_iter-1)*nonmem_options$num_steps$indiv+1):(ind_iter*nonmem_options$num_steps$indiv+1)
    label_opt1="indiv"
    ind_indiv<- type_mcmc
    x0<- generate_indiv_vector (xtot,main_struct,ind_indiv)
    Tann<-compl_struct$Tann[,ind_indiv]
    num_steps<-nonmem_options$num_steps$indiv
    proposed_var<- compl_struct$proposed_var[[ind_indiv]]
    candapp <- compl_struct$candapp[[ind_indiv]][ind_runi,]
    results <- compl_struct$results[ind_runi,ind_indiv]
    accepted <- compl_struct$accepted[ind_runi,ind_indiv]
    prevres <- compl_struct$prevres[[ind_indiv]]
  }
  x<-x0
  #if('sd0'%in%names(nonmem_options))
  #{
  #  nonmem_options$sd <- nonmem_options$sd0/num_var
  #  nonmem_options$update_tann<-1.1
  #  nonmem_options$Tanninit<-nonmem_options$sd
 # }else{
 #   nonmem_options$update_tann<-1
 #   nonmem_options$Tanninit <- 1
 # }
  if('sd0'%in%names(nonmem_options))
  {
    nonmem_options$sd <- (nonmem_options$sd0/length(x))^0.5#Chol is equivalent to square root!!!
    nonmem_options$update_tann<-1
    nonmem_options$Tanninit<-nonmem_options$sd
  }else{
    nonmem_options$update_tann<-1.1
    nonmem_options$Tanninit <- 1
  }
  candapp[1,]<-x
  Tann[1:length(Tann)] <- nonmem_options$Tanninit
  accepted[1]<- 1
  #y<-func_pop_thrombo(x,xtot,main_struct)
  
  ygoalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt=label_opt1,
                       nLL_indiv=nonmem_options$nLL_indiv,ind_indiv1= ind_indiv)
  y<- (-ygoalf$nLL_tot)
  results[1] <- y
  out_$nLL_indiv<-  ygoalf$nLL_indiv
  for(j in 2:num_steps)
  {
    #####
    
    
    ############
    #if(j==2)
    #{
      x1<-ProposalDistInit(x,proposed_var*Tann[j])
    #}else{
    #  x1<-ProposalDistrMultiv(prev_chain=candapp[1:j-1,],omega_=proposed_var*Tann[j],nonmem_options$ignore,opt_=j-10*floor(j/10))
    #}
     
    #cand[j,]<-x1
    #y1<-func_pop_thrombo(x1,xtot,main_struct)
    xtot1 <- update_listindivvect(x1,xtot,main_struct,ind_indiv)#generate_listvect(x1,xtot,main_struct)
    #xtot1<-InitUpdate(xtot1,main_struct)
    if(main_struct$num_estim_indiv[ind_indiv]==0)
    {  
      xtot1<-InitUpdate(xtot1,main_struct)
    }else{
      xtot1<-InitUpdateMultiple(xtot1,ind_indiv,main_struct)
    }
    
    
    ygoalf1<-GoalFuncComp(xtot=xtot1,main_struct=main_struct,opt_mod= "new",label_opt=label_opt1,
                          nLL_indiv=nonmem_options$nLL_indiv,ind_indiv1=ind_indiv)
    y1<-  (- ygoalf1$nLL_tot)
    y2<- (- ygoalf1$quadrresid)
    alpha_ <- min(1,(exp(y1-y)))#min(1,(exp(-y1+y))^(1/T_ann(j)))
    print('y,y1,alpha:')
    print(c(y,y1))
    print(c(y1-y,alpha_)) 
    
    if(is.na(alpha_))
    {
      print(alpha_)
    }
    if(alpha_==1)
    {
      candapp[j,]<-x1
      results[j]<-y1
      accepted[j]<-1
      out_$nLL_indiv <-  ygoalf1$nLL_indiv
      residapp[j,] <- y2
      residapp_cumul[j,]<- ygoalf1$quadrresid_cumul
    }else{
      add_guess <- runif(1)
      #generate uniformly distributed random number add_guess
      #("additional guess") and accept the candidate if it exeeds
      #alpha_; otherwise reject it and use the preceeding value.
      #Update the fields.
      if(add_guess<= alpha_)
      {
        candapp[j,]<-x1
        results[j]<-y1
        accepted[j]<-1
        out_$nLL_indiv <-  ygoalf1$nLL_indiv
        residapp[j,] <- y2
        residapp_cumul[j,]<- ygoalf1$quadrresid_cumul
      }else{
        candapp[j,]<-candapp[j-1,]
        results[j]<-results[j-1]
        accepted[j]<-0
        residapp[j,] <- residapp[j-1,]
        residapp_cumul[j,] <- residapp_cumul[j-1,]
      }
    }
    x <- candapp[j,]#for the next iteration!!!!
    y <- results[j]
    if(j>nonmem_options$temp_check)
    {
      print("annealing criteria:")
      print(mean(accepted[(j-nonmem_options$temp_check+1):j]))
      if((mean(accepted[(j-nonmem_options$temp_check+1):j])<nonmem_options$acceptl)&(Tann[j]>nonmem_options$minvar))
      {
        Tann[(j+1):num_steps]<-Tann[j]/nonmem_options$update_tann#1.1#1.2
      }
      if(mean(accepted[(j-nonmem_options$temp_check+1):j])>nonmem_options$acceptu)
      {
        Tann[(j+1):num_steps]<-Tann[j]*nonmem_options$update_tann#1.1#1.2
      }
      
    }
    print("results:")
    print(c(j,Tann[j+1],results[j],accepted[j],alpha_,mean(accepted[1:j])))
 ###   if(nonmem_options$write == 1)
 ###  {
  ###    if(j==round(j/20)*20)
  ###    {
  ###      print("writing")
  ###      #WriteMCMCResults(output_mcmc=outpop,main_struct,dir_name=nonmem_options$dir_name)
 ###       WriteMCMCIntermResults(j,results,accepted,Tann,candapp,main_struct,dir_name=nonmem_options$dir_name)
 ###     }
  ###  } 
    if((j==round(j/nonmem_options$jump)*nonmem_options$jump)&(j>1*nonmem_options$jump))
    {
      subseti<-SubsetVect(1:j,step_=nonmem_options$jump)
      newres<-rbind(prevres,candapp[subseti,])
      if(('mem'%in%names(nonmem_options))&((dim(newres)[1])-mem_step+1>0))#Corrected 23.04.2020!!!!
      {
        mem_step <- ceiling(nonmem_options$mem/nonmem_options$jump)
        newres<-newres[((dim(newres)[1])-mem_step+1):dim(newres)[1],]
        
      }
      #upd_cov<-FastCovAdv(arr=candapp[subseti,],prevcov=proposed_var,steps_=nonmem_options$step_init)
      #nonmem_options$step_init <- nonmem_options$step_init+1
      proposed_var<- t(chol(cov(newres)+diag(nonmem_options$eps_chol,dim(cov(newres)))))
    }
  }
  
  #out_$cand <- cand
  out_$results <-results
  out_$accepted <- accepted
  out_$Tann <- Tann
  out_$candapp <- candapp
  out_$nLL_indiv
  out_$residapp <- residapp
  out_$residapp_cumul <- residapp_cumul
  ###
  out_
}

GenerateMCMCSpecIndivStruct<-function(main_struct,nonmem_options)#proposed_var is a matriX!!!!
{
  compl_struct <- list()
  
  compl_struct$proposed_var<-list()
  
  
  Tann <- rep.int(x=1,times=nonmem_options$m+1)
  compl_struct$candapp <- list()
  compl_struct$Tann <- array(dim=c(nonmem_options$num_steps$indiv*nonmem_options$num_indiv_runs,main_struct$numindiv))
  compl_struct$results <- array(dim=c(nonmem_options$num_steps$indiv*nonmem_options$num_indiv_runs,main_struct$numindiv))
  compl_struct$accepted <- array(dim=c(nonmem_options$num_steps$indiv*nonmem_options$num_indiv_runs,main_struct$numindiv))
  compl_struct$prevres <- list()
  #compl_struct <- list()
  #compl_struct <- list()
  sum_i<-main_struct$numpoppar
  for(indi in 1:main_struct$numindiv)
  {
    ind_relevi <-(1+sum_i):(sum_i+main_struct$num_estim_indiv[indi])
    #compl_struct$proposed_var[[indi]]<- nonmem_options$proposed_var[ind_relevi,ind_relevi]
    compl_struct$prevres[[indi]] <- nonmem_options$prevres[,ind_relevi]
    compl_struct$proposed_var[[indi]]<- t(chol(cov(compl_struct$prevres[[indi]])+diag(nonmem_options$eps_chol,dim(cov(compl_struct$prevres[[indi]])))))
    compl_struct$candapp[[indi]] <- array(dim=c(nonmem_options$num_steps$indiv*nonmem_options$num_indiv_runs,main_struct$num_estim_indiv[indi]))
    sum_i <- sum_i+main_struct$num_estim_indiv[indi]
  }
  ##
  compl_struct
}
UpdateSAEMOptons<-function(nonmem_options)#proposed_var is a matriX!!!!
{
  nonmem_options$K1 <- 10
  nonmem_options$m <- 200
  nonmem_options$Ktot <- 100
  nonmem_options$gam<-vector(length = nonmem_options$Ktot)
  for(k in 1:nonmem_options$Ktot)
  {
    nonmem_options$gam[k]<-1/(1+max(0,k-nonmem_options$K1))
  }
  ###Output:
  nonmem_options
}
SAEMforResid<- function(xtot,ind_indiv,main_struct,nonmem_options)
{
  subseti<-SubsetVect(1:nonmem_options$m,step_=nonmem_options$jump)
  num_eit<-length(subseti)
  acurr_acc<-array(dim=c(nonmem_options$Ktot,length(main_struct$acurr)))
  acurr_acc[1,] <- main_struct$acurr
  all_goals <- vector(length = num_eit*(nonmem_options$Ktot-1))
  all_samp <- array(dim=c(num_eit*(nonmem_options$Ktot-1),main_struct$numvar-1))
  for(k in 2:nonmem_options$Ktot)
  {
    acurr_acc[k,] <- acurr_acc[k-1,]
    outpop <-MCMCDomesticHaarGeneral(xtot,main_struct,nonmem_options)
    if(sum(outpop$results==max(outpop$results))>1)
    {  
      xres<-outpop$candapp[outpop$results==max(outpop$results),][1,]
    }else
    {
      xres<-outpop$candapp[outpop$results==max(outpop$results),]
    }
    xtot<-InitParam(x=xres,xtot,stoch_obj0,main_struct)#outpop3$candapp[2000,]
    #param_<-ConstructParam(x=xtot$phiopt,ind_indiv=1,main_struct)
    goals<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
                        nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)#
    
    for(indsubs in 1: length(subseti))
    {
      acurr_acc[k,]<- acurr_acc[k,] + nonmem_options$gam[k]*(outpop$residapp[subseti[indsubs],]-acurr_acc[k,])
    }
    all_goals[(num_eit*(k-2)+1):(num_eit*(k-1))] <- outpop$results[subseti]
    all_samp[(num_eit*(k-2)+1):(num_eit*(k-1)),] <- outpop$candapp[subseti,]
    main_struct$acurr<- acurr_acc[k,]
    nonmem_options$prevres<-rbind(nonmem_options$prevres,outpop$candapp[subseti,])
  }
  out_ <- list()
  out_$acurr_acc <- acurr_acc
  out_$acurr <- acurr_acc[nonmem_options$Ktot,]
  out_$all_goals <- all_goals
  out_$all_samp <- all_samp
    ###
    out_
}


SAEM_SMStep<- function(gam_,outpopi,prev_pop,nonmem_options,relev_ind,num_par,xtot,main_struct,stoch_objp)#notreat_ind,
{
  curr_pop <- prev_pop
  burn<-nonmem_options$burn
  num_sim<- dim(outpopi[[relev_ind[1]]]$candapp)[1]#nonmem_options$m
  if(gam_==1)
  {  
    curr_pop$s1<-array(dim=c(max(relev_ind),num_par))
    curr_pop$s2 <- 0
    curr_pop$s5 <- 0
    curr_pop$s51 <- 0
  }else{
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
    curr_pop$s1[ind_indiv,] <- StochUpd(x=prev_pop$s1[ind_indiv,],m=outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],gam_)
    #curr_pop$s1[ind_indiv,] <- prev_pop$s1[ind_indiv,]+gam_*(apply(outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],2,mean)-prev_pop$s1[ind_indiv,])
    #ind_measi <- (!is.na(stoch_objp[[ind_indiv]]$statesm))
    #num_measi <- sum(ind_measi)
    #stoch_obj0 <-CreateStochProtocolDataPrimitive(main_struct,ind_indiv)
    
    for(ind_pat in (burn+1):num_sim)
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
   help_s5 <- help_s5/(num_sim-burn)
   help_s51 <- help_s51/(num_sim-burn)
  #approximation step:
  
  curr_pop$s5<- curr_pop$s5+gam_*(help_s5-curr_pop$s5)
  curr_pop$s51<- curr_pop$s51+gam_*(help_s51-curr_pop$s51)
  if(num_indiv>1)
  { 
    curr_pop$mu<-apply(curr_pop$s1[relev_ind,],2,mean)
    ind_indiv <- relev_ind[1]
    help1<-AddMatrVectRows(m=outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],x= - curr_pop$mu)^2
    for(ind_indiv in relev_ind[2:length( relev_ind)])
    {
      help1 <- help1 + AddMatrVectRows(m=outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],x= - curr_pop$mu)^2
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
SAEMforResidAdv<- function(gam_,outpopi,prev_pop,nonmem_options,relev_ind,num_par,xtot,main_struct,stoch_objp)#notreat_ind,
{
  curr_pop <- prev_pop
  burn<-nonmem_options$burn
  num_sim<- dim(outpopi[[relev_ind[1]]]$candapp)[1]#nonmem_options$m
  if(gam_==1)
  {  
    curr_pop$s1<-array(dim=c(max(relev_ind),num_par))
    curr_pop$s2 <- 0
    curr_pop$s5 <- 0
  }else{
    curr_pop$s5 <- prev_pop$s5
  }
  
  help_s5 <- array(dim=c(max(relev_ind),main_struct$YTYPES$num))
  
  
  res_simstochp<-list()
  num_steps_tot<-0
  num_steps_no_treat<-0
  num_steps_meas<-0
  num_indiv <- length(relev_ind)
  num_measi <- array(dim=c(max(relev_ind),main_struct$YTYPES$num))
  for(ind_indiv in relev_ind)
  {
    help_s5[ind_indiv,] <- rep.int(x=0, times = main_struct$YTYPES$num)
    curr_pop$s1[ind_indiv,] <- StochUpd(x=prev_pop$s1[ind_indiv,],m=outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],gam_)
    #curr_pop$s1[ind_indiv,] <- prev_pop$s1[ind_indiv,]+gam_*(apply(outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],2,mean)-prev_pop$s1[ind_indiv,])
    #ind_measi <- (!is.na(stoch_objp[[ind_indiv]]$statesm))
    #num_measi <- sum(ind_measi)
    #stoch_obj0 <-CreateStochProtocolDataPrimitive(main_struct,ind_indiv)
    
    for(ind_pat in (burn+1):num_sim)
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
      #goals<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod=opt_mod,label_opt=label_opt,
      #                    nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)#  
    }
    num_measi[ind_indiv,] <- as.numeric(main_struct$measurements[[ind_indiv]]$length_data)#temp_resid$numobs
    #num_steps_meas <- num_steps_meas+num_measi#population residual error is calculated only for untreated subjects!!! 
    
  }
  #average on chains:
  help_s5 <- help_s5/(num_sim-burn)
  
  #approximation step:
  
  curr_pop$s5<- curr_pop$s5+gam_*(help_s5-curr_pop$s5)
  if(num_indiv>1)
  { 
    curr_pop$mu<-apply(curr_pop$s1[relev_ind,],2,mean)
    ind_indiv <- relev_ind[1]
    help1<-AddMatrVectRows(m=outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],x= - curr_pop$mu)^2
    for(ind_indiv in relev_ind[2:length( relev_ind)])
    {
      help1 <- help1 + AddMatrVectRows(m=outpopi[[ind_indiv]]$candapp[(burn+1):num_sim,1:num_par],x= - curr_pop$mu)^2
    }
    curr_pop$s2 <- StochUpd(x=curr_pop$s2,m=help1,gam_)#mean on chains
    curr_pop$om<- (curr_pop$s2/num_indiv)^0.5
  }
  #if(sum(relev_ind%in%notreat_ind)>0)
  #{  
  curr_pop$sig_resid <- (curr_pop$s5/num_measi)^0.5#num_steps_meas
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
StochUpd <-function(x,m,gam_)#update by the given chain
{
  if(gam_<1)
  {
    out_ <- x+gam_*(apply(m,2,mean)-x)
  }else{
    out_ <- apply(m,2,mean)
  }
  
  return(out_)
}
AddMatrVectRows<-function(m,x)
{
  out_ <- m
  for(indi in 1:dim(m)[1])
  {
    out_[indi,]<-out_[indi,]+x
  }
  return(out_)
}
MetricsMultivar <- function(ve1,ve2,mat)
{
  deltamu<- ve1-ve2
  return((deltamu)%*%solve(mat)%*%t(t(deltamu)))
}