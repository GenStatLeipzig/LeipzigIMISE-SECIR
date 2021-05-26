ProposalDistrMultiv<-function(prev_chain,omega_,ignore_,opt_)
{
  dchain<-dim(prev_chain)
  if(opt_<=9)
  {
    center_<-prev_chain[dchain[1],]
  }
  if(opt_==9)
  {
    if(dchain[1]>ignore_)
    {
      permut<-round((dchain[1]-ignore_)*runif(1),0)
      center_<-prev_chain[permut+dchain[1]-(dchain[1]-ignore_),]
    }else{
      center_<-prev_chain[dchain[1],]
    }
  }
  if(is.array(omega_))
  {
    out_<-center_+(omega_)%*%rnorm(dchain[2])
  }else
  {
    out_<-center_+omega_*rnorm(dchain[2])
  }
  
  
  #####
  out_
}
ProposalDistInit<-function(x,omega_)
{
  center_<-x
  
  #out_<-center_+omega_*rnorm(length(x))
  if(is.array(omega_))
  {
    out_<-center_+(omega_)%*%rnorm(length(x))
  }else
  {
    out_<-center_+omega_*rnorm(length(x))
  }
  #####
  out_
}
MCMCDomestic<-function(x,xtot,main_struct,nonmem_options)
{
  ##num_chains is currently 1
  #xx<-generate_pop_vector(xtot,main_struct1)
  cand <- array(dim=c(nonmem_options$m,length(x)))
  candapp <- array(dim=c(nonmem_options$m,length(x)))
  results <- vector(length = nonmem_options$m)
  accepted <- vector(length = nonmem_options$m)
  Tann <- rep.int(x=1,times=nonmem_options$m+1)
  print(length(Tann))
  cand[1,]<-x
  candapp[1,]<-x
  Tann[1] <- 1
  accepted[1]<- 1
  y<-func_pop_thrombo(x,xtot,main_struct)
  results[1] <- y
  
  
  for(j in 2:nonmem_options$m)
  {
    
    if(j==2)
    {
      x1<-ProposalDistInit(x,nonmem_options$proposed_var*Tann[j])
    }else{
      x1<-ProposalDistrMultiv(prev_chain=candapp[1:j-1,],omega_=nonmem_options$proposed_var*Tann[j],nonmem_options$ignore,opt_=j-10*floor(j/10))
    }
    
    cand[j,]<-x1
    y1<-func_pop_thrombo(x1,xtot,main_struct)
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
      }else{
        candapp[j,]<-candapp[j-1,]
        results[j]<-results[j-1]
        accepted[j]<-0
      }
    }
    y <- results[j]
    if(j>nonmem_options$temp_check)
    {
      print("annealing criteria:")
      print(mean(accepted[(j-nonmem_options$temp_check+1):j]))
      if((mean(accepted[(j-nonmem_options$temp_check+1):j])<nonmem_options$acceptl)&(Tann[j]>nonmem_options$minvar))
      {
        Tann[(j+1):nonmem_options$m]<-Tann[j]/1.1#1.2
      }
      if(mean(accepted[(j-nonmem_options$temp_check+1):j])>nonmem_options$acceptu)
      {
        Tann[(j+1):nonmem_options$m]<-Tann[j]*1.1#1.2
      }
    }
    print("results:")
    print(c(j,Tann[j+1],results[j],accepted[j],alpha_,mean(accepted[1:j])))
    if(nonmem_options$write == 1)
    {
      if(j==round(j/20)*20)
      {
        print("writing")
        #WriteMCMCResults(output_mcmc=outpop,main_struct,dir_name=nonmem_options$dir_name)
        WriteMCMCIntermResults(j,results,accepted,Tann,candapp,main_struct,dir_name=nonmem_options$dir_name)
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
  #out_$ <-
  ###
  out_
  
}


MCMCDomesticMulticyclic<-function(x,xtot,main_struct,nonmem_options)
{
  outglob<-list()
  outpopi <- MCMCDomestic(x,xtot,main_struct,nonmem_options)
  subseti<-SubsetVect(1:nonmem_options$m,step_=nonmem_options$jump)
  ni<-length(subseti)
  for(namei in names(outpopi))
  {
    if(is.vector(outpopi[[namei]]))
    {
      outglob[[namei]]<-vector(length = length(subseti)*nonmem_options$nglob)#subseti[length(subseti)]
      outglob[[namei]][1:ni]<-outpopi[[namei]][subseti]
    }
    if(is.array(outpopi[[namei]]))
    {
      outglob[[namei]]<-array(dim=c(length(subseti)*nonmem_options$nglob,dim(outpopi[[namei]])[2]))
      outglob[[namei]][1:ni,]<-outpopi[[namei]][subseti,]
    }
  }
  xx<-outpopi$candapp[nonmem_options$m,]
  nonmem_options$proposed_var <- nonmem_options$proposed_var*outpopi$Tann[nonmem_options$m]
  for(ind_g in 2:nonmem_options$nglob)
  {
    print("global index = ")
    print("global index = ")
    print("global index = ")
    print(ind_g)
    outpopi <- MCMCDomestic(xx,xtot,main_struct,nonmem_options)
    for(namei in names(outpopi))
    {
      if(is.vector(outpopi[[namei]]))
      {
        outglob[[namei]][1:ni+ni*(ind_g-1)]<-outpopi[[namei]][subseti]
      }
      if(is.array(outpopi[[namei]]))
      {
        outglob[[namei]][1:ni+ni*(ind_g-1),]<-outpopi[[namei]][subseti,]
      }
    }
    xx<-outpopi$candapp[nonmem_options$m,]
    nonmem_options$proposed_var <- nonmem_options$proposed_var*outpopi$Tann[nonmem_options$m]
  }
  ###output:
  outglob
} 
MCMCDomesticHaar<-function(x,xtot,main_struct,nonmem_options)
{
  ##num_chains is currently 1
  #xx<-generate_pop_vector(xtot,main_struct1)
  cand <- array(dim=c(nonmem_options$m,length(x)))
  candapp <- array(dim=c(nonmem_options$m,length(x)))
  results <- vector(length = nonmem_options$m)
  accepted <- vector(length = nonmem_options$m)
  Tann <- rep.int(x=1,times=nonmem_options$m+1)
  print(length(Tann))
  cand[1,]<-x
  candapp[1,]<-x
  Tann[1] <- 1
  accepted[1]<- 1
  y<-func_pop_thrombo(x,xtot,main_struct)
  results[1] <- y
  
  
  for(j in 2:nonmem_options$m)
  {
    
    if(j==2)
    {
      x1<-ProposalDistInit(x,nonmem_options$proposed_var*Tann[j])
    }else{
      x1<-ProposalDistrMultiv(prev_chain=candapp[1:j-1,],omega_=nonmem_options$proposed_var*Tann[j],nonmem_options$ignore,opt_=j-10*floor(j/10))
    }
    
    cand[j,]<-x1
    y1<-func_pop_thrombo(x1,xtot,main_struct)
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
      }else{
        candapp[j,]<-candapp[j-1,]
        results[j]<-results[j-1]
        accepted[j]<-0
      }
    }
    y <- results[j]
    if(j>nonmem_options$temp_check)
    {
      print("annealing criteria:")
      print(mean(accepted[(j-nonmem_options$temp_check+1):j]))
      if((mean(accepted[(j-nonmem_options$temp_check+1):j])<nonmem_options$acceptl)&(Tann[j]>nonmem_options$minvar))
      {
        Tann[(j+1):nonmem_options$m]<-Tann[j]/1.1#1.2
      }
      if(mean(accepted[(j-nonmem_options$temp_check+1):j])>nonmem_options$acceptu)
      {
        Tann[(j+1):nonmem_options$m]<-Tann[j]*1.1#1.2
      }
      
    }
    print("results:")
    print(c(j,Tann[j+1],results[j],accepted[j],alpha_,mean(accepted[1:j])))
    if(nonmem_options$write == 1)
    {
      if(j==round(j/20)*20)
      {
        print("writing")
        #WriteMCMCResults(output_mcmc=outpop,main_struct,dir_name=nonmem_options$dir_name)
        WriteMCMCIntermResults(j,results,accepted,Tann,candapp,main_struct,dir_name=nonmem_options$dir_name)
      }
    } 
    if((j==round(j/nonmem_options$jump)*nonmem_options$jump)&(j>10*nonmem_options$jump))
    {
      subseti<-SubsetVect(1:j,step_=nonmem_options$jump)
      newres<-rbind(nonmem_options$prevres,candapp[subseti,])
      #upd_cov<-FastCovAdv(arr=candapp[subseti,],prevcov=nonmem_options$proposed_var,steps_=nonmem_options$step_init)
      #nonmem_options$step_init <- nonmem_options$step_init+1
      nonmem_options$proposed_var<- t(chol(cov(newres)+diag(0.000000001,dim(cov(newres)))))
    }
  }
   
  ###output:
  out_ <- list()
  out_$cand <- cand
  out_$results <-results
  out_$accepted <- accepted
  out_$Tann <- Tann
  out_$candapp <- candapp
  #out_$ <-
  ###
  out_
  
}

SubsetVect <- function(x,step_)
{
  ind_ <- step_*(1:floor(length(x)/step_))
  ###
  x[ind_]
}
FastCov <- function(arr,prevcov)
{
  kk<-dim(arr)[1]
  uu<-kk-1
  av_arr<-apply(arr,2,mean)
  if(kk==2)
  {
    out_<-(arr[1,] %*%t(arr[1,])+arr[2,] %*%t(arr[2,])-2*av_arr%*%t(av_arr))#(1/(kk-1))* (kk+1)
  }
  else
  {
    av_arrpr<-apply(arr[1:(kk-1),],2,mean)
    out_<-((uu-1)/uu)*prevcov+av_arrpr%*%t(av_arrpr) + (1/uu)*(-(uu+1)*av_arr%*%t(av_arr)+arr[kk,] %*%t(arr[kk,]))
  }
  out_
}
FastCovAdv <- function(arr,prevcov,prevav,steps_)
{

  uu<-steps_-1
  av_arr<-apply(arr,2,mean)
  
    av_arrpr<-apply(arr[1:(steps_-1),],2,mean)
    out_<-((uu-1)/uu)*prevcov+av_arrpr%*%t(av_arrpr) + (1/uu)*(-(uu+1)*av_arr%*%t(av_arr)+arr[steps_,] %*%t(arr[steps_,]))
 
  out_
}
