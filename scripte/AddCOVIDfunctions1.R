AlternativeshutdownsCollectPredictions <- function(predt,xtot,outpop,main_struct,stoch_obj0,ind_indiv,param_,schutd_sc,alpha_,ind_blockeff)#alpha_=0.05
{
  xtot$num_steps<-xtot$num_steps+predt
  stoch_obj0$num_steps <- xtot$num_steps
  xtot$stox <- vect2arr(x = rep.int(x=0,times = main_struct$numvar*stoch_obj0$num_steps),nrows = stoch_obj0$num_steps,ncols=main_struct$numvar)
  out_str <- list()
  out_str$altern_sc<-list()
  if(!is.null(schutd_sc))
  {
    num_alt_scenarios<-length(schutd_sc)
  }else{
    num_alt_scenarios<-1
  }
  num_sol<-sum(outpop$accepted)*num_alt_scenarios
  out_str$accepted<-outpop$accepted
  if(main_struct$label_treat_info==1)
  {  
    out_str$arr<-array(dim=c(xtot$num_steps,(main_struct$numvar+3),num_sol))
  }else{
    {  
      out_str$arr<-array(dim=c(xtot$num_steps,(main_struct$numvar+1),num_sol))
    }
  }
  ind_ <- 0
  if(main_struct$label_treat_info==1)
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
  if(main_struct$label_treat_info==1)
  {
    sim_res_best<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
  }else{
    sim_res_best<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
  }
  #sim_res_best<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
  sim_res_best<-NameOutput(sim_res_best,main_struct)
  xtot_alt_sc<-xtot
  
  #if(main_struct$names[main_struct$subs_popi_num[[ind_indiv]]])
  
  #sim_res_full_shutd<-SimCOVIDStochastic(xtot,main_struct,ind_indiv=1)
  #sim_res<-NameOutput(sim_res,main_struct)
  #fits_res$xtot
  num_accepted <- sum(outpop$accepted)
  if(num_alt_scenarios>0)
  {
    for(ind_sc in 1:num_alt_scenarios)
    {
      ind_ <- 0
      main_struct$init_mean[main_struct$names==paste('blockeff',ind_blockeff,sep='')]<-schutd_sc[ind_sc]
      for(indi in 1:dim(outpop$candapp)[1] )
      {
        if(outpop$accepted[indi]==1)
        {
          ind_<-ind_+1
          xtot<-InitParam(x=outpop$candapp[indi,],ind_indiv,xtot,stoch_obj0,main_struct)
          #xtot$num_steps<-xtot$num_steps+predt
          if(main_struct$label_treat_info==1)
          {
            sim_res<-SimCOVIDStochasticAdvTreat(xtot,main_struct,ind_indiv)
          }else{
            sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)#ind_indiv=1
          }
          #sim_res<-SimCOVIDStochastic(xtot,main_struct,ind_indiv)
          sim_res<-NameOutput(sim_res,main_struct)
          sim_res <- sim_res[,1:(main_struct$numvar+3)]
          sim_resp<-sim_res
          out_str$arr[,,ind_+(ind_sc-1)*num_accepted]<-sim_res
          
          sn<-sn+sim_res+((indi*sim_res-mn)^2)/(indi*(indi+1))
          mn<-mn+sim_res
          
          for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
          {
            if(ind_ytype=="Death")
            {
              sim_resy[[ind_ytype]]<-Deathreportscumul(S=sim_res[,],param_)
            }else{
              sim_resy[[ind_ytype]]<- CalculateOutput(sim_res,ind_ytype,param_)
            }
            out_str[[ind_ytype]]$sn <- out_str[[ind_ytype]]$sn + sim_resy[[ind_ytype]] + ((indi*sim_resy[[ind_ytype]]-out_str[[ind_ytype]]$mn)^2)/(indi*(indi+1))
            out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mn+sim_resy[[ind_ytype]]
          }
          sim_resyp <- sim_resy
        }else{
          sn<-sn+sim_resp+((indi*sim_resp-mn)^2)/(indi*(indi+1))
          mn<-mn+sim_resp
          for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
          {
            out_str[[ind_ytype]]$sn <- out_str[[ind_ytype]]$sn  + sim_resyp[[ind_ytype]] + ((indi*sim_resyp[[ind_ytype]]-out_str[[ind_ytype]]$mn)^2)/(indi*(indi+1))
            out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mn+sim_resyp[[ind_ytype]]
          }
          #out_str$arr[,,ind_]<-out_str$arr[,,ind_-1]#!!!!
        }
        print(c(ind_,indi))
      }
    }
  }
 
  out_str$mean <- mn/dim(outpop$candapp)[1]
  out_str$sd <- (sn/(dim(outpop$candapp)[1]+1))^0.5
  for(ind_ytype in main_struct$YTYPE[[ind_indiv]]$list)
  {
    out_str[[ind_ytype]]$sn <- (out_str[[ind_ytype]]$sn/(dim(outpop$candapp)[1]+1))^0.5
    out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mn/dim(outpop$candapp)[1]
    ###Ovewrite!!!
    
    arr3d0<-out_str$arr
    print('times of CalcCondifenceFun, ind_ytype = ')
    print(ind_ytype)
    print(', cumul:')
    ptm <- proc.time()
    aaa<-CalcConfidenceFunAdv(outpop=outpop,arr3d=arr3d0,fun_=sum,arr_names=colnames(sim_res),
                           ind_ytype=ind_ytype,alpha_,param_,label_cummul=1,ind_indiv,main_struct,num_alt_scenarios)
    print((proc.time() - ptm))
    print('times of CalcCondifenceFun, ind_ytype = ')
    print(ind_ytype)
    print(', daily:')
    ptm <- proc.time()
    arr3d1 <- arr3d0#DiffCalclArr(out_str$arr) !!!10.04.2020, complicated, redo!!!!
    bbb<-CalcConfidenceFunAdv(outpop=outpop,arr3d=arr3d1,fun_=sum,arr_names=colnames(sim_res),
                           ind_ytype=ind_ytype,alpha_,param_,label_cummul=0,ind_indiv,main_struct,num_alt_scenarios)
    print((proc.time() - ptm))
    out_str[[ind_ytype]]$mean <- aaa$mn#NameOutput(out_str$mean,main_struct)
    out_str[[ind_ytype]]$sd <-  aaa$sd#NameOutput(out_str$sd,main_struct)
    
    out_str[[ind_ytype]]$mn <- out_str[[ind_ytype]]$mean 
    
    out_str[[ind_ytype]]$PerDay$mean <- bbb$mn#NameOutput(out_str$mean,main_struct)
    out_str[[ind_ytype]]$PerDay$sd <-  bbb$sd#NameOutput(out_str$sd,main_struct)
    if(length(alpha_)==1)
    {  
      out_str[[ind_ytype]]$PerDay$lb <-  bbb$lb
      out_str[[ind_ytype]]$PerDay$ub <-  bbb$ub
      out_str[[ind_ytype]]$lb <-  aaa$lb
      out_str[[ind_ytype]]$ub <-  aaa$ub
    }else{
      out_str[[ind_ytype]]$PerDay$quantile <-  bbb$quantile
      out_str[[ind_ytype]]$quantile <-  aaa$quantile
    }
  }
  
  
  out_str
}
CalcConfidenceFunAdv<-function(outpop,arr3d,fun_,arr_names,ind_ytype,alpha_,param_,label_cummul,ind_indiv,main_struct,num_alt_scenarios)
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
  #help_arr<-array(dim=c(length(outpop$accepted),length(sub_arr)))
  help_arr <- array(dim=c(dim(arr3d)[1],num_alt_scenarios*length(outpop$accepted)))
  #help_arr[1]<-vecti[1]
  #ind_<-1
  ind_t <- 0
  for(ind_sc in 1:num_alt_scenarios)
  {
    ####
    for(indi in 1:dim(outpop$candapp)[1] )
    {
      if(outpop$accepted[indi]==1)
      {
        ind_t<-ind_t+1
      }
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
      help_arr[,indi+(ind_sc-1)*dim(outpop$candapp)[1]] <- vecti#[ind_]
      #if(outpop$accepted[indi]==1)
      #{
      #  print(vecti)
      #}
    }
  }
    
    for(ind_time in 1:dim(arr3d)[1])
    {
      help_vect <-help_arr[ind_time,]
      out_$mn[ind_time] <-mean(help_vect)
      out_$sd[ind_time] <-sd(help_vect)
      if(length(alpha_)==1)
      {  
        out_$lb[ind_time] <- sort(help_vect)[round((alpha_/2)*dim(help_arr)[2],digits = 0)]
        out_$ub[ind_time] <- sort(help_vect)[round((1-alpha_/2)*dim(help_arr)[2],digits = 0)]
      }else{
        for(ind_q in 1: length(alpha_))
        {
          out_$quantile[ind_time,ind_q ] <- sort(help_vect)[round((alpha_[ind_q])*dim(help_arr)[2],digits = 0)]
        }
      }
      #Older:
      #if(length(alpha_)==1)
      #{  
      #  out_$lb[ind_time] <- sort(help_vect)[round((alpha_/2)*dim(outpop$candapp)[1],digits = 0)]
      #  out_$ub[ind_time] <- sort(help_vect)[round((1-alpha_/2)*dim(outpop$candapp)[1],digits = 0)]
      #}else{
      #  for(ind_q in 1: length(alpha_))
     #   {
     #     out_$quantile[ind_time,ind_q ] <- sort(help_vect)[round((alpha_[ind_q])*dim(outpop$candapp)[1],digits = 0)]
      #  }
    #  }
    }
    ###
  
   
  
  
  
  out_
}
AdvCalcConfidenceFunNext<-function(outpop,arr3d,fun_,arr_names,ind_ytype,alpha_,param_,label_cummul,ind_indiv,main_struct)
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
  #help_arr<-array(dim=c(length(outpop$accepted),length(sub_arr)))
  help_arr <- array(dim=dim(arr3d)[c(1,3)])#array(dim=c(dim(arr3d)[1],length(outpop$accepted)))
  #help_arr[1]<-vecti[1]
  #ind_<-1
  ind_t <- 0
  ###
  num_samp<-dim(outpopi[[ind_indiv]]$candapp)[1]
  ind_accepted<-(1:num_samp)[outpopi[[ind_indiv]]$accepted==1]
  num_accepted<-DiffCalcl(ind_accepted)
  N_accepted <- sum(outpopi[[ind_indiv]]$accepted)
  if(max(ind_accepted)==num_samp)
  {
    num_accepted<-c(num_accepted,1)
  }else{
    num_accepted<-c(num_accepted,num_samp-max(ind_accepted)+1)
  }
  
  prob_accepted<-num_accepted/sum(num_accepted)
  #out_str[[paste('sol',ind_ytype,sep='')]] <- array(dim=dim(out_str$arr)[c(1,3)])
  for(ind_s in 1:N_accepted)
  {
    s_i<-arr3d[,,ind_s]
    colnames(s_i) <- arr_names
    help_arr[,ind_s]<- CalculateOutput(S=s_i,ind_ytype,param_)
    help_arr[,ind_s] <- SimpDailyFromCumulative(x = help_arr[,ind_s],label_cummul,ind_ytype,av_num = param_$DelayData)
    
  }
  #test_ <- sample(x=out_str[[paste('sol',ind_ytype,sep='')]][298,], size=100000, replace = TRUE, prob = num_accepted/sum(num_accepted))
  
  ###
  num_samp_all <- 100000
  for(ind_time in 1:dim(arr3d)[1])
  {
    help_vect0 <- sample(x=help_arr[ind_time,], size=num_samp_all, replace = TRUE, prob = num_accepted/sum(num_accepted))
    pert_ <- rnorm(num_samp_all,mean=0,  sd= sd_rand)
    help_vectp<- pmax(TransformSol(y=help_vect0,label_transf ,main_struct,ind_ytype,ind_ytypenum,label_cummul)+pert_,0)
    help_vect<- RevTransformSol(y=help_vectp,label_transf ,main_struct,ind_ytype,ind_ytypenum,label_cummul)
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
      }
    }
  }
  
  
  
  out_
}