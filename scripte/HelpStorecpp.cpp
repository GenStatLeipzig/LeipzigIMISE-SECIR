if(indi>0)
{
  treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff1, param_.Lock1,param_.delaylock);
  if(num_locks>1)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff2, param_.Lock2,param_.delaylock);
  }
  if(num_locks>2)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff3, param_.Lock3,param_.delaylock);
  }
  if(num_locks>3)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff4, param_.Lock4,param_.delaylock);
  }
  if(num_locks>4)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff5, param_.Lock5,param_.delaylock);
  }
  if(num_locks>5)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff6, param_.Lock6,param_.delaylock);
  }
  if(num_locks>6)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.blockeff7, param_.Lock7,param_.delaylock);
  }
  treat_arr(indi,1) = treat_arr(indi,1)*DelayEff(indi,1,param_.unlockeff1,param_.Unlock1,param_.delayunlock);
  
  if(num_unlocks>1)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff2, param_.Unlock2,param_.delayunlock);
  }
  if(num_unlocks>2)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff3, param_.Unlock3,param_.delayunlock);
  }
  if(num_unlocks>3)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff4, param_.Unlock4,param_.delayunlock);
  }
  if(num_unlocks>4)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff5, param_.Unlock5,param_.delayunlock);
  }
  if(num_unlocks>5)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff6, param_.Unlock6,param_.delayunlock);
  }
  if(num_unlocks>6)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff7, param_.Unlock7,param_.delayunlock);
  }
  if(num_unlocks>7)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff8, param_.Unlock8,param_.delayunlock);
  }
  if(num_unlocks>8)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1, param_.unlockeff9, param_.Unlock9,param_.delayunlock);
  }
  /*for(int ind_locks= 1; ind_locks <num_locks; ind_locks++)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,1,par2=param_[[paste('blockeff',ind_locks,sep='')]],
                                                  t1=param_[[paste('Lock',ind_locks,sep='')]],
                                                  param_,delaylock);
  }
  for(int ind_unlocks= 1; ind_unlocks <num_locks; ind_unlocks++)
  {
    treat_arr(indi,1)= treat_arr(indi,1)*DelayEff(indi,par1=1,par2=param_[[paste('unlockeff',ind_unlocks,sep='')]],
                                                  t1= param_[[paste('Unlock',ind_unlocks,sep='')]],
                                                  delaylock=param_$delayunlock);
  }*/
    
    
    if(param_.parsim_treat==0){
      treat_arr(indi,2) = treat_arr(indi,2)*DelayEff(indi,1,param_.blockeff_r2,param_.Lock2,
                                                     param_.delaylock);
    }
}