
func_indiv_thrombo<-function(x,xtot,ind_indiv,main_struct)#activation function of the first compartment
{
  print("start")
  xtot$indiv[[ind_indiv]]<-x
  #goalinew<-GoalFuncIndiv(xtot=xtot,ind_indiv,main_struct=main_struct1,opt_mod= "new")
  print("x = ")
  print(x)
  goalinew<-GoalFuncIndiv(xtot=xtot,ind_indiv,main_struct,opt_mod= "new")
  print(-goalinew$tot)
  ###output:
  -goalinew$tot
}
generate_indiv_vector <- function(xtot,main_struct,ind_indiv)
{
  x<- xtot$indiv[[ind_indiv]]
  ###
  x
}
generate_pop_vector <- function(xtot,main_struct)
{
  x<-vector(length=main_struct$numpoppar+sum(main_struct$num_estim_indiv))
  x[1:main_struct$numpoppar] <- xtot$pop
  sum_i<-main_struct$numpoppar
  for(indi in 1:main_struct$numindiv)
  {
    if(main_struct$num_estim_indiv[indi]>0)
    {   
      x[(1+sum_i):(sum_i+main_struct$num_estim_indiv[indi])] <- xtot$indiv[[indi]]
      sum_i <- sum_i+main_struct$num_estim_indiv[indi]
    }
  }
  ###
  x
}
generate_names_vecot <- function(xtot,main_struct)
{
  x<-vector(length=main_struct$numpoppar+sum(main_struct$num_estim_indiv))
  x[1:main_struct$numpoppar] <- main_struct$names[main_struct1$subs_pop_mu_opt==1]
  sum_i<-main_struct$numpoppar
  for(indi in 1:main_struct$numindiv)
  {
    x[(1+sum_i):(sum_i+main_struct$num_estim_indiv[indi])] <- paste(main_struct$names[main_struct$subs_indivi[indi,]==1],indi,sep='_')
    sum_i <- sum_i+main_struct$num_estim_indiv[indi]
  }
  ###
  x
}
generate_listvect<-function(x,xtot,main_struct)#activation function of the first compartment
{
  
  xtot$pop<-x[1:main_struct$numpoppar]
  sum_i<-main_struct$numpoppar
  for(indi in 1:main_struct$numindiv)
  {
    if(main_struct$num_estim_indiv[indi]>0)
    {  
      xtot$indiv[[indi]]<-x[(1+sum_i):(sum_i+main_struct$num_estim_indiv[indi])] 
      sum_i <- sum_i+main_struct$num_estim_indiv[indi]
    }  
  }
  xtot
}
update_listindivvect<-function(x,xtot,main_struct,ind_indiv)#activation function of the first compartment
{
  
  #sum_i<-main_struct$numpoppar
  
    xtot$indiv[[ind_indiv]]<-x#[(1+sum_i):(sum_i+main_struct$num_estim_indiv[indi])] 
    #sum_i <- sum_i+main_struct$num_estim_indiv[indi]
  
  xtot
}
func_pop_thrombo<-function(x,xtot,main_struct)#activation function of the first compartment
{
  
  xtot<-generate_listvect(x,xtot,main_struct)
  #print("start")
  
  #goalinew<-GoalFuncIndiv(xtot=xtot,ind_indiv,main_struct=main_struct1,opt_mod= "new")
  #print("x = ")
  #print(x)
  tot_goalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt="pop",
                          nLL_indiv=rep.int(x=0,times=main_struct$numindiv),ind_indiv1=ind_indiv)
  print("goalf:")
  print(-tot_goalf$nLL_tot)
  ###output:
  -tot_goalf$nLL_tot
}
func_general_thrombo<-function(x,xtot,ind_indiv,main_struct,optim_options,label_opt)#activation function of the first compartment
{
  
  xtot<-generate_listvect(x,xtot,main_struct)
  #print("start")
  
  #goalinew<-GoalFuncIndiv(xtot=xtot,ind_indiv,main_struct=main_struct1,opt_mod= "new")
  #print("x = ")
  #print(x)
  tot_goalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct,opt_mod= "new",label_opt=label_opt,#="pop",
                          nLL_indiv=optim_options$nLL_indiv_min,ind_indiv1=ind_indiv)
  print("goalf:")
  print(-tot_goalf$nLL_tot)
  ###output:
  -tot_goalf$nLL_tot
}
WriteMCMCResults<-function(output_mcmc,main_struct,dir_name)
{
  resultspath<-paste(main_struct$resultspath,dir_name,'/',sep='')
  if(!dir.exists(resultspath))
  {
    dir.create(resultspath)
  }
  running_inf<-cbind(outpop1$results,outpop1$accepted,outpop1$Tann[1:nonmem_options$m])
  colnames(running_inf)<-c('results','accepted','om_adjust')
  write.table(x=output_mcmc$candapp,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(resultspath,"accepted_phi.csv",sep=''))
  write.table(x=running_inf,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(resultspath,"running_info.csv",sep=''))
}
WriteMCMCIntermResults<-function(j,results,accepted,Tann,candapp,main_struct,dir_name)
{
  resultspath<-paste(main_struct$resultspath,dir_name,'/',sep='')
  if(!dir.exists(resultspath))
  {
    dir.create(resultspath)
  }
  running_inf<-cbind(results[1:j],accepted[1:j],Tann[1:j])
  colnames(running_inf)<-c('results','accepted','om_adjust')
  write.table(x=candapp[1:j,],sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(resultspath,"accepted_phi.csv",sep=''))
  write.table(x=running_inf,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(resultspath,"running_info.csv",sep=''))
}
UpdatePrevResFromResults<-function(main_struct, res_dir)#,relev_row
{
  allfiles <-list.files(res_dir)#main_struct$resultspath
  # for(ind_f in allfiles)
  # {
  #   file.copy(paste(old_dir,ind_f,sep =""), paste(res_dir,ind_f,sep =""))
  # }
  
  indiv_data<-read.csv(file=paste(res_dir,'indivpar.csv',sep=""), sep = ";",dec = main_struct$dec)
  resid_data<-read.csv(file=paste(res_dir,'resid.csv',sep=""), sep = ";",dec = main_struct$dec)
  
  error_data <- read.csv(file=paste(res_dir,'error.csv',sep=""), sep = ";",dec = main_struct$dec)
  num_rows<-sum(!is.na(error_data))
  num_columns<- main_struct$numpoppar+sum(main_struct$num_estim_indiv)
  out_ <- list()
  out_$phi <-array( dim = c(num_rows,num_columns) )
  out_$psi <-array( dim = c(num_rows,num_columns) )
  colnames(out_$phi) <- rep.int(x='a',times=num_columns )
  colnames(out_$psi) <- rep.int(x='a',times=num_columns )
  ###############Read results until prescribed row:
  sum_i<-main_struct$numpoppar
  for(ind_indiv in 1: main_struct$numindiv)
  {
    tempori <- read.csv(file=paste(res_dir,'maini',ind_indiv,'.csv',sep=""), sep = ";",dec = main_struct$dec)
    temporresidi <- read.csv(file=paste(res_dir,'residi',ind_indiv,'.csv',sep=""), sep = ";",dec = main_struct$dec)
    for(relev_row in 1:num_rows)
    {
      colnames(out_$psi)[1:main_struct$numpoppar]<-main_struct1$names[main_struct1$subs_pop_mu_opt==1]
      colnames(out_$phi)[1:main_struct$numpoppar]<-main_struct1$names[main_struct1$subs_pop_mu_opt==1]
      out_$psi[relev_row,1:main_struct$numpoppar]<-unlist(tempori[relev_row,main_struct1$subs_pop_mu_opt==1])
      
      #for(indi in 1:main_struct$numindiv)
      #{
      if(main_struct$num_estim_indiv[ind_indiv]>0)
      {   
        colnames(out_$psi)[(1+sum_i):(sum_i+main_struct$num_estim_indiv[ind_indiv])]<-main_struct1$names[main_struct1$subs_indivi_num[[ind_indiv]]]
        colnames(out_$phi)[(1+sum_i):(sum_i+main_struct$num_estim_indiv[ind_indiv])]<-main_struct1$names[main_struct1$subs_indivi_num[[ind_indiv]]]
        out_$psi[relev_row,(1+sum_i):(sum_i+main_struct$num_estim_indiv[ind_indiv])] <- unlist(tempori[relev_row,main_struct1$subs_indivi_num[[ind_indiv]]])
        
      }
      #}
      #indiv_data[ind_indiv,2:dim(indiv_data)[2]] <- tempori[relev_row,]
      #resid_data[ind_indiv,2:dim(resid_data)[2]] <- temporresidi[relev_row,]
    }
    sum_i <- sum_i+main_struct$num_estim_indiv[ind_indiv]
    
  }
  for(col_ind in 1:num_columns)
  {
    for(row_ind in 1:num_rows)
    {
      out_$phi[row_ind,col_ind]<-RevTransformParameters(x=out_$psi[row_ind,col_ind],main_struct,name=colnames(out_$phi)[col_ind])
    }
  }
  
  
  #write.table(x = indiv_data, file = paste(res_dir,'indivpar.csv',sep=""), sep = ";" ,  dec = main_struct$dec, row.names=F)
  #write.table(x = resid_data, file = paste(res_dir,'resid.csv',sep=""), sep = ";" ,  dec = main_struct$dec, row.names=F)
  #main_data<-read.csv(file=paste(res_dir,'main.csv',sep=""), sep = ";")#!!!!!!res_dir
  ##main_data$init_mean[2:(dim(indiv_data)[2]-1)] <- as.numeric(indiv_data[1,2:dim(indiv_data)[2]])#Any is good
  #main_data$init_mean <- as.numeric(indiv_data[1,2:dim(indiv_data)[2]])#Any is good
  #main_data$init_mean[main_data$names%in%c("xleuko0","xloc0","xanc0","xplc0","xtpo0")]<-1
  #write.table(x =main_data, file = paste(res_dir,'main.csv',sep=""), sep = ";" ,  dec = main_struct$dec, row.names=F)
  ###output:
  out_#res_dir
} 