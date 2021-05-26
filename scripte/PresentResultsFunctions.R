PlotCovidFits <- function(sim_res,main_struct,ind_indiv,param_,label_save,log_option,date_option,label_cummul)
{
  #if(main_struct0$numindiv==1)
  #{  
  #  delaydatai<-main_struct$covariates[ind_indiv,'DelayData']
 # }else{
  #  delaydatai<-main_struct$covariates$DelayData[ind_indiv]
  #}
  #delaydatai <- param_$DelayData
  #ind_indiv<-1
  if(label_save==1)
  {
    if(log_option==1)
    {
      add_name<-'log_scale'
    }else{
      add_name<-'normal_scale'
    }
    if(label_cummul==1)
    {
      add_name<-paste(add_name,"Cumulative",sep="")
    }
    else{
      add_name<-paste(add_name,"Daily",sep="")
    }
    add_name<-paste(add_name,main_struct$covariates$Land[ind_indiv])
    #tiff(filename = paste(main_struct$resultspath,"PlotFits",add_name,".tif",sep=''),
    #     width = 2250, height = 1300, units = "px", pointsize = 12, res= 300)#px 2250 2625
    png(filename = paste(main_struct$resultspath,"PlotFits",add_name,".png",sep=''),
        width = 2250/2, height = 1300/2, units = "px", pointsize = 12, res= 300/2) 
        #width = 2250, height = 1300, units = "px", pointsize = 12, res= 300)#px 2250 2625
  }#bmp
  set_plot_dimensions(5,5)
  num_plots_in_row<- main_struct$YTYPE[[ind_indiv]]$num#2
  par(mfrow=c(1,main_struct$YTYPE[[ind_indiv]]$num+1)) #+1
  #all_times<-unique(main_struct$measurements[[ind_indiv]]$data[1,])
  if(date_option==0)
  {  
    all_times<-unique(main_struct$measurements[[ind_indiv]]$data[1,])
  }else{
    all_times<-unique(main_struct$measurements[[ind_indiv]]$date)
  }
  for (ind_ytype0 in 1:main_struct$YTYPE[[ind_indiv]]$num)
  {#main_struct$YTYPES$num#main_struct$YTYPE[[ind_indiv]]$num
    ind_ytype=main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0];
    ind_ytypenum=main_struct$YTYPE[[ind_indiv]]$listnum[ind_ytype0];
    print(ind_ytype)
    #print(c(main_struct$YTYPE[[ind_indiv]]$size[[ind_ytypenum]],))
    currentylab<-Deriveylabels(ind_ytype,label_cummul)
    if(main_struct$YTYPE[[ind_indiv]]$size[[ind_ytypenum]]>0){
      if(date_option==0)
      {
        xx <- main_struct$measurements[[ind_indiv]]$data[1,main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
        xx_comp<-xx
      }else{
        xx <- main_struct$measurements[[ind_indiv]]$date[main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
        xx_comp<-CompleteDates(all_times)
      }
      
      #xx <- main_struct$measurements[[ind_indiv]]$data[1,main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
      yy <- main_struct$measurements[[ind_indiv]]$data[2,main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
      relev_indi<-(all_times%in%xx)
      zzdense<-CalculateOutput(S=sim_res[,],ind_ytype,param_)#sim_res[relev_indi,]
      
      zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul,ind_ytype,av_num = param_$DelayData)
      if((label_cummul==0)&(ind_ytype!="Critical"))
      {
        yy <- main_struct$Daily$measurements[[ind_indiv]]$data[2,main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
       
        
      }
      #if(ind_ytype=='Total')
      #{  
      #    zzdense<-apply(sim_res[relev_indi,c('Cc','Dc', 'Rc1')],1,sum)#+param_$procentmeas*sim_res[relev_indi,'ISc']
      # } 
      # if(ind_ytype=='Death')
      # {  
      #   zzdense<- sim_res[relev_indi,'Dc']
      # }
      # if(ind_ytype=='Recovery1')
      # {  
      #   zzdense<- sim_res[relev_indi,'Rc1']
      # }
      # if(ind_ytype=='Critical')
      # {  
      #   zzdense<- sim_res[relev_indi,'Cc']
      # }
      ylims <-c(min(min(zzdense,na.rm =T),min(yy,na.rm =T)),max(max(zzdense,na.rm =T),max(yy,na.rm =T)))
      xlims<-c(min(all_times),max(all_times))#xlims <- c(min(xx),max(xx))#_comp !!!!
      
      #if(ind_ytype0==1)
      #{  
      #  xlims1<-c(xlims[1],xlims[2]+1.5*(xlims[2]-xlims[1]))
      #  ylims1<-c(ylims[1],ylims[2]+0.3*(ylims[2]-ylims[1]))
      #}else{
      xlims1 <- xlims
      ylims1 <- ylims
      #}
      
      if(log_option==1)
      {  
        ylims[1]<-max(1,ylims[1])
        plot(xx,yy,col = "blue",pch = 1,ylim=ylims,xlim=xlims1,xlab="Time, days",cex=2,lwd=2,log= "y",  
             cex.lab=1.5, cex.axis=1.2,mex=1.5,ylab=currentylab)# , ylab =currentylab mex=1.5mex=2pch = ind_ytype0
      }else{
        plot(xx,yy,col = "blue",pch = 1,ylim=ylims1,xlim=xlims1,xlab="Time, days", ylab =currentylab,cex=2,lwd=2,
             cex.lab=1.5, cex.axis=1.2,mex=1.5)#,log= "y")#pch = ind_ytype0
      }
      #mtext(currentylab,side = 2, padj = 1.5)#2
      #points(xx,zzdense,col = rgb(1,0.3,0),pch = 1,cex=2,lwd=2) #, padj = 1.5
      ##not necessary:
     # if((label_cummul==0)&(ind_ytype!="Critical"))
     # {
      #  lines(xx,AverSim(yy,av_num= param_$DelayData),col='red',lwd=3)#c(0,)
     # }
      if(ind_ytype=="Critical")
      {  
          lines(xx_comp,zzdense,col = rgb(0,1,0),lwd=3)#rgb(1,0.3,0)
      }
      if(ind_ytype=="Death")
      {
        if(label_cummul==1)
        {  
          #zzdensem<-ComplexCumulDeathreportsTransform(S=sim_res[,],param_)
          zzdensem<-Deathreportscumul(S=sim_res[,],param_)
          zzreported <- Deathreportsonlycumul(S=sim_res[,],param_)
          
        }else{  
          zzdensem<-DeathreportsTransform(S=sim_res[,],param_)
          zzreported <- Deathreportsonly(S=sim_res[,],param_)
           
        }
        #lines(xx_comp,AverSim(zzdensem,av_num= param_$DelayData),col = 'magenta',lwd=3)#rgb(1,0.3,0)
        lines(xx_comp,zzdensem,col = rgb(0,1,0),lwd=3)#rgb(1,0.3,0) # 'magenta
        ##not necessary:
        #lines(xx_comp,zzreported,col = 'black',lwd=3)#rgb(1,0.3,0)
        
        #res_[[ind_ytypenum]]<-DeathreportsTransform(S,param_)
      }
      if(ind_ytype=="Total")
      {
        if(label_cummul==1)
        {  
          #zzdensem<-ComplexCumulDeathreportsTransform(S=sim_res[,],param_)
          zzdensem<-NewCasesreportscumul(S=sim_res[,],param_)
          #zzreported <- Deathreportsonlycumul(S=sim_res[,],param_)
          
        }else{  
          zzdensem<-NewCasesreportsTransformDaily(S=sim_res[,],param_)
         # zzreported <- Deathreportsonly(S=sim_res[,],param_)
          
        }
        #lines(xx_comp,AverSim(zzdensem,av_num= param_$DelayData),col = 'magenta',lwd=3)#rgb(1,0.3,0)
        lines(xx_comp,zzdensem,col  = rgb(0,1,0),lwd=3)#rgb(1,0.3,0)# magenta
        ##not necessary:
        #lines(xx_comp,zzreported,col = 'black',lwd=3)#rgb(1,0.3,0)
        
        #res_[[ind_ytypenum]]<-DeathreportsTransform(S,param_)
      }
      if((main_struct$label_treat_info==1)||(main_struct$label_treat_info==2))
      {
        if(main_struct$label_treat_info==2)
        {
          main_struct$num_locks <- param_$num_locks 
          main_struct$num_unlocks <- param_$num_unlocks 
        }
        for(ind_locks in 1:main_struct$num_locks)
        {
          restricti <- min(all_times)+param_[[paste('Lock',ind_locks,sep='')]]
          if(ind_locks==2)
          {
            lines(c(restricti,restricti),ylims,col='black',lwd=1,lty=3)#lwd =3
          }else{
            lines(c(restricti,restricti),ylims,col='black',lwd=1,lty=1)#lwd =3
          }
          #if(!main_struct$covariates[[paste('DateLock',ind_locks,sep='')]][ind_indiv]==restricti)
          #{
           # print("Error tretament einformation!!!")
          #}
        }
        for(ind_unlocks in 1:main_struct$num_unlocks)
        {
          restricti <- min(all_times)+param_[[paste('Unlock',ind_unlocks,sep='')]]
          
            lines(c(restricti,restricti),ylims,col='red',lwd=1,lty=1)#lwd =3
          
          
         # if(!main_struct$covariates[[paste('Unlock',ind_unlocks,sep='')]][ind_indiv]==restricti)
          #{
         #   print("Error tretament einformation!!!")
         # }
          
          
        }
        
        par(new=TRUE)
        #lines(xx_comp,sim_res[,'vr1']*ylims1[2]/max(max(sim_res[,'vr1']),max(sim_res[,'vr2'])),col='cyan',lwd=3,lty=3)
        #lines(xx_comp,sim_res[,'vr2']*ylims1[2]/max(max(sim_res[,'vr1']),max(sim_res[,'vr2'])),col='brown',lwd=3,lty=3)
        plot(xx_comp,sim_res[,'vr1']/max(max(sim_res[,'vr1']),max(sim_res[,'vr2'])),col='cyan',lwd=3,lty=3, ylim=c(0,1),
             axes = FALSE,xlab = '',ylab = '',type='l')#pch=NA_integer_
        lines(xx_comp,sim_res[,'vr2']/max(max(sim_res[,'vr1']),max(sim_res[,'vr2'])),col='brown',lwd=3,lty=3)
        if(ind_ytype=='Critical')
        {
          lines(xx_comp,sim_res[,'pcrit']/max(sim_res[,'pcrit']),col='yellow',lwd=3,lty=3)
        }
        if(ind_ytype=="Death")
        {  
          lines(xx_comp,sim_res[,'pdeath']/max(sim_res[,'pdeath']),col='violet',lwd=3,lty=3)
        }
        axis(side=4,ylim=c(0.01,1),at = c(0,0.5,1))#,at = pretty(range(0:1),2)
        #if(ind_ytype=="Total")
       # {  
          mtext('Relative infection rates',side = 4, padj = 0.7,adj=0.2)#2
       # }
        #treat_arr[indi,'r2']<- treat_arr[indi,'r2']*DelayEff(t0=indi,par1=1,par2=param_$blockeff_r2,t1=param_$Lock2,
         #                                                    delaylock=param_$delaylock)
         
      }
      else
      {  
        restrict1 <- min(all_times)+param_$Lock1
        restrict2 <- min(all_times)+param_$Lock2
        if(main_struct$label_treat_info==0)
        {  
          restrict3 <- min(all_times)+param_$Unlock
        }else{
          restrict3 <- min(all_times)+param_$Unlock1
        }
        if(min(xx)<restrict1)
        {
          lines(c(restrict1,restrict1),ylims,col='black',lwd=1,lty=3)#lwd =3
          
        }
        if(min(xx)<restrict2)
        {
          
          lines(c(restrict2,restrict2),ylims,col='black',lwd=1,lty=1)#lwd =3
        }
        if(min(xx)<restrict3)
        {
          
          lines(c(restrict3,restrict3),ylims,col='red',lwd=1,lty=1)#lwd =3
        }
        
        # if(ind_ytype0==1)
        # {
        #   AddLegend(label_grades,label_barplots,xlims,ylims1)
        # }
        }
    }
  }
  ###legends:
  #plot.new()
  plot(0,type='n',axes=FALSE,ann=FALSE)
  #par(mar=rep(0,4),mai=0.3*c(0,0,0,0),pin=0.5*c(1,2),fin=1*c(1,2))#,pin=c(5,5) (1,2),mar = c(bottom=0, left=0, top=0, right=0),
  
  legend("topleft",#xlims[2]-0.5*(xlims[2]-xlims[1]),ylims1[1],#"topleft",##"topright",inset=0,
         box.lty = 0,inset=0,# cex=1,
         legend=c("Data","Fitted","Shut-down events","Release events", 
                  'Relative dynamic b1',  'Relative dynamic b2','Relative pcrit', 'Relative pdeath'),
         col = c('blue','green','black','red','cyan','brown','yellow','violet'),lty = c(0,1,1,1,3,3,3,3),
         lwd=c(3,3,1,1,3,3,3,3),cex = 0.8,#c(1,1,1,1,1,1,1,1)#rep.int(x=0.6,times=8)
         pch = c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_,NA_integer_,NA_integer_,NA_integer_))
  if(label_save==1)
  {
    dev.off()
  }
}
set_plot_dimensions <- function(width_choice, height_choice) {
  options(repr.plot.width=width_choice, repr.plot.height=height_choice)
}
####
PlotCovidPredictions <- function(coll_pred,main_struct,ind_indiv,param_,label_save,log_option,date_option,
                                 label_all,label_cummul,label_altern,predt0)#,add_name
{
  if(main_struct0$numindiv==1)
  {  
    delaydatai<-main_struct$covariates[ind_indiv,'DelayData']
  }else{
    delaydatai<-main_struct$covariates$DelayData[ind_indiv]
  }
  if(label_cummul==0){
  add_name1<-paste('Daily')
  }else{
  add_name1<-paste('Cumulative')
  }
  if(label_save==1)
  {
    if(log_option==1)
    {
      add_name<-paste('For',predt0,'Days','log_scale',sep='')
    }else{
      add_name<-paste('For',predt0,'Days','normal_scale',sep='')
    }
    
    if(label_cummul==1)
    {
      add_name<-paste(add_name,"Cumulative",sep="")
    }
    else{
      add_name<-paste(add_name,"Daily",sep="")
    }
   # tiff(filename = paste(main_struct$resultspath,"PlotPredictions",add_name,".tif",sep=''),
   #      width = 2250/2, height = 1300/2, units = "px", pointsize = 12, res= 300/2)#px 2250 2625
    png(filename = paste(main_struct$resultspath,"PlotPredictions",add_name,".png",sep=''),
        width = 2250/2, height = 1300/2, units = "px", pointsize = 12, res= 300/2)#px 2250 2625 res= 300
  }#bmp
  ####
  #ind_indiv<-1
  num_plots_in_row<- main_struct$YTYPE[[ind_indiv]]$num#2
  par(mfrow=c(1,main_struct$YTYPE[[ind_indiv]]$num)) #+1
  if(date_option==0)
  {  
    all_times<-unique(main_struct$measurements[[ind_indiv]]$data[1,])
  }else{
    all_times<-unique(main_struct$measurements[[ind_indiv]]$date)
  }
  
  all_times<-c(min(all_times),min(all_times+1):max(all_times))#!!!!
  
  for (ind_ytype0 in 1:main_struct$YTYPE[[ind_indiv]]$num)
  {#main_struct$YTYPES$num#main_struct$YTYPE[[ind_indiv]]$num
    ind_ytype=main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0];
    ind_ytypenum=main_struct$YTYPE[[ind_indiv]]$listnum[ind_ytype0];
    print(ind_ytype)
    #print(c(main_struct$YTYPE[[ind_indiv]]$size[[ind_ytypenum]],))
    currentylab<-Deriveylabels(ind_ytype,label_cummul)
    if(length(coll_pred$altern_sc)>0)
    {
      help_v <- SimpDailyFromCumulative(x = coll_pred$altern_sc[[ind_ytype]]$best[,1],
                                   label_cummul,ind_ytype,av_num = param_$DelayData)
      for(ind_sc in 2: dim(coll_pred$altern_sc[[ind_ytype]]$best)[2])
      {
        help_vi <- SimpDailyFromCumulative(x = coll_pred$altern_sc[[ind_ytype]]$best[,ind_sc],
                                      label_cummul,ind_ytype,av_num = param_$DelayData)
        help_v<- cbind(help_v,help_vi)
      }
    }
    if(main_struct$YTYPE[[ind_indiv]]$size[[ind_ytypenum]]>0){
      if(date_option==0)
      {
        xx <- main_struct$measurements[[ind_indiv]]$data[1,main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
       # xx0<-c(xx,(max(xx)+1):(max(xx)+predt0))
        xlims <- c(min(xx0),max(xx0))
      }else{
        xx <- main_struct$measurements[[ind_indiv]]$date[main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
        xx0<-c(all_times,(max(all_times)+1):(max(all_times)+predt0))#xx0<-c(xx,(max(xx)+1):(max(xx)+predt0))
        xlims <- c(min(xx0),max(xx0))
      }
      yy <- main_struct$measurements[[ind_indiv]]$data[2,main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
      #relev_indi<-c((all_times%in%xx),rep.int(x=TRUE,times = predt0))
      #yy <- SimpDailyFromCumulative(x = yy,label_cummul,ind_ytype,av_num = delaydatai)##!!!!
      if((label_cummul==0)&(ind_ytype!="Critical"))
      {
        yy <- main_struct$Daily$measurements[[ind_indiv]]$data[2,main_struct$measurements[[ind_indiv]]$datatypes==main_struct$YTYPE[[ind_indiv]]$list[ind_ytype0]]
        
        
      }
      
      ylims <-c(min(yy),max(yy))
     # ylims <- c(min(ylims[1],apply(help_v,2,min)),max(ylims[1],apply(help_v,2,max)))
      #if(label_all==1)
       
      for(ind_sol in 1:dim(coll_pred$arr)[3])
      {
        sim_res<-coll_pred$arr[,,ind_sol]
        sim_res<-NameOutput(sim_res,main_struct)
        #colnames(sim_res)<- main_struct$colnames_sim_res#colnames(coll_pred$mean)#NameOutput(sim_res,main_struct)
        
        zzdense<-CalculateOutput(S=sim_res,ind_ytype,param_)#sim_res[relev_indi,]
        zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul,ind_ytype,av_num = delaydatai)
        
        
        if((min(zzdense)>=0)&(label_all==1))
        {  
          ylims <-c(min(min(zzdense),ylims[1]),max(max(zzdense),ylims[2]))
        }
      }
      #temp_pred<-cbind(xx0,coll_pred[[ind_ytype]]$mn[relev_indi],coll_pred[[ind_ytype]]$sn[relev_indi])
      temp_pred<- list()
      temp_pred$Date<-as.Date(xx0)
      temp_pred$lb <- min(coll_pred[[ind_ytype]]$quantile)#coll_pred[[ind_ytype]]$lb#temp_pred$mean_predict <- coll_pred[[ind_ytype]]$mn#[relev_indi]
      temp_pred$ub <- max(coll_pred[[ind_ytype]]$quantile)#coll_pred[[ind_ytype]]$ub#temp_pred$std_predict <- coll_pred[[ind_ytype]]$sn#[relev_indi]
      temp_pred$mean_predict <- coll_pred[[ind_ytype]]$mean
      if(label_cummul==0)
      {
        temp_pred$lb <- min(coll_pred[[ind_ytype]]$PerDay$quantile)#coll_pred[[ind_ytype]]$PerDay$lb
        temp_pred$ub <- max(coll_pred[[ind_ytype]]$PerDay$quantile)#coll_pred[[ind_ytype]]$PerDay$ub
        temp_pred$mean_predict <- coll_pred[[ind_ytype]]$PerDay$mean
      }
      #temp_pred$lb <- DailyFromCumulative(x = temp_pred$lb,label_cummul,ind_ytype,av_num = delaydatai)
      #temp_pred$ub <- DailyFromCumulative(x = temp_pred$ub,label_cummul,ind_ytype,av_num = delaydatai)
     # temp_pred$mean_predict <- DailyFromCumulative(x = temp_pred$mean_predict,label_cummul,ind_ytype,av_num = delaydatai)
      
        
      #colnames(temp_pred)<-paste(c('Date','mean_predict','std_predict'),ind_ytype,sep='')
      if(main_struct$write_csv_from_plots==1)
      {
        write.table(x=temp_pred,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"Simulations",ind_ytype,add_name1,".csv",sep=''))
      }
      if((label_all==0))
      {
        #ylims <-c(min(0.7*min(coll_pred[[ind_ytype]]$mn[relev_indi]-coll_pred[[ind_ytype]]$sn[relev_indi]),ylims[1]),max(1.2*max(coll_pred[[ind_ytype]]$mn[relev_indi]+coll_pred[[ind_ytype]]$sn[relev_indi]),ylims[2]))
        #min(coll_pred[[ind_ytype]]$mn[relev_indi]-coll_pred[[ind_ytype]]$sn[relev_indi])
        if((label_cummul==1)||(ind_ytype!="Critical"))
        {
          ylims <- c(min(0.7*min(temp_pred$lb),ylims[1],na.rm=T),max(1.2*max(temp_pred$ub,ylims[2],na.rm=T)))
        }else{
          ylims <- c(min(0.7*min(temp_pred$lb),ylims[1],na.rm=T),max(1.2*max(temp_pred$ub,ylims[2],na.rm=T)))
        }
        #min(temp_pred$mean_predict-temp_pred$std_predict)
      }
      if((label_altern==1)&(length(coll_pred$altern_sc)>0))
      {
        if(log_option==1)#when log option, upper limit is srawn
        {
          ylims <- c(min(0.7*min(help_v),ylims[1],na.rm=T),max(1.2*max(ylims[2],help_v,na.rm=T),na.rm=T))
        }else{
          ylims <- c(min(0.7*min(help_v[,1]),ylims[1],na.rm=T),max(1.2*max(ylims[2],help_v[,1],na.rm=T),na.rm=T))
        }
        #ylims <- c(min(0.7*min(help_v),ylims[1]),max(1.2*max(help_v)))
      }
      if(log_option==1)
      {  
        ylims[1]<-max(1,ylims[1])
        plot(xx,yy,col = "blue",pch = 1,ylim=ylims,xlim=xlims,xlab="Time, days", ylab =currentylab,cex=2,cex.lab=1.5, cex.axis=1.2,lwd=2,log= "y")#pch = ind_ytype0
        
      }else{
        plot(xx,yy,col = "blue",pch = 1,ylim=ylims,xlim=xlims,xlab="Time, days", ylab =currentylab,cex=2,cex.lab=1.5, cex.axis=1.2,lwd=2)#pch = ind_ytype0
      }
      if(label_all==1)
      {  
        for(ind_sol in 1:dim(coll_pred$arr)[3])
        {
          sim_res<-coll_pred$arr[,,ind_sol]
          sim_res<-NameOutput(sim_res,main_struct)#colnames(sim_res)<-main_struct$colnames_sim_res##NameOutput(sim_res,main_struct)
          
          zzdense<-CalculateOutput(S=sim_res,ind_ytype,param_)#sim_res[relev_indi,]
          zzdense <- SimpDailyFromCumulative(x = zzdense,label_cummul,ind_ytype,av_num = delaydatai)
          
          if(min(zzdense>=0))
          {  
            lines(xx0,zzdense,col = rgb(1,0.3,0),lwd=1)
          }
        }
      }
      points(xx,yy,col = "blue",pch = 1,ylim=ylims,xlim=xlims,cex=2,lwd=2)#,log= "y")#pch = ind_ytype0
      ##
    #  if((label_cummul==1)||((label_cummul==1)&(ind_ytype!="Critical")))
    #  {
        if((label_cummul==0))#&(ind_ytype!="Critical"))
        {
          AddCIBarstoPlot(t_points=xx0,l_points= coll_pred[[ind_ytype]]$PerDay$quantile[,1],#temp_pred$lb,#coll_pred[[ind_ytype]]$PerDay$lb,#[relev_indi]
                          u_points= coll_pred[[ind_ytype]]$PerDay$quantile[,dim(coll_pred[[ind_ytype]]$PerDay$quantile)[2]],col_='black')# temp_pred$ub[relev_indi]
          if(main_struct$write_csv_from_plots==1)
          {
            temp_CI<-cbind(temp_pred$lb,temp_pred$ub)
            colnames(temp_CI)<-paste(c('LB','UB'),'_',ind_ytype,sep='')
            write.table(x=temp_CI,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"ConfidenceIntervalsDaily",ind_ytype,add_name1,".csv",sep=''))
          }
         # AddErrorBarstoPlot(t_points=xx0,x_points=coll_pred[[ind_ytype]]$mn[relev_indi],
          #                   bar_points=coll_pred[[ind_ytype]]$sn[relev_indi],col_='black')
        }else
        {  
          #AddErrorBarstoPlot(t_points=xx0,x_points=coll_pred[[ind_ytype]]$mn[relev_indi],
          #               bar_points=coll_pred[[ind_ytype]]$sn[relev_indi],col_='black')
          AddCIBarstoPlot(t_points=xx0,l_points= coll_pred[[ind_ytype]]$quantile[,1],#temp_pred$lb,#coll_pred[[ind_ytype]]$PerDay$lb,#[relev_indi]
                          u_points= coll_pred[[ind_ytype]]$quantile[,dim(coll_pred[[ind_ytype]]$quantile)[2]],col_='black')# temp_pred$ub[relev_indi]
          
          
          #AddCIBarstoPlot(t_points=xx0,l_points=coll_pred[[ind_ytype]]$lb,#[relev_indi]
           #               u_points=coll_pred[[ind_ytype]]$ub,col_='black')#[relev_indi]
          if(main_struct$write_csv_from_plots==1)
          {
            temp_CI<-cbind(coll_pred[[ind_ytype]]$lb,coll_pred[[ind_ytype]]$ub)
            colnames(temp_CI)<-paste(c('LB','UB'),'_',ind_ytype,sep='')
            write.table(x=temp_CI,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"ConfidenceIntervals",ind_ytype,add_name1,".csv",sep=''))
          }
        }
     # }
      ##
      lines(xx0,temp_pred$mean_predict,col = rgb(0.6,0.6,0),lwd=2)#coll_pred[[ind_ytype]]$mn[relev_indi]
      ind_pred<-xx0>max(xx)
      lines(xx0[ind_pred],temp_pred$mean_predict[ind_pred],col = rgb(0,1,0),lwd=3,lty=3)#coll_pred[[ind_ytype]]$mn[relev_indi][ind_pred]
      
      if(label_altern==1)
      {
        if(length(coll_pred$altern_sc)>0)
        {  
          lines(xx0,help_v[,1],col='red',lwd=3,lty=1)#coll_pred$altern_sc[[ind_ytype]]$best[,1]
          if(log_option==1)#when log option, upper limit is srawn
          {#coll_pred$altern_sc[[ind_ytype]]$best[,2]
        ##    lines(xx0,help_v[,2],col='red',lwd=3,lty=1)
          }
          if(dim(help_v)[2]>2)
          {
            for(ind_sc in 3:dim(help_v)[2])
            {#coll_pred$altern_sc[[ind_ytype]]$best[,ind_sc]
              lines(xx0,help_v[,ind_sc],col='green',lwd=3,lty=1)
            }
            
          }
        }
        if((main_struct$write_csv_from_plots==1)&(length(coll_pred$altern_sc)>0))#&(dim(help_v)[2]>2))
        {
          temp_alt<- help_v[,1]#help_v[,1:dim(help_v)[2]]
          #colnames(temp_alt)<-paste(main_struct$val_altern[1:dim(help_v)[2]],'_',ind_ytype,sep='')
          #write.table(x=temp_alt,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"AlternativeShutdown",ind_ytype,add_name1,".csv",sep=''))
          write.table(x=temp_alt,sep = ";" ,  dec = main_struct$dec, row.names=F,file=paste(main_struct$resultspath,"MostProbable",ind_ytype,add_name1,".csv",sep=''))
        }
      }
      restrict1 <- min(all_times)+param_$Lock1
      restrict2 <- min(all_times)+param_$Lock2
      if(main_struct$label_treat_info==0)
      {  
        restrict3 <- min(all_times)+param_$Unlock
      }else{
        restrict3 <- min(all_times)+param_$Unlock1
      }

      #restrict3 <- min(all_times)+param_$Lock3
      if(min(xx)<restrict1)
      {
        lines(c(restrict1,restrict1),ylims,col='black',lwd=3,lty=3)
        
      }
      if(min(xx)<restrict2)
      {
        
        lines(c(restrict2,restrict2),ylims,col='black',lwd=3,lty=1)
      }
      if(min(xx)<restrict3)
      {
      #  
        lines(c(restrict3,restrict3),ylims,col='red',lwd=3,lty=1)
      }
      points(xx,yy,col = "blue",pch = 1,ylim=ylims,xlim=xlims,xlab="Time, days", ylab =currentylab,cex=2,cex.lab=1.5, cex.axis=1.2,lwd=2)
      min(all_times)+param_$Lock1
      
    }
  }
  ####
  if(label_save==1)
  {
    dev.off()
  }
}
DiffCalcl<-function(x)
{
  x[2:length(x)]-x[1:(length(x)-1)]
}
SumCalcl<-function(x)
{
  le<-length(x)
  out_<-vector(length= le)
  out_[1]<-x[1]
  for(indi in 2:le)
  {
    out_[indi] <- out_[indi-1]+x[indi]
  }
  return(out_)
}
DiffCalcl1<-function(x)
{
  c(0,x[2:length(x)]-x[1:(length(x)-1)])
}
DiffCalclAdv<-function(x,av_num)
{
  le<-length(x)
  out_<-vector(length= le)
  x_d<-DiffCalcl1(x)
  if(av_num>0)
  {  
    for(indi in 1:le)
    {
      out_[indi]<-mean(x_d[max(1,indi-av_num):min(indi+av_num,le)])#mean(x_d[max(1,indi-av_num):indi])
    }
  }else{
    out_<-x_d
  }
  out_
}
AverSim<-function(x,av_num)
{
  le<-length(x)
  out_<-vector(length= le)
  for(indi in 1:le)
  {
    out_[indi]<-mean(x[max(1,indi-av_num):min(indi+av_num,le)])
  }
  out_
}
AverAdv<-function(x,av_num)
{
  le<-length(x)
  if(av_num>0)
  {
    out_<-vector(length= le)
    for(indi in 1:le)
    {
      out_[indi]<-mean(x[max(1,indi-av_num):indi])
    }
  }else{
    out_<-x
  }
  out_
}
DiffCalclArr<-function(x)
{
  dimi<-dim(x)
  dimi1<-dim(x)
  dimi1[1] <-  dimi1[1]-1
  out_<-array(dim=dimi1)
  if(length(dimi1)==2)
  {
    out_ <- x[2:dimi[1],]-x[1:(dimi[1]-1),]
  }
  if(length(dimi1)==3)
  {
    out_ <- x[2:dimi[1],,]-x[1:(dimi[1]-1),,]
  }
  out_
}
AddCIBarstoPlot<- function(t_points,l_points,u_points,col_)#,ylims)
{
  for(indi in 1: length(t_points))
  {
    x_temp<-c(t_points[indi],t_points[indi])
    y_temp<-c(l_points[indi],u_points[indi])
    #ylims <- c(min(y_temp,ylims[1]),max(y_temp,ylims[2]))
    lines(x_temp,y_temp,col=col_,lwd=1.5)#,ylim=ylims)
  }
  
}

