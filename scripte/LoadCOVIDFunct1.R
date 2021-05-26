AdvancedLoad <- function(path_code,root_path,label_data, label_clean, ind_indiv, inputdec )
{
  print("start") 
  data_name <- '040_datint_rki_meldedat_kreise_bl_d_2020-10-05.txt'
  #path1<-paste(root_path,'/StochODEMS_YK/',sep='')
  if(label_data=='Land')
  {
    path_model<-paste(path_code,'AllLandsTogether/',sep='')
  }else if(label_data=='Kreis')
  {
    path_model<-paste(path_code,'AllKreise/',sep='')
  }else if(label_data=='LandPrior')
  {
    path_model<-paste(path_code,'AllLAndsTogetherPrior/',sep='')
  }else if (label_data=='LandAdv')
  {
    path_model<-paste(path_code,'AdvancedPrior/',sep='')
  }
  else if (label_data== 'LandKI')
  {
    path_model<-paste(path_code,'KIPrior/',sep='')
  }
  else if (label_data== 'LandKIAdv')
  {
    if(label_clean==1)
    {  
      path_model<-paste(path_code,'KICleaned/',sep='')
    }else{
      path_model<-paste(path_code,'KIAdvanced/',sep='')
    }
  }
  else if (label_data== 'LandMuKIAdv')
  {
    if(label_clean==1)
    {  
      path_model<-paste(path_code,'KIMutantCleaned/',sep='')
    }else{
      path_model<-paste(path_code,'KIMutant/',sep='')
    }
  }
  else if (label_data== 'LandMuPolyKIAdv')
  {
    #if(label_clean==1)
    path_model<-paste(path_code,'KIMutPoly/',sep='')
  }else if(label_data=='LandMuPolyKIAdvDiag')
  {
    path_model<-paste(path_code,'KIMuPolyDiag/',sep='')
  }
  setwd(path_code)
  source('PredictionFunction.R')
  source("GeneralFunctions/MCMCMonolix.R")
  source("GeneralFunctions/HookJeevesLib.R")
  source("GeneralFunctions/MCMCFunctionsPoiesis.R") 
  
  source("PresentResultsFunctions.R")
  #source("LoadCOVIDFunct1.R")#source("LoadCOVIDFunct.R") 
  #source("LoadFunct.R")
  source("AuxFunct.R") 
  source("CreateUpdateParameters.R")
  source("COVIDfunctions.R")
  source("MCMCMonolixSpecialCovid.R")
  source("AddCOVIDfunctions.R")
  source("AddCOVIDfunctionsTreatments.R")
  source("AddCOVIDfunctions1.R")
  source("MutantModels.R")
  source('AddCOVIDfunctionsTreatments.R')
  source("SimulateCovid.R")
  source("MCMCAlternativesMixture.R")
  source("SAEMMixtureFunctions.R")
  #'040_datint_rki_meldedat_kreise_bl_d_2020-05-18.txt'
  
  main_struct <- UploadDataDefinitionsCOVID(ind_indiv,path_model,root_path,
                                            file_data_name = data_name,label_data)#,#"D:/2003_covid_modell"
  
  main_struct$label_cummul <- 0#0
  main_struct$write_csv_from_plots <- 0#0#1
  main_struct$LogSc <-0#1
  ####advanced:
  relev_file_name <- FindActualData(ind_indiv,root_path,label_clean,label_age)#relev_file_name <- FindActualData(ind_indiv,root_path,label_clean=0,label_age)
  main_struct$weekly <- T#F
  main_struct$inputdec <- inputdec#'.'#','
  #main_struct$dec <- ','
  main_struct <- LoadECDCData(root_path,ecdc_file_data_name= relev_file_name,main_struct,ind_indiv)
 
  
  #main_struct <- LoadECDCData(root_path,ecdc_file_data_name='040_datint_ecdc_saxony_2020-09-09_1.txt',main_struct)
  main_struct$acurr_cumul <-main_struct$acurr#*10
  if(!'label_treat_info'%in%names(main_struct))
  {
    main_struct$label_treat_info<-0
  }
  main_struct$write_csv_from_plots<-0
  main_struct <- UploadAgeStruct(ind_indiv,root_path,main_struct)###Problems!!!???
  main_struct$week_label <- 0#!!!!1#0
  main_struct$label_branch<- 2#1#2#0
  main_struct$delay_report <-1
  #4#3
  ###################
  phiopt<-GeneratePhiopt(main_struct)
  
  param_<-ConstructParam(x=phiopt,ind_indiv,main_struct) 
  
  main_struct <-UpdateMainstuctVparams (param_,main_struct,label_data)
  stoch_obj0 <-CreateStochProtocolDataPrimitive(main_struct,ind_indiv)
  xtot<-InitParam(x=phiopt$indiv[[ind_indiv]],ind_indiv ,phiopt,stoch_obj0,main_struct)#phiopt$pop
  xtot$phiopt <-phiopt
  main_struct$speclength <- xtot$num_steps
  param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct) 
  main_struct$param_ <- param_
  out_<-list()
  out_$main_struct <- main_struct
  
  out_$stoch_obj0 <- stoch_obj0
  out_$xtot <- xtot
  out_$param_ <- param_
  ###
  return(out_)
  #return(main_struct)
}
SimpleLoad <- function(path_code,root_path,data_name,label_data, label_clean)
{
  #path1<-paste(root_path,'/StochODEMS_YK/',sep='')
  if(label_data=='Land')
  {
    path_model<-paste(path_code,'AllLandsTogether/',sep='')
  }else if(label_data=='Kreis')
  {
    path_model<-paste(path_code,'AllKreise/',sep='')
  }else if(label_data=='LandPrior')
  {
    path_model<-paste(path_code,'AllLAndsTogetherPrior/',sep='')
  }else if (label_data=='LandAdv')
  {
    path_model<-paste(path_code,'AdvancedPrior/',sep='')
  }
  else if (label_data== 'LandKI')
  {
    path_model<-paste(path_code,'KIPrior/',sep='')
  }
  else if (label_data== 'LandKIAdv')
  {
    if(label_clean==1)
    {  
      path_model<-paste(path_code,'KICleaned/',sep='')
    }else{
      path_model<-paste(path_code,'KIAdvanced/',sep='')
    }
  }
  else if (label_data== 'LandMuKIAdv')
  {
    if(label_clean==1)
    {  
      path_model<-paste(path_code,'KIMutantCleaned/',sep='')
    }else{
      path_model<-paste(path_code,'KIMutant/',sep='')
    }
  }
  else if (label_data== 'LandMuPolyKIAdv')
  {
    #if(label_clean==1)
    path_model<-paste(path_code,'KIMutPoly/',sep='')
  }else if(label_data=='LandMuPolyKIAdvDiag')
  {
    path_model<-paste(path_code,'KIMuPolyDiag/',sep='')
  }
  setwd(path_code)
  source("GeneralFunctions/MCMCMonolix.R")
  source("GeneralFunctions/HookJeevesLib.R")
  source("GeneralFunctions/MCMCFunctionsPoiesis.R") 
  
  source("PresentResultsFunctions.R")
  #source("LoadCOVIDFunct1.R")#source("LoadCOVIDFunct.R") 
  #source("LoadFunct.R")
  source("AuxFunct.R") 
  source("CreateUpdateParameters.R")
  source("COVIDfunctions.R")
  source("MCMCMonolixSpecialCovid.R")
  source("AddCOVIDfunctions.R")
  source("AddCOVIDfunctionsTreatments.R")
  source("AddCOVIDfunctions1.R")
  source("MutantModels.R")
  #'040_datint_rki_meldedat_kreise_bl_d_2020-05-18.txt'
  
  main_struct <- UploadDataDefinitionsCOVID(ind_indiv,path_model,root_path,
                                             file_data_name = data_name,label_data)#,#"D:/2003_covid_modell"
  
  main_struct$label_cummul <- 0#0
  main_struct$write_csv_from_plots <- 1#0#1
  main_struct$LogSc <-0#1
  if(!dir.exists(main_struct$resultspath))
  {
    dir.create(main_struct$resultspath)
  }
  #4#3
  main_struct$model_opt <- 31#32
  main_struct$treatpardir <- 2
  sourceCpp('solvecovidmodelMuKI.cpp')
  return(main_struct)
} 
UploadDataDefinitionsCOVID<-function(ind_indiv,path_model,root_path,file_data_name,label_data)#root_gen_funct,path_critical
{ 
  path_critical<-paste(root_path,'/data/',sep ='')#"R:/modellclub/2003_covid_modell/data/"
  path_data<-paste(root_path,'/data/',sep ='')
  log_option <- 1#1#0
  date_option <- 1#0#0 1
  label_all <-0
  main_struct <- ReadDataAndDefinitionsGUI(path_model,root_path,file_data_name,label_data)#(path_model=path1)#path_data for options main_struct<-LoadModelPerPat(path0 =path0,path1 = path1, newpat = 0)#main_struct basic definitions without measurements, compartments and treatments also to the case of biological data
  #main_struct <- ReadCOVIDData(path_model=path1, path_data= path_data,# "D:/2003_covid_modell/data/", 
  #                             path_critical, 
  ##                            file_data_name ,# "covid_data_sachsen_sms_200403.txt",
  #                             main_struct=main_struct)
  if((label_data== 'LandMuKIAdv')||(label_data== "LandMuPolyKIAdv")||(label_data== "LandMuPolyKIAdvDiag"))
  {
    main_struct$model_opt<- 31
  }else{
    main_struct$model_opt<- 3
  }
  main_struct <- AssignCompartmentsSimple(ind_indiv,path_model=path_model, main_struct)#,path_data=path1
  
  name_dir <- sub(x=file_data_name,pattern = '.txt',replacement = '')
  main_struct$resultspath<-paste(root_path,'/PredictionsValidations','/',name_dir,'/',sep='')
  
  if(!dir.exists(main_struct$resultspath))
  {
    dir.create(main_struct$resultspath)
  }
  main_struct$dec<-','
  #main_struct$model_opt<- 3#!!!!2#1
  #main_struct$model_opt<-2
  ind_indiv <- 1
  if(main_struct$model_opt==1)
  {
    main_struct$numvar<-7
  }else{
    main_struct$numvar<- main_struct$compartments[[ind_indiv]]$NumofVariables#8
  }
  #if(main_struct$model_opt==2)#Many compartments
  #{
  #  main_struct$numvar<- main_struct$compartments[[ind_indiv]]$NumofVariables#8
 # }
 # if(main_struct$model_opt==3)
 # {
 #   main_struct$numvar<- main_struct$compartments[[ind_indiv]]$NumofVariables
 # }
  #main_struct$weights<- rep.int(x = 1, times= main_struct$YTYPES$num)
  #main_struct <- UploadAgeStruct(root_path,main_struct)
  ####output:
  return(main_struct)
  
}

ReadDataAndDefinitionsGUI<-function(path_model,root_path,file_data_name,label_data){#(path0, main_struct) ,path_data
  #setwd(path_model)#path0
  if((label_data=='Land')||(label_data=='LandPrior')||(label_data=='LandAdv')||(label_data=='LandKI')||(label_data=='LandKIAdv')||(label_data=="LandMuKIAdv")||(label_data=='LandMuPolyKIAdv')||(label_data=='LandMuPolyKIAdvDiag'))
  {
    main_struct<-ReadCovariatesCOVID(root_path,file_data_name,path_model)#list()#initialisation, in matlab version:  handles.main_struct
  }else if(label_data=='Kreis')
  {
    main_struct<-ReadKreiseCovariatesCOVID(root_path,file_data_name)#list()#initialisation, in matlab version:  handles.main_struct
  }
  #####
  optim_options <- read.csv(paste(path_model,'simulationoptionshookjeeves.csv',sep=''), sep = ";",dec = ',')
  main_struct$optim_options<-list()
  for(opt_name in as.character(optim_options$Name))
  {
    main_struct$optim_options[[opt_name]]<-optim_options$Val[as.character(optim_options$Name)==opt_name]
  }
  main_struct$nonmem_options0<-read.csv(paste(path_model,'simulationoptionsmcmc.csv',sep=''), sep = ";",dec = ',')
  
  
  main_struct<-LoadIntegralModelCOVID(path_model, main_struct)#(path0, main_struct)
  if("weights.csv"%in%list.files(path = (path_model)))
  {
    temp_weights<-read.csv(paste(path_model,'weights.csv',sep = ''), sep = ";",dec = ',')
    if(dim(temp_weights)[2]-1==main_struct$numoutcomesnames){
      for (ind_names in 1:main_struct$numoutcomesnames)
      {
        main_struct$badoutcomes<-main_struct$badoutcomes+1-(main_struct$outcomesnames[ind_names]==names(temp_weights)[ind_names+1])
      }
    }else{
      main_struct$badoutcomes<-1
    }
    main_struct$weights<-temp_weights[,2:(main_struct$numoutcomesnames+1)]
  }else{
    main_struct$weights<- rep.int(x = 1, times= main_struct$YTYPES$num)
  }
  
  
  main_struct
}
ReadCovariatesCOVID<-function(root_path,file_data_name,path_model)
{  

  #'040_datint_rki_meldedat_kreise_bl_d_2020-05-18.txt'
  #'040_datint_rki_meldedat_kreise_bl_d_20-05-14.txt'
  #''040_datint_rki_meldedat_kreise_bl_d_2020-05-18.txt'
  together_data<-read.csv(paste(root_path,'/data/',file_data_name,sep =''), sep = "\t",dec = ',',header =T)
  #together_data$CountryExp
  #together_data$DateRep <- as.Date(as.character(temp_critical$DateRep),format="%Y-%m-%d")
  #ind_critical<-!is.na(together_data$covid_inICU_upscaled)
  temp_measurements0 <- list()
  measurements_daily <- list()
  measurements <- list()
  #measurements$RelevClinData <-list()
  #unique(together_data$CountryExp)[260:270]
  
  bundeslands<- as.character(unique(together_data$CountryExp[together_data$level=='Bundesland']))
  #main_struct$covariates <- data.frame("ID"= numeric(), "Land" = character(), "Sc0" = numeric(),
  #                                       "DateLock1" = as.Date(character()))
  num_lands <- length(bundeslands)+1
  ####
  kreise<- as.character(unique(together_data$CountryExp[together_data$level=='Kreis']))
  #main_struct$covariates <- data.frame("ID"= numeric(), "Land" = character(), "Sc0" = numeric(),
  #                                       "DateLock1" = as.Date(character()))
  num_kreise <- length(kreise)
  kreise_land<-kreise#initialization
  for(ind_kreise in 1:num_kreise)
  {
    ind_relev_land<-(together_data$level=='Kreis')&(together_data$CountryExp==kreise[ind_kreise])
    together_data$CountryExp[(together_data$level=='Kreis')]
    kreise_land[ind_kreise]<-as.character(together_data$Bundesland[ind_relev_land][1])
  }
  ####
  main_struct <- list()
  main_struct$measurements <- list()
  main_struct$completedates <- list()
  main_struct$completeweekdates <- list()
  
  main_struct$kreise_land <- kreise_land
  
  main_struct$numindiv <- num_lands
  main_struct$erstsonntagind <- vector(length = main_struct$numindiv)
  main_struct$letztsammstagind <- vector(length = main_struct$numindiv)
  
  #the following is predefined, raather than red from error data
  alloutcomes<- c("Total","Death","Critical")#as.character(temp_measurements$YTYPE)
  main_struct$outcomesnames <- alloutcomes
  alloutcomesnum<-vector(length = length(alloutcomes))
  main_struct$numoutcomesnames <- length(alloutcomes)
  listhelp=1:length(alloutcomes)
  for (ind_names in 1:main_struct$numoutcomesnames)
  {
    alloutcomesnum[listhelp[(main_struct$outcomesnames[ind_names]==alloutcomes)]]<-ind_names; 
  }
  main_struct$outcomesnumlist=1:main_struct$numoutcomesnames
  main_struct$YTYPES<-list()
  main_struct$YTYPES$listnum<-unique(sort(alloutcomesnum) )#unique(sort(alloutcomesnum[index_id]) )!!!
  main_struct$YTYPES$list<-main_struct$outcomesnames[main_struct$YTYPES$listnum]
  main_struct$YTYPES$num <- main_struct$numoutcomesnames
  main_struct$covariates <- list()#list("ID", "Land", "Sc0", "DateLock1", "DateLock2", "DateLock3",
  
  main_struct$covariates$ID <- vector(length = num_lands)
  main_struct$covariates$Land <- character(length = num_lands)
  main_struct$covariates$Sc0 <- vector(length = num_lands)
  together_data<-read.csv(paste(root_path,'/data/',file_data_name,sep =''), sep = "\t",dec = ',',header =T)
  #if('DateLock.csv'%in%list.files(path= paste(root_path,'/StochODEMS_YK/AllLAndsTogetherPrior',sep ='')))#!!!!
  if('DateLock.csv'%in%list.files(path= path_model))#!!!!
      
  {
    main_struct$label_treat_info <- 1
   # DateLockTemp<-read.csv(paste(root_path,'/StochODEMS_YK/AllLAndsTogetherPrior/','DateLock.csv',sep =''), sep = "\t",dec = ',',header =T)
   # DateUnlockTemp<-read.csv(paste(root_path,'/StochODEMS_YK/AllLAndsTogetherPrior/','DateUnlock.csv',sep =''), sep = "\t",dec = ',',header =T)
    
    DateLockTemp<-read.csv(paste(path_model,'DateLock.csv',sep =''), sep = "\t",dec = ',',header =T)
    DateUnlockTemp<-read.csv(paste(path_model,'DateUnlock.csv',sep =''), sep = "\t",dec = ',',header =T)
    
    
    main_struct$num_locks <- length(DateLockTemp$DateLock)
    main_struct$num_unlocks <- length(DateUnlockTemp$DateUnlock)#as.Date(DateLockTemp)
    for(ind_locki in 1:main_struct$num_locks)
    {
      main_struct$covariates[[paste('DateLock',ind_locki,sep='')]] <- rep.int(x = as.character(DateLockTemp$DateLock[ind_locki]),times = num_lands)#character(length = num_lands)
      main_struct$covariates[[paste('DateLock',ind_locki,sep='')]] <- as.Date(main_struct$covariates[[paste('DateLock',ind_locki,sep='')]],tryFormats = c("%d.%m.%Y"))  
    }
    for(ind_locki in 1:main_struct$num_unlocks)
    {
      main_struct$covariates[[paste('Unlock',ind_locki,sep='')]] <- rep.int(x = as.character(DateUnlockTemp$DateUnlock[ind_locki]),times = num_lands)#character(length = num_lands)
      main_struct$covariates[[paste('Unlock',ind_locki,sep='')]] <- as.Date(main_struct$covariates[[paste('Unlock',ind_locki,sep='')]],tryFormats = c("%d.%m.%Y"))
    }
  }else if('DatesUncertainty.csv'%in%list.files(path= path_model))
  {
    main_struct$label_treat_info <- 2
  }
    else{
    main_struct$label_treat_info <- 0
    main_struct$covariates$DateLock1 <- character(length = num_lands)
    main_struct$covariates$DateLock2 <- character(length = num_lands)
    main_struct$covariates$DateLock3 <- character(length = num_lands)
    main_struct$covariates$DateStart <- character(length = num_lands)
    #main_struct$covariates$DelayData <- vector(length = num_lands)
    main_struct$covariates$Unlock <- character(length = num_lands)
    main_struct$covariates$Unlock1 <- character(length = num_lands)
  }
  ####step functions pcrit pdeaths
  if('DateStepsPDeath.csv'%in%list.files(path= path_model))#!!!!
  {
    main_struct$covariates$DateStepsPDeath <- list()
    main_struct$covariates$DateStepsPCrit <- list()
    main_struct$label_step_death <- 1
    DatePdeathTemp <- read.csv(paste(path_model,'DateStepsPDeath.csv',sep =''), sep = "\t",dec = ',',header =T)
    DatePcritTemp <-read.csv(paste(path_model,'DateStepsPCrit.csv',sep =''), sep = "\t",dec = ',',header =T)
    main_struct$num_pdeath <- length(DatePdeathTemp$DateStepsPDeath)
    main_struct$num_pcrit <- length(DatePcritTemp$DateStepsPCrit)#as.Date(DateLockTemp)
    for(indi in 1: num_lands)
    {  
      main_struct$covariates$DateStepsPDeath[[indi]] <- as.Date(DatePdeathTemp$DateStepsPDeath,tryFormats = c("%d.%m.%Y")) 
      main_struct$covariates$DateStepsPCrit[[indi]] <-  as.Date(DatePcritTemp$DateStepsPCrit,tryFormats = c("%d.%m.%Y")) 
    }
  }
  ####
  main_struct$Daily<-list()
  #ind_<-1
  
    ind_ <- 0
  for(landi in c('Land',bundeslands))
  {
    ind_ <- ind_ +1
    if(ind_==1)
    {
      ind_relev <- (together_data$level== 'Land')
      main_struct$covariates$Land[ind_] <- "Deutschland"
    }else{
      ind_relev <- (together_data$level== 'Bundesland')&(together_data$CountryExp== landi)
      main_struct$covariates$Land[ind_] <- landi
    }
    
    main_struct$covariates$ID[ind_] <- ind_
    
    main_struct$covariates$Sc0[ind_] <- unique(as.numeric(as.character(together_data$einwohnergesamt[ind_relev])))
    for(nami in names(together_data))
    {
      temp_measurements0[[nami]]<-together_data[[nami]][ind_relev]
    }
    temp_measurements0$Meldedatum <- as.Date(together_data$DateRep[ind_relev])
    temp_measurements0$DateRep<- as.Date(together_data$DateRep[ind_relev])
    temp_measurements0$covid_inICU_upscaled<-as.numeric(as.character(temp_measurements0$covid_inICU_upscaled))
    date_start <- max(temp_measurements0$Meldedatum[together_data$AllConfCases[ind_relev]==min(together_data$AllConfCases[ind_relev],na.rm=T)],na.rm=T)
    
    
    main_struct$covariates$DateStart[ind_] <- as.character(date_start)
    if(main_struct$label_treat_info==0)
    {  
      main_struct$covariates$DateLock1[ind_] <- "2020-03-10"#as.character(main_struct$covariates[1,'DateLock1'])
      main_struct$covariates$DateLock2[ind_] <- "2020-03-16"#as.character(main_struct$covariates[1,'DateLock2'])
      main_struct$covariates$DateLock3[ind_] <- "2020-12-12"#as.character(main_struct$covariates[1,'DateLock3'])
      main_struct$covariates$Unlock[ind_] <- "2020-04-30"# und "2020-06-01" "2020-04-28"#"2020-05-01"#as.character(main_struct$covariates[1,'Unlock'])
      main_struct$covariates$Unlock1[ind_] <- "2020-06-01"
    }
    #main_struct$covariates$DelayData[ind_] <- 2# To special parameters!!!!!#main_struct$covariates[1,'DelayData']
    
    #######
    #####
    ind_allcases<-!is.na(temp_measurements0$AllConfCases)
    ind_deaths<-!is.na(temp_measurements0$AllDeaths)
    ind_cases_daily<-!is.na(temp_measurements0$NewConfCases)
    ind_deaths_daily<-!is.na(temp_measurements0$NewDeaths)
    ind_time_relev <- (temp_measurements0$DateRep)>= as.Date(main_struct$covariates$DateStart[ind_])
    ind_critical <- !is.na(temp_measurements0$covid_inICU_upscaled)
    ind_critical <- ind_critical&(ind_time_relev)
    
    
    ind_allcases <- ind_allcases&(ind_time_relev)
    ind_deaths <- ind_deaths&(ind_time_relev)
    
    ind_cases_daily <- ind_cases_daily&(ind_time_relev)
    ind_deaths_daily <- ind_deaths_daily&(ind_time_relev)
    
    dates_critical<-temp_measurements0$DateRep[ind_critical]
    num_critical<-temp_measurements0$covid_inICU_upscaled[ind_critical]
    num_all1<-temp_measurements0$AllConfCases[ind_critical]
    numindiv_obs<-sum(ind_deaths)+sum(ind_allcases)+sum(ind_critical)
    numindiv_obs_daily<-sum(ind_deaths_daily)+sum(ind_cases_daily)+sum(ind_critical)
    all_cases_date<-temp_measurements0$Meldedatum[ind_allcases]
    all_cases<-temp_measurements0$AllConfCases[ind_allcases]
    
    cases_date_daily<-temp_measurements0$Meldedatum[ind_cases_daily]#-temp_cov[['DelayData']]
    #deaths_date_daily<-temp_measurements0$Meldedatum[ind_deaths_daily]-temp_cov[['DelayData']]
    
    cases_daily<-temp_measurements0$NewConfCases[ind_cases_daily]
    if(min(temp_measurements0$DateRep)>main_struct$covariates$DateStart[ind_])
    {
      numindiv_obs <- numindiv_obs+1
      all_cases<-c(0,all_cases)
      all_cases_date<-c(main_struct$covariates$DateStart[ind_],all_cases_date)
    }
    if(min(cases_date_daily)>main_struct$covariates$DateStart[ind_])
    {
      numindiv_obs_daily <- numindiv_obs_daily+1
      cases_daily<-c(0,cases_daily)
      cases_date_daily<-c(main_struct$covariates$DateStart[ind_],cases_date_daily)
    }
    
    measurements$namesmeasdata<-c("T","Y")
    measurements$date<-c(all_cases_date,#temp_measurements0$Meldedatum[ind_allcases],
                         temp_measurements0$Meldedatum[ind_deaths],
                         dates_critical)
    measurements$data <- rbind(as.numeric(measurements$date-min(measurements$date)),
                               c(all_cases,#temp_measurements0$AllConfCases[ind_allcases]
                                 temp_measurements0$AllDeaths[ind_deaths],
                                 temp_measurements0$covid_inICU_upscaled[ind_critical]))
    measurements$datatypes <- c(rep.int(x="Total",times=length(all_cases)),#sum(ind_allcases)
                                rep.int(x="Death",times=sum(ind_deaths)),
                                rep.int(x="Critical",times=sum(ind_critical)))
    measurements$scope<- as.numeric(max(measurements$date)-min(measurements$date))
    #measurements$RelevClinData <- (dates_critical>"2020-04-19")
    ind_ytypes1<-0
    measurements$length_data<- vector(length = length(unique(measurements$datatypes)))
    for(ind_ytype in unique(measurements$datatypes))
    {
      ind_ytypes1 <- ind_ytypes1 +1
      observedi<-measurements$data[measurements$namesmeasdata=="Y",measurements$datatypes==ind_ytype]
      measurements$length_data[ind_ytypes1]<- sum(!is.na(observedi))
    }   
    names(measurements$length_data) <- unique(measurements$datatypes)
    main_struct$measurements[[ind_]]<-measurements
    main_struct$completedates[[ind_]]<-CompleteDates(unique(main_struct$measurements[[ind_]]$date))
    #if('DateStepsPDeath'%in%names(main_struct$covariates))
   # {
   #   main_struct$covariates$DateStepsPDeath[[ind_]]<-c(main_struct$covariates$DateStepsPDeath[[ind_]],max(main_struct$completedates[[ind_]]))
    
   #   main_struct$covariates$DateStepsPCrit[[ind_]]<-c(main_struct$covariates$DateStepsPCrit[[ind_]],max(main_struct$completedates[[ind_]]))
   # }
      ###Daily:
     
    measurements_daily$namesmeasdata<-c("T","Y")
    measurements_daily$date<-c(cases_date_daily,#temp_measurements0$Meldedatum[ind_allcases],
                               temp_measurements0$Meldedatum[ind_deaths_daily],
                               dates_critical)
    measurements_daily$data <- rbind(as.numeric(measurements_daily$date-min(measurements_daily$date)),
                                     c(cases_daily,#temp_measurements0$AllConfCases[ind_allcases]
                                       temp_measurements0$NewDeaths[ind_deaths_daily],
                                       temp_measurements0$covid_inICU_upscaled[ind_critical]))
    measurements_daily$datatypes <- c(rep.int(x="Total",times=length(cases_daily)),#sum(ind_allcases)
                                      rep.int(x="Death",times=sum(ind_deaths_daily)),
                                      rep.int(x="Critical",times=sum(ind_critical)))
    measurements_daily$RelevClinData <- (dates_critical>"2020-04-19")
    
    main_struct$Daily$measurements[[ind_]]<-measurements_daily
    
    
    main_struct$YTYPE[[ind_]]<-list()
    main_struct$YTYPE[[ind_]]$listnum<-unique(sort(alloutcomesnum) )#unique(sort(alloutcomesnum[index_id]) )!!!
    main_struct$YTYPE[[ind_]]$list<-main_struct$outcomesnames[main_struct$YTYPE[[ind_]]$listnum]
    main_struct$YTYPE[[ind_]]$num<-length(main_struct$YTYPE[[ind_]]$list)
    main_struct$YTYPE[[ind_]]$size<-list()
    main_struct$YTYPE[[ind_]]$index<-list()
    main_struct$observed[[ind_]]<-list()
    main_struct$prierr_indiv[[ind_]]<-list()
    main_struct$initname[[ind_]]<-list()
    main_struct$initval[[ind_]]<-vector(length = max(main_struct$YTYPE[[ind_]]$listnum) )
    
    for(ind_ytype1 in 1:main_struct$YTYPES$num)
    {
      ind_ytype<-main_struct$YTYPES$listnum[ind_ytype1]
      ind_ytype_symb<-main_struct$YTYPES$list[ind_ytype1]
      main_struct$YTYPE[[ind_]]$size[ind_ytype]<-sum(alloutcomesnum==main_struct$YTYPES$listnum[ind_ytype1]) #alloutcomesnum[index_id]
      main_struct$YTYPE[[ind_]]$index[[ind_ytype]]<-(alloutcomesnum==main_struct$YTYPES$listnum[ind_ytype1])#alloutcomesnum[index_id]
      #index_id_y_ytype is unnecessary, use main_struct$YTYPE[[ind_]]$size[ind_ytype] instead!!!!!
      ##obsolete and not correct!!!:
    #  if(sum(main_struct$YTYPE[[ind_]]$index[[ind_ytype]],na.rm = T)>0)
    #  {
    #    main_struct$observed[[ind_]][[ind_ytype]]<-main_struct$measurements[[ind_]]$data[2,][main_struct$YTYPE[[ind_]]$index[[ind_ytype]]]#main_struct$Y[[ind_]] #alternative to matlab definition !!!! 
          
    #  }
    }
    #####
    #######
  }
  main_struct$indiv <- list()
  main_struct$indiv$ID <- 1:main_struct$numindiv
  main_struct$ID <- main_struct$covariates$ID
  if(main_struct$label_treat_info==0)
  {  
    main_struct$covariates$DateLock1 <- as.Date(main_struct$covariates$DateLock1)
    main_struct$covariates$DateLock2 <- as.Date(main_struct$covariates$DateLock2)
    main_struct$covariates$DateLock3 <- as.Date(main_struct$covariates$DateLock3)
    main_struct$covariates$DateStart <- as.Date(main_struct$covariates$DateStart)
    main_struct$covariates$Unlock <- as.Date(main_struct$covariates$Unlock)
    main_struct$covariates$Unlock1 <- as.Date(main_struct$covariates$Unlock1)
    main_struct$covariates$Unlock2 <- main_struct$covariates$Unlock1#initialization
  } 
  
  return(main_struct)
}
####
ReadKreiseCovariatesCOVID<-function(root_path,file_data_name)
{  
  
  #'040_datint_rki_meldedat_kreise_bl_d_2020-05-18.txt'
  #'040_datint_rki_meldedat_kreise_bl_d_20-05-14.txt'
  #''040_datint_rki_meldedat_kreise_bl_d_2020-05-18.txt'
  together_data<-read.csv(paste(root_path,'/data/',file_data_name,sep =''), sep = "\t",dec = ',',header =T)
  #together_data$CountryExp
  #together_data$DateRep <- as.Date(as.character(temp_critical$DateRep),format="%Y-%m-%d")
  #ind_critical<-!is.na(together_data$covid_inICU_upscaled)
  temp_measurements0 <- list()
  measurements_daily <- list()
  measurements <- list()
  #measurements$RelevClinData <-list()
  #unique(together_data$CountryExp)[260:270]
  
  bundeslands<- as.character(unique(together_data$CountryExp[together_data$level=='Bundesland']))
  #main_struct$covariates <- data.frame("ID"= numeric(), "Land" = character(), "Sc0" = numeric(),
  #                                       "DateLock1" = as.Date(character()))
  num_lands <- length(bundeslands)+1
  ####
  kreise<- as.character(unique(together_data$CountryExp[together_data$level=='Kreis']))
  #main_struct$covariates <- data.frame("ID"= numeric(), "Land" = character(), "Sc0" = numeric(),
  #                                       "DateLock1" = as.Date(character()))
  num_kreise <- length(kreise)
  kreise_land<-kreise#initialization
  for(ind_kreise in 1:num_kreise)
  {
    ind_relev_land<-(together_data$level=='Kreis')&(together_data$CountryExp==kreise[ind_kreise])
    together_data$CountryExp[(together_data$level=='Kreis')]
    kreise_land[ind_kreise]<-as.character(together_data$Bundesland[ind_relev_land][1])
  }
  ####
  main_struct <- list()
  main_struct$measurements <- list()
  main_struct$kreise_land <- kreise_land
  
  main_struct$numindiv <- num_kreise
  #the following is predefined, raather than red from error data
  alloutcomes<- c("Total","Death","Critical")#as.character(temp_measurements$YTYPE)
  main_struct$outcomesnames <- alloutcomes
  alloutcomesnum<-vector(length = length(alloutcomes))
  main_struct$numoutcomesnames <- length(alloutcomes)
  listhelp=1:length(alloutcomes)
  for (ind_names in 1:main_struct$numoutcomesnames)
  {
    alloutcomesnum[listhelp[(main_struct$outcomesnames[ind_names]==alloutcomes)]]<-ind_names; 
  }
  main_struct$outcomesnumlist=1:main_struct$numoutcomesnames
  main_struct$YTYPES<-list()
  main_struct$YTYPES$listnum<-unique(sort(alloutcomesnum) )#unique(sort(alloutcomesnum[index_id]) )!!!
  main_struct$YTYPES$list<-main_struct$outcomesnames[main_struct$YTYPES$listnum]
  main_struct$YTYPES$num <- main_struct$numoutcomesnames
  main_struct$covariates <- list()#list("ID", "Land", "Sc0", "DateLock1", "DateLock2", "DateLock3",
  
  main_struct$covariates$ID <- vector(length = num_kreise)
  main_struct$covariates$Land <- character(length = num_kreise)
  main_struct$covariates$Sc0 <- vector(length = num_kreise)
  main_struct$covariates$DateLock1 <- character(length = num_kreise)
  main_struct$covariates$DateLock2 <- character(length = num_kreise)
  main_struct$covariates$DateLock3 <- character(length = num_kreise)
  main_struct$covariates$DateStart <- character(length = num_kreise)
  #main_struct$covariates$DelayData <- vector(length = num_kreise)
  main_struct$covariates$Unlock <- character(length = num_kreise)
  main_struct$covariates$Unlock1 <- character(length = num_kreise)
  main_struct$Daily<-list()
  #ind_<-1
  
  ind_kreise <- 0
  for(landi in kreise)#c('Land',bundeslands))
  {
    ind_kreise <- ind_kreise +1
  
    ind_relev <- (together_data$level=='Kreis')&(together_data$CountryExp==kreise[ind_kreise])#(together_data$level== 'Bundesland')&(together_data$CountryExp== landi)
    main_struct$covariates$Land[ind_kreise] <- landi
    
    
    main_struct$covariates$ID[ind_kreise] <- ind_kreise
    
    main_struct$covariates$Sc0[ind_kreise] <- unique(as.numeric(as.character(together_data$einwohnergesamt[ind_relev])))
    for(nami in names(together_data))
    {
      temp_measurements0[[nami]]<-together_data[[nami]][ind_relev]
    }
    temp_measurements0$Meldedatum <- as.Date(together_data$DateRep[ind_relev])
    temp_measurements0$DateRep<- as.Date(together_data$DateRep[ind_relev])
    temp_measurements0$covid_inICU_upscaled<-as.numeric(as.character(temp_measurements0$covid_inICU_upscaled))
    date_start <- max(temp_measurements0$Meldedatum[together_data$AllConfCases[ind_relev]==min(together_data$AllConfCases[ind_relev],na.rm=T)],na.rm=T)
    main_struct$covariates$DateStart[ind_kreise] <- as.character(date_start)
    main_struct$covariates$DateLock1[ind_kreise] <- "2020-03-10"#as.character(main_struct$covariates[1,'DateLock1'])
    main_struct$covariates$DateLock2[ind_kreise] <- "2020-03-16"#as.character(main_struct$covariates[1,'DateLock2'])
    main_struct$covariates$DateLock3[ind_kreise] <- "2020-12-12"#as.character(main_struct$covariates[1,'DateLock3'])
    main_struct$covariates$Unlock[ind_kreise] <- "2020-04-30"# und "2020-06-01" "2020-04-28"#"2020-05-01"#as.character(main_struct$covariates[1,'Unlock'])
    main_struct$covariates$Unlock1[ind_kreise] <- "2020-06-01"
    #main_struct$covariates$DelayData[ind_kreise] <- 2# To special parameters!!!!!#main_struct$covariates[1,'DelayData']
    
    #######
    #####
    ind_allcases<-!is.na(temp_measurements0$AllConfCases)
    ind_deaths<-!is.na(temp_measurements0$AllDeaths)
    ind_cases_daily<-!is.na(temp_measurements0$NewConfCases)
    ind_deaths_daily<-!is.na(temp_measurements0$NewDeaths)
    ind_time_relev <- (temp_measurements0$DateRep)>= as.Date(main_struct$covariates$DateStart[ind_kreise])
    ind_critical <- !is.na(temp_measurements0$covid_inICU_upscaled)
    ind_critical <- ind_critical&(ind_time_relev)
    
    
    ind_allcases <- ind_allcases&(ind_time_relev)
    ind_deaths <- ind_deaths&(ind_time_relev)
    
    ind_cases_daily <- ind_cases_daily&(ind_time_relev)
    ind_deaths_daily <- ind_deaths_daily&(ind_time_relev)
    
    dates_critical<-temp_measurements0$DateRep[ind_critical]
    num_critical<-temp_measurements0$covid_inICU_upscaled[ind_critical]
    num_all1<-temp_measurements0$AllConfCases[ind_critical]
    numindiv_obs<-sum(ind_deaths)+sum(ind_allcases)+sum(ind_critical)
    numindiv_obs_daily<-sum(ind_deaths_daily)+sum(ind_cases_daily)+sum(ind_critical)
    all_cases_date<-temp_measurements0$Meldedatum[ind_allcases]
    all_cases<-temp_measurements0$AllConfCases[ind_allcases]
    
    cases_date_daily<-temp_measurements0$Meldedatum[ind_cases_daily]#-temp_cov[['DelayData']]
    #deaths_date_daily<-temp_measurements0$Meldedatum[ind_deaths_daily]-temp_cov[['DelayData']]
    
    cases_daily<-temp_measurements0$NewConfCases[ind_cases_daily]
    if(min(temp_measurements0$DateRep)>main_struct$covariates$DateStart[ind_kreise])
    {
      numindiv_obs <- numindiv_obs+1
      all_cases<-c(0,all_cases)
      all_cases_date<-c(main_struct$covariates$DateStart[ind_kreise],all_cases_date)
    }
    if(min(cases_date_daily)>main_struct$covariates$DateStart[ind_kreise])
    {
      numindiv_obs_daily <- numindiv_obs_daily+1
      cases_daily<-c(0,cases_daily)
      cases_date_daily<-c(main_struct$covariates$DateStart[ind_kreise],cases_date_daily)
    }
    
    measurements$namesmeasdata<-c("T","Y")
    measurements$date<-c(all_cases_date,#temp_measurements0$Meldedatum[ind_allcases],
                         temp_measurements0$Meldedatum[ind_deaths],
                         dates_critical)
    measurements$data <- rbind(as.numeric(measurements$date-min(measurements$date)),
                               c(all_cases,#temp_measurements0$AllConfCases[ind_allcases]
                                 temp_measurements0$AllDeaths[ind_deaths],
                                 temp_measurements0$covid_inICU_upscaled[ind_critical]))
    measurements$datatypes <- c(rep.int(x="Total",times=length(all_cases)),#sum(ind_allcases)
                                rep.int(x="Death",times=sum(ind_deaths)),
                                rep.int(x="Critical",times=sum(ind_critical)))
    measurements$scope<- as.numeric(max(measurements$date)-min(measurements$date))
    #measurements$RelevClinData <- (dates_critical>"2020-04-19")
    ind_ytypes1<-0
    measurements$length_data<- vector(length = length(unique(measurements$datatypes)))
    for(ind_ytype in unique(measurements$datatypes))
    {
      ind_ytypes1 <- ind_ytypes1 +1
      observedi<-measurements$data[measurements$namesmeasdata=="Y",measurements$datatypes==ind_ytype]
      measurements$length_data[ind_ytypes1]<- sum(!is.na(observedi))
    }   
    names(measurements$length_data) <- unique(measurements$datatypes)
    main_struct$measurements[[ind_kreise]]<-measurements
    
    ###Daily:
    
    measurements_daily$namesmeasdata<-c("T","Y")
    measurements_daily$date<-c(cases_date_daily,#temp_measurements0$Meldedatum[ind_allcases],
                               temp_measurements0$Meldedatum[ind_deaths_daily],
                               dates_critical)
    measurements_daily$data <- rbind(as.numeric(measurements_daily$date-min(measurements_daily$date)),
                                     c(cases_daily,#temp_measurements0$AllConfCases[ind_allcases]
                                       temp_measurements0$NewDeaths[ind_deaths_daily],
                                       temp_measurements0$covid_inICU_upscaled[ind_critical]))
    measurements_daily$datatypes <- c(rep.int(x="Total",times=length(cases_daily)),#sum(ind_allcases)
                                      rep.int(x="Death",times=sum(ind_deaths_daily)),
                                      rep.int(x="Critical",times=sum(ind_critical)))
    measurements_daily$RelevClinData <- (dates_critical>"2020-04-19")
    
    main_struct$Daily$measurements[[ind_kreise]]<-measurements_daily
    
    
    main_struct$YTYPE[[ind_kreise]]<-list()
    main_struct$YTYPE[[ind_kreise]]$listnum<-unique(sort(alloutcomesnum) )#unique(sort(alloutcomesnum[index_id]) )!!!
    main_struct$YTYPE[[ind_kreise]]$list<-main_struct$outcomesnames[main_struct$YTYPE[[ind_kreise]]$listnum]
    main_struct$YTYPE[[ind_kreise]]$num<-length(main_struct$YTYPE[[ind_kreise]]$list)
    main_struct$YTYPE[[ind_kreise]]$size<-list()
    main_struct$YTYPE[[ind_kreise]]$index<-list()
    main_struct$observed[[ind_kreise]]<-list()
    main_struct$prierr_indiv[[ind_kreise]]<-list()
    main_struct$initname[[ind_kreise]]<-list()
    main_struct$initval[[ind_kreise]]<-vector(length = max(main_struct$YTYPE[[ind_kreise]]$listnum) )
    for(ind_ytype1 in 1:main_struct$YTYPES$num)
    {
      ind_ytype<-main_struct$YTYPES$listnum[ind_ytype1]
      ind_ytype_symb<-main_struct$YTYPES$list[ind_ytype1]
      main_struct$YTYPE[[ind_kreise]]$size[ind_ytype]<-sum(alloutcomesnum==main_struct$YTYPES$listnum[ind_ytype1]) #alloutcomesnum[index_id]
      main_struct$YTYPE[[ind_kreise]]$index[[ind_ytype]]<-(alloutcomesnum==main_struct$YTYPES$listnum[ind_ytype1])#alloutcomesnum[index_id]
      #index_id_y_ytype is unnecessary, use main_struct$YTYPE[[ind_]]$size[ind_ytype] instead!!!!!
      if(sum(main_struct$YTYPE[[ind_kreise]]$index[[ind_ytype]],na.rm = T)>0)
      {
        main_struct$observed[[ind_kreise]][[ind_ytype]]<-main_struct$measurements[[ind_kreise]]$data[2,][main_struct$YTYPE[[ind_kreise]]$index[[ind_ytype]]]#main_struct$Y[[ind_]] #alternative to matlab definition !!!! 
        
      }
    }
    #####
    #######
  }
  main_struct$indiv <- list()
  main_struct$indiv$ID <- 1:main_struct$numindiv
  main_struct$ID <- main_struct$covariates$ID
  main_struct$covariates$DateLock1 <- as.Date(main_struct$covariates$DateLock1)
  main_struct$covariates$DateLock2 <- as.Date(main_struct$covariates$DateLock2)
  main_struct$covariates$DateLock3 <- as.Date(main_struct$covariates$DateLock3)
  main_struct$covariates$DateStart <- as.Date(main_struct$covariates$DateStart)
  main_struct$covariates$Unlock <- as.Date(main_struct$covariates$Unlock)
  main_struct$covariates$Unlock1 <- as.Date(main_struct$covariates$Unlock1)
  main_struct$covariates$Unlock2 <- main_struct$covariates$Unlock1#initialization
  return(main_struct)
}
####
LoadIntegralModelCOVID<-function(path_model, main_struct)
{
  #setwd(path_model)
  #####
  temp_main<-read.csv(paste(path_model,'main.csv',sep=''), sep = ";",dec = ',')
  main_struct$numallpar<-dim(temp_main)[1]
  for(ind_names in names(temp_main))
  {
    #print(ind_names)
    main_struct[[ind_names]]<-temp_main[[ind_names]]
  }
  main_struct$names<-as.character(main_struct$names)
  #colnames(main_struct$init_mean) <- main_struct$names
  main_struct$togethertransform_label_LB_UB<-cbind(temp_main$transform_label,temp_main$LB,temp_main$UB)
  main_struct$poptransform_label_LB_UB <- main_struct$togethertransform_label_LB_UB[main_struct$subs_pop_mu_opt==1,]
  #!!! for every subject separately!main_struct$num_estim_indiv<-sum(main_struct$subs_indiv)
  temp_allpar<-1:main_struct$numallpar
  if('age_spec'%in%names(temp_main))
  {
    main_struct$label_poly_age <- 1
    main_struct$age_spec <- temp_main$age_spec*temp_main$subs_indiv
    main_struct$age_spec_param_num<- sum(temp_main$age_spec*temp_main$subs_indiv)
    if(main_struct$age_spec_param_num==0)
    {
      main_struct$label_poly_age <- 0
    }
    
  }else{
    main_struct$label_poly_age <- 0
  }
  if(main_struct$label_poly_age==1)
  {
    file_names <- list.files(path_model)
    num_files <- length(file_names)
    ind_age_speci <- 0
    ind_age_speci <- 0
    main_struct$age_spec_param<-list()
    main_struct$age_spec_phi<-list()
    main_struct$age_spec_param_names<-main_struct$names[main_struct$age_spec==1]
    main_struct$age_spec_id<-NULL
    for(indi in 1:num_files)
    {
      if(length(grep(pattern = 'agespecpar_',x=file_names[[indi]]))>0)
      {
        ind_age_speci<- ind_age_speci+1
        temp_age<-read.csv(paste(path_model,file_names[[indi]],sep=''), sep = ";",dec = ',')
        #print(indi)
        
        help_1 <-sub(pattern = 'agespecpar_',x=file_names[[indi]],replacement='')
        main_struct$age_spec_id<-c(main_struct$age_spec_id,as.numeric(sub(pattern = '.csv',x=help_1,replacement='')))
        main_struct$age_spec_dim <- c(length(temp_age$age),main_struct$age_spec_param_num)
        main_struct$age_spec_param[[as.numeric(sub(pattern = '.csv',x=help_1,replacement=''))]]<-array(dim=main_struct$age_spec_dim)
        main_struct$age_spec_phi[[as.numeric(sub(pattern = '.csv',x=help_1,replacement=''))]]<-array(dim=main_struct$age_spec_dim)
        for(ind_poly in 1:main_struct$age_spec_param_num)
        {  
          main_struct$age_spec_param[[as.numeric(sub(pattern = '.csv',x=help_1,replacement=''))]][,ind_poly]<- temp_age[[main_struct$age_spec_param_names[ind_poly]]]
          for(ind_agei in 1:length(temp_age$age))
          {
            main_struct$age_spec_phi[[as.numeric(sub(pattern = '.csv',x=help_1,replacement=''))]][ind_agei,ind_poly]<-RevTransformParameters(x=main_struct$age_spec_param[[as.numeric(sub(pattern = '.csv',x=help_1,replacement=''))]][ind_agei,ind_poly],
            main_struct,name=main_struct$age_spec_param_names[ind_poly])
          }
        }
        
      }
    }
    
    main_struct$subs_indiv_num<-temp_allpar[main_struct$subs_indiv==1]
  }else{
    main_struct$subs_indiv_num<-temp_allpar[main_struct$subs_indiv==1]
  }
  #main_struct$subs_indiv_num<-temp_allpar[main_struct$subs_indiv==1]
  
  main_struct$init_mu<-RevTransformParametersAll(main_struct$init_mean,main_struct)
  main_struct$priormuall<-RevTransformParametersAll(main_struct$priorval,main_struct)
  main_struct$ind_prior<-list()
  
  main_struct$init_om<-main_struct$init_sd
  ######
  #errors:
  if("error.csv"%in%list.files(path = (path_model)))
  {
    main_struct<-LoadErrorsModel(path_model, main_struct)
  }
  ########################
  main_struct$subs_indiv <- (main_struct$subs_indiv==1)#label of subset of parameters which are individualized, nonzero omega
  
  main_struct$indivoptions<- array(dim = c(main_struct$numindiv,main_struct$numallpar))
  for(ind_r in 1:main_struct$numindiv)
  {
    main_struct$indivoptions[ind_r,] <- main_struct$subs_indiv##the new assumption, every subject has the same estimated parameters!!!???
  }
  
  main_struct$subs_pop_opt <- main_struct$subs_pop_mu_opt==1;#label ofsubset of parameters, whose population value is optimized 
  main_struct$ind_opt_parameters <- (1:main_struct$numallpar)[main_struct$subs_pop_opt>0]
  main_struct$names_pop_opt_parameters <- main_struct$names[main_struct$ind_opt_parameters]
  if(sum(main_struct$subs_indiv*main_struct$subs_pop_mu_opt)>0)
  {
    'inconsistent estimated individual and population values, redundant population values are cancelled '
    temp_par <- 1:main_struct$numallpar
    main_struct$subs_pop_mu_opt[temp_par[main_struct$subs_indiv*main_struct$subs_pop_mu_opt>0]] <- 0
    main_struct$subs_pop_opt[temp_par[main_struct$subs_indiv*main_struct$subs_pop_mu_opt>0]] <- 0
  }#end
  main_struct$mucurr<-list()
  main_struct$mucurr$pop <- main_struct$init_mu[main_struct$subs_pop_opt>0]#
  main_struct$omcurr$pop <- main_struct$init_om[main_struct$subs_pop_opt>0]#allom(main_struct$subs_pop_mu_opt_gen);
  temp_priormu_pop <- main_struct$priormuall[main_struct$subs_pop_opt>0]
  temp_priorom_pop <- main_struct$priorom[main_struct$subs_pop_opt>0]
  #main_struct$prior1$pop=temp_prior_pop;
  temp_prior_pop <- main_struct$prior[main_struct$subs_pop_opt>0]
  main_struct$numprior <- list()
  main_struct$numprior$pop <- sum(temp_prior_pop)
  main_struct$ind_prior$pop <- main_struct$prior[main_struct$subs_pop_opt>0]
  if(main_struct$numprior$pop>0)
  {
    main_struct$xprior<-list()
    main_struct$omprior<-list()
    main_struct$xprior$pop <- temp_priormu_pop[temp_prior_pop==1]
    main_struct$omprior$pop=temp_priorom_pop[temp_prior_pop==1]
  }
  temp_allpar <- 1:main_struct$numallpar
  main_struct$subs_indiv_num <- temp_allpar[main_struct$subs_indiv==1]
  main_struct$subs_indivi<-array(dim=c(main_struct$numindiv,main_struct$numallpar))
  main_struct$subs_indivi_num <- list()
  main_struct$names_indiv_opt_parameters <- list()
  main_struct$num_estim_indiv<- vector(length = main_struct$numindiv )
  main_struct$indivtransform_label_LB_UB <- list()
  
  main_struct$subs_popi<- array(dim=c(main_struct$numindiv,main_struct$numallpar))
  main_struct$num_estim_pop <- vector(length = main_struct$numindiv )
  main_struct$poptransform_label_LB_UB <- vector(length = main_struct$numindiv )
  main_struct$subs_popi_num <- list()
  main_struct$mucurr$indiv <- list()
  main_struct$omcurr$indiv <- list()
  temp_priormu_indiv <- list() 
  temp_priorom_indiv <- list()
  
  main_struct$ind_prior$indiv <- list()
  main_struct$ind_prior$indivspec <- list()
  main_struct$numprior$indiv <- list()
  main_struct$avprior <- list()#obsolete?!!!
  main_struct$avprior$indiv <- list()
  main_struct$priornumeric <- list()
  main_struct$priornumeric$indiv  <- list()
  main_struct$subs_indivi<-array(dim=c(main_struct$numindiv,main_struct$numallpar))
  for (Row in 1:main_struct$numindiv)#size(main_struct$indivoptions,1)
  {
    ##this is redundant:
    main_struct$subs_indivi[Row,] <- as.vector(main_struct$indivoptions[Row,]*main_struct$subs_indiv)#!!!not *
    
    main_struct$num_estim_indiv[Row] <- sum(main_struct$subs_indivi[Row,])
    main_struct$subs_indivi_num[[Row]] <- temp_allpar[main_struct$subs_indivi[Row,]==1]
    main_struct$names_indiv_opt_parameters[[Row]] <- main_struct$names[main_struct$subs_indivi_num[[Row]]]
    main_struct$indivtransform_label_LB_UB[[Row]] <- main_struct$togethertransform_label_LB_UB[main_struct$subs_indivi[Row,]==1,]
    main_struct$subs_popi[Row,] <- main_struct$indivoptions[Row,]*main_struct$subs_pop_mu_opt
    main_struct$num_estim_pop[Row] <- sum(main_struct$subs_popi[Row,])
    if(sum(main_struct$subs_popi[Row,]==1)>0)
    {
      main_struct$subs_popi_num[[Row]] <- temp_allpar[main_struct$subs_popi[Row,]==1]
    }
    
    
    main_struct$mucurr$indiv[[Row]] <- main_struct$init_mu[main_struct$subs_indivi[Row,]>0]#
    main_struct$omcurr$indiv[[Row]] <- main_struct$init_om[main_struct$subs_indivi[Row,]>0]#allom(main_struct$subs_pop_mu_opt_gen)
    
    temp_priormu_indiv[[Row]] <- main_struct$priormuall[main_struct$subs_indivi[Row,]>0]
    temp_priorom_indiv[[Row]] <- main_struct$priorom[main_struct$subs_indivi[Row,]>0]
    main_struct$ind_prior$indiv[[Row]] <- main_struct$prior[main_struct$subs_indivi[Row,]>0]
    main_struct$ind_prior$indivspec[[Row]] <- main_struct$ind_prior$indiv[[Row]]*main_struct$subs_special[main_struct$subs_indivi[Row,]>0]
    
    main_struct$numprior$indiv[[Row]] <- sum(main_struct$ind_prior$indiv[[Row]])
    temp_ind <- 0
    main_struct$avprior$indiv[[Row]]<-vector(length=main_struct$num_estim_indiv[Row])#!!!!
    for (ind_indiv_par in 1:main_struct$num_estim_indiv[Row])
    {
      if(main_struct$num_estim_indiv[Row]>0)
      {
        if(main_struct$ind_prior$indiv[[Row]][ind_indiv_par]==1)
        {
          temp_ind <- temp_ind+1
          if(main_struct$ind_prior$indivspec[[Row]][ind_indiv_par]==1)
          {
            main_struct$avprior$indiv[[Row]][temp_ind] <- 1
          }
          else{
            main_struct$avprior$indiv[[Row]][temp_ind] <- 0
          }#end
          #main_struct$avprior.indiv[[Row]](temp_ind) <- main_struct$subs_indivi_num[[Row]][ind_indiv_par]
        }#end
      }
    }#end
    #!!!! per indiv!!!!
    
    
    if(main_struct$numprior$indiv[[Row]]>0)
    {
      main_struct$xprior$indiv[[Row]] <- temp_priormu_indiv[[Row]][main_struct$ind_prior$indiv[[Row]]==1]
      main_struct$omprior$indiv[[Row]] <- temp_priorom_indiv[[Row]][main_struct$ind_prior$indiv[[Row]]==1]  
      main_struct$priornumeric$indiv[[Row]] <- main_struct$subs_indivi_num[[Row]][main_struct$ind_prior$indiv[[Row]]==1]
    }
  }#end
  ###
  
  
  ###???end
  main_struct$numindivpar1 <- sum(main_struct$subs_indiv)
  main_struct$numpoppar <- sum(main_struct$subs_pop_opt)
  main_struct$subs_popisubs <- array(dim =c(dim(main_struct$subs_popi)[1],sum(main_struct$subs_pop_mu_opt)))
  main_struct$subs_popisubs[,] <- main_struct$subs_popi[,main_struct$subs_pop_mu_opt==1]
  #main_struct$relev_input_indiv=1:main_struct$numindiv#initalization
  #################
  
  ###indivpar:
  if("indivpar.csv"%in%list.files(path = (path_model)))
  {  
    temp_indiv<-read.csv(paste(path_model,'indivpar.csv',sep = ''), sep = ";",dec = ',')
    num_cols <- dim(temp_indiv)[2]-1
    temp_indiv<- temp_indiv[,2:dim(temp_indiv)[2]]
  }else{
    temp_indiv<- array(dim = c(main_struct$numindiv,main_struct$numallpar))
    num_cols <- dim(temp_indiv)[2]
    for(ind_r in 1:main_struct$numindiv)
    {
      temp_indiv[ind_r,] <- main_struct$init_mean
    }
  }
  if("indivconst.csv"%in%list.files(path = (path_model)))
  {  
    temp_indiv_const<-read.csv(paste(path_model,'indivconst.csv',sep = ''), sep = ";",dec = ',')
    num_cols_const <- dim(temp_indiv_const)[2]-1
    temp_indiv_const<- temp_indiv_const[,2:dim(temp_indiv_const)[2]]
    main_struct$indiv$const<-temp_indiv_const
  }
  if("DatesUncertainty.csv"%in%list.files(path = (path_model)))
  {
    ##treatments
    temp_date<-read.csv(paste(path_model,'DatesUncertainty.csv',sep = ''), sep = ";",dec = ',')
    main_struct$dates_treatments <- temp_date
  }
  if("resid.csv"%in%list.files(path = (path_model)))
  {  
    temp_resid<-read.csv(paste(path_model,'resid.csv', sep = ''), sep = ";",dec = ',')
  }
  main_struct$numtypesresidinit<-dim(temp_resid)[2]-1#size( ddataresid,2)-1;
  main_struct$typesresidinit<-colnames(temp_resid)[2:(main_struct$numtypesresidinit+1)]#=dtextresid(2:end);###dnumresid(1,2:end);
  main_struct$indiv$psi<-array(dim=c(dim(temp_indiv)[1],num_cols))
  main_struct$indiv$phi<-array(dim=c(dim(temp_indiv)[1],num_cols))
  colnames(main_struct$indiv$psi) <- main_struct$names
  colnames(main_struct$indiv$phi) <- main_struct$names
  main_struct$indiv$phiopt<-list()
  main_struct$subs_indivfixed <-array(dim=c(dim(temp_indiv)[1],num_cols))
  main_struct$num_estim_indivfixed <- vector(length = dim(temp_indiv)[1])###Corrected!!!! 20.11.2018
  main_struct$mucurr$indivfixed <-list()
  main_struct$omcurr$indivfixed <- list()
  temp_priormu_indivfixed <- list()
  temp_priorom_indivfixed <- list()
  main_struct$ind_prior$indivfixed <- list()
  main_struct$numprior$indivfixed <- list()
  main_struct$xprior$indivfixed <- list()
  main_struct$omprior$indivfixed  <- list()
  main_struct$residinit <- array(dim=c(dim(temp_resid)[1],sum(main_struct$act_err)))
  main_struct$allresidinit <- array(dim=c(dim(temp_resid)[1],dim(temp_resid)[2]-1))
  main_struct$allresidnum<-dim(temp_resid)[2]-1
  colnames(main_struct$allresidinit)<-colnames(temp_resid[2:dim(temp_resid)[2]])
  main_struct$numindivpar2 <- vector(length=main_struct$numindiv) 
  #########################
  for (Row in 1:main_struct$numindiv)#numnotnanind#size(dnum,1)
  {
    #   Row
    ###17.04.2019:!!!??? May be unnecessary!!!???
    if(main_struct$numtypesresidinit>0)
    {
      tempresidi<- as.matrix(temp_resid)[Row,2:(dim(temp_resid)[2])]#tempresid=dnumresid(Row,2:end);
      if(sum(names(main_struct)=='act_err')>0)#if(sum(strcmp(fields(main_struct),'act_err'))>0)
      {
        main_struct$residinit[Row,] <- tempresidi[main_struct$act_err==1]
      }else
      {
        main_struct$residinit[Row,] <- tempresidi[main_struct$opt_err==1]#seqrelev[Row]
      }
      main_struct$allresidinit[Row,] <- tempresidi
    }
    
    for (Col in 1:(num_cols))#(Col in 1:(dim(temp_indiv)[2]-1)) !!!!
    {
                #print(c(Row,Col))  
      main_struct$indiv$psi[Row,Col] <- temp_indiv[Row,Col]#!!!+1#eval('main_struct$indiv$psi(Row,Col)= ddataindiv{Row+1,Col+1};');#seqrelev[Row] dnum(seqrelev[Row],Col+1)ddata{Row,Col+1}
    }
    
    main_struct$indiv$phi[Row,] <- RevTransformParametersAll(main_struct$indiv$psi[Row,],main_struct);
    #main_struct$indiv$phiopt.pop[[Row]]=main_struct$indiv$phi(Row,main_struct$subs_pop_mu_opti[Row,]==1);
    
    main_struct$indiv$phiopt[[Row]]<-array(dim=c(dim(temp_resid)[1],main_struct$numindivpar1))
    main_struct$indiv$phiopt[[Row]] <- main_struct$indiv$phi[Row,main_struct$subs_indivi[Row,]>0]
    if(main_struct$label_poly_age>0)
    {
      if(Row%in%main_struct$age_spec_id)
      {
        main_struct$indiv$phiopt[[Row]][main_struct$names_indiv_opt_parameters[[Row]]%in%main_struct$age_spec_param_names] <- main_struct$age_spec_phi[[Row]][1,]
        main_struct$indiv$phi[Row,main_struct$subs_indivi[Row,]>0][main_struct$names_indiv_opt_parameters[[Row]]%in%main_struct$age_spec_param_names] <- main_struct$age_spec_phi[[Row]][1,]
        main_struct$indiv$psi[Row,main_struct$subs_indivi[Row,]>0][main_struct$names_indiv_opt_parameters[[Row]]%in%main_struct$age_spec_param_names] <- main_struct$age_spec_param[[Row]][1,]
        
        main_struct$indiv$phi
        for(ind_agei in 2:dim(main_struct$age_spec_phi[[Row]])[1])
        {
          main_struct$indiv$phiopt[[Row]] <- c(main_struct$indiv$phiopt[[Row]],main_struct$age_spec_phi[[Row]][ind_agei,])
        }
      }
    }
    main_struct$numindivpar2[Row]<-length(main_struct$indiv$phiopt[[Row]] )
    #fixed individual do not participate in prior!!!! Later to fix
    main_struct$subs_indivfixed[Row,] <- ((!main_struct$indivoptions[Row,])*main_struct$subs_indiv)
    main_struct$num_estim_indivfixed[Row] <- sum(main_struct$subs_indivfixed[Row,]);
    if(main_struct$num_estim_indivfixed[Row]>0)
    {
      main_struct$mucurr$indivfixed[[Row]] <- main_struct$indiv$phi[Row,main_struct$subs_indivfixed[Row,]>0]#
      main_struct$omcurr$indivfixed[[Row]] <- main_struct$init_om[main_struct$subs_indivfixed[Row,]>0]#allom(main_struct$subs_pop_mu_opt_gen);
      temp_priormu_indivfixed[[Row]] <- main_struct$priormuall[main_struct$subs_indivfixed[Row,]>0]
      temp_priorom_indivfixed[[Row]] <- main_struct$priorom[main_struct$subs_indivfixed[Row,]>0]
      
      #temp.prior.indiv[[Row]]=main_struct$prior(main_struct$subs_indivi[Row,]>0);
      main_struct$ind_prior$indivfixed[[Row]] <- main_struct$prior[main_struct$subs_indivfixed[Row,]>0]
      main_struct$numprior$indivfixed[[Row]] <- sum(main_struct$ind_prior$indivfixed[[Row]]);
      if(main_struct$numprior$indivfixed[[Row]]>0)
      {
        main_struct$xprior$indivfixed[[Row]] <- temp_priormu_indivfixed[[Row]][main_struct$ind_prior$indivfixed[[Row]]==1]
        main_struct$omprior$indivfixed[[Row]] <- temp_priorom_indivfixed[[Row]][main_struct$ind_prior$indivfixed[[Row]]==1]
      }#end;
    }#end
    ###
  }#end
  ###############################################
  
  temp_specialpar<-read.csv(paste(path_model,'specialpar.csv',sep = ''), sep = ";",dec = ',')
  main_struct$specialpar<-list()
  main_struct$specialpar$names<-temp_specialpar[,1]
  main_struct$specialpar$val<-temp_specialpar[,2]
  main_struct$numspecialpar<-length(main_struct$specialpar$val)
  
  if("ModelClass.csv"%in%list.files(path = (path_model)))
  {
    temp_modelproperties<-read.csv(paste(path_model,'ModelClass.csv',sep = ''), sep = ";",dec = ',')
    main_struct$modelproperties<-temp_modelproperties
  }
  main_struct
}

###################
####General project-non-specific functions
###From LoadFunctions:
AssignCompartmentsSimple<-function(ind_indiv,path_model,main_struct)#,path_data
{
  main_struct$compartments<-list()
  #temp_measurements <- read.csv(paste(path_data,'measurements.csv',sep=''), sep = ";",dec = ',')
  temp_compartments<-read.csv(paste(path_model,'compartments.csv',sep=''), sep = ";",dec = ',')
  NumOfRows<-dim(temp_compartments)[1]
  basic_vect<-as.character(temp_compartments$type)=='basic'
  treat_vect<-as.character(temp_compartments$type)=='treatment'
  all_id <- unique(main_struct$covariates$ID)
  main_struct$ProtocolNew <- list()
  
  for(i0 in 1: main_struct$numindiv)#length(main_struct$ID))
  {
    i<-main_struct$ID[i0]#i <- all_id[i0]# main_struct$ID[i0]
    if(i %in% all_id)#for the case when certain patient has no treatment
    {  
      #index_id<-temp_measurements$ID==i
      #index_id_meas<-temp_measurements$ID==i
      main_struct$compartments[[i0]]<-list()
      # print(i0)
      temp_comp <- basic_vect
      # print(temp_comp)
      #temp_treat_list <- unique(as.character(temp_treatments$CMT)[index_id])#unique(main_struct$Protocol[[i0]]$cmt1)#main_struct$Protocol$cmt1[[i0]]
      
      NumOfCompartmentsi <- sum(temp_comp)
      main_struct$compartments[[i0]]$NumofVariables<-sum(temp_compartments$num_var[temp_comp==1])
      main_struct$compartments[[i0]]$check<-0*array(dim=c(main_struct$compartments[[i0]]$NumofVariables,1))
      main_struct$compartments[[i0]]$num_ind<-list()
      main_struct$compartments[[i0]]$compind<-array(dim = c(NumOfCompartmentsi,max(temp_compartments$num_var)))
      main_struct$compartments[[i0]]$compnum<-vector(length = NumOfCompartmentsi)
      main_struct$compartments[[i0]]$names<-vector(length = NumOfCompartmentsi)
      ind_ <- 0
      ind_i<-0
      #print(NumOfRows)
      for (Row in 1:NumOfRows)#(Row in 2:(NumOfRows+1))
      {
        if(temp_comp[Row]==1)#Row-1!!!!
        {
          ind_i<-ind_i+1
          main_struct$compartments[[i0]]$names[ind_i]=as.character(temp_compartments$names[Row])#datacompartments{Row,1};not main_struct$compartments[[i0]]$names[Row]!!!!
          main_struct$compartments[[i0]][[paste(as.character(temp_compartments$names[Row]),'num',sep="")]]<-temp_compartments$num_var[Row]
          main_struct$compartments[[i0]][[paste(as.character(temp_compartments$names[Row]),'_ind',sep="")]]<-(ind_+1):(ind_+temp_compartments$num_var[Row])
          main_struct$compartments[[i0]]$compind[ind_i,1:temp_compartments$num_var[Row]]<-(ind_+1):(ind_+temp_compartments$num_var[Row])
          main_struct$compartments[[i0]]$compnum[ind_i]<-temp_compartments$num_var[Row]
          #print(paste(as.character(temp_compartments$names[Row]),'num',sep=""))
          # print(names(main_struct$compartments[[i0]]))
          #eval(['main_struct$compartments[[i0]]$' datacompartments{Row,1} 'num= datacompartments{Row,2};']);
          
          #    eval(['main_struct$compartments[[i0]]$' datacompartments{Row,1}  '_ind= [ind_+1:(ind_+datacompartments{Row,2})];']);
          main_struct$compartments[[i0]]$num_ind[[Row]]<-(ind_+1):(ind_+temp_compartments$num_var[Row])
          if(Row==1){
            main_struct$compartments[[i0]]$check<-rep.int(x=1,times=temp_compartments$num_var[Row])
          }else
          {
            main_struct$compartments[[i0]]$check<-c(main_struct$compartments[[i0]]$check,rep.int(x=1,times=temp_compartments$num_var[Row]))
          }
          #main_struct$compartments[[i0]]$check(ind_+1:(ind_+datacompartments{Row,2}))=1;
          #    main_struct$compartments[[i0]]$num_ind{Row-1}= [ind_+1:(ind_+datacompartments{Row,2})];
          ind_<-ind_+temp_compartments$num_var[Row]#datacompartments{Row,2};
          main_struct$compartments[[i0]]$NumOfCompartments=NumOfCompartmentsi
        }
        
      }
      
      # main_struct$compartmentsall$structure=datacompartments;
      main_struct$compartments[[i0]]$baddefined=sum(main_struct$compartments[[i0]]$check==0);
      
    } 
  }
  ###########
  main_struct$colnames_sim_res <- CreateNameOutput(ind_indiv,main_struct)
  return(main_struct)
}
#########################################
#########################################
UpdateStructures<-function(main_struct)
{
  main_struct$num_totobstype<-rep.int(x=0, times=max(main_struct$YTYPES$listnum))
  main_struct$num_totobstypew<-rep.int(x=0, times=max(main_struct$YTYPES$listnum))
  main_struct$num_sepobstype<-array(dim = c(main_struct$numindiv,max(main_struct$YTYPES$listnum)))
  main_struct$num_sepobstypew<-array(dim = c(main_struct$numindiv,max(main_struct$YTYPES$listnum)))
  for( i in 1: main_struct$numindiv)
  {
    for (ind_ytype0 in 1:main_struct$YTYPES$num)
    {
      ind_ytype=main_struct$YTYPES$listnum[ind_ytype0]
      
      if(main_struct$YTYPE[[i]]$size[[ind_ytype]]>0)
      { 
        main_struct$num_totobstype[ind_ytype] <- main_struct$num_totobstype[ind_ytype]+length(main_struct$observed[[i]][[ind_ytype]])
        main_struct$num_totobstypew[ind_ytype] <- main_struct$num_totobstypew[ind_ytype]+main_struct$weights[i,ind_ytype]*length(main_struct$observed[[i]][[ind_ytype]])
        #}!!!! in R ?????
        #print(c(i,ind_ytype))
        main_struct$num_sepobstype[i,ind_ytype] <- length(main_struct$observed[[i]][[ind_ytype]])
        main_struct$num_sepobstypew[i,ind_ytype] <- main_struct$weights[i,ind_ytype]*length(main_struct$observed[[i]][[ind_ytype]])
      }
    }
    #print(i)
  }
  ####
  main_struct
}
#############################
#############################
GeneratePhiopt<-function(main_struct)
{
  phiopt<-list()
  phiopt$pop<-main_struct$mucurr$pop#=main_struct.mucurr.pop;
  phiopt$indiv<-main_struct$indiv$phiopt
  phiopt$resid <- main_struct$residinit
  phiopt
}
LoadErrorsModel<-function(path_model, main_struct)
{
  #setwd(path_model)
  temp_error<-read.csv(paste(path_model,'error.csv',sep = ''), sep = ";",dec = ',')
  main_struct$acurr <- as.numeric(temp_error[1,2:dim(temp_error)[2]])
  main_struct$opt_err <- as.numeric(temp_error[2,2:dim(temp_error)[2]])
  main_struct$logres <- as.numeric(temp_error[6,2:dim(temp_error)[2]])
  main_struct$prierre <- as.numeric(temp_error[7,2:dim(temp_error)[2]])
  main_struct$transformres <- as.numeric(temp_error[3,2:dim(temp_error)[2]])
  main_struct$expdist <- as.numeric(temp_error[4,2:dim(temp_error)[2]])
  main_struct$scalingini <- as.numeric(temp_error[5,2:dim(temp_error)[2]])
  main_struct$act_err <- as.numeric(temp_error[8,2:dim(temp_error)[2]] )
  main_struct$act_err_num<-(1:(dim(temp_error)[2]-1))[main_struct$act_err==1]# the field act_err exists always!!!!
  if(dim(temp_error)[1]>8)#now is a default!!!
  {
    main_struct$expsens<-as.numeric(temp_error[9,2:dim(temp_error)[2]])
    if(dim(temp_error)[1]>9)
    {  
      #main_struct$grade_resid_err<-as.numeric(temp_error[10,2:dim(temp_error)[2]])
      #if(dim(temp_error)[1]>10)#now is a default!!!
      #{
        main_struct$cumulexpsens<-as.numeric(temp_error[10,2:dim(temp_error)[2]])
        main_struct$grade_resid_err<-rep.int(x = 0,dim(temp_error)[2])
      #}
    }else
    {
      main_struct$grade_resid_err<-rep.int(x = 0,dim(temp_error)[2])#dim(temp=temp_error)[2]
    }
  }else{
    main_struct$expsens<-rep.int(x = 0,dim(temp_error)[2])
    main_struct$grade_resid_err<-rep.int(x = 0,dim(temp_error)[2])
  }
  main_struct$error_structure <- temp_error
  ##the follows are pre-defined in the beginning
  #main_struct$outcomesnames<-colnames(temp_error)[2:dim(temp_error)[2]]
  #main_struct$numoutcomesnames=dim(temp_error)[2]-1
  #main_struct$outcomesnumlist=1:main_struct$numoutcomesnames
  ####output:
  main_struct
}
ChangeDec<-function(path_dir, main_struct)
{
  #dir1 <-"R:/Blutmodelle/Yuri/daten/Ulrike_Klotz/Aufgabe 1/Morstyn/Bilderdaten/4 mit SE/"
  #dir.create(paste(dir1,"/corrected/",sep = ""))
  #dir2<-paste(dir1,"/corrected/",sep = "")
  file_names <- list.files(path_dir)
  num_files <- length(file_names)
  for(indi in 1:num_files)
  {
    if(length(grep(pattern = 'csv',x=file_names[[indi]]))>0)
    {
      #print(file_names[[indi]])
      temp_dat<-read.csv(paste(path_dir,'/',file_names[[indi]], sep = ""),sep=";",dec=main_struct$decold)
      if(file_names[[indi]]=='treatments.csv')
      {
        temp_dat<-read.csv(paste(path_dir,'/','treatments.csv', sep = ""),sep=";",dec=main_struct1$decold,colClasses = c("integer","numeric","numeric","character","character","character","integer","integer","integer","numeric","character","numeric"))
        temp_dat$TINF<-as.numeric(temp_dat$TINF)
      }
      write.table(x=temp_dat,sep = ";" ,  dec = main_struct$dec, row.names=F,col.names=T,file=paste(path_dir,'/',file_names[[indi]],sep = ''),na ='.')
    }
  }
}


ChangeDecSimp<-function(path_dir, decnew, decold)
{
  #dir1 <-"R:/Blutmodelle/Yuri/daten/Ulrike_Klotz/Aufgabe 1/Morstyn/Bilderdaten/4 mit SE/"
  #dir.create(paste(dir1,"/corrected/",sep = ""))
  #dir2<-paste(dir1,"/corrected/",sep = "")
  file_names <- list.files(path_dir)
  num_files <- length(file_names)
  for(indi in 1:num_files)
  {
    if(length(grep(pattern = 'csv',x=file_names[[indi]]))>0)
    {
      #print(file_names[[indi]])
      temp_dat<-read.csv(paste(path_dir,'/',file_names[[indi]], sep = ""),sep=";",dec=decold)
      if(file_names[[indi]]=='treatments.csv')
      {
        temp_dat<-read.csv(paste(path_dir,'/','treatments.csv', sep = ""),sep=";",dec=decold,colClasses = c("integer","numeric","numeric","character","character","character","integer","integer","integer","numeric","character","numeric"))
        temp_dat$TINF<-as.numeric(temp_dat$TINF)
      }
      write.table(x=temp_dat,sep = ";" ,  dec = decnew, row.names=F,col.names=T,file=paste(path_dir,'/',file_names[[indi]],sep = ''),na ='.')
    }
  }
}

ChangeDecSimpPerfile<-function(path_dir, file_name,decnew, decold)
{
  #dir1 <-"R:/Blutmodelle/Yuri/daten/Ulrike_Klotz/Aufgabe 1/Morstyn/Bilderdaten/4 mit SE/"
  #dir.create(paste(dir1,"/corrected/",sep = ""))
  #dir2<-paste(dir1,"/corrected/",sep = "")
  #file_names <- list.files(path_dir)
  #num_files <- length(file_names)
  
  
  #print(file_names[[indi]])
  temp_dat<-read.csv(paste(path_dir,'/',file_name, sep = ""),sep=";",dec=decold)
  if(file_name=='treatments.csv')
  {
    temp_dat<-read.csv(paste(path_dir,'/','treatments.csv', sep = ""),sep=";",dec=decold,colClasses = c("integer","numeric","numeric","character","character","character","integer","integer","integer","numeric","character","numeric"))
    temp_dat$TINF<-as.numeric(temp_dat$TINF)
  }
  write.table(x=temp_dat,sep = ";" ,  dec = decnew, row.names=F,col.names=T,file=paste(path_dir,'/',file_name,sep = ''),na ='.')
  
  
}
#####
LoadMCMCParam <- function(main_struct,ind_indiv)
{
  
  nonmem_options<-list()
  nonmem_options$m<-5000#10000#2000
  nonmem_options$mem <- 400# 1000#400
  nonmem_options$ignore<- 0#20
  nonmem_options$acceptl<-0.2
  nonmem_options$acceptu<-0.4
  nonmem_options$temp_check<-10
  nonmem_options$minvar <-0.00001#0.00001  #0.001
  nonmem_options$write <-0#1
  nonmem_options$dir_name <-'MCMC'
  
  nonmem_options$jump<-2
  nonmem_options$step_init<-10
  nonmem_options$num_indiv_runs <- 0#20
  nonmem_options$num_steps<-list()
  nonmem_options$num_steps$indiv<- 10
  nonmem_options$nLL_indiv <- 0#tot_goalf0$nLL_indiv
  #nonmem_options$prevres <-res_phi_arr[subseti,]
  nonmem_options$num_indiv_runs <-0
  nonmem_options$sd0<- (2.4)^2
  nonmem_options$eps_chol <- 0.00001#0.000001#0.00001# 0.000000001
  #nonmem_options$proposed_var<- diag(c(0.001,0.01,rep.int(x=0.01,times=main_struct$numpoppar-2)))#diag(0.001,c(main_struct$numpoppar,main_struct$numpoppar))
  if(main_struct$numpoppar>0)
  {
     nonmem_options$proposed_var<- diag(rep.int(x=0.01,times=main_struct$numpoppar))#
  }else{
    nonmem_options$proposed_var<- diag(rep.int(x=0.01,times=main_struct$num_estim_indiv[ind_indiv]))#
  } 
  ####
  nonmem_options
}


UploadAgeStructVariant <- function(ind_indiv,root_path,main_struct)
{  
  age_spec_incend<-read.csv(paste(root_path,'/data/','200709_age_specific_weekly_incidences_deaths.txt',sep =''), sep = "\t",dec = ',',header =T)
  age_spec_prop<-read.csv(paste(root_path,'/data/','200709_death_rate_week16.txt',sep =''), sep = "\t",dec = ',',header =T)
  age_spec_ICE<-read.csv(paste(root_path,'/data/','guess2_icu_age_distrib_germany_2020-07-14.txt',sep =''), sep = "\t",dec = ',',header =T)
  #guess_icu_age_distrib_germany_2020-07-13.txt'
  main_struct$age_struct<-list()
  
  main_struct$age_struct$ages<-levels(age_spec_incend$Altersgruppe)[1:6]
  
  main_struct$age_struct$TotalDeathRate<-as.numeric(as.character(age_spec_prop$TotalDeathRate))[age_spec_prop$Altersgruppe%in%main_struct$age_struct$ages]
  str(as.Date(age_spec_incend$Meldedatum_woche))
  age_spec_incend$Meldedatum_woche <- (as.Date(age_spec_incend$Meldedatum_woche))
  
  all_times<-unique(main_struct$measurements[[ind_indiv]]$date)
  xx_comp<-CompleteDates(all_times)
  #xx_comp%in%unique(sort(age_spec_incend$Meldedatum_woche))
  #xx_comp[!xx_comp%in%unique(sort(age_spec_incend$Meldedatum_woche))]
  num_weeks<-max(age_spec_incend$woche)                                  
  relev_weeks <- unique(sort(age_spec_incend$Meldedatum_woche))[unique(sort(age_spec_incend$Meldedatum_woche))>=(min(xx_comp)-7)]
  num_weeks<-length(relev_weeks)
  main_struct$age_struct$New_cases_week_distr <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  main_struct$age_struct$Death_cases_week_distr <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  main_struct$age_struct$Empir_Death_rate_week <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  main_struct$age_struct$Empir_Death_rate_week_del <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  
  colnames(main_struct$age_struct$New_cases_week_distr)<- c(main_struct$age_struct$ages,'Num Cases','ExpecDeathRate')
  colnames(main_struct$age_struct$Death_cases_week_distr)<- c(main_struct$age_struct$ages,'Num Cases','ExpecDeathRate')
  colnames(main_struct$age_struct$Empir_Death_rate_week)<- c(main_struct$age_struct$ages,'Mean Death','ExpecDeathRate')
  colnames(main_struct$age_struct$Empir_Death_rate_week_del)<- c(main_struct$age_struct$ages,'Mean Death','Mean Death1')#'ExpecCritRate'
  death_delay<-1#1#2
  ind_day_good<-16
  for(ind_w in 1:num_weeks)
  {
    relev_ind<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
    num_all_casesi <- sum(age_spec_incend$NewConfCases[ relev_ind])
    num_all_deathi <- sum(age_spec_incend$NewDeaths[ relev_ind])
    main_struct$age_struct$New_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+1] <- num_all_casesi
    main_struct$age_struct$Death_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+1] <- num_all_deathi
    for(ind_a in 1:length(main_struct$age_struct$ages))
    {
      relev_indi <- relev_ind&(age_spec_incend$Altersgruppe==main_struct$age_struct$ages[ind_a])
      main_struct$age_struct$Empir_Death_rate_week[ind_w,ind_a] <- age_spec_incend$AllDeaths[ relev_indi]/age_spec_incend$AllConfCases[ relev_indi]  
    }
    #main_struct$age_struct$Empir_Death_rate_week[ind_w,1:length(main_struct$age_struct$ages)] <- age_spec_incend$AllDeaths[ relev_ind]/age_spec_incend$AllConfCases[ relev_ind]
    main_struct$age_struct$Empir_Death_rate_week[ind_w,1+length(main_struct$age_struct$ages)] <- sum(age_spec_incend$AllDeaths[ relev_ind])/sum(age_spec_incend$AllConfCases[ relev_ind])
    if(ind_w>death_delay)
    {
      relev_ind_prev<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w-death_delay])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
      main_struct$age_struct$Empir_Death_rate_week_del[ind_w,1:length(main_struct$age_struct$ages)] <- age_spec_incend$AllDeaths[ relev_ind]/age_spec_incend$AllConfCases[ relev_ind_prev]
      main_struct$age_struct$Empir_Death_rate_week_del[ind_w,1+length(main_struct$age_struct$ages)] <- sum(age_spec_incend$AllDeaths[ relev_ind])/sum(age_spec_incend$AllConfCases[ relev_ind_prev])
    }
  }
  relev_ind1<- ((age_spec_incend$Meldedatum_woche==relev_weeks[1])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
  #num_all_casesi1 <- sum(age_spec_incend$NewConfCases[ relev_ind1])
  #num_all_deathi1 <- sum(age_spec_incend$NewDeaths[ relev_ind1])
  for(ind_w in 1:num_weeks)
  {  
    #if(ind_w<3)
    #{
    #  relev_ind<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
    #}else{
      relev_ind<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
    #}
    
    num_all_casesi <- sum(age_spec_incend$NewConfCases[ relev_ind])
    num_all_deathi <- sum(age_spec_incend$NewDeaths[ relev_ind])
    for(ind_a in 1:length(main_struct$age_struct$ages))
    {
      relev_indi <- relev_ind&(age_spec_incend$Altersgruppe==main_struct$age_struct$ages[ind_a])
      #relev_indi1 <- relev_ind1&(age_spec_incend$Altersgruppe==main_struct$age_struct$ages[ind_a])
      if(num_all_casesi>0)
      {
        
        main_struct$age_struct$New_cases_week_distr[ind_w,ind_a]<- age_spec_incend$NewConfCases[relev_indi]/num_all_casesi
        #if(ind_w==2)
        #{
        #  main_struct$age_struct$New_cases_week_distr[ind_w,ind_a]<- (age_spec_incend$NewConfCases[relev_indi]+age_spec_incend$NewConfCases[relev_indi1])/(num_all_casesi+num_all_casesi1)
         # main_struct$age_struct$New_cases_week_distr[1,ind_a]<-main_struct$age_struct$New_cases_week_distr[ind_w,ind_a]
        #}
      }else{
        main_struct$age_struct$New_cases_week_distr[ind_w,ind_a]<- 0
      }
      if(num_all_deathi>0)
      {
        main_struct$age_struct$Death_cases_week_distr[ind_w,ind_a]<- age_spec_incend$NewDeaths[relev_indi]/num_all_deathi
      }else{
        main_struct$age_struct$Death_cases_week_distr[ind_w,ind_a]<- 0
      }
    }
  } 
  #Correct for stability:
 # for(ind_w in 1:7)
 # {  
 #   main_struct$age_struct$Death_cases_week_distr[ind_w,1:6]<-main_struct$age_struct$Death_cases_week_distr[8,1:6]
 #   main_struct$age_struct$New_cases_week_distr[ind_w,1:6] <-main_struct$age_struct$New_cases_week_distr[8,1:6] 
 # }
  for(ind_w in 1:num_weeks)
  {
    #main_struct$age_struct$New_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$New_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(as.character(age_spec_prop$TotalDeathRate[age_spec_prop$Altersgruppe%in%main_struct$age_struct$ages])))
    #main_struct$age_struct$Death_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$Death_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(as.character(age_spec_prop$TotalDeathRate[age_spec_prop$Altersgruppe%in%main_struct$age_struct$ages])))
    main_struct$age_struct$New_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$New_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(main_struct$age_struct$Empir_Death_rate_week_del[ind_day_good,1:length(main_struct$age_struct$ages)]))
    main_struct$age_struct$Death_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$Death_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(main_struct$age_struct$Empir_Death_rate_week_del[ind_day_good,1:length(main_struct$age_struct$ages)]))
    # older for ICU:main_struct$age_struct$Empir_Death_rate_week_del[ind_w,2+length(main_struct$age_struct$ages)] <- sum(main_struct$age_struct$New_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(as.character(age_spec_ICE$icu_guessed)))
    if(ind_w>(death_delay))#+2
    {
      #relev_ind_prev<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w-death_delay])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
      main_struct$age_struct$Empir_Death_rate_week_del[ind_w,2+length(main_struct$age_struct$ages)] <- sum(main_struct$age_struct$Empir_Death_rate_week[16,1:length(main_struct$age_struct$ages)]*main_struct$age_struct$New_cases_week_distr[ind_w-death_delay,1:length(main_struct$age_struct$ages)])#*main_struct$age_struct$New_cases_week_distr[ind_w-death_delay,length(main_struct$age_struct$ages)+1])
    }
  }
  
  apply(main_struct$age_struct$New_cases_week_distr[,1:6],1,sum)
  main_struct$age_struct$death_exp_w <-main_struct$age_struct$New_cases_week_distr[,length(main_struct$age_struct$ages)+2]
  main_struct$age_struct$death_dir_w <-main_struct$age_struct$Death_cases_week_distr[,length(main_struct$age_struct$ages)+2]
  main_struct$age_struct$relev_weeks <- relev_weeks
  
  death_exp_d <- vector(length = length(xx_comp))
  #critical_exp_d <- vector(length = length(xx_comp))
  for(ind_d in 1:length(xx_comp))
  {
    #relev_ind <- max(max(which(xx_comp[ind_d]>relev_weeks)-death_delay),1)#length(relev_weeks)
    #death_exp_d[ind_d]<-main_struct$age_struct$death_exp_w[relev_ind]
    relev_ind <- max(max(which(xx_comp[ind_d]>relev_weeks)),1)#length(relev_weeks)
    death_exp_d[ind_d]<-main_struct$age_struct$Empir_Death_rate_week_del[,1+length(main_struct$age_struct$ages)][relev_ind] #1+length(main_struct$age_struct$ages
    #critical_exp_d[ind_d] <-main_struct$age_struct$Empir_Death_rate_week_del[,2+length(main_struct$age_struct$ages)][relev_ind] 
    
  }
  help_v <- main_struct$age_struct$Empir_Death_rate_week_del[,1+length(main_struct$age_struct$ages)]
  ind_interp<-max(which(is.na(help_v)|help_v==0))+1
  death_exp_d[is.na(death_exp_d)|(death_exp_d==0)]<- main_struct$age_struct$Empir_Death_rate_week_del[ind_interp,1+length(main_struct$age_struct$ages)]#0#NANcase
  ##check:
  #par(mfrow=c(1,1)) 
  #plot(all_times,death_exp_d)
  #lines(relev_weeks,main_struct$age_struct$death_exp_w)
  main_struct$age_struct$death_exp_d <- death_exp_d
  main_struct$age_struct$norm_death_exp_d <- main_struct$age_struct$death_exp_d/mean(main_struct$age_struct$death_exp_d )
  
  #main_struct$age_struct$critical_exp_d <- critical_exp_d
  main_struct$age_struct$norm_critical_exp_d <- main_struct$age_struct$critical_exp_d/mean(main_struct$age_struct$critical_exp_d )
  
  return(main_struct)
}


UploadAgeStruct <- function(ind_indiv,root_path,main_struct)
{  
  age_spec_incend<-read.csv(paste(root_path,'/data/','200709_age_specific_weekly_incidences_deaths.txt',sep =''), sep = "\t",dec = ',',header =T)
  age_spec_prop<-read.csv(paste(root_path,'/data/','200709_death_rate_week16.txt',sep =''), sep = "\t",dec = ',',header =T)
  age_spec_ICE<-read.csv(paste(root_path,'/data/','guess2_icu_age_distrib_germany_2020-07-14.txt',sep =''), sep = "\t",dec = ',',header =T)
  #guess_icu_age_distrib_germany_2020-07-13.txt'
  main_struct$age_struct<-list()
  
  main_struct$age_struct$ages<-levels(age_spec_incend$Altersgruppe)[1:6]
  
  main_struct$age_struct$TotalDeathRate<-as.numeric(as.character(age_spec_prop$TotalDeathRate))[age_spec_prop$Altersgruppe%in%main_struct$age_struct$ages]
  str(as.Date(age_spec_incend$Meldedatum_woche))
  age_spec_incend$Meldedatum_woche <- (as.Date(age_spec_incend$Meldedatum_woche))
  
  all_times<-unique(main_struct$measurements[[ind_indiv]]$date)
  xx_comp<-CompleteDates(all_times)
  
  #xx_comp%in%unique(sort(age_spec_incend$Meldedatum_woche))
  #xx_comp[!xx_comp%in%unique(sort(age_spec_incend$Meldedatum_woche))]
  num_weeks<-max(age_spec_incend$woche)                                  
  relev_weeks <- unique(sort(age_spec_incend$Meldedatum_woche))[unique(sort(age_spec_incend$Meldedatum_woche))>=(min(xx_comp)-7)]
  num_weeks<-length(relev_weeks)
  main_struct$age_struct$New_cases_week_distr <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  main_struct$age_struct$Death_cases_week_distr <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  main_struct$age_struct$Empir_Death_rate_week <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  main_struct$age_struct$Empir_Death_rate_week_del <- array(dim=c(num_weeks,length(main_struct$age_struct$ages)+2))
  
  colnames(main_struct$age_struct$New_cases_week_distr)<- c(main_struct$age_struct$ages,'Num Cases','ExpecDeathRate')
  colnames(main_struct$age_struct$Death_cases_week_distr)<- c(main_struct$age_struct$ages,'Num Cases','ExpecDeathRate')
  colnames(main_struct$age_struct$Empir_Death_rate_week)<- c(main_struct$age_struct$ages,'Mean Death','ExpecDeathRate')
  colnames(main_struct$age_struct$Empir_Death_rate_week_del)<- c(main_struct$age_struct$ages,'Mean Death','ExpecCritRate')
  death_delay<-2#1 2
  ind_day_good<-16
  for(ind_w in 1:num_weeks)
  {
    relev_ind<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
    num_all_casesi <- sum(age_spec_incend$NewConfCases[ relev_ind])
    num_all_deathi <- sum(age_spec_incend$NewDeaths[ relev_ind])
    main_struct$age_struct$New_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+1] <- num_all_casesi
    main_struct$age_struct$Death_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+1] <- num_all_deathi
    main_struct$age_struct$Empir_Death_rate_week[ind_w,1:length(main_struct$age_struct$ages)] <- age_spec_incend$AllDeaths[ relev_ind]/age_spec_incend$AllConfCases[ relev_ind]
    main_struct$age_struct$Empir_Death_rate_week[ind_w,1+length(main_struct$age_struct$ages)] <- sum(age_spec_incend$AllDeaths[ relev_ind])/sum(age_spec_incend$AllConfCases[ relev_ind])
    if(ind_w>death_delay)
    {
      relev_ind_prev<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w-death_delay])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
      main_struct$age_struct$Empir_Death_rate_week_del[ind_w,1:length(main_struct$age_struct$ages)] <- age_spec_incend$AllDeaths[ relev_ind]/age_spec_incend$AllConfCases[ relev_ind_prev]
      main_struct$age_struct$Empir_Death_rate_week_del[ind_w,1+length(main_struct$age_struct$ages)] <- sum(age_spec_incend$AllDeaths[ relev_ind])/sum(age_spec_incend$AllConfCases[ relev_ind_prev])
    }
  }
  age_mortality<-as.numeric(main_struct$age_struct$Empir_Death_rate_week_del[ind_day_good,1:length(main_struct$age_struct$ages)])
  for(ind_w in 1:num_weeks)
  {  
    relev_ind<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
    num_all_casesi <- sum(age_spec_incend$NewConfCases[ relev_ind])
    num_all_deathi <- sum(age_spec_incend$NewDeaths[ relev_ind])
    for(ind_a in 1:length(main_struct$age_struct$ages))
    {
      relev_indi <- relev_ind&(age_spec_incend$Altersgruppe==main_struct$age_struct$ages[ind_a])
      if(num_all_casesi>0)
      {
        main_struct$age_struct$New_cases_week_distr[ind_w,ind_a]<- age_spec_incend$NewConfCases[relev_indi]/num_all_casesi
      }else{
        main_struct$age_struct$New_cases_week_distr[ind_w,ind_a]<- 0
      }
      if(num_all_deathi>0)
      {
        main_struct$age_struct$Death_cases_week_distr[ind_w,ind_a]<- age_spec_incend$NewDeaths[relev_indi]/num_all_deathi
      }else{
        main_struct$age_struct$Death_cases_week_distr[ind_w,ind_a]<- 0
      }
    }
    #main_struct$age_struct$New_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$New_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(as.character(age_spec_prop$TotalDeathRate[age_spec_prop$Altersgruppe%in%main_struct$age_struct$ages])))
    #main_struct$age_struct$Death_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$Death_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(as.character(age_spec_prop$TotalDeathRate[age_spec_prop$Altersgruppe%in%main_struct$age_struct$ages])))
    main_struct$age_struct$New_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$New_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(main_struct$age_struct$Empir_Death_rate_week_del[ind_day_good,1:length(main_struct$age_struct$ages)]))
    main_struct$age_struct$Death_cases_week_distr[ind_w,length(main_struct$age_struct$ages)+2] <- sum(main_struct$age_struct$Death_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(main_struct$age_struct$Empir_Death_rate_week_del[ind_day_good,1:length(main_struct$age_struct$ages)]))
    main_struct$age_struct$Empir_Death_rate_week_del[ind_w,2+length(main_struct$age_struct$ages)] <- sum(main_struct$age_struct$New_cases_week_distr[ind_w,1:length(main_struct$age_struct$ages)]*as.numeric(as.character(age_spec_ICE$icu_guessed)))
  }
  
  apply(main_struct$age_struct$New_cases_week_distr[,1:6],1,sum)
  main_struct$age_struct$death_exp_w <-main_struct$age_struct$New_cases_week_distr[,length(main_struct$age_struct$ages)+2]
  main_struct$age_struct$death_dir_w <-main_struct$age_struct$Death_cases_week_distr[,length(main_struct$age_struct$ages)+2]
  main_struct$age_struct$relev_weeks <- relev_weeks
  
  death_exp_d <- vector(length = length(xx_comp))
  critical_exp_d <- vector(length = length(xx_comp))
  for(ind_d in 1:length(xx_comp))
  {
    #relev_ind <- max(max(which(xx_comp[ind_d]>relev_weeks)-death_delay),1)#length(relev_weeks)
    #death_exp_d[ind_d]<-main_struct$age_struct$death_exp_w[relev_ind]
    relev_ind <- max(max(which(xx_comp[ind_d]>relev_weeks)),1)#length(relev_weeks)
    death_exp_d[ind_d]<-main_struct$age_struct$Empir_Death_rate_week_del[,1+length(main_struct$age_struct$ages)][relev_ind] 
    critical_exp_d[ind_d] <-main_struct$age_struct$Empir_Death_rate_week_del[,2+length(main_struct$age_struct$ages)][relev_ind] 
    
  }
  help_v <- main_struct$age_struct$Empir_Death_rate_week_del[,1+length(main_struct$age_struct$ages)]
  ind_interp<-max(which(is.na(help_v)|help_v==0))+1
  death_exp_d[is.na(death_exp_d)|(death_exp_d==0)]<- main_struct$age_struct$Empir_Death_rate_week_del[ind_interp,1+length(main_struct$age_struct$ages)]#0#NANcase
  ##check:
  #par(mfrow=c(1,1)) 
  #plot(all_times,death_exp_d)
  #lines(relev_weeks,main_struct$age_struct$death_exp_w)
  main_struct$age_struct$death_exp_d <- death_exp_d
  main_struct$age_struct$norm_death_exp_d <- main_struct$age_struct$death_exp_d/mean(main_struct$age_struct$death_exp_d )
  
  main_struct$age_struct$critical_exp_d <- critical_exp_d
  main_struct$age_struct$norm_critical_exp_d <- main_struct$age_struct$critical_exp_d/mean(main_struct$age_struct$critical_exp_d )
  
  return(main_struct)
}

UploadAgeStructNext <- function(ind_indiv,root_path,main_struct)
{
  age_spec_incend<-read.csv(paste(root_path,'/data/','200709_age_specific_weekly_incidences_deaths.txt',sep =''), sep = "\t",dec = ',',header =T)
  age_spec_prop<-read.csv(paste(root_path,'/data/','200709_death_rate_week16.txt',sep =''), sep = "\t",dec = ',',header =T)
  age_spec_ICE<-read.csv(paste(root_path,'/data/','guess2_icu_age_distrib_germany_2020-07-14.txt',sep =''), sep = "\t",dec = ',',header =T)
  ####
  ####
  death_delay<-2#1 2
  ind_day_good<-16
  for(ind_w in 1:num_weeks)
  {
    relev_ind<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
    if(ind_w>death_delay)
    {
      relev_ind_prev<- ((age_spec_incend$Meldedatum_woche==relev_weeks[ind_w-death_delay])&(age_spec_incend$Altersgruppe%in%main_struct$age_struct$ages))
      main_struct$age_struct$Empir_Death_rate_week_del[ind_w,1:length(main_struct$age_struct$ages)] <- age_spec_incend$AllDeaths[ relev_ind]/age_spec_incend$AllConfCases[ relev_ind_prev]
      main_struct$age_struct$Empir_Death_rate_week_del[ind_w,1+length(main_struct$age_struct$ages)] <- sum(age_spec_incend$AllDeaths[ relev_ind])/sum(age_spec_incend$AllConfCases[ relev_ind_prev])
    }
  }
  
  
  apply(main_struct$age_struct$New_cases_week_distr[,1:6],1,sum)
  main_struct$age_struct$death_exp_w <-main_struct$age_struct$New_cases_week_distr[,length(main_struct$age_struct$ages)+2]
  main_struct$age_struct$death_dir_w <-main_struct$age_struct$Death_cases_week_distr[,length(main_struct$age_struct$ages)+2]
  main_struct$age_struct$relev_weeks <- relev_weeks
  
  death_exp_d <- vector(length = length(xx_comp))
  critical_exp_d <- vector(length = length(xx_comp))
  for(ind_d in 1:length(xx_comp))
  {
    #relev_ind <- max(max(which(xx_comp[ind_d]>relev_weeks)-death_delay),1)#length(relev_weeks)
    #death_exp_d[ind_d]<-main_struct$age_struct$death_exp_w[relev_ind]
    relev_ind <- max(max(which(xx_comp[ind_d]>relev_weeks)),1)#length(relev_weeks)
    death_exp_d[ind_d]<-main_struct$age_struct$Empir_Death_rate_week_del[,1+length(main_struct$age_struct$ages)][relev_ind] 
    critical_exp_d[ind_d] <-main_struct$age_struct$Empir_Death_rate_week_del[,2+length(main_struct$age_struct$ages)][relev_ind] 
    
  }
  help_v <- main_struct$age_struct$Empir_Death_rate_week_del[,1+length(main_struct$age_struct$ages)]
  ind_interp<-max(which(is.na(help_v)|help_v==0))+1
  death_exp_d[is.na(death_exp_d)|(death_exp_d==0)]<- main_struct$age_struct$Empir_Death_rate_week_del[ind_interp,1+length(main_struct$age_struct$ages)]#0#NANcase
  ##check:
  #par(mfrow=c(1,1)) 
  #plot(all_times,death_exp_d)
  #lines(relev_weeks,main_struct$age_struct$death_exp_w)
  main_struct$age_struct$death_exp_d <- death_exp_d
  main_struct$age_struct$norm_death_exp_d <- main_struct$age_struct$death_exp_d/mean(main_struct$age_struct$death_exp_d )
  
  main_struct$age_struct$critical_exp_d <- critical_exp_d
  main_struct$age_struct$norm_critical_exp_d <- main_struct$age_struct$critical_exp_d/mean(main_struct$age_struct$critical_exp_d )
  
  return(main_struct)
  
}
##################
##################


LoadECDCData<-function(root_path,ecdc_file_data_name,main_struct,ind_indiv)
{  
  ind_ <- ind_indiv#1
  measurements <- main_struct$measurements[[ind_]]
  measurements_daily <- main_struct$Daily$measurements[[ind_]]
  together_data<-read.csv(paste(root_path,'/data/',ecdc_file_data_name,sep =''), sep = "\t",dec = main_struct$inputdec,header =T)
  if(grep(pattern = "agestrat", x =ecdc_file_data_name)==1)#not relev_file_name !!!!21.05.2021
  {
    temp_measurements0 <- together_data[together_data$Altersgruppe=='all',]
    temp_measurements0_14 <- together_data[together_data$Altersgruppe=='A00-A14',]
    temp_measurements15_59 <- together_data[together_data$Altersgruppe=='A15-A59',]
    temp_measurements60_79 <- together_data[together_data$Altersgruppe=='A60-A79',]
    temp_measurements80 <- together_data[together_data$Altersgruppe=='A80+',]
    
  }else{
    temp_measurements0 <- together_data
  }
  
  if(main_struct$inputdec == ',')
  {  
    temp_measurements0$Meldedatum <- as.Date(temp_measurements0$DateRep,format = "%d.%m.%Y")#tryFormats = c("%Y-%m-%d","%d.%m.%Y"))#"%d.%m.%Y"
    temp_measurements0$DateRep<- as.Date(temp_measurements0$DateRep,format = "%d.%m.%Y")#, tryFormats = c("%Y-%m-%d","%d.%m.%Y"))#forma="%d.%m.%Y"
    temp_measurements0$covid_inICU_upscaled<-as.numeric(as.character(temp_measurements0$covid_inICU_upscaled))
  }
  else if(main_struct$inputdec == '.')
  {  
    temp_measurements0$Meldedatum <- as.Date(temp_measurements0$DateRep,format = "%Y-%m-%d")#tryFormats = c("%Y-%m-%d","%d.%m.%Y"))#"%d.%m.%Y"
    temp_measurements0$DateRep<- as.Date(temp_measurements0$DateRep,format = "%Y-%m-%d")#, tryFormats = c("%Y-%m-%d","%d.%m.%Y"))#forma="%d.%m.%Y"
    temp_measurements0$covid_inICU_upscaled<-as.numeric(as.character(temp_measurements0$covid_inICU_upscaled))
    
    temp_measurements0$NewConfCases7[1:3] <-  temp_measurements0$NewConfCases7[4]
    temp_measurements0$NewDeaths7[1:3] <- as.numeric(as.character(temp_measurements0$NewDeaths7[4]))
    temp_measurements0$NewConfCases7 <- as.numeric(as.character(temp_measurements0$NewConfCases7))
    bad_index <- is.na(temp_measurements0$NewConfCases7)
    temp_measurements0<-  temp_measurements0[!bad_index,]
  }  
  date_start <- min(temp_measurements0$Meldedatum,na.rm=T)#max(temp_measurements0$Meldedatum[temp_measurements0$AllConfCases==min(temp_measurements0$AllConfCases,na.rm=T)],na.rm=T)
  
  #######
  #####
  if(main_struct$weekly==T)
  {
    death_daily_col <- 'NewDeaths7'
    conf_daily_col <- 'NewConfCases7'
  }else{
    death_daily_col <- 'NewDeaths'
    conf_daily_col <- 'NewConfCases'
  }
  ind_allcases<-!is.na(temp_measurements0$AllConfCases)
  ind_deaths<-!is.na(temp_measurements0$AllDeaths)
  ind_cases_daily<-!is.na(temp_measurements0[[conf_daily_col]])
  ind_deaths_daily<-!is.na(temp_measurements0[[death_daily_col]])
  ind_time_relev <- (temp_measurements0$DateRep)>= as.Date(main_struct$covariates$DateStart[ind_])
  ind_critical <- !is.na(temp_measurements0$covid_inICU_upscaled)
  ind_critical <- ind_critical&(ind_time_relev)
  ind_allcases <- ind_allcases&(ind_time_relev)
  ind_deaths <- ind_deaths&(ind_time_relev)
  ind_cases_daily <- ind_cases_daily&(ind_time_relev)
  ind_deaths_daily <- ind_deaths_daily&(ind_time_relev)
  
  dates_critical<-temp_measurements0$DateRep[ind_critical]
  num_critical<-temp_measurements0$covid_inICU_upscaled[ind_critical]
  num_all1<-temp_measurements0$AllConfCases[ind_critical]
  
  numindiv_obs<-sum(ind_deaths)+sum(ind_allcases)+sum(ind_critical)
  numindiv_obs_daily<-sum(ind_deaths_daily)+sum(ind_cases_daily)+sum(ind_critical)
  all_cases_date<-temp_measurements0$Meldedatum[ind_allcases]
  all_cases<-temp_measurements0$AllConfCases[ind_allcases]
  cases_date_daily<-temp_measurements0$Meldedatum[ind_cases_daily]#-temp_cov[['DelayData']]
  #deaths_date_daily<-temp_measurements0$Meldedatum[ind_deaths_daily]-temp_cov[['DelayData']]
  cases_daily<-temp_measurements0[[conf_daily_col]][ind_cases_daily]
  #????:
 # if(min(temp_measurements0$DateRep)>main_struct$covariates$DateStart[ind_])
 # {
 #   numindiv_obs <- numindiv_obs+1
 #   all_cases<-c(0,all_cases)
 #   all_cases_date<-c(main_struct$covariates$DateStart[ind_],all_cases_date)
 # }
 # if(min(cases_date_daily)>main_struct$covariates$DateStart[ind_])
 # {
 #   numindiv_obs_daily <- numindiv_obs_daily+1
 #   cases_daily<-c(0,cases_daily)
 #   cases_date_daily<-c(main_struct$covariates$DateStart[ind_],cases_date_daily)
 # }
  measurements$namesmeasdata<-c("T","Y")
  measurements$date<-c(all_cases_date,#temp_measurements0$Meldedatum[ind_allcases],
                       temp_measurements0$Meldedatum[ind_deaths],
                       dates_critical)
  measurements$data <- rbind(as.numeric(measurements$date-min(measurements$date)),
                             c(all_cases,#temp_measurements0$AllConfCases[ind_allcases]
                               temp_measurements0$AllDeaths[ind_deaths],
                               temp_measurements0$covid_inICU_upscaled[ind_critical]))
  measurements$datatypes <- c(rep.int(x="Total",times=length(all_cases)),#sum(ind_allcases)
                              rep.int(x="Death",times=sum(ind_deaths)),
                              rep.int(x="Critical",times=sum(ind_critical)))
  measurements$scope<- as.numeric(max(measurements$date)-min(measurements$date))
  #measurements$RelevClinData <- (dates_critical>"2020-04-19")
  ind_ytypes1<-0
  measurements$length_data<- vector(length = length(unique(measurements$datatypes)))
  for(ind_ytype in unique(measurements$datatypes))
  {
    ind_ytypes1 <- ind_ytypes1 +1
    observedi<-measurements$data[measurements$namesmeasdata=="Y",measurements$datatypes==ind_ytype]
    measurements$length_data[ind_ytypes1]<- sum(!is.na(observedi))
  }   
  names(measurements$length_data) <- unique(measurements$datatypes)
  main_struct$measurements[[ind_]]<-measurements
  ###Daily:
  measurements_daily$namesmeasdata<-c("T","Y")
  measurements_daily$date<-c(cases_date_daily,#temp_measurements0$Meldedatum[ind_allcases],
                             temp_measurements0$Meldedatum[ind_deaths_daily],
                             dates_critical)
  measurements_daily$data <- rbind(as.numeric(measurements_daily$date-min(measurements_daily$date)),
                                   c(cases_daily,#temp_measurements0$AllConfCases[ind_allcases]
                                     temp_measurements0[[death_daily_col]][ind_deaths_daily],
                                     temp_measurements0$covid_inICU_upscaled[ind_critical]))
  measurements_daily$datatypes <- c(rep.int(x="Total",times=length(cases_daily)),#sum(ind_allcases)
                                    rep.int(x="Death",times=sum(ind_deaths_daily)),
                                    rep.int(x="Critical",times=sum(ind_critical)))
  measurements_daily$RelevClinData <- (dates_critical>"2020-04-19")
  main_struct$Daily$measurements[[ind_]]<-measurements_daily
  #main_struct$observed[[ind_]]<-list()
  #main_struct$initval[[ind_]]<-vector(length = max(main_struct$YTYPE[[ind_]]$listnum) )
  for(ind_ytype1 in 1:main_struct$YTYPES$num)
  {
    ind_ytype<-main_struct$YTYPES$listnum[ind_ytype1]
    ind_ytype_symb<-main_struct$YTYPES$list[ind_ytype1]
    #main_struct$YTYPE[[ind_]]$size[ind_ytype]<-sum(alloutcomesnum==main_struct$YTYPES$listnum[ind_ytype1]) #alloutcomesnum[index_id]
    #main_struct$YTYPE[[ind_]]$index[[ind_ytype]]<-(alloutcomesnum==main_struct$YTYPES$listnum[ind_ytype1])#alloutcomesnum[index_id]
    #index_id_y_ytype is unnecessary, use main_struct$YTYPE[[ind_]]$size[ind_ytype] instead!!!!!
   #obsolete and not correct!!!
    # if(sum(main_struct$YTYPE[[ind_]]$index[[ind_ytype]],na.rm = T)>0)
    #{
    #  main_struct$observed[[ind_]][[ind_ytype]]<-main_struct$measurements[[ind_]]$data[2,][main_struct$YTYPE[[ind_]]$index[[ind_ytype]]]#main_struct$Y[[ind_]] #alternative to matlab definition !!!! 
      
    #}
  }
  main_struct$completedates[[ind_]]<-CompleteDates(unique(main_struct$measurements[[ind_]]$date))
  main_struct$completeweekdates[[ind_]]<-weekdays(main_struct$completedates[[ind_indiv]])
  main_struct$erstsonntagind[ind_]<-min((1:length(main_struct$completeweekdates[[ind_]])) [main_struct$completeweekdates[[ind_]]=="Sonntag"])
  main_struct$letztsammstagind[ind_]<-max((1:length(main_struct$completeweekdates[[ind_]])) [main_struct$completeweekdates[[ind_]]=="Samstag"])
  main_struct$num_weeks <- (main_struct$letztsammstagind[[ind_]]-main_struct$erstsonntagind[[ind_]]+1)/7
  if('DateStepsPDeath'%in%names(main_struct$covariates))
  {
    main_struct$covariates$DateStepsPDeath[[ind_]]<-c(main_struct$covariates$DateStepsPDeath[[ind_]],max(main_struct$completedates[[ind_]]))
    
    main_struct$covariates$DateStepsPCrit[[ind_]]<-c(main_struct$covariates$DateStepsPCrit[[ind_]],max(main_struct$completedates[[ind_]]))
  }
  #####
  temp_name <- gsub(pattern="Truncated", replacement='', x=ecdc_file_data_name, ignore.case = FALSE, perl = FALSE,
                    fixed = FALSE, useBytes = FALSE)
  temp_name <- gsub(pattern=".txt", replacement='', x=temp_name, ignore.case = FALSE, perl = FALSE,
                    fixed = FALSE, useBytes = FALSE)
  main_struct$resultspath<-paste(root_path,'/PredictionsValidations','/',temp_name,'/',sep='')
  if('dates_treatments'%in%names(main_struct))
  {
    dim_treatments <- dim(main_struct$dates_treatment)
    Row <- (1:dim_treatments[1])[as.numeric(as.character(main_struct$dates_treatment$ID))==ind_indiv]
    #relev_indiv_names<-main_struct$names[main_struct$subs_indivi[ind_indiv,]==1]#all subjects have the same individual parameters!
    for (indpari in 1:(dim_treatments[2]-1))#(Col in 1:(dim(temp_indiv)[2]-1)) !!!!
    {
      #print(c(Row,Col)) 
      #Col <- (1:main_struct$numindivpar1)[names(main_struct$dates_treatment)[[indpari+1]]==relev_indiv_names]
      main_struct$indiv$psi[Row,names(main_struct$dates_treatment)[[indpari+1]]] <- as.numeric(as.Date(main_struct$dates_treatment[Row,indpari+1],tryFormats = "%d.%m.%Y")-min(main_struct$measurements[[ind_indiv]]$date))#as.Date(main_struct$covariates$DateStart[ind_indiv]))
      main_struct$init_mean[names(main_struct$dates_treatment)[[indpari+1]]==main_struct$names] <- main_struct$indiv$psi[Row,names(main_struct$dates_treatment)[[indpari+1]]]#correct dates for population parameters!
    }
    
    main_struct$indiv$phi[Row,] <- RevTransformParametersAll(main_struct$indiv$psi[Row,],main_struct)
    main_struct$indiv$phiopt[[Row]] <- main_struct$indiv$phi[Row,main_struct$subs_indivi[Row,]>0]
    main_struct$treat_par_ind<-setdiff(grep(pattern ='treat', x = main_struct$names),grep(pattern ='date_treat', x = main_struct$names))
    main_struct$treat_par_names <- main_struct$names[main_struct$treat_par_ind]
    main_struct$treat_par_num <- length(main_struct$treat_par_names)
    
    main_struct$treat_date_ind<-grep(pattern ='date_treat', x = main_struct$names)
    main_struct$treat_date_names <- main_struct$names[main_struct$treat_date_ind]
    #main_struct$treat_par_num <- length(main_struct$treat_par_names)
  }
  return(main_struct) 
}  

UploadTestStrategy <- function(ind_indiv,root_path,main_struct)
{
  main_struct$TestStrategy<- list()
  all_times<-unique(main_struct$measurements[[ind_indiv]]$date)
  xx_comp<-CompleteDates(all_times)
  test_struct<-read.csv(paste(root_path,'/data/','tests_covid_until_2020-09-13_1.txt',sep =''), sep = "\t",dec = ',',header =T)
  test_struct$week.start<-as.Date(as.character(test_struct$week.start), tryFormats = c("%d.%m.%Y"))
  test_struct$week.end<-as.Date(as.character(test_struct$week.end), tryFormats = c("%d.%m.%Y"))
  num_days <- length(xx_comp)
  num_weaks<-dim(test_struct)[1]
  main_struct$TestStrategy$Daily_tests<-vector(length=num_days)
  main_struct$TestStrategy$PrPos_tests<-vector(length=num_days)
  main_struct$TestStrategy$Pos_tests<-vector(length=num_days)
  for(ind_d in 1:length(xx_comp))
  {
    relev_ind <- (xx_comp[ind_d]>test_struct$week.start)&(xx_comp[ind_d]<=test_struct$week.end)
    main_struct$TestStrategy$Daily_tests[ind_d] <- test_struct$Anzahl_Testungen[relev_ind]/7
    main_struct$TestStrategy$PrPos_tests[ind_d] <- test_struct$Positivenquote_proz[relev_ind]
    main_struct$TestStrategy$Pos_tests[ind_d] <- test_struct$Positiv_getestet[relev_ind] 
  }
  return(main_struct)
  
}
CreateNameOutput<-function(ind_indiv,main_struct)
{
  
 
 # if(main_struct$model_opt==3)
 # {
    #colnames(sim_res)<-c('Sc','IAc','ISc','ISc2', 'Cc','Dc', 'Rc1','Rc2','Rc0') 
    compartments<-main_struct$compartments[[1]]
    colnames_sim_res <-vector(length=compartments$NumofVariables+1)
    colnames_sim_res[compartments$Sc_ind]<- 'Sc'
    colnames_sim_res[compartments$Ec_ind]<- 'Ec'
    colnames_sim_res[compartments$IAc_ind] <- c('IAc1','IAc2','IAc3')
    colnames_sim_res[compartments$ISc_ind] <- c('ISc1','ISc2','ISc3')
    colnames_sim_res[compartments$Cc_ind] <- paste('Cc',1:compartments$Ccnum,sep='')
    colnames_sim_res[compartments$Dc_ind] <- 'Dc'
    colnames_sim_res[compartments$Rc_ind] <- 'Rc'
    colnames_sim_res[compartments$StCare_ind] <- 'StCare'
    colnames_sim_res[compartments$NumofVariables+1]<-'SumIsc'
    if(main_struct$label_treat_info==1)
    {
      colnames_sim_res<-c(colnames_sim_res,c('vr1','vr2','pcrit','pdeath'))
    }  
    #colnames(sim_res) <- colnames_sim_res
 # }
    
    
    
    return(colnames_sim_res)
}
FindActualData<-function(ind_indiv,root_path,label_clean,label_age)
{
  all_file_names <- list.files(paste(root_path,'/data/',sep =''))
  if(ind_indiv==1)
  {
    help_pattern <- 'germany'
  }else if(ind_indiv==9)
  {
    help_pattern <- 'saxony'
  }
  relev_f_ind <- grep(pattern = paste('042_datint_ecdc_',help_pattern,'_', sep =''), x=all_file_names, ignore.case = FALSE, perl = FALSE, value = FALSE,
                      fixed = FALSE, useBytes = FALSE, invert = FALSE)
  bad_f_ind <- grep(pattern = 'Truncated', x=all_file_names, ignore.case = FALSE, perl = FALSE, value = FALSE,
                    fixed = FALSE, useBytes = FALSE, invert = FALSE)
 
  relev_f_ind <- setdiff(relev_f_ind, bad_f_ind)
  if(label_clean==1)
  {  
    clean_ind <- grep(pattern = 'cleaned', x=all_file_names, ignore.case = FALSE, perl = FALSE, value = FALSE,
                      fixed = FALSE, useBytes = FALSE, invert = FALSE)
    relev_f_ind <- intersect(relev_f_ind,clean_ind)
  }
  if(label_age==1)
  {
    age_ind <- grep(pattern = 'agestrat', x=all_file_names, ignore.case = FALSE, perl = FALSE, value = FALSE,
                      fixed = FALSE, useBytes = FALSE, invert = FALSE)
    relev_f_ind <- intersect(relev_f_ind,age_ind)
  }
  sort(all_file_names[relev_f_ind])
  relev_datesch <- sub(pattern = paste('042_datint_ecdc_',help_pattern,'_', sep=''), replacement ='',x=all_file_names[relev_f_ind])
  relev_datesch <- sub(pattern = '_v3.txt', replacement ='',x=relev_datesch)
  #which.max(as.Date(relev_datesch, format = "%Y-%m-%d" ))
  return(all_file_names[relev_f_ind][which.max(as.Date(relev_datesch, format = "%Y-%m-%d" ))])
}
  