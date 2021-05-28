rm(list = ls())
library(Rcpp)
#root_path <- "F:/2003_covid_modell"
#root_path <- "D:/2003_covid_modell"
#root_path <- "R:/modellclub/2003_covid_modell"
root_path <- "R:/modellclub/2003_covid_modell/github/LeipzigIMISE-SECIR"
#path_code<-paste(root_path,'/StochODEMS_YK/',sep='')#working path, yuri
#path_code<-paste(root_path,'/github/LeipzigIMISE-SECIR/scripte/',sep='')#working path Sandra
path_code<-paste(root_path,'/scripte/',sep='')#working path Sandra
setwd(path_code)
source("LoadCOVIDFunct1.R")
label_data <- 'LandMuKIAdv'#'LandMuPolyKIAdv'
label_clean <- 1#1#0
label_age <- 1#0
ind_indiv <- 1#9  #1-Deutschland, 9 - Sachsen
#main_struct0 <- SimpleLoad(root_path,data_name = '040_datint_rki_meldedat_kreise_bl_d_2020-10-05.txt',label_data,label_clean)
helpstr <- AdvancedLoad(path_code,root_path,label_data, label_clean,  ind_indiv = ind_indiv, inputdec = '.' )
###main structure, governing simulations
main_struct0 <- helpstr$main_struct
main_struct0$resultspath<- path_code<-paste(root_path,'/results/',sep='')
### "stochastic object2- structure, dwetermines time steps. Our model can also be stochastic (not now, in perspective)
stoch_obj0 <- helpstr$stoch_obj0

###all model's parameters:
param_ <- helpstr$param_
##all individual transformed parameters, initial conditions and stochastic perturbations (currently zero)
xtot <- helpstr$xtot
####individual index:
#ind_indiv <- 1#1#9#10
##1- Germany, 9-Saxony
main_struct0$treatpardir <- 2
sourceCpp('solvecovidmodelMuKI.cpp')
##updating of model's parameters
param_<-ConstructParam(x=xtot$phiopt,ind_indiv,main_struct0) 
####Simulation results:
sims<-SimulCppCovidModel(param_,xtot,ind_indiv,main_struct0)
sims<-NameOutput(sims,main_struct0)#naming all columns


###Plots of results
#label_cummul - label of cumulative or daily simulations/measurements: 1 or 0, respectively
#date_option=1 time is given in the date format. currently is the single option
#log_option: normal or logarithmic transformation of Y axis: 0 or 1, respectively
#label_save if 1- save results
PlotCovidFits(sims,main_struct0,ind_indiv,param_,label_save=0,log_option=1,date_option=1,label_cummul=0)
PlotCovidFits(sims,main_struct0,ind_indiv,param_,label_save=0,log_option=1,date_option=1,label_cummul=0)
PlotCovidFits(sims,main_struct0,ind_indiv,param_,label_save=0,log_option=0,date_option=1,label_cummul=0)
PlotCovidFits(sims,main_struct0,ind_indiv,param_,label_save=0,log_option=1,date_option=1,label_cummul=1)
PlotCovidFits(sims,main_struct0,ind_indiv,param_,label_save=0,log_option=0,date_option=1,label_cummul=1)

goalf<-GoalFuncComp(xtot=xtot,main_struct=main_struct0,opt_mod=opt_mod,label_opt=label_opt,
                                                  nLL_indiv=rep.int(x=0,times=main_struct0$numindiv),ind_indiv1=ind_indiv)#
#liminfo shows how close are the parameters' values to their respective values
liminfo <- CompareLimits(x=xtot$phiopt$indiv[[ind_indiv]],ind_indiv,main_struct0)
#Sensitivanalysis:
#Provide sensitivity analysis for given perturbation delt
#rel_opt whether relative perturbation is proveded (1- yes, 0- no)
#transform_opt- whether transformed parameter is perturbed (1 yes, 0 no)
#for example transform_opt=0,rel_opt=1 means, that absolute parameter's value is multiplied by delt. If we want to check +-2.5% we should use delt= 1.025 or 0.975 in this case
#transform_opt=0,rel_opt=1 must not go beyond the respective limits, which can be determined by liminfo

sens_results1_1 <- Sensitivanalysis(xtot,delt=1.1,ind_indiv,main_struct0,transform_opt=0,rel_opt=1)
sens_results1_025 <- Sensitivanalysis(xtot,delt=1.025,ind_indiv,main_struct0,transform_opt=0,rel_opt=1)
sens_results0_975 <- Sensitivanalysis(xtot,delt=0.975,ind_indiv,main_struct0,transform_opt=0,rel_opt=1)
sens_results1_01 <- Sensitivanalysis(xtot,delt=1.01,ind_indiv,main_struct0,transform_opt=0,rel_opt=1)
sens_results1_001 <- Sensitivanalysis(xtot,delt=1.001,ind_indiv,main_struct0,transform_opt=0,rel_opt=1)
sens_results1_0001 <- Sensitivanalysis(xtot,delt=1.0001,ind_indiv,main_struct0,transform_opt=0,rel_opt=1)

#Sensitivanalysis is based on the function PerturbPar, which perturbs given parameter ( parameter's name par_name) by perturbation delt
#rel_opt whether relative perturbation is proveded (1- yes, 0- no)
#transform_opt- whether transformed parameter is perturbed (1 yes, 0 no)

#liminfo, PerturbPar and Sensitivanalysis are in COVIDfunctions.R 