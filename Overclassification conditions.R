##############################################################
#
#       RESPONSE-TIME THRESHOLD PROCEDURE COMPARISONS        #
#
##############################################################

#clearing memory
rm(list=ls())

######### Specifying packages needed for analyses ############

library(dplyr) #combines the simulation results into a condition matrix
library(SimDesign) #calculate RMSE
library(mirt) #estimating IRT item parameters
#install.packages("robustbase",lib='C:\\R-4.0.0\\library\\')
#library(robustbase,lib='C:\\R-4.0.0\\library\\')
#library(RcppZiggurat,lib='C:\\R-4.0.0\\library\\')
#library(Rfast,lib='C:\\R-4.0.0\\library\\') #package needed to calculate column medians
#library(sirt,lib='C:\\R-4.0.0\\library\\') #package for estimating HG and MW models

##################### Variables that won't be changed ##############
N<-1000 #total sample size

reps<-100 #number of reps for study
reps.EM.imputed<-100 #number of reps for EM-imputation procedure

#################### creating conditions matrix ####################


# level codings will be specified in if else statements


#MISCLASSIFICATION PERCENTAGE
misclassify.percent <-c(.02,.04,.06,.08,.12) #corresponds to underclassifying 1, 2, 3, 4, and 6 items on average per misclassified examinee

#TEST DIFFICULTY
#there are two levels: easy [mean b = -1] and moderate [mean b = 0]
test.difficulty <- c("Easy", "Moderate")

# We need to get a conditions matrix for a fully crossed design. In order to do this

con <- expand.grid(misclassify.percent)

##################### Setting working directory ####################

setwd("G:\\My Drive\\UMN\\Research 2021-2022\\Comparison of RT Threshold Procedures to RG Misclassifications\\Simulation\\")

#reading in item parameters
item.3PL.pars<-as.matrix(read.csv("Generating 3pl item parameters (50 items, moderate).csv"))
I<-nrow(item.3PL.pars) #number of items
a.3PL.true<-as.matrix(item.3PL.pars[,1],ncol=1)
b.3PL.true<-as.matrix(item.3PL.pars[,2],ncol=1)
c.3PL.true<-as.matrix(item.3PL.pars[,3],ncol=1)


#creating matrix for overall results by condition

overall.results<-matrix(-999,nrow=nrow(con),ncol=51) #creating a matrix to place results in by condition, which will be used for the overall descriptive results

colnames(overall.results)<-c(
  
  "Bias.a.ML","Bias.a.penalized","Bias.a.EM","Bias.a.EM.imputed","Bias.a.HG",
  "Bias.b.ML","Bias.b.penalized","Bias.b.EM","Bias.b.EM.imputed","Bias.b.HG",
  # "Bias.c.ML","Bias.c.penalized","Bias.c.EM","Bias.c.EM.imputed","Bias.c.HG",
  "Bias.theta.ML","Bias.theta.penalized","Bias.theta.EM","Bias.theta.EM.imputed","Bias.theta.HG",
  
  #"SE.a.ML","SE.a.penalized","SE.a.EM","SE.a.EM.imputed","SE.a.HG",
  #"SE.b.ML","SE.b.penalized","SE.b.EM","SE.b.EM.imputed","SE.b.HG",
  #"SE.c.ML","SE.c.penalized","SE.c.EM","SE.c.EM.imputed","SE.c.HG",
  #"SE.theta.ML","SE.theta.penalized","SE.theta.EM","SE.theta.EM.imputed","SE.theta.HG",
  
  "RMSE.a.ML","RMSE.a.penalized","RMSE.a.EM","RMSE.a.EM.imputed","RMSE.a.HG",
  "RMSE.b.ML","RMSE.b.penalized","RMSE.b.EM","RMSE.b.EM.imputed","RMSE.b.HG",
  #"RMSE.c.ML","RMSE.c.penalized","RMSE.c.EM","RMSE.c.EM.imputed","RMSE.c.HG",
  "RMSE.theta.ML","RMSE.theta.penalized","RMSE.theta.EM","RMSE.theta.EM.imputed","RMSE.theta.HG",
  "BIC.ML","BIC.penalized","BIC.EM","BIC.EM.imputed","BIC.HG","r.theta.RG.HG",
  "Non.Converge.ML","Non.Converge.penalized","Non.Converge.EM","Non.Converge.EM.imputed","Non.Converge.HG",
  
  "Bias.unmot.theta.ML","Bias.unmot.theta.penalized","Bias.unmot.theta.EM","Bias.unmot.theta.EM.imputed","Bias.unmot.theta.HG",
  "RMSE.unmot.theta.ML","RMSE.unmot.theta.penalized","RMSE.unmot.theta.EM","RMSE.unmot.theta.EM.imputed","RMSE.unmot.theta.HG"
  
)



#condition loop starts here


for(count.con in 1:nrow(con)) { 
  
  
  results.condition<-matrix(-999,nrow=reps,ncol=51) #creating a matrix to place results in by condition, which will be used for the overall descriptive results
  
  
  
  ######### REPLICATION LOOP STARTS HERE
  for (r in 1:reps){
    #setting seed
    set.seed(123+(r*count.con)+r) #rep to this counter
    ##################### Sampling ability parameters ##############
    
    theta<-matrix(rnorm(N,0,1),ncol=1) #motivated ability, standard normal
    
    
    ##################### Generating response probabilities ##############
    
    #creating function to compute 3pl logistic model
    PL3<-function(thetas,a.MC,b.MC,c.MC,N,I.MC){
      alphas.betas <-matrix(a.MC * b.MC,ncol=1)
      alphas.thetas <-matrix(thetas,ncol=1) %*% matrix(a.MC,nrow=1)
      logit <-alphas.thetas - t(matrix(alphas.betas,I.MC,N))
      c.par <- t(matrix(c.MC,I.MC,N))
      probs <- (c.par+((1-c.par)*(1/(1+exp(-logit)))))
      return(probs)}
    
    probs<-PL3(theta,a.3PL.true,b.3PL.true,c.3PL.true,N,I) #generating response probabilities for entire sample
    theta.probs<-cbind(theta,probs)
    probs.high <- theta.probs[ which((theta.probs[,1])>quantile(theta.probs[,1],.69)), ][,-1] #SAMPLING EXAMINEES WITH 0.5 SD above average (top 31% of examinees)
    theta.low<-matrix(theta.probs[ which((theta.probs[,1])<=quantile(theta.probs[,1],.69)), ][,1],ncol=1)
    theta.high <- matrix(theta.probs[ which((theta.probs[,1])>quantile(theta.probs[,1],.69)), ][,1],ncol=1) #SAMPLING EXAMINEES WITH 0.5 SD above average (top 31% of examinees)
    probs.low<-theta.probs[ which((theta.probs[,1])<=quantile(theta.probs[,1],.69)), ][,-1]
    theta.combined<-rbind(theta.high,theta.low)
    N.over<-nrow(probs.high) #sample size of examinees overclassified
    N.low<-N-N.over #sample size of examinees overclassified
    
    ######################### SAMPLING RANDOM NUMBER OF RGs FOR EACH UNMOTIVATED SIMULEE #######################
    #this is the random number of RGs for condition with 10% RGs
    
    aa<-N.over*I #total number of item responses (N x I)
    bb<-aa*con[count.con,1] #number of RG responses to sample
    prob.sample<-1/N.over #probability of sampling each simulee
    samples<-1 #number of samples of random numbers to pull 
    
    RG.1<-matrix(rmultinom(samples, size = bb,prob=c(rep(prob.sample,N.over))),ncol=samples,nrow=N.over)
    #sampling number of RG responses per simulee

    
    ########################  IMPUTING RG MISCLASSIFICATIONS ##############################
    
    
    prob.unmot.NA<-matrix(-999,nrow=N.over,ncol=I)
    prob.unmot.zero<-matrix(-999,nrow=N.over,ncol=I)
    
    RG.matrix.mot<-matrix(0,nrow=N.low,ncol=I) #CREATING RG Matrix for Motivated Sample (all motivated)
    
    RG.matrix.unmot<-matrix(0,nrow=N.over,ncol=I) #CREATING RG Matrix for Sample with Misclassifications
    
    for(p in 1:N.over){
      
      
      ############## SELECTING ITEMS TO RECEIVE RG ##################

        
        #rank ordering the true probability of success from lowest to highest; ties are randomly ordered 
        rank.prob.success.guesser<-order( probs.high[p,],decreasing=TRUE) #order in descending order the true probability of success from lowest to highest; ties are randomly ordered 
        item.id<-rank.prob.success.guesser[c(1:RG.1[p,])]
        
     
      ######################### UPDATING GUESSING PROBABILITIES & RG MATRIX ############################

      #replacing 1 with 0 to represent RG responses
      RG.matrix.unmot[p,item.id]<-1
      
      
      #print(p)
    }
    ############## NEED TO RECODE PROBABILITIES BASED ON UPDATED RG MATRIX
    for(p in 1:N.over){
      for (i in 1:I){
        ifelse(RG.matrix.unmot[p,i]==1,prob.unmot.NA[p,i]<-NA,prob.unmot.NA[p,i]<-probs.high[p,i])
        ifelse(RG.matrix.unmot[p,i]==1,prob.unmot.zero[p,i]<-0,prob.unmot.zero[p,i]<-probs.high[p,i])
        
      }
    }
    
    
    RG.combined<-rbind(RG.matrix.unmot,RG.matrix.mot) #combining RG classifications
    
    RG.combined[N,]<-1 #adding in case here to make sure that I can estimate model (need at least two response options for each item)
    
    #combining response probabilities for unmotivated and motivated samples
    probs.combined.ML<-rbind(probs.high,probs.low)
    probs.combined.NA<-rbind(prob.unmot.NA,probs.low)
    probs.combined.zero<-rbind(prob.unmot.zero,probs.low)
    
    ########################################################################################################
    # ASSIGNING 1'S AND 0'S USING RANDOM NUMBERS #
    ########################################################################################################
    #getting random probabilities
    
    
    rnd.unif <- matrix(runif(N*I, 0, 1), nrow = N, ncol = I)
    
    #transforming probs to 0/1 for when noneffortful responses are in the data matrix
    IR.combined.ML <- (ifelse(probs.combined.ML > rnd.unif,1,0)) #Coding responses as 0/1 for IRT software#
    IR.combined.NA <- (ifelse(probs.combined.NA > rnd.unif,1,0)) #Coding responses as 0/1 for IRT software#
    IR.combined.zero <- (ifelse(probs.combined.zero > rnd.unif,1,0)) #Coding responses as 0/1 for IRT software#
    
    #adding column names for use in MIRT Package
    colnames(IR.combined.ML)<-c(1:I)
    colnames(IR.combined.NA)<-c(1:I)
    colnames(IR.combined.zero)<-c(1:I)
    
    
    #############################################################################
    #                                                                           #
    #                   ESTIMATING ITEM PARAMETERS & ABILITIES                  #
    #                                                                           #
    #############################################################################
    
    ###### ML scoring (ignoring RG responses) ########
    
    ML<-mirt(IR.combined.ML, 1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000))
    ipars.ML<-coef(ML,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
    
    
    results.condition[r,37]<-ifelse(extract.mirt(ML,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(results.condition[r,37]==0){ 
      a.ML<-ipars.ML[c(seq(1,200,4))] 
      b.ML<-ipars.ML[c(seq(2,200,4))] 
     #c.ML<-ipars.ML[c(seq(3,200,4))] 
      theta.ML<-as.matrix(fscores(ML,method='ML',max_theta=4))
    } else {
      a.ML<-rep(NA,I) 
      b.ML<-rep(NA,I) 
      #c.ML<-rep(NA,I) 
      theta.ML<-rep(NA,N) 
    }

    
    ###### penalized scoring ######
    
    penalized<-mirt(IR.combined.zero,  1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000))
    ipars.penalized<-coef(penalized,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
    
    results.condition[r,38]<-ifelse(extract.mirt(penalized,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(results.condition[r,38]==0){ 
      a.penalized<-ipars.penalized[c(seq(1,200,4))] 
      b.penalized<-ipars.penalized[c(seq(2,200,4))] 
      #c.ML<-ipars.ML[c(seq(3,200,4))] 
      theta.penalized<-as.matrix(fscores(penalized,method='ML',max_theta=4))
    } else {
      a.penalized<-rep(NA,I) 
      b.penalized<-rep(NA,I) 
      #c.ML<-rep(NA,I) 
      theta.penalized<-rep(NA,N) 
    }
    
    
    
    ###### EM scoring ######
    EM<-mirt(IR.combined.NA, 1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000))
    ipars.EM<-coef(EM,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
    
    results.condition[r,39]<-ifelse(extract.mirt(EM,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(results.condition[r,39]==0){ 
      a.EM<-ipars.EM[c(seq(1,200,4))] 
      b.EM<-ipars.EM[c(seq(2,200,4))] 
      c.EM<-ipars.ML[c(seq(3,200,4))] 
      theta.EM<-as.matrix(fscores(EM,method='ML',max_theta=4))
    } else {
      a.EM<-rep(NA,I) 
      b.EM<-rep(NA,I) 
      c.EM<-rep(NA,I) 
      theta.EM<-rep(NA,N) 
    }
    
    
    ##### EM-IMPUTATION scoring #####
    
    #imputing probabilities
    imputed.probabilities<-PL3(theta.EM,a.EM,b.EM,c.EM,N,I) #based on the theta and item parameter estimates from EM scoring
    
    
    a.EM.imputed<-matrix(-999,nrow=I,ncol=reps.EM.imputed)
    b.EM.imputed<-matrix(-999,nrow=I,ncol=reps.EM.imputed)
    #c.EM.imputed<-matrix(-999,nrow=I,ncol=reps.EM.imputed)
    converged.EM.imputed<-matrix(-999,nrow=1,ncol=reps.EM.imputed)
    theta.EM.imputed<-matrix(-999,nrow=N,ncol=reps.EM.imputed)
    
    
    if(results.condition[r,39]==1){
      
      a.EM.imputed<-matrix(rep(matrix(rep(NA,I),ncol=1),reps.EM.imputed),ncol=reps.EM.imputed) 
      b.EM.imputed<-matrix(rep(matrix(rep(NA,I),ncol=1),reps.EM.imputed),ncol=reps.EM.imputed) 
      #c.ML<-rep(NA,I) 
      theta.EM.imputed<-matrix(rep(matrix(rep(NA,N),ncol=1),reps.EM.imputed),ncol=reps.EM.imputed)
      converged.EM.imputed<-matrix(1,nrow=1,ncol=reps.EM.imputed)
        
    } else if (results.condition[r,39]==0){
    
    
    for(r.EM.impute in 1:reps.EM.imputed){
      
      rnd.unif.impute <- matrix(runif(N*I, 0, 1), nrow = N, ncol = I)
      
      #transforming probs to 0/1 for when noneffortful responses are in the data matrix
      IR.EM.imputed <- (ifelse(imputed.probabilities > rnd.unif.impute,1,0)) #Coding responses as 0/1 for IRT software#
      colnames(IR.EM.imputed)<-rep(1:I)
      
      EM.imputed<-mirt(IR.EM.imputed,  1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000))
      ipars.EM.imputed<-coef(EM.imputed,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
      
      converged.EM.imputed[,r.EM.impute]<-ifelse(extract.mirt(EM.imputed,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
      
      #if model converged, take parameter estimates. if not, impute missing values
      if(extract.mirt(EM.imputed,'converged')==TRUE){ 
        a.EM.imputed[,r.EM.impute]<-ipars.EM.imputed[c(seq(1,200,4))] 
        b.EM.imputed[,r.EM.impute]<-ipars.EM.imputed[c(seq(2,200,4))] 
        #c.ML<-ipars.ML[c(seq(3,200,4))] 
        theta.EM.imputed[,r.EM.impute]<-as.matrix(fscores(EM.imputed,method='ML',max_theta=4))
      } else if (extract.mirt(EM.imputed,'converged')==FALSE){
        a.EM.imputed[,r.EM.impute]<-matrix(rep(NA,I),ncol=1) 
        b.EM.imputed[,r.EM.impute]<-rep(NA,I) 
        #c.ML<-rep(NA,I) 
        theta.EM.imputed[,r.EM.impute]<-rep(NA,N) 
      }
      
    }
    
    } #closes if/else statement above
    
    
    a.EM.imputed<-matrix(rowMeans(a.EM.imputed,na.rm=TRUE),ncol=1)
    b.EM.imputed<-matrix(rowMeans(b.EM.imputed,na.rm=TRUE),ncol=1)
    #c.EM.imputed<-matrix(rowMeans(c.EM.imputed),ncol=1)
    theta.EM.imputed<-matrix(rowMeans(theta.EM.imputed,na.rm=TRUE),ncol=1)
    
    results.condition[r,40]<-matrix(rowMeans(converged.EM.imputed,na.rm=TRUE),ncol=1) #average reps nonconverged
    
    ##### HG scoring ########      
    
    #COMBINING IR AND RG MATRICES TO ESTIMATE 2-FACTOR CORRELATED-TRAITS MODEL#
    
    
    #combining data for both factors into single dataset
    data.factor<-data.matrix(cbind(IR.combined.NA,RG.combined))
    
    colnames(data.factor)<-c(1:(2*I)) #adding in column names to estimate data in mirt r package
    
    #estimating multidimensional 2pl modelwith mirt.model definition
    model <- 'F1 = 1-50
    F2 = 51-100
    COV = F1*F2'
    
    item.type<-c(rep('2PL',50),rep('Rasch',50))
    HG<-mirt(data.factor, model,itemtype =item.type,guess=as.vector(c(rep(.25,50),rep(0,50))), TOL = .0001, technical = list(NCYCLES = 10000))
    ipars.HG<-coef(HG,as.data.frame=TRUE)
    
    results.condition[r,41]<-ifelse(extract.mirt(HG,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(results.condition[r,41]==0){ 
      #need to transform parameter estimates onto theta scale
      a.HG<-ipars.HG[c(seq(1,250,5))] 
      b.HG<--(ipars.HG[c(seq(3,250,5))])/a.HG # b=-threshold/slope
      #c.HG<-ipars.HG[c(seq(4,250,5))] 
      results.condition[r,36]<-ipars.HG[504] #Covariance between theta and RG propensity
      theta.HG<-matrix(fscores(HG,method='ML',max_theta=4)[,1],ncol=1) #setting upper and lower limits of theta at 4 and -4
      
    } else if (results.condition[r,41]==0){
      a.HG<-rep(NA,2*I) 
      b.HG<-rep(NA,2*I)  #I am not sure if there is a problem; trying to convert cfa parameters to irt b=threshold/slope
      #c.HG<-ipars.HG[c(seq(4,250,5))] 
      results.condition[r,36]<-NA #Covariance between theta and RG propensity
      theta.HG<-rep(NA,N) #setting upper and lower limits of theta at 4 and -4
      
    }
    
   
    ###################################
    #                                 #
    # CALCULATE DEPENDENT VARIABLES   #
    #                                 #
    ###################################
    #have to drop missing data to run simdesign functions
    
    #dropping missing cases for a.ML
    a.ML.data<-cbind(a.ML,a.3PL.true)
    a.ML.recoded<-matrix(a.ML.data[complete.cases(a.ML.data), ],ncol=2)
    colnames(a.ML.recoded)<-c("a.ML","a.True")
    
    #dropping missing cases for b.ML
    b.ML.data<-cbind(b.ML,b.3PL.true)
    b.ML.recoded<-matrix(b.ML.data[complete.cases(b.ML.data), ],ncol=2)
    colnames(b.ML.recoded)<-c("b.ML","b.True")
    
    #dropping missing cases for theta.ML
    theta.ML.data<-cbind(theta.ML,theta.combined)
    theta.ML.recoded<-matrix(theta.ML.data[complete.cases(theta.ML.data), ],ncol=2)
    theta.ML.recoded <- theta.ML.recoded[!is.infinite(rowSums(theta.ML.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.ML.recoded)<-c("theta.ML","theta.True")
    
 
    #dropping missing cases for a.penalized
    a.penalized.data<-cbind(a.penalized,a.3PL.true)
    a.penalized.recoded<-matrix(a.penalized.data[complete.cases(a.penalized.data), ],ncol=2)
    colnames(a.penalized.recoded)<-c("a.penalized","a.True")
    
    #dropping missing cases for b.penalized
    b.penalized.data<-cbind(b.penalized,b.3PL.true)
    b.penalized.recoded<-matrix(b.penalized.data[complete.cases(b.penalized.data), ],ncol=2)
    colnames(b.penalized.recoded)<-c("b.penalized","b.True")
    
    #dropping missing cases for theta.penalized
    theta.penalized.data<-cbind(theta.penalized,theta.combined)
    theta.penalized.recoded<-matrix(theta.penalized.data[complete.cases(theta.penalized.data), ],ncol=2)
    theta.penalized.recoded <- theta.penalized.recoded[!is.infinite(rowSums(theta.penalized.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.penalized.recoded)<-c("theta.penalized","theta.True")
    
    
    #dropping missing cases for a.EM
    a.EM.data<-cbind(a.EM,a.3PL.true)
    a.EM.recoded<-matrix(a.EM.data[complete.cases(a.EM.data), ],ncol=2)
    colnames(a.EM.recoded)<-c("a.EM","a.True")
    
    #dropping missing cases for b.EM
    b.EM.data<-cbind(b.EM,b.3PL.true)
    b.EM.recoded<-matrix(b.EM.data[complete.cases(b.EM.data), ],ncol=2)
    colnames(b.EM.recoded)<-c("b.EM","b.True")
    
    #dropping missing cases for theta.EM
    theta.EM.data<-cbind(theta.EM,theta.combined)
    theta.EM.recoded<-matrix(theta.EM.data[complete.cases(theta.EM.data), ],ncol=2)
    theta.EM.recoded <- theta.EM.recoded[!is.infinite(rowSums(theta.EM.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.EM.recoded)<-c("theta.EM","theta.True")
    
    
    #dropping missing cases for a.EM.imputed
    a.EM.imputed.data<-cbind(a.EM.imputed,a.3PL.true)
    a.EM.imputed.recoded<-matrix(a.EM.imputed.data[complete.cases(a.EM.imputed.data), ],ncol=2)
    colnames(a.EM.imputed.recoded)<-c("a.EM.imputed","a.True")
    
    #dropping missing cases for b.EM.imputed
    b.EM.imputed.data<-cbind(b.EM.imputed,b.3PL.true)
    b.EM.imputed.recoded<-matrix(b.EM.imputed.data[complete.cases(b.EM.imputed.data), ],ncol=2)
    colnames(b.EM.imputed.recoded)<-c("b.EM.imputed","b.True")
    
    #dropping missing cases for theta.EM.imputed
    theta.EM.imputed.data<-cbind(theta.EM.imputed,theta.combined)
    theta.EM.imputed.recoded<-matrix(theta.EM.imputed.data[complete.cases(theta.EM.imputed.data), ],ncol=2)
    theta.EM.imputed.recoded <- theta.EM.imputed.recoded[!is.infinite(rowSums(theta.EM.imputed.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.EM.imputed.recoded)<-c("theta.EM.imputed","theta.True")
    
    
    #dropping missing cases for a.HG
    a.HG.data<-cbind(a.HG,a.3PL.true)
    a.HG.recoded<-matrix(a.HG.data[complete.cases(a.HG.data), ],ncol=2)
    colnames(a.HG.recoded)<-c("a.HG","a.True")
    
    #dropping missing cases for b.HG
    b.HG.data<-cbind(b.HG,b.3PL.true)
    b.HG.recoded<-matrix(b.HG.data[complete.cases(b.HG.data), ],ncol=2)
    colnames(b.HG.recoded)<-c("b.HG","b.True")
    
    #dropping missing cases for theta.HG
    theta.HG.data<-cbind(theta.HG,theta.combined)
    theta.HG.recoded<-matrix(theta.HG.data[complete.cases(theta.HG.data), ],ncol=2)
    theta.HG.recoded <- theta.HG.recoded[!is.infinite(rowSums(theta.HG.recoded)),] #removing any cases with inifite theta estimates
    colnames(theta.HG.recoded)<-c("theta.HG","theta.True")
    
    
    
    ####################################
    # calculating overall descriptives #
    ####################################
    #mean bias
    results.condition[r,1]<-bias(a.ML.recoded[,"a.ML"],a.ML.recoded[,"a.True"],type="bias")
    results.condition[r,2]<-bias(a.penalized.recoded[,"a.penalized"],a.penalized.recoded[,"a.True"],type="bias")
    results.condition[r,3]<-bias(a.EM.recoded[,"a.EM"],a.EM.recoded[,"a.True"],type="bias")
    results.condition[r,4]<-bias(a.EM.imputed.recoded[,"a.EM.imputed"],a.EM.imputed.recoded[,"a.True"],type="bias")
    results.condition[r,5]<-bias(a.HG.recoded[,"a.HG"],a.HG.recoded[,"a.True"],type="bias")
    
    results.condition[r,6]<-bias(b.ML.recoded[,"b.ML"],b.ML.recoded[,"b.True"],type="bias")
    results.condition[r,7]<-bias(b.penalized.recoded[,"b.penalized"],b.penalized.recoded[,"b.True"],type="bias")
    results.condition[r,8]<-bias(b.EM.recoded[,"b.EM"],b.EM.recoded[,"b.True"],type="bias")
    results.condition[r,9]<-bias(b.EM.imputed.recoded[,"b.EM.imputed"],b.EM.imputed.recoded[,"b.True"],type="bias")
    results.condition[r,10]<-bias(b.HG.recoded[,"b.HG"],b.HG.recoded[,"b.True"],type="bias")
    
    
    results.condition[r,11]<-bias(theta.ML.recoded[,"theta.ML"],theta.ML.recoded[,"theta.True"],type="bias")
    results.condition[r,12]<-bias(theta.penalized.recoded[,"theta.penalized"],theta.penalized.recoded[,"theta.True"],type="bias")
    results.condition[r,13]<-bias(theta.EM.recoded[,"theta.EM"],theta.EM.recoded[,"theta.True"],type="bias")
    results.condition[r,14]<-bias(theta.EM.imputed.recoded[,"theta.EM.imputed"],theta.EM.imputed.recoded[,"theta.True"],type="bias")
    results.condition[r,15]<-bias(theta.HG.recoded[,"theta.HG"],theta.HG.recoded[,"theta.True"],type="bias")
    
    
    #results.condition[r,11]<-bias(as.vector(c.ML),as.vector(c.3PL.true),type="bias")
    #results.condition[r,12]<-bias(as.vector(c.penalized),as.vector(c.3PL.true),type="bias")
    #results.condition[r,13]<-bias(as.vector(c.EM),as.vector(c.3PL.true),type="bias")
    #results.condition[r,14]<-bias(as.vector(c.EM.imputed),as.vector(c.3PL.true),type="bias")
    #results.condition[r,15]<-bias(as.vector(c.HG),as.vector(c.3PL.true),type="bias")
    
  
    
    #mean se
    #results.condition[r,21]<- sqrt(sum((a.ML-mean(a.ML))^2)/I)
    #results.condition[r,22]<- sqrt(sum((a.penalized-mean(a.penalized))^2)/I)
    #results.condition[r,23]<- sqrt(sum((a.EM-mean(a.EM))^2)/I)
    #results.condition[r,24]<- sqrt(sum((a.EM.imputed-mean(a.EM.imputed))^2)/I)
    #results.condition[r,25]<- sqrt(sum((a.HG-mean(a.HG))^2)/I)
    
    #results.condition[r,26]<- sqrt(sum((b.ML-mean(b.ML))^2)/I)
    #results.condition[r,27]<- sqrt(sum((b.penalized-mean(b.penalized))^2)/I)
    #results.condition[r,28]<- sqrt(sum((b.EM-mean(b.EM))^2)/I)
    #results.condition[r,29]<- sqrt(sum((b.EM.imputed-mean(b.EM.imputed))^2)/I)
    #results.condition[r,30]<- sqrt(sum((b.HG-mean(b.HG))^2)/I)
    
    #results.condition[r,31]<- sqrt(sum((c.ML-mean(c.ML))^2)/I)
    #results.condition[r,32]<- sqrt(sum((c.penalized-mean(c.penalized))^2)/I)
    #results.condition[r,33]<- sqrt(sum((c.EM-mean(c.EM))^2)/I)
    #results.condition[r,34]<- sqrt(sum((c.EM.imputed-mean(c.EM.imputed))^2)/I)
    #results.condition[r,35]<- sqrt(sum((c.HG-mean(c.HG))^2)/I)
    
    #results.condition[r,36]<- sqrt(sum((theta.ML[1:N.unmot,]-mean(theta.unmot))^2)/I)
    #results.condition[r,37]<- sqrt(sum((theta.penalized[1:N.unmot,]-mean(theta.unmot))^2)/I)
    #results.condition[r,38]<- sqrt(sum((theta.EM[1:N.unmot,]-mean(theta.unmot))^2)/I)
    #results.condition[r,39]<- sqrt(sum((theta.EM.imputed[1:N.unmot,]-mean(theta.unmot))^2)/I)
    #results.condition[r,40]<- sqrt(sum((theta.HG[1:N.unmot,]-mean(theta.unmot))^2)/I)
    
    #RMSE
    
    
    #mean RMSE
    results.condition[r,16]<-RMSE(a.ML.recoded[,"a.ML"],a.ML.recoded[,"a.True"],type="RMSE")
    results.condition[r,17]<-RMSE(a.penalized.recoded[,"a.penalized"],a.penalized.recoded[,"a.True"],type="RMSE")
    results.condition[r,18]<-RMSE(a.EM.recoded[,"a.EM"],a.EM.recoded[,"a.True"],type="RMSE")
    results.condition[r,19]<-RMSE(a.EM.imputed.recoded[,"a.EM.imputed"],a.EM.imputed.recoded[,"a.True"],type="RMSE")
    results.condition[r,20]<-RMSE(a.HG.recoded[,"a.HG"],a.HG.recoded[,"a.True"],type="RMSE")
    
    results.condition[r,21]<-RMSE(b.ML.recoded[,"b.ML"],b.ML.recoded[,"b.True"],type="RMSE")
    results.condition[r,22]<-RMSE(b.penalized.recoded[,"b.penalized"],b.penalized.recoded[,"b.True"],type="RMSE")
    results.condition[r,23]<-RMSE(b.EM.recoded[,"b.EM"],b.EM.recoded[,"b.True"],type="RMSE")
    results.condition[r,24]<-RMSE(b.EM.imputed.recoded[,"b.EM.imputed"],b.EM.imputed.recoded[,"b.True"],type="RMSE")
    results.condition[r,25]<-RMSE(b.HG.recoded[,"b.HG"],b.HG.recoded[,"b.True"],type="RMSE")
    
    
    results.condition[r,26]<-RMSE(theta.ML.recoded[,"theta.ML"],theta.ML.recoded[,"theta.True"],type="RMSE")
    results.condition[r,27]<-RMSE(theta.penalized.recoded[,"theta.penalized"],theta.penalized.recoded[,"theta.True"],type="RMSE")
    results.condition[r,28]<-RMSE(theta.EM.recoded[,"theta.EM"],theta.EM.recoded[,"theta.True"],type="RMSE")
    results.condition[r,29]<-RMSE(theta.EM.imputed.recoded[,"theta.EM.imputed"],theta.EM.imputed.recoded[,"theta.True"],type="RMSE")
    results.condition[r,30]<-RMSE(theta.HG.recoded[,"theta.HG"],theta.HG.recoded[,"theta.True"],type="RMSE")
    
    
    
    
    #results.condition[r,51]<-RMSE(as.vector(c.ML),as.vector(c.3PL.true),type="RMSE")
    #results.condition[r,52]<-RMSE(as.vector(c.penalized),as.vector(c.3PL.true),type="RMSE")
    #results.condition[r,53]<-RMSE(as.vector(c.EM),as.vector(c.3PL.true),type="RMSE")
    #results.condition[r,54]<-RMSE(as.vector(c.EM.imputed),as.vector(c.3PL.true),type="RMSE")
    #results.condition[r,55]<-RMSE(as.vector(c.HG),as.vector(c.3PL.true),type="RMSE")
    
    
    
    
    #bias and RMSE for misclassified examinees only
    if(nrow(theta.ML.recoded)>0){
    results.condition[r,42]<-bias(theta.ML.recoded[1:N.over,"theta.ML"],theta.ML.recoded[1:N.over,"theta.True"],type="bias")
    results.condition[r,47]<-RMSE(theta.ML.recoded[1:N.over,"theta.ML"],theta.ML.recoded[1:N.over,"theta.True"],type="RMSE")
    } else if (nrow(theta.ML.recoded)==0){
      results.condition[r,42]<-NA
        results.condition[r,47]<-NA
    }
    
    if(nrow(theta.penalized.recoded)>0){
    results.condition[r,43]<-bias(theta.penalized.recoded[1:N.over,"theta.penalized"],theta.penalized.recoded[1:N.over,"theta.True"],type="bias")
    results.condition[r,48]<-RMSE(theta.penalized.recoded[1:N.over,"theta.penalized"],theta.penalized.recoded[1:N.over,"theta.True"],type="RMSE")
    } else if (nrow(theta.penalized.recoded)==0){
      results.condition[r,43]<-NA
      results.condition[r,48]<-NA
    }
    
    
    if(nrow(theta.EM.recoded)>0){
    results.condition[r,44]<-bias(theta.EM.recoded[1:N.over,"theta.EM"],theta.EM.recoded[1:N.over,"theta.True"],type="bias")
    results.condition[r,49]<-RMSE(theta.EM.recoded[1:N.over,"theta.EM"],theta.EM.recoded[1:N.over,"theta.True"],type="RMSE")
    } else if (nrow(theta.EM.recoded)==0){
      results.condition[r,44]<-NA
      results.condition[r,49]<-NA
    }
    
    
    if(nrow(theta.EM.imputed.recoded)>0){
    results.condition[r,45]<-bias(theta.EM.imputed.recoded[1:N.over,"theta.EM.imputed"],theta.EM.imputed.recoded[1:N.over,"theta.True"],type="bias")
    results.condition[r,50]<-RMSE(theta.EM.imputed.recoded[1:N.over,"theta.EM.imputed"],theta.EM.imputed.recoded[1:N.over,"theta.True"],type="RMSE")
    } else if (nrow(theta.EM.imputed.recoded)==0){
      results.condition[r,45]<-NA
      results.condition[r,50]<-NA
    }
    
    if(nrow(theta.HG.recoded)>0){
    
    results.condition[r,46]<-bias(theta.HG.recoded[1:N.over,"theta.HG"],theta.HG.recoded[1:N.over,"theta.True"],type="bias")
    results.condition[r,51]<-RMSE(theta.HG.recoded[1:N.over,"theta.HG"],theta.HG.recoded[1:N.over,"theta.True"],type="RMSE")
    
   } else if (nrow(theta.HG.recoded)==0){
      results.condition[r,46]<-NA
      results.condition[r,51]<-NA
    }
    
    
    #BIC VALUES
    if(results.condition[r,37]==1){ 
      results.condition[r,31]<-NA
      
    } else if (results.condition[r,37]==0){
      results.condition[r,31]<-anova(ML)$BIC
     }
    
    
    #Penalized scoring
    
    if(results.condition[r,38]==1){
      results.condition[r,32]<-NA
     
    } else if (results.condition[r,38]==0){
      results.condition[r,32]<-anova(penalized)$BIC
       }
    
    #EM 
    if(results.condition[r,39]==1){
      results.condition[r,33]<-NA
      
    } else if (results.condition[r,39]==0){
      results.condition[r,33]<-anova(EM)$BIC
     
    }
  
    #EM-IMPUTED
    if(results.condition[r,40]==1|(is.nan(results.condition[r,40])==TRUE)){
      results.condition[r,34]<-NA
     
    } else {
      results.condition[r,34]<-anova(EM.imputed)$BIC
    
    }
    
    
    #HG scoring
    if(results.condition[r,41]==1){
      results.condition[r,35]<-NA
     
      } else if (results.condition[r,41]==0){
        results.condition[r,35]<-anova(HG)$BIC
        }
    
  } #closes the rep loop
  
  #place descriptive results for each condition into an overall matrix that will be printed out
  overall.results[count.con,1:51]<-colMeans(results.condition[,1:51],na.rm = TRUE)
  
  print(count.con) #printing number of condition
  #count.con = count.con + 1
  
  
} #closes the condition loop; the parenthesis around the bracket will calculate the run time

##### WRITING OUT RESULTS #####
overall.results.output<-data.frame(cbind(matrix(con[,1],ncol=1),
                                         overall.results))

names(overall.results.output)[names(overall.results.output) == "V1"]<-"Overclassification Percentage"
               

write.csv(overall.results.output,"Overall Results - Overclassification Conditions 2022_03_31.csv",row.names = FALSE)
