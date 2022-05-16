##############################################################
#
#       RESPONSE-TIME THRESHOLD PROCEDURE COMPARISONS        #
#
##############################################################

#clearing memory
rm(list=ls())

##################### Setting working directory ####################

setwd("G:\\My Drive\\Research\\Projects\\2021-2022\\Two-stage Scoring\\HG-MI")


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
Unmot.perc<-.28 #percentage of sample that engage in RG
N.unmot<-N*Unmot.perc #size of unmotivated sample
N.mot<-N*(1-Unmot.perc) #size of motivated sample
reps<-100 #number of reps for study
reps.EM.imputed<-100 #number of reps for EM-imputation procedure
reps.HG.imputed <- 100

#################### creating conditions matrix ####################


# level codings will be specified in if else statements


#RG PATTERN
#there are two levels: 1 = difficulty based; 2 = progressive. (for v2, reviewers wanted me to drop this factor) 
rg.type <- c("Difficulty","Idiosyncratic")

#PERCENTAGE OF RG RESPONSES FOR UNMOTIVATED SIMULEES
#there are three levels: 10%; 20% 
rg.percent <- c(.1,.2)

#MISCLASSIFICATION PERCENTAGE
misclassify.percent <-0.6

#TEST DIFFICULTY
#there are two levels: easy [mean b = -1] and moderate [mean b = 0]
test.difficulty <- c("Easy", "Moderate")


# We need to get a conditions matrix for a fully crossed design. In order to do this
# efficiently we can use the expand grid function which will give us a 48x5 matrix
# 2x2x4 = 16 levels crossed with the 4 variables: impact, NER relationship with ability, NER type, 
#percent unmotivated simulees in focal group, percent NER for unmotivated simulees in focal group
con <- expand.grid(rg.type,rg.percent,misclassify.percent, test.difficulty)


#reading in item parameters
item.3PL.pars.easy<-as.matrix(read.csv("Generating 3pl item parameters (50 items, easy).csv"))
item.3PL.pars.mod <- as.matrix(read.csv("Generating 3pl item parameters (50 items, moderate).csv"))
I<-50 #number of items
a.3PL.true.easy<-as.matrix(item.3PL.pars.easy[,1],ncol=1)
b.3PL.true.easy<-as.matrix(item.3PL.pars.easy[,2],ncol=1)
c.3PL.true.easy<-as.matrix(item.3PL.pars.easy[,3],ncol=1)

a.3PL.true.mod<-as.matrix(item.3PL.pars.mod[,1],ncol=1)
b.3PL.true.mod<-as.matrix(item.3PL.pars.mod[,2],ncol=1)
c.3PL.true.mod<-as.matrix(item.3PL.pars.mod[,3],ncol=1)


#creating matrix for overall results by condition
overall.results <- matrix(-999, nrow = nrow(con), ncol = 20)
colnames(overall.results) <- c(
  "ML.all","EM.all","EM-I.all","HG.all","HG-MI.all",
  "ML.low", "EM.low", "EM-I.low", "HG.low", "HG-MI.low",
  "ML.medium", "EM.medium", "EM-I.medium", "HG.medium", "HG-MI.medium",
  "ML.high", "EM.high","EM-I.high","HG.high", "HG-MI.medium"
)


#condition loop starts here


for(count.con in 1:nrow(con)) { 
  
  
  results.condition<-matrix(-999,nrow=reps,ncol=20) #creating a matrix to place results in by condition, which will be used for the overall descriptive results
  
  
  
  ######### REPLICATION LOOP STARTS HERE
  for (r in 1:reps){
    #setting seed
    set.seed(123+(r*count.con)+r) #rep to this counter
    ##################### Sampling ability parameters ##############
    
    theta.dif<- -0.5 #difference between effortful and noneffortful simulees
    
    theta.unmot<-matrix(rnorm(N.unmot,theta.dif,1),ncol=1) #sampling ability of unmotivated is -0.5 sd lower than motivated
    theta.mot<-matrix(rnorm(N.mot,0,1),ncol=1) #motivated ability, standard normal
    theta.combined<-rbind(theta.unmot,theta.mot)
    
    #breakdown unmotivated simulees into low, medium, and high ability subgroups by their true theta
    if(con[count.con,1]=="Difficulty"){
      id.high <- which(theta.unmot >0 | theta.unmot <0.5)
      id.medium <- which(theta.unmot >=-1 | theta.unmot <=0)
      id.low <- which(theta.unmot < -1)
    }else if(con[count.con,1]=="Idiosyncratic"){
      id.high <- which(theta.unmot >0.5)
      id.medium <- which(theta.unmot >=-0.5 | theta.unmot <=0.5)
      id.low <- which(theta.unmot < -0.5)
    }
    
    ##################### Generating response probabilities ##############
    
    #creating function to compute 3pl 
    PL3<-function(thetas,a.MC,b.MC,c.MC,N,I.MC){
      alphas.betas <-matrix(a.MC * b.MC,ncol=1)
      alphas.thetas <-matrix(thetas,ncol=1) %*% matrix(a.MC,nrow=1)
      logit <-alphas.thetas - t(matrix(alphas.betas,I.MC,N))
      c.par <- t(matrix(c.MC,I.MC,N))
      probs <- (c.par+((1-c.par)*(1/(1+exp(-logit)))))
      return(probs)}
    
    if(con[count.con,4]=="Easy"){
      Unmot.prob<-PL3(theta.unmot,a.3PL.true.easy,b.3PL.true.easy,c.3PL.true.easy,N.unmot,I) #generating unmotivated response probabilities
      Mot.prob<-PL3(theta.mot,a.3PL.true.easy,b.3PL.true.easy,c.3PL.true.easy,N.mot,I) #generating motivated response probabilities
    }else if(con[count.con,4]=="Moderate"){
      Unmot.prob<-PL3(theta.unmot,a.3PL.true.mod,b.3PL.true.mod,c.3PL.true.mod,N.unmot,I) #generating unmotivated response probabilities
      Mot.prob<-PL3(theta.mot,a.3PL.true.mod,b.3PL.true.mod,c.3PL.true.mod,N.mot,I) #generating motivated response probabilities
    }
    
    
    ######################### SAMPLING RANDOM NUMBER OF RGs FOR EACH UNMOTIVATED SIMULEE #######################
    #this is the random number of RGs for condition with 10% RGs
    
    aa<-N*I #total number of item responses (N x I)
    bb<-aa*con[count.con,2] #number of RG responses to sample
    prob.sample<-1/N.unmot #probability of sampling each simulee
    samples<-1 #number of samples of random numbers to pull 
    
    RG.1<-matrix(rmultinom(samples, size = bb,prob=c(rep(prob.sample,N))),ncol=samples,nrow=N)
    #sampling number of RG responses per simulee
    
    
    
    ########################  CREATING RG ##############################
    
    
    prob.unmot.NA<-matrix(-999,nrow=N.unmot,ncol=I)
    prob.unmot.zero<-matrix(-999,nrow=N.unmot,ncol=I)
    
    RG.matrix.mot<-matrix(0,nrow=N.mot,ncol=I) #CREATING RG Matrix for Motivated Sample (all motivated)
    
    RG.matrix.unmot<-matrix(0,nrow=N.unmot,ncol=I) #CREATING RG Matrix for Unmotivated Sample 
    
    for(p in 1:N.unmot){
      
      
      ############## SELECTING ITEMS TO RECEIVE RG ##################
      if(con[count.con,1]=="Difficulty"){
        
        #rank ordering the true probability of success from lowest to highest; ties are randomly ordered 
        rank.prob.success.guesser<-order(Unmot.prob[p,],decreasing=TRUE) #order in descending order the true probability of success from lowest to highest; ties are randomly ordered 
        item.id<-rank.prob.success.guesser[c(1:RG.1[p,])]
        
      } else if(con[count.con,1]=="Idiosyncratic"){
        
        #idiosyncratic RG
        item.id<-sample(matrix(1:I,ncol=1),RG.1[p,]) #sampling a random item number#
      }
      
      ######################### UPDATING GUESSING PROBABILITIES & RG MATRIX ############################
      Unmot.prob[p,item.id]<-.25 #replacing original probability with .25 for the random item number #
      
      #replacing 1 with 0 to represent RG responses
      RG.matrix.unmot[p,item.id]<-1
      
      
      ############################## ADDING RG MISCLASSIFICATIONS ###############################
      
      if(con[count.con,3]==0){
        RG.matrix.unmot[p,]<-RG.matrix.unmot[p,] #NO MISCLASSIFICATION CONDITIONS
      } else if (con[count.con,3]!=0){ #MISCLASSIFICATION CONDITIONS
        
        #BELOW REPRESENTS UNDERCLASSIFICATION
        #HERE I NEED TO SAMPLE ITEM.IDs TO MISCLASSIFY RGs IDENTIFY ITEMS
        #2 STEP PROCESS: (1) IDENTIFY HOW MANY RGs TO MISCLASSIFY; (2) RANDOMLY MISCLASSIFY
        #the line below is random sampling items that have been identified as RGs
        #this code is written for one unmotivated simulee at a time
        
        
        #UNDERCLASSIFICATION CONDITION
        if(round((con[count.con,3]*nrow(matrix(item.id,ncol=1))))>0){ #if the number of underclassified RG responses is greater than 0
          
          #random sampling number of RG responses to underclassify for each simulee
          RG.U.item.id<-sample(item.id,round((con[count.con,3]*nrow(matrix(item.id,ncol=1))))) 
          
          #replacing 1 with 0 to represent underclassified RG responses
          
          RG.matrix.unmot[p,RG.U.item.id]<-0
        }else{
          RG.matrix.unmot[p,]<-RG.matrix.unmot[p,] #if zero misclassifications are specified, then true RG classification is unchanged
        }
      }
      #print(p)
    }
    ############## NEED TO RECODE PROBABILITIES BASED ON UPDATED RG MATRIX
    for(p in 1:N.unmot){
      for (i in 1:I){
        ifelse(RG.matrix.unmot[p,i]==1,prob.unmot.NA[p,i]<-NA,prob.unmot.NA[p,i]<-Unmot.prob[p,i])
        ifelse(RG.matrix.unmot[p,i]==1,prob.unmot.zero[p,i]<-0,prob.unmot.zero[p,i]<-Unmot.prob[p,i])
        
      }
    }
    
    
    RG.combined<-rbind(RG.matrix.unmot,RG.matrix.mot) #combining RG classifications
    
    RG.combined[N,]<-1 #adding in case here to make sure that I can estimate model (need at least two response options for each item)
    
    #combining response probabilities for unmotivated and motivated samples
    probs.combined.ML<-rbind(Unmot.prob,Mot.prob)
    probs.combined.NA<-rbind(prob.unmot.NA,Mot.prob)
    probs.combined.zero<-rbind(prob.unmot.zero,Mot.prob)
    
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
    #ipars.ML<-coef(ML,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
    
    
   # results.condition[r,37]<-ifelse(extract.mirt(ML,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(extract.mirt(ML,'converged')==TRUE){ 
      theta.ML<-as.matrix(fscores(ML,method='ML',max_theta=4))
      
    } else {
      theta.ML<-rep(NA,N)
    }
    
    
   
    ###### EM scoring ######
    EM<-mirt(IR.combined.NA, 1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000))
   # ipars.EM<-coef(EM,IRTpars=TRUE,as.data.frame=TRUE) #need to read in item parameters
    
    results.condition[r,39]<-ifelse(extract.mirt(EM,'converged')==TRUE,0,1) #if model failed to converge, give 1; otherwise, 0
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(extract.mirt(EM,'converged')==TRUE){ 
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
    
    
    theta.EM.imputed<-matrix(-999,nrow=N,ncol=reps.EM.imputed)
    
    
    if(extract.mirt(EM,'converged')==FALSE){ 
      theta.EM.imputed<-matrix(rep(matrix(rep(NA,N),ncol=1),reps.EM.imputed),ncol=reps.EM.imputed)
      
    } else if (extract.mirt(EM,'converged')==TRUE){
      
      for(r.EM.impute in 1:reps.EM.imputed){
        
        rnd.unif.impute <- matrix(runif(N*I, 0, 1), nrow = N, ncol = I)
        
        #transforming probs to 0/1 for when noneffortful responses are in the data matrix
        IR.EM.imputed <- (ifelse(imputed.probabilities > rnd.unif.impute,1,0)) #Coding responses as 0/1 for IRT software#
        colnames(IR.EM.imputed)<-rep(1:I)
        
        EM.imputed<-mirt(IR.EM.imputed,  1, itemtype ='2PL',guess=.25, TOL = .0001, technical = list(NCYCLES = 10000))
        
        
        #if model converged, take parameter estimates. if not, impute missing values
        if(extract.mirt(EM.imputed,'converged')==TRUE){ 
          theta.EM.imputed[,r.EM.impute]<-as.matrix(fscores(EM.imputed,method='ML',max_theta=4))
          
        } else if (extract.mirt(EM.imputed,'converged')==FALSE){
          theta.EM.imputed[,r.EM.impute]<-rep(NA,N) 
        }
        
      }
      
    } #closes if/else statement above
    
    
    theta.EM.imputed<-matrix(rowMeans(theta.EM.imputed,na.rm=TRUE),ncol=1)
    
    
    
   
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
    
    
    #if model converged, take parameter estimates. if not, impute missing values
    if(extract.mirt(HG,'converged')==TRUE){ 
      #need to transform parameter estimates onto theta scale
      a.HG<-ipars.HG[c(seq(1,250,5))] 
      b.HG<- -(ipars.HG[c(seq(3,250,5))])/a.HG #I am not sure if there is a problem; trying to convert cfa parameters to irt b=threshold/slope
      theta.HG<-matrix(fscores(HG,method='ML',max_theta=4)[,1],ncol=1) #setting upper and lower limits of theta at 4 and -4
      
    } else {
      theta.HG<-rep(NA,N) #setting upper and lower limits of theta at 4 and -4
      
    }
    
    
    ##### HG-MI scoring #####
    
    # parameters obtatined from HG scoring: a.HG, b.HG, d.HG, theta.HG
    
    #imputing probabilities
    HG.imputed.probabilities<-PL3(theta.HG,a.HG,b.HG,c.HG,N,I) #based on the theta and item parameter estimates from HG scoring
    
    
    theta.HG.imputed<-matrix(-999,nrow=N,ncol=reps.HG.imputed)
    
    
    if(extract.mirt(HG,'converged')==FALSE){
      
      theta.HG.imputed<-matrix(rep(matrix(rep(NA,N),ncol=1),reps.HG.imputed),ncol=reps.HG.imputed)
      
    } else if (extract.mirt(HG,'converged')==TRUE){
      
      for(r.HG.impute in 1:reps.HG.imputed){
        
        rnd.unif.impute <- matrix(runif(N*I, 0, 1), nrow = N, ncol = I)
        
        #transforming probs to 0/1 for when noneffortful responses are in the data matrix
        IR.HG.imputed <- (ifelse(HG.imputed.probabilities > rnd.unif.impute,1,0)) #Coding responses as 0/1 for IRT software#
        colnames(IR.HG.imputed)<-rep(1:I)
        
        #fix item parameters
        sv <-  mirt(IR.HG.imputed, 1, itemtype = "2PL", guess = .25, TOL = .0001, technical = list(NCYCLES = 10000), pars = 'values')
        sv$value[sv$name == 'a1'] <- a.HG
        sv$value[sv$name =='d'] <- d.HG
        sv$est <- FALSE
        
        HG.imputed<-mirt(IR.HG.imputed, 1, itemtype = "2PL", guess = .25, TOL = .0001, technical = list(NCYCLES = 10000), pars = sv)
        
        #if model converged, take parameter estimates. if not, impute missing values
        if(extract.mirt(HG.imputed,'converged')==TRUE){ 
          theta.HG.imputed[,r.HG.impute]<-as.matrix(fscores(HG.imputed,method='ML',max_theta=4))
          
        } else if (extract.mirt(HG.imputed,'converged')==FALSE){
          theta.HG.imputed[,r.HG.impute]<-rep(NA,N) 
          
        }
        
      }
      
    } #closes if/else statement above
    
    
    theta.HG.imputed<-matrix(rowMeans(theta.HG.imputed,na.rm=TRUE),ncol=1)  
    
    
    ###################################
    #                                 #
    # CALCULATE DEPENDENT VARIABLES   #
    #                                 #
    ###################################
    #have to drop missing data to run simdesign functions
    
    #dropping missing cases for theta.ML
    theta.ML.data<-cbind(theta.ML,theta.combined)
    theta.ML.recoded<-matrix(theta.ML.data[complete.cases(theta.ML.data), ],ncol=2)
    theta.ML.recoded <- theta.ML.recoded[!is.infinite(rowSums(theta.ML.recoded)),] #removing any cases with inifite theta estimates
    theta.ML.recoded <- theta.ML.recoded[1:N.unmot]
    colnames(theta.ML.recoded)<-c("theta.ML","theta.True")
    
    #dropping missing cases for theta.EM
    theta.EM.data<-cbind(theta.EM,theta.combined)
    theta.EM.recoded<-matrix(theta.EM.data[complete.cases(theta.EM.data), ],ncol=2)
    theta.EM.recoded <- theta.EM.recoded[!is.infinite(rowSums(theta.EM.recoded)),] #removing any cases with inifite theta estimates
    theta.EM.recoded <- theta.EM.recoded[1:N.unmot]
    colnames(theta.EM.recoded)<-c("theta.EM","theta.True")
    
    #dropping missing cases for theta.EM.imputed
    theta.EM.imputed.data<-cbind(theta.EM.imputed,theta.combined)
    theta.EM.imputed.recoded<-matrix(theta.EM.imputed.data[complete.cases(theta.EM.imputed.data), ],ncol=2)
    theta.EM.imputed.recoded <- theta.EM.imputed.recoded[!is.infinite(rowSums(theta.EM.imputed.recoded)),] #removing any cases with inifite theta estimates
    theta.EM.imputed.recoded <- theta.EM.imputed.recoded[1:N.unmot]
    colnames(theta.EM.imputed.recoded)<-c("theta.EM.imputed","theta.True")
    
    #dropping missing cases for theta.HG
    theta.HG.data<-cbind(theta.HG,theta.combined)
    theta.HG.recoded<-matrix(theta.HG.data[complete.cases(theta.HG.data), ],ncol=2)
    theta.HG.recoded <- theta.HG.recoded[!is.infinite(rowSums(theta.HG.recoded)),] #removing any cases with inifite theta estimates
    theta.HG.recoded <- theta.HG.recoded[1:N.unmot]
    colnames(theta.HG.recoded)<-c("theta.HG","theta.True")
    
    #dropping missing cases for theta.HG.imputed
    theta.HG.imputed.data<-cbind(theta.HG.imputed,theta.combined)
    theta.HG.imputed.recoded<-matrix(theta.HG.imputed.data[complete.cases(theta.HG.imputed.data), ],ncol=2)
    theta.HG.imputed.recoded <- theta.HG.imputed.recoded[!is.infinite(rowSums(theta.HG.imputed.recoded)),] #removing any cases with inifite theta estimates
    theta.HG.imputed.recoded <- theta.HG.imputed.recoded[1:N.unmot]
    colnames(theta.HG.imputed.recoded)<-c("theta.HG.Imputed","theta.True")
    
    
    
    
    ################ CALCULATING DEPENDENT VARIABLES ###########
    
    #ML
    if(nrow(theta.ML.recoded)>0){
      results.condition[r,1]<-RMSE(theta.ML.recoded[,"theta.ML"],theta.ML.recoded[,"theta.True"],type="RMSE")#all
      results.condition[r,6]<-RMSE(theta.ML.recoded[id.low,"theta.ML"],theta.ML.recoded[id.low,"theta.True"],type="RMSE")#low
      results.condition[r,11]<-RMSE(theta.ML.recoded[id.medium,"theta.ML"],theta.ML.recoded[id.medium,"theta.True"],type="RMSE")#medium
      results.condition[r,16]<-RMSE(theta.ML.recoded[id.high,"theta.ML"],theta.ML.recoded[id.high,"theta.True"],type="RMSE")#high
      
    } else {
      results.condition[r,1]<-NA
      esults.condition[r,6]<-NA
      results.condition[r,11]<-NA
      results.condition[r,16]<-NA
    }
    

    
    #EM
    if(nrow(theta.EM.recoded)>0){
      results.condition[r,2]<-RMSE(theta.EM.recoded[,"theta.ML"],theta.EM.recoded[,"theta.True"],type="RMSE")#all
      results.condition[r,7]<-RMSE(theta.EM.recoded[id.low,"theta.ML"],theta.EM.recoded[id.low,"theta.True"],type="RMSE")#low
      results.condition[r,12]<-RMSE(theta.EM.recoded[id.medium,"theta.ML"],theta.EM.recoded[id.medium,"theta.True"],type="RMSE")#medium
      results.condition[r,17]<-RMSE(theta.EM.recoded[id.high,"theta.ML"],theta.EM.recoded[id.high,"theta.True"],type="RMSE")#high
      
    } else {
      results.condition[r,2]<-NA
      esults.condition[r,7]<-NA
      results.condition[r,12]<-NA
      results.condition[r,17]<-NA
    }
    
    
    #EM Imputation
    if(nrow(theta.EM.imputed.recoded)>0){
      results.condition[r,3]<-RMSE(theta.EM.imputed.recoded[,"theta.ML"],theta.EM.imputed.recoded[,"theta.True"],type="RMSE")#all
      results.condition[r,8]<-RMSE(theta.EM.imputed.recoded[id.low,"theta.ML"],theta.EM.imputed.recoded[id.low,"theta.True"],type="RMSE")#low
      results.condition[r,13]<-RMSE(theta.EM.imputed.recoded[id.medium,"theta.ML"],theta.EM.imputed.recoded[id.medium,"theta.True"],type="RMSE")#medium
      results.condition[r,18]<-RMSE(theta.EM.imputed.recoded[id.high,"theta.ML"],theta.EM.imputed.recoded[id.high,"theta.True"],type="RMSE")#high
      
    } else {
      results.condition[r,3]<-NA
      esults.condition[r,8]<-NA
      results.condition[r,13]<-NA
      results.condition[r,18]<-NA
    }
    
    
    
    #HG
    if(nrow(theta.HG.recoded)>0){
      results.condition[r,4]<-RMSE(theta.HG.recoded[,"theta.ML"],theta.HG.recoded[,"theta.True"],type="RMSE")#all
      results.condition[r,19]<-RMSE(theta.HG.recoded[id.low,"theta.ML"],theta.HG.recoded[id.low,"theta.True"],type="RMSE")#low
      results.condition[r,14]<-RMSE(theta.HG.recoded[id.medium,"theta.ML"],theta.HG.recoded[id.medium,"theta.True"],type="RMSE")#medium
      results.condition[r,19]<-RMSE(theta.HG.recoded[id.high,"theta.ML"],theta.HG.recoded[id.high,"theta.True"],type="RMSE")#high
      
    } else {
      results.condition[r,4]<-NA
      esults.condition[r,9]<-NA
      results.condition[r,14]<-NA
      results.condition[r,19]<-NA
    }
      
   
     
    
    
  } #closes the rep loop
  
  #place descriptive results for each condition into an overall matrix that will be printed out
  overall.results[count.con,1:20]<-colMeans(results.condition[,1:20],na.rm = TRUE)
  
  print(count.con) #printing number of condition
  #count.con = count.con + 1
  
  
} #closes the condition loop; the parenthesis around the bracket will calculate the run time

##### WRITING OUT RESULTS #####
overall.results.output<-data.frame(cbind(matrix(con[,1],ncol=1),
                                         matrix(con[,2],ncol=1),
                                         matrix(con[,3],ncol=1),
                                         matrix(con[,4],ncol = 1),
                                         overall.results))

names(overall.results.output)[names(overall.results.output) == "V1"]<-"RG Type"
names(overall.results.output)[names(overall.results.output) == "V2"]<-"RG Percent"                 
names(overall.results.output)[names(overall.results.output) == "V3"]<-"Misclassify Percent"    
names(overall.results.output)[names(overall.results.output) == "V4"]<-"Test Difficulty"

write.csv(overall.results.output,"Overall Results - 60% Overclassifications 2022_05_16.csv",row.names = FALSE)
