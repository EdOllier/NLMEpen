library(Rcpp)
library(RcppArmadillo)

sourceCpp("RCPP_SAPG_GAUSS_FAMILY_StatEx_ModifiedCholesky.cpp")


###################################
### PARTICLE SWARM OPTIMIZATION ###
###################################


PSO_PARALLEL_SAPG_ADAGRAD_R <- function(OBS,OBS_TIME,Design_Ind,MU_INIT,DELTA_INIT,GAMMA_INIT,ERROR_INIT,LAMBDA_INIT,LAMBDA_MAX,N_ITER_PerSEQ,AdaptWEIGHTS,N_ITER_PSO,CONTROL=list(N2_ITER,N_MCMC , beta, delta_0, Nbr_Iter_Prox, step_adagrad)){
    
    library(MASS)
    library(parallel)
    
    NbrSubject = dim(OBS)[1]
    NbrTimes = dim(OBS)[2]
    NbrParam = dim(GAMMA_INIT)[1]
    NbrCov = length(MU_INIT) - NbrParam
    NbrParamCov = length(MU_INIT)
    
    N_PARTICLE = dim(LAMBDA_INIT)[1]
    N_ITER_FIRST = 500
    N_ITER_PerSEQ_R = 250
    Nbr_Iter_Prox_R = 20
    ### PARTICLE SWARM OPTIMIZATION PARAMETER ###
    Velocity = matrix(0,N_PARTICLE,2)
    Velocity[,1] = runif(N_PARTICLE, min = -LAMBDA_MAX[1]*0.01, max = LAMBDA_MAX[1]*0.01)
    Velocity[,2] = runif(N_PARTICLE, min = -LAMBDA_MAX[2]*0.01, max = LAMBDA_MAX[2]*0.01)
    
    Vmax = LAMBDA_MAX*0.2
    
    INIERTIA_MIN = 0.9
    INIERTIA_MAX = 0.2
    C_AMP_1 = 2
    C_AMP_2 = 2
    
    ### SPARSITY PARAMETER ###
    LAMBDA = LAMBDA_INIT
    Nbr_Iter_Prox = CONTROL$Nbr_Iter_Prox
    
    AdaptWEIGHTS_MU = AdaptWEIGHTS[[1]]
    AdaptWEIGHTS_COV = AdaptWEIGHTS[[2]]
    
    ### PARAM GRADIENT DESCENT ###
    N_ITER_N = CONTROL$N_ITER_N
    StepAdapt = CONTROL$StepAdapt
    eps_cst = 1e-8
    step_adagrad = CONTROL$step_adagrad
    
    ### SOCHASTIC APPROXIMATION ###
    C_SA = CONTROL$delta_0
    beta = CONTROL$beta
    
    ### MCMC ###
    NbrIterN2 = CONTROL$N2_ITER
    N_ITER_MCMC = CONTROL$N_MCMC
    
    ### INITIALISATION ###
    STAT_EX_INITind0 =list(0*MU_PI_INIT, lapply(1:NbrSubject, function(x) matrix(0, nrow=NbrParam, ncol=NbrParam)) ,matrix(0,NbrParam,NbrParam),matrix(0,NbrParam*NbrParam,NbrParam*NbrParam),0,0)
    STAT_EX_INIT = lapply(1:N_PARTICLE, FUN= function(x){ list(0*MU_PI_INIT, lapply(1:NbrSubject, function(x) matrix(0, nrow=NbrParam, ncol=NbrParam)) ,matrix(0,NbrParam,NbrParam),matrix(0,NbrParam*NbrParam,NbrParam*NbrParam),0,0)} )
    PHI_INIT = lapply(1:N_PARTICLE, FUN= function(x){MU_PI_INIT} )
    SD_STEP_INIT = lapply(1:N_PARTICLE, FUN= function(x){rep(1,NbrParam)} )
    G_SQUARED_INIT = lapply(1:N_PARTICLE, FUN= function(x){ list(as.matrix(rep(0,NbrParamCov)), matrix(0,NbParam,NbParam),0) } )
    N_ITER_PREV =  rep(0,N_PARTICLE)
    
    DELTA_INIT_PART = lapply(1:N_PARTICLE, FUN= function(x){DELTA_INIT})
    
    MU_P = lapply(1:N_PARTICLE, FUN= function(x){MU_INIT})
    DELTA_P = lapply(1:N_PARTICLE, FUN= function(x){DELTA_INIT})
    GAMMA_P = lapply(1:N_PARTICLE, FUN= function(x){GAMMA_INIT})
    ERROR_P = lapply(1:N_PARTICLE, FUN= function(x){ERROR_INIT})
    
    MU_P_R = lapply(1:N_PARTICLE, FUN= function(x){MU_INIT})
    DELTA_P_R = lapply(1:N_PARTICLE, FUN= function(x){DELTA_INIT})
    GAMMA_P_R = lapply(1:N_PARTICLE, FUN= function(x){GAMMA_INIT})
    ERROR_P_R = lapply(1:N_PARTICLE, FUN= function(x){ERROR_INIT})
    
    MU_P_R_SEQ = list(); MU_P_R_SEQ[[1]] = MU_P_R
    DELTA_P_R_SEQ = list(); DELTA_P_R_SEQ[[1]] = DELTA_P_R
    GAMMA_P_R_SEQ = list(); GAMMA_P_R_SEQ[[1]] = GAMMA_P_R
    ERROR_P_R_SEQ = list(); ERROR_P_R_SEQ[[1]] = ERROR_P_R
    
    MU_P_SEQ = list(); MU_P_SEQ[[1]] = MU_P
    DELTA_P_SEQ = list(); DELTA_P_SEQ[[1]] = DELTA_P
    GAMMA_P_SEQ = list(); GAMMA_P_SEQ[[1]] = GAMMA_P
    ERROR_P_SEQ = list(); ERROR_P_SEQ[[1]] = ERROR_P
    
    EBIC_SEQ = matrix(0,N_ITER_PSO,N_PARTICLE)
    EBIC_R_SEQ = matrix(0,N_ITER_PSO,N_PARTICLE)
    
    PersonalBest = matrix(0,N_PARTICLE,2)
    GlobalBest = c(0,0)
    
    LAMBDA_SEQ_PSO = list()
    LAMBDA_SEQ_PSO[[1]] = LAMBDA
    
    PersonalBest_VALUE = rep(10000000,N_PARTICLE);
    GlobalBest_VALUE = 10000000
    
    MU_G_SEQ = list();
    DELTA_G_SEQ = list();
    GAMMA_G_SEQ = list();
    ERROR_G_SEQ = list();
    
    #########################
    ### Boucle Principale ###
    #########################
    
    RUN_SAPG <- function(x){ set.seed(123456);return(SAPG_ADAGRAD_OneRun(OBS,OBS_TIME, Design_Ind, MU_P[[x]]*c(rep(1,NbrParam),rep(0,NbrCov)),GAMMA_INIT,DELTA_INIT_PART[[x]],ERROR_P[[x]], list(LAMBDA[x,1],LAMBDA[x,2]), N_ITER_PerSEQ + (k==1)*N_ITER_FIRST, AdaptWEIGHTS, NbrIterN2, N_ITER_MCMC, beta , C_SA, Nbr_Iter_Prox, step_adagrad, STAT_EX_INIT[[x]], PHI_INIT[[x]], SD_STEP_INIT[[x]], G_SQUARED_INIT[[x]], N_ITER_PREV[x] ) ) }
    
    
    RUN_SAPG_RELAX <- function(x){ set.seed(123456);
        
        PHI_RELAX = 0
        
        AdaptWEIGHTS_R = list()
        AdaptWEIGHTS_R[[1]] = as.matrix(c(rep(0,NbParam),rep(1,NbrCOV*NbParam))*(1/abs(SAPG_OneRun[[x]][[1]])))
        AdaptWEIGHTS_R[[1]][which(!is.infinite(AdaptWEIGHTS_R[[1]]))] = PHI_RELAX*LAMBDA[x,1]
        AdaptWEIGHTS_R[[1]][which(is.infinite(AdaptWEIGHTS_R[[1]]))] = 1e4
        
        AdaptWEIGHTS_R[[2]] = (1/abs(SAPG_OneRun[[x]][[22]]))
        diag(AdaptWEIGHTS_R[[2]])=0
        AdaptWEIGHTS_R[[2]][which(!is.infinite(AdaptWEIGHTS_R[[2]]))] = PHI_RELAX*LAMBDA[x,2]
        AdaptWEIGHTS_R[[2]][which(is.infinite(AdaptWEIGHTS_R[[2]]))] = 1e4
        
        
        return(SAPG_ADAGRAD_OneRun(OBS,OBS_TIME, Design_Ind, SAPG_OneRun[[x]][[1]]*c(rep(1,NbrParam),rep(1,NbrCov)),SAPG_OneRun[[x]][[22]],SAPG_OneRun[[x]][[23]],SAPG_OneRun[[x]][[3]], list(1,1), N_ITER_PerSEQ_R , AdaptWEIGHTS_R, NbrIterN2, N_ITER_MCMC, beta, C_SA, Nbr_Iter_Prox_R, step_adagrad, SAPG_OneRun[[x]][[7]], SAPG_OneRun[[x]][[7]][[1]], SD_STEP_INIT[[x]], G_SQUARED_INIT[[x]], 0 ) )
        
    }
    
    for (k in 1:N_ITER_PSO){
        set.seed(k*100)
        RR = runif(2)
        INIERTIA = (INIERTIA_MIN - INIERTIA_MAX)*(N_ITER_PSO - k)/N_ITER_PSO + INIERTIA_MAX*4*RR*(1-RR)
        #browser()
        ### Evaluation des LAMBDA ###
        #browser()
        SAPG_OneRun <-mclapply(1:N_PARTICLE,RUN_SAPG,mc.cores=3)
        
        SAPG_OneRunR <-mclapply(1:N_PARTICLE,RUN_SAPG_RELAX,mc.cores=3)
        
        for (i in 1:N_PARTICLE){
            
            MU_P[[i]] = SAPG_OneRun[[i]][[1]]
            GAMMA_P[[i]] = SAPG_OneRun[[i]][[22]]
            DELTA_P[[i]] = SAPG_OneRun[[i]][[23]]
            ERROR_P[[i]] = SAPG_OneRun[[i]][[3]]
            
            MU_P_R[[i]] = SAPG_OneRunR[[i]][[1]]
            GAMMA_P_R[[i]] = SAPG_OneRunR[[i]][[22]]
            DELTA_P_R[[i]] = SAPG_OneRunR[[i]][[23]]
            ERROR_P_R[[i]] = SAPG_OneRunR[[i]][[3]]
            
            EBIC_SEQ[k,i] = SAPG_OneRun[[i]][[15]]
            EBIC_R_SEQ[k,i] = SAPG_OneRunR[[i]][[15]]
            
            ### Mise à jour des paramètres d'initialisation ###
            STAT_EX_INIT[[i]] = SAPG_OneRun[[i]][[7]]
            PHI_INIT[[i]] = SAPG_OneRun[[i]][[7]][[1]]#SAPG_OneRun[[i]][[8]] #MU_PI_MiseForme(MU_P[[i]],Design_Ind) #
            SD_STEP_INIT[[i]] = rep(1,NbParam) #SAPG_OneRun[[i]][[9]]
            G_SQUARED_INIT[[i]] = list(as.matrix(rep(0,dim(MU_INIT)[1])), matrix(0,NbParam,NbParam),0) #SAPG_OneRun[[i]][[10]]
            N_ITER_PREV[i] = SAPG_OneRun[[i]][[11]]
            
            DELTA_INIT_PART[[i]] = diag(diag(chol(var(SAPG_OneRun[[i]][[7]][[1]]))))
            
            if ((EBIC_R_SEQ[k,i] <= PersonalBest_VALUE[i]) & !is.na(EBIC_R_SEQ[k,i]) ){
                PersonalBest[i,] = LAMBDA[i,]
                PersonalBest_VALUE[i] = EBIC_R_SEQ[k,i]
            }
            
            if ( (EBIC_R_SEQ[k,i] <= GlobalBest_VALUE) & !is.na(EBIC_R_SEQ[k,i]) ){
                GlobalBest = LAMBDA[i,]
                GlobalBest_VALUE = EBIC_R_SEQ[k,i]
                GlobalBest_ID = i
            }
            
        }
        #browser()
        MU_P_SEQ[[1+k]] = MU_P
        GAMMA_P_SEQ[[1+k]] = GAMMA_P
        DELTA_P_SEQ[[1+k]] = DELTA_P
        ERROR_P_SEQ[[1+k]] = ERROR_P
        
        MU_P_R_SEQ[[1+k]] = MU_P_R
        GAMMA_P_R_SEQ[[1+k]] = GAMMA_P_R
        DELTA_P_R_SEQ[[1+k]] = DELTA_P_R
        ERROR_P_R_SEQ[[1+k]] = ERROR_P_R
        
        if (sum(EBIC_R_SEQ[k,]==GlobalBest_VALUE)!=0){
            MU_G_SEQ[[k]] = SAPG_OneRunR[[GlobalBest_ID]][[1]]
            GAMMA_G_SEQ[[k]] = SAPG_OneRunR[[GlobalBest_ID]][[22]]
            DELTA_G_SEQ[[k]] = SAPG_OneRunR[[GlobalBest_ID]][[23]]
            ERROR_G_SEQ[[k]] = SAPG_OneRunR[[GlobalBest_ID]][[3]]
        }else{
            MU_G_SEQ[[k]] = MU_G_SEQ[[k-1]]
            GAMMA_G_SEQ[[k]] = GAMMA_G_SEQ[[k-1]]
            DELTA_G_SEQ[[k]] = DELTA_G_SEQ[[k-1]]
            ERROR_G_SEQ[[k]] = ERROR_G_SEQ[[k-1]]
        }
        
        
        
        for (i in 1:N_PARTICLE){
            
            ### Mise à jour des LAMBDA ###
            set.seed(k*100*i)
            R1 = runif(2)
            R2 = runif(2)
            Velocity[i,] = INIERTIA*Velocity[i,] + C_AMP_1*(PersonalBest[i,] - LAMBDA[i,])*R1 + C_AMP_2*(GlobalBest - LAMBDA[i,])*R2
            Velocity[i,] = pmin(Velocity[i,],Vmax)
            Velocity[i,] = pmax(Velocity[i,],-Vmax)
            
            LAMBDA[i,] = LAMBDA[i,] + Velocity[i,]
            
            LAMBDA[i,1] = max(LAMBDA[i,1] , 0)
            if(LAMBDA[i,1] ==0 ){ LAMBDA[i,1] = LAMBDA[i,1] - runif(1)*Velocity[i,1]  }
            LAMBDA[i,1] = min(LAMBDA[i,1] , LAMBDA_MAX[1])
            
            LAMBDA[i,2] = max(LAMBDA[i,2] , 0)
            if(LAMBDA[i,2] ==0 ){ LAMBDA[i,2] = LAMBDA[i,2] - runif(1)*Velocity[i,2]  }
            LAMBDA[i,2] = min(LAMBDA[i,2] , LAMBDA_MAX[2])
            
        }
        
        LAMBDA_SEQ_PSO[[k + 1]] = LAMBDA
        
    }
    
    return(list(LAMBDA_SEQ_PSO,EBIC_SEQ,MU_P_SEQ,GAMMA_P_SEQ,DELTA_P_SEQ,ERROR_P_SEQ,EBIC_R_SEQ,MU_P_R_SEQ,GAMMA_P_R_SEQ,DELTA_P_R_SEQ,ERROR_P_R_SEQ,MU_G_SEQ,GAMMA_G_SEQ,DELTA_G_SEQ,ERROR_G_SEQ,GlobalBest))
    
}


