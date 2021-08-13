

source('CodeR_SPAG_WR_PSO.R')

library(lhs)
library(MASS)

###########
### TEST ###
###########
set.seed(1234)

library(MASS)
set.seed(500)
NbrSubject = 60
NbParam = 4
NbrCOV = 50

Time = matrix(rep(c(0.1,0.33,0.75,1,2,4,8),times=NbrSubject),NbrSubject,7,byrow=T)


mu_Vc = log(6.17)
mu_Vp = log(9.56)
mu_Cl = log(5.31)
mu_Q = log(22.3)

MU_TH =  c(mu_Vc,mu_Vp,mu_Q,mu_Cl)

CARACT_MAT = rbind(CARACT_MAT,data.frame(rep=rep_pso,NbrCOV=NbrCOV))

set.seed((rep_pso)*100)
PI_TH = rbind(diag(c(0,0,0,0)),matrix(0,NbrCOV-4,NbParam))
PI_TH[2,1] = c(0.4)
#PI_TH[2,2] = c(0.4)
PI_TH[4,4] = c(0.4)

DESIGN_COV = matrix(rnorm(NbrSubject*NbrCOV),nrow=NbrSubject ,ncol=NbrCOV)

Design_Ind = list()
for (i in 1:NbrSubject){
    MatD =  cbind(diag(rep(1,NbParam)), matrix(0,NbParam,NbParam*NbrCOV))
    for (j in 1:NbParam){
        MatD[j,(NbParam+1+(j-1)*NbrCOV):(NbParam + j*NbrCOV )] = DESIGN_COV[i,]
    }
    Design_Ind[[i]] = MatD
}

MU_PI_TH = matrix(rep(MU_TH,times=NbrSubject),byrow=T,ncol=4) + DESIGN_COV%*%PI_TH

sigma_Vc = 0.4
sigma_Vp = 0.3
sigma_Q = 0.3
sigma_Cl = 0.4

DELTA_TH = diag(c(sigma_Vc,sigma_Vp,sigma_Q,sigma_Cl))
GAMMA_TH = diag(rep(1,NbParam))
GAMMA_TH[4,1]  = 1.5

OMEGA_TH = DELTA_TH%*%GAMMA_TH%*%t(GAMMA_TH)%*%DELTA_TH
invOMEGA_TH = ginv(OMEGA_TH)

chol.OMEGA_TH <- t(GAMMA_TH)%*%DELTA_TH

B_TH = matrix(rnorm(NbrSubject * NbParam), ncol = NbParam) %*% chol.OMEGA_TH
PHI_TH = MU_PI_TH + B_TH

PredModel_TH <- A1C1_MODEL(PHI_TH,Time)

ERROR= 5
Obs = PredModel_TH + ERROR*matrix(rnorm(prod(dim(PredModel_TH))), dim(PredModel_TH)[1] , dim(PredModel_TH)[2] )

################
### RUN  PSO ###
################

MU_INIT =  as.matrix(c(log(5),log(8),log(10),log(1),rep(0.2,NbrCOV*NbParam)))
MU_PI_INIT = MU_PI_MiseForme(MU_INIT,Design_Ind)
DELTA_INIT = diag(rep(1,4))
GAMMA_INIT = diag(rep(1,4))
GAMMA_INIT[which(upper.tri(GAMMA_INIT)==TRUE)]=0
diag(GAMMA_INIT) = 1
invOMEGA_INIT = diag(rep(1,4))
ERROR_INIT = 25


delta_CSA = 0.9
NiterPSO = 10
N_P = 25

LAMBDA_MAX = c(80,40)
LAMBDA_INIT = (as.matrix(expand.grid(seq(0.1,(LAMBDA_MAX[1]),length.out=round(sqrt(N_P))) , seq(0.1,(LAMBDA_MAX[2]),length.out=round(sqrt(N_P))) ))) #exp(as.matrix(expand.grid(seq(-0.1,log(LAMBDA_MAX[1]),length.out=round(sqrt(N_P))) , seq(-0.1,log(LAMBDA_MAX[2]),length.out=round(sqrt(N_P))) )))

AdaptWEIGHTS = list()
AdaptWEIGHTS[[1]] = as.matrix(c(rep(0,NbParam),rep(1,NbrCOV*NbParam)))
AdaptWEIGHTS[[2]] = matrix(1,NbParam,NbParam)
AdaptWEIGHTS[[2]][which(upper.tri(AdaptWEIGHTS[[2]])==TRUE)]=0
diag(AdaptWEIGHTS[[2]])=0

############
### SAPG ###
############

set.seed(1234)
LAMBDA = c(20,6)
N_ITER = 10000
NbrIterN2 = 2
N_ITER_MCMC = 2
beta = 0.75
delta_0 = 0.9
Nbr_Iter_Prox = 20
STAT_EX_INIT = list(0*MU_PI_INIT, lapply(1:NbrSubject, function(x) matrix(0, nrow=NbParam, ncol=NbParam)) ,matrix(0,NbParam,NbParam),matrix(0,NbParam*NbParam,NbParam*NbParam),0,0)
PHI_INIT = MU_PI_INIT
SD_STEP_INIT = rep(1,4)
G_SQUARED_INIT = list(as.matrix(rep(0,dim(MU_INIT)[1])), matrix(0,NbParam,NbParam),0)
N_ITER_PREV = 0
step_adagrad = 0.2

test_PG_SA_OneRun <- SAPG_ADAGRAD_OneRun(Obs,Time, Design_Ind, MU_INIT,GAMMA_INIT,DELTA_INIT,ERROR_INIT, LAMBDA, N_ITER, AdaptWEIGHTS, NbrIterN2, N_ITER_MCMC, beta, delta_0, Nbr_Iter_Prox, step_adagrad, STAT_EX_INIT, PHI_INIT, SD_STEP_INIT, G_SQUARED_INIT, N_ITER_PREV )

# Estimated Beta
test_PG_SA_OneRun[[1]]

# Estimated Gamma
test_PG_SA_OneRun[[22]]

###########
### PSO ###
###########

set.seed(1234)
test_PSO <- PSO_PARALLEL_SAPG_ADAGRAD_R(Obs,Time,Design_Ind,MU_INIT,DELTA_INIT,GAMMA_INIT,ERROR_INIT,LAMBDA_INIT,LAMBDA_MAX,N_ITER_PerSEQ=200,AdaptWEIGHTS,NiterPSO,CONTROL=list(N2_ITER=c(2),N_MCMC = c(2) , beta = 0.75, delta_0 = delta_CSA, Nbr_Iter_Prox = 20, step_adagrad=0.2))

# Selected beta
BETA_SELECT = test_PSO[[12]][[NiterPSO]][-c(1:4)]
BETA_TRUE = as.numeric(PI_TH)
cbind(BETA_TRUE,BETA_SELECT)

# Selected Gamma
GAMMA_SELECT = test_PSO[[13]][[NiterPSO]]
GAMMA_TRUE = GAMMA_TH

