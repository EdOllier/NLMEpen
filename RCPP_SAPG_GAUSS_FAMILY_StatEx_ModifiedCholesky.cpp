#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat A1C1_MODEL(arma::mat PHI, arma::mat TIME) {
    
    int nsubj = PHI.n_rows, ntimes = TIME.n_cols;
    arma::mat out(nsubj, ntimes);
    
    double dose = 1000;
    /*double TINF = 8;*/
    
    for(int i = 0; i < nsubj; ++i) {
        
        double Vc = exp(PHI(i,0));
        double Vp = exp(PHI(i,1));
        double Q = exp(PHI(i,2));
        double Cl = exp(PHI(i,3));
        
        double Beta = ( (Q + Cl)/Vc + Q/Vp - sqrt( pow((Q + Cl)/Vc + Q/Vp,2.0) - 4*(Q*Cl/(Vc*Vp)) ) )/2 ;
        double Alpha = (Q*Cl)/(Vc*Vp*Beta) ;
        double Amc = ( ( Alpha -(Q/Vp) )/( Vc*( Alpha - Beta ) ) );
        double Bmc = ( ( Beta - (Q/Vp) )/( Vc*( Beta - Alpha ) ) );
        
        
        for(int j = 0; j < ntimes; ++j) {
            
            out(i,j) =  dose*( Amc*exp(-Alpha*TIME(i,j)) + Bmc*exp(-Beta*TIME(i,j))  );
            
        }
        
    }
    
    return out;
    
}


// [[Rcpp::export]]
arma::mat GRAD_LLcomp_MU(arma::mat S_1 ,Rcpp::List Design_Ind, arma::mat MU_PI, arma::mat invOMEGA ){
	
	arma::mat D0 = Design_Ind[0];
	int n_param = D0.n_cols;
	int n_subject = S_1.n_rows;
	
	arma::mat GRAD_MU = arma::zeros(n_param,1);
	
	for(int i=0; i<n_subject; i++){ 
		arma::mat DI = Design_Ind[i];
		GRAD_MU += (DI.t())*invOMEGA*(S_1.row(i).t() - MU_PI.row(i).t());
	}
	
	return(GRAD_MU);
    
}


// [[Rcpp::export]]
arma::mat GRAD_LLcomp_GAMMA(arma::mat Sigma, arma::mat DELTA, arma::mat GAMMA ,double n_subject){
    
    int n_param = GAMMA.n_cols;
    arma::mat invDELTA = inv(DELTA);
    
    arma::mat GRAD_GAMMA =   (inv(GAMMA*GAMMA.t()))*(invDELTA*Sigma*invDELTA + (invDELTA*Sigma*invDELTA).t())*inv(GAMMA*GAMMA.t())*GAMMA/2;
    
    GRAD_GAMMA = trimatl(GRAD_GAMMA,-1);
    
    return(GRAD_GAMMA);
}


// [[Rcpp::export]]
arma::mat GRAD_LLcomp_DELTA(arma::mat Sigma,arma::mat DELTA,arma::mat GAMMA,double n_subject){
    
    int n_param = DELTA.n_cols;
    arma::mat invDELTA = inv(DELTA);
    arma::mat invGAMMA = inv(GAMMA);
    arma::mat invSigma = inv(Sigma);
    
    arma::mat GRAD_DELTA =  -n_subject*invDELTA + ( inv(DELTA)*inv(GAMMA*GAMMA.t())*inv(DELTA*invSigma*DELTA) );
    
    GRAD_DELTA =  GRAD_DELTA%arma::eye(n_param,n_param);
    
    return(GRAD_DELTA);
}

// [[Rcpp::export]]
double GRAD_LLcomp_ERROR(double S_3,double ERR, double n_times, double n_subject){
	
	double GRAD_ERROR = -(n_times*n_subject)/ERR + ( S_3/pow(ERR,3.0) );
	
	return(GRAD_ERROR);
}



// [[Rcpp::export]]
Rcpp::List Stat_Ex_SA(Rcpp::List PHI_SIM, Rcpp::List PredModel_SIM, arma::mat Obs, arma::mat beta, Rcpp::List Design_Ind, double delta_SA, Rcpp::List Stat_Ex_Old){
    
    int n_sim = PHI_SIM.size();
    int n_subject = Obs.n_rows;
    arma::mat PHI0 = PHI_SIM[0];
    int n_param = PHI0.n_cols;
    
    arma::mat  S1_Old = Stat_Ex_Old[0];
    Rcpp::List S2_Old = Stat_Ex_Old[1];
    arma::mat  Sigma_Old = Stat_Ex_Old[2];
    arma::mat  Sigma_Squared_Old = Stat_Ex_Old[3];
    double S3_Old = Stat_Ex_Old[4];
    double S3_Squared_Old = Stat_Ex_Old[5];
    
    /* Initialisation */
    arma::mat S_1 = arma::zeros(n_subject,n_param);
    Rcpp::List S_2(n_subject);
    for(int j=0; j<n_subject; j++){
        S_2[j] =  arma::zeros(n_param,n_param);
    }
    arma::mat Sigma = arma::zeros(n_param,n_param);
    arma::mat Sigma_Squared = arma::zeros(n_param*n_param,n_param*n_param);
    double S_3 = 0;
    double S3_Squared = 0;
    
    for(int i=0; i<n_sim; i++){
        
        /* S1 */
        arma::mat PHItemp = PHI_SIM[i];
        S_1 += PHItemp;
        
        /* S2, Sigma et Sigma squared */
        arma::mat Sigmatemp = arma::zeros(n_param,n_param);
        
        for(int j=0; j<n_subject; j++){
            arma::mat S2temp = ((PHItemp.row(j)).t())*(PHItemp.row(j));
            arma::mat Design_Ind_j = Design_Ind[j];
            Sigmatemp += ( S2temp - (Design_Ind_j*beta)*(PHItemp.row(j)) -  ((PHItemp.row(j)).t())*((Design_Ind_j*beta).t()) + (Design_Ind_j*beta)*((Design_Ind_j*beta).t()) );
            arma::mat S2J = S_2[j];
            S_2[j] =  S2J + S2temp;
        }
        
        Sigma += Sigmatemp;
        Sigma_Squared += kron(Sigmatemp,Sigmatemp);
        
        /* S3 et S3 squared */
        arma::mat pred_mat = PredModel_SIM[i];
        double S3temp = arma::accu(arma::pow(Obs-pred_mat,2.0));
        S_3 += S3temp;
        S3_Squared += S3temp*S3temp;
        
    }
    
    for(int j=0; j<n_subject; j++){
        arma::mat S2Imat = S2_Old[j];
        arma::mat S2J = S_2[j];
        S_2[j] = S2Imat + delta_SA*(S2J/n_sim - S2Imat);
    }
    
    S_1 = S1_Old + delta_SA*(S_1/n_sim - S1_Old);
    Sigma = Sigma_Old + delta_SA*(Sigma/n_sim - Sigma_Old);
    Sigma_Squared = Sigma_Squared_Old + delta_SA*(Sigma_Squared/n_sim - Sigma_Squared_Old);
    S_3 = S3_Old + delta_SA*(S_3/n_sim - S3_Old);
    S3_Squared = S3_Squared_Old + delta_SA*(S3_Squared/n_sim - S3_Squared_Old);
    
    Rcpp::List out(6);
    out[0] = S_1;
    out[1] = S_2;
    out[2] = Sigma;
    out[3] = Sigma_Squared;
    out[4] = S_3;
    out[5] = S3_Squared;
    
    return(out);
    
}


// [[Rcpp::export]]
arma::mat LLcompIS(arma::mat Obs,arma::mat PredModel,arma::mat B, arma::mat OMEGA, arma::mat invOMEGA, double ERR){
	
	int n_obs = Obs.n_cols;
	int n_param = B.n_cols;
	int n_subject = Obs.n_rows;

	
	arma::mat out = arma::zeros(n_subject,1);
	
	double Cst = -n_obs*log(ERR) - 0.5*(n_obs + n_param)*log(2*arma::datum::pi) - 0.5*log(arma::det(OMEGA));
	
	for(int i=0; i<n_subject; i++){
		
		out(i,0) = Cst - 0.5*sum(arma::pow((Obs.row(i)-PredModel.row(i))/ERR,2.0)) - 0.5*sum(B.row(i) % (B.row(i) *  invOMEGA) );
		
	}
	
	return(out);
}

// [[Rcpp::export]]
double LIK_IS(arma::mat OBS_init,arma::mat OBS_TIME_init,arma::mat OMEGA, arma::mat invOMEGA, double ERR, arma::mat MU_PI_init, arma::mat MeanPhiPost_init , arma::mat VarPhiPost_init){ /*,Rcpp::List PHI_SIM_MCMC, int NBURN*/
    
    double nuis = 4;
    int MM = 100;
    int NbrIterIS = 10000;
    int NbrSubject = OBS_init.n_rows;
    int NbrParam = OMEGA.n_rows;
    int NITER = round(NbrIterIS/MM);
    /*int NSIM = PHI_SIM_MCMC.size();*/
    arma::vec LLIK = arma::zeros(NITER);
    arma::mat llik = arma::zeros(NbrSubject,1);

    arma::mat SdPhiPost_init = arma::sqrt(VarPhiPost_init);
    
    arma::mat MeanPhiPost = repmat(MeanPhiPost_init,MM,1);
    arma::mat VarPhiPost = repmat(VarPhiPost_init,MM,1);
    arma::mat SdPhiPost = repmat(SdPhiPost_init,MM,1);
    arma::mat OBS = repmat(OBS_init,MM,1);
    arma::mat OBS_TIME = repmat(OBS_TIME_init,MM,1);
    arma::mat MU_PI = repmat(MU_PI_init,MM,1);
    
    arma::mat rnd = arma::zeros((NbrSubject*MM),NbrParam);
    arma::mat rnd_pdf = arma::zeros((NbrSubject*MM),NbrParam);
    arma::mat llik_new = arma::zeros((NbrSubject*MM),1);
    arma::mat llik_new_vect = arma::zeros(NbrSubject,1);
    arma::mat rnd_prob = arma::zeros((NbrSubject*MM),1);
    arma::mat var_phi_cst = arma::sum(arma::log(VarPhiPost),1)/2;
    
    for (int iter = 0; iter<NITER; iter++){
        
        for (int i = 0; i<(NbrSubject*MM); i++){
            for (int j = 0; j<NbrParam; j++){
                rnd(i,j) =  R::rt(nuis);
                rnd_pdf(i,j) =  R::dt(rnd(i,j),nuis,true);
            }
        }
        
        rnd_prob = arma::sum(rnd_pdf,1) - var_phi_cst;
   
        arma::mat PHI = MeanPhiPost + SdPhiPost%rnd;
        arma::mat B = PHI - MU_PI;
        
        arma::mat Pred = A1C1_MODEL(PHI,OBS_TIME);
        
        llik_new = exp(LLcompIS(OBS,Pred,B,OMEGA,invOMEGA,ERR)  - rnd_prob);
        
        llik_new.reshape(NbrSubject,MM);
        llik_new = arma::max(llik_new,2.225074e-308*arma::ones(NbrSubject,MM));
        /*llik_new = arma::mean(llik_new, 1);*/
        
        for (int i = 0; i<NbrSubject; i++){
            arma::mat llik_new_ind = llik_new.row(i);
            double M_ll_ind = arma::mean( llik_new_ind.elem( find_finite(llik_new_ind) ) );
            llik_new_vect(i,0) = M_ll_ind;
        }
        
        
        llik = llik + (llik_new_vect - llik)/(iter+1);
        
        LLIK[iter] = arma::accu(arma::log(llik));
    
    }

    return(LLIK[NITER-1]);
}



// [[Rcpp::export]]
arma::mat SOFT_TRESH( arma::mat X, arma::mat lambda){
    
    return( arma::sign(X) % arma::max(arma::abs(X) - lambda, arma::zeros( size(X) ) )  );
    
}


// [[Rcpp::export]]
arma::mat DEFPOS_PROJ( arma::mat X){
    
    double tresh = 0.0000001;
    
    arma::mat eigvec;
    arma::vec eigval;
    arma::eig_sym(eigval,eigvec,X);
    
    /*eigvec = arma::sign(eigvec) % arma::max(arma::abs(eigvec) - 0.00000000001, arma::zeros( size(X) ) ) ;*/
    arma::vec eigval_tresh = arma::max(eigval, tresh*arma::ones( size(eigval) ) ) ;
    
    arma::mat X_tresh =  eigvec*( arma::diagmat(eigval_tresh) )*(eigvec.t());
    X_tresh = arma::sign(X_tresh) % arma::max(arma::abs(X_tresh) - 0.000000001, arma::zeros( size(X) ) );
    
    return( X_tresh  );
    
}

// [[Rcpp::export]]
arma::mat DIAG_POS_PROJ( arma::mat X){
    
    double tresh = 0.00001;
    
    arma::vec diag_tresh = arma::max(diagvec(X), tresh*arma::ones( size(diagvec(X)) ) ) ;
    
    arma::mat X_tresh =   arma::diagmat(diag_tresh);
    
    return( X_tresh  );
    
}

// [[Rcpp::export]]
Rcpp::List MCMC_SIM(arma::mat PHI_INIT,int N_ITER_MCMC,int NbrIterN2,arma::vec SD_STEP_INIT, arma::mat OMEGA,arma::mat invOMEGA,arma::mat MU_PI,double ERR,arma::mat OBS, arma::mat OBS_TIME) {


	double AcRateTheo = 0.3;
	double RW = 0.4;

	int NbParam = MU_PI.n_cols;
	int NbSujet = OBS.n_rows;

	arma::mat chol_omega = arma::chol(OMEGA);
	
	arma::mat B = PHI_INIT - MU_PI;
	arma::mat PHI = PHI_INIT;

	arma::mat PredModel = A1C1_MODEL(PHI,OBS_TIME);
	/* Rcpp::List Grad_PredModel = A1C1_MODEL_GRAD(PHI,OBS_TIME); */
	arma::vec Prob = sum( 0.5 * pow((OBS - PredModel)/ERR,2.0),1);

	Rcpp::List B_SIM((NbrIterN2+1)*N_ITER_MCMC);
    Rcpp::List PHI_SIM((NbrIterN2+1)*N_ITER_MCMC);
	Rcpp::List PredModel_SIM((NbrIterN2+1)*N_ITER_MCMC);
	/* Rcpp::List Grad_PredModel_SIM((NbrIterN2+1)*N_ITER_MCMC); */
	
	arma::vec NbrAccept = arma::zeros( NbParam);
	arma::vec SD_STEP = SD_STEP_INIT;

	int cursor = 0;

	for(int iter=0; iter<N_ITER_MCMC; iter++){
		
		arma::mat B_NEW = arma::randn(NbSujet,NbParam) * chol_omega;
		PHI = MU_PI+ B_NEW;
		arma::mat PredModel_NEW = A1C1_MODEL(PHI,OBS_TIME);
		
		arma::vec Prob_NEW = sum( 0.5 * pow((OBS - PredModel_NEW)/ERR,2.0),1);
        arma::vec DELTA = Prob_NEW  -  Prob + log(arma::randu(NbSujet));
        
		for (int i=0; i<NbSujet; i ++){
			
			if (DELTA[i]<0){
				B.row(i) = B_NEW.row(i);
				Prob[i] = Prob_NEW[i];
				PredModel.row(i) = PredModel_NEW.row(i);
			}
		
			
		}
		
		PHI = MU_PI+ B;
        /* Grad_PredModel = A1C1_MODEL_GRAD(PHI,OBS_TIME); */
		
		B_SIM[cursor] =  B;
        PHI_SIM[cursor] =  PHI;
        PredModel_SIM[cursor] =  PredModel;
        /* Grad_PredModel_SIM[cursor] = Grad_PredModel; */
		
		cursor += 1;
		
		NbrAccept = arma::zeros( NbParam);
		arma::vec Prob_B =  0.5*sum(B %( B * invOMEGA), 1);
		
		for(int iter_N2=0; iter_N2<NbrIterN2; iter_N2++){
			
			for(int id_par=0; id_par<NbParam; id_par++){
			
				B_NEW = B;
				B_NEW.col(id_par) = B.col(id_par) + arma::randn(NbSujet) * SD_STEP[id_par];
				
				PHI = MU_PI + B_NEW;
			 	PredModel_NEW = A1C1_MODEL(PHI,OBS_TIME);
			 	
				Prob_NEW = sum( 0.5 * pow((OBS - PredModel_NEW)/ERR,2.0),1);
    			arma::vec Prob_B_NEW = 0.5*sum(B_NEW %( B_NEW * invOMEGA), 1);
    			
    			DELTA = Prob_NEW  -  Prob + Prob_B_NEW - Prob_B + log(arma::randu(NbSujet));
			
				for (int i=0; i<NbSujet; i ++){
			
					if ( DELTA[i] < 0 ){
						B.row(i) = B_NEW.row(i);
						Prob[i] = Prob_NEW[i];
						PredModel.row(i) = PredModel_NEW.row(i);
						Prob_B[i]= Prob_B_NEW[i];
						NbrAccept[id_par] += 1; 
					}
		
			
				}
			
				
			}
			
			PHI = MU_PI+ B;
        	/* Grad_PredModel = A1C1_MODEL_GRAD(PHI,OBS_TIME); */
        	
        	B_SIM[cursor] =  B;
            PHI_SIM[cursor] =  PHI;
        	PredModel_SIM[cursor] =  PredModel;
        	/* Grad_PredModel_SIM[cursor] = Grad_PredModel; */
		
			cursor += 1;
		
		
		}
		
		SD_STEP += RW*(SD_STEP%(NbrAccept/(NbSujet*NbrIterN2) - AcRateTheo) );
		
	}

	Rcpp::List out(5);

	out[0] = PHI_SIM;
	out[1] = PredModel_SIM;
	out[2] = SD_STEP;
    out[3] = B_SIM;
	out[4] = NbrAccept;

	return(out);

}

// [[Rcpp::export]]
arma::mat MU_PI_MiseForme(arma::mat MU, Rcpp::List DesignInd){
    
    arma::mat DD = DesignInd[0];
    int NP = DD.n_rows;
    int NS = DesignInd.size();
    
    arma::mat MU_PI = arma::zeros(NS,NP);
    
    for (int i = 0; i<NS; i++){
        arma::mat DI = DesignInd[i];
        MU_PI.row(i) = (DI*MU).t();
    }
    
    return(MU_PI);
    
}


// [[Rcpp::export]]
Rcpp::List SAPG_ADAGRAD_OneRun(arma::mat OBS,arma::mat OBS_TIME, Rcpp::List Design_Ind, arma::mat MU_INIT, arma::mat GAMMA_INIT, arma::mat DELTA_INIT,double ERROR_INIT, Rcpp::List LAMBDA, int N_ITER, Rcpp::List AdaptWEIGHTS, int NbrIterN2, int N_ITER_MCMC, double beta, double delta_0, int Nbr_Iter_Prox, double step_adagrad, Rcpp::List STAT_EX_INIT, arma::mat PHI_INIT, arma::vec SD_STEP_INIT, Rcpp::List G_SQUARED_INIT, int N_ITER_PREV ){
    
    
    int NbrSubject = OBS.n_rows;
    int NbrTimes = OBS.n_cols;
    int NbrParam = GAMMA_INIT.n_rows;
    int NbrParamCov = MU_INIT.n_rows;
    
    double L1 = LAMBDA[0];
    double L2 = LAMBDA[1];
    
    int Nbr_ITER_MAX = N_ITER;
    
    /*### SPARSITY PARAMETER ###*/
    
    arma::mat AdaptWEIGHTS_MU = AdaptWEIGHTS[0];
    arma::mat AdaptWEIGHTS_GAMMA = AdaptWEIGHTS[1];
    
    /*### PARAM GRADIENT DESCENT ###*/
    double eps_cst = 0.00000001;
    
    /*### SOCHASTIC APPROXIMATION ###*/
    double C_SA = delta_0;
    
    arma::vec C_SA2 = arma::zeros(Nbr_ITER_MAX);
    
    for (int i = 0; i<Nbr_ITER_MAX; i++){
        /*if ( i < N_ITER_PREV ){
            C_SA2[i] = C_SA;
        }else{*/
            C_SA2[i] = C_SA/( pow( i + N_ITER_PREV + 1, beta ) );
       /* }*/
    }
    
    /*### PARAM MCMC ###*/
    
    int NbrSIM = N_ITER_MCMC*(NbrIterN2 + 1);
    arma::vec SD_STEP = SD_STEP_INIT;
    
    arma::mat MU = MU_INIT;
    arma::mat MU_PI = MU_PI_MiseForme(MU,Design_Ind);
    arma::mat GAMMA = GAMMA_INIT;
    arma::mat DELTA = DELTA_INIT;
    arma::mat OMEGA = DELTA*GAMMA*GAMMA.t()*DELTA;
    arma::mat invOMEGA = inv(OMEGA);
    double ERRORsd = ERROR_INIT;
    
    /*##########################*/
    /*### BOUCLE PRINCIPALE ###*/
    /*#########################*/
    
    arma::vec IC_Q_SEQ = arma::zeros(Nbr_ITER_MAX*Nbr_Iter_Prox);
    arma::vec BIC_IC_Q_SEQ = arma::zeros(Nbr_ITER_MAX*Nbr_Iter_Prox);
    arma::vec EBIC_IC_Q_SEQ = arma::zeros(Nbr_ITER_MAX*Nbr_Iter_Prox);
    
    Rcpp::List MU_SEQ(Nbr_ITER_MAX+1);
    MU_SEQ[0] = MU;
    Rcpp::List GAMMA_SEQ(Nbr_ITER_MAX+1);
    GAMMA_SEQ[0] = GAMMA;
    Rcpp::List DELTA_SEQ(Nbr_ITER_MAX+1);
    DELTA_SEQ[0] = DELTA;
    Rcpp::List OMEGA_SEQ(Nbr_ITER_MAX+1);
    OMEGA_SEQ[0] = OMEGA;
    Rcpp::List ERROR_SEQ(Nbr_ITER_MAX+1);
    ERROR_SEQ[0] = ERRORsd;
    
    Rcpp::List STEP_MU_SEQ(Nbr_ITER_MAX+1);
    Rcpp::List STEP_GAMMA_SEQ(Nbr_ITER_MAX+1);
    Rcpp::List STEP_DELTA_SEQ(Nbr_ITER_MAX+1);
    Rcpp::List STEP_ERROR_SEQ(Nbr_ITER_MAX+1);

    
    /*Rcpp::List Stat_Ex_SEQ(Nbr_ITER_MAX);*/
    
    Rcpp::List Stat_Ex(6);
    Stat_Ex[0] = STAT_EX_INIT[0];
    Stat_Ex[1] = STAT_EX_INIT[1];
    Stat_Ex[2] = STAT_EX_INIT[2];
    Stat_Ex[3] = STAT_EX_INIT[3];
    Stat_Ex[4] = STAT_EX_INIT[4];
    Stat_Ex[5] = STAT_EX_INIT[5];
    
    Rcpp::List Stat_Ex_old = Stat_Ex;
    
    arma::mat GRAD_MU = arma::zeros(NbrParamCov,1);
    arma::mat GRAD_GAMMA = arma::zeros(NbrParam,NbrParam);
    arma::mat GRAD_DELTA = arma::zeros(NbrParam,NbrParam);
    double GRAD_ERROR = 0;
    
    arma::mat HESS_MU = arma::zeros(NbrParamCov,1);
    arma::mat HESS_GAMMA = arma::zeros(NbrParam,NbrParam);
    arma::mat HESS_DELTA = arma::zeros(NbrParam,NbrParam);
    double HESS_ERROR = 0;
    
    arma::vec SUP_SIZE_MU_SEQ = arma::zeros(Nbr_ITER_MAX);
    arma::vec SUP_SIZE_OMEGA_SEQ = arma::zeros(Nbr_ITER_MAX);
    
    int ddd = 0;
    
    for (int iter = 0; iter<Nbr_ITER_MAX; iter ++){
        
        /*###########################*/
        /*### SIMULATION PAR MCMC ###*/
        /*###########################*/
        
        Rcpp::List SIM = MCMC_SIM(PHI_INIT,N_ITER_MCMC, NbrIterN2,SD_STEP,OMEGA,invOMEGA,MU_PI,ERRORsd,OBS,OBS_TIME);
        
        arma::mat SD_STEP_temp = SIM[2];
        SD_STEP = SD_STEP_temp;
        
        Rcpp::List PHI_SIM = SIM[0];
        arma::mat PHI_INIT_temp = PHI_SIM[NbrSIM-1];
        PHI_INIT = PHI_INIT_temp;
        
        /*##############################################*/
        /*### APPROXIMATION STATISTIQUES EXHAUSTIVES ###*/
        /*##############################################*/
        
        Stat_Ex_old = Stat_Ex;
        Stat_Ex = Stat_Ex_SA(SIM[0], SIM[1], OBS, MU, Design_Ind, C_SA2[iter],  Stat_Ex_old);
        /*Stat_Ex_SEQ[iter] = Stat_Ex;*/
        
        for (int iter_PROX = 0; iter_PROX<Nbr_Iter_Prox; iter_PROX ++){
            
            /*############################*/
            /*### GRADIENT COMPUTATION ###*/
            /*############################*/
            
            GRAD_MU = GRAD_LLcomp_MU(Stat_Ex[0] , Design_Ind, MU_PI, invOMEGA);
            GRAD_GAMMA = GRAD_LLcomp_GAMMA(Stat_Ex[2], DELTA, GAMMA, NbrSubject);
            GRAD_DELTA = GRAD_LLcomp_DELTA(Stat_Ex[2], DELTA, GAMMA, NbrSubject);
            GRAD_ERROR = GRAD_LLcomp_ERROR(Stat_Ex[4], ERRORsd, NbrTimes, NbrSubject);
            
            /*###################################*/
            /*### GRADIENT CUM L2 COMPUTATION ###*/
            /*###################################*/
            
            HESS_MU += arma::pow(GRAD_MU,2.0);
            HESS_GAMMA +=  arma::pow(GRAD_GAMMA,2.0);
            HESS_DELTA +=  arma::pow(GRAD_DELTA,2.0);
            HESS_ERROR +=  pow(GRAD_ERROR,2.0);
            
            /*##############################*/
            /*### Proximal Gradient Step ###*/
            /*##############################*/
            
            /*### MU ###*/
            arma::mat STEP_MU = step_adagrad/(arma::sqrt(HESS_MU + eps_cst ));
            arma::mat MU_NEW = SOFT_TRESH( MU + (STEP_MU % GRAD_MU),L1*(STEP_MU % AdaptWEIGHTS_MU) );
            arma::mat MU_PI_NEW = MU_PI_MiseForme(MU_NEW,Design_Ind);
            
            /*### GAMMA ###*/
            arma::mat STEP_GAMMA = step_adagrad/(arma::sqrt(HESS_GAMMA + eps_cst ));
            arma::mat GAMMA_NEW = SOFT_TRESH(GAMMA + (STEP_GAMMA % GRAD_GAMMA),L2*(STEP_GAMMA % AdaptWEIGHTS_GAMMA));
            
            /*### DELTA ###*/
            arma::mat STEP_DELTA = step_adagrad/(arma::sqrt(HESS_DELTA + eps_cst ));
            arma::mat DELTA_NEW = DELTA + (STEP_DELTA % GRAD_DELTA);
            DELTA_NEW = DIAG_POS_PROJ(DELTA_NEW);
            
            /*### OMEGA ###*/
            arma::mat OMEGA_NEW = DELTA_NEW*GAMMA_NEW*GAMMA_NEW.t()*DELTA_NEW;
            arma::mat invOMEGA_NEW = inv(OMEGA_NEW);
            
            /*### ERROR ###*/
            
            double STEP_ERROR = step_adagrad/(sqrt(HESS_ERROR + eps_cst ));
            double ERROR_NEW = ERRORsd + STEP_ERROR*GRAD_ERROR;
            if (ERROR_NEW<=0.00005){ERROR_NEW=0.00005;}
            
            
            /*### UPDATE CURRENT ITERATES ###*/
            
            MU = MU_NEW;
            MU_PI = MU_PI_NEW;
            GAMMA = GAMMA_NEW;
            DELTA = DELTA_NEW;
            OMEGA = OMEGA_NEW;
            invOMEGA= invOMEGA_NEW;
            ERRORsd = ERROR_NEW;
            
            /*if(iter_PROX==0){
                STEP_MU_SEQ[iter+1] = STEP_MU;
                STEP_GAMMA_SEQ[iter+1] =  STEP_GAMMA;
                STEP_DELTA_SEQ[iter+1] =  STEP_DELTA;
                STEP_ERROR_SEQ[iter+1] = STEP_ERROR;
            }*/
            
            /*IC_Q_SEQ[ddd] = -2*IC_Q(Stat_Ex[2], Stat_Ex[4], OMEGA, invOMEGA, ERRORsd, NbrTimes, NbrSubject);
            BIC_IC_Q_SEQ[ddd] = IC_Q_SEQ[ddd] + L1*accu(abs(AdaptWEIGHTS_MU%MU)) + L2*accu(abs(AdaptWEIGHTS_GAMMA%GAMMA));
            ddd += 1;*/
        }
        
        /*MU_SEQ[iter+1] = MU;
        GAMMA_SEQ[iter+1] =  GAMMA;
        DELTA_SEQ[iter+1] =  DELTA;
        OMEGA_SEQ[iter+1] =  OMEGA;
        ERROR_SEQ[iter+1] = ERRORsd;*/
        
        /*MU_SEQ[iter+1] = MU;
        GAMMA_SEQ[iter+1] =  GAMMA;
        DELTA_SEQ[iter+1] =  DELTA;
        OMEGA_SEQ[iter+1] =  OMEGA;
        ERROR_SEQ[iter+1] = ERRORsd;*/
        
        
        
        
    }
    
    /*MU = median(MU_SEQ.cols(Nbr_ITER_MAX-100,Nbr_ITER_MAX),1);
    MU_PI = MU_PI_MiseForme(MU,Design_Ind);
    arma::mat OMEGAtt = (arma::median(OMEGA_SEQ.cols(Nbr_ITER_MAX-100, Nbr_ITER_MAX),1));
    OMEGAtt.reshape(NbrParam,NbrParam);
    OMEGA = OMEGAtt;
    invOMEGA = inv(OMEGA);*/
    
    arma::uvec NZ_MU = find(MU);
    int SUP_SIZE_MU = NZ_MU.size();
    
    arma::uvec NZ_GAMMA = find(GAMMA);
    int SUP_SIZE_GAMMA = (NZ_GAMMA.size() - NbrParam);
    
    /*arma::uvec NZ_OMEGA = find(OMEGA);
    int SUP_SIZE_OMEGA = (NZ_OMEGA.size() - NbrParam)/2;*/
    
    SD_STEP = arma::ones(NbrParam);
    /*Rcpp::List SIM = MCMC_SIM(PHI_INIT,30, 100,SD_STEP,OMEGA,invOMEGA,MU_PI,ERRORsd,OBS,OBS_TIME);
     Rcpp::List SIM_PHI = SIM[0];*/
    arma::mat COND_MEAN_PHI = Stat_Ex[0];
    arma::mat COND_VAR_PHI = arma::zeros(NbrSubject,NbrParam);
    Rcpp::List S2_LIST = Stat_Ex[1];
    for(int j=0; j<NbrSubject; j++){
        arma::mat S2_IND_MAT = S2_LIST[j];
        arma::vec S2_IND_MAT_DIAG = S2_IND_MAT.diag();
        COND_VAR_PHI.row(j) = S2_IND_MAT_DIAG.t();
    }
    COND_VAR_PHI = COND_VAR_PHI - arma::pow(COND_MEAN_PHI,2.0);
    double LIK = -2*LIK_IS(OBS, OBS_TIME, OMEGA, invOMEGA, ERRORsd, MU_PI,COND_MEAN_PHI , COND_VAR_PHI);/*SIM_PHI,1000);*/
    int NparamTotSel = NbrParamCov + NbrParam*(NbrParam-1)/2;
    int NparamSel = SUP_SIZE_MU + SUP_SIZE_GAMMA;
    double BIC = LIK + ( NparamSel + NbrParam )*log(NbrSubject);
    double EBIC = LIK + (NparamSel + NbrParam )*log(NbrSubject);
    
    Rcpp::List HESS_LIST(4);
    HESS_LIST[0] = HESS_MU;
    HESS_LIST[1] = HESS_GAMMA;
    HESS_LIST[2] = HESS_DELTA;
    HESS_LIST[3] = HESS_ERROR;
    
    Rcpp::List out(27);
    
    out[0] = MU;
    out[1] = OMEGA;
    out[2] = ERRORsd;
    out[3] = IC_Q_SEQ;
    out[4] = BIC_IC_Q_SEQ;
    out[5] = EBIC_IC_Q_SEQ;
    out[6] = Stat_Ex;
    out[7] = PHI_INIT;
    out[8] = SD_STEP;
    out[9] = HESS_LIST;
    out[10] = N_ITER_PREV+Nbr_ITER_MAX;
    out[11] = SUP_SIZE_MU_SEQ;
    out[12] = SUP_SIZE_OMEGA_SEQ;
    out[13] = LIK;
    out[14] = BIC;
    out[15] = EBIC;
    out[16] = MU_SEQ;
    out[17] = OMEGA_SEQ;
    out[18] = ERROR_SEQ;
    out[19] = GAMMA_SEQ;
    out[20] = DELTA_SEQ;
    out[21] = GAMMA;
    out[22] = DELTA;
    out[23] = STEP_MU_SEQ;
    out[24] = STEP_GAMMA_SEQ;
    out[25] = STEP_DELTA_SEQ;
    out[26] = STEP_ERROR_SEQ;
    
    
    
    return(out);
    
}





