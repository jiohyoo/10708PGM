Description 
--  This package provides the Matlab implementation of learning and inference methods
--  for sparse conditional Gaussian graphical models (sCGGMs) as described in the following paper:
-- 
--  L. Zhang & S. Kim (2014). "Learning Gene Networks under SNP Perturbations Using eQTL Datasets", 
--  PLoS Computational Biology. 

Software URL
--  http://www.cs.cmu.edu/~sssykim/softwares/softwares.html#scggm


////////////////////////////////////////////////////
Matlab routines included in this package
////////////////////////////////////////////////////

------------------------------------------------------------
scggm.m		fit an sCGGM given regularization parameters
------------------------------------------------------------
Description
--  The main function that estimates sCGGM parameters given 
--  regularization parameters lambda_1 and lambda_2. 

Usage
--  "OPT = scggm( x, y, lambda_1, lambda_2, option)"

Details
--  sCGGM parameters, Theta.xy and Theta.yy, are estimated using training data. 
--  If option is not supplied, all fields will be set to default values. User 
--  can also specify only some of the fields in the option variable. The list 
--  of fields is below in "Arguments -- option". 


Arguments
--  x - training data in N by J matrix for genotypes, where N is the number of training samples and J is 
    		the number of SNPs.
--  y - training data in N by K matrix for gene expressions, where K is the number of gene-expression
    		traits.
    (x and y are not necessarily centered)
--  lambda_1/lambda_2 - regularization parameters. lambda_1 controls the sparsity level of Theta.xy and lambda_2 for Theta.yy. 
		Both are positive real numbers. 
--  option (optional) - a structure that specifies options to run the algorithm.
		It contains the following fields:
--  option.maxiter - the maximum number of iterations in accelerated proximal gradient descent. Default is 1000. 
--  option.ifrefit - a flag to specify whether to perform a refitting procedure or not.
		If option.ifrefit=false, the routine will return the parameter estimates from optimizaing L1-regularized
		data loglikelihood. If option.ifrefit=true, the routine will perform an additional refitting step
		of re-estimating the non-zero entries of the estimated parameters without L1 penalty. The purpose
		of this re-fitting step is to remove the effect of bias introduced by L1 penalty. The default value
		is option.ifrefit=true.
--  option.centered_input - a true/false flag to specify whether the input dataset has been mean-centered. The default
		value is false and the genotype data will be mean-centered. 
--  option.tol - tolerance for convergence of the optimization algorithm. Convgence is monitored
                 through the relative change of the absolute values of objective function 
                 over 10 iterations. Default value is 10^-7. 
--  option.verbose - a true/false flag that specifies whether to print out intermediate results to STDOUT. Default is false. 
--  option.eta - line search rate within the optimization algorithm. It must be a real number greater than 1. 
                 Large values of option.eta lead to fast line search. Default is 1.5.
--  option.Theta0 - initial value of parameters for the iterative optimization algorithm. This is a structure of two fields:
                 Theta0.xy (a J by K matrix) and Theta0.yy (a K by K matrix). Theta0.yy must be symmetric positive definite.  
		 Theta.xy and Theta.yy will be randomly initialized as sparse matrix and sparse positive definite matrix. 
                 For large datasets (e.g., K>1000), initializing Theta0.yy with graphical Lasso estimator 
		 (Friedman et al, 2008) may speed up convergence. 

Return Variable
--  OPT - a structure containing the following fields:
--  OPT.Theta - estimated sCGGM parameters, a structure of two fields: OPT.Theta.xy and OPT.Theta.yy. OPT.Theta.xy is 
		a J by K matrix that contains the estimated direct eQTL effects, and OPT.Theta.yy is a K by K matrix that 
		contains the network. 
--  OPT.intercept - the intercept term is a length K vector that corresponds to the intercept in the standard regression
		model. It represents the part of gene-expression level not explained by SNPs. Together with Theta.xy and Theta.yy, 
		OPT.intercept is used to predict gene-expression levels based on SNP data. 

------------------------------------------------- 
scggm_cv.m 	sparse CGGM with cross validation
------------------------------------------------- 

Description
--  Estimate a sparse CGGM with cross validation to select optimal
--  regularization parameters lambda1 and lambda2. 

Details
--  Regularization parameter lambda1 controls the amount of sparsity on
--  Theta.xy, and lambda2 controls the sparsity of Theta.yy. Optimal lambdas
--  are selected based on smallest cross-validation error. 

Usage
-- "OPT = scggm_cv(x, y, kcv, lambda1_seq, lambda2_seq, option)"

Arguments
--  x - N by J matrix of genotype data. 
--  y - N by K matrix of gene-expression data. 
--  kcv - the number of folds for cross-validation. Minimum value allowed is 3. In each fold, 
		(kcv-1)/kcv of the samples in the dataset is used as a training dataset to fit an sCGGM, 
		and the remaining 1/kcv of the samples is used as a validation dataset to calculate
          	cross-validation error. 
--  lambda<1/2>_seq (optional) - the vectors of regularization parameter lambda1 and lambda2's values. 
		Two vectors can have different lengths. Cross-validation will search all the 
		possible combinations of lambda1 and lambda2 values. 

		Default values of both vectors are [2^5, 2^4, ..., 2^1, 2^0] * 0.01. 
--  option (optional) - a structure that contains the parameter settings to fit an sCGGM, which is the
		same as in scggm.m above. If option.verbose=true, this routine will print 
		cross-validation errors to STDOUT.



Return Variable
--  OPT - a structure containing the following fields: 
--  OPT.Theta, OPT.intercept - sCGGM parameters, as in scggm.m return values.
--  OPT.lambdas - a length 2 vector that holds the optimal lambda1/lambda2 selected
                  by cross-validation. 


//////////////////////////////////////////////
Functions to Perform Inferences given an sCGGM
//////////////////////////////////////////////

--  Given a fitted sCGGM model, there are several inference schemes to
--  characterize SNP perturbations of gene-expression network in more detail,
--  in addition to direct SNP perturbations and gene-expression network.

--  1) overall indirect SNP perturbations - is characterized by a J by K matrix
--     Beta, obtained by routine scggm_indirect_SNP_overall().
--  2) decomposition of SNP perturbation effects - are summarized in a J by K 
--     matrix, obtained on a per-gene basis by routine
--     scggm_indirect_SNP_decompose(). 
--  3) decomposition of gene-expression covariance - contains three K by K
--     matrices, computed by scggm_cov_decompose().

------------------------------------------------------------
scggm_indirect_SNP_overall  perform inference given an sCGGM
------------------------------------------------------------
Description 
--  Perform inference of the graphical model given sCGGM parameters
--  to obtain the overall indirect SNP perturbations on gene-expression. 

Usage
--  "Beta = scggm_indirect_SNP_overall( Theta )"

Arguments
--  Theta - a fitted sCGGM parameter containing Theta.xy and Theta.yy. 

Return Variable
--  Beta - the overall indirect SNP perturbations, J by K matrix. 

------------------------------------------------------------------------------
scggm_indirect_SNP_decompose      decompose overall indirect SNP perturbations
------------------------------------------------------------------------------
Description
--  Compute decomposition of indirect SNP perturbations for a particular gene 
--  given an sCGGM model. 

Usage
--  "Beta_k = scggm_indirect_SNP_decompose(Theta, k)"

Arguments
--  Theta - sCGGM parameters Theta.xy and Theta.yy. 
--  k - the index of the gene that passes on the indrect effects. 

Return Variable
--  Beta_k - the indirect SNP effects passed on by the k-th gene, a 
             J by K matrix. 

-------------------------------------------------------------
scggm_cov_decompose      decompose gene-expression covariance
-------------------------------------------------------------
Description
--  Compute covariance decomposition of gene-expression given an sCGGM model. 

Usage
--  "Cov = scggm_cov_decompose(Theta, x, y, centered_input)"

Arguments
--  Theta - a fitted sCGGM parameter containing Theta.xy and Theta.yy. 
--  x - SNP data, N by J matrix. 
--  y - gene-expression data, N by K matrix.  
--  centered_input (optional) - indicates if the data is centered.
                                Default is false.

Return Variable
--  Cov - A structure containing three K by K matrix: Cov.Overall is the
--  overall gene-expression covariance; Cov.Network_Induced means the part
--  of covariance that is explained by the internal network structure
--  Theta_yy, and Cov.SNP_Induced means the part of covariance that is 
--  explained by SNP perturbations. 




///////////////////////////////////
Using an sCGGM model for prediction
///////////////////////////////////

----------------------------------------------------------------------------
scggm_predict.m 	predict gene-expression level with SNP data given an sCGGM
----------------------------------------------------------------------------
Description
--  Given an estimated sCGGM, predict gene-expression levels based on SNP data. 

Usage
--  "[Ey_ts, prederr] = scggm_predict(Theta, Intercept, x_ts, y_ts)"

Arguments
--  Theta - sCGGM parameter containing Theta.xy and Theta.yy. 
--  Intercept - sCGGM parameter for a length K vector of intercepts. 
--  x_ts - N_ts by J matrix of genotype data, where N_ts is the number of samples. 
		Gene-expression levels are predicted based on the SNPs in x_ts. 
--  y_ts (optional) - N_ts by K matrix of gene-expression levels. If y_ts is provided,
		it is used to compute prediction error between observed values y_ts
		and gene-expression levels predicted by the sCGGM. 

Return Variable
--  Ey_ts -- predicted gene-expression level. 
--  prederr -- prediction error. Computed only when y_ts is supplied. 




////////////////////////////////////////////////////////////////
Example code (demo) for using sCGGM package
////////////////////////////////////////////////////////////////

-------------------------------------------------------------------------------------------------
demo.m 		sample code for learning sparse CGGM given a pre-specified regularization parameter
		(without cross-validation for selecting regularization parameters), and performing inference 
-------------------------------------------------------------------------------------------------
Description
--  A demo code for running scggm.m (described above) to learn sCGGM parameters with a given
--  regularization parameter (without cross-validation for selecting the optimal regularization parameters).
--  The demo code also performs inference on the sCGGM to obtain the indirect SNP perturbation effects, 
--  the decomposition of indirect SNP perturbations, and the decomposition of gene-expression covariance matrix. 

Details
--  For the purpose of this demo, regularization parameters lambda_1 and lambda_2 were set to 0.1. 
--  Input data are read from ./data directory and results are saved to ./results/demo/ directory. 

Usage
--  ">demo": type "demo" on Matlab command prompt.

Input Directory Files
--  ./data/ 
    	1) xtrain.txt - N x J matrix of SNP data, where N is the number of samples, J 
          the number of SNPs. 
   	2) ytrain.txt - N x K matrix, gene-expression data, where 
	  K is the number of gene-expression traits.

Output Directory Files
--  ./results/demo/
	1) optimal_Theta_xy.txt - contains the estimated Theta_xy (J by K matrix). 
	2) optimal_Theta_yy.txt - contains the estimated Theta_yy (K by K matrix). 
	3) optimal_intercept.txt - estimated intercepts (a length K vector). 
	4) Beta.txt - indirect SNP perturbations (J by K matrix). 
	5) Beta_2.txt - decomposition of overall indirect SNP effects. The results are saved 
	  only for the component corresponding to the second gene-expression trait that passes on SNP 
	  perturbation effects. (J by K matrix). 
	6) Cov_Overall.txt - overall gene-expression covariance. 
	7) Cov_Network_Induced.txt - the component of gene-expression covariance 
	  induced by gene-expression network Theta.yy, after a decomposition of gene-expression covariance.
	8) Cov_SNP_Induced.txt - the component of gene-expression covariance induced
	  by SNP perturbations, after a decomposition of gene-expression covariance. 

----------------------------------------------------------------------
demo_cv.m 	sample code for learning sCGGM with cross-validation to select the optimal
	regularization parameters and performing inference on the estimated sCGGM
----------------------------------------------------------------------
Description
--  An example script for running sparse CGGM with five fold cross-validation
--  using scggm_cv.m (as described above), computing prediction error for test set, performing
--  inference on the sCGGM to obtain the indirect SNP perturbations, a decomposition of
--  indirect SNP perturbations, and a decomposition of gene-expression covariance.

Details
--  5-fold cross validation is used with the dataset of sample size N=300. Test data is used to compute 
--  prediction error. Results are saved to ./results/demo_cv directory. 

Usage
--  ">demo_cv": type "demo_cv" on Matlab command prompt.

Input Directory
--  ./data/ 
    	1) xtrain.txt - N x J matrix of SNP data for training set 
   	2) ytrain.txt - N x K matrix of gene-expression traits for training set
    	3) xtest.txt - N_ts x J matrix of SNP data for test set
    	4) ytest.txt - N_ts x K matrix of gene-expression data for test set, where
		N: the number of samples in the training set 
	  	N_ts: the number of samples in the test set
		J: the number of SNPs 
		K: the number of gene-expression traits

Output Directory
--  ./results/demo_cv/
	1) optimal_<Theta_xy, Theta_yy, intercept>.txt - sCGGM parameters, using the same format
           as in demo.m
        2) optimal_lambdas.txt - optimal regularization parameters lambda1 and 
            lambda2. 
	3) Beta.txt - overall indirect SNP perturbations, 
	4) Beta_2.txt - SNP perturbations after decomposition with respect to each gene-expression trait
	  that passes on the SNP perturbation effects. Results are saved only for the second gene-expression trait,
	  using the same format as in demo.m
	4) Cov_<Overall, Network_Induced, SNP_Induced>.txt - decomposition of
	   gene-expression covariance, using the same format as in demo.m 

See Also
--  demo.m



