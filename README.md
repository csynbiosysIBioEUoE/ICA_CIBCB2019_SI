# Information content analysis reveals desirable aspects of in vivo experiments of a synthetic circuit

Here the scripts used to generate, analyse and visualise the data presented in the paper are made available. 
The use of the scripts requires the RStan Package, available at the link:
https://cran.r-project.org/web/packages/rstan/index.html 
The data is organised in the following subfolders:

-	Inference
o	ODE_Model.stan, stan statistical model script used to perform Bayesian Inference of the experimental data from [1]. 
o	MultiExtractExp.R, script designed to access the experimental data and experimental schemes from [1] to generate an appropriate list of objects to be passed to the stan model to perform the inference. The csv files are generated using the script DataExtraction.m.
o	DataExtraction.m, script to extract the desired experimental data and experimental profiles from [1]. 
o	masterRun.R, script to perform inference through RStan using the designed model ODE_Model.stan and the list of data extracted from MultiExtractExp.R. The script allows to perform inference on single datasets in series or on the combined set. 
o	masterRunOptim.R, script designed as the masterRun.R script but including an initial optimisation process for the initialisation of the 4 MCMC chains used in the inference. 
-	PriorDefinition
o	ExtractingInitialPriorsLugagneLog.mat, matlab script to compute the mean and standard deviation of our priors (10 lognormal and 4 normal distributions) based on the results of the fit obtained in [1].
-	InformationContent
o	VarCovarMatrix.R, function that extracts the post-warmup samples from a list of selected stanfit objects and computes the covariance matrix of the posterior distribution of parameters, saving the determinant of the matrices of interest in a CSV file. 
o	RelativeEntropyFunction.R, function used to approximate the joint posterior from the samples of a stanfit object using Gaussian Mixtures and compute prior and approximate posterior entropy as presented in [2]. The function saves the Gaussian Mixture results as RDS objects and the summary of the entropies as a CSV file. 
o	runRelativeEntropy.R, script to run RelativeEntropyFunction.R for a set of stanfit object results selected in series.
-	Predictions&Analysis
o	ODE_Model_Function.stan, stan script containing the proposed ODE system and the implementation of the event-based representation of the inputs to simulate the response to a selected input (processed with the function MultiExtractExp.R).
o	PostPredCheckSimul.R, function to simulate the ODEs for a determined experimental profile selected using all the MCMC samples from a stanfit object result selected and save the results in CSV format. 
o	ConfidenceIntervalsParam.R, function to generate a plot for each parameter of the model with the 95% confidence interval from the MCMC samples for a set of stanfit objects.
o	ConfidenceIntervalPlotsFunction.R, function to extract all the MCMC samples from a stanfit object and simulate a determined experimental profile, obtaining the 95% confidence intervals and saving the plot of the simulation for checks. 
o	AccuracyPredictions.R, functions to compute the nRMSE distributions from the PostPredCheckSimul.R function and compute the Bhattacharyya distance between pairs nRMSEii and nRMSEij saving the results as a CSV matrix. 
o	SampleExtract.R, function to extract all the MCMC samples from a stanfit object result for a user-specified parameter, saving the results in CSV files. 
-	Sensitivity
o	SensitivityAnalyses, PDF document specifying where to locate the sensitivity index calculation function. 
o	SenseFileConv.R, function to convert the desired CSV prediction files obtained from PostPredCheckSimul.R and the MCMC samples from the associated simulation into compatible .mat files to be used by the function described in SensitivityAnalyses.
o	SensitivitiesSingle.mat, script to run the sensitivity index calculation for our set of MCMC samples and predictions extracted from SenseFileConv.R. 

To run the scripts plase place all the scripts and data required in the R working directory. 

The data associated with these scripts can be found at:

References:
[1] Jean-Baptiste Lugagne, Sebastián Sosa Carrillo, Melanie Kirch, Agnes Köhler, Gregory Batt & Pascal Hersen, 2017. Balancing a genetic toggle switch by real-time feedback control and periodic forcing. Nature Communications, 8 (1671), pp. 1-7.
[2] Huber, M. F., Bailey, T., Durrant-whyte, H. & Hanebeck, U. D., 2008. On Entropy Approxilllation for Gaussian Mixture Randolll Vectors. 2008 IEEE Int. Conf. Multisens. Fusion Integr. Intell. Syst. pp. 181–188.
