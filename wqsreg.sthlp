{smcl}
{* *! version 6.0 September 2025}{...}

{cmd:help wqsreg}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{hi:wqsreg} {hline 2}}Weighted Quantile Sum (WQS) regression{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 13 2}
{cmd:wqsreg} {it:yvar} {it:expvars} [{it:cvars}] , {opt mixture(varlist)} {opt boot(#)} [ {it:options} ]


{phang}
{it:yvar} is the name of the outcome variable. 

{phang}
{it:expvars} is the varlist containing the names of the variables representing the exposure.  

{phang}
[{it:cvars}] is the varlist containing the names of the variables representing the confounders.  


{synoptset 31 tabbed}{...}
{synopthdr :options}
{synoptline}
{p2coldent: * {opt mixture(varlist)}}specify the list of exposure variables{p_end}
{p2coldent: * {opt boot(#)}}set the number of bootstrap samples (>0); setting boot = 1 disables bootstrapping{p_end}

{synopt :{opt validation(integer)}}specify the validation set percentage (default: 0 = no split; range: [0, 100)){p_end}
{synopt :{opt q(integer)}}define the quantiles; the default is 4 (quartiles){p_end}
{synopt :{opt b1_neg(integer)}}set uni-directionality constraint: 0 = positive (default); 1 = negative{p_end}
{synopt :{opt cvar(varlist)}}specify the list of confounders{p_end}
{synopt :{opt seed(integer)}}set the random seed; default is 0{p_end}
{synopt :{opt conv_maxiter(integer)}}set the maximum number of iterations to be performed before optimization; default is 2000{p_end}
{synopt :{opt conv_vtol(real)}}set the tolerance; default is 0.000000001{p_end}
{synopt :{opt technique(string)}}choose the optimization method: 'bfgs' (Broyden–Fletcher–Goldfarb–Shanno, default) or 'nr' (Modified Newton–Raphson){p_end}
{synopt :{opt model_fam(string)}}choose the model: 'Linear', 'Poisson' or 'Logistic'{p_end}
{synopt :{opt saveWQSindex(integer)}}specify whether to save a dataset containing the WQS index for each observation: use 1 to save; the default is 0 (do not save). This option is not allowed for repeated holdout validation{p_end}
{synopt :{opt saveWeights(integer)}}specify whether to save a dataset containing the weights: use 1 to save; the default is 0 (do not save){p_end}
{synopt :{opt datasetWQSindexName(string)}}specify the name of the new dataset that will include the WQS index{p_end}
{synopt :{opt datasetWeightsName(string)}}specify the name of the new dataset that will include the weights{p_end}
{synopt :{opt figureName(string)}}specify the filename for the plot of weights, only if you want to save the figure{p_end}
{synopt :{opt id(string)}}specify the variable that identifies the observations{p_end}
{synopt :{opt rh_rep(integer)}}set the number of repetitions for repeated holdout validation; 1 is the default (no repeated holdout validation){p_end}
{p2colreset}{...}

{p 4 6 2}* are required.{p_end}


{title:Description}

{pstd}
Weighted Quantile Sum (WQS) regression is a flexible statistical method for quantifying the association between a set of possibly correlated predictors and a health outcome.
It allows estimating the overall effect of complex sets of exposures and the specific contributions of each factor.

{pstd}
{cmd:wqsreg} – enables users to fit WQS regression while allowing for the several flexible components of this framework
(e.g., confounders adjustment, direction specification, bootstrapping, training-validation data splitting and repeated holdout validation).

{pstd}
{cmd:wqsreg} returns the estimates from WQS regression, generates plots of the estimated weights, and saves additional related information. It requires Stata version 11 or higher.



{title:References}

{phang} Bellavia A., Statistical Methods for Environmental Mixtures, Springer, 2025. doi: 10.1007/978-3-031-78987-8

{phang} Carrico C, Gennings C, Wheeler DC, Factor-Litvak P. Characterization of Weighted Quantile Sum Regression for Highly Correlated Data in a Risk Analysis Setting.
J Agric Biol Environ Stat. 2015 Mar;20(1):100-120. doi: 10.1007/s13253-014-0180-3.
Epub 2014 Dec 24. PMID: 30505142; PMCID: PMC6261506.

{phang} Czarnota J, Gennings C, Wheeler DC. Assessment of weighted quantile sum regression for modeling chemical mixtures and cancer risk.
Cancer Inform. 2015 May 13;14(Suppl 2):159-71. doi: 10.4137/CIN.S17295.
PMID: 26005323; PMCID: PMC4431483.

{phang} Renzetti, S., Curtin, P., Just, A. C., Bello, G., Gennings, C., Renzetti, M. S., & Rsolnp, I. (2021). Package ‘gWQS’.


{title:Authors}

{pstd}Marta Ponzano [1][2], Stefano Renzetti [3], Andrea Bellavia [4]{p_end}

{pstd}[1] {it:Department of Health Sciences, University of Genoa, Genoa, Italy}{p_end}
{pstd}[1] {it:Department of Life Science, Health, and Health Professions, Link Campus University, Rome, Italy}{p_end}
{pstd}[3] {it:Department of Medicine and Surgery, University of Parma, Parma, Italy}{p_end}
{pstd}[4] {it:Department of Environmental Health, Harvard T.H. Chan School of Public Health}{p_end}
