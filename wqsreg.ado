
/*
wqsreg: fitting Weighted Quantile Sum (WQS) regression models for continuous outcomes
*/

capture program drop wqsreg 
program define wqsreg  , rclass

 syntax varlist(max=1 numeric) ///
 [, validation(integer 60) q(integer 4) b1_neg(integer 0) cvar(varlist) seed(integer 0) conv_maxiter(integer 2000) conv_vtol(real 0.000000001) technique(string) model_fam(string) savingWQSindex(string)  savingWeights(string) figureName(string) id(string) rh_rep(integer 1)  boot(integer 100)] ///
 mixture(varlist)   
 
 local n_mixt : word count `mixture'
 scalar n_mixt = `n_mixt'

   
 if (`rh_rep'<1 ) {
  di as err "Error, please set the number of repetitions for repeated holdout validation (rh_rep); 1 is the default (no repeated holdout validation)"
  exit 198
  } 
  
   if (`rh_rep'==1 ) {
   di as txt "Note: for a more robust estimation the use of repeated holdouts is recommended"
  } 
  
  if (`rh_rep'!=1 & "`savingWQSindex'"!="") {
								  	        di as err "Error, savingWQSindex is not allowed for repeated holdout validation (rh_rep)"
                                            exit 198
								                         }
  
  if (`rh_rep'==1 | `rh_rep'==. ) {
  wqsreg_no_rh `varlist'  , validation(`validation') q(`q') b1_neg(`b1_neg') cvar(`cvar') /// 
                                       seed(`seed') conv_maxiter(`conv_maxiter') conv_vtol(`conv_vtol') ///
                                       technique(`technique') model_fam(`model_fam') ///
                                        savingWeights(`savingWeights') savingWQSindex(`savingWQSindex')   ///
                                       figureName(`figureName') id(`id')  ///
                                       mixture(`mixture') boot(`boot')
  } 
  
   if (`rh_rep'>1) {
   	
     local start = `seed'
     local end   = `seed' + `rh_rep'-1
     
	 matrix coef_WQS_indexM = J(`rh_rep', 1, .)
	 matrix WeightsM = J(`rh_rep',`n_mixt', .)
     forvalues s = `start'/`end' {	   	
										
	 	                               local index = `s' - `start' + 1	
									   
									    
	                                   di "......repeated holdout `index'......."
									   capture noisily wqsreg_no_rh `varlist'  , validation(`validation') q(`q') b1_neg(`b1_neg') cvar(`cvar') /// 
                                                                 seed(`s') conv_maxiter(`conv_maxiter') conv_vtol(`conv_vtol') ///
                                                                 technique(`technique') model_fam(`model_fam') ///
                                                                 savingWQSindex(`savingWQSindex') savingWeights(`savingWeights') ///
                                                                 figureName(`figureName') id(`id')  ///
                                                                 mixture(`mixture') boot(`boot')
										
																				
										if (_rc==1) {   // user break
                                                         di as err "Interrupted by the user"
                                                         exit 198
                                                        }
                                      else if (_rc!=0) {
											         di as txt "Warning - repeated holdout `index' not used"
										            }						 
									   
									  if (_rc==0)   {
									  	matrix coef_WQS_indexM[`s'-`start'+1,1]=coef_WQS_index
									 									    
									   forvalues i = 1/`n_mixt' {
									                             scalar val`i' = mean_cols[1,`i']
									                             matrix WeightsM[`s'-`start'+1,`i']=val`i'
									                            }     
																} 
                                  }
								  

								  capture frame drop fTemp
								  frame create fTemp
                                  frame change fTemp
								  
								  qui clear 
								  qui set obs `rh_rep'
								  								
								  local mylist "`mixture'"
                                  forvalues j = 1/`n_mixt'{
								                        qui gen Weights`j' = .
								                        forvalues i = 1/`rh_rep' {
                                                                                 qui replace Weights`j' = WeightsM[`i', `j'] in `i' 
																				 if ("`savingWeights'"!=""){
																				 	                          savingWeights, filename("`savingWeights'")
																					                       }		 
								                                                 }
														 }
														 														 					 				     
														
														 
														 gen rep = _n
														 qui reshape long Weights, i(rep) j(mixt)														 
														 
														 local i = 1
                                                         foreach lbl in `mixture' {
                                                                                   label define mixtlbl `i' "`lbl'", add
                                                                                   local ++i
                                                                                     }		 
														 label values mixt mixtlbl
														 
														 local threshold = 1/`n_mixt'
														 
														 qui graph box Weights 	, ///
                                                                   over(mixt, sort(1)  descending label(angle(90))) title("Weights") ytitle("") box(1, lcolor(gs9))   ///
                                                                   yline(`threshold', lcolor(red) lpattern(dash)) 
																   
														 qui graph save box_Weights.gph, replace 
								 
								  
								  frame change default
								  
								  
								  capture frame drop fTemp
								  frame create fTemp
                                  frame change fTemp
								  
								  qui clear 
								  qui set obs `rh_rep'
								  
								  qui gen coef_WQS_index_var = .
								  forvalues i = 1/`rh_rep' {
                                                             qui replace coef_WQS_index_var = coef_WQS_indexM[`i', 1] in `i'			                             
								                          }
								 qui histogram coef_WQS_index_var ,   frequency bin(10) title("WQS index, coefficient")	color(gs9) xtitle("")
								 
								 qui graph save coef_WQS_index_var.gph, replace
								 
								 display ""
								 display ""
                                 display  "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
								 display ""
								 display "Repeated Holdout validation - WQS index coefficient: "
								 sum 	coef_WQS_index_var	
								 local n = r(N)
								 display ""
								 display ""
								 display "				95%-sd CI: [" r(mean)-1.96*r(sd) ", " r(mean)+1.96*r(sd) "]"
								 qui centile coef_WQS_index_var, centile(2.5 97.5)
                                 local p2_5 = r(c_1)
                                 local p97_5 = r(c_2) 
								 display "				2.5-97.5th : ["`p2_5' ", "`p97_5' "]"
								 
								  
								 frame change default
								 
								 graph combine coef_WQS_index_var.gph box_Weights.gph, title("rh - `n' repetitions")

								 if "`figureName'" != ""  {
		                                                   qui graph save "`figureName'.gph" , replace
	                                                      }
								 								 							 
								 erase box_Weights.gph
								 erase coef_WQS_index_var.gph			  
			} 			  
  end

/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*--*-*-*-*-*-*/
  
capture program drop wqsreg_no_rh
program define wqsreg_no_rh , rclass
 
 version 11
 
 syntax varlist ///
 [, validation(integer 60) q(integer 4) b1_neg(integer 0) cvar(varlist) seed(integer 0) conv_maxiter(integer 2000) conv_vtol(real 0.000000001) technique(string) model_fam(string) savingWQSindex(string) savingWeights(string)   figureName(string) id(string) boot(integer 100) ] ///
 mixture(varlist)  
 
 qui frame copy default fWQS, replace /*create & copy current data to new frame*/
 frame change fWQS                // switch to the new frame
 
  /*conv_maxiter*/
  if (`conv_maxiter'<=0) {
     di as err "Error, maximum number of iterations to be performed before optimization must be positive"
     exit 198
  } 
  /*conv_vtol*/
   if (`conv_vtol'<=0) {
     di as err "Error, conv_vtol must be positive"
     exit 198
  } 
  /*technique*/
  if ("`technique'" != "nr" & "`technique'" != "bfgs" & "`technique'" != "") {
     di as err "Error, please choose a valid optimization technique — either Broyden–Fletcher–Goldfarb–Shanno (bfgs) or modified Newton–Raphson(nr)"
     exit 198
  } 
  if ("`technique'" == "") {
     local technique "bfgs"
                           }
  /*model_fam*/					   
  if ("`model_fam'" != "Linear" & "`model_fam'" != "Poisson" & "`model_fam'" != "Logistic" & "`model_fam'" != "") {
       di as err "Error, please choose a valid model — either Linear or Poisson or Logistic)"
       exit 198
  } 
  if ("`model_fam'" == "") {
        local model_fam "Linear"
                           }						   
 
  foreach var of varlist `varlist' {
       qui drop if `var'==.
  }
   if ("`cvar'"!="") {
	foreach var of varlist `cvar' {
       qui drop if `var'==.
  }
  }
    foreach var of varlist `mixture' {
       qui drop if `var'==.
  }
  
  /*b1_neg*/
  if (`b1_neg'!=0 & `b1_neg'!=1) {
        di as err "Error, b1_neg can be 0 (positive b1, the default) or 1 (negative b1)"
        exit 198
  } 
  
  /* savingWQSindex*/
if ("`savingWQSindex'"!="" & "`id'"=="") {
        di as err "Error, id is required when savingWQSindex is reported"
        exit 198
  } 
  
  if ("`id'"!="") {
        confirm variable `id'
  } 
 
  /*datasetWeightsName; datasetWQSindexName*/
  if ("`savingWeights'"=="`savingWQSindex'" & "`savingWeights'"!="") {
        di as err "Error, please provide different names for savingWQSindex and savingWeights"
        exit 198
  } 
 
  /*boot*/
  if (`boot'<=0) {
           di as err "Error, please insert a postive number of bootstrap samples. Note that boot=1 means no bootstrapping"
           exit 198
  } 

  /*validation*/
  if (`validation'<0 | `validation'>=100) {
            di as err "Error, Validation must be numeric in [0; 100)"
            exit 198
  } 
  else {
            set seed `seed'
            quietly gen Validation=rbinomial(1,`validation'/100)
  }
  
  qui gen Y_wqs = `varlist'

  local n_mixt : word count `mixture'
  scalar n_mixt = `n_mixt'

  local n_conf : word count `cvar'
  scalar n_conf = `n_conf'
  
  if n_conf>0   {
                 quietly: ds `cvar'  
                 local adj = r(varlist)
                 foreach var of varlist `adj' {
                                           qui gen Z_conf_wqs_`var' = `var' 
                                           }
                 }
  if n_conf == 0 {
                  quietly gen Z_conf = . 
                 }

  
  foreach var of varlist `mixture' {
                                    xtile Q_wqs_`var' = `var', nq(`q')
                                   }   
  
  matrix W_boot = J(1, `n_mixt', .)
 
  if `boot'==1 { 
                
                preserve
				 
				
                quietly  drop if Validation==1 
  
                
                capture noisily  WQS_single_boot Y_wqs Q_wqs_* `cvar',   b1_neg(`b1_neg')  n_conf(`n_conf') conv_maxiter(`conv_maxiter') conv_vtol(`conv_vtol') technique(`technique') model_fam(`model_fam')
		        
		        if (_rc!=0) {
											         di as txt "Warning - bootstrap sample not used"
													  }
				
				else {
					if( /*p_WQS<0.05 &*/  b_WQS>0 & b1_neg==0){
			                                          matrix W_boot =  W
		                                            }
		        if( /*p_WQS<0.05 &*/  b_WQS<0 & b1_neg==1){
			                                         matrix W_boot =  W
		                                           }	
												   }
               restore
               }
  else         {
              	forvalues i = 1/`boot' {
					
                                        preserve
										
                                        quietly  drop if Validation==1 
                                        bsample 
                                        capture noisily  WQS_single_boot  Y_wqs Q_wqs_* `cvar',  b1_neg(`b1_neg')  n_conf(`n_conf') conv_maxiter(`conv_maxiter') conv_vtol(`conv_vtol') technique(`technique')  model_fam(`model_fam')
										
										if (_rc!=0) {
											         di as txt "Warning - bootstrap sample `i' not used"
										            }
										
                                        else {
											  
										if( /*p_WQS<0.05 &*/  b_WQS>0 & b1_neg==0){
             	             	                                              matrix W`i' = W
             	             	                                              matrix W_boot = W_boot \ W`i'
             	             	                                              }
                                        if(/*p_WQS<0.05 &*/ b_WQS<0 & b1_neg==1){
             	             	                                              matrix W`i' = W
             	             	                                              matrix W_boot = W_boot \ W`i'
             	             	                                              }
										}
                                        restore
										
										if _rc == 1 {   
                                                     di as err "Interrupted by the user"
                                                      exit 1	
													}
	                                    }
              }
			  
   local rows = rowsof(W_boot)
   local cols = colsof(W_boot)
   
   if `rows'>1 {
	            matrix W_boot = W_boot[2..`rows', 1..`cols']
               }
   
   else	if (`rows'==1) {
  	                    scalar x = W_boot[1,1]
		                if (x==.) {
		                             
		                            if( `b1_neg'==1){
			                                        di as err "Error: There are no negative b1"	
		                                            exit 198
		                                            }
		                            if( `b1_neg'==0){
			                                        di as err "Error: There are no positive b1"
			                                        exit 198
			                                        }
		                          }
	                    }
						
  
  local rows = rowsof(W_boot)
  qui gen WQS_index = .
  
  mata{	
       W_boot = st_matrix("W_boot");
       mean_cols = mean(W_boot);
       st_matrix("mean_cols", mean_cols) ;
	   
       Q = st_data(., "Q_wqs*");
       Y_wqs = st_data(., "Y_wqs");
       N = st_nobs()

       WQS_index=Q*mean_cols';
	   
       st_store(., "WQS_index", WQS_index);
	}
	
    preserve
    clear
	
    qui svmat mean_cols, names(col)
	local mixture "`mixture'"
	local i=1
	foreach component of local mixture {
	                                 	rename c`i'   `component'
										local ++i
		                                }									
									
	if ("`savingWeights'"!=""){
		                    savingWeights, filename("`savingWeights'")
	                      }									

  restore
  
  display ""
  display ""
  display  "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
  display "N observations used - Total: " _N
  di      "Number of bootstrap samples used: "     "`rows'"
  
  preserve
  
  if (`validation'!=0) {
	                   quietly  keep if Validation==1
                       }
					   
  if n_conf>0  {
  	            if "`model_fam'"=="Linear" {
										  	qui regress Y_wqs WQS_index  Z_conf_wqs_*
											}
				if "`model_fam'"=="Poisson" {
										  	qui poisson Y_wqs WQS_index  Z_conf_wqs_*
											}
				if "`model_fam'"=="Logistic" {
										  	qui logit   Y_wqs WQS_index  Z_conf_wqs_*
											}
               } 
  if n_conf==0  {
  	            if "`model_fam'"=="Linear" {
        	                                qui regress Y_wqs WQS_index  
											}
  	            if "`model_fam'"=="Poisson" {
        	                                qui poisson Y_wqs WQS_index  
											}
  	            if "`model_fam'"=="Logistic" {
        	                                qui logit   Y_wqs WQS_index  
											}
                }
			 
  matrix tbl_results = r(table)
  
  scalar coef_WQS_index  = tbl_results[1,1]
  scalar se_WQS_index    = tbl_results[2,1]
  scalar pval_WQS_index  = tbl_results[4,1]
  scalar ci_lo_WQS_index = tbl_results[5,1]
  scalar ci_hi_WQS_index = tbl_results[6,1]
  scalar Nobs = e(N)

			 
  display "N observations used - Validation: " Nobs
  display "WQS index Coef: " coef_WQS_index
  display "Std. Err.: " se_WQS_index
  display "p-value: " pval_WQS_index
  display "95% CI: [" ci_lo_WQS_index ", " ci_hi_WQS_index "]" 
  display "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
			 
  restore

  /*Plot weights*/
  preserve
  clear
  qui set obs `n_mixt'
  
  local threshold = 1/`n_mixt'

  local mylist "`mixture'"
  qui  gen varname = ""
  
  local i = 1
  foreach component of local mylist {
                                     qui  replace varname = "`component'" in `i'
                                     local ++i
                                    }

   
  qui  gen weight = .
  qui  forvalues i=1/`n_mixt' {
                               replace weight = mean_cols[1, `i']  in `i'
                               }
 
  graph hbar weight,   over(varname, sort(1) descending label(angle(0))) bar(1, color(gs9)) title("Weights") ytitle("") /*blabel(bar, format(%5.2f))*/ yline(`threshold', lcolor(gs9) lpattern(dash))
  
  if "`figureName'" != ""  {
		                     qui graph save "`figureName'.gph" , replace
	                       }
  restore

  if "`savingWQSindex'"!="" {
  	                       savingWQSindex , filename("`savingWQSindex'")    id("`id'") 
                     }
 
  frame change default  // switch back to previous frame
  frame drop fWQS       // drop new frame
  
end

/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

capture program drop WQS_single_boot
program define WQS_single_boot , rclass

  version 11

  syntax varlist [, b1_neg(integer 0) conv_maxiter(integer 2000) conv_vtol(real 0.000000001) technique(string) model_fam(string) ]  n_conf(integer)  
 
  global  model_fam  "`model_fam'"
 
  quietly gen wqs=0 
  
  if n_conf>0 {
  	           if "`model_fam'"=="Linear" {
                                            qui regress Y_wqs Z_conf_wqs_* wqs
										   }  
  	           if "`model_fam'"=="Poisson" {
                                            qui poisson Y_wqs Z_conf_wqs_* wqs
										   }
  	           if "`model_fam'"=="Logistic" {
                                            qui logit   Y_wqs Z_conf_wqs_* wqs
										   }
               local model=r(varlist)    
              }
  if n_conf  == 0 {
                   if "`model_fam'"=="Linear" {
                                               qui regress  Y_wqs   wqs 
											   } 
  	               if "`model_fam'"=="Poisson" {
                                               qui poisson  Y_wqs   wqs 
											   } 
  	               if "`model_fam'"=="Logistic" {
                                                qui logit    Y_wqs   wqs 
											    }  
                   local model=r(varlist)     
                  }
  matrix init_reg=e(b)
	
  if (`b1_neg'==0)     {
                         matrix init_reg[1,n_conf+1]=0.0001
                       } 
  else if (`b1_neg'==1){
                         matrix init_reg[1,n_conf+1]=-0.0001
                       }
  scalar b1_neg = `b1_neg'
  scalar conv_maxiter=`conv_maxiter'
  scalar conv_vtol=real("`conv_vtol'")
   
  scalar b_cons = init_reg[1,n_conf+2]
  scalar b_wqs  = init_reg[1,n_conf+1]
  if n_conf>0  {
  	            matrix b_confounders= init_reg[1,1..n_conf]
               }
  if n_conf==0 {
  	            matrix b_confounders=0
               }
  
  quietly gen WQS_index_var = .
 
  quietly mata{
	
                S  = optimize_init()  ;
                 
                n_mixt = st_numscalar("n_mixt");
                n_conf = st_numscalar("n_conf");
				N = st_nobs()
				
				conv_maxiter = st_numscalar("conv_maxiter");
				conv_vtol = st_numscalar("conv_vtol");
				technique=st_local("technique")

                b_cons = st_numscalar("b_cons");
                b_wqs = st_numscalar("b_wqs");
                b_confounders=st_matrix("b_confounders")

                initial_w=J(1, n_mixt, 1/sqrt(n_mixt)) 
                params = (b_wqs,b_cons, initial_w, b_confounders)
  
                if (n_conf == 0) {
	                              params = (b_wqs,b_cons, initial_w)
                                 }

                optimize_init_which(S, "min");
                optimize_init_technique(S, technique);
                optimize_init_evaluator(S, &Optim());
				optimize_init_trace_params(S);
                optimize_init_params(S, (params));
                optimize_init_conv_maxiter(S, conv_maxiter);
                optimize_init_conv_vtol(S, conv_vtol);
                optimize_init_conv_warning(S, "off" );
			
                bh = optimize(S); 
				

                bh = bh'
	            			
				den = sum(bh[3::(n_mixt+2), 1] :^ 2)
				W = (bh[3::(n_mixt+2), 1] :^ 2) / den
                W = W'   
                    
	            st_matrix("W", W)
                Q = st_data(., "Q_wqs_*")
	            Y = st_data(., "Y_wqs")
		
	            WQS_index = Q*W'
				
	             st_store(., "WQS_index_var", WQS_index)
	           }
	
	if (n_conf >0) {
		            if "`model_fam'"=="Linear" {
					                              qui regress Y_wqs WQS_index Z_conf_wqs_*
						                         }
					if "`model_fam'"=="Poisson" {
					                              qui poisson Y_wqs WQS_index Z_conf_wqs_*
						                         }		
					if "`model_fam'"=="Logistic" {
					                              qui logit Y_wqs WQS_index Z_conf_wqs_*
						                         }
                     }
	if (n_conf == 0) {
		              if "`model_fam'"=="Linear" {
					  	                           qui regress Y_wqs WQS_index
						                         }
		              if "`model_fam'"=="Poisson" {
					  	                           qui poisson Y_wqs WQS_index
						                         }
		              if "`model_fam'"=="Logistic" {
					  	                           qui logit Y_wqs WQS_index
						                         }
                     }
	
	matrix model=e(b) 
	scalar b_WQS = model[1,1]
	
	matrix p_table = r(table)
    scalar p_WQS = p_table[4, 1]
end  

/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/

mata: mata clear
mata

  void Optim(real scalar todo, real vector b,    ///
               val, grad, hess)
    {
	  real vector y, xb, Interc, wqs_, diff2, w
	  real matrix Q , Z 
	 
	  model_fam = st_global("model_fam")

	  y = st_data(., "Y_wqs")
	  Q = st_data(., "Q_wqs_*")
	  Z = st_data(., "Z_conf*")
	  
	  
	  n_mixt = st_numscalar("n_mixt")
	  n_conf = st_numscalar("n_conf")
	  
	  N = st_nobs()
	
	  Interc = J(N, 1, 1)
	  
	  real scalar den  
	  						
	  w = b[1, 3::(3+n_mixt-1)]:^2 / sum(b[1, 3::(3+n_mixt-1)]:^2)
	  wqs_ = Q * w'						
				
	  xb=Interc*b[1,2]+(wqs_)*b[1,1] 
	  
	
	  if (n_conf>0) {
	                  xb=xb+Z*(b[1,n_mixt+2+1..2+n_mixt+n_conf])'
	                }
	
     mu = J(N, 1, .)
	 if (model_fam == "Logistic") {  
	 	                             mu=1:/(1:+exp(-xb))
                                     for (j=1; j<=N; j++) { 
									 	                   if (mu[j]<0.000000000001) {
																                       mu[j]=0.000000000001
															                          } 
															 if (mu[j]>1-0.000000000001) {
																                       mu[j]=1-0.000000000001
															                          } 
									                       }    
	                                } 
									
	 lambda = J(N, 1, .)
	 if (model_fam == "Poisson") {	
	 	                           lambda=exp(xb)
	                             }
	 diff2 = J(N, 1, 0)
      
	 if (model_fam == "Linear") {
							 	  diff2=(xb:-y):^2 
							      } 
	
	 if (model_fam == "Logistic") {
	  	
		                            diff2= (y:==0) :* (-(2 :* ( log(1:-mu)))) + ///
							 	           (y:==1) :* (-(2 :* (log(mu))))
							        } 	
									
	  if (model_fam == "Poisson") {
							 	    diff2=-(-lambda + y:*ln(lambda) - ln(factorial(y)))
							       } 	
                        				   
     val = sum(diff2) 

	}
end

capture program drop savingWQSindex
program define savingWQSindex 

        syntax [varlist] , filename(string) id(string)
		
		preserve
                       keep `id' Validation WQS_index 	
					   qui save "`filename'.dta"  , replace      

        restore
end

capture program drop savingWeights
program define savingWeights 

        syntax [varlist] , filename(string) 
					   qui save "`filename'.dta"  , replace      
end





/*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/