/*This ado file executes non-linear interval regressions where the error term is 
distributed in the GB2 or SGT family tree

Author--Jacob Orchard
v 1.2
Update--6/8/2016*/



*program drop intreg2
program intreg2, eclass
version 13.0
	if replay() {
		display "Replay not implemented"
	}
	else {
		set more off
		syntax varlist(min=2 fv ts)  [aw fw pw iw] [if] [in] ///
		[, DISTribution(string) /// 
		sigma(varlist) ///
		lambda(varlist) ///
		p(varlist) ///
		q(varlist) ///
		b(varlist) ///
		beta(varlist)  ///
		INITial(numlist) ///
		vce(passthru)  ///
		eyx(string) ///
		CONSTraints(passthru) DIFficult TECHnique(passthru) ITERate(passthru)  /// 
		nolog TRace GRADient showstep HESSian SHOWTOLerance TOLerance(passthru) ///
		LTOLerance(passthru) NRTOLerance(passthru) robust cluster(passthru) ///
		svy SHOWConstonly] 
		
		*Defines Independent and Dependent Variables
		local depvar1: word 1 of `varlist'
		local depvar2: word 2 of `varlist'
		local tempregs: list varlist - depvar1 
		local regs: list tempregs - depvar2
				
		*Defines variables for other parameters
		if "`sigma;" != ""{
			local sigmavars `sigma'
			}
			
		if "`lambda;" != ""{
			local lambdavars `lambda'
			}
		if "`p;" != ""{
			local pvars `p'
			}
		if "`q;" != ""{
			local qvars `q'
			}
		
			
		local nregs: word count `regs'
		local nsigma: word count `sigmavars'
		local nlambda: word count `lambdavars'
		local np: word count `pvars'
		local nq: word count `qvars'
		
		
		*Displays error if using the wrong parameter with chosen distribution
		if  (`nlambda' > 0) & (("`distribution'" != "sgt") & ("`distribution'" != "sged")){
				di as err "Lambda is not a parameter of the chosen distribution"  
				exit 498 
			}
			
		if `np' > 0 & ("`distribution'" != "sgt" & "`distribution'" != "gb2" & "`distribution'" ///
								!= "gg" & "`distribution'" != "sged") {
					di as err "p is not a parameter of the chosen distribution"  
					exit 498
				}
		if  `nq' > 0 &  ("`distribution'" != "sgt" & "`distribution'" != "gb2") {

						di as err "q is not a parameter of the chosen distribution"
						exit 498 
				}
		
		*Defines titles used when running the program
	    local gb2title "Interval Regression with GB2 Distribution"
		local ggtitle "Interval Regression with Generalized Gamma Distribution"
	    local lntitle "Interval Regression with Log-Normal Distribution"
		local normaltitle "Interval Regression with Normal Distribution"
		local sgttitle "Interval Regression with SGT Distribution"
		
		*Warns users specifying inital values if too few inital values are specified.
		if "`initial'" != ""{			
			local initiallen: word count `initial'
			
			if ("`distribution'" == "gb2"  & `initiallen' != 4) {
				di as err "initial does not have the correct amount of numbers."
				di as err "You must have 4 starting guesses in initial for the GB2 distribution."
				exit 503
			}
			else if ("`distribution'" == "gg"  & `initiallen' != 3) {
				di as err "initial does not have the correct amount of numbers."
				di as err "You must have 3 starting guesses in initial for the Generalized Gamma distribution."
				exit 503
			}
			else if ("`distribution'" == "ln" | "`distribution'" == "lnormal"  & `initiallen' != 3) {
				di as err "initial does not have the correct amount of numbers."
				di as err "You must have 2 starting guesses in initial for the Log-Normal distribution."
				exit 503
			}
			else if ("`distribution'" == "sgt"   & `initiallen' != 5) {
				di as err "initial does not have the correct amount of numbers."
				di as err "You must have 5 starting guesses in initial for the SGT distribution."
				exit 503
			}
			else if ("`distribution'" == "sged"   & `initiallen' != 4) {
				di as err "initial does not have the correct amount of numbers."
				di as err "You must have 4 starting guesses in initial for the SGED distribution."
				exit 503
			}
			else if ("`distribution'" == "" | "`distribution'" == "normal" & `initiallen' != 2) {
				di as err "initial does not have the correct amount of numbers."
				di as err "You must have 2 starting guesses in initial for the Normal distribution."
				exit 503
			}
		}
		
		*Decides which observations to use in analysis.
		
		marksample touse, nov
		
		foreach i in  `regs' `sigmavars' `pvars' `qvars'{
			qui replace `touse' = 0 if `i' ==.
			}
		qui replace `touse' = 0 if `depvar1' == `depvar2' == .
		
		*Gets rid of uncensored observations with a non-positive dependent
		*variable if user is using a positive distribution.
	
		if "`distribution'" == "lnormal" | "`distribution'" == "ln" | ///
		 "`distribution'" == "gg" | "`distribution'" == "gb2"{
			quietly{ 
			  count if `depvar1' < 0 & `touse' & `depvar1' == `depvar2'
			  local n =  r(N) 
			  if `n' > 0 {
				noi di " "
				if `n' == 1{
					noi di as txt " {res:`depvar1'} has `n' uncensored value < 0;" _c
					noi di as text " not used in calculations"
				}
				else{
					noi di as txt " {res:`depvar1'} has `n' uncensored values < 0;" _c
					noi di as text " not used in calculations"
					}
				}

			  count if `depvar1' == 0 & `touse' & `depvar1' == `depvar2'
			  local n =  r(N) 
			  if `n' > 0 {
				noi di " "
				noi di as txt " {res:`depvar1'} has `n' uncensored values = 0;" _c
				noi di as text " not used in calculations"
				}
				
			  count if `depvar1' <= 0 & `depvar2' <= 0 & `touse' & `depvar1' != `depvar2'
			  local n =  r(N) 
			  if `n' > 0 {
				noi di " "
				noi di as txt " {res:`depvar1'} has `n' intervals < 0;" _c
				noi di as text " not used in calculations"
				}
				

		  replace `touse' = 0 if `depvar1' <= 0 & `depvar2' <= 0
		  }
		}
		
		*Counts the number of each type of interval
		quietly{
			count if `depvar1' != . & `depvar2' != . & `depvar1' == `depvar2'  /// 
			& `touse' == 1
			local nuncensored = r(N)
			count if `depvar1' != . & `depvar2' != . & `depvar1' != `depvar2'  ///
			& `touse' == 1
			local ninterval = r(N)
			count if `depvar1' != . & `depvar2' == .  & `touse' == 1
			local nright = r(N)
			count if `depvar1' == . & `depvar2' != . & `touse' == 1
			local nleft = r(N)
			count if `depvar1' == . & `depvar2' ==. & `touse' == 1
			local nnoobs = r(N)
			
		}
		
		*Evaluates model
		
		if "`distribution'" == "normal"{
			
			local evaluator intllf_normal
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2'=)  (sigma: ///
			) if `touse' ==1 , missing search(on)   maximize /// 
			initial(`initial') `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy'	 
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (sigma: ///
			`sigmavars') if `touse' ==1 , missing search(on) continue  maximize /// 
			initial(`initial') `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy'
			
			ml display
			
		}
		else if "`distribution'" == "lnormal" | "`distribution'" == "ln" {
			
			local evaluator intllf_ln
			
			di " "
			di as txt "Fitting constant only model:"
				
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = )  /// 
			(sigma: ) if `touse' ==1 , missing search(on)  /// 
			maximize initial(`initial') `constraints' `technique'  `difficult'  ///
			`iterate' `log' `trace' `gradient' `showstep' `hessian'  ///
			`showtolerance' `tolerance' `ltolerance' `nrtolerance' title(`lntitle') ///
			`vce'  `robust' `cluster' `svy' 
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  /// 
			(sigma: `sigmavars' ) if `touse' ==1 , missing search(on) continue /// 
			maximize initial(`initial') `constraints' `technique'  `difficult'  ///
			`iterate' `log' `trace' `gradient' `showstep' `hessian'  ///
			`showtolerance' `tolerance' `ltolerance' `nrtolerance' title(`lntitle') ///
			`vce'  `robust' `cluster' `svy'
			
			ml display
			
		}
		else if "`distribution'" == "sgt"{
			
			local evaluator intllf_sgt
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = )   ///
			(lambda: ) (sigma:  ) (p: ) (q:  ///
			) if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy'
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')   ///
			(lambda: `lambdavars') (sigma: `sigmavars' ) (p: `pvars') (q:  ///
			`qvars') if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy'
			
			ml display
			
		}
		else if "`distribution'" == "sged"{
			
			local evaluator intllf_sged
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (m: `depvar1' `depvar2' = )   ///
			(lambda: ) (sigma:  ) (p: )   ///
			 if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy'
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (m: `depvar1' `depvar2' = `regs')   ///
			(lambda: `lambdavars') (sigma: `sigmavars' ) (p: `pvars') [`weight'`exp']   ///
			 if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy'
			
			ml display
			
		}
		else if "`distribution'" == "gb2"{
			
			local evaluator intllf_gb2sigma
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = )   ///
			(sigma:  ) (p: ) (q:  ///
			) [`weight'`exp'] if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy'
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(sigma: `sigmavars') (p: `pvars') (q:  ///
			`qvars') [`weight'`exp'] if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy'
			
			ml display, plus
			
		}
		else if "`distribution'" == "gg"{
			
			local evaluator intllf_ggsigma
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = )   ///
			(sigma: ) (p: ) [`weight'`exp']  ///
			 if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///
			`robust' `cluster' `svy'
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (delta: `depvar1' `depvar2' = `regs')   ///
			(sigma: `sigmavars') (p: `pvars') [`weight'`exp']   ///
			 if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///
			`robust' `cluster' `svy'
			
			ml display
			
		}
		else{
			
			local evaluator intllf_normal
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2'=)  (sigma: ///
			)  [`weight'`exp']  if `touse' ==1 , missing search(on)   maximize /// 
			initial(`initial') `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy'	
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (mu: `depvar1' `depvar2' = `regs')  (sigma: ///
			`sigmavars')  [`weight'`exp']  if `touse' ==1 , missing search(on) continue  maximize /// 
			initial(`initial') `constraints' `technique'  `difficult' `iterate' ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`normaltitle') `vce' ///
			`robust' `cluster' `svy'
			
			ml display
		}
		
		mat betas = e(b) //coefficient matrix
		
		*Find the Conditional expected value at specified level
		if "`distribution'" == "gb2" | "`distribution'" == "gg" | "`distribution'" ///
							== "ln" | "`distribution'" == "lnormal"		{
			
			mat mid_Xs = 1
			if "`eyx'" == ""{
				local eyx "mean"
				di "{res:`eyx'}         {c |}"
			}
			else if "`eyx'" == "mean" {
			di "{res:`eyx'}         {c |}" 
			}
			else if "`eyx'" == "p50" | "`eyx'" == "p10" | "`eyx'" == "p25" | ///
			        "`eyx'" == "p75" | "`eyx'" == "p90" | "`eyx'" == "p95" | ///
					"`eyx'" == "p99" | "`eyx'" == "min" | "`eyx'" == "max" {		
			
					di "{res:`eyx'}          {c |}"
				}
			else if "`eyx'" == "p1" | "`eyx'" == "p5" {
				di "{res:`eyx'}           {c |}"
				}
			else{
				di as err "Not a valid option for eyx"
				exit 498
				}
			
			
			quietly foreach x in `regs' {
				sum `x', detail
				scalar mid_ = r(`eyx')
				mat mid_Xs = mid_Xs, mid_
			}
			mat sigma = betas[1,"sigma:_cons"]
			scalar sigma = sigma[1,1]
			
			if "`distribution'" == "gb2"{
			
				mat deltas = betas[1,"delta:"]
				mat deltas = deltas'
				mata: st_matrix("deltas", flipud(st_matrix("deltas"))) //flips matrix around
		                        									// to conform with Xs									
				mat p = betas[1,"p:_cons"]
				scalar p = p[1,1]
				mat q = betas[1,"q:_cons"]
				scalar q = q[1,1]
				mat xbeta = mid_Xs*deltas
				scalar xbeta = xbeta[1,1]
				mat expected = exp(xbeta)*( (exp(lngamma(p+sigma))*exp(lngamma(q-sigma)))/  ///
											( exp(lngamma(p))*exp(lngamma(q))))
			}
			
			if "`distribution'" == "gg"{
			
				mat deltas = betas[1,"delta:"]
				mat deltas = deltas'
				mata: st_matrix("deltas", flipud(st_matrix("deltas"))) //flips matrix around
																		// to conform with Xs										
				mat p = betas[1,"p:_cons"]
				scalar p = p[1,1]
				mat xbeta = mid_Xs*deltas
				scalar xbeta = xbeta[1,1]
				mat expected = exp(xbeta)*( (exp(lngamma(p+sigma)))/  ///
											( exp(lngamma(p))))
			}
			
			if "`distribution'" == "ln" | "`distribution'" == "lnormal" {
				
				mat mu = betas[1,"mu:"]
				mat mu = mu'
				mata: st_matrix("mu", flipud(st_matrix("mu"))) //flips matrix around
																		// to conform with Xs										
				mat xbeta = mid_Xs*mu
				scalar xbeta = xbeta[1,1]
				mat expected = exp((xbeta + sigma^2)/2)
			}
					
			scalar eyx = expected[1,1]
			table_line "E[Y|X]" eyx 
			di as text "{hline 13}{c BT}{hline 64}"
		
		}
		
		*Observation type count
		noi di " "
		if `nleft' != 1{
			noi di as txt " {res:`nleft'} left-censored observations" 
		}
		if `nleft' == 1{
			noi di as txt " {res:`nleft'} left-censored observation" 
		}
		if `nuncensored' != 1{
			noi di as txt " {res: `nuncensored'} uncensored observations" 
		}
		if `nuncensored' == 1{
			noi di as txt " {res:`nuncensored'} uncensored observation" 
		}
		if `nright' != 1{
			noi di as txt " {res:`nright'} right-censored observations" 
		}
		if `nright' == 1{
			noi di as txt " {res:`nright'} right-censored observation" 
		}
		if `ninterval' != 1{
			noi di as txt " {res:`ninterval'} interval observations" 
		}
		if `ninterval' == 1{
			noi di as txt " {res:`ninterval'} interval observation" 
		}
		
		qui ereturn list
		}
end

*program drop table_line
program table_line
	args vname coef se z p 95l 95h
	if (c(linesize) >= 100){
		local abname = "`vname'"
		}
	else if (c(linesize) > 80){
	local abname = abbrev("`vname'", 12+(c(linesize)-80))
	}
	else{
	local abname = abbrev("`vname'", 12)
	}
	local abname = abbrev("`vname'",12)
	display as text %12s "`abname'" " { c |}" /*
	*/ as result /*
	*/ " " %8.0g `coef' " " /*
	*/ %9.0g `se' " " %9.0g `z' " " /*
	*/ %9.0g `p' "  " %9.0g `95l' "   " /*
	*/ %9.0g `95h'
end

*mata: mata drop delta_sigma()
version 13
mata: 
	void delta_sigma( ) //Computes delta and sigma from a and
											// b for the gb2 distribution
	{
		as = st_matrix("As")
		bs = st_matrix("Bs")
		u = J(1,cols(as),1)
		Sigmas = u :/ as
		Deltas = log(bs)
		st_matrix("Sigmas",Sigmas)
		st_matrix("Deltas",Deltas)
	}
end

version 13.0
mata:
matrix function flipud(matrix X)
{
return(rows(X)>1 ? X[rows(X)..1,.] : X)
}
end
