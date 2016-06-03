/*This ado file executes non-linear interval regressions where the error term is 
distributed in the GB2 or SGT family tree

Author--Jacob Orchard
v 1
Update--5/25/2016*/




program intreg2, eclass
version 13.0
	if replay() {
		display "Replay not implemented"
	}
	else {
		set more off
		syntax varlist(min=2 fv ts) [if] [in] [, DISTribution(string) /// 
		sigma(varlist) ///
		lambda(varlist) p(varlist) q(varlist) b(varlist) beta(varlist)  ///
		INITial(numlist) vce(passthru)  ///
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
		if "`b;" != ""{
			local bvars `b'
			}
		if "`beta;" != ""{
			local bvars `beta'
			}
			
		local nregs: word count `regs'
		local nsigma: word count `sigmavars'
		local np: word count `pvars'
		local nq: word count `qvars'
		
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
			
			
		}
		else if "`distribution'" == "sgt"{
			
			local evaluator intllf_sgt
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (m: `depvar1' `depvar2' = )   ///
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
			
			ml model lf `evaluator' (m: `depvar1' `depvar2' = `regs')   ///
			(lambda: `lambdavars') (sigma: `sigmavars' ) (p: `pvars') (q:  ///
			`qvars') if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy'
			
			
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
			(lambda: `lambdavars') (sigma: `sigmavars' ) (p: `pvars')   ///
			 if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`sgttitle') `vce'  ///
			`robust' `cluster' `svy'
			
			
		}
		else if "`distribution'" == "gb2"{
			
			local evaluator intllf_gb2
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (a: `depvar1' `depvar2' = )   ///
			(b: ) (p: ) (q:  ///
			) if `touse' ==1 , missing search(on) maximize  ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy'
			
			if "`showconstonly'" != ""{
				ml display, showeqns //Shows constant only model
			}
			
			di " "
			di as txt "Fitting Full model:"
			
			ml model lf `evaluator' (a: `depvar1' `depvar2' = `regs')   ///
			(b: `bvars') (p: `pvars') (q:  ///
			`qvars') if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`gb2title') `vce'  ///
			`robust' `cluster' `svy'
			
			
		}
		else if "`distribution'" == "gg"{
			
			local evaluator intllf_gg
			
			di " "
			di as txt "Fitting constant only model:"
			
			ml model lf `evaluator' (a: `depvar1' `depvar2' = )   ///
			(beta: ) (p: )   ///
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
			
			ml model lf `evaluator' (a: `depvar1' `depvar2' = `regs')   ///
			(beta: `betavars') (p: `pvars')   ///
			 if `touse' ==1 , missing search(on) maximize continue ///
			initial(`initial') `constraints' `technique'  `difficult' `iterate'  ///
			`log' `trace' `gradient' `showstep' `hessian' `showtolerance'  ///
			`tolerance' `ltolerance' `nrtolerance' title(`ggtitle') `vce'  ///
			`robust' `cluster' `svy'
			
			
		}
		else{
			
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
			
		}
		
		ml display
		
		*Additional output
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
