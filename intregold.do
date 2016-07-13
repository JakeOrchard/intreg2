*Old Stuff for intreg2

*Additional output
		
		if "`distribution'" == "gb2"{
		
			*Extract parameter estimates
			mat As = betas[1,"a:"]
			mat Bs = betas[1,"b:"]
			mat Ps = betas[1,"p:"]
			mat Qs = betas[1,"q:"]

			
			
			if `nsigma' == 0 & `nregs' == 0{
				
				mata: As = st_matrix("As")
				mata: Bs = st_matrix("Bs")
				mata: delta_sigma() 
				
				mat colnames Deltas = delta:
				mat colnames Deltas = _cons 
				mat colnames Sigmas = sigma:
				mat colnames Sigmas = _cons 
				
				*Initial value matrix
				mat new_initial = Deltas,Sigmas,Ps,Qs
				mat rownames new_initial = y1
				
				*Compute Standard errors for Delta and Sigma
				ml model lf intllf_gb2sigma (delta: `depvar1' `depvar2' = `regs')  (sigma: ///
				`sigmavars') (p: `pvars') (q: `qvars')  [`weight'`exp'] if `touse' == 1  ///
				, missing  `constraints' 

				qui ml check
				
				qui ml search

				qui ml init new_initial
				
				ml maximize, nowarn  /// 
				 `vce' 	`robust' `cluster' `svy'
			
			ml display
			
			}		
		
		
		}
		
		
		
		
		
		
		
		
		Conditional mean for positive distributions
		if "`distribution'" == "gb2"{
			
			*A parameter
			mat As = betas[1,"a:"]
			mat out_table = r(table)
			mat Ase = out_table["se","a:"]
			mat A95l = out_table["ll","a:"]
			mat A95h = out_table["ul","a:"]
			scalar asize = colsof(As)
			global anum = asize
			
			*B parameter
			mat Bs = betas[1,"b:"]
			mat Bse = out_table["se","b:"]
			mat B95l = out_table["ll","b:"]
			mat B95h = out_table["ul","b:"]
			scalar bsize = colsof(Bs)
			
			*CRIT
			mat crit = out_table["crit",1]
			local crit = crit[1,1]
			
			if `nregs' == 0{
				di "{res: delta}" "       {c |}"
				
				local deltaval: di %7.5g = log(Bs[1,1])
				local bse = Bse[1,1]
				*fix standard errors
				local deltase: di %7.5g = `bse' * `deltaval'
				local deltall:  di %7.5g = `deltaval' - `crit'*`deltase'
				local deltaul: di %7.5g = `deltaval' + `crit'*`deltase'
					
				table_line _cons `deltaval' `deltase' . . `deltall' `deltaul' 
					
				
				di as text "{hline 13}{c +}{hline 64}"
			}
			
			if `nsigma' == 0{
				di "{res: sigma}" "       {c |}"
				
				local sigmaval: di %7.5g = 1/As[1,1]
				local ase = Ase[1,1]
				*Fix standard errors
				local sigmase: di %7.5g = `ase' * `sigmaval'
				local sigmall:  di %7.5g = `sigmaval' - `crit'*`sigmase'
				local sigmaul: di %7.5g = `sigmaval' + `crit'*`sigmase'
					
				table_line _cons `sigmaval' `sigmase' . . `sigmall' `sigmaul' 
					
				
				di as text "{hline 13}{c BT}{hline 64}"
			}
			
			
			
			
			
		}
		
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
