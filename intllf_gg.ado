/*This ado file gives the log likelihood function used in interval regressions
for the GG distribution.
It works with intreg2.ado
v 1
Author--Jacob Orchard
Update--5/24/2016*/


program intllf_gg
version 13
		args lnf a beta p 
		tempvar Fu Fl zu zl 
		qui gen double `Fu' = .
		qui gen double `Fl' = .
		qui gen double `zu' = . 
		qui gen double `zl' = .
		
		*Point data
			tempvar x w y z
			
			qui gen double `x' = log( abs(`a')) + (`a'*`p' - 1)*log($ML_y1) ///
								if $ML_y1 != . & $ML_y2 != . & $ML_y1 == $ML_y2
								
			qui gen double `w' = ($ML_y1/`beta')^`a'
								
			qui gen double `y' = (`a'*`p')*log(`beta')  if $ML_y1 != . & $ML_y2 != . ///
								& $ML_y1 == $ML_y2
								
			qui gen double `z' = lngamma(`p') if $ML_y1 != . ///
									& $ML_y2 != . & $ML_y1 == $ML_y2
									
			qui replace `lnf' = `x' - `w' - `y' - `z' if $ML_y1 != . & $ML_y2 != .  ///
													& $ML_y1 == $ML_y2
			
		
		*Interval data
			qui replace `zu' = ($ML_y2/`beta')^`a' if $ML_y1 != . & $ML_y2 != . & ///
															$ML_y1 != $ML_y2
								
			qui replace `Fu' = gammap(`p',`zu') if $ML_y1 != . & $ML_y2 != . ///
								&  $ML_y1 != $ML_y2
								
			qui replace `zl' = ($ML_y1/`beta')^`a' if $ML_y1 != . & $ML_y2 != .  ///
												&  $ML_y1 != $ML_y2
								
			qui replace `Fl' = gammap(`p',`zl') if $ML_y1 != . & $ML_y2 != . ///
										&  $ML_y1 != $ML_y2
										
			qui replace `lnf' = log(`Fu' -`Fl') if $ML_y1 != . & $ML_y2 != . ///
								&  $ML_y1 != $ML_y2
		
		
		*Bottom coded data
			qui replace `zl' = ($ML_y1/`beta')^`a' if $ML_y1 != . & $ML_y2 == .
									
			qui replace `Fl' = gammap(`p',`zl') if $ML_y1 != . & $ML_y2 == .
							
			qui replace `lnf' = log(1-`Fl') if $ML_y1 != . & $ML_y2 == .
		
		*Top coded data
		
			qui replace `zu' = ($ML_y2/`beta')^`a'	if $ML_y2 != . & $ML_y1 == .
								
			qui replace `Fu' = gammap(`p',`zu') if $ML_y2 != . & $ML_y1 == .
									
			qui replace `lnf' = log(`Fu') if $ML_y2 != . & $ML_y1 == .
			
		*Missing values
			qui replace `lnf' = 0 if $ML_y2 == . & $ML_y1 == .
				
		
end		
