/*This ado file gives the log likelihood function used in interval regressions
for the GB2 distribution.
It works with intreg2.ado
v 1.1
Author--Jacob Orchard
Update--6/1/2016*/


program intllf_gb2
version 13
		args lnf a b p q
		tempvar Fu Fl zu zl 
		qui gen double `Fu' = .
		qui gen double `Fl' = .
		qui gen double `zu' = . 
		qui gen double `zl' = .
		
		*Point data
			tempvar x y z
			
			qui gen double `x' = log( abs(`a')) + (`a'*`p' - 1)*log($ML_y1) ///
								if $ML_y1 != . & $ML_y2 != . & $ML_y1 == $ML_y2
								
			qui gen double `y' = (`a'*`p')*log(`b') + (lngamma(`p') + lngamma(`q') /// 
								- lngamma(`p' + `q')) if $ML_y1 != . & $ML_y2 != . ///
								& $ML_y1 == $ML_y2
								
			qui gen double `z' = (`p' + `q')*log(1+($ML_y1/`b')^`a') if $ML_y1 != . ///
									& $ML_y2 != . & $ML_y1 == $ML_y2
									
			qui replace `lnf' = `x' - `y' - `z' if $ML_y1 != . & $ML_y2 != .  ///
													& $ML_y1 == $ML_y2
			
		
		*Interval data
			qui replace `zu' = ($ML_y2/`b')^`a'/(1+($ML_y2/`b')^`a') if ///
								$ML_y1 != . & $ML_y2 != . &  $ML_y1 != $ML_y2
								
			qui replace `Fu' = exp(lngamma(`p')+lngamma(`q')-lngamma(`p'+`q'))* /// 
								ibeta(`p',`q',`zu') if $ML_y1 != . & $ML_y2 != . ///
								&  $ML_y1 != $ML_y2
								
			qui replace `zl' = ($ML_y1/`b')^`a'/(1+($ML_y1/`b')^`a') if /// 
								$ML_y1 != . & $ML_y2 != . &  $ML_y1 != $ML_y2
								
			qui replace `Fl' = exp(lngamma(`p')+lngamma(`q')-lngamma(`p'+`q'))* ///
								ibeta(`p',`q',`zl') if $ML_y1 != . & $ML_y2 != . ///
										&  $ML_y1 != $ML_y2
										
			qui replace `lnf' = log(`Fu' -`Fl') if $ML_y1 != . & $ML_y2 != . ///
								&  $ML_y1 != $ML_y2
		
		
		*Bottom coded data
			qui replace `zl' = ($ML_y1/`b')^`a'/(1+($ML_y1/`b')^`a') if ///
									$ML_y1 != . & $ML_y2 == .
									
			qui replace `Fl' = exp(lngamma(`p')+lngamma(`q')-lngamma(`p'+  ///
							`q'))*ibeta(`p',`q',`zl') if $ML_y1 != . & $ML_y2 == .
							
			qui replace `lnf' = log(1-`Fl') if $ML_y1 != . & $ML_y2 == .
		
		*Top coded data
		
			qui replace `zu' = ($ML_y2/`b')^`a'/(1+($ML_y2/`b')^`a') ///
								if $ML_y2 != . & $ML_y1 == .
								
			qui replace `Fu' = exp(lngamma(`p')+lngamma(`q')-lngamma(`p'+ ///
								`q'))*ibeta(`p',`q',`zu') if $ML_y2 != . & ///
									$ML_y1 == .
									
			qui replace `lnf' = log(`Fu') if $ML_y2 != . & $ML_y1 == .
			
		*Missing values
			qui replace `lnf' = 0 if $ML_y2 == . & $ML_y1 == .
		
		
		
end		
