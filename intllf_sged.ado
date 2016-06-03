/*This ado file gives the log likelihood function used in interval regressions
for the SGT distribution.
It works with intreg2.ado
v 1
Author--Jacob Orchard
Update--6/1/2016*/


program intllf_sged
version 13
		args lnf m lambda sigma p 
		tempvar Fu Fl zu zl 
		qui gen double `Fu' = .
		qui gen double `Fl' = .
		qui gen double `zu' = . 
		qui gen double `zl' = .
		
		*Point data
			tempvar x s l 
			qui gen double `x' = $ML_y1 - (`m') if $ML_y1 != . & $ML_y2 != . ///
												& $ML_y1 == $ML_y2
												
			qui gen double `s' = log(`p') - (abs(`x')^`p'/(`sigma'*(1+`lambda'*sign(`x')))^`p') ///
								if $ML_y1 != . & $ML_y2 != .  & $ML_y1 == $ML_y2
												
			qui gen double `l' = log(2) + log(`sigma') + lngamma((1/`p')) ///
								if $ML_y1 != . & $ML_y2 != . & $ML_y1 == $ML_y2
																					
			qui replace `lnf' = `s' - `l' if $ML_y1 != . & $ML_y2 != . & ///
								$ML_y1 == $ML_y2
		
		
		*Interval data
			qui replace `zu' = (abs($ML_y2 - `m')^`p')/( ///
								(`sigma'^`p')*(1+`lambda'*sign($ML_y2 -`m'))^`p') ///
								if $ML_y1 != . & $ML_y2 != . &  $ML_y1 != $ML_y2
								
			qui replace `Fu' = .5*(1-`lambda') + .5*(1+`lambda'*sign($ML_y2- ///
								`m'))*sign($ML_y2 - `m')*gammap(`zu', (1/`p')) ///
								if $ML_y1 != . & $ML_y2 != . &  $ML_y1 != $ML_y2
								
			qui replace `zl' = (abs($ML_y1 - `m')^`p')/( ///
								(`sigma'^`p')*(1+`lambda'*sign($ML_y1 -`m'))^`p') ///
								if $ML_y1 != . & $ML_y2 != . &  $ML_y1 != $ML_y2
								
			qui replace `Fl' = .5*(1-`lambda') + .5*(1+`lambda'*sign($ML_y1- ///
								`m'))*sign($ML_y1 - `m')*gammap(`zl', (1/`p'))  ///
								if $ML_y1 != . & $ML_y2 != . &  $ML_y1 != $ML_y2
								
			qui replace `lnf' = log(`Fu' -`Fl') if $ML_y1 != . & $ML_y2 != . &  ///
														$ML_y1 != $ML_y2
		
		*Bottom coded data
			qui replace `zl' = (abs($ML_y1 - `m')^`p')/( ///
								(`sigma'^`p')*(1+`lambda'*sign($ML_y1 -`m'))^`p') ///
								if $ML_y1 != . & $ML_y2 == .
								
			qui replace `Fl' = .5*(1-`lambda') + .5*(1+`lambda'*sign($ML_y1- ///
								`m'))*sign($ML_y1 - `m')*gammap(`zl', (1/`p')) ///
								if $ML_y1 != . & $ML_y2 == .
								
			qui replace `lnf' = log(1-`Fl') if $ML_y1 != . & $ML_y2 == .
		
		*Top coded data
			qui replace `zu' = (abs($ML_y2 - `m')^`p')/( ///
								(`sigma'^`p')*(1+`lambda'*sign($ML_y2 -`m'))^`p') ///
								if $ML_y2 != . & $ML_y1 == .
								
			qui replace `Fu' = .5*(1-`lambda') + .5*(1+`lambda'*sign($ML_y2- ///
								`m'))*sign($ML_y2 - `m')*gammap(`zu', (1/`p')) ///
								if $ML_y2 != . & $ML_y1 == .
								
			qui replace `lnf' = log(`Fu') if $ML_y2 != . & $ML_y1 == .
		
		*Missing values
			qui replace `lnf' = 0 if $ML_y2 == . & $ML_y1 == .
		
		
		
end		