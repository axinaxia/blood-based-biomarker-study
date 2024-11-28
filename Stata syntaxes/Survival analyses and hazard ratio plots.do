/*Syntax for the study
"Blood-based biomarkers of Alzheimer's disease and incident dementia in the community: a 16-year-long study"

In this study, we had two main objectives:

To estimate the hazard of all-cause and Alzheimer's disease (AD) dementia in relation to baseline levels of AD biomarkers.
To evaluate the predictive performance of these biomarkers for all-cause and AD dementia within a 10-year timeframe.

AIM 1: We conducted Cox proportional hazards regression models using calendar time as the time scale. These analyses were performed using Stata, version 17.

AIM 2: The search for optimal cut-offs and the evaluation of biomarker performance using these cut-offs were done in R. For further details, please refer to the R script provided.

The syntax below outlines the analysis for p-tau217 toward all-cause dementia (dem). To analyze other biomarkers, replace P-tau217 with the relevant biomarker name. To analyze AD dementia replace time_to_dementia with time_to_AD and dem with AD_dem. 
*/

* Survival analyses 

* Figure 3

global confounder "age educ3trnew AnyE4_2014 sex BLAnemia BLAtrial_fibrillation BLCerebrovascular_dis BLChronic_kidney_dis BLHeart_failure BLHypertension BLIschemic_heart_dis BLObesity"
su nfl ttau  gfap ptau181 abetaratio ptau217
global biomarkers "nfl gfap ttau ptau181 ptau217 abetaratio" 
tabstat $biomarkers, stat(min mean sd max )
local plot ""
local j = 1
foreach b of global biomarkers {
	su `b', d 
	mkspline `b's = `b', nk(3) cubic displayknots 
	mat knot_`b' = r(knots)
	stcox  `b's1  `b's2 $confounder

	qui su `b's1 
	scalar ref1 = r(min) // pick the minimum value as referent 
	qui su `b's2 if `b's1 == scalar(ref1)
	scalar ref2 = r(mean)
	
	predictnl loghrs_`b' = _b[`b's1]*(`b's1-scalar(ref1)) + _b[`b's2]*(`b's2-scalar(ref2)) , ci(lb_`b' ub_`b')
	gen hrs_`b' = exp(loghrs_`b')
	gen lbs_`b' = exp(lb_`b')
	gen ubs_`b' = exp(ub_`b')
	qui su  `b', d
	local min = r(min)
	local max = r(p95)
	*local plot "`plot' (line hrs_`b' `b' if inrange(`b', `min', `max') , sort lw(thick) ) "
	*local legend "`legend' label(`j' `b')"
	local j = `j' + 1
}

capture drop y_cases
gen y_cases = 32+3
gen pipe = "|"

su ptau217 ,d
scalar min = r(min)
scalar max = r(p95)
set scheme s2color
twoway ///
(kdensity ptau217, lc(green%10) fc(green%10) yaxis(2) recast(area)) ///
(line hrs_ptau217 ptau217 , sort lw(thick) lc(black)) ///
(rarea lbs_ptau217 ubs_ptau217 ptau217, sort fc(grey%10) lc(grey%10)) ///
(scatter y_cases ptau217 if _d == 1, jitter(5) ms(i) mlabel(pipe) mlabcolor(black%50) ) ///
if inrange(ptau217, min, max) ,  yscale(log) ///
yscale(axis(2) off) yscale(axis(1) alt) ///
ylabel(0.5 1 2 4 8 16 32 64, format(%3.1fc))  ///
xtitle("Ptau-217") legend(off)  yscale(range(0.5 16))  ///
plotregion(fcolor(white)) graphregion(fcolor(white)) ///
xlabel(#10) ///
ytitle("Adjusted HR") yline(1) ysize(4) xsize(4)    name(f_ptau217, replace )

graph export nfl_dem_white.jpg, as(jpg) name("f_ptau217") quality(100) replace

* Table 2

xtile ptau217_cat = z_ptau217, nq(4)
stset time_to_dementia, failure(dem) id(Lopnr)
stcox i.ptau217_cat age sex educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis AnyE4_2014 
stptime , by(ptau217_cat) per(100)
stcox i.ptau217_cat age sex educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis AnyE4_2014 if mmse >=27
stcox i.ptau217_cat age sex educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis AnyE4_2014 if agegr2 ==0
stcox i.ptau217_cat age sex educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis AnyE4_2014 if agegr2 ==1
stcox i.ptau217_cat age  educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis AnyE4_2014 if sex ==0
stcox i.ptau217_cat age  educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis AnyE4_2014 if sex ==1
stcox i.ptau217_cat age sex educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis  if AnyE4_2014 ==0
stcox i.ptau217_cat age sex educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis  if AnyE4_2014 ==1
stcox i.ptau217_cat age sex educ3trnew BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLAnemia BLObesity BLChronic_kidney_dis BLCerebrovascular_dis  if subj_complaints ==1

* Inverse probability weighting 
stset time_to_dementia, failure(dem) id(Lopnr)
gen miss = _st==0

//estimating predictive probability of missing 

ttest age, by(miss)
tab sex miss,chi col
tab educ3trnew miss,chi col

*"Chron_num" is the number of chronic diseases
oneway Chron_num  miss, t 
ttest Chron_num, by(miss)

logistic miss  age sex i.educ3trnew Chron_num     

predict pmiss, pr // Predicted probability of missing
gen weight_censor = 1/(1 - pmiss) if miss == 0 //censoring weight 

stset time_to_dementia [pweight=weight_censor], failure(dem) id(Lopnr) 

stcox i.ptau217_cat age sex educ3trnew  BLCerebrovascular_dis BLHeart_failure BLIschemic_heart_dis BLAtrial_fibrillation BLHypertension BLAnemia BLObesity BLChronic_kidney_dis  AnyE4_2014 
