/***************************************************************************************************/
/***************************************************************************************************/
/****************************  Make Data for Programs and Analysis  ********************************/
/***************************************************************************************************/
/***************************************************************************************************/

*  Notes
*  For each year, we make 4 files:  one for all cohorts, one for young cohort, one for middle cohort, and one for older cohort
*  For each year-cohort, there are three files.  
*  The reason we need three files is we need to allocate part time workers to both the home sector with 0.5 weight and their respective market occ with 0.5 weight.
*  To fascilitate this, we made two occ_codes in the above file:  occ_code and occ_code2
*  occ_code2 puts the people working part time into the home sector (with 0.5 weight)
*  occ_code  puts the people working part time into the market sector (with 0.5 weight)
*  The third file for each year-cohort is merging these two back together
*  In total, there will be 12 files for each year (three files for each cohort (including the "all age" cohort)).
*  Finally, for each cohort, we will transform the composite data into longform for exporting. 
*  After these steps, we simply combine the cohort-year data into one file for analysis 

	
/*************************************************************************************************/
/****************************************** Cohort All *******************************************/
/*************************************************************************************************/
	
	#delimit ;
	
	foreach i in 2012 2000 1990 1980 1970 1960 {

	/* Step 1:  Compute the averages for all occupations - include the part time people as 0.5 in their paid occupation (using occ_code)*/ ; 									

	cd c:\ErikMain\discrimination_growth;
	use `i'_extract_composite_main, clear ; 	

	/*   A:  Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 
			  
	keep if cohort == 1 | cohort == 2 | cohort == 3 ;  /*  This is the line that changes across the cohort files */
	
	/*  B:  Define variables for counts, earnings, education and wages by occupation*group */
	
	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;	
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt], by(occ_code) ;																	
								
	drop if occ_code == 0;
								
	save `i'_extract_cohort_all_step1, replace ;								
										

	/* Step 2:  Compute the averages for all occupations - include the part time people as 0.5 in home sector (using occ_code2)*/ ; 									
									
	#delimit ;
	
	use `i'_extract_composite_main, clear ; 	

	/*   A:  Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 

	keep if cohort == 1 | cohort == 2 | cohort == 3 ;  /*  This is the line that changes across the cohort files */
	
	
	/*  B:  Define variables for counts, earnings, education and wages by occupation*group */

	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt],  by(occ_code2) ;					

	keep if occ_code2 == 0;   /*  This is key difference from step 1 file - we only keep the home sector here */
	rename occ_code2 occ_code ; 		 

	save `i'_extract_cohort_all_step2, replace ;	

	
	/* Step 3:  Merge step1 and step2 files.  Then extract key variables.  Reshape data for easier use. */ ; 									
	
	use 			`i'_extract_cohort_all_step1; 
	append using	`i'_extract_cohort_all_step2;

	/* Define variables */
	
	*  Note "0" refers to all groups
	*  Note "1"  refers to white men
	*  Note "2"  refers to white women	
	*  Note "3"  refers to white black men	
	*  Note "4"  refers to white black women	
	*  These are the code that chad used for running our analysis file
	
	*  Note "occ_income" is average earnings for full time workers in a given occupation-group
	*  Note "occ_grade" is average years of schooling for workers in a given occupation-group
	*  Note "wage" is the average hourly wage for full time workers in a given occupation-group
	*  Note "num" is the number of people in a given occupation-group. 
	
	*  Note:  The program is set to make both geometric and arithmatic means of variables ;

	gen year = `i' ;
	gen cohort = 0 ; 
	
	gen num0 	= count_num_adj;
	gen num1  	= white_man_adj;
	gen num2	= white_woman_adj;
	gen num3  	= black_man_adj;
	gen num4	= black_woman_adj;
	
	gen occ_income0 = income_full_all/working_full;
	gen occ_income1 = income_full_wm/working_full_wm;
	gen occ_income2 = income_full_ww/working_full_ww;
	gen occ_income3 = income_full_bm/working_full_bm;
	gen occ_income4 = income_full_bw/working_full_bw;

	gen occ_ln_income_arith0 = ln(occ_income0);
	gen occ_ln_income_arith1 = ln(occ_income1);
	gen occ_ln_income_arith2 = ln(occ_income2);
	gen occ_ln_income_arith3 = ln(occ_income3);
	gen occ_ln_income_arith4 = ln(occ_income4);
	
	gen occ_ln_income_geo0 = ln_income_full_all/working_full;
	gen occ_ln_income_geo1 = ln_income_full_wm/working_full_wm;
	gen occ_ln_income_geo2 = ln_income_full_ww/working_full_ww;
	gen occ_ln_income_geo3 = ln_income_full_bm/working_full_bm;
	gen occ_ln_income_geo4 = ln_income_full_bw/working_full_bw;
	
	gen occ_grade0 = highgrade_adj/count_num_adj;
	gen occ_grade1 = highgrade_adj_wm/white_man_adj;
	gen occ_grade2 = highgrade_adj_ww/white_woman_adj;
	gen occ_grade3 = highgrade_adj_bm/black_man_adj;
	gen occ_grade4 = highgrade_adj_bw/black_woman_adj;

	gen occ_wage0 = . ;
	gen occ_wage1 = wage_wm/working_full_wm;
	gen occ_wage2 = wage_ww/working_full_ww;
	gen occ_wage3 = wage_bm/working_full_bm;
	gen occ_wage4 = wage_bw/working_full_bw;

	gen occ_ln_wage_arith0 = ln(occ_wage0);
	gen occ_ln_wage_arith1 = ln(occ_wage1);
	gen occ_ln_wage_arith2 = ln(occ_wage2);
	gen occ_ln_wage_arith3 = ln(occ_wage3);
	gen occ_ln_wage_arith4 = ln(occ_wage4);
	
	gen occ_ln_wage_geo0 = . ;
	gen occ_ln_wage_geo1 = ln_wage_wm/working_full_wm;
	gen occ_ln_wage_geo2 = ln_wage_ww/working_full_ww;
	gen occ_ln_wage_geo3 = ln_wage_bm/working_full_bm;
	gen occ_ln_wage_geo4 = ln_wage_bw/working_full_bw;
	
	gen occ_var_ln_income0 = (ln_income_full_all_sd)^2;
	gen occ_var_ln_income1 = (ln_income_full_wm_sd)^2;
	gen occ_var_ln_income2 = (ln_income_full_ww_sd)^2;
	gen occ_var_ln_income3 = (ln_income_full_bm_sd)^2;
	gen occ_var_ln_income4 = (ln_income_full_bw_sd)^2;	

	gen occ_var_ln_wage0 = . ;
	gen occ_var_ln_wage1 = (ln_wage_wm_sd)^2;
	gen occ_var_ln_wage2 = (ln_wage_ww_sd)^2;
	gen occ_var_ln_wage3 = (ln_wage_bm_sd)^2;
	gen occ_var_ln_wage4 = (ln_wage_bw_sd)^2;
	
	
	/*  Keep the main variables for our analysis file*/
	
	keep 	year cohort occ_code num0-num4 occ_income0-occ_income4 occ_grade0-occ_grade4  occ_wage0-occ_wage4 
			occ_ln_income_arith0-occ_ln_income_arith4 occ_ln_income_geo0-occ_ln_income_geo4
			occ_ln_wage_arith0-occ_ln_wage_arith4 occ_ln_wage_geo0-occ_ln_wage_geo4
			occ_var_ln_income0-occ_var_ln_income4  occ_var_ln_wage0-occ_var_ln_wage4 ;

	reshape long num occ_income occ_grade occ_wage 
			occ_ln_income_arith occ_ln_income_geo occ_ln_wage_arith occ_ln_wage_geo
			occ_var_ln_income occ_var_ln_wage, i(occ_code) j(group);
	
	sort group occ_code ; 
	
	save `i'_extract_cohort_all_final, replace ; 

	
	
 
 
/*************************************************************************************************/
/****************************************** Cohort Young *****************************************/
/*************************************************************************************************/
								

	/* Step 1:  Compute the averages for all occupations - include the part time people as 0.5 in their paid occupation (using occ_code)*/ ; 									

	#delimit ;  
	cd c:\ErikMain\discrimination_growth;
	use `i'_extract_composite_main, clear ; 	

	/*  A. Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 
			  
	keep if cohort == 1 ;  /*  This is the line that changes across the cohort files */
	
	
	/*  B. Define variables for counts, earnings, education and wages by occupation*group */
	
	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;	
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt], by(occ_code) ;																					
								
	drop if occ_code == 0;
								
	save `i'_extract_cohort_young_step1, replace ;								
										

	/* Step 2:  Compute the averages for all occupations - include the part time people as 0.5 in home sector (using occ_code2)*/ ; 									
									
	#delimit ;
	
	use `i'_extract_composite_main, clear ; 	

	/*   A. Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 

	keep if cohort == 1 ;  	/*  This is the line that changes across the cohort files */
	
	
	/*  B. Define variables for counts, earnings, education and wages by occupation*group */

	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt],  by(occ_code2) ;																

	keep if occ_code2 == 0;   /*  This is key difference from step 1 file - we only keep the home sector here */
	rename occ_code2 occ_code ; 		 

	save `i'_extract_cohort_young_step2, replace ;	

	
	/* Step 3:  Merge step1 and step2 files.  Then extract key variables.  Reshape data for easier use. */ ; 									

	use 			`i'_extract_cohort_young_step1; 
	append using	`i'_extract_cohort_young_step2;

	/* Define variables */
	
	*  Note "0" refers to all groups
	*  Note "1"  refers to white men
	*  Note "2"  refers to white women	
	*  Note "3"  refers to white black men	
	*  Note "4"  refers to white black women	
	*  These are the code that chad used for running our analysis file
	
	*  Note "occ_income" is average earnings for full time workers in a given occupation-group
	*  Note "occ_grade" is average years of schooling for workers in a given occupation-group
	*  Note "wage" is the average hourly wage for full time workers in a given occupation-group
	*  Note "num" is the number of people in a given occupation-group. 
	
	*  Note:  The program is set to make both geometric and arithmatic means of variables ;
	
	gen year = `i' ;
	
	gen 	cohort = . ; 
	replace cohort = 1 if `i' == 2012;
	replace cohort = 2 if `i' == 2000;
	replace cohort = 3 if `i' == 1990;
	replace cohort = 4 if `i' == 1980;
	replace cohort = 5 if `i' == 1970;
	replace cohort = 6 if `i' == 1960;
	
	gen num0 	= count_num_adj;
	gen num1  	= white_man_adj;
	gen num2	= white_woman_adj;
	gen num3  	= black_man_adj;
	gen num4	= black_woman_adj;
	
	gen occ_income0 = income_full_all/working_full;
	gen occ_income1 = income_full_wm/working_full_wm;
	gen occ_income2 = income_full_ww/working_full_ww;
	gen occ_income3 = income_full_bm/working_full_bm;
	gen occ_income4 = income_full_bw/working_full_bw;

	gen occ_ln_income_arith0 = ln(occ_income0);
	gen occ_ln_income_arith1 = ln(occ_income1);
	gen occ_ln_income_arith2 = ln(occ_income2);
	gen occ_ln_income_arith3 = ln(occ_income3);
	gen occ_ln_income_arith4 = ln(occ_income4);
	
	gen occ_ln_income_geo0 = ln_income_full_all/working_full;
	gen occ_ln_income_geo1 = ln_income_full_wm/working_full_wm;
	gen occ_ln_income_geo2 = ln_income_full_ww/working_full_ww;
	gen occ_ln_income_geo3 = ln_income_full_bm/working_full_bm;
	gen occ_ln_income_geo4 = ln_income_full_bw/working_full_bw;
	
	gen occ_grade0 = highgrade_adj/count_num_adj;
	gen occ_grade1 = highgrade_adj_wm/white_man_adj;
	gen occ_grade2 = highgrade_adj_ww/white_woman_adj;
	gen occ_grade3 = highgrade_adj_bm/black_man_adj;
	gen occ_grade4 = highgrade_adj_bw/black_woman_adj;

	gen occ_wage0 = . ;
	gen occ_wage1 = wage_wm/working_full_wm;
	gen occ_wage2 = wage_ww/working_full_ww;
	gen occ_wage3 = wage_bm/working_full_bm;
	gen occ_wage4 = wage_bw/working_full_bw;

	gen occ_ln_wage_arith0 = ln(occ_wage0);
	gen occ_ln_wage_arith1 = ln(occ_wage1);
	gen occ_ln_wage_arith2 = ln(occ_wage2);
	gen occ_ln_wage_arith3 = ln(occ_wage3);
	gen occ_ln_wage_arith4 = ln(occ_wage4);
	
	gen occ_ln_wage_geo0 = . ;
	gen occ_ln_wage_geo1 = ln_wage_wm/working_full_wm;
	gen occ_ln_wage_geo2 = ln_wage_ww/working_full_ww;
	gen occ_ln_wage_geo3 = ln_wage_bm/working_full_bm;
	gen occ_ln_wage_geo4 = ln_wage_bw/working_full_bw;
	
	gen occ_var_ln_income0 = (ln_income_full_all_sd)^2;
	gen occ_var_ln_income1 = (ln_income_full_wm_sd)^2;
	gen occ_var_ln_income2 = (ln_income_full_ww_sd)^2;
	gen occ_var_ln_income3 = (ln_income_full_bm_sd)^2;
	gen occ_var_ln_income4 = (ln_income_full_bw_sd)^2;	

	gen occ_var_ln_wage0 = . ;
	gen occ_var_ln_wage1 = (ln_wage_wm_sd)^2;
	gen occ_var_ln_wage2 = (ln_wage_ww_sd)^2;
	gen occ_var_ln_wage3 = (ln_wage_bm_sd)^2;
	gen occ_var_ln_wage4 = (ln_wage_bw_sd)^2;
	
	
	/*  Keep the main variables for our analysis file*/
	
	keep 	year cohort occ_code num0-num4 occ_income0-occ_income4 occ_grade0-occ_grade4  occ_wage0-occ_wage4 
			occ_ln_income_arith0-occ_ln_income_arith4 occ_ln_income_geo0-occ_ln_income_geo4
			occ_ln_wage_arith0-occ_ln_wage_arith4 occ_ln_wage_geo0-occ_ln_wage_geo4
			occ_var_ln_income0-occ_var_ln_income4  occ_var_ln_wage0-occ_var_ln_wage4 ;

	reshape long num occ_income occ_grade occ_wage 
			occ_ln_income_arith occ_ln_income_geo occ_ln_wage_arith occ_ln_wage_geo
			occ_var_ln_income occ_var_ln_wage, i(occ_code) j(group);
	
	sort group occ_code ; 
	
	save `i'_extract_cohort_young_final, replace ; 

	
	
/*************************************************************************************************/
/****************************************** Cohort Mid *******************************************/
/*************************************************************************************************/
								

	/* Step 1:  Compute the averages for all occupations - include the part time people as 0.5 in their paid occupation (using occ_code)*/ ; 									

	#delimit ;  
	cd c:\ErikMain\discrimination_growth;
	use `i'_extract_composite_main, clear ; 	

	/*   A. Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 
			  
	keep if cohort == 2 ;  /*  This is the line that changes across the cohort files */
	
	/*  B. Define variables for counts, earnings, education and wages by occupation*group */
	
	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;	
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt], by(occ_code) ;																							
								
	drop if occ_code == 0;
								
	save `i'_extract_cohort_mid_step1, replace ;								
										

	/* Step 2:  Compute the averages for all occupations - include the part time people as 0.5 in home sector (using occ_code2)*/ ; 									
									
	#delimit ;
	
	use `i'_extract_composite_main, clear ; 	

	/*   A. Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 

	keep if cohort == 2 ;  	/*  This is the line that changes across the cohort files */
	
	/*  B. Define variables for counts, earnings, education and wages by occupation*group */

	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt],  by(occ_code2) ;														

	keep if occ_code2 == 0;   /*  This is key difference from step 1 file - we only keep the home sector here */
	rename occ_code2 occ_code ; 		 

	save `i'_extract_cohort_mid_step2, replace ;	

	
	/* Step 3:  Merge step1 and step2 files.  Then extract key variables.  Reshape data for easier use. */ ; 									

	use 			`i'_extract_cohort_mid_step1; 
	append using	`i'_extract_cohort_mid_step2;

	/* Define variables */
	
	*  Note "0" refers to all groups
	*  Note "1"  refers to white men
	*  Note "2"  refers to white women	
	*  Note "3"  refers to white black men	
	*  Note "4"  refers to white black women	
	*  These are the code that chad used for running our analysis file
	
	*  Note "occ_income" is average earnings for full time workers in a given occupation-group
	*  Note "occ_grade" is average years of schooling for workers in a given occupation-group
	*  Note "wage" is the average hourly wage for full time workers in a given occupation-group
	*  Note "num" is the number of people in a given occupation-group. 
	
	*  Note:  The program is set to make both geometric and arithmatic means of variables ;

	gen year = `i' ;

	gen 	cohort = . ; 
	replace cohort = 2 if `i' == 2012;
	replace cohort = 3 if `i' == 2000;
	replace cohort = 4 if `i' == 1990;
	replace cohort = 5 if `i' == 1980;
	replace cohort = 6 if `i' == 1970;
	replace cohort = 7 if `i' == 1960;
	
	gen num0 	= count_num_adj;
	gen num1  	= white_man_adj;
	gen num2	= white_woman_adj;
	gen num3  	= black_man_adj;
	gen num4	= black_woman_adj;
	
	gen occ_income0 = income_full_all/working_full;
	gen occ_income1 = income_full_wm/working_full_wm;
	gen occ_income2 = income_full_ww/working_full_ww;
	gen occ_income3 = income_full_bm/working_full_bm;
	gen occ_income4 = income_full_bw/working_full_bw;

	gen occ_ln_income_arith0 = ln(occ_income0);
	gen occ_ln_income_arith1 = ln(occ_income1);
	gen occ_ln_income_arith2 = ln(occ_income2);
	gen occ_ln_income_arith3 = ln(occ_income3);
	gen occ_ln_income_arith4 = ln(occ_income4);
	
	gen occ_ln_income_geo0 = ln_income_full_all/working_full;
	gen occ_ln_income_geo1 = ln_income_full_wm/working_full_wm;
	gen occ_ln_income_geo2 = ln_income_full_ww/working_full_ww;
	gen occ_ln_income_geo3 = ln_income_full_bm/working_full_bm;
	gen occ_ln_income_geo4 = ln_income_full_bw/working_full_bw;
	
	gen occ_grade0 = highgrade_adj/count_num_adj;
	gen occ_grade1 = highgrade_adj_wm/white_man_adj;
	gen occ_grade2 = highgrade_adj_ww/white_woman_adj;
	gen occ_grade3 = highgrade_adj_bm/black_man_adj;
	gen occ_grade4 = highgrade_adj_bw/black_woman_adj;

	gen occ_wage0 = . ;
	gen occ_wage1 = wage_wm/working_full_wm;
	gen occ_wage2 = wage_ww/working_full_ww;
	gen occ_wage3 = wage_bm/working_full_bm;
	gen occ_wage4 = wage_bw/working_full_bw;

	gen occ_ln_wage_arith0 = ln(occ_wage0);
	gen occ_ln_wage_arith1 = ln(occ_wage1);
	gen occ_ln_wage_arith2 = ln(occ_wage2);
	gen occ_ln_wage_arith3 = ln(occ_wage3);
	gen occ_ln_wage_arith4 = ln(occ_wage4);
	
	gen occ_ln_wage_geo0 = . ;
	gen occ_ln_wage_geo1 = ln_wage_wm/working_full_wm;
	gen occ_ln_wage_geo2 = ln_wage_ww/working_full_ww;
	gen occ_ln_wage_geo3 = ln_wage_bm/working_full_bm;
	gen occ_ln_wage_geo4 = ln_wage_bw/working_full_bw;
	
	gen occ_var_ln_income0 = (ln_income_full_all_sd)^2;
	gen occ_var_ln_income1 = (ln_income_full_wm_sd)^2;
	gen occ_var_ln_income2 = (ln_income_full_ww_sd)^2;
	gen occ_var_ln_income3 = (ln_income_full_bm_sd)^2;
	gen occ_var_ln_income4 = (ln_income_full_bw_sd)^2;	

	gen occ_var_ln_wage0 = . ;
	gen occ_var_ln_wage1 = (ln_wage_wm_sd)^2;
	gen occ_var_ln_wage2 = (ln_wage_ww_sd)^2;
	gen occ_var_ln_wage3 = (ln_wage_bm_sd)^2;
	gen occ_var_ln_wage4 = (ln_wage_bw_sd)^2;
	
	
	/*  Keep the main variables for our analysis file*/
	
	keep 	year cohort occ_code num0-num4 occ_income0-occ_income4 occ_grade0-occ_grade4  occ_wage0-occ_wage4 
			occ_ln_income_arith0-occ_ln_income_arith4 occ_ln_income_geo0-occ_ln_income_geo4
			occ_ln_wage_arith0-occ_ln_wage_arith4 occ_ln_wage_geo0-occ_ln_wage_geo4
			occ_var_ln_income0-occ_var_ln_income4  occ_var_ln_wage0-occ_var_ln_wage4 ;

	reshape long num occ_income occ_grade occ_wage 
			occ_ln_income_arith occ_ln_income_geo occ_ln_wage_arith occ_ln_wage_geo
			occ_var_ln_income occ_var_ln_wage, i(occ_code) j(group);
	
	sort group occ_code ; 
	
	save `i'_extract_cohort_mid_final, replace ; 


/*************************************************************************************************/
/****************************************** Cohort Old *******************************************/
/*************************************************************************************************/
								

	/* Step 1:  Compute the averages for all occupations - include the part time people as 0.5 in their paid occupation (using occ_code)*/ ; 									

	#delimit ;  
	cd c:\ErikMain\discrimination_growth;
	use `i'_extract_composite_main, clear ; 	

	/*   A. Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 
			  
	keep if cohort == 3 ;  /*  This is the line that changes across the cohort files */
	
	
	/*  B. Define variables for counts, earnings, education and wages by occupation*group */
	
	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;	
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt], by(occ_code) ;																								
								
	drop if occ_code == 0;
								
	save `i'_extract_cohort_old_step1, replace ;								
										

	/* Step 2:  Compute the averages for all occupations - include the part time people as 0.5 in home sector (using occ_code2)*/ ; 									
									
	#delimit ;
	
	use `i'_extract_composite_main, clear ; 	

	/*   A.  Define Cohorts */ 

	gen     cohort = 1 if age >= 25 & age <= 34 ;
	replace cohort = 2 if age >= 35 & age <= 44 ;
	replace cohort = 3 if age >= 45 & age <= 54 ; 

	keep if cohort == 3 ;  	/*  This is the line that changes across the cohort files */

	
	/*  B.  Define variables for counts, earnings, education and wages by occupation*group */

	gen working_full = emp_full_lastyear == 1 & emp_full_adj == 1 & incwage_full ~=. ;	
	
	gen count_num = 1;

	gen count_num_adj = count_num * person_adj ; 

	gen white_man_adj 	= white_man * person_adj ;
	gen white_woman_adj = white_woman * person_adj ;
	gen black_man_adj 	= black_man * person_adj ;
	gen black_woman_adj = black_woman * person_adj ;
	
	gen highgrade_adj 	= highgrade * person_adj ;

	gen income_full_all = incwage_full * working_full ;
	gen income_full_wm = white_man * incwage_full * working_full;
	gen income_full_ww = white_woman * incwage_full * working_full;
	gen income_full_bm = black_man * incwage_full * working_full;
	gen income_full_bw = black_woman * incwage_full * working_full;

	gen highgrade_adj_wm = highgrade * white_man ;
	gen highgrade_adj_ww = highgrade * white_woman ;
	gen highgrade_adj_bm = highgrade * black_man ;
	gen highgrade_adj_bw = highgrade * black_woman ;

	gen wage_wm = white_man * wage * working_full;
	gen wage_ww = white_woman * wage * working_full;
	gen wage_bm = black_man * wage * working_full;
	gen wage_bw = black_woman * wage * working_full;

	gen working_full_wm = white_man * working_full ;
	gen working_full_ww = white_woman * working_full ;
	gen working_full_bm = black_man * working_full ;
	gen working_full_bw = black_woman * working_full ;
	
	gen ln_income_full_all = ln(income_full_all);
	gen ln_income_full_wm = ln(income_full_wm) ;
	gen ln_income_full_ww = ln(income_full_ww) ;
	gen ln_income_full_bm = ln(income_full_bm) ;
	gen ln_income_full_bw = ln(income_full_bw) ;
	
	gen ln_wage_wm = ln(wage_wm) ;
	gen ln_wage_ww = ln(wage_ww) ;
	gen ln_wage_bm = ln(wage_bm) ;
	gen ln_wage_bw = ln(wage_bw) ;
	
	gen ln_income_full_all_sd = ln(income_full_all) ;
	gen ln_income_full_wm_sd = ln(income_full_wm) ;
	gen ln_income_full_ww_sd = ln(income_full_ww) ;
	gen ln_income_full_bm_sd = ln(income_full_bm) ;
	gen ln_income_full_bw_sd = ln(income_full_bw) ;
	
	gen ln_wage_wm_sd = ln(wage_wm) ;
	gen ln_wage_ww_sd = ln(wage_ww) ;
	gen ln_wage_bm_sd = ln(wage_bm) ;
	gen ln_wage_bw_sd = ln(wage_bw) ;

	format white_man_adj %12.1f ;
	format black_man_adj %12.1f ;
	format black_woman_adj %12.1f ;
	format white_woman_adj %12.1f ;
	format count_num %18.1f ;
	format count_num_adj %18.1f ; 
	format highgrade_adj %12.1f ; 
	format incwage_full	%18.1f ;					
	format income_full_all %18.1f ;
	format income_full_wm %18.1f;
	format income_full_bm %18.1f;
	format income_full_ww %18.1f;
	format income_full_bw %18.1f;
	format highgrade_adj_wm %12.1f;
	format highgrade_adj_ww %12.1f;
	format highgrade_adj_bm %12.1f;
	format highgrade_adj_bw %12.1f;

	format wage_wm %18.1f;
	format wage_bm %18.1f;
	format wage_ww %18.1f;
	format wage_bw %18.1f;

	format working_full % 12.1f;
	format working_full_wm %12.1f;
	format working_full_ww %12.1f;
	format working_full_bm %12.1f;
	format working_full_bw %12.1f;

	format ln_income_full_all %18.1f;
	format ln_income_full_wm %18.1f;
	format ln_income_full_bm %18.1f;
	format ln_income_full_ww %18.1f;
	format ln_income_full_bw %18.1f;
	
	format ln_wage_wm %18.1f;
	format ln_wage_bm %18.1f;
	format ln_wage_ww %18.1f;
	format ln_wage_bw %18.1f;	

	format ln_income_full_all_sd %18.1f;
	format ln_income_full_wm_sd %18.1f;
	format ln_income_full_bm_sd %18.1f;
	format ln_income_full_ww_sd %18.1f;
	format ln_income_full_bw_sd %18.1f;
	
	format ln_wage_wm_sd %18.1f;
	format ln_wage_bm_sd %18.1f;
	format ln_wage_ww_sd %18.1f;
	format ln_wage_bw_sd %18.1f;
	
	collapse 	(count) count_num  
				(sum)   incwage_full ln_incwage_full count_num_adj black_man_adj white_man_adj white_woman_adj black_woman_adj highgrade_adj
						income_full_all income_full_wm income_full_bm income_full_ww income_full_bw 
						highgrade_adj_wm highgrade_adj_ww highgrade_adj_bm highgrade_adj_bw 
						wage_wm wage_bm wage_ww wage_bw  working_full working_full_wm working_full_ww working_full_bm working_full_bw 
						ln_income_full_all ln_income_full_wm ln_income_full_ww ln_income_full_bm ln_income_full_bw 
						ln_wage_wm ln_wage_ww ln_wage_bm ln_wage_bw 
				(sd)	ln_income_full_all_sd ln_income_full_wm_sd ln_income_full_ww_sd ln_income_full_bm_sd ln_income_full_bw_sd 
						ln_wage_wm_sd ln_wage_ww_sd ln_wage_bm_sd ln_wage_bw_sd 
						[fw=perwt],  by(occ_code2) ;						
	
	keep if occ_code2 == 0;   /*  This is key difference from step 1 file - we only keep the home sector here */
	rename occ_code2 occ_code ; 		 

	save `i'_extract_cohort_old_step2, replace ;	

	
	/* Step 3:  Merge step1 and step2 files.  Then extract key variables.  Reshape data for easier use. */ ; 									

	use 			`i'_extract_cohort_old_step1; 
	append using	`i'_extract_cohort_old_step2;

	/* Define variables */
	
	*  Note "0" refers to all groups
	*  Note "1"  refers to white men
	*  Note "2"  refers to white women	
	*  Note "3"  refers to white black men	
	*  Note "4"  refers to white black women	
	*  These are the code that chad used for running our analysis file
	
	*  Note "occ_income" is average earnings for full time workers in a given occupation-group
	*  Note "occ_grade" is average years of schooling for workers in a given occupation-group
	*  Note "wage" is the average hourly wage for full time workers in a given occupation-group
	*  Note "num" is the number of people in a given occupation-group. 
	
	*  Note:  The program is set to make both geometric and arithmatic means of variables ;

	gen year = `i' ;
	
	gen 	cohort = . ; 
	replace cohort = 3 if `i' == 2012;
	replace cohort = 4 if `i' == 2000;
	replace cohort = 5 if `i' == 1990;
	replace cohort = 6 if `i' == 1980;
	replace cohort = 7 if `i' == 1970;
	replace cohort = 8 if `i' == 1960;
	
	gen num0 	= count_num_adj;
	gen num1  	= white_man_adj;
	gen num2	= white_woman_adj;
	gen num3  	= black_man_adj;
	gen num4	= black_woman_adj;
	
	gen occ_income0 = income_full_all/working_full;
	gen occ_income1 = income_full_wm/working_full_wm;
	gen occ_income2 = income_full_ww/working_full_ww;
	gen occ_income3 = income_full_bm/working_full_bm;
	gen occ_income4 = income_full_bw/working_full_bw;

	gen occ_ln_income_arith0 = ln(occ_income0);
	gen occ_ln_income_arith1 = ln(occ_income1);
	gen occ_ln_income_arith2 = ln(occ_income2);
	gen occ_ln_income_arith3 = ln(occ_income3);
	gen occ_ln_income_arith4 = ln(occ_income4);
	
	gen occ_ln_income_geo0 = ln_income_full_all/working_full;
	gen occ_ln_income_geo1 = ln_income_full_wm/working_full_wm;
	gen occ_ln_income_geo2 = ln_income_full_ww/working_full_ww;
	gen occ_ln_income_geo3 = ln_income_full_bm/working_full_bm;
	gen occ_ln_income_geo4 = ln_income_full_bw/working_full_bw;
	
	gen occ_grade0 = highgrade_adj/count_num_adj;
	gen occ_grade1 = highgrade_adj_wm/white_man_adj;
	gen occ_grade2 = highgrade_adj_ww/white_woman_adj;
	gen occ_grade3 = highgrade_adj_bm/black_man_adj;
	gen occ_grade4 = highgrade_adj_bw/black_woman_adj;

	gen occ_wage0 = . ;
	gen occ_wage1 = wage_wm/working_full_wm;
	gen occ_wage2 = wage_ww/working_full_ww;
	gen occ_wage3 = wage_bm/working_full_bm;
	gen occ_wage4 = wage_bw/working_full_bw;

	gen occ_ln_wage_arith0 = ln(occ_wage0);
	gen occ_ln_wage_arith1 = ln(occ_wage1);
	gen occ_ln_wage_arith2 = ln(occ_wage2);
	gen occ_ln_wage_arith3 = ln(occ_wage3);
	gen occ_ln_wage_arith4 = ln(occ_wage4);
	
	gen occ_ln_wage_geo0 = . ;
	gen occ_ln_wage_geo1 = ln_wage_wm/working_full_wm;
	gen occ_ln_wage_geo2 = ln_wage_ww/working_full_ww;
	gen occ_ln_wage_geo3 = ln_wage_bm/working_full_bm;
	gen occ_ln_wage_geo4 = ln_wage_bw/working_full_bw;
	
	gen occ_var_ln_income0 = (ln_income_full_all_sd)^2;
	gen occ_var_ln_income1 = (ln_income_full_wm_sd)^2;
	gen occ_var_ln_income2 = (ln_income_full_ww_sd)^2;
	gen occ_var_ln_income3 = (ln_income_full_bm_sd)^2;
	gen occ_var_ln_income4 = (ln_income_full_bw_sd)^2;	

	gen occ_var_ln_wage0 = . ;
	gen occ_var_ln_wage1 = (ln_wage_wm_sd)^2;
	gen occ_var_ln_wage2 = (ln_wage_ww_sd)^2;
	gen occ_var_ln_wage3 = (ln_wage_bm_sd)^2;
	gen occ_var_ln_wage4 = (ln_wage_bw_sd)^2;
	
	
	/*  Keep the main variables for our analysis file*/
	
	keep 	year cohort occ_code num0-num4 occ_income0-occ_income4 occ_grade0-occ_grade4  occ_wage0-occ_wage4 
			occ_ln_income_arith0-occ_ln_income_arith4 occ_ln_income_geo0-occ_ln_income_geo4
			occ_ln_wage_arith0-occ_ln_wage_arith4 occ_ln_wage_geo0-occ_ln_wage_geo4
			occ_var_ln_income0-occ_var_ln_income4  occ_var_ln_wage0-occ_var_ln_wage4 ;

	reshape long num occ_income occ_grade occ_wage 
			occ_ln_income_arith occ_ln_income_geo occ_ln_wage_arith occ_ln_wage_geo
			occ_var_ln_income occ_var_ln_wage, i(occ_code) j(group);
	
	sort group occ_code ; 
	
	save `i'_extract_cohort_old_final, replace ;

	
/*****************************************************************************************************/
/**************************  Combine the Cohort Files For Year ***************************************/
/*****************************************************************************************************/
	
	#delimit ;
	
	use 		  `i'_extract_cohort_all_final, clear ;
	append using  `i'_extract_cohort_young_final	;	
	append using  `i'_extract_cohort_mid_final	;
	append using  `i'_extract_cohort_old_final ; 
	
	sort cohort group occ_code ; 
	
	save `i'_extract_cohort_all_final, replace ;
	
	};
	

	
/*****************************************************************************************************/
/**************************  Combine the Cohort Files For Year ***************************************/
/*****************************************************************************************************/

	* Note:  The above part of the program is one big loop
	* Note:  After running that loop, need to manually combine the files ;
	
	#delimit ;
	
	use 		  2012_extract_cohort_all_final, clear ;
	append using  2000_extract_cohort_all_final	;
	append using  1990_extract_cohort_all_final ; 
	append using  1980_extract_cohort_all_final ; 
	append using  1970_extract_cohort_all_final ; 
	append using  1960_extract_cohort_all_final ; 
	
	gen occ_code_num = occ_code + 0 ;
	
	outsheet year group cohort occ_code num occ_income occ_grade occ_wage occ_ln_income_arith occ_ln_income_geo occ_ln_wage_arith occ_ln_wage_geo occ_var_ln_income occ_var_ln_wage using chad_output_file.csv, comma nolabel replace ;
	
	save 1960_2012_extract_cohort_allyears_final, replace ;	
				