/* Now with matrix notation and PROC IML */
PROC IML;
	/* create the design matrix */	
	X={
		1 135 3,
		1 120 4,
		1 100 3,
		1 105 2,
		1 130 4,
		1 125 5,
		1 125 2,
		1 105 3,
		1 120 5,
		1  90 4,
		1 120 2,
		1  95 3,
		1 120 3,
		1 150 4,
		1 160 3,
		1 125 3};
	/* create the outcome vector */
	Y = {89,90,83,77,92,98,82,85,96,95,80,79,86,97,92,88};

	/* X'X, X'Y, Y'Y, (X'X)^{-1} */
	/* RECALL: SAS does not consider capitalization */
	xpx = x`*x;
	xpy = x`*y;
	ypy = y`*y;
	xpxi = INV(xpx); 

	/* Calculate regression components, sums of squares */
	betahat = xpxi*xpy;
	ssm=betahat`*x`*y-SUM(Y)**2/NROW(y); 
	sse=(y-x*betahat)`*(y-x*betahat);
	sst=ssm+sse;
	sighat2 = sse/(NROW(y)-NROW(betahat));
	varb = xpxi*sighat2; 

	/* Calculate beta_1, variance, t-stat, p-val */
	cb1 = {0 1 0};
	b1 = cb1*betahat; 
	varb1=cb1*xpxi*cb1`*sighat2;
	seb1=SQRT(varb1);
	tb1 = cb1*betahat/seb1;
	pb1 = 2*PROBT(-ABS(tb1),(NROW(y)-NROW(betahat)));

	/* Calculate beta_2, variance, t-stat, p-val */
	cb2 = {0 0 1};
	b2 = cb2*betahat;
	varb2=cb2*xpxi*cb2`*sighat2;
	seb2=SQRT(varb2);
	tb2 = cb2*betahat/seb2;
	pb2 = 2*PROBT(-ABS(tb2),(NROW(y)-NROW(betahat)));

	/* Print various summary objects to review */
	PRINT y x;
	PRINT xpx xpy ypy;
	PRINT xpxi;
	PRINT betahat;
	PRINT ssm;
	PRINT sse;
	PRINT sst;
	PRINT sighat2;
	PRINT varb;
	PRINT b1 seb1 tb1 pb1;
	PRINT b2 seb2 tb2 pb2; 
	/* NOTE: You can change the field width for printing numeric values with the RESET command */
	RESET fw=12;
	PRINT xpxi;
	RESET fw=4;
	PRINT b2 seb2 tb2 pb2;
QUIT;

/* Check our matrix-derived answers with PROC REG */
DATA birth;
INPUT sbp weight age;
LABEL 	sbp 	 = "SBP (mmHg)"
		weight = "Weight (oz)"
		age    = "Age (days)";
DATALINES;
89 135 3
90 120 4
83 100 3
77 105 2
92 130 4
98 125 5
82 125 2
85 105 3
96 120 5
95  90 4 
80 120 2
79  95 3
86 120 3
97 150 4
92 160 3
88 125 3
;
RUN;

PROC REG DATA=birth;
	MODEL sbp = weight age / clb covb xpx; /* xpx prints matrix crossproducts */
RUN;
