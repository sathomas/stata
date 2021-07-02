// This do-file is an example script demonstrating various approaches for
// finding parameters of arbitrary models in Stata.
//
// In general the following code tries to stick with common Stata coding
// conventions. The major exception is abbreviations; few abbreviations
// are found in the lines below as they can make code harder to read.

version 17          // Stata version used to create/execute this script
set more off        // disable pause of output display every full screen
clear all           // ensure clean slate
capture log close   // close any pending logs
set seed 123456789  // consistent random number generation

**# ########## Observations ##########

// Data is based on an example is from P. M. Berthouex and L. C. Brown,
// Statistics for Environmental Engineers, CRC Press, 2002.
//
// The model is Monod growth:
//
// y = θ₁ t / (θ₂ + t) + ε, ε ~ N(0, Iσ²)
//
// x: 28    55    83    110   138   225   375
// y: 0.053 0.060 0.112 0.105 0.099 0.122 0.125

quietly {
	input x y
		 28 0.053
		 55 0.060
		 83 0.112
		110 0.105
		138 0.099
		225 0.122
		375 0.125
	end
}

// Initial values for model parameters, obtained from external analysis

local theta1 = 0.15
local theta2 = 50
local sigma  = 0.01

**# ########## Non-Linear Least Squares ##########

// A simple approach for finding parameter values is non-linear least
// squares. It may be helpful to use "nl" as a first step even when
// using other approaches, as the "nl" results can provide good initial
// values and, e.g., estimates for the error variance.

// The "nl" command can be passed an expression directly, but for this
// example we'll use a program to demonstrate the more general approach.

// Stata requires programs used for non-linear least squares to have a
// name that starts with "nl". That prefix is then omitted when the
// program is referenced in the actual "nl" command.

capture program drop nlmonod
program define nlmonod
	syntax varlist(min=2 max=2) if, at(name)
	local y: word 1 of `varlist'
	local x: word 2 of `varlist'
	tempname theta1 theta2
	scalar `theta1' = `at'[1, 1]
	scalar `theta2' = `at'[1, 2]
	replace `y' = (`theta1' * `x') / (`theta2' + `x') `if'
end

nl monod @ y x, ///
    parameters(theta1 theta2) ///
	initial(theta1 `theta1' theta2 `theta2')
estimates store nl

local theta1_nl = _b[/theta1]
local theta2_nl = _b[/theta2]

// Now we have a reasonable estimate of the error variance. This can be
// helpful in other methods (e.g. maximum likelihood and bayes)
// σ² ≈ MSE = RSS / (n - p)
local sigma = sqrt(e(rss) / (_N - 2))

display "NL Estimates: θ₁ `theta1' θ₂ `theta2' σ `sigma' RSS " e(rss)

// One benefit of "nl" is the ability to use "predictnl" to generate
// new variables based on the parameters. Note, though, that we can
// only do this if the model can be written as an expression in addition
// or or instead of a program. That's the case here. As an option, we
// can also ask for confidence intervals.
quietly predictnl nl_pred = (_b[/theta1] * x) / (_b[/theta2] + x), ///
	ci(nl_low nl_high)

**# ########## Maximum likelihood (ml) estimation ##########

// Custom program to evaluate the log-likelihood function. We'll use the
// linear form for the maximum likelihood estimation, so this program is
// a "method-lf evaluator."
//
// For this problem we could get away with "mlexp" and use an expression
// instead of a program. As with non-linear least squares, though, we'll
// go with a full program for maximum generality.
//
// Method-lf evaluators are required to evaluate the observation-by-
// observation log likelihood ln(Lj), j = 1, …, N. The subscript j 
// indexes the observations.

capture program drop monod_lfeval
program define monod_lfeval
    // Arguments:
    //   lnfj: variable to be filled in with observation-by-observation
    //         values of ln(Lj)
    //   theta1x: variable containing evaluation of first equation
    //            θ1j = x1j b1
    //   theta2x: variable containing evaluation of second equation
    //            θ2j = x2j b2
    //   sigma: variable containing free parameter for standard deviation
    //          of errors
    //
    // Equations are defined such that
	//    theta1x = theta1 ⨉ x   [note multiplication, no constant term]
	//    theta2x = theta2 + x   [note addition, only a constant term]
	args lnfj theta1x theta2x sigma

    // Get access to the dependent variable (y). This isn't strictly
    // necessary (and probably degrades performance), but it makes the
    // expressions below easier to read.
	local y "$MH_y"
	
	// For convenience define a temporary variable to store the residuals.
	// Again, this probably degrades performance slightly but impreoves
	// readability.
	tempvar residuals
	
	// Calculate log-likelihood for each observation
	// ln(Lj) = -ln(2ᴨ)/2 - ln(σ) - (yj - θ₁ xj / (θ₂ + xj))² / 2σ²
	//
	// Note that $ML_samp isn't necessary for method-lf evaluators, but
	// it doesn't hurt and may improve performance. The performance
	// improvement is especially significant for other evaluators as it
	// allows the use of the "nopreserve" option with the "ml" command.
	quietly {
    	generate double `residuals' = `y' - `theta1x' / `theta2x' if $ML_samp == 1
		replace `lnfj' = -0.5 * ln(2 * _pi) - ln(`sigma') - ///
			0.5 * (`residuals' / `sigma')^2 if $ML_samp == 1
	}
end

// Calculate the maximum likelihood estimate using three "equations"
//    θ₁x
//    θ₂ + x
//    σ
// This is non-interactive mode, so the "maximize" option is specified. We
// also specify "nopreserve" as good practice even though it has no effect
// for linearform (lf) methods. In general, it tells Stata to assume that the
// evaluator is correctly using $ML_samp; without that assumption Stata takes
// precautions that degrade performance.
ml model linearform monod_lfeval ///
	(theta1x: y = x, noconstant) ///
	(theta2x:, offset(x)) ///
	(sigma:, freeparm) ///
	, maximize search(off) technique(nr) difficult nopreserve ///
	init(theta1x:x = `theta1' theta2x:_cons = `theta2' sigma = `sigma')
estimates store ml

// Save the parameter estimates
local theta1_ml = _b[theta1x:x]
local theta2_ml = _b[theta2x:_cons]
local sigma  = e(b)[1, colnumb(e(b), "sigma")]

// Update sigma to account for sample bias
quietly generate ml_pred = (`theta1_ml' * x) / (`theta2_ml' + x)
quietly generate ml_sqerr = (y - ml_pred) ^ 2
summarize ml_sqerr, meanonly
local s = sqrt(r(sum) / (_N - 2))
display "ML Estimates: θ₁ `theta1' θ₂ `theta2' σ `sigma' s `s'"

**# ########## Bayesian (MCMC) ##########

// Bayesian version using a log-likelihood evaluator. We can use a
// log-likelihood evaluator instead of a log-posterior evaluator because
// we're okay using built-in distributions as priors for the paraemters.

capture program drop monod_lleval
program define monod_lleval

	// We don't have any linear predictors, so the arguments are
	//   lnden: the name of a temporary scalar to be filled in with an overall
	//          log-likelihood value
	//   theta1: first parameter
	//   theta2: second parameter
	args lnden theta1 theta2 

	// We do need access to the unmodified (by a regression coefficient)
	// indepedent variable x; it's passed as an extra variable. In addition
	// the standard deviation is passed through, and we create a convenient
	// reference to the dependent variable y.
	local x "$MH_extravars"
	local y "$MH_y"
	local sigma "$MH_passthruopts"

	// Calculate squared error at each point; ensure double precision for
	// maximum accuracy. Technically we don't need $MH_touse since we always
	// use all observations, but it's safer in general to include it.
	tempvar sq_err
	generate double `sq_err' = (`y' - (`theta1' * `x' / (`theta2' + `x'))) ^ 2 ///
		if $MH_touse

	// Find sum of squared errors (i.e. SSE)
	summarize `sq_err', meanonly

	// If there was a problem, don't report an actual result but return a
	// missing value instead. This is the Stata convention. To detect problems
	// we ensure that the number of generated squared error values is the same
	// as the number of observations. (The number of squared error values is
	// reported in r(N) by the summarize command above.)
	if r(N) < $MH_n {
		scalar `lnden' = .
		exit
	}
	
	// No problems, so sum of squared errors is in r(sum) reported by the
	// summarize command.

	// Negative log likelihood ∝ sum of squared errors. Stata allows us to
	// ignore any constant terms, but it does warn that "some of the reported
	// statistics such as DIC and log marginal-likelihood may not be
	// applicable."
	scalar `lnden' = -1 / (2 * `sigma'^2) * r(sum)
end

// Bayes results include the MCMC posteriors and require significantly more
// storage space. Use a temporary file to hold them.
tempfile bayes_file

// Note that we're passing through the estimate for σ obtained from the
// maximum likelihood estimate, and we're assuming it is constant. It
// would definitely be better to adjust this value dynamically, but--at
// least so far--I haven't figured out a way to do that in Stata.
local mcmc_samples = 100000
bayesmh y, noconstant ///
	llevaluator( ///
		monod_lleval, ///
		parameters({theta1} {theta2}) ///
		extravars(x) ///
		passthruopts(`sigma') ///
	) ///
	initial({theta1} `theta1' {theta2} `theta2') ///
	prior({theta1}, flat) ///
	prior({theta2}, flat) ///
	mcmcsize(`mcmc_samples') burnin(25000) ///
	saving(`bayes_file')
estimates store bayes

// For a simple comparison, use the mean values from the MCMC posterior
local theta1_b = e(mean)[1, colnumb(e(mean), "theta1")]
local theta2_b = e(mean)[1, colnumb(e(mean), "theta2")]
quietly generate bayes_pred = (`theta1_b' * x) / (`theta2_b' + x)

display "Bayes Estimates: θ₁ `theta1' θ₂ `theta2'"

// For more extensive analysis load the posterior data (in a separate
// data frame). The resulting dataset will contain the parameter values
// realized in the MCMC chains. With a small bit of housekeeping, we can
// randomly sample from those posterior values.

tempname posterior
frame create `posterior'
frame `posterior': use `bayes_file'
quietly frame `posterior' {
	expand _frequency
	keep eq*
	rename eq0_p1 theta1
	rename eq0_p2 theta2
	generate double randu = runiform()
	isid randu
	sort randu
}

**# ########## Compare Results ##########

set scheme s1color

// One plot for each method separately to highlight that method's unique
// features, and then a combined plot to allow comparison of the results.

// Non-linear Least Squares - include confidence interval
twoway ///
	(rarea nl_low nl_high x, fcolor(gs10%25) lwidth(none)) ///
	(function y = ((`theta1_nl' * x) / (`theta2_nl' + x)),  range(x)) ///
	(scatter y x), ///
	legend(off) ///
	subtitle("Non-linear Least Squares", size(8pt) position(12) ring(0)) ///
	xtitle("") xlabel(none) ytitle("") ylabel(none) ///
	name(nl, replace) nodraw

// Maximum Likelihood - not really any options to addition
twoway ///
	(function y = ((`theta1_ml' * x) / (`theta2_ml' + x)),  range(x)) ///
	(scatter y x), ///
	legend(off) ///
	subtitle("Maximim Likelihood", size(8pt) position(12) ring(0)) ///
	xtitle("") xlabel(none) ytitle("") ylabel(none) ///
	name(ml, replace) nodraw

// Bayesian - show plausible trajectories by sampling from the posterior
forvalues sample = 1/500 {
	local i = runiformint(1, `mcmc_samples')
	local t1 = _frval(`posterior', theta1, `i')
	local t2 = _frval(`posterior', theta2, `i')
	local graph = "`graph' (function y = ((`t1' * x) / (`t2' + x)),"
	local graph = "`graph' range(x) lwidth(thin) lcolor(gs10%10))"
}

twoway ///
	`graph' ///
	(function y = ((`theta1_b' * x) / (`theta2_b' + x)),  range(x)) ///
	(scatter y x), ///
	legend(off) ///
	subtitle("Bayesian (MCMC)", size(8pt) position(12) ring(0)) ///
	xtitle("") xlabel(none) ytitle("") ylabel(none) ///
	name(bayes, replace) nodraw

// Combined results on single plot - using mean values for parameters
twoway ///
	(function y = ((`theta1_nl' * x) / (`theta2_nl' + x)),  range(x)) ///
	(function y = ((`theta1_ml' * x) / (`theta2_ml' + x)),  range(x)) ///
	(function y = ((`theta1_b' * x) / (`theta2_b' + x)),  range(x)) ///
	(scatter y x), ///
	subtitle("Mean Estimates", size(8pt) position(12) ring(0)) ///
	legend(order( ///
		4 "Observations" ///
		1 "Non-linear Least Squares" ///
		2 "Maximum Likelihood" ///
		3 "Bayesian (MCMC)" ///
	) ///
	rowgap(0) keygap(1) symxsize(9) size(6pt) region(margin(small)) ///
	cols(1) position(5) ring(0)) ///
	xtitle("") xlabel(none) ytitle("") ylabel(none) ///
	name(all, replace) nodraw
	
graph combine nl ml bayes all, cols(2) name(models, replace)
graph export "models.png", as(png) name(models) replace

