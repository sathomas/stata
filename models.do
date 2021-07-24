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

// Data is based on an example from OpenBUGS.
//
//     y ~ ɑ - β λˣ + N(0, Iσ)
//

/*

R-language nls fit

       Estimate Std. Error t value Pr(>|t|)
alpha   2.65807    0.06151   43.21  < 2e-16
beta    0.96352    0.06968   13.83  6.3e-13
lamda   0.87146    0.02460   35.42  < 2e-16


Stan Model:

    data {
        int<lower=0> N;
        array[N] real<lower=0> x;
        array[N] real<lower=0> y;
    }
    parameters {
        real alpha;
        real beta;
        real<lower=.5, upper= 1> lamda;
        real<lower=0> sigma;
    } 
    transformed parameters {
        array[N] real m;
        for (i in 1:N)
            m[i] = alpha - beta * pow(lamda, x[i]);
    }
    model {
        // priors
        alpha ~ normal(0.0, 1000);
        beta  ~ normal(0.0, 1000);
        lamda ~ uniform(0.5, 1);
        sigma ~ inv_gamma(0.01, 0.01);
    
        // likelihood
        y ~ normal(m, sigma);
    }


Stan Summary:

                 Mean     MCSE   StdDev     5%    50%   95%

alpha             2.7  3.5e-03    0.069    2.6    2.6   2.8 
beta             0.97  3.2e-03    0.072   0.85   0.97   1.1
lamda            0.86  1.6e-03    0.031   0.81   0.87  0.91
sigma           0.099  7.2e-04    0.015  0.078  0.097  0.12

*/


quietly {
    input x y
         1    1.8
         1.5  1.85
         1.5  1.87
         1.5  1.77
         2.5  2.02
         4    2.27
         5    2.15
         5    2.26
         7    2.47
         8    2.19
         8.5  2.26
         9    2.4
         9.5  2.39
         9.5  2.41 
        10    2.5
        12    2.32 
        12    2.32  
        13    2.43
        13    2.47 
        14.5  2.56
        15.5  2.65   
        15.5  2.47
        16.5  2.64
        17    2.56
        22.5  2.7
        29    2.72
        31.5  2.57     
    end
}

// Initial values for model parameters, obtained from external analysis

local alpha = 1
local beta  = 1
local lamda = 0.9

**# ########## Non-Linear Least Squares ##########

// A simple approach for finding parameter values is non-linear least
// squares. It may be helpful to use "nl" as a first step even when
// using other approaches, as the "nl" results can provide good initial
// values and, e.g., estimates for the error variance.

// The simplest way to use "nl" is with expressions.

nl (y = {alpha} - {beta} * {lamda} ^ x), ///
    initial(alpha `alpha' beta `beta' lamda `lamda') ///
    nolog

// Although not necessary for this example, more complicated models may
// be implemented as "programs." 

// Stata requires programs used for non-linear least squares to have a
// name that starts with "nl". That prefix is then omitted when the
// program is referenced in the actual "nl" command.

capture program drop nlgrowth
program define nlgrowth
    syntax varlist(min=2 max=2) if, at(name)
    local y: word 1 of `varlist'
    local x: word 2 of `varlist'
    tempname alpha beta lamda
    scalar `alpha' = `at'[1, 1]
    scalar `beta'  = `at'[1, 2]
    scalar `lamda' = `at'[1, 3]
    replace `y' = (`alpha' - `beta' * `lamda' ^ `x') `if'
end

nl growth @ y x, ///
    parameters(alpha beta lamda) ///
    initial(alpha `alpha' beta `beta' lamda `lamda') ///
    nolog

local alpha_nl = _b[/alpha]
local beta_nl  = _b[/beta]
local lamda_nl = _b[/lamda]

// Least squares provides reasonable estimate of the error variance.
// This estimate can be helpful in other methods (e.g. maximum likelihood
// and bayes)
//
// σ² ≈ MSE = RSS / (n - p)
local sigma = sqrt(e(rss) / (_N - 2))

// One benefit of "nl" is the ability to use "predictnl" to generate
// new variables based on the parameters. Note, though, that we can
// only do this if the model can be written as an expression in addition
// or or instead of a program. That's the case here. As an option, we
// can also ask for confidence intervals.
quietly predictnl pred_nl = (_b[/alpha] - _b[/beta] * _b[/lamda] ^ x), ///
    ci(low_nl high_nl)

**# ########## Bayesian (MCMC) ##########

// The simplest approach uses substitutable expressions.

bayesmh y = ({alpha} - {beta} * {lamda} ^ x), ///
    likelihood(normal({var})) ///
    prior({alpha beta}, normal(0, 1000)) ///
    prior({lamda}, uniform(0.5, 1)) ///
    prior({var}, igamma(0.1, 0.1)) ///
    initial({alpha} `alpha' {beta} `beta' {lamda} `lamda' {var} `sigma'^2)

// For more complicated models, we can use a custom "program evaluator."
// We can use a log-likelihood evaluator instead of a log-posterior
// evaluator because we're okay using built-in distributions as priors
// for the paraemters.

// Our model has no linear combinations for the independent variable
// x, so ideally we would just define a constraint that fixed the slope.
// Currently the "bayesmh" command doesn't allow constraints with program
// evaluators, though. To work around that and avoid any spurious warnings
// for variables that we don't use, we pass x as an extra variable instead.
capture program drop growth_lleval
program define growth_lleval
    // Arguments:
    //   lnden: the name of a temporary scalar to be filled in with an overall
    //          log-likelihood value
    //   alpha, beta, lamda, sigma: model parameters (not linear equations)
    args lnden alpha beta lamda sigma 

    // Make the code a little bit more readable with clearer references
    local x "$MH_extravars"
    local y "$MH_y"

    tempvar lnfj
    quietly generate double `lnfj' = ///
        lnnormalden(`y', (`alpha' - `beta' * `lamda' ^ `x'), `sigma') ///
        if $MH_touse
    
    summarize `lnfj', meanonly

    // If there was a problem, don't report an actual result but return a
    // missing value instead. This is the Stata convention. To detect problems
    // we ensure that the number of generated values is the same as the number
    // of observations.
    if r(N) < $MH_n {
        scalar `lnden' = .
        exit
    }

    scalar `lnden' = r(sum)
end

// Bayes results include the MCMC posteriors and require significantly more
// storage space. Use a temporary file to hold them.
tempfile bayes_file

bayesmh y, noconstant ///
    llevaluator( ///
        growth_lleval, ///
        extravars(x) ///
        parameters({alpha} {beta} {lamda} {sigma}) ///
    ) ///
    prior({alpha beta}, flat) ///
    prior({lamda}, uniform(0.5, 1)) ///
    prior({sigma}, igamma(0.01, 0.01)) ///
    initial({alpha} `alpha' {beta} `beta' {lamda} `lamda' {sigma} `sigma') ///
    saving(`bayes_file')

// For a simple comparison, use the mean values from the MCMC posterior
local alpha_bs = e(mean)[1, "alpha"]
local beta_bs  = e(mean)[1, "beta"]
local lamda_bs = e(mean)[1, "lamda"]
local sigma_bs = e(mean)[1, "sigma"]
local mcmc_samples = e(mcmcsize)

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
    rename eq0_p1 alpha
    rename eq0_p2 beta
    rename eq0_p3 lamda
    rename eq0_p4 sigma
    generate double randu = runiform()
    isid randu
    sort randu
}

**# ########## Maximum likelihood (ml) estimation ##########

// The simplest way to use maximum likelihood is also with a
// substitutable expression.

mlexp (lnnormalden(y, {alpha} - {beta} * {lamda} ^ x, {sigma})), ///
    from(alpha = `alpha' beta = `beta' lamda = `lamda' sigma = `sigma') ///
    difficult nolog

// The "predictnl" command also works after "mlexp"
quietly predictnl pred_ml = (_b[/alpha] - _b[/beta] * _b[/lamda] ^ x), ///
    ci(low_ml high_ml)

// More complex problems may benefit from a complete "program."  We'll use
// the linear form for the maximum likelihood estimation, so this program is
// a "method-lf evaluator."
//
// Method-lf evaluators are required to evaluate the observation-by-
// observation log likelihood ln(Lj), j = 1, …, N. The subscript j 
// indexes the observations.

// Because the model doesn't have any parameters that create linear
// functions of x (i.e. no scale or location) we constrain the linear
// coefficient to be 1. We'll also force the constant term to be zero
// when we execute the "ml" command
constraint define 1 [xb]:x = 1

capture program drop growth_lfeval
program define growth_lfeval
    // Arguments:
    //   lnfj: variable to be filled in with observation-by-observation
    //         values of ln(Lj)
    //   x: independent variable obtained via a linear equation where the
    //      slope is constrained to be 1 (above) and no constant is allowed
    //      (below)
    //   alpha, beta, lamda, sigma: model parameters (not linear equations)
    args lnfj x alpha beta lamda sigma

    // Get access to the dependent variable (y). This isn't strictly
    // necessary (and possibly degrades performance), but it makes the
    // expressions easier to read.
    local y $ML_y1

    quietly replace `lnfj' = ///
        lnnormalden(`y', (`alpha' - `beta' * `lamda' ^ `x'), `sigma') ///
        if $ML_samp == 1
end

// This is non-interactive mode, so the "maximize" option is specified. We
// also specify "nopreserve" as good practice even though it has no effect
// for linearform (lf) methods. In general, it tells Stata to assume that the
// evaluator is correctly using $ML_samp; without that assumption Stata takes
// precautions that degrade performance.
ml model linearform growth_lfeval ///
    (xb: y = x, noconstant) ///
    (alpha:, freeparm) ///
    (beta:,  freeparm) ///
    (lamda:, freeparm) ///
    (sigma:, freeparm) ///
    , constraint(1) ///
    init(xb:x = 1 alpha = `alpha' beta = `beta' lamda = `lamda' sigma = `sigma') ///
    maximize search(off) difficult nolog nopreserve

local alpha_ml = _b[/alpha]
local beta_ml  = _b[/beta]
local lamda_ml = _b[/lamda]
local sigma_ml = _b[/sigma]

**# ########## Compare Results ##########

set scheme s1color

// One plot for each method separately to highlight that method's unique
// features, and then a combined plot to allow comparison of the results.

// Non-linear Least Squares - include confidence interval
twoway ///
    (rarea low_nl high_nl x, fcolor(gs10%25) lwidth(none)) ///
    (function y = (`alpha_nl' - `beta_nl' * `lamda_nl' ^ x),  range(x)) ///
    (scatter y x), ///
    legend(off) ///
    title("Non-linear Least Squares", size(8pt) position(12) ring(0)) ///
    subtitle("w/ 95% CI", size(6pt) position(5) ring(0)) ///
    xtitle("") xlabel(none) ytitle("") ylabel(none) ///
    name(nl, replace) nodraw

// Maximum Likelihood - also with confidence interval
twoway ///
    (rarea low_ml high_ml x, fcolor(gs10%25) lwidth(none)) ///
    (function y = (`alpha_ml' - `beta_ml' * `lamda_ml' ^ x),  range(x)) ///
    (scatter y x), ///
    legend(off) ///
    title("Maximim Likelihood", size(8pt) position(12) ring(0)) ///
    subtitle("w/ 95% CI", size(6pt) position(5) ring(0)) ///
    xtitle("") xlabel(none) ytitle("") ylabel(none) ///
    name(ml, replace) nodraw

// Bayesian - show plausible trajectories by sampling from the posterior
forvalues sample = 1/500 {
    local i = runiformint(1, `mcmc_samples')
    local a = _frval(`posterior', alpha, `i')
    local b = _frval(`posterior', beta,  `i')
    local l = _frval(`posterior', lamda, `i')
    local graph = "`graph' (function y = (`a' - `b' * `l' ^ x),"
    local graph = "`graph' range(x) lwidth(thin) lcolor(gs12%10))"
}

twoway ///
    `graph' ///
    (function y = (`alpha_bs' - `beta_bs' * `lamda_bs' ^ x),  range(x)) ///
    (scatter y x), ///
    legend(off) ///
    title("Bayesian (MCMC)", size(8pt) position(12) ring(0)) ///
    subtitle("w/ 500 posterior samples", size(6pt) position(5) ring(0)) ///
    xtitle("") xlabel(none) ytitle("") ylabel(none) ///
    name(bayes, replace) nodraw

// Combined results on single plot - using mean values for parameters
twoway ///
    (function y = (`alpha_nl' - `beta_nl' * `lamda_nl' ^ x),  range(x)) ///
    (function y = (`alpha_ml' - `beta_ml' * `lamda_ml' ^ x),  range(x)) ///
    (function y = (`alpha_bs' - `beta_bs' * `lamda_bs' ^ x),  range(x)) ///
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
    
graph combine nl ml bayes all, cols(2) name(models, replace) ///
    title("Stata Model Fitting Approaches", size(8pt))
graph export "models.png", as(png) name(models) replace
