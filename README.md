
# Code supplement for "Optimal post-selection inference" by Woody, Padilla & Scott

This code replicates the simulation and real data examples from the
paper ["Optimal post-selection inference"][saFAB paper] by Woody, Padilla,
and Scott (2020).

## Install saFAB R package

First, install the R package [`saFAB`][saFAB package] with the
following command:

```r
## install.packages("devtools")
devtools::install_github("spencerwoody/saFAB")
```

## Contents

There are three main applications in the paper. Here is a list of
these applications, along with the corresponding scripts.

1. The toy example Section 2.3:
    - `toy-length-comparison.R` compares the lengths of UMAU and saFAB
      intervals
    - `toy-coverage-plot.R` compares thes coverage properties of saFAB
      intervals to Bayesian selection-adjusted credible intervals and
      non-selection-adjusted UMAU intervals.
2. The simulation examples in Sections 4.1--4.3:
    - `sims-create-intervals.R` is a wrapper function for calculating
	intervals for a set of selected observed y's, given a spending
	function.  This makes the other scripts much more compact.
    - The other files `sims-*.R` create the simulation examples,
      comparing the efficiency and coverage of intervals created by
      three versions of saFAB: (i) oracle, where the prior is known
      exactly, (ii) parametric empircal Bayes, where we assume the
      form of the prior and estimate hyperparameters through maximum
      marginal likelihood, and (iii) nonparametric empirical Bayes.
      These are all compared against the UMAU intervals. These scripts
      use the parameters defined in the corresponding scripts
      `setpars-*.R`.  
	  There are five simulation settings:
        + The prior is known exactly, with selection joint on (y,
          theta) (`sims-nonpar-joint.R`)
        + The prior is known, with selection conditional on theta
          (`sims-nonpar-conditional.R`)
        + The prior for the paramatric saFAB method is misspecified,
          with the true prior having a bimodal density (joint
          selection only) (`sims-misspec-bimodal.R`)
        + The prior for the paramatric saFAB method is misspecified,
          with the true prior having a skewed density (joint selection
          only) (`sims-misspec-skew.R`)
        + The selection rule is data-dependent through the
          Benjamini-Hochberg procedure (`sims-bh-joint.R`)
3. The neural synchrony data analysis in Section 4.4:
    - `cortex-analysis.R` compares the lengths of saFAB intervals
      versus UMAU intervals for the neural synchrony data.

<!-- Links -->

[saFAB paper]: https://arxiv.org/abs/1810.11042
[saFAB package]: https://github.com/spencerwoody/safab


