
# cQL

`cQL` implements the proposed clustered Q-learning procedure with the
M-out-of-N cluster bootstrap for making inference on tailoring variables
in a two-stage clustered SMART (cSMART).

The package is designed for one completed cSMART dataset at a time.
Given a user-supplied dataset and user-specified stage-1 and stage-2
Q-functions, `cQL()`:

1.  fits the stage-2 Q-function;
2.  estimates the degree of non-regularity;
3.  selects the stage-1 resample size `M`;
4.  constructs the stage-1 pseudo-outcome; and
5.  returns bootstrap-based inference for both stages.

The output includes:

- regression coefficient estimates;
- 95% confidence intervals;
- bootstrap standard errors;
- bootstrap p-values;
- significance stars; and
- a stage summary showing `N`, `N_rand`, and `M` for each stage.

## Installation

For local development from the `clusterQ/cQL` folder:

``` r
devtools::install(".")
```

Then load the package:

``` r
library(cQL)
```

## Data format required by `cQL()`

Before using the package, prepare the cSMART data in the following
format.

1.  The dataset must be at the individual level, so each row is one
    individual.
2.  The dataset must contain a cluster id column.
3.  The dataset must contain one stage-1 treatment column and one
    stage-2 treatment column.
4.  If the cSMART is a Design II or III trial with limited second-stage
    re-randomization, then the stage-2 treatment must be `NA` for
    clusters that were not re-randomized at stage 2.
5.  Cluster-level variables such as the treatments and candidate
    tailoring variables should be repeated across all individuals within
    the same cluster.

In other words:

- for Design I cSMART data, the stage-2 treatment column is fully
  observed;
- for Design II or III cSMART data, the stage-2 treatment column
  contains `NA` for non-re-randomized clusters.

The package can work with any binary coding of the two treatment
columns. If the treatments are not already coded as `-1` and `1`,
`cQL()` internally recodes them and stores the mapping in the fitted
object.

## Recommended workflow

The recommended workflow for your own data is:

1.  prepare the data in the required format;
2.  specify the stage-2 and stage-1 Q-functions;
3.  fit `cQL()`;
4.  inspect `fit$stage_summary` to see how many clusters are used at
    each stage and what `M` was selected; and
5.  inspect `fit$stage2` and `fit$stage1` to interpret the
    stage-specific tailoring effects.

## Example 1: Design I cSMART with full second-stage re-randomization

Below we generate a toy Design I cSMART dataset. In this design, all
clusters are re-randomized at stage 2, so the stage-2 treatment column
has no missing values.

### Step 1: generate a full re-randomization dataset

``` r
design1_data <- simulate_csmart_data(
  n_clusters = 40,
  cluster_size = 12,
  rerandomization = "full",
  seed = 101
)

head(design1_data)
#>   cluster_id patient_id X1 A1 X2 A2 response_status rerandomized          Y
#> 1          1          1 -1  1 -1  1              NA            1  0.3080466
#> 2          1          2 -1  1 -1  1              NA            1  0.2029207
#> 3          1          3 -1  1 -1  1              NA            1 -0.9791916
#> 4          1          4 -1  1 -1  1              NA            1  0.4294969
#> 5          1          5 -1  1 -1  1              NA            1  1.3481200
#> 6          1          6 -1  1 -1  1              NA            1  0.2041069
table(is.na(design1_data$A2))
#> 
#> FALSE 
#>   480
```

The table above should show that `A2` is never missing, which is what we
expect for Design I.

### Step 2: specify the Q-functions

``` r
s2_formula <- Y ~ X1 * A1 + A1 * A2 + A2:X2
s1_formula <- Y1 ~ X1 * A1
```

In this example:

- the stage-2 candidate tailoring variables are `A1` and `X2`, because
  they interact with `A2`;
- the stage-1 candidate tailoring variable is `X1`, because it interacts
  with `A1`.

### Step 3: fit the clustered Q-learning analysis

``` r
fit_design1 <- cQL(
  data = design1_data,
  stage2_formula = s2_formula,
  stage1_formula = s1_formula,
  cluster = "cluster_id",
  stage1_treat = "A1",
  stage2_treat = "A2",
  stage2_tailoring_vars = c("A1", "X2"),
  n_boot = 150,
  seed = 1,
  verbose = FALSE
)
```

### Step 4: inspect how the algorithm works internally

``` r
fit_design1$stage_summary
#>     stage  N N_rand  M                                        bootstrap
#> 1 Stage 2 40     40 40 Cluster bootstrap on stage-2 randomized clusters
#> 2 Stage 1 40     40 39                     M-out-of-N cluster bootstrap
```

How to read this output:

- `N` is the total number of clusters in the dataset.
- `N_rand` is the number of clusters randomized at that stage.
- `M` is the bootstrap resample size used at that stage.

For Design I:

- at stage 2, all clusters are randomized, so `N_rand = N`, and stage 2
  uses a cluster bootstrap over those randomized clusters, so
  `M = N_rand`;
- at stage 1, all clusters are also randomized, but the manuscript’s
  M-out-of-N rule may choose `M < N` when non-regularity is present.

### Step 5: inspect the stage-2 inference

``` r
fit_design1$stage2
#>          term  estimate   conf.low conf.high  std.error    p.value significance
#> 1 (Intercept) 0.4901491 0.33786123 0.5965443 0.06855097 0.01324503            *
#> 2          X1 0.5631818 0.40286217 0.6943769 0.07309110 0.01324503            *
#> 3          A1 0.3193648 0.18084590 0.4739967 0.06928196 0.01324503            *
#> 4          A2 0.4062106 0.28155217 0.5113478 0.06493599 0.01324503            *
#> 5       X1:A1 0.2638915 0.11258743 0.3863879 0.07140048 0.01324503            *
#> 6       A1:A2 0.2509546 0.12880095 0.3824527 0.06372556 0.01324503            *
#> 7       A2:X2 0.1649829 0.03323396 0.2893730 0.06508515 0.03973510            *
```

To interpret the stage-2 table, focus especially on the terms involving
`A2`. If an interaction involving `A2` has a confidence interval that
excludes zero and a small p-value, that suggests the corresponding
variable may be useful as a stage-2 tailoring variable.

### Step 6: inspect the stage-1 inference

``` r
fit_design1$stage1
#>          term  estimate  conf.low conf.high  std.error    p.value significance
#> 1 (Intercept) 0.8953649 0.7289788 1.0271973 0.07851363 0.01324503            *
#> 2          X1 0.5455391 0.3849705 0.6552291 0.07228551 0.01324503            *
#> 3          A1 0.5948832 0.4592119 0.7633633 0.07888810 0.01324503            *
#> 4       X1:A1 0.2815342 0.1199928 0.4338442 0.08139256 0.01324503            *
```

The stage-1 table is built after constructing the stage-1 pseudo-outcome
from the fitted stage-2 model. Terms involving `A1` are the stage-1
tailoring effects of interest.

## Example 2: Design II or III cSMART with limited second-stage re-randomization

Now we generate a toy cSMART dataset in which only a subset of clusters
is re-randomized at stage 2. This mimics the data format required for
Design II or III cSMARTs.

### Step 1: generate a partial re-randomization dataset

``` r
design23_data <- simulate_csmart_data(
  n_clusters = 48,
  cluster_size = 12,
  rerandomization = "responder",
  p_rerand = 0.65,
  seed = 2026
)

head(design23_data)
#>   cluster_id patient_id X1 A1 X2 A2 response_status rerandomized            Y
#> 1          1          1 -1  1 -1  1               1            1  1.603547763
#> 2          1          2 -1  1 -1  1               1            1  1.139607586
#> 3          1          3 -1  1 -1  1               1            1  1.212212930
#> 4          1          4 -1  1 -1  1               1            1  1.748961883
#> 5          1          5 -1  1 -1  1               1            1 -0.006233853
#> 6          1          6 -1  1 -1  1               1            1  0.078092873
table(is.na(design23_data$A2))
#> 
#> FALSE  TRUE 
#>   192   384
```

Here the stage-2 treatment column contains `NA` for clusters that were
not re-randomized at stage 2. This is the key formatting rule for Design
II or III data.

The package does not need separate arguments telling it whether the
re-randomized clusters were responders or non-responders. The crucial
input is simply:

- one row per individual; and
- `NA` in the stage-2 treatment column for clusters not re-randomized.

### Step 2: use the same Q-functions

``` r
s2_formula
#> Y ~ X1 * A1 + A1 * A2 + A2:X2
s1_formula
#> Y1 ~ X1 * A1
```

### Step 3: fit the clustered Q-learning analysis

``` r
fit_design23 <- cQL(
  data = design23_data,
  stage2_formula = s2_formula,
  stage1_formula = s1_formula,
  cluster = "cluster_id",
  stage1_treat = "A1",
  stage2_treat = "A2",
  stage2_tailoring_vars = c("A1", "X2"),
  n_boot = 150,
  seed = 11,
  verbose = FALSE
)
```

### Step 4: inspect the stage summary

``` r
fit_design23$stage_summary
#>     stage  N N_rand  M                                        bootstrap
#> 1 Stage 2 48     16 16 Cluster bootstrap on stage-2 randomized clusters
#> 2 Stage 1 48     48 46                     M-out-of-N cluster bootstrap
```

This output is often the easiest way to understand what the algorithm is
doing internally.

For a partial re-randomization design:

- the stage-2 row should have `N_rand < N`, because only some clusters
  are randomized at stage 2;
- the stage-2 bootstrap resamples those stage-2 randomized clusters, so
  the stage-2 row has `M = N_rand`;
- the stage-1 row reports the manuscript-selected stage-1 `M` for the
  M-out-of-N bootstrap.

### Step 5: inspect the stage-2 inference

``` r
fit_design23$stage2
#>          term   estimate    conf.low conf.high std.error    p.value
#> 1 (Intercept) 0.82219945  0.41575672 1.1887897 0.2309933 0.01324503
#> 2          X1 0.48175821  0.07423095 0.7917144 0.2018352 0.02649007
#> 3          A1 0.13125283 -0.22044677 0.4586660 0.2026440 0.76821192
#> 4          A2 0.25409355 -0.16644961 0.6783350 0.2306893 0.15894040
#> 5       X1:A1 0.09933682 -0.26642517 0.5806213 0.2183360 0.56953642
#> 6       A1:A2 0.30850107 -0.04729216 0.7043165 0.2110339 0.13245033
#> 7       A2:X2 0.53327211  0.11894888 1.2115914 0.2467119 0.01324503
#>   significance
#> 1            *
#> 2            *
#> 3             
#> 4             
#> 5             
#> 6             
#> 7            *
```

The interpretation is the same as before, but now the stage-2 regression
is fit only to the clusters that were actually re-randomized at stage 2.

### Step 6: inspect the stage-1 inference

``` r
fit_design23$stage1
#>          term  estimate  conf.low conf.high std.error    p.value significance
#> 1 (Intercept) 1.0270226 0.6847385 1.1649680 0.1289050 0.01324503            *
#> 2          X1 0.5926596 0.4492606 0.8767815 0.1150225 0.01324503            *
#> 3          A1 0.6173133 0.3963612 0.7854266 0.1018887 0.01324503            *
#> 4       X1:A1 0.2962095 0.1051471 0.6409662 0.1369957 0.01324503            *
```

For limited second-stage re-randomization, the manuscript’s stage-1
pseudo-outcome rule is used:

- for clusters re-randomized at stage 2, the pseudo-outcome uses the
  fitted stage-2 Q-function;
- for clusters not re-randomized at stage 2, the pseudo-outcome equals
  the observed outcome.

## Applying `cQL()` to your own data

After your data are prepared, the analysis for your own cSMART dataset
will look like this:

``` r
my_fit <- cQL(
  data = my_csmart_data,
  stage2_formula = Y ~ X1 * A1 + A1 * A2 + A2:X2,
  stage1_formula = Y1 ~ X1 * A1,
  cluster = "cluster_id",
  stage1_treat = "A1",
  stage2_treat = "A2",
  stage2_tailoring_vars = c("A1", "X2"),
  n_boot = 1000,
  fixed_xi = 0.025
)

my_fit$stage_summary
my_fit$stage2
my_fit$stage1
```

In practice:

- start with `my_fit$stage_summary` to see the cluster counts used by
  the algorithm;
- then inspect `my_fit$stage2` to evaluate candidate stage-2 tailoring
  variables;
- then inspect `my_fit$stage1` to evaluate candidate stage-1 tailoring
  variables.

## Notes

- The current package targets the two-stage cSMART setting developed in
  the manuscript.
- The stage-2 formula must include the main effect of the stage-2
  treatment.
- Every interaction involving the stage-2 treatment must be a two-way
  interaction with one of the variables listed in
  `stage2_tailoring_vars`.
- In the printed stage summary, the stage-2 row has `M = N_rand` because
  stage 2 uses the cluster bootstrap over the stage-2 randomized
  clusters, whereas the stage-1 row uses the manuscript-selected
  M-out-of-N resample size.
