
# cQL

This R package `cQL` implements a novel clustered Q-learning method with
the $M$-out-of-$N$ cluster bootstrap for making inference on tailoring
variables in a two-stage clustered SMART (cSMART).

The package is designed for analyzing one completed cSMART dataset at a
time. Given a user-supplied dataset and user-specified stage-1 and
stage-2 Q-functions, the main function `cQL()` could:

1.  fits the stage-2 Q-function;
2.  estimates the degree of non-regularity;
3.  selects the stage-1 resample size `M`;
4.  constructs the stage-1 pseudo-outcome; and
5.  returns bootstrap-based inference for both stages.

The output includes:

- regression coefficient estimates;
- confidence intervals at the requested level (95% by default);
- bootstrap standard errors;
- bootstrap p-values;
- significance stars; and
- a stage summary showing `N`, `N_rand`, and `M` for each stage.

## Installation

To install this package from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("SelinaSong0412/cQL")
```

For local development:

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

1.  The dataset should be at the individual level, so each row is one
    individual.
2.  The dataset should contain a cluster id column.
3.  The dataset should contain one stage-1 treatment column and one
    stage-2 treatment column.
4.  The dataset should contain the final observed outcome column `Y`.
5.  If the dataset also contains an observed stage-1 intermediate
    outcome, that column should be named `Y1`.
6.  If the cSMART is a Design II or III trial with limited second-stage
    re-randomization, then the stage-2 treatment should be `NA` for
    clusters that were not re-randomized at stage 2.
7.  Cluster-level variables such as the treatments and candidate
    tailoring variables should be repeated across all individuals within
    the same cluster.

In other words:

- for Design I cSMART data, the stage-2 treatment column is fully
  observed;
- for Design II or III cSMART data, the stage-2 treatment column
  contains `NA` for non-re-randomized clusters.
- if the user only has the final observed outcome `Y`, then the stage-1
  pseudo-outcome is built from the stage-2 pseudo-outcome alone;
- if the user has both `Y1` and `Y`, then the stage-1 pseudo-outcome is
  built as `Y1 +` the stage-2 pseudo-outcome.

The package can work with any binary coding of the two treatment
columns. If the treatments are not already coded as `-1` and `1`,
`cQL()` internally recodes them and stores the mapping in the fitted
object.

The stage-1 formula should always use `Y1_tilde` as its response, for
example `Y1_tilde ~ X1 * A1`. If `Y1` is available, the stage-2 formula
may include `Y1` as a predictor.

The function also lets the user choose the significance level through
`alpha`. If `alpha = NULL`, the package uses the default `alpha = 0.05`,
which gives 95% confidence intervals. For example, set `alpha = 0.10`
for 90% confidence intervals.

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

Below, we generate a toy Design I cSMART dataset. In this design, all
clusters are re-randomized at stage 2, so the stage-2 treatment column
has no missing values.

### Step 1: generate a full re-randomization dataset

``` r
design1_data <- simulate_csmart_data(
  n_clusters = 40,
  cluster_size = 20,
  rerandomization = "full",
  seed = 111
)

head(design1_data)
#>   cluster_id patient_id X1 A1 X2 A2 response_status rerandomized        Y
#> 1          1          1  1  1  1  1              NA            1 2.706349
#> 2          1          2  1  1  1  1              NA            1 4.285709
#> 3          1          3  1  1  1  1              NA            1 1.089771
#> 4          1          4  1  1  1  1              NA            1 2.328823
#> 5          1          5  1  1  1  1              NA            1 3.537776
#> 6          1          6  1  1  1  1              NA            1 1.249130
table(is.na(design1_data$A2))
#> 
#> FALSE 
#>   800
```

The table above should show that `A2` is never missing, which is what we
expect for Design I.

### Step 2: specify the Q-functions

``` r
s2_formula <- Y ~ X1 * A1 + A1 * A2 + A2:X2
s1_formula <- Y1_tilde ~ X1 * A1
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
  working_correlation = "exchangeable",
  alpha = NULL,
  n_boot = 150,
  seed = 412,
  verbose = FALSE
)
```

Here we use the exchangeable working correlation model by specifying
`working_correlation = "exchangeable"`, which is also the default. If
desired, the user can instead set
`working_correlation = "independence"`.

### Step 4: inspect how the algorithm works internally

``` r
fit_design1$stage_summary
#>     stage  N N_rand  M                                        bootstrap
#> 1 Stage 2 40     40 40 Cluster bootstrap on stage-2 randomized clusters
#> 2 Stage 1 40     40 38                     M-out-of-N cluster bootstrap
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
#>          term   estimate    conf.low conf.high  std.error    p.value
#> 1 (Intercept) 0.50123442  0.32607911 0.6501724 0.08184193 0.01324503
#> 2          X1 0.67502287  0.52848726 0.7825857 0.06418361 0.01324503
#> 3          A1 0.25006349  0.10566007 0.3666382 0.07306235 0.01324503
#> 4          A2 0.27970228  0.15297888 0.4273272 0.06991907 0.01324503
#> 5       X1:A1 0.06841392 -0.05378174 0.2066801 0.06925594 0.33112583
#> 6       A1:A2 0.30846504  0.17474942 0.4430473 0.07198258 0.01324503
#> 7       A2:X2 0.34480828  0.21617948 0.5095903 0.06857789 0.01324503
#>   significance
#> 1            *
#> 2            *
#> 3            *
#> 4            *
#> 5             
#> 6            *
#> 7            *
```

To interpret the stage-2 table, focus especially on the terms involving
`A2`. If an interaction involving `A2` has a confidence interval that
excludes zero and a small p-value, that suggests the corresponding
variable may be useful as a stage-2 tailoring variable.

### Step 6: inspect the stage-1 inference

``` r
fit_design1$stage1
#>          term   estimate   conf.low conf.high  std.error    p.value
#> 1 (Intercept) 0.90998435  0.6874282 1.0947726 0.11091811 0.01324503
#> 2          X1 0.66366705  0.5121467 0.7849744 0.06881373 0.01324503
#> 3          A1 0.30633507  0.1244976 0.4751951 0.08885414 0.01324503
#> 4       X1:A1 0.05514058 -0.0937298 0.1616316 0.07016592 0.47682119
#>   significance
#> 1            *
#> 2            *
#> 3            *
#> 4
```

The stage-1 table is built after constructing the stage-1 pseudo-outcome
from the fitted stage-2 model. In this example there is no observed
`Y1`, so `Y1_tilde` is the stage-2 pseudo-outcome itself. Terms
involving `A1` are the stage-1 tailoring effects of interest.

## Example 2: Design II or III cSMART with limited second-stage re-randomization

Now we generate a toy cSMART dataset in which only a subset of clusters
is re-randomized at stage 2. This mimics the data format required for
Design II or III cSMARTs.

### Step 1: generate a partial re-randomization dataset

``` r
design23_data <- simulate_csmart_data(
  n_clusters = 40,
  cluster_size = 20,
  rerandomization = "nonresponder",
  p_rerand = 0.7,
  seed = 222
)

head(design23_data)
#>   cluster_id patient_id X1 A1 X2 A2 response_status rerandomized          Y
#> 1          1          1 -1 -1 -1 -1               0            1  0.5567023
#> 2          1          2 -1 -1 -1 -1               0            1  1.0253915
#> 3          1          3 -1 -1 -1 -1               0            1  0.5940704
#> 4          1          4 -1 -1 -1 -1               0            1 -0.4244146
#> 5          1          5 -1 -1 -1 -1               0            1 -1.3192431
#> 6          1          6 -1 -1 -1 -1               0            1 -0.9299379
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
#> Y1_tilde ~ X1 * A1
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
  working_correlation = "exchangeable",
  alpha = NULL,
  n_boot = 150,
  seed = 412,
  verbose = FALSE
)
```

### Step 4: inspect the stage summary

``` r
fit_design23$stage_summary
#>     stage  N N_rand  M                                        bootstrap
#> 1 Stage 2 40     26 26 Cluster bootstrap on stage-2 randomized clusters
#> 2 Stage 1 40     40 39                     M-out-of-N cluster bootstrap
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
#>          term  estimate     conf.low conf.high  std.error    p.value
#> 1 (Intercept) 0.5829939  0.355488360 0.8466682 0.11767260 0.01324503
#> 2          X1 0.6734640  0.440746446 0.8941773 0.11612442 0.01324503
#> 3          A1 0.1853665 -0.067189568 0.4243392 0.13019465 0.17218543
#> 4          A2 0.5662042  0.425349181 0.6897801 0.07342242 0.01324503
#> 5       X1:A1 0.1912644  0.001678435 0.4007254 0.10940519 0.06622517
#> 6       A1:A2 0.1014308 -0.074627778 0.2570886 0.08823046 0.29139073
#> 7       A2:X2 0.3173520  0.048673612 0.5304483 0.12687627 0.03973510
#>   significance
#> 1            *
#> 2            *
#> 3             
#> 4            *
#> 5             
#> 6             
#> 7            *
```

The interpretation is the same as before, but now the stage-2 regression
is fit only to the clusters that were actually re-randomized at stage 2.

### Step 6: inspect the stage-1 inference

``` r
fit_design23$stage1
#>          term  estimate   conf.low conf.high  std.error    p.value significance
#> 1 (Intercept) 0.8873112  0.6761015 1.0409854 0.09915426 0.01324503            *
#> 2          X1 0.7327109  0.5535310 0.9471056 0.10558044 0.01324503            *
#> 3          A1 0.4909462  0.2092183 0.7740152 0.14766397 0.01324503            *
#> 4       X1:A1 0.1349813 -0.0460383 0.3416661 0.09649030 0.19867550
```

For limited second-stage re-randomization, the manuscript’s stage-1
pseudo-outcome rule is used:

- for clusters re-randomized at stage 2, the pseudo-outcome uses the
  fitted stage-2 Q-function;
- for clusters not re-randomized at stage 2, the pseudo-outcome equals
  the observed outcome.

## Example 3: Design II cSMART with an observed stage-1 outcome `Y1`

This example uses limited second-stage re-randomization again, but now
the input data include both an observed stage-1 intermediate outcome
`Y1` and the final outcome `Y`.

### Step 1: create a Design II-style dataset with both `Y1` and `Y`

``` r
design23_y1_data <- simulate_csmart_data(
  n_clusters = 40,
  cluster_size = 20,
  rerandomization = "nonresponder",
  p_rerand = 0.7,
  seed = 333
)

set.seed(412)
cluster_effect_y1 <- stats::rnorm(
  length(unique(design23_y1_data$cluster_id)),
  sd = 0.25
)
names(cluster_effect_y1) <- as.character(unique(design23_y1_data$cluster_id))

design23_y1_data$Y1 <- with(
  design23_y1_data,
  0.4 +
    0.5 * X1 +
    0.3 * A1 +
    0.2 * X1 * A1 +
    cluster_effect_y1[as.character(cluster_id)] +
    stats::rnorm(nrow(design23_y1_data), sd = 0.6)
)

design23_y1_data$Y <- design23_y1_data$Y + 0.35 * design23_y1_data$Y1

head(design23_y1_data[c("cluster_id", "patient_id", "X1", "A1", "X2", "A2", "Y1", "Y")])
#>   cluster_id patient_id X1 A1 X2 A2          Y1         Y
#> 1          1          1  1 -1 -1  1  0.13237210  1.106925
#> 2          1          2  1 -1 -1  1  0.74401271 -1.149350
#> 3          1          3  1 -1 -1  1  0.76141055  1.883481
#> 4          1          4  1 -1 -1  1  0.85311440  1.452291
#> 5          1          5  1 -1 -1  1  0.01521741  1.556063
#> 6          1          6  1 -1 -1  1 -0.37251960  1.418785
```

### Step 2: specify the Q-functions

``` r
s2_formula_y1 <- Y ~ Y1 + X1 * A1 + A1 * A2 + A2:X2
s1_formula_y1 <- Y1_tilde ~ X1 * A1
```

Here `Y1` is allowed in the stage-2 formula, but the stage-1 formula
still uses `Y1_tilde` as its response.

### Step 3: fit the clustered Q-learning analysis

``` r
fit_design23_y1 <- cQL(
  data = design23_y1_data,
  stage2_formula = s2_formula_y1,
  stage1_formula = s1_formula_y1,
  cluster = "cluster_id",
  stage1_treat = "A1",
  stage2_treat = "A2",
  stage2_tailoring_vars = c("A1", "X2"),
  working_correlation = "exchangeable",
  alpha = NULL,
  n_boot = 150,
  seed = 412,
  verbose = FALSE
)
```

### Step 4: inspect the results

``` r
fit_design23_y1$stage_summary
#>     stage  N N_rand  M                                        bootstrap
#> 1 Stage 2 40     26 26 Cluster bootstrap on stage-2 randomized clusters
#> 2 Stage 1 40     40 38                     M-out-of-N cluster bootstrap
fit_design23_y1$stage2
#>          term  estimate    conf.low conf.high  std.error    p.value
#> 1 (Intercept) 0.6969863  0.44862376 0.9197475 0.12875327 0.01324503
#> 2          Y1 0.2530020  0.06700422 0.4010286 0.09103561 0.02649007
#> 3          X1 0.5902799  0.35923198 0.8030698 0.11440995 0.01324503
#> 4          A1 0.4161941  0.16176055 0.6595230 0.13818738 0.02649007
#> 5          A2 0.4068680  0.15723999 0.6569425 0.13294951 0.02649007
#> 6       X1:A1 0.2474096 -0.01332455 0.4632225 0.12647576 0.09271523
#> 7       A1:A2 0.4300958  0.16377534 0.7011124 0.14156313 0.01324503
#> 8       A2:X2 0.3728801  0.10805505 0.5271479 0.10716127 0.01324503
#>   significance
#> 1            *
#> 2            *
#> 3            *
#> 4            *
#> 5            *
#> 6             
#> 7            *
#> 8            *
fit_design23_y1$stage1
#>          term  estimate  conf.low conf.high std.error    p.value significance
#> 1 (Intercept) 1.5749700 1.2395897 1.8775642 0.1589357 0.01324503            *
#> 2          X1 1.2922897 0.9842671 1.5013229 0.1274876 0.01324503            *
#> 3          A1 0.9664480 0.6947004 1.3110249 0.1511807 0.01324503            *
#> 4       X1:A1 0.5808454 0.2285540 0.8005064 0.1425428 0.01324503            *
```

In this scenario the stage-1 pseudo-outcome is:

- `Y1 +` the fitted stage-2 pseudo-outcome for clusters re-randomized at
  stage 2;
- `Y1 + Y` for clusters not re-randomized at stage 2.

## Applying `cQL()` to your own data

After your data are prepared, the analysis for your own cSMART dataset
will look like this:

``` r
my_fit <- cQL(
  data = my_csmart_data,
  stage2_formula = Y ~ X1 * A1 + A1 * A2 + A2:X2,
  stage1_formula = Y1_tilde ~ X1 * A1,
  cluster = "cluster_id",
  stage1_treat = "A1",
  stage2_treat = "A2",
  stage2_tailoring_vars = c("A1", "X2"),
  working_correlation = "exchangeable",
  alpha = NULL,
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

If your data contain an observed stage-1 outcome `Y1`, keep the stage-1
formula as `Y1_tilde ~ ...` and optionally include `Y1` in
`stage2_formula`.

If you want a confidence level other than 95%, specify `alpha` directly.
For example, use `alpha = 0.10` to request 90% confidence intervals.

## Notes

- The current package targets the two-stage cSMART setting developed in
  the manuscript.
- The stage-2 formula must include the main effect of the stage-2
  treatment.
- Every interaction involving the stage-2 treatment must be a two-way
  interaction with one of the variables listed in
  `stage2_tailoring_vars`.
- In the printed stage summary, the stage-2 row has `M = N_rand` because
  stage 2 uses the full cluster bootstrap over the stage-2 randomized
  clusters, whereas the stage-1 row uses the selected M-out-of-N
  resample size. See more technical detail about the algorithm design in
  [the original manuscript](https://arxiv.org/abs/2505.00822)
