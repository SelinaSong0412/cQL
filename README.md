
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
- 95% confidence intervals;
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
4.  If the cSMART is a Design II or III trial with limited second-stage
    re-randomization, then the stage-2 treatment should be `NA` for
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

Below, we generate a toy Design I cSMART dataset. In this design, all
clusters are re-randomized at stage 2, so the stage-2 treatment column
has no missing values.

### Step 1: generate a full re-randomization dataset

``` r
design1_data <- simulate_csmart_data(
  n_clusters = 40,
  cluster_size = 12,
  rerandomization = "full",
  seed = 412
)

head(design1_data)
#>   cluster_id patient_id X1 A1 X2 A2 response_status rerandomized          Y
#> 1          1          1 -1 -1 -1  1              NA            1 -0.2008901
#> 2          1          2 -1 -1 -1  1              NA            1 -2.0253075
#> 3          1          3 -1 -1 -1  1              NA            1 -3.0686516
#> 4          1          4 -1 -1 -1  1              NA            1 -0.2653847
#> 5          1          5 -1 -1 -1  1              NA            1 -1.1638457
#> 6          1          6 -1 -1 -1  1              NA            1 -1.4434476
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
  seed = 412,
  verbose = FALSE
)
```

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
#>          term  estimate   conf.low conf.high  std.error    p.value significance
#> 1 (Intercept) 0.4205342 0.28617359 0.5716164 0.07163976 0.01324503            *
#> 2          X1 0.5414709 0.39956939 0.6995669 0.07377064 0.01324503            *
#> 3          A1 0.3171185 0.17776460 0.4375005 0.06732913 0.01324503            *
#> 4          A2 0.4927690 0.36868516 0.6154612 0.06530733 0.01324503            *
#> 5       X1:A1 0.1663990 0.02967405 0.3076275 0.07203819 0.03973510            *
#> 6       A1:A2 0.3002660 0.16664787 0.4224097 0.06740352 0.01324503            *
#> 7       A2:X2 0.2711164 0.13636781 0.4161889 0.07049529 0.01324503            *
```

To interpret the stage-2 table, focus especially on the terms involving
`A2`. If an interaction involving `A2` has a confidence interval that
excludes zero and a small p-value, that suggests the corresponding
variable may be useful as a stage-2 tailoring variable.

### Step 6: inspect the stage-1 inference

``` r
fit_design1$stage1
#>          term  estimate   conf.low conf.high  std.error    p.value significance
#> 1 (Intercept) 0.9505510 0.75090683 1.1733883 0.10211132 0.01324503            *
#> 2          X1 0.6211356 0.42074962 0.7883756 0.09538244 0.01324503            *
#> 3          A1 0.6275822 0.43279404 0.8157717 0.10339072 0.01324503            *
#> 4       X1:A1 0.2426261 0.08182366 0.4074709 0.08635830 0.01324503            *
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
  n_clusters = 40,
  cluster_size = 20,
  rerandomization = "responder",
  p_rerand = 0.65,
  seed = 2026
)

head(design23_data)
#>   cluster_id patient_id X1 A1 X2 A2 response_status rerandomized         Y
#> 1          1          1 -1  1  1 NA               0            0 1.6633267
#> 2          1          2 -1  1  1 NA               0            0 0.7023146
#> 3          1          3 -1  1  1 NA               0            0 1.0616128
#> 4          1          4 -1  1  1 NA               0            0 1.7423697
#> 5          1          5 -1  1  1 NA               0            0 2.1988416
#> 6          1          6 -1  1  1 NA               0            0 1.1747294
table(is.na(design23_data$A2))
#> 
#> FALSE  TRUE 
#>   200   600
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
#> <environment: 0x1582e7278>
s1_formula
#> Y1 ~ X1 * A1
#> <environment: 0x1582e7278>
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
  seed = 2026,
  verbose = FALSE
)
#> Warning: Only 140 successful bootstrap refits were obtained out of 150.
```

### Step 4: inspect the stage summary

``` r
fit_design23$stage_summary
#>     stage  N N_rand  M                                        bootstrap
#> 1 Stage 2 40     10 10 Cluster bootstrap on stage-2 randomized clusters
#> 2 Stage 1 40     40 40                     M-out-of-N cluster bootstrap
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
#>          term  estimate    conf.low conf.high  std.error   p.value significance
#> 1 (Intercept) 0.6867430  0.53043776 0.8528338 0.09988243 0.0141844            *
#> 2          X1 0.4600669  0.39629319 0.5140551 0.04822348 0.0141844            *
#> 3          A1 0.1378116 -0.02827919 0.2941169 0.10485366 0.2553191             
#> 4          A2 0.4738526  0.25885386 0.6939315 0.11981833 0.0141844            *
#> 5       X1:A1 0.2078228  0.15383462 0.2715965 0.04822348 0.0141844            *
#> 6       A1:A2 0.1920257  0.03572043 0.3581165 0.10156664 0.0141844            *
#> 7       A2:X2 0.2350227  0.17124895 0.4067727 0.10126986 0.0141844            *
```

The interpretation is the same as before, but now the stage-2 regression
is fit only to the clusters that were actually re-randomized at stage 2.

### Step 6: inspect the stage-1 inference

``` r
fit_design23$stage1
#>          term  estimate   conf.low conf.high  std.error    p.value significance
#> 1 (Intercept) 0.9279898 0.71522973 1.1234748 0.10765712 0.01324503            *
#> 2          X1 0.4250531 0.24635415 0.5871342 0.08953314 0.01324503            *
#> 3          A1 0.5746654 0.41975261 0.7790176 0.09203841 0.01324503            *
#> 4       X1:A1 0.2405181 0.07771541 0.4642341 0.09766087 0.01324503            *
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
  stage 2 uses the full cluster bootstrap over the stage-2 randomized
  clusters, whereas the stage-1 row uses the selected M-out-of-N
  resample size. See more technical detail about the algorithm design in
  [the original manuscript](https://arxiv.org/abs/2505.00822)
