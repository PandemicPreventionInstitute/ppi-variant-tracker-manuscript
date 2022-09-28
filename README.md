# ppi-variant-tracker-manuscript
This repo provides the code needed to generate the figures and outputs in
the manuscript. The code calls the GISAID SARS-COv-2 metadata, obtained at 6 
different time points between April 28th, 2022 and July 1st, 2022, for the 
bulk of the analysis. A one time pull of the GISAID flu metadata is also used.
This data is available for authenticated GISAID users, and can be regenerated
as it is used here by filtering by submission date in the line list data. 
We thank GISAID for making this data publicly available. 

## Pipeline for historical validation 
1. Run `code/auto_extract_gisaid.R` to get updated gisaid line list metadata in the 
   form of a csv. For the 6 reference datasets, we manually time stamped these 
   corresponding to the day they were generated. These are saved into 
   `../data/raw/{reference_date}_metadata.csv`. Alternatively, generate these from 
   filtering the current metadata by submission date, for each reference date. 
2. Run `historical_validation_code/pre_processing_HV.R`, which loads in the
   `../data/raw/{reference_date}_metadata.csv` data, finds a consensus set
   of lineages across all 6 reference datasets, and saves timestamped aggregated
   outputs of the counts of sequences of a particular lineage in each country on each
   day as `../data/processed/validation_data/lineage_t_{reference_date}.csv
3. Run `historical_validation_code/prepare_data_for_cmdstan_HV.R` which loads in
   the lineage-country-timepoint data from each reference date 
   `../data/processed/validation_data/lineage_t_{reference_date} and loops through
   and transforms it into `json` object for cmdstan as 
   `../data/processed/validation_data/data_for_cmdstan_{reference_date}.json`. It also
   generates the csv data needed for the single country MLE estimation model at each 
   reference data as `../data/processed/validation_data/data_for_nnet_{reference_date}.csv
4. Run the multicountry model using `historical_validation_code/run_model_using_
   shell_scripts.R` which reads in the reference date json objects `../data/processed
   /validation_data/data_for_cmdstan_{reference_date}.json`, resets the directory
   to the ppi-variant-tracker-manuscript folder, runs `run_cmdstan_HV.sh` which runs
   each model (takes about 4 hours on Mac OS XXX for each reference dataset). See the 
   next section for more information on how to get the model set up in 
   `cmdstan`. We **strongly** recommend using `cmdstan` or, if not on an M1 
   Mac `CmdStanR`, over RStan because of some improvements in Stan's HMC 
   sampler that are slow to update in RStan/CRAN.
   - If the first time the model is being run, it needs to be compiled.  See 
     the next section for an explanation of how to compile a stan model with 
     `cmdstan`. 
     The script then runs `process_output_HV.sh` which generates mode output files in
   `../data/output/multicountry_output/validation/processed_output_{reference_date}.csv`.
   Within this R script, the `$p_{hat}$` (model estimated prevalence) and `$\tilde{Y}$`
   (model estimated observed sequences of each lineage) get summarized for each reference 
   date. Those are saved here: `../data/output/multicountry_output/validation/
   p_hat_{reference_date}.csv` and `../data/output/multicountry_output/validation/
   Y_tilde_{reference_date}.csv`
 5. Run `historical_validation_code/process_results_from_stansummary_HV.R which loads in
   the lineage-country-timepoints from each reference date, and the summarized model
   estimated variant dynamics, joining them together so that the countries and lineages
   are properly matched. It generates a time stamped `../data/processed/validation_data/
   clean_global_df_{reference_date}.csv` that contains the model estimates summarized lineage
   prevalences alongside the observed data.
 6. Run `historical_validation_code/process_results_from_cmdstan_HV.R` which loads
   in the lineage-country-timepoints from each reference date and the raw model output. 
   This script computes the Briar score for each country and each country-lineage based on the 
   predicted variant prevalence ($p_{hat}$) and the observed variant prevalence. These are
   saved as `../data/processed/validation_data/country_metrics_{reference_date}.csv` and
   `../data/processed/validation_data/country_lineage_metrics_{reference_date}.csv`. 
   This script also generates the summary stats and full posterior distributions of the 
   relative fitness advantage estimates at the country (`../data/processed/validation_data/
   `r_summary_{reference_date}.csv` and `../data/processed/validation_data/r_distrib_
   {reference_date}.csv`) and at the global level (`../data/processed/validation_data/
   mu_hat_{reference_date}.csv and `../data/processed/validation_data/mu_distrib_
   {reference_date}.csv`).
 7. Run `historical_validation_code/evaluate_output_HV.R` which loads in all of the time stamped
   model outputs, concatenates them together and combines them with the final one 
   (from July 1st, 2022) for direct comparison. This generates a single dataset for each 
   output type, that contains a column for the reference date. These outputs are:
   `../data/output/validation/r_summary.csv`, `../data/output/validation/clean_global_df.csv`,
    `../data/output/validation/mu_hat.csv`, `../data/output/validation/r_distrib.csv`, 
    `../data/output/validation/country_metrics.csv`, `../data/output/validation/country_
    lineage_metrics.csv`

These are loaded into `code/figures.R` to generate Figure 4 in the main text of the manuscript and
multiple supplemental figures. 


## Setting up the Stan model

We strongly recommend running the model using `cmdstan`. The RStan sampler 
will probably work, but it's several generations behind the version released 
for `cmdstan` due to issues getting the RStan updates approved for CRAN.  
Although `CmdStanR` is an excellent workaround for R users, it seems to have 
a serious bug for M1 Macs (like those used by the PPI Data team). As a result, 
we use `cmdstan` in our workflow and will recommend that here.

We recommend _carefully_ following the [setup instructions for 
cmdstan](https://mc-stan.org/docs/2_29/cmdstan-guide/cmdstan-installation.html).  
However, do not use the Anaconda installation approach! We've had serious 
issues getting it set up through conda and has lead to subsequent issues with 
`clang` compilers. Instead, we strongly suggest using the Source installation 
suggested in section 1.2.

Once the model is set up as described in section 1.2 and validated as 
described in section 1.3, `cmdstan` is ready to build the model. This is 
clearly described in [section 3.1 of the cmdstan 
guide](https://mc-stan.org/docs/2_29/cmdstan-guide/compiling-a-stan-program.html#invoking-the-make-utility).  
In brief, the model can be generated using the `make` program from the 
`cmdstan` directory. In other words, you can cd to the `cmdstan` directory and 
run `make 
path/to/ppi-variant-tracker/code/stancode/multivariate_variant_multinomial_ncp`.  
Note, the `.stan` extension is omitted.

This should produce two objects in the `stancode/` directory, the 
`multivariate_variant_multinomial_ncp` executable and the 
`multivariate_variant_multinomial_ncp.hpp` header file. These indicate that 
the model has successfully compiled and is ready to be run.

## Variant tracking model

The variant model uses SARS-CoV-2 sequence metadata from GISAID to produce 
estimates of relative variant growth rates and true proportion of total cases 
in a given country. It does this with a hierarchical Bayesian approach, 
treating the observed sequence data as multinomial and the observed variant 
dynamics in each country as related to those in other countries. This approach 
has a bit of the flavor of a Bradley-Terry model or the model presented on 
pages 423-425 of BDA3. 

Considering first the case for a variants within a single country, if we 
consider all sequences observed over the past 6 months (6 months ago 
$t = 0$; $t \in \mathbb{W}$), on day $t$ we observe $N_t$ sequences total 
($N_t \in \mathbb{W}$), of which $n_{it}$ sequences are of variant lineage $i$ 
($n_{it} \in \mathbb{W}; N_t = \sum_i n_{it}$). Then the set of observed 
variant lineages on that day is $i_t \in \{1, ..., I\}$ where I is the number 
of observed lineages over the entire time period.  Although this quantity 
changes over time due to evolution, we consider it here to be both fixed and 
known (i.e. the number of PANGO lineages known today). If
$p_{it} = \frac{n_{it}}{N_t}$, then the set of observed variant proportions on 
day $t$ is $p_t = \{p_{it}, ..., p_{It}\}$, $p \in [0, 1]$. Assuming that the 
random draw $n_{i, t}$ is independent (albeit correlated) with the random draw 
$n_{i, t+1}$, this scenario describes a multinomial distribution, which has 
probability mass function (PMF):

$$f(n_{it}, ..., n_{It}; N_t, p_{it}, ..., p_{It}) = \frac{N_t!}{n_{it}! ...  
n_{It}!} p_{it} \times ... \times p_{It}$$

However, this PMF conditions on the _known_ true probabilities of drawing $i$; these probabilities are of course unknown. With vector of lineages 
$\{i, ..., j, I\}$, where $I$ is the dominant local variant, this PMF and its unknowns motivate the regression model formulation

$$ln\left(\frac{P(X = i)}{P(X = I)} \right)= \beta_{0i} + \beta_{1i} t$$
$$\vdots$$
$$ln\left(\frac{P(X = j)}{P(X = I)} \right)= \beta_{0j} + \beta_{1j} t$$

for each category $i$ except the base case $I$, producing $I-1$ cases. 
Classically, the intercepts ($\beta_{0i}$) in a multinomial logistic 
regression are interpreted as the average difference in log odds being in 
category $i$ than category $I$ when the coefficients are 0 (i.e. $t = 0$). The 
slope, $\beta_{1i}$ is interpreted as the average change in log odds of being 
in category $i$ relative to category $I$ for a unit change in the covariate 
$t$. Fit to observed variant proportions, the slope coefficients describe the 
difference in intrinsic growth rates for the variant in the numerator and the 
variant in the denominator (i.e. $\beta_{1i} = r_i - r_I$ where $r_I$ is the 
intrinsic growth rate of the dominant variant). However, this model does not 
produce an estimate of the actual intrinsic growth rates (i.e. $r_i$ or $r_I$) 
because fitting $I$ cases would make the model overdetermined.

Despite the ability of the multinomial likelihood to account for changes in 
sample size over time (i.e., $N_t$), the $\beta_{1i}$ coefficients are still 
likely to be biased at small sample sizes or when first observed. Many of 
these variants will be just as fit or roughly similar to the current dominant 
variants (i.e. neutral evolution), but we are only able to observe the mutant 
lineages that grow more quickly than the dominant lineage -- regardless of 
whether this apparent growth advantage is truly real or simply stochastic and 
transient. As a result, of the variants that are growing quickest at time $t$, 
many will, on average, have MLE $\hat \beta_{1i}$ coefficients stochastically 
higher than their true average growth rate $\beta_{1i}$. If they were not 
growing more quickly at the current moment, we wouldn't be able to observe 
them. In other words, our observations are a biased subset of the entire 
population -- but this bias is transient and will decrease with time and 
number of samples. It's a function of our observation process. If we just fit 
a multinomial to these samples, we're almost certainly going to see our 
estimates decline over time (regress to the mean) as we get more data because 
we're inducing a selection bias in the data we're considering.


We can decrease this bias by specifying a model that shrinks these small 
sample size, uncertain $\hat \beta_{1i}$ coefficients away from more extreme 
estimates and towards values we think are more likely. However, in order to 
shrink the inflated $\hat \beta_{1i}$ coefficients to smaller values, we need 
to specify: 1) what we are shrinking these coefficients towards and 2) the 
strength of the shrinkage. In other words, we need specify what it means to be 
extreme and what it means to reasonable. In a Bayesian model, we perform this 
shrinkage through two separate, but related mechanisms:

1. Specifying a candidate generative process through which the data are 
   observed (AKA the model or, more formally, the likelihood)
2. Specifying our prior beliefs (more formally, the prior likelihood). 

One obvious choice for a "reasonable" model of $\beta_{1i}$ that would perform 
shrinkage is one in which most new lineages are expected to be basically 
similar to the lineages that have already occurred, unless the data strongly 
suggest otherwise.  Biologically this model would make sense; most new 
lineages have a few genetic mutations that differentiate them from their 
parent lineages, but otherwise remain the same.  

It seems reasonable, therefore, that these new lineages should have growth 
rates similar to their parent lineages most of the time, sometimes growing 
little bit quicker and sometimes a little bit slower, but occasionally 
a particularly impactful constellation of mutations will lead to a radical 
change. More technically, this description is of a generative process in which 
most $\beta_{1i}$ values are similar to the overall mean of the growth rates, 
$\mu_{\beta_1}$. This is the value towards which we would shrink $\beta_{1i}$ 
values. The strength of the shrinkage is controlled by $\sigma^2_{\beta_1}$. 
the variance of the $\beta_{1i}$ values. If most $\beta_{1i}$ values are 
tightly clustered around $\mu_{\beta_1}$, their variance would be low and we 
would be surprised to see a very high $\beta_{1i}$ value, so we would see more 
shrinkage of that outlier. On the other hand, if there's high variance in 
$\beta_{1i}$ values, we would see less shrinkage because an outlier would be 
more expected. Even more technically, this process suggests a unimodal, 
symmetric distribution over the $\beta_{1i}$ values.  In building this model 
we should evaluate a range of distributions that fit this description (i.e.  
Student's $t$ distributions with varying degrees of freedom), but for the rest 
of this document I'll use a Normal distribution for notational simplicity and 
clarity.  While this solution is imperfect because it does not model the 
actual Bernoulli observation process (which would require modeling both the 
observed $\hat \beta_{1i}$ and $\beta_{1i}$ as well as the probability of 
observing $i$ as an increasing function of $\hat \beta_{1i}$), it will help to 
somewhat reduce the bias.

If we consider relative variant fitness as drawn from a distribution 
$\beta_{1i} \sim N(\mu_{\beta_1}, \sigma^2_{\beta_1})$, new variants with 
small sample sizes will be pulled toward $\mu_{\beta_1}$ with the strength of 
the shrinkage controlled by $\sigma^2_{\beta_1}$. The $\sigma^2_{\beta_1}$ 
parameter in the model will be informed by both observed data and our prior 
expectation of its value. Since these differences in growth rates are relative 
to the dominant variant as the base case, a $N(0, 0.5)$ prior distribution is 
likely appropriate for $N(\mu_{\beta_1}, \sigma^2_{\beta_1})$, which assumes 
that a new variant's intrinsic growth is 50\% likely to be better than the 
base case and has a 95% probability to be $\pm$ 1 from the mean intrinsic 
growth rate. If the data are looser or wider than this prior value, the 
$\sigma^2_{\beta_1}$ parameter will be updated during the model by the 
information present in the data.

## Multi-country extension

Comparison of intrinsic growth rates from logistic growth models across time 
or location is complicated by the dependence of these growth rates on local 
conditions. One strain may be fitter in the presence of particular conditions, 
driven by the contact or immunity landscape, but these conditions may be 
transient or particular to that location. A strain that is more fit in one 
setting is not necessarily more fit in another.

However, if one is willing to relax the interpretability of their model, one 
could extend these multinomial models across countries -- borrowing 
information from conditions in one country to inform the model about those of 
another. This approach can help handle problems of limited sampling in one 
country and, vice versa, help increase the uncertainty about relative fitness 
to our true level of skepticism when it is artificially deflated by very high 
sampling in one or two locations.

This can be accomplished by generalizing the 
$\beta_{1i} \sim N(\mu_{\beta_1}, \sigma^2_{\beta_1})$ distribution into 
a multi-country extension: 
$\beta_{1i}^K \sim MVN(\mu_{\beta_1}^K, \Sigma_{\beta_1}^{K \times K})$.
Here, each variant lineage $I$ has a single 
$K$-length vector $\beta_{1i}^K$ as a single realization from the multivariate 
normal distribution, in which each element $\beta_{1i}$ is one variant lineage 
in country $k$ (from the set of all countries $\{k, ..., K\}$, so there are 
$I$ of these $K$ length vectors).  The mean over all lineages of the relative 
growth rate of variants in that country ($E[\beta_{1}|k, \{i, ..., I\}]$) is 
described by the $K$-length vector $\mu_{\beta_1}^K$, in which each element is 
the expected relative growth rate of a lineage in that country. Rather than 
a single $\sigma^2_{\beta_1}$ parameter, we can produce a $K \times K$ 
$\Sigma_{\beta_1}$ covariance matrix, describing the similarity of conditions 
in one country with those of another. The $\Sigma_{\beta_1}$ matrix can be 
decomposed into scale parameters $\tau$ and a correlation matrix $\Omega$ with 
a more interpretable LKJ prior. The expected variability of the relative 
growth rates of variants within the country is described by 
diag$(\Sigma_{\beta_1}^{K \times K})$, the same as the $\sigma^2_{\beta_1}$ 
values in the original model. The off-diagonal elements of
$\Omega^{K x K}$ describe the expected correlation of the relative growth in 
country $k$ 
with that of another country -- in other words the amount of information 
sharing across countries. If the observed relative growth rate of variant $i$ 
is higher in country $k$, the the off-diagonal elements of 
$\Omega^{K \times K}$ determine how strongly we should expect to be higher in other
countries too. 

In order to maintain the ability to share information across variants that 
have only been sparingly observed, the mean differences in growth rates are 
themselves drawn from a hierarchical distribution: 
$\mu_k \sim N(\mu_{\mu}, \sigma_{\mu})$. This hierarchical distribution is 
similar to the hierarchical 
distribution that implemented shrinkage on the growth rates in the single 
country model, but in this case it is instead on the mean growth rates 
averaged over countries. The intercepts, $\beta_{0, k}$ are each drawn from 
lineage specific normal distributions, with no correlation matrix across 
countries: $\beta_{0, i, k} \sim N(\mu_{0i}, \sigma_{0i})$. This approach 
allows for different timing, emergence, and patterns of spread for each 
variant while still enforcing a roughly similar introduction time for 
different variants across countries.


Less technically, we would be jointly estimating a vector of global mean 
growth rates for each variant ($\mu_{\beta_1}^K$) over all countries. The 
individual observations of the growth rates in each country would be shrunk 
towards the overall global mean with strength of shrinkage controlled by the 
country-specific variances across the countries (diag
$(\Sigma_{\beta_1}^{K \times K})$). Information sharing about variants is then 
also possible across countries, with the correlations between countries 
describing how strongly we expect individual countries to be like each other.  
If the correlation between countries $k$ and $m$ is high and if $\beta_{1i}$ 
is observed higher than $\mu_{\beta_1, k}$ in country $k$, then we would 
likewise expect it to be higher than $\mu_{\beta_1, m}$ in country $m$. We 
would be assuming that each variant (in every country) is a single observation 
from this multivariate normal distribution with the similarities across 
countries following the pattern described by its correlation/covariance 
matrix.

All the variant growth rates in a country are treated as a single draw from 
this distribution, with each country treated as a conditionally independent 
draw. The vector of (variant-specific) means of the multivariate normal 
distribution drawn from yet another hierarchical normal distribution.

$$\mu_i \sim N(\mu_{hierarchical}, \sigma_{hierarchical})$$

This approach shrinks the estimated growth rates of new variant with few known 
sequences to the overall mean. It has the effect of regularizing estimates of 
new variants, which are consistently overestimated -- likely because of 
a combination of selection bias and transient effects supporting their 
emergence before regression to the mean.

## Rt estimation

To estimate the $R_v(t$) of a variant in a country, we start by inferring the 
cases with a variant in each country. To do this, we take the daily reported 
cases in a country and multiply by each draw from the posterior mean 
prevalence estimate of that variant in that country. For each case trajectory, 
we follow the methods described in (Liu et al. 2020), focusing only on locally 
acquired infections. Briefly, this approach assumes that the number of locally 
acquired infections at time t depends on the transmissibility at time t (i.e.  
$R(t)$) and on the number of cases at any time before t, weighted by the 
distribution of the generation time. We assume that cases of a variant at time 
t ($C_v(t)$) are generated via a Poisson process according to the renewal 
equation:

$$C_v(t) \approx Pois(R_v(t) \sum_{s=1}^{t}\phi(s)C_v(t-s))$$

Where $\phi(s)$ is the generation time distribution and $R_v(t)$ is the 
effective reproductive number of the variant in that country at time t. For 
the generation time distribution, we assume this follows a gamma distribution 
with shape = 5.8 and rate = 0.43 [UPDATE & CITE]. The likelihood of the 
observed time seres of cases from day 1 to T can be written as:

$$L = \prod_{t=1}^{T} P(C_v(t), \lambda = R_v(t) \sum_{s=1}^{t} \phi(s)C_v(t-s))$$

Where $P(x,\lambda)$ is the Poisson density distribution of observing $x$ 
events given $lambda$ parameter. For each time point of the case trajectory, 
we independently estimate R_v(t). We then take the 7 day rolling average of 
the estimated R_v(t)s. We repeat the estimation procedure for 100 samples from 
the posterior estimated prevalence described above, and combine the posterior 
distributions from each case trajectory to capture the full uncertainty.

Liu Q-H, Bento AI, Yang, K, Zhang H, Yang X, Merler S, et al. (2020) The 
COVID-19 outbreak in Sichuan, China: Epidemiology and impact of interventions.  
PLoS Comput Biol 16(12):e1008467 
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008467
