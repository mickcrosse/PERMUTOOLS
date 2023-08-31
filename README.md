# <img src="docs/permutools_logo.png">

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

PERMUTOOLS is a statistical software package for multivariate permutation testing in MATLAB. Permutation tests based on the t-statistic can be performed on one-sample, paired-sample and independent two-sample datasets. For two-sample situations, a permutation test based on the F-statistic can be performed to determine variance equivalence. By comparing the magnitude of the test statistic of interest with those obtained using permutations of the data, it provides powerful, distribution-free hypothesis testing. Family-wise error rate (FWER) is controlled using the maximum statistic method (Blair et al., 1994), making it suitable for multivariate or multiple permutation tests. PERMUTOOLS offers a range of test statistics including the t-statistic (paired, unpaired), the F-statistic (unpaired), the correlation coefficient (Pearson, Spearman, rankit), as well as measures of effect size with bootstrapped confidence intervals (Cohen's d, Hedges' g, Glass' delta).

- [Installation](#installation)
- [Documentation](#documentation)
- [Max Statistic Correction](#max-statistic-correction)
- [Contents](#contents)
- [Examples](#examples)
- [License](#license)

## Installation

Download and unzip PERMUTOOLS to a local directory, then in the MATLAB/GNU Octave command window enter:

```matlab
addpath 'directory/PERMUTOOLS'
savepath
```

## Documentation

For documentation and citation, please refer to the [PERMUTOOLS paper](docs/Crosse_etal_2018.pdf).

For usage, please see [examples](#examples) and [example M-files](examples).

## Max Statistic Correction

The maximum statistic method provides strong control of FWER, even for small sample sizes, and is much more powerful than traditional correction methods (Gondan, 2010; Groppe et al., 2011a). For unpaired testing it is also rather insensitive to differences in population variance when samples of equal size are used (Groppe et al., 2011b). For samples of unequal size or variance, Welch's t-statistic may be used as it is less sensitive to differences in variance (but also less sensitive to differences in means). For nonlinear correlations, the raw data may be transformed to rank orders using a Spearman's or a rankit transformation (Bishara & Hittner, 2012).

## Contents

* `permuttest()` - one-sample or paired-sample permutation test with tmax correction
* `permuttest2()` - unpaired two-sample permutation test with tmax correction
* `permuvartest2()` - permutation-based F-test with max statistic correction
* `permucorr()` - permutation-based correlation with max statistic correction
* `deffectsize()` - effect size measure with bootstrapped confidence intervals

## Examples

### Paired and one-sample permutation testing

Here, we generate multivariate random data for 2 conditions, each with 20 variables and 30 observations and perform paired permutation tests between the corresponding variables of each condition.

```matlab
% Generate random data
x = randn(30,20);
y = randn(30,20);

% Make certain variables have different means
y(:,1:8) = y(:,1:8)-1;

% Compute adjusted test statistics for paired-sample test
[t,stats] = permuttest(x,y)

% Compute adjusted test statistics for one-sample test
[t,stats] = permuttest(x)
```

Here, we generate univariate random data for 5 conditions, each with 30 observations and perform paired permutation tests between every pair of conditions (5 conditions = 10 comparisons).

```matlab
% Generate random data
x = randn(30,5);

% Make certain variables have different means
x(:,3:5) = x(:,3:5)-1;

% Compute adjusted test statistics for paired-sample test
[t,stats] = permuttest(x,[],'sample','paired')
```

### Unpaired permutation testing

Here, generate multivariate random data for 2 samples, each with 20 variables and 30 observations and perform unpaired permutation tests between the corresponding variables of each sample.

```matlab
% Generate random data
x = randn(30,20);
y = randn(30,20);

% Make certain variables have different means
y(:,1:8) = y(:,1:8)-1;

% Compute adjusted test statistics for unpaired test
[t,stats] = permuttest2(x,y)
```

Here, we generate univariate random data for 5 samples, each with 30 observations and perform unpaired permutation tests between every sample (5 samples = 10 comparisons). Note that each column of X represents an independent sample and may contain NaNs for samples with a smaller number of observations.

```matlab
% Generate random data
x = randn(30,5);

% Make certain variables have different means
x(:,3:5) = x(:,3:5)-1;

% Compute adjusted test statistics for unpaired test
[t,stats] = permuttest2(x)
```

### Correlation testing

Here, we generate multivariate random data for 2 conditions, each with 20 variables and 30 observations and calculate the correlation between the corresponding variables of each condition.

```matlab
% Generate random data
x = randn(30,20);
y = randn(30,20);

% Make certain variables more correlated
y(:,1:5) = y(:,1:5)+0.5*x(:,1:5);
y(:,15:20) = y(:,15:20)-x(:,15:20);

% Compute correlation and adjusted test statistics
[r,stats] = permucorr(x,y)
```

Here, we generate univariate random data for 5 conditions, each with 1 variable and 30 observations and calculate the correlation between every pair of conditions (5 conditions = 10 correlations).

```matlab
% Generate random data
x = randn(30,5);

% Make certain variables more correlated
x(:,1:2) = x(:,1:2)-0.5*x(:,4:5);

% Compute correlation and adjusted test statistics
[r,stats] = permucorr(x)
```

## License

[BSD 3-Clause License](LICENSE)