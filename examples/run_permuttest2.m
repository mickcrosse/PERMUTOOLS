function run_permuttest2
%RUN_PERMUTTEST2  Run unpaired permutation test examples.
%   Example 1: generate multivariate data for 2 samples, each with 20
%   variables and 30 observations and perform unpaired permutation tests
%   between the corresponding variables of each sample.
%
%   Example 2: generate univariate data for 5 samples, each with 30
%   observations and perform unpaired permutation tests between every
%   sample (5 samples = 10 comparisons). Note that each column of X
%   represents an independent sample and may contain NaNs for samples with
%   smaller number of observations.
%
%   See also PERMUTTEST2 RUN_PERMUTTEST RUN_PERMUTCORR.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

% Multivariate Data

% Generate random data
x = randn(30,20);
y = randn(30,20);

% Make certain variables have different means
y(:,1:8) = y(:,1:8)-1;

% Compute adjusted test statistics for unpaired test
[t,stats] = permuttest2(x,y)

% Display results
disp('Multivariate Data')
disp('t:')
disp(t)
disp('stats:')
disp(stats)

% Permutation Test Matrix

% Generate random data
x = randn(30,5);

% Make certain variables have different means
x(:,3:5) = x(:,3:5)-1;

% Compute adjusted test statistics for unpaired test
[t,stats] = permuttest2(x)

% Display results
disp('Permutation Test Matrix')
disp('t:')
disp(t)
disp('stats:')
disp(stats)