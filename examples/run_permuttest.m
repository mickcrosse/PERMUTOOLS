function run_permuttest
%RUN_PERMUTTEST  Run paired permutation test examples.
%   Example 1: generate multivariate data for 2 conditions, each with 20
%   variables and 30 observations and perform paired permutation tests
%   between the corresponding variables of each condition.
%
%   Example 2: generate univariate data for 5 conditions, each with 30
%   observations and perform paired permutation tests between every pair of
%   conditions (5 conditions = 10 comparisons).
%
%   See also PERMUTTEST RUN_PERMUTTEST2 RUN_PERMUTCORR.
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

% Compute adjusted test statistics for paired-sample test
[t,stats] = permuttest(x,y)

% Compute adjusted test statistics for one-sample test
[t,stats] = permuttest(x)

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

% Compute adjusted test statistics for paired-sample test
[t,stats] = permuttest(x,[],'sample','paired')

% Display results
disp('Permutation Test Matrix')
disp('t:')
disp(t)
disp('stats:')
disp(stats)