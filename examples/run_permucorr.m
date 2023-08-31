function run_permucorr
%RUN_PERMUCORR  Run permutation-based correlation examples.
%   Example 1: generate multivariate data for 2 conditions, each with 20
%   variables and 30 observations and calculate the correlation between the
%   corresponding variables of each condition.
%
%   Example 2: generate univariate data for 5 conditions, each with 1
%   variable and 30 observations and calculate the correlation between
%   every pair of conditions (5 conditions = 10 correlations).
%
%   See also PERMUCORR RUN_PERMUTTEST RUN_PERMUTTEST2.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

% Multivariate Data

% Generate random data
x = randn(30,20);
y = randn(30,20);

% Make certain variables more correlated
y(:,1:5) = y(:,1:5)+0.5*x(:,1:5);
y(:,15:20) = y(:,15:20)-x(:,15:20);

% Compute correlation and adjusted test statistics
[r,stats] = permucorr(x,y);

% Display results
disp('Multivariate Data')
disp('r:')
disp(r)
disp('stats:')
disp(stats)

% Correlation Matrix

% Generate random data
x = randn(30,5);

% Make certain variables more correlated
x(:,1:2) = x(:,1:2)-0.5*x(:,4:5);

% Compute correlation and adjusted test statistics
[r,stats] = permucorr(x)

% Display results
disp('Correlation Matrix')
disp('r:')
disp(r)
disp('stats:')
disp(stats)