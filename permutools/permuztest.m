function [z,p,ci,stats,dist] = permuztest(x,m,sigma,varargin)
%PERMUZTEST  One-sample permutation-based Z-test.
%   Z = PERMUZTEST(X,M,SIGMA) performs a one-sample permutation test based
%   on the Z-statistic of the null hypothesis that the data in X come from
%   a distribution with mean M, and returns the test statistic. M and SIGMA
%   must be scalars. If X is a matrix, separate permutation tests are
%   performed along each column of X, and a vector of results is returned.
%
%   PERMUZTEST treats NaNs as missing values, and ignores them.
%
%   [Z,P] = PERMUZTEST(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by permuting the data,
%   no assumption is made about the shape of the distribution that the data
%   come from, except that the standard deviation is SIGMA. P-values are
%   automatically adjusted for multiple comparisons using the max
%   correction method.
%
%   [Z,P,CI] = PERMUZTEST(...) returns a 100*(1-ALPHA)% confidence interval
%   (CI) for the true mean. CIs are also adjusted for multiple comparisons
%   using the max correction method.
%
%   [Z,P,CI,STATS] = PERMUZTEST(...) returns a structure with the following
%   fields:
%       'sd'        -- the estimated population standard deviation of X
%       'mu'        -- the estimated population mean of X
%
%   [Z,P,CI,STATS,DIST] = PERMUZTEST(...) returns the permuted sampling
%   distribution of the test statistic.
%
%   [...] = PERMUZTEST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows.
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      mean is not M (two-tailed, default)
%                       'right'     mean is greater than M (right-tailed)
%                       'left'      mean is less than M (left-tailed)
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether
%                   to control FWER using max correction (default=true).
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator. By default, the
%                   generator is initialised based on the current time,
%                   resulting in a different permutation on each call.
%       'verbose'   A numeric or logical specifying whether to execute in
%                   verbose mode: pass in 1 for verbose mode (default), or
%                   0 for non-verbose mode.
%
%   See also ZTEST PERMUTTEST BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%       [2] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

narginchk(3,Inf);

% Parse input arguments
arg = ptparsevarargin(varargin);

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Orient data column-wise
if arg.dim==2 || isrow(x)
    x = x';
end

% Use only rows with no NaN values if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(x),2),:);
end

% Get data dimensions, ignoring NaNs
[maxnobs,nvar] = size(x);
nobs = sum(~isnan(x));

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute mean value
mu = sum(x,nanflag)./nobs;

% Compute test statistic
se = sigma./sqrt(nobs);
z = (mu-m)./se;

if nargout > 1

    % Generate random permutations
    rng(arg.seed);
    signx = sign(rand(maxnobs,arg.nperm)-0.5);

    % Estimate sampling distribution
    diffxm = x-m;
    sen = se.*nobs;
    dist = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        xp = diffxm.*repmat(signx(:,i),1,nvar);
        smx = sum(xp,nanflag);
        dist(i,:) = smx./sen;
    end

    % Apply max correction if specified
    if arg.correct
        dist = max(abs(dist),[],2);
    end

    % Add negative values
    dist(arg.nperm+1:2*arg.nperm,:) = -dist;
    arg.nperm = 2*arg.nperm;
    if arg.verbose
        fprintf('Adding negative of values to permutation distribution.\n')
        fprintf('Number of permutations used: %d\n',arg.nperm)
    end

    % Compute p-value & CI
    switch arg.tail
        case 'both'
            p = 2*(sum(abs(z)<=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha/2)).*se;
                ci = [mu-crit;mu+crit];
            end
        case 'right'
            p = (sum(z<=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha)).*se;
                ci = [mu-crit;Inf(1,nvar)];
            end
        case 'left'
            p = (sum(z>=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha)).*se;
                ci = [-Inf(1,nvar);mu+crit];
            end
    end

end

% Store statistics in a structure
if nargout > 3
    stats.sd = std(x,nanflag);
    stats.mu = mu;
end