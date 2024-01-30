function [t,p,ci,stats,dist] = permuttest(x,m,varargin)
%PERMUTTEST  One-sample and paired-sample permutation-based t-test.
%   T = PERMUTTEST(X) performs a one-sample permutation test based on the
%   t-statistic of the hypothesis that the data in X come from a
%   distribution with mean zero, and returns the test statistic. If X is a
%   matrix, separate permutation tests are performed along each column of
%   X, and a vector of results is returned. If the 'compare' parameter is
%   set to 'pairwise', two-tailed permutation tests between every pair of
%   columns in X are performed, and a matrix of results is returned.
%
%   PERMUTTEST treats NaNs as missing values, and ignores them.
%
%   T = PERMUTTEST(X,M) returns the results of a one-sample permutation
%   test of the hypothesis that the data in X come from a distribution with
%   mean M. M must be a scalar.
%
%   T = PERMUTTEST(X,Y) returns the results of a paired-sample permutation
%   test between dependent samples X and Y of the hypothesis that the data
%   in X and Y come from distributions with equal means. X and Y must have
%   the same length. If X and Y are matrices, separate permutation tests
%   are performed between each corresponding pair of columns in X and Y,
%   and a vector of results is returned.
%
%   [T,P] = PERMUTTEST(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by permuting the data,
%   no assumption is made about the shape of the distribution that the data
%   come from. P-values are automatically adjusted for multiple comparisons
%   using the max correction method.
%
%   [T,P,CI] = PERMUTTEST(...) returns a 100*(1-ALPHA)% confidence interval
%   (CI) for the true mean of X, or of X-Y for a paired test. CIs are also
%   adjusted for multiple comparisons using the max correction method.
%
%   [T,P,CI,STATS] = PERMUTTEST(...) returns a structure with the following
%   fields:
%       'df'        -- the degrees of freedom of each test
%       'sd'        -- the estimated population standard deviation of X, or
%                      of X-Y for a paired test
%       'mu'        -- the estimated population mean of X, or of X-Y for a
%                      paired test
%
%   [T,P,CI,STATS,DIST] = PERMUTTEST(...) returns the permuted sampling
%   distribution of the test statistic.
%
%   [...] = PERMUTTEST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      mean is not M (two-tailed, default)
%                       'right'     mean is greater than M (right-tailed)
%                       'left'      mean is less than M (left-tailed)
%       'compare'   A string specifying what to compare each variable to
%                   when only X is entered:
%                       'zero'      compare each column of X to zero and
%                                   return a vector of results (default)
%                       'pairwise'  compare every pair of columns in X to
%                                   each other using two-tailed tests and
%                                   return a matrix of results
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
%   See also TTEST PERMUTTEST2 BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%       [2] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [3] Gondan M (2010) A permutation test for the race model
%           inequality. Behav Res Methods, 42(1):23-28.
%       [4] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(m)
    y = [];
elseif isscalar(m)
    y = repmat(m,size(x));
elseif ismatrix(m)
    y = m;
end

% Parse input arguments
arg = ptparsevarargin(varargin);

% Validate input parameters
ptvalidateparamin(x,y,arg)

% Orient data column-wise
if arg.dim==2 || isrow(x)
    x = x';
end
if ~isempty(y) && (arg.dim==2 || isrow(y))
    y = y';
end

% Set up comparison
if isempty(y)
    switch arg.compare
        case 'zero'
            y = zeros(size(x));
        case 'pairwise'
            warning('Comparing all columns of X using two-tailed test...')
            [x,y] = ptpaircols(x);
            arg.tail = 'both';
            arg.mat = true;
    end
else
    switch arg.compare
        case 'pairwise'
            error('The PAIRWISE option can only be used with one sample.')
    end
end
if size(x)~=size(y)
    error('X and Y must be the same size.')
end

% Compute difference between samples
x = x-y;

% Use only rows with no NaN values if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(x),2),:);
end

% Get data dimensions, ignoring NaNs
[maxnobs,nvar] = size(x);
nobs = sum(~isnan(x));

% Compute degrees of freedom
df = nobs-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute standard deviation
sd = std(x,nanflag);

% Compute mean difference
mu = sum(x,nanflag)./nobs;

% Compute test statistic
se = sd./sqrt(nobs);
t = mu./se;

if nargout > 1

    % Generate random permutations
    rng(arg.seed);
    signx = sign(rand(maxnobs,arg.nperm)-0.5);

    % Estimate sampling distribution
    sqrtn = sqrt(nobs.*df);
    dist = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        xp = x.*repmat(signx(:,i),1,nvar);
        smx = sum(xp,nanflag);
        dist(i,:) = smx./nobs./(sqrt(sum(xp.^2)-(smx.^2)./nobs)./sqrtn);
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
            p = 2*(sum(abs(t)<=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha/2)).*se;
                ci = [mu-crit;mu+crit];
            end
        case 'right'
            p = (sum(t<=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha)).*se;
                ci = [mu-crit;Inf(1,nvar)];
            end
        case 'left'
            p = (sum(t>=dist)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(dist,100*(1-arg.alpha)).*se;
                ci = [-Inf(1,nvar);mu+crit];
            end
    end

end

% Arrange results in a matrix if specified
if arg.mat
    t = ptvec2mat(t);
    if nargout > 1
        p = ptvec2mat(p);
    end
    if nargout > 2
        ciLwr = ptvec2mat(ci(1,:));
        ciUpr = ptvec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
    end
    if nargout > 3
        df = ptvec2mat(df);
        sd = ptvec2mat(sd);
        mu = ptvec2mat(mu);
    end
end

% Store statistics in a structure
if nargout > 3
    stats.df = df;
    stats.sd = sd;
    stats.mu = mu;
end
