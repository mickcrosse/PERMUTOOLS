function [h,p,ci,stats] = permuttest(x,y,varargin)
%PERMUTTEST  One-sample and paired-sample permutation-based t-test.
%   H = PERMUTTEST(X,Y) returns the results of a paired-sample permutation
%   test between X and Y based on the t-statistic. H=0 indicates that the
%   null hypothesis that X and Y have equal means cannot be rejected at the
%   5% significance level, wheras H=1 indicates that the null hypothesis
%   can be rejected at the 5% significance level. As the null distribution
%   is generated empirically by permuting the data, no assumption is made
%   about the shape of the distribution that the data come from. X and Y
%   must have the same length.
%
%   If X and Y are matrices, multiple permutation tests are performed
%   simultaneously between each corresponding pair of columns in X and Y,
%   and a vector of results is returned. Family-wise error rate (FWER) is
%   controlled for multiple permutation tests using the maximum statistic
%   correction method (Blair et al., 1994). This method provides strong
%   control of FWER, even for small sample sizes, and is much more powerful
%   than traditional correction methods (Gondan, 2010; Groppe et al., 2011).
%
%   For one-sample permutation tests, enter the data column-wise in X and
%   leave Y empty. Here, the null hypothesis is that X has a mean of zero.
%
%   PERMUTTEST treats NaNs as missing values, and ignores them.
%
%   [H,P] = PERMUTTEST(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%
%   [H,P,CI] = PERMUTTEST(...) returns a 100*(1-ALPHA)% confidence interval
%   for the true mean of X, or of X-Y for a paired test.
%
%   [H,P,CI,STATS] = PERMUTTEST(...) returns a structure with the following
%   fields:
%       'tstat'     -- the value of the test statistic
%       'df'        -- the degrees of freedom of each test
%       'sd'    	-- the estimated population standard deviation of X, or
%                      of X-Y for a paired test
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
%       'm'         A scalar or row vector specifying the mean of the null
%                   hypothesis for each variable (default=0).
%       'sample'    A string specifying whether to perform a one-sample
%                   test or a paired-sample test when only X is entered:
%                       'one'       compare each column of X to zero and
%                                   store the results in a vector (default)
%                       'paired'    compare each pair of columns in X and
%                                   store the results in a matrix
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000, or all possible permutations for less
%                   than 14 observations).
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
%
%   See also TTEST PERMUTTEST2 BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Gondan M (2010) A permutation test for the race model
%           inequality. Behav Res Methods, 42(1):23-28.
%       [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

if nargin<2
    y = [];
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

% Set up permutation test
switch arg.sample
    case 'paired'
        if isempty(y)
            warning('Comparing all columns of X using two-tailed test...')
            [x,y] = ptpaircols(x);
            arg.tail = 'both';
            arg.mat = true;
        else
            error('The paired-sample test option only applies to X.')
        end
end
if isempty(y)
    y = 0;
elseif size(x)~=size(y)
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
dfp = sqrt(nobs.*df);

% Remove mean of null hypothesis from data
if isscalar(arg.m)
    x = x-arg.m;
else
    x = x-repmat(arg.m,maxnobs,1);
end

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute test statistic
sd = std(x,nanflag);
mx = sum(x,nanflag)./nobs;
se = sd./sqrt(nobs);
tstat = mx./se;

% Use all possible permutations if less than 14 observations
if min(nobs) < 14
    warning('Computing all possible permutations due to small N.')
    arg.nperm = 2^min(nobs);
end

% Generate permutation distribution
rng(arg.seed);
signx = sign(rand(maxnobs,arg.nperm)-0.5);
pd = zeros(arg.nperm,nvar);
for i = 1:arg.nperm
    xp = x.*repmat(signx(:,i),1,nvar);
    sm = sum(xp,nanflag);
    pd(i,:) = sm./nobs./(sqrt(sum(xp.^2)-(sm.^2)./nobs)./dfp);
end

% Apply max correction if specified
if arg.correct
    pd = max(abs(pd),[],2);
end

% Add negative values
pd(arg.nperm+1:2*arg.nperm,:) = -pd;
arg.nperm = 2*arg.nperm;

% Compute test statistics
switch arg.tail
    case 'both'
        p = 2*(sum(abs(tstat)<=pd)+1)/(arg.nperm+1);
        crit = prctile(pd,100*(1-arg.alpha/2)).*se;
        ci = [mx-crit;mx+crit];
    case 'right'
        p = (sum(tstat<=pd)+1)/(arg.nperm+1);
        crit = prctile(pd,100*(1-arg.alpha)).*se;
        ci = [mx-crit;Inf(1,nvar)];
    case 'left'
        p = (sum(tstat>=pd)+1)/(arg.nperm+1);
        crit = prctile(pd,100*(1-arg.alpha)).*se;
        ci = [-Inf(1,nvar);mx+crit];
end

% Determine if p-values exceed alpha level
h = cast(p<=arg.alpha,'like',p);
h(isnan(p)) = NaN;

% Arrange test results in a matrix if specified
if arg.mat
    h = ptvec2mat(h);
    p = ptvec2mat(p);
    ciLwr = ptvec2mat(ci(1,:));
    ciUpr = ptvec2mat(ci(2,:));
    ci = cat(3,ciLwr,ciUpr);
    ci = permute(ci,[3,1,2]);
    tstat = ptvec2mat(tstat);
    df = ptvec2mat(df);
    sd = ptvec2mat(sd);
end

% Store test statistics in a structure
stats = struct('tstat',tstat,'df',df,'sd',sd);