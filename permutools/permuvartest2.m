function [h,p,ci,stats] = permuvartest2(x,y,varargin)
%PERMUVARTEST2  Unpaired two-sample permutation-based F-test.
%   H = PERMUVARTEST2(X,Y) returns the results of a two-sample permutation
%   test between X and Y based on the F-statistic. H=0 indicates that the
%   null hypothesis that X and Y have equal variances cannot be rejected at
%   the 5% significance level, wheras H=1 indicates that the null
%   hypothesis can be rejected at the 5% significance level. As the null
%   distribution is generated empirically by permuting the data, no
%   assumption is made about the shape of the distribution that the data
%   come from, except that the variance is equal. X and Y can have
%   different lengths.
%
%   If X and Y are matrices, multiple permutation tests are performed
%   simultaneously between each corresponding pair of columns in X and Y,
%   and a vector of results is returned. Family-wise error rate (FWER) is
%   controlled for multiple permutation tests using the maximum statistic
%   correction method (Blair et al., 1994). This method provides strong
%   control of FWER, even for small sample sizes, and is much more powerful
%   than traditional correction methods (Groppe et al., 2011).
%
%   If Y is empty, permutation tests between every pair of columns in X are
%   performed and a matrix of results is returned.
%
%   PERMUVARTEST2 treats NaNs as missing values, and ignores them.
%
%   [H,P] = PERMUVARTEST2(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%
%   [H,P,CI] = PERMUVARTEST2(...) returns a 100*(1-ALPHA)% confidence
%   interval for the true ratio of sample variances.
%
%   [H,P,CI,STATS] = PERMUVARTEST2(...) returns a structure with the
%   following fields:
%       'fstat'     -- the value of the test statistic
%       'df1'       -- the numerator degrees of freedom of each test
%       'df2'    	-- the denominator degrees of freedom of each test
%
%   [...] = PERMUVARTEST2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%                       'both'      means are not equal (default)
%                       'right'     mean of X is greater than mean of Y
%                       'left'      mean of X is less than mean of Y
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
%
%   See also VARTEST2 PERMUTTEST2 BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
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
if isempty(y)
    warning('Comparing all columns of X using two-tailed test...')
    [x,y] = ptpaircols(x);
    arg.tail = 'both';
    arg.mat = true;
end
if size(x,2)~=size(y,2)
    error('X and Y must have the same number of variables.')
end

% Use only rows with no NaN values if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(x),2),:);
        y = y(~any(isnan(y),2),:);
end

% Get data dimensions, ignoring NaNs
maxnobsx = size(x,1);
nobsx = sum(~isnan(x));
nobsy = sum(~isnan(y));

% Compute degrees of freedom
df1 = nobsx-1;
df2 = nobsy-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:))) || any(isnan(y(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute test statistic
varx = (sum(x.^2,nanflag)-(sum(x,nanflag).^2)./nobsx)./df1;
vary = (sum(y.^2,nanflag)-(sum(y,nanflag).^2)./nobsy)./df2;
fstat = varx./vary;

% Concatenate data
x = [x;y];
[maxnobs,nvar] = size(x);

% Generate permutation distribution
rng(arg.seed);
[~,idx] = sort(rand(maxnobs,arg.nperm));
ix = idx(1:maxnobsx,:);
iy = idx(maxnobsx+1:maxnobs,:);
pd = zeros(arg.nperm,nvar);
for i = 1:arg.nperm
    xp = x(ix(:,i),:);
    yp = x(iy(:,i),:);
    varx = (sum(xp.^2,nanflag)-(sum(xp,nanflag).^2)./nobsx)./df1;
    vary = (sum(yp.^2,nanflag)-(sum(yp,nanflag).^2)./nobsy)./df2;
    pd(i,:) = varx./vary;
end

% Apply max correction if specified
if arg.correct
    [~,imax] = max(pd,[],2);
    [~,imin] = min(pd,[],2);
    csvar = [0;cumsum(ones(arg.nperm-1,1)*nvar)];
    pd = pd';
    pdmax = pd(imax+csvar);
    pdmin = pd(imin+csvar);
    k = 1;
else
    pdmax = pd;
    pdmin = pd;
    k = 2;
end

% Compute corrected test statistics using max statistic correction
switch arg.tail
    case 'both'
        p = k*(min(sum(fstat<=pdmax),sum(fstat>=pdmin))+1)/(arg.nperm+1);
        if nargout > 2
            crit = [prctile(pdmin,100*arg.alpha/2);...
                prctile(pdmax,100*(1-arg.alpha/2))];
            ci = fstat./crit;
        end
    case 'right'
        p = (sum(fstat<=pdmax)+1)/(arg.nperm+1);
        if nargout > 2
            crit = prctile(pdmax,100*(1-arg.alpha));
            ci = [fstat./crit;Inf(1,nvar)];
        end
    case 'left'
        p = (sum(fstat>=pdmin)+1)/(arg.nperm+1);
        if nargout > 2
            crit = prctile(pdmin,100*arg.alpha);
            ci = [zeros(1,nvar);fstat./crit];
        end
end

% Determine if p-values exceed alpha level
h = cast(p<=arg.alpha,'like',p);
h(isnan(p)) = NaN;

% Arrange test results in a matrix if specified
if arg.mat
    h = ptvec2mat(h);
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
        fstat = ptvec2mat(fstat);
        df1 = ptvec2mat(df1);
        df2 = ptvec2mat(df2);
    end
end

% Store test statistics in a structure
if nargout > 3
    stats = struct('fstat',fstat,'df1',df1,'df2',df2);
end