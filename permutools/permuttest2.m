function [h,p,ci,stats] = permuttest2(x,y,varargin)
%PERMUTTEST2  Unpaired two-sample permutation-based t-test.
%   H = PERMUTTEST2(X,Y) returns the results of a two-sample permutation
%   test between X and Y based on the t-statistic. H=0 indicates that the
%   null hypothesis that X and Y have equal means cannot be rejected at the
%   5% significance level, wheras H=1 indicates that the null hypothesis
%   can be rejected at the 5% significance level. As the null distribution
%   is generated empirically by permuting the data, no assumption is made
%   about the shape of the distribution that the data come from, except
%   that the variance is equal. X and Y can have different lengths.
%
%   If X and Y are matrices, multiple permutation tests are performed
%   simultaneously between each corresponding pair of columns in X and Y,
%   and a vector of results is returned. Family-wise error rate (FWER) is
%   controlled for multiple permutation tests using the maximum statistic
%   correction method (Blair et al., 1994). This method provides strong
%   control of FWER, even for small sample sizes, and is much more powerful
%   than traditional correction methods (Groppe et al., 2011). It is also
%   rather insensitive to differences in population variance when samples
%   of equal size are used (Groppe et al., 2011b). For samples of unequal
%   size or variance, Welch's t-statistic may be used as it is less
%   sensitive to differences in variance (but also less sensitive to
%   differences in means).
%
%   If Y is empty, permutation tests between every pair of columns in X are
%   performed and a matrix of results is returned.
%
%   PERMUTTEST2 treats NaNs as missing values, and ignores them.
%
%   [H,P] = PERMUTTEST2(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%
%   [H,P,CI] = PERMUTTEST2(...) returns a 100*(1-ALPHA)% confidence
%   interval for the true difference of population means (Groppe, 2016).
%
%   [H,P,CI,STATS] = PERMUTTEST2(...) returns a structure with the
%   following fields:
%       'tstat'     -- the value of the test statistic
%       'df'        -- the degrees of freedom of each test
%       'sd'    	-- the estimated population standard deviation of X, or
%                      of X-Y for a paired test
%
%   [...] = PERMUTTEST2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%       'vartype'   A string specifying the variance equivalence of X and Y
%                   to determine the SD and t-statistic estimation method:
%                       'equal'   	assume equal variances (default)
%                       'unequal' 	assume unequal variances
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
%   See also TTEST2 PERMUTTEST PERMUVARTEST2 BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Groppe DM, Urbach TP, Kutas M (2011a) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.
%       [3] Groppe DM, Urbach TP, Kutas M (2011b) Mass univariate analysis
%           of event-related brain potentials/fields II: Simulation studies
%           Psychophysiology, 48(12):1726-1737.
%       [4] Groppe DM (2016) Combating the scientific decline effect with
%           confidence (intervals). Psychophysiology, 54(1):139-145.

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

% Use only rows with no NaNs if specified
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
dfx = nobsx-1;
dfy = nobsy-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:))) || any(isnan(y(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute sample variance using fast algo
smx = sum(x,nanflag);
smy = sum(y,nanflag);
varx = (sum(x.^2,nanflag)-(smx.^2)./nobsx)./dfx;
vary = (sum(y.^2,nanflag)-(smy.^2)./nobsy)./dfy;

% Concatenate samples
x = [x;y];
[maxnobs,nvar] = size(x);
nobs = sum(~isnan(x));
sqrtn = sqrt(nobs./(nobsx.*nobsy));

% Compute test statistic
switch arg.vartype
    case 'equal'
        df = nobs-2;
        sd = sqrt((dfx.*varx+dfy.*vary)./df);
        se = sd.*sqrtn;
    case 'unequal'
        se2x = varx./nobsx;
        se2y = vary./nobsy;
        df = (se2x+se2y).^2./(se2x.^2./dfx+se2y.^2./dfy);
        sd = sqrt([varx;vary]);
        se = sqrt(se2x+se2y);
end
mx = smx./nobsx-smy./nobsy;
tstat = mx./se;

% Generate permutation distribution
rng(arg.seed);
[~,idx] = sort(rand(maxnobs,arg.nperm));
ix = idx(1:maxnobsx,:);
iy = idx(maxnobsx+1:maxnobs,:);
pd = zeros(arg.nperm,nvar);
for i = 1:arg.nperm
    xp = x(ix(:,i),:);
    yp = x(iy(:,i),:);
    smx = sum(xp,nanflag);
    smy = sum(yp,nanflag);
    s2x = (sum(xp.^2,nanflag)-(smx.^2)./nobsx)./dfx;
    s2y = (sum(yp.^2,nanflag)-(smy.^2)./nobsy)./dfy;
    switch arg.vartype
        case 'equal'
            sep = sqrt((dfx.*s2x+dfy.*s2y)./df).*sqrtn;
        case 'unequal'
            sep = sqrt(s2x./nobsx+s2y./nobsy);
    end
    pd(i,:) = (smx./nobsx-smy./nobsy)./sep;
end

% Apply max correction if specified
if arg.correct
    [~,idx] = max(abs(pd),[],2);
    csvar = [0;cumsum(ones(arg.nperm-1,1)*nvar)];
    pd = pd';
    pd = pd(idx+csvar);
end

% Compute test statistics
switch arg.tail
    case 'both'
        pd = abs(pd);
        p = (sum(abs(tstat)<=pd)+1)/(arg.nperm+1);
        crit = prctile(pd,100*(1-arg.alpha)).*se;
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