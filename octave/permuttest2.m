function [t,p,ci,stats,dist] = permuttest2(x,y,varargin)
%PERMUTTEST2  Unpaired two-sample permutation-based t-test.
%   T = PERMUTTEST2(X,Y) performs a two-sample permutation test based on
%   the t-statistic of the hypothesis that the data in X and Y come from
%   distributions with equal means, and returns the test statistic. If X
%   and Y are matrices, separate permutation tests are performed between
%   each corresponding pair of columns in X and Y, and a vector of results
%   is returned. If Y is empty, two-tailed permutation tests between every
%   pair of columns in X are performed, and a matrix of results is
%   returned. X and Y can have different lengths.
%
%   For samples of unequal size or variance, Welch's t-statistic may be
%   used by setting the 'vartype' parameter to 'unequal' as it is less
%   sensitive to differences in variance (but also less sensitive to
%   differences in means).
%
%   PERMUTTEST2 treats NaNs as missing values, and ignores them.
%
%   [T,P] = PERMUTTEST2(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data, no assumption is made about the shape of the distributions that
%   the data come from, except that they have equal variances. P-values are
%   automatically adjusted for multiple comparisons using the max
%   correction method.
%
%   [T,P,CI] = PERMUTTEST2(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for the true difference of population means. CIs are also
%   adjusted for multiple comparisons using the max correction method.
%
%   [T,P,CI,STATS] = PERMUTTEST2(...) returns a structure with the
%   following fields:
%       'df'        -- the degrees of freedom of each test
%       'sd'        -- the pooled estimate of the population standard
%                      deviation (equal variances) or a vector containing
%                      the unpooled estimates of the population standard
%                      deviations (unequal variance)
%       'mu'        -- the estimated population mean of X-Y
%
%   [T,P,CI,STATS,DIST] = PERMUTTEST2(...) returns the permuted sampling
%   distribution of the test statistic.
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
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%       [2] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.
%       [4] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields II: Simulation
%           studies. Psychophysiology, 48(12):1726-1737.
%       [5] Groppe DM (2016) Combating the scientific decline effect with
%           confidence (intervals). Psychophysiology, 54(1):139-145.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

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

% Set up comparison
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
[maxnobsx,nvar] = size(x);
nobsx = sum(~isnan(x));
nobsy = sum(~isnan(y));

% Compute degrees of freedom
dfx = nobsx-1;
dfy = nobsy-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:))) || any(isnan(y(:)))
    nanflag = true;
else
    nanflag = false;
end

% Compute sample variance using fast algo
if nanflag
    smx = nansum(x);
    smy = nansum(y);
    varx = (nansum(x.^2)-(smx.^2)./nobsx)./dfx;
    vary = (nansum(y.^2)-(smy.^2)./nobsy)./dfy;
else
    smx = sum(x);
    smy = sum(y);
    varx = (sum(x.^2)-(smx.^2)./nobsx)./dfx;
    vary = (sum(y.^2)-(smy.^2)./nobsy)./dfy;
end

% Concatenate samples
x = [x;y];
nobs = sum(~isnan(x));
sqrtn = sqrt(nobs./(nobsx.*nobsy));

% Compute standard error
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

% Compute mean difference
mu = smx./nobsx-smy./nobsy;

% Compute test statistic
t = mu./se;

if nargout > 1

    % Generate random permutations
    rng(arg.seed);
    maxnobs = size(x,1);
    [~,idx] = sort(rand(maxnobs,arg.nperm));
    i1 = idx(1:maxnobsx,:);
    i2 = idx(maxnobsx+1:maxnobs,:);

    % Estimate sampling distribution
    dist = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        x1 = x(i1(:,i),:);
        x2 = x(i2(:,i),:);
        if nanflag
            sm1 = nansum(x1);
            sm2 = nansum(x2);
            var1 = (nansum(x1.^2)-(sm1.^2)./nobsx)./dfx;
            var2 = (nansum(x2.^2)-(sm2.^2)./nobsy)./dfy;
        else
            sm1 = sum(x1);
            sm2 = sum(x2);
            var1 = (sum(x1.^2)-(sm1.^2)./nobsx)./dfx;
            var2 = (sum(x2.^2)-(sm2.^2)./nobsy)./dfy;
        end
        switch arg.vartype
            case 'equal'
                sep = sqrt((dfx.*var1+dfy.*var2)./df).*sqrtn;
            case 'unequal'
                sep = sqrt(var1./nobsx+var2./nobsy);
        end
        dist(i,:) = (sm1./nobsx-sm2./nobsy)./sep;
    end

    % Apply max correction if specified
    if arg.correct
        [~,idx] = max(abs(dist),[],2);
        csvar = [0;cumsum(ones(arg.nperm-1,1)*nvar)];
        dist = dist';
        dist = dist(idx+csvar);
    end

    % Compute p-value & CI
    switch arg.tail
        case 'both'
            pdabs = abs(dist);
            p = (sum(abs(t)<=pdabs)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(pdabs,100*(1-arg.alpha)).*se;
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