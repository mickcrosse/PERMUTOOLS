function [chisq,p,ci,stats,dist] = permuvartest(x,v,varargin)
%PERMUVARTEST  Bootstrap-based one-sample test of variance.
%   CHISQ = PERMUVARTEST(X,V) performs a one-sample bootstrap test based on
%   the Chi-squared statistic of the the null hypothesis that the data in X
%   come from a distribution with variance V, and returns the test
%   statistic. If X is a matrix, separate bootstrap tests are performed
%   along each column of X, and a vector of results is returned. V must be
%   a scalar.
%
%   PERMUVARTEST treats NaNs as missing values, and ignores them.
%
%   [CHISQ,P] = PERMUVARTEST(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by bootstrapping the
%   data, no assumption is made about the shape of the distributions that
%   the data come from. P-values are automatically adjusted for multiple
%   comparisons using the max correction method.
%
%   [CHISQ,P,CI] = PERMUVARTEST(...) returns a 100*(1-ALPHA)% confidence
%   interval for the true variance. CIs are also adjusted for multiple
%   comparisons using the max correction method.
%
%   [CHISQ,P,CI,STATS] = PERMUVARTEST(...) returns a structure with the
%   following fields:
%       'chisqstat' -- the value of the test statistic
%       'df'        -- the degrees of freedom of each test
%       'var'       -- the estimated population variance
%       'method'    -- the statistical method used
%
%   [CHISQ,P,CI,STATS,DIST] = PERMUVARTEST(...) returns the bootstrapped
%   sampling distribution of the test statistic.
%
%   [...] = PERMUVARTEST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%                       'both'      variance is not V (default)
%                       'right'     variance is greater than V
%                       'left'      variance is less than V
%       'nboot'     An integer scalar specifying the number of bootstraps
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   control FWER using max correction (default=1).
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the bootstrap generator (default=shuffle).
%
%   See also VARTEST PERMUVARTEST2 PERMUEFFECTSIZE.
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

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

narginchk(2,Inf);

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

% Compute degrees of freedom
df = nobs-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitnan';
else
    nanflag = 'includenan';
end

% Compute observed statistics
mu = sum(x,nanflag)./nobs;
xdm = x-mu;
sumsq = sum(xdm.^2,nanflag);
if v > 0
    chisq = sumsq./v;
else
    chisq = Inf(size(1,nvar),'like',sumsq);
    chisq(sumsq==0) = NaN;
end

if nargout > 1

    rng(arg.seed);

    % Compute sample variance using fast algo and scale to V
    varx = (sum(x.^2,nanflag)-(sum(x,nanflag).^2)./nobs)./df;
    xnull = xdm.*sqrt(v./varx);

    % Estimate sampling distribution
    dist = zeros(arg.nboot,nvar);
    try
        idx = randi(maxnobs,maxnobs,arg.nboot,'uint32');
        for i = 1:nvar
            xcol = xnull(:,i);
            xp = xcol(idx);
            nobsp = sum(~isnan(xp),1);
            mup = sum(xp,1,nanflag)./nobsp;
            sumsqp = sum((xp-mup).^2,1,nanflag);
            dist(:,i) = sumsqp'./v;
        end
    catch
        for i = 1:arg.nboot
            idx = randi(maxnobs,maxnobs,1,'uint32');
            xp = xnull(idx,:);
            nobsp = sum(~isnan(xp),1);
            mup = sum(xp,1,nanflag)./nobsp;
            sumsqp = sum((xp-mup).^2,1,nanflag);
            dist(i,:) = sumsqp./v;
        end
    end

    % Apply max correction if specified
    if arg.correct
        pdmax = max(dist,[],2);
        pdmin = min(dist,[],2);
    else
        pdmax = dist;
        pdmin = dist;
    end

    % Compute p-value
    switch arg.tail
        case {'both','two'}
            p = min(1,2*(min(sum(chisq<=pdmax),...
                sum(chisq>=pdmin))+1)/(arg.nboot+1));
        case 'right'
            p = (sum(chisq<=pdmax)+1)/(arg.nboot+1);
        case 'left'
            p = (sum(chisq>=pdmin)+1)/(arg.nboot+1);
    end

end

% Compute confidence interval
if nargout > 2
    switch arg.tail
        case {'both','two'}
            crit = [prctile(pdmax,100*(1-arg.alpha/2));...
                prctile(pdmin,100*arg.alpha/2)];
            ci = sumsq./crit;
        case 'right'
            crit = prctile(pdmax,100*(1-arg.alpha));
            ci = [sumsq./crit;Inf(1,nvar)];
        case 'left'
            crit = prctile(pdmin,100*arg.alpha);
            ci = [zeros(1,nvar);sumsq./crit];
    end
end

% Store statistics in a structure
if nargout > 3
    stats.chisqstat = chisq;
    stats.df = df;
    stats.var = varx;
    stats.method = 'chisqtest';
end