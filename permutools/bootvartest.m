function [chi2,p,ci,stats,dist] = bootvartest(x,v,varargin)
%BOOTVARTEST  One-sample bootstrapped chi-square test.
%   CHI2 = BOOTVARTEST(X,V) performs a one-sample bootstrapped test based
%   on the chi-squared statistic of the the null hypothesis that the data
%   in X come from a distribution with variance V, and returns the test
%   statistic. If X is a matrix, separate permutation tests are performed
%   along each column of X, and a vector of results is returned. V must be
%   a scalar.
%
%   BOOTVARTEST treats NaNs as missing values, and ignores them.
%
%   [CHI2,P] = BOOTVARTEST(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by bootstrapping the
%   data, no assumption is made about the shape of the distributions that
%   the data come from. P-values are automatically adjusted for multiple
%   comparisons using the max correction method.
%
%   [CHI2,P,CI] = BOOTVARTEST(...) returns a 100*(1-ALPHA)% confidence
%   interval for the true variance. CIs are also adjusted for multiple
%   comparisons using the max correction method.
%
%   [CHI2,P,CI,STATS] = BOOTVARTEST(...) returns a structure with the
%   following fields:
%       'df'        -- the degrees of freedom of each test
%       'varx'      -- the estimated population variance of X
%
%   [CHI2,P,CI,STATS,DIST] = BOOTVARTEST(...) returns the bootstrapped
%   sampling distribution of the test statistic.
%
%   [...] = BOOTVARTEST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%       'nboot'     An integer scalar specifying the number of bootstraps
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
%   See also VARTEST BOOTVARTEST2 BOOTEFFECTSIZE.
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

% Compute sum of squares
mu = sum(x,nanflag)./nobs;
xdm = x-mu;
sumsq = sum(xdm.^2,nanflag);

% Compute test statistic
if v > 0
    chi2 = sumsq./v;
else
    chi2 = Inf(size(1,nvar),'like',sumsq);
    chi2(sumsq==0) = NaN;
end

if nargout > 1

    % Compute sample variance using fast algo and scale to V
    varx = (sum(x.^2,nanflag)-(sum(x,nanflag).^2)./nobs)./df;
    xnull = xdm.*sqrt(v./varx);

    % Generate random bootstraps
    rng(arg.seed);

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
        [~,imax] = max(dist,[],2);
        [~,imin] = min(dist,[],2);
        csvar = [0;cumsum(ones(arg.nboot-1,1)*nvar)];
        dist = dist';
        pdmax = dist(imax+csvar);
        pdmin = dist(imin+csvar);
        k = 1;
    else
        pdmax = dist;
        pdmin = dist;
        k = 2;
    end

    % Compute p-value & CI
    switch arg.tail
        case 'both'
            p = k*(min(sum(chi2<=pdmax),...
                sum(chi2>=pdmin))+1)/(arg.nboot+1);
            if nargout > 2
                crit = [prctile(pdmin,100*arg.alpha/2);...
                    prctile(pdmax,100*(1-arg.alpha/2))];
                ci = sumsq./crit;
            end
        case 'right'
            p = (sum(chi2<=pdmax)+1)/(arg.nboot+1);
            if nargout > 2
                crit = prctile(pdmax,100*(1-arg.alpha));
                ci = [sumsq./crit;Inf(1,nvar)];
            end
        case 'left'
            p = (sum(chi2>=pdmin)+1)/(arg.nboot+1);
            if nargout > 2
                crit = prctile(pdmin,100*arg.alpha);
                ci = [zeros(1,nvar);sumsq./crit];
            end
    end

end

if nargout > 3
    varx = var(x,nanflag);
end

% Arrange results in a matrix if specified
if arg.mat
    chi2 = ptvec2mat(chi2);
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
        varx = ptvec2mat(varx);
    end
end

% Store statistics in a structure
if nargout > 3
    stats.df = df;
    stats.varx = varx;
end