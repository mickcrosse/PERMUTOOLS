function [stat,p,ci,stats,dist] = permuvartest2(x,y,varargin)
%PERMUVARTEST2  Permutation-based two-sample test of variance and Conover squared ranks test.
%   STAT = PERMUVARTEST2(X,Y) performs a two-sample permutation test based
%   on the F-statistic of the null hypothesis that the data in X and Y come
%   from distributions with equal variances, and returns the test
%   statistic. X and Y can have different lengths. If X and Y are matrices,
%   separate permutation tests are performed between each corresponding
%   pair of columns in X and Y, and a vector of results is returned. If Y
%   is empty, two-tailed permutation tests between every pair of columns in
%   X are performed, and a matrix of results is returned.
%
%   For non-normally distributed samples, the raw data may be transformed
%   to squared ranks in order to compute a Conover squared ranks test by
%   setting the 'type' parameter to 'squarerank' or 'rank'.
%
%   PERMUVARTEST2 treats NaNs as missing values, and ignores them.
%
%   [STAT,P] = PERMUVARTEST2(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by permuting the data,
%   no assumption is made about the shape of the distributions that the
%   data come from. P-values are automatically adjusted for multiple
%   comparisons using the max correction method.
%
%   [STAT,P,CI] = PERMUVARTEST2(...) returns a 100*(1-ALPHA)% confidence
%   interval for the true ratio of population variances. CIs are also
%   adjusted for multiple comparisons using the max correction method.
%
%   [STAT,P,CI,STATS] = PERMUVARTEST2(...) returns a structure with the
%   following fields:
%       'fstat'     -- the value of the test statistic
%       'df1'       -- the numerator degrees of freedom of each test
%       'df2'       -- the denominator degrees of freedom of each test
%       'var'       -- the estimated population variance
%       'method'    -- the statistical method used
%
%   For rank-based tests ('type','squarerank'), STATS contains:
%       'tstat'     -- the value of the test statistic
%       'median'    -- the estimated population median
%       'iqr'       -- the estimated population interquartile range
%       'method'    -- the statistical method used
%
%   [STAT,P,CI,STATS,DIST] = PERMUVARTEST2(...) returns the permuted
%   sampling distribution of the test statistic.
%
%   [...] = PERMUVARTEST2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of test to perform:
%                       'ftest'         two-sample F-test of variance for
%                                       continuous unpaired data (default)
%                       'squarerank'    Conover squared ranks test of scale
%                                       for continuous heavy-tailed data
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'  variances are not equal (default)
%                       'right' variance of X is greater than variance of Y
%                       'left'  variance of X is less than variance of Y
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   control FWER using max correction (default=1).
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'matrix'    A numeric scalar (0,1) or logical indicating whether to
%                   return results as a matrix (default=0).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator (default=shuffle).
%
%   See also VARTEST2 PERMUTTEST2 BOOTEFFECTSIZE.
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

if nargin < 2
    y = [];
end

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'ftest';
elseif strcmpi(arg.type,'rank')
    arg.type = 'squarerank';
end

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
    warning('Comparing all columns of X using two-tailed tests...')
    [x,y] = ptpaircols(x);
    arg.tail = 'both';
    arg.matrix = true;
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
[maxnobsx,nvar] = size(x);
nobsx = sum(~isnan(x),1);
nobsy = sum(~isnan(y),1);
df1 = nobsx-1;
df2 = nobsy-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:))) || any(isnan(y(:)))
    nanflag = 'omitnan';
else
    nanflag = 'includenan';
end

% Compute observed statistics
switch arg.type
    case 'ftest'
        mux = sum(x,1,nanflag)./nobsx;
        muy = sum(y,1,nanflag)./nobsy;
        varx = (sum(x.^2,1,nanflag)-(sum(x,1,nanflag).^2)./nobsx)./df1;
        vary = (sum(y.^2,1,nanflag)-(sum(y,1,nanflag).^2)./nobsy)./df2;
        stat = varx./vary;
    case 'squarerank'
        if nargout > 3
            iqrxy = [iqr(x);iqr(y)];
        end
        medx = median(x,1,nanflag);
        medy = median(y,1,nanflag);
        devxy = tiedrank(abs([x-medx;y-medy]));
        x = devxy(1:size(x,1),:).^2;
        y = devxy(size(x,1)+1:end,:).^2;
        mux = sum(x,1,nanflag)./nobsx;
        muy = sum(y,1,nanflag)./nobsy;
        varx = (sum(x.^2,1,nanflag)-(sum(x,1,nanflag).^2)./nobsx)./df1;
        vary = (sum(y.^2,1,nanflag)-(sum(y,1,nanflag).^2)./nobsy)./df2;
        sp2 = (df1.*varx+df2.*vary)./(nobsx + nobsy - 2);
        se = sqrt(sp2.*(1./nobsx+1./nobsy));
        stat = (mux-muy)./se;
    otherwise
        error('The TYPE parameter value must be FTEST or SQUARERANK.')
end

if nargout > 1

    rng(arg.seed);

    % Demean and concatenate samples
    x = x-mux;
    y = y-muy;
    x = [x;y];

    % Generate random permutations
    maxnobs = size(x,1);
    [~,idx] = sort(rand(maxnobs,arg.nperm));
    idx1 = idx(1:maxnobsx,:);
    idx2 = idx(maxnobsx+1:maxnobs,:);

    % Estimate sampling distribution
    dist = zeros(arg.nperm,nvar);
    switch nanflag
        case 'omitnan'
            % Dynamic N-tracking for missing data
            switch arg.type
                case 'ftest'
                    for i = 1:arg.nperm
                        x1 = x(idx1(:,i),:);
                        x2 = x(idx2(:,i),:);
                        nobs1 = sum(~isnan(x1),1);
                        nobs2 = sum(~isnan(x2),1);
                        var1 = (sum(x1.^2,1,nanflag)-...
                            (sum(x1,1,nanflag).^2)./nobs1)./(nobs1-1);
                        var2 = (sum(x2.^2,1,nanflag)-...
                            (sum(x2,1,nanflag).^2)./nobs2)./(nobs2-1);
                        dist(i,:) = var1./var2;
                    end
                case 'squarerank'
                    for i = 1:arg.nperm
                        x1 = x(idx1(:,i),:);
                        x2 = x(idx2(:,i),:);
                        nobs1 = sum(~isnan(x1),1);
                        nobs2 = sum(~isnan(x2),1);
                        df1p = nobs1-1;
                        df2p = nobs2-1;
                        mu1 = sum(x1,1,nanflag)./nobs1;
                        mu2 = sum(x2,1,nanflag)./nobs2;
                        var1 = (sum(x1.^2,1,nanflag)-...
                            (sum(x1,1,nanflag).^2)./nobs1)./df1p;
                        var2 = (sum(x2.^2,1,nanflag)-...
                            (sum(x2,1,nanflag).^2)./nobs2)./df2p;
                        sp2 = (df1p.*var1+df2p.*var2)./(nobs1+nobs2-2);
                        se = sqrt(sp2.*(1./nobs1+1./nobs2));
                        dist(i,:) = (mu1-mu2)./se;
                    end
            end
        case 'includenan'
            % Fast vectorized calculation for complete data
            I = zeros(maxnobs,arg.nperm);
            colidx = repmat(1:arg.nperm,maxnobsx,1);
            linidx = sub2ind([maxnobs,arg.nperm],idx1,colidx);
            I(linidx) = 1;
            sum1 = I'*x;
            sum2 = sum(x,1)-sum1;
            sumsq1 = I'*(x.^2);
            sumsq2 = sum(x.^2,1)-sumsq1;
            var1 = (sumsq1-(sum1.^2)./nobsx)./df1;
            var2 = (sumsq2-(sum2.^2)./nobsy)./df2;
            switch arg.type
                case 'ftest'
                    dist = var1./var2;
                case 'squarerank'
                    mu1 = sum1./nobsx;
                    mu2 = sum2./nobsy;
                    sp2 = (df1.*var1+df2.*var2)./(nobsx+nobsy-2);
                    se = sqrt(sp2.*(1./nobsx+1./nobsy));
                    dist = (mu1-mu2)./se;
            end
    end

    % Apply max correction if specified
    if arg.correct
        distmax = max(dist,[],2);
        distmin = min(dist,[],2);
    else
        distmax = dist;
        distmin = dist;
    end

    % Compute p-value
    switch arg.tail
        case 'both'
            p = min(1,2*(min(sum(stat>=distmin,1),...
                sum(stat<=distmax,1))+1)/(arg.nperm+1));
        case 'right'
            p = (sum(stat<=distmax,1)+1)/(arg.nperm+1);
        case 'left'
            p = (sum(stat>=distmin,1)+1)/(arg.nperm+1);
    end

end

% Compute confidence interval
if nargout > 2
    switch arg.type
        case 'ftest'
            switch arg.tail
                case 'both'
                    crit = [prctile(distmin,100*arg.alpha/2);...
                        prctile(distmax,100*(1-arg.alpha/2))];
                    ci = stat./crit;
                case 'right'
                    crit = prctile(distmax,100*(1-arg.alpha));
                    ci = [stat./crit;Inf(1,nvar)];
                case 'left'
                    crit = prctile(distmin,100*arg.alpha);
                    ci = [zeros(1,nvar);stat./crit];
            end
        case 'squarerank'
            ci = NaN(2,nvar);
    end
end

% Format descriptive statistics
if nargout > 3
    switch arg.type
        case 'ftest'
            stats.fstat = stat;
            stats.df1 = df1;
            stats.df2 = df2;
            stats.var = [varx;vary];
        case 'squarerank'
            stats.tstat = stat;
            stats.median = [medx;medy];
            stats.iqr = iqrxy;
    end
    stats.method = arg.type;
end

% Arrange results in a matrix if specified
if arg.matrix
    stat = ptvec2mat(stat);
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
        switch arg.type
            case 'ftest'
                stats.fstat = stat;
                stats.df1 = ptvec2mat(df1);
                stats.df2 = ptvec2mat(df2);
                stats.var = ptvec2mat([varx;vary]);
            case 'squarerank'
                stats.tstat = stat;
                stats.median = ptvec2mat([medx;medy]);
                stats.iqr = ptvec2mat(iqrxy);
        end
    end
end