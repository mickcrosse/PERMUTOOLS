function [stat,p,ci,stats,dist] = permuttest2(x,y,varargin)
%PERMUTTEST2  Permutation-based two-sample t-test and Wilcoxon rank-sum test.
%   STAT = PERMUTTEST2(X,Y) performs a two-sample permutation test based on
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
%   For non-normally distributed samples, the raw data may be transformed
%   to rank orders in order to compute a Wilcoxon rank-sum / Mann-Whitney U
%   test by setting the 'type' parameter to 'ranksum' or 'rank'.
%
%   PERMUTTEST2 treats NaNs as missing values, and ignores them.
%
%   [STAT,P] = PERMUTTEST2(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data, no assumption is made about the shape of the distributions that
%   the data come from, except that they have equal variances. P-values are
%   automatically adjusted for multiple comparisons using the max
%   correction method.
%
%   [STAT,P,CI] = PERMUTTEST2(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for the true difference of population means. For rank-
%   based tests ('type','ranksum'), CI returns NaNs. CIs are also adjusted
%   for multiple comparisons using the max correction method.
%
%   [STAT,P,CI,STATS] = PERMUTTEST2(...) returns a structure with the
%   following fields for standard tests ('type','ttest2'):
%       'tstat'     -- the value of the test statistic
%       'df'        -- the degrees of freedom of each test
%       'sd'        -- the pooled estimate of the population standard
%                      deviation (equal variances) or a vector containing
%                      the unpooled estimates of the population standard
%                      deviations (unequal variance)
%       'mean'      -- the estimated population mean
%       'method'    -- the statistical method used
%
%   For rank-based tests ('type','ranksum'), STATS contains:
%       'zstat'     -- the value of the test statistic
%       'iqr'       -- the pooled estimate of the population interquartile
%                      range (equal variances) or a vector containing the
%                      unpooled estimates of the population interquartile
%                      ranges (unequal variance)
%       'median'    -- the estimated population median
%       'method'    -- the statistical method used
%
%   [STAT,P,CI,STATS,DIST] = PERMUTTEST2(...) returns the permuted sampling
%   distribution of the test statistic.
%
%   [...] = PERMUTTEST2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'ttest2'    two-sample t-test (default)
%                       'ranksum'   Wilcoxon rank-sum / Mann-Whitney U test
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
%   See also TTEST2 RANKSUM PERMUTTEST PERMUVARTEST2 BOOTEFFECTSIZE.
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

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2
    y = [];
end

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'ttest2';
elseif strcmpi(arg.type,'rank')
    arg.type = 'ranksum';
end
if strcmpi(arg.tail,'two')
    arg.tail = 'both';
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
if size(x,2) ~= size(y,2)
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
    nanflag = 'omitnan';
else
    nanflag = 'includenan';
end

% Compute variance using fast algo
switch arg.type
    case 'ttest2'
        sumx = sum(x,nanflag); sumy = sum(y,nanflag);
        varx = (sum(x.^2,nanflag)-(sumx.^2)./nobsx)./dfx;
        vary = (sum(y.^2,nanflag)-(sumy.^2)./nobsy)./dfy;
    case 'ranksum'
        if nargout > 3
            iqrx = [iqr(x);iqr(y)];
            medx = [median(x,nanflag);median(y,nanflag)];
        end
end

% Concatenate samples
x = [x;y];
nobs = sum(~isnan(x));
sqrtn = sqrt(nobs./(nobsx.*nobsy));

% Compute observed statistics
switch arg.type
    case 'ttest2'
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
        mu = sumx./nobsx-sumy./nobsy;
        stat = mu./se;
    case 'ranksum'
        df = nobs-2;
        [r,tieadj] = tiedrank(x);
        rx = r(1:maxnobsx,:);
        w = sum(rx,1,nanflag);
        meanw = (nobsx.*(nobsx+nobsy+1))./2;
        varw = (nobsx.*nobsy)./12.*((nobsx+nobsy+1)-tieadj./...
            ((nobsx+nobsy).*(nobsx+nobsy-1)));
        stat = (w-meanw)./sqrt(varw);
    otherwise
        error('The TYPE parameter value must be TTEST2, or RANKSUM.')
end

if nargout > 1

    rng(arg.seed);

    % Generate random permutations
    maxnobs = size(x,1);
    [~,idx] = sort(rand(maxnobs,arg.nperm));
    i1 = idx(1:maxnobsx,:);
    i2 = idx(maxnobsx+1:maxnobs,:);

    % Estimate sampling distribution
    switch arg.type
        case 'ttest2'
            switch nanflag
                case 'omitnan'
                    % Dynamic N-tracking for missing data
                    dist = zeros(arg.nperm,nvar);
                    for i = 1:arg.nperm
                        x1 = x(i1(:,i),:);
                        x2 = x(i2(:,i),:);
                        nobs1 = sum(~isnan(x1));
                        nobs2 = sum(~isnan(x2));
                        df1 = nobs1-1;
                        df2 = nobs2-1;
                        sum1 = sum(x1,nanflag);
                        sum2 = sum(x2,nanflag);
                        var1 = (sum(x1.^2,nanflag)-(sum1.^2)./nobs1)./df1;
                        var2 = (sum(x2.^2,nanflag)-(sum2.^2)./nobs2)./df2;
                        switch arg.vartype
                            case 'equal'
                                sep = sqrt((df1.*var1+df2.*var2)./...
                                    (nobs1+nobs2-2)).*...
                                    sqrt((nobs1+nobs2)./(nobs1.*nobs2));
                            case 'unequal'
                                sep = sqrt(var1./nobs1+var2./nobs2);
                        end
                        dist(i,:) = (sum1./nobs1-sum2./nobs2)./sep;
                    end
                case 'includenan'
                    % Fast vectorized calculation for complete data
                    I = zeros(maxnobs,arg.nperm);
                    colidx = repmat(1:arg.nperm,maxnobsx,1);
                    linidx = sub2ind([maxnobs,arg.nperm],i1,colidx);
                    I(linidx) = 1;
                    sum1 = I'*x;
                    sum2 = sum(x,1)-sum1;
                    sumsq1 = I'*(x.^2);
                    sumsq2 = sum(x.^2,1)-sumsq1;
                    var1 = (sumsq1-(sum1.^2)./nobsx)./dfx;
                    var2 = (sumsq2-(sum2.^2)./nobsy)./dfy;
                    switch arg.vartype
                        case 'equal'
                            sep = sqrt((dfx.*var1+dfy.*var2)./df).*sqrtn;
                        case 'unequal'
                            sep = sqrt(var1./nobsx+var2./nobsy);
                    end
                    dist = (sum1./nobsx-sum2./nobsy)./sep;
            end
        case 'ranksum'
            switch nanflag
                case 'omitnan'
                    % Dynamic N-tracking for missing data
                    dist = zeros(arg.nperm,nvar);
                    for i = 1:arg.nperm
                        r1 = r(i1(:,i),:);
                        wp = sum(r1,1,nanflag);
                        nobs1p = sum(~isnan(r1));
                        nobs2p = nobs-nobs1p;
                        meanwp = nobs1p.*(nobs+1)./2;
                        varwp = (nobs1p.*nobs2p)./12.*((nobs+1)-tieadj./...
                            (nobs.*(nobs-1)));
                        dist(i,:) = (wp-meanwp)./sqrt(varwp);
                    end
                case 'includenan'
                    % Fast vectorized calculation for complete data
                    I = zeros(maxnobs,arg.nperm);
                    colidx = repmat(1:arg.nperm,maxnobsx,1);
                    linidx = sub2ind([maxnobs,arg.nperm],i1,colidx);
                    I(linidx) = 1;
                    wp = I'*r;
                    dist = (wp-meanw)./sqrt(varw);
            end
    end

    % Apply max correction if specified
    if arg.correct
        switch arg.tail
            case 'both'
                dist = max(abs(dist),[],2);
            case 'right'
                dist = max(dist,[],2);
            case 'left'
                dist = min(dist,[],2);
        end
    else
        switch arg.tail
            case 'both'
                dist = abs(dist);
        end
    end

    % Compute p-value
    switch arg.tail
        case 'both'
            p = (sum(abs(stat)<=dist)+1)/(arg.nperm+1);
        case 'right'
            p = (sum(stat<=dist)+1)/(arg.nperm+1);
        case 'left'
            p = (sum(stat>=dist)+1)/(arg.nperm+1);
    end

end

% Compute confidence interval
if nargout > 2
    switch arg.type
        case 'ttest2'
            switch arg.tail
                case 'both'
                    crit = prctile(dist,100*(1-arg.alpha)).*se;
                    ci = [mu-crit;mu+crit];
                case 'right'
                    crit = prctile(dist,100*(1-arg.alpha)).*se;
                    ci = [mu-crit;Inf(1,nvar)];
                case 'left'
                    crit = prctile(-dist,100*(1-arg.alpha)).*se;
                    ci = [-Inf(1,nvar);mu+crit];
            end
        case 'ranksum'
            ci = NaN(2,nvar);
    end
end

% Store statistics in a structure
if nargout > 3
    switch arg.type
        case 'ttest2'
            stats.tstat = stat;
            stats.df = df;
            stats.sd = sd;
            stats.mean = mu;
        case 'ranksum'
            stats.zstat = stat;
            stats.iqr = iqrx;
            stats.median = medx;
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
            case 'ttest2'
                stats.tstat = stat;
                stats.df = ptvec2mat(df);
                stats.sd = ptvec2mat(sd);
                stats.mean = ptvec2mat(mu);
            case 'ranksum'
                stats.zstat = stat;
                stats.iqr = ptvec2mat(iqrx);
                stats.median = ptvec2mat(medx);
        end
    end
end