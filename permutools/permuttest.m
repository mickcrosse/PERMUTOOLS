function [stat,p,ci,stats,dist] = permuttest(x,m,varargin)
%PERMUTTEST  Permutation-based one/paired-sample t-test and Wilcoxon signed-rank test.
%   STAT = PERMUTTEST(X) performs a one-sample permutation test based on
%   the t-statistic of the hypothesis that the data in X come from a
%   distribution with mean zero, and returns the test statistic. If X is a
%   matrix, separate permutation tests are performed along each column of
%   X, and a vector of results is returned. If the 'compare' parameter is
%   set to 'pairwise', two-tailed permutation tests between every pair of
%   columns in X are performed, and a matrix of results is returned.
%
%   For non-normally distributed samples, the raw data may be transformed
%   to signed ranks in order to compute a Wilcoxon signed-rank test by
%   setting the 'type' parameter to 'signrank' or 'rank'.
%
%   PERMUTTEST treats NaNs as missing values, and ignores them.
%
%   STAT = PERMUTTEST(X,M) returns the results of a one-sample permutation
%   test of the hypothesis that the data in X come from a distribution with
%   mean M. M must be a scalar.
%
%   STAT = PERMUTTEST(X,Y) returns the results of a paired-sample
%   permutation test between dependent samples X and Y of the hypothesis
%   that the data in X and Y come from distributions with equal means. X
%   and Y must have the same length. If X and Y are matrices, separate
%   permutation tests are performed between each corresponding pair of
%   columns in X and Y, and a vector of results is returned.
%
%   [STAT,P] = PERMUTTEST(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by permuting the data,
%   no assumption is made about the shape of the distribution that the data
%   come from. P-values are automatically adjusted for multiple comparisons
%   using the max correction method.
%
%   [STAT,P,CI] = PERMUTTEST(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for the true mean of X, or of X-Y for a paired test.
%   For rank-based tests ('type','signrank'), CI returns NaNs. CIs are also
%   adjusted for multiple comparisons using the max correction method.
%
%   [STAT,P,CI,STATS] = PERMUTTEST(...) returns a structure with the
%   following fields for standard tests ('type','ttest'):
%       'tstat'     -- the value of the test statistic
%       'df'        -- the degrees of freedom of each test
%       'sd'        -- the estimated population standard deviation based on
%                      X, or on X-Y for a paired test
%       'mean'      -- the estimated population mean based on X, or on X-Y
%                      for a paired test
%       'method'    -- the statistical method used
%
%   For rank-based tests ('type','signrank'), STATS contains:
%       'zstat'     -- the value of the test statistic
%       'iqr'       -- the estimated population interquartile range based
%                      on X, or on X-Y for a paired test
%       'median'    -- the estimated population median based on X, or on
%                      X-Y for a paired test
%       'method'    -- the statistical method used
%
%   [STAT,P,CI,STATS,DIST] = PERMUTTEST(...) returns the permuted sampling
%   distribution of the test statistic.
%
%   [...] = PERMUTTEST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'ttest'     one/paired-sample t-test (default)
%                       'signrank'  Wilcoxon signed-rank test
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      mean is not M (default)
%                       'right'     mean is greater than M
%                       'left'      mean is less than M
%       'compare'   A string specifying what to compare each variable to
%                   when only X is entered:
%                       'zero'      compare each column of X to zero and
%                                   return a vector of results (default)
%                       'pairwise'  compare every pair of columns in X to
%                                   each other using two-tailed tests and
%                                   return a matrix of results
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
%       'verbose'   A numeric scalar (0,1) or logical indicating whether to
%                   run in verbose mode (default=1).
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

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
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
if isempty(arg.type)
    arg.type = 'ttest';
elseif strcmpi(arg.type,'rank')
    arg.type = 'signrank';
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
    switch arg.compare
        case 'zero'
            y = zeros(size(x));
        case 'pairwise'
            warning('Comparing all columns of X using two-tailed tests...')
            [x,y] = ptpaircols(x);
            arg.tail = 'both';
            arg.matrix = true;
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

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitnan';
else
    nanflag = 'includenan';
end

% Compute observed statistics
switch arg.type
    case 'ttest'
        df = nobs-1;
        sd = std(x,nanflag);
        mu = sum(x,nanflag)./nobs;
        se = sd./sqrt(nobs);
        stat = mu./se;
    case 'signrank'
        if nargout > 3
            iqrx = iqr(x);
            medx = median(x,nanflag);
        end
        [r,tieadj] = tiedrank(abs(x));
        x = sign(x).*r;
        w = sum(x,1,nanflag);
        varw = (nobs.*(nobs+1).*(2.*nobs+1))/6-(tieadj/6);
        stat = w./sqrt(varw);
    otherwise
        error('The TYPE parameter value must be TTEST, or SIGNRANK.')
end

if nargout > 1

    rng(arg.seed);

    % Generate random permutations
    signp = double(randi([0,1],arg.nperm,maxnobs,'int8')*2-1);

    % Estimate sampling distribution
    switch nanflag
        case 'omitnan'
            x(isnan(x)) = 0;
    end
    sump = signp*x;
    switch arg.type
        case 'ttest'
            sumsq = sum(x.^2,nanflag);
            varp = (sumsq-(sump.^2)./nobs)./df;
            dist = (sump./nobs)./sqrt(varp./nobs);
        case 'signrank'
            dist = sump./sqrt(varw);
    end

    % Add negative values
    dist(arg.nperm+1:2*arg.nperm,:) = -dist;
    arg.nperm = 2*arg.nperm;
    if arg.verbose
        fprintf('Adding negative of values to permutation distribution.\n')
        fprintf('Number of permutations used: %d\n',arg.nperm)
    end

    % Apply max correction if specified
    if arg.correct
        switch arg.tail
            case 'both'
                refdist = max(abs(dist),[],2);
            case 'right'
                refdist = max(dist,[],2);
            case 'left'
                refdist = min(dist,[],2);
        end
    else
        refdist = dist;
    end

    % Compute p-value
    switch arg.tail
        case 'both'
            p = (sum(abs(stat)<=abs(refdist))+1)/(arg.nperm+1);
        case 'right'
            p = (sum(stat<=refdist)+1)/(arg.nperm+1);
        case 'left'
            p = (sum(stat>=refdist)+1)/(arg.nperm+1);
    end

end

% Compute confidence interval
if nargout > 2
    switch arg.type
        case 'ttest'
            switch arg.tail
                case 'both'
                    crit = prctile(abs(refdist),100*(1-arg.alpha)).*se;
                    ci = [mu-crit;mu+crit];
                case 'right'
                    crit = prctile(refdist,100*(1-arg.alpha)).*se;
                    ci = [mu-crit;Inf(1,nvar)];
                case 'left'
                    crit = prctile(-refdist,100*(1-arg.alpha)).*se;
                    ci = [-Inf(1,nvar);mu+crit];
            end
        case 'signrank'
            ci = NaN(2,nvar);
    end
end

% Format descriptive statistics
if nargout > 3
    switch arg.type
        case 'ttest'
            stats.tstat = stat;
            stats.df = df;
            stats.sd = sd;
            stats.mean = mu;
        case 'signrank'
            stats.zstat = stat;
            stats.iqr = iqrx;
            stats.median = medx;
    end
    stats.method = arg.type;
end

% Format sampling distribution
if nargout > 4
    if arg.correct
        dist = refdist;
    end
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
            case 'ttest'
                stats.tstat = stat;
                stats.df = ptvec2mat(df);
                stats.sd = ptvec2mat(sd);
                stats.mean = ptvec2mat(mu);
            case 'signrank'
                stats.zstat = stat;
                stats.iqr = ptvec2mat(iqrx);
                stats.median = ptvec2mat(medx);
        end
    end
end