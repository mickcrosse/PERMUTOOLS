function [stat,p,ci,stats,dist] = permucorr(x,varargin)
%PERMUCORR  Permutation-based Pearson, Spearman's rank, and rankit correlation coefficient.
%   STAT = PERMUCORR(X) returns a matrix containing the pairwise linear
%   correlation coefficients between each pair of columns in X based on
%   Pearson's r.
%
%   For nonlinear correlations, the raw data may be transformed to rank
%   orders in order to compute a Spearman's rank correlation by setting
%   the 'type' parameter to 'spearman' or 'rank', or a rankit correlation
%   by setting the 'type' parameter to 'rankit'.
%
%   PERMUCORR treats NaNs as missing values, and ignores them.
%
%   STAT = PERMUCORR(X,Y) returns the pairwise correlation coefficient
%   between vectors X and Y. X and Y must have the same length. If X and Y
%   are matrices, the correlation coefficients between each corresponding
%   pair of columns in X and Y are calculated, and a vector of results is
%   returned.
%
%   [STAT,P] = PERMUCORR(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data, no assumption is made about the shape of the distribution that
%   the data come from. When only one sample is entered in X, two-tailed
%   permutation tests are automatically used. P-values are automatically
%   adjusted for multiple comparisons using the max correction method.
%
%   [STAT,P,CI] = PERMUCORR(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for each coefficient. CIs are also adjusted for multiple
%   comparisons using the max correction method.
%
%   [STAT,P,CI,STATS] = PERMUCORR(...) returns a structure with the
%   following fields:
%       'df'        -- the degrees of freedom of each measure
%       'method'    -- the correlation method used
%
%   [STAT,P,CI,STATS,DIST] = PERMUCORR(...) returns the permuted sampling
%   distribution of the test statistic.
%
%   [...] = PERMUCORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of correlation coefficient
%                   to compute:
%                       'pearson'   Pearson correlation coefficient
%                                   (default)
%                       'spearman'  Spearman's rank correlation coefficient
%                       'rankit'    rankit correlation coefficient
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      correlation is not zero (default)
%                       'right'     correlation is greater than zero
%                       'left'      correlation is less than zero
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000 or all possible permutations for less
%                   than 8 observations).
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
%   See also CORR CORRCOEF PARTIALCORR TIEDRANK.
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
%       [4] Bishara AJ, Hittner JB, (2012) Testing the Significance of a
%           Correlation With Nonnormal Data: Comparison of Pearson,
%           Spearman, Transformation, and Resampling Approaches. Psychol
%           Methods, 17(3):399-417.
%       [5] Bishara AJ, Hittner JB, (2017) Confidence intervals for
%           correlations when data are not normal. Behav Res, 49:294-309.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin<2 || ischar(varargin{1})
    y = [];
else
    y = varargin{1};
    varargin = varargin(2:end);
end

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'pearson';
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
    warning('Comparing all columns of X using two-tailed correlations...')
    [x,y] = ptpaircols(x);
    arg.tail = 'both';
    arg.matrix = true;
end
if size(x)~=size(y)
    error('X and Y must be the same size.')
end

% Use only rows with no NaN values if specified
switch arg.rows
    case 'all'
        x(isnan(y)) = NaN;
        y(isnan(x)) = NaN;
    case 'complete'
        valid_rows = ~any(isnan(x),2) & ~any(isnan(y),2);
        x = x(valid_rows,:);
        y = y(valid_rows,:);
end

% Get data dimensions
[maxnobs,nvar] = size(x);
nobs = sum(~isnan(x));

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:))) || any(isnan(y(:)))
    nanflag = 'omitnan';
else
    nanflag = 'includenan';
end

% Transform raw data to rank orders if specified
switch arg.type
    case 'pearson'
    case {'spearman','rank'}
        x = tiedrank(x);
        y = tiedrank(y);
    case 'rankit'
        x = norminv((tiedrank(x)-0.5)./nobs);
        y = norminv((tiedrank(y)-0.5)./nobs);
    otherwise
        error('The TYPE parameter value must be PEARSON, SPEARMAN, or RANKIT.')
end

% Compute observed statistics
sdxy = sqrt((sum(x.^2,nanflag)-(sum(x,nanflag).^2)./nobs)...
    .*(sum(y.^2,nanflag)-(sum(y,nanflag).^2)./nobs));
mu = sum(x,nanflag).*sum(y,nanflag)./nobs;
stat = (sum(x.*y,nanflag)-mu)./sdxy;

if nargout > 1

    rng(arg.seed);

    % Generate random permutations
    if maxnobs < 8
        arg.nperm = factorial(maxnobs);
        idx = perms(1:maxnobs)';
        if arg.verbose
            warning('Computing all possible permutations due to small N.')
            fprintf('Number of permutations used: %d\n',arg.nperm)
        end
    else
        [~,idx] = sort(rand(maxnobs,arg.nperm));
    end

    % Estimate sampling distribution
    dist = zeros(arg.nperm,nvar);
    switch nanflag
        case 'omitnan'
            % Exact pair-wise calculation for missing data
            for i = 1:arg.nperm
                xp = x(idx(:,i),:);
                valid = ~isnan(xp) & valid_y;
                xp(~valid) = 0;
                yp = y;
                yp(~valid) = 0;
                nobsp = sum(valid,1);
                sumx = sum(xp,1);
                sumy = sum(yp,1);
                sdxyp = sqrt((sum(xp.^2,1)-(sumx.^2)./nobsp).*...
                    (sum(yp.^2,1)-(sumy.^2)./nobsp));
                mup = sumx.*sumy./nobsp;
                dist(i,:) = (sum(xp.*yp,1)-mup)./sdxyp;
            end
        case 'includenan'
            % Fast vectorized calculation for complete data
            for j = 1:nvar
                xp = x(:,j);
                yp = y(:,j);
                sumxy = yp'*xp(idx);
                dist(:,j) = (sumxy'-mu(j))./sdxy(j);
            end
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
    zstat = atanh(stat);
    zdist = atanh(refdist);
    switch arg.tail
        case 'both'
            crit = prctile(abs(zdist), 100*(1-arg.alpha));
            ci = [tanh(zstat - crit); tanh(zstat + crit)];
        case 'right'
            crit = prctile(zdist, 100*(1-arg.alpha));
            ci = [tanh(zstat - crit); Inf(1,nvar)];
        case 'left'
            crit = prctile(-zdist, 100*(1-arg.alpha));
            ci = [-Inf(1,nvar); tanh(zstat + crit)];
    end
end

% Format descriptive statistics
if nargout > 3
    df = nobs-2;
    stats.df = df;
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
        stats.df = ptvec2mat(df);
    end
end