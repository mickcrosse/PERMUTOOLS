function [r,p,ci] = permucorr(x,y,varargin)
%PERMUCORR  Linear or rank permutation-based correlation.
%   R = PERMUCORR(X,Y) returns the pairwise linear correlation coefficient
%   based on Pearson's r. For nonlinear correlations, the raw data may be
%   transformed to rank orders using the Spearman's or a rankit method via
%   the 'TYPE' parameter (Bishara & Hittner, 2012). X and Y must have the
%   same length.
%
%   If X and Y are matrices, the correlation between each corresponding
%   pair of columns in X and Y is calculated, and a vector of results is
%   returned.
%
%   If Y is not entered, the correlation between each pair of columns in X
%   is computed and a correlation matrix is returned.
%
%   [R,P] = PERMUCORR(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data, no assumption is made about the shape of the distribution that
%   the data come from. For multivariate data, Family-wise error rate
%   (FWER) is controlled using the maximum statistic correction method
%   (Blair et al., 1994). This method provides strong control of FWER, even
%   for small sample sizes, and is much more powerful than traditional
%   correction methods (Groppe et al., 2011).
%
%   [R,P,CI] = PERMUCORR(...) returns a 100*(1-ALPHA)% confidence interval
%   for each coefficient. For Pearson and rankit correlations, the Fisher
%   z' method is used and for Spearman rank correlations, Fieller et al.'s
%   estimate is used (Bishara & Hittner, 2017).
%
%   [...] = PERMUCORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%                       'both'      correlation is not zero (default)
%                       'right'     correlation is greater than zero
%                       'left'      correlation is less than zero
%       'type'      A string specifying the type of correlation measure:
%                       'Pearson'   Pearson's correlation coefficient (def)
%                       'Spearman'  Spearman's rank correlation coefficient
%                       'Rankit'    Bliss's rankit correlation coefficient
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000 or all possible permutations for less
%                   than 14 observations).
%       'correct'   A numeric scalar (0,1) or logical indicating whether
%                   to control FWER using rmax correction (default=true).
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator. By default, the
%                   generator is initialised based on the current time,
%                   resulting in a different permutation on each call.
%
%   See also CORR CORRCOEF PARTIALCORR TIEDRANK.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Bishara AJ, Hittner JB, (2012) Testing the Significance of a
%           Correlation With Nonnormal Data: Comparison of Pearson,
%           Spearman, Transformation, and Resampling Approaches. Psychol
%           Methods, 17(3):399-417.
%       [2] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.
%       [4] Bishara AJ, Hittner JB, (2017) Confidence intervals for
%           correlations when data are not normal. Behav Res, 49:294-309.

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
    warning('Comparing all columns of X in a correlation matrix...')
    [x,y] = ptpaircols(x);
    arg.tail = 'both';
    arg.mat = true;
end
if size(x)~=size(y)
    error('X and Y must be the same size.')
end

% Use only rows with no NaN values if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(y),2),:);
        y = y(~any(isnan(y),2),:);
        y = y(~any(isnan(x),2),:);
        x = x(~any(isnan(x),2),:);
    case 'all'
        if any(isnan(x(:))) || any(isnan(y(:)))
            error('X or Y is missing values. Set ROWS to ''complete''.')
        end
end

% Get data dimensions
[nobs,nvar] = size(x);

% Transform raw data to rank-orders if specified
switch arg.type
    case 'Rankit'
        x = norminv((tiedrank(x)-0.5)/nobs);
        y = norminv((tiedrank(y)-0.5)/nobs);
    case 'Spearman'
        x = tiedrank(x);
        y = tiedrank(y);
end

% Compute test statistic
muxy = sum(x).*sum(y)/nobs;
sdxy = sqrt((sum(x.^2)-(sum(x).^2)/nobs).*(sum(y.^2)-(sum(y).^2)/nobs));
r = (sum(x.*y)-muxy)./sdxy;

if nargout > 1

    % Use all possible permutations if less than 8 observations
    if nobs < 8
        warning('Computing all possible permutations due to small N.')
        arg.nperm = factorial(nobs);
        idx = perms(1:nobs)';
    else
        rng(arg.seed);
        [~,idx] = sort(rand(nobs,arg.nperm));
    end

    % Generate permutation distribution
    pd = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        pd(i,:) = (sum(x(idx(:,i),:).*y)-muxy)./sdxy;
    end

    % Apply max correction if specified
    if arg.correct
        switch arg.tail
            case 'both'
                [~,idx] = max(abs(pd),[],2);
                csvar = [0;cumsum(ones(arg.nperm-1,1)*nvar)];
                pd = pd';
                pd = pd(idx+csvar);
            case 'right'
                pd = max(pd,[],2);
            case 'left'
                pd = min(pd,[],2);
        end
    end

    % Compute unadjusted test statistics
    switch arg.tail
        case 'both'
            p = 2*(min(sum(r<=pd),sum(r>=pd))+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(pd,100*(1-arg.alpha/2));
                ci = [max(-1,r-crit);min(1,r+crit)];
            end
        case 'right'
            p = (sum(r<=pd)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(pd,100*(1-arg.alpha));
                ci = [max(-1,r-crit);Inf(1,nvar)];
            end
        case 'left'
            p = (sum(r>=pd)+1)/(arg.nperm+1);
            if nargout > 2
                crit = prctile(-pd,100*(1-arg.alpha));
                ci = [-Inf(1,nvar);min(1,r+crit)];
            end
    end

end

% Arrange test results in a matrix if specified
if arg.mat
    r = ptvec2mat(r);
    if nargout > 1
        p = ptvec2mat(p);
    end
    if nargout > 2
        ciLwr = ptvec2mat(ci(1,:));
        ciUpr = ptvec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
    end
end