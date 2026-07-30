function [chi2,p,stats,tbl,dist] = permuchi2(x,varargin)
%PERMUCHI2  Permutation-based Chi-square test of independence.
%   CHI2 = PERMUCHI2(X) performs a permutation test of independence based
%   on the Chi-square statistic for the contingency table X, and returns
%   the test statistic. If X is a 3D matrix (Rows x Columns x Variables),
%   the test evaluates all variables simultaneously, and a vector of
%   results is returned.
%
%   PERMUCHI2 serves as a robust and computationally efficient alternative
%   to Fisher's exact test for sparse or complex contingency tables because
%   permuting the raw data preserves the table's marginal totals to
%   empirically model the exact probability space.
%
%   [CHI2,P] = PERMUCHI2(...) returns the probabilities (i.e. p-values)
%   of observing the given results by chance if the null hypothesis is
%   true. As the null distribution is generated empirically by permuting
%   the data, no assumption is made about the distribution that the data
%   come from. P-values are automatically adjusted for multiple comparisons
%   using the max correction method.
%
%   [CHI2,P,STATS] = PERMUCHI2(...) returns a structure with the following
%   fields:
%       'method'    -- the statistical method used
%       'df'        -- the degrees of freedom
%       'expfreq'   -- the expected frequencies matrix
%
%   [CHI2,P,STATS,TBL] = PERMUCHI2(...) returns the table contents as a
%   cell array.
%
%   [CHI2,P,STATS,TBL,DIST] = PERMUCHI2(...) returns the permuted sampling
%   distribution of the test statistic.
%
%   [...] = PERMUCHI2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      association is not equal to chance
%                       'right'     association is greater than chance (default)
%                       'left'      association is less than chance
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   control FWER using max correction (default=1).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator (default=shuffle).
%
%   See also CROSSTAB CHI2GOF PERMUANOVA1 PERMUTTEST.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Parse input arguments
arg = ptparsevarargin(varargin);
if ~any(strcmpi(varargin,'tail'))
    arg.tail = 'right';
end

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Handle explicit NaNs
if any(isnan(x(:)))
    x(isnan(x)) = 0;
    nanflag = 'omitnan';
else
    nanflag = 'includenan';
end

% Get data dimensions
[r,c,v] = size(x);

% Compute degrees of freedom
df = (r-1)*(c-1);

% Validate sample sizes
nobsall = squeeze(sum(x,[1,2]));
if any(nobsall ~= nobsall(1))
    error('The total number of observations must be equal across all variables for multivariate permutations.')
end
nobs = nobsall(1);

% Compute marginal sums
rowsums = sum(x,2);
colsums = sum(x,1);

% Compute expected frequencies
expfreq = (rowsums.*colsums)./nobs;

% Compute test statistic
chi2 = sum((x-expfreq).^2./expfreq,[1,2],nanflag);
chi2 = chi2(:)';

if nargout > 1

    % Generate random permutations
    rng(arg.seed);

    % Vectorized reconstruction of raw categorical vectors
    vara = zeros(nobs,v);
    varb = zeros(nobs,v);
    [I,J] = ndgrid(1:r,1:c);
    for k = 1:v
        xk = x(:,:,k);
        vara(:,k) = repelem(I(:),xk(:));
        varb(:,k) = repelem(J(:),xk(:));
    end

    % Precompute subscript arrays for 3D accumarray
    varK = repmat(1:v,nobs,1);
    subs13 = [vara(:),varK(:)];

    % Estimate sampling distribution
    dist = zeros(arg.nperm,v);
    for i = 1:arg.nperm
        randidx = randperm(nobs);
        varbperm = varb(randidx,:);
        xperm = accumarray([subs13(:,1),varbperm(:),subs13(:,2)],1,[r,c,v]);
        chi2perm = sum((xperm-expfreq).^2./expfreq,[1,2],nanflag);
        dist(i,:) = chi2perm(:)';
    end

    % Apply max-correction if specified
    if arg.correct
        distmax = max(dist,[],2);
        distmin = min(dist,[],2);
    else
        distmax = dist;
        distmin = dist;
    end

    % Compute p-values
    switch arg.tail
        case {'both','two'}
            p = min(1,2*(min(sum(chi2>=distmin,1),...
                sum(chi2<=distmax,1))+1)/(arg.nperm+1));
        case 'right'
            p = (sum(chi2<=distmax,1)+1)/(arg.nperm+1);
        case 'left'
            p = (sum(chi2>=distmin,1)+1)/(arg.nperm+1);
    end

end

% Store statistics in a structure
if nargout > 2
    stats.method = arg.type;
    stats.df = df;
    stats.expfreq = expfreq;
end

% Create Chi-square table
if nargout > 3
    if v == 1
        tbl = {
            'Statistic','df','Value','Prob>Chi2';
            'Pearson Chi2',df,chi2(1),p(1)};
    else
        tbl = cell(v+1,4);
        tbl(1,:) = {'Statistic','df','Value','Prob>Chi2'};
        for k = 1:v
            tbl(k+1,:) = {['Variable ',num2str(k)],df,chi2(k),p(k)};
        end
    end
end