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
%                       'right'     right-tailed test (default)
%                       'left'      left-tailed test
%                       'both'      two-tailed test
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric or logical value specifying whether to apply
%                   max-correction for multiple comparisons (default=1).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator.
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

% Get data dimensions
[r,c,v] = size(x);

% Compute degrees of freedom
df = (r-1)*(c-1);

% Validate max-correction sample sizes
if arg.correct
    nobsall = squeeze(sum(x,[1,2]));
    if any(nobsall ~= nobsall(1))
        error('For max-correction, the total number of observations must be equal across all variables.')
    end
end
nobs = sum(x(:,:,1),'all');

% Compute marginal sums
rowsums = sum(x,2);
colsums = sum(x,1);

% Compute expected frequencies
expfreq = (rowsums.*colsums)./nobs;

% Compute test statistic
chi2 = squeeze(sum((x-expfreq).^2./expfreq,[1,2]))';

if nargout > 1

    % Generate random permutations
    rng(arg.seed);

    % Reconstruct raw categorical vectors for all variables
    vara = zeros(nobs,v);
    varb = zeros(nobs,v);
    for k = 1:v
        idx = 1;
        for i = 1:r
            for j = 1:c
                count = x(i,j,k);
                if count > 0
                    vara(idx:idx+count-1,k) = i;
                    varb(idx:idx+count-1,k) = j;
                    idx = idx+count;
                end
            end
        end
    end

    % Estimate uncorrected sampling distribution
    dist = zeros(arg.nperm,v);
    for i = 1:arg.nperm
        randidx = randperm(nobs);
        varbperm = varb(randidx,:);
        for k = 1:v
            xperm = accumarray([vara(:,k),varbperm(:,k)],1,[r,c]);
            dist(i,k) = sum((xperm-expfreq(:,:,k)).^2./expfreq(:,:,k),'all');
        end
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
        case 'both'
            p = min(1,2*(min(sum(chi2>=distmin),...
                sum(chi2<=distmax))+1)/(arg.nperm+1));
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