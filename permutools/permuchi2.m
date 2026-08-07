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
%       'chi2stat'              -- the value of the test statistic
%       'df'                    -- the degrees of freedom
%       'O'                     -- the observed count in each cell
%       'E'                     -- the expected count in each cell
%       'cramersv'              -- Cramer's V effect size (strength of association)
%       'oddsratio'             -- the odds ratio (for 2x2 tables only)
%       'oddsratioci'           -- the asymptotic confidence interval for
%                                  the odds ratio (for 2x2 tables only)
%       'method'                -- the statistical method used
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

% Replace NaNs with zeros
if any(isnan(x(:)))
    x(isnan(x)) = 0;
    warning('Replacing NaN values with 0 assuming they represent empty categories.');
end

% Get data dimensions
[nrow,ncol,nvar] = size(x);

% Compute degrees of freedom
df = (nrow-1)*(ncol-1);

% Validate sample sizes
nobsall = squeeze(sum(x,[1,2]));
if any(nobsall ~= nobsall(1))
    error('The total number of observations must be equal across all variables for multivariate permutations.')
end
nobs = nobsall(1);

% Compute observed statistics
rowsums = sum(x,2);
colsums = sum(x,1);
E = (rowsums.*colsums)./nobs;
chi2 = sum((x-E).^2./E,[1,2]);
chi2 = chi2(:)';

if nargout > 1

    rng(arg.seed);

    % Vectorized reconstruction of raw categorical vectors
    vara = zeros(nobs,nvar);
    varb = zeros(nobs,nvar);
    [I,J] = ndgrid(1:nrow,1:ncol);
    for k = 1:nvar
        xk = x(:,:,k);
        vara(:,k) = repelem(I(:),xk(:));
        varb(:,k) = repelem(J(:),xk(:));
    end

    % Precompute subscript arrays for 3D accumarray
    vark = repmat(1:nvar,nobs,1);
    subs13 = [vara(:),vark(:)];

    % Estimate sampling distribution
    dist = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        randidx = randperm(nobs);
        varbperm = varb(randidx,:);
        xperm = accumarray([subs13(:,1),varbperm(:),subs13(:,2)],1,...
            [nrow,ncol,nvar]);
        chi2perm = sum((xperm-E).^2./E,[1,2],nanflag);
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
    stats.chi2stat = chi2;
    stats.df = df;
    stats.O = x;
    stats.E = E;
    stats.cramersv = sqrt(chi2./(nobs*(min(nrow,ncol)-1)));
    stats.oddsratio = NaN(1,nvar);
    stats.oddsratioci = NaN(2,nvar);
    stats.method = 'chi2test';
    if nrow == 2 && ncol == 2
        stats.oddsratio = squeeze((x(1,1,:).*x(2,2,:))./(x(1,2,:).*...
            x(2,1,:)))';
        z = norminv(1-arg.alpha/2);
        for k = 1:nvar
            if any(x(:,:,k) == 0, 'all')
                stats.oddsratioci(:, k) = [-Inf;Inf];
            else
                selogor = sqrt(sum(1./x(:,:,k),'all'));
                logor = log(stats.oddsratio(k));
                stats.oddsratioci(1,k) = exp(logor-z*selogor);
                stats.oddsratioci(2,k) = exp(logor+z*selogor);
            end
        end
    end
end

% Create Chi-square table
if nargout > 3
    if nvar == 1
        tbl = {
            'Statistic','df','Value','Prob>Chi2';
            'Pearson Chi2',df,chi2(1),p(1)};
    else
        tbl = cell(nvar+1,4);
        tbl(1,:) = {'Statistic','df','Value','Prob>Chi2'};
        for k = 1:nvar
            tbl(k+1,:) = {['Variable ',num2str(k)],df,chi2(k),p(k)};
        end
    end
end