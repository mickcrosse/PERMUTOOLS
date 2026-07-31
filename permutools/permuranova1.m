function [f,p,ci,stats,tbl,dist] = permuranova1(x,varargin)
%PERMURANOVA1  Permutation-based one-way repeated-measures ANOVA and Friedman test.
%   F = PERMURANOVA1(X) performs a balanced permutation-based one-way 
%   repeated-measures analysis of variance (ANOVA) for comparing the means 
%   of two or more columns of data in matrix X (where rows represent 
%   subjects and columns represent conditions), and returns the test 
%   statistic.
%
%   For non-normally distributed data, the raw data may be transformed
%   to rank orders in order to compute a one-way Friedman test by setting 
%   the 'type' parameter to 'friedman' or 'rank'.
%
%   PERMURANOVA1 requires a fully balanced design and does not support 
%   missing values (NaNs).
%
%   [F,P] = PERMURANOVA1(...) returns the probabilities (i.e. p-values) of
%   observing the given results by chance if the null hypothesis (that the
%   means of the groups are equal) is true. As the null distribution is
%   generated empirically by permuting the data within each subject, no
%   assumption is made about the shape of the distributions.
%
%   [F,P,CI] = PERMURANOVA1(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for the true difference of population means.
%
%   [F,P,CI,STATS] = PERMURANOVA1(...) returns a structure with the
%   following fields:
%       'source'    -- the function used to compute the ANOVA
%       'method'    -- the statistical method used
%       'sigmasq'   -- the error mean square
%       'colmeans'  -- the column means
%       'coln'      -- the column sample sizes (subjects)
%       'df'        -- the error degrees of freedom
%
%   [F,P,CI,STATS,TBL] = PERMURANOVA1(...) returns the ANOVA table contents
%   as a cell array.
%
%   [F,P,CI,STATS,TBL,DIST] = PERMURANOVA1(...) returns the permuted
%   sampling distributions of the test statistics.
%
%   [...] = PERMURANOVA1(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'ranova1'   one-way repeated-measures ANOVA (default)
%                       'friedman'  one-way Friedman rank test
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows.
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator (default=shuffle).
%
%   See also RANOVA PERMUANOVA1 PERMUANOVA2 PERMUTTEST BOOTEFFECTSIZE.
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
if isempty(arg.type)
    arg.type = 'ranova1';
end

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Orient data column-wise
if arg.dim == 2
    x = x';
end

% Balanced design check
if any(isnan(x(:)))
    error('Repeated-measures ANOVA requires a perfectly balanced design. Data must not contain NaNs.')
end

% Transform raw data to rank orders if specified
switch arg.type
    case 'ranova1'
    case {'friedman','rank'}
        x = tiedrank(x')';
    otherwise
        error('The TYPE parameter value must be RANOVA1, or FRIEDMAN.')
end

% Get data dimensions
[n,k] = size(x);
nobs = n*k;

% Compute grand mean
gm = sum(x,'all')/nobs;

% Compute sum of squares
colmeans = mean(x,1);
rowmeans = mean(x,2);
sst = sum((x-gm).^2,'all');
ssb = k*sum((rowmeans-gm).^2,'all');
ssw = sst-ssb;
ssc = n*sum((colmeans-gm).^2,'all');
sse = ssw-ssc;

% Compute degrees of freedom
dfc = k-1;
dfe = (n-1)*(k-1);

% Compute mean squares
msc = ssc/dfc;
mse = sse/dfe;

% Compute F-statistic
f = msc/mse;

if nargout > 1

    % Generate random permutations
    rng(arg.seed);
    rowidx = repmat((1:n)',1,k);

    % Estimate sampling distribution
    dist = zeros(arg.nperm,1);
    for i = 1:arg.nperm
        [~,idx] = sort(rand(n,k),2);
        linidx = rowidx+(idx-1)*n;
        xp = x(linidx);
        colmeansp = mean(xp,1);
        sscp = n*sum((colmeansp-gm).^2,'all');
        ssep = ssw-sscp;
        mscp = sscp/dfc;
        msep = ssep/dfe;
        dist(i) = mscp/msep;
    end

    % Compute p-value
    p = (sum(f<=dist)+1)/(arg.nperm+1);

end

% Compute confidence intervals
if nargout > 2
    crit = prctile(dist,100*(1-arg.alpha));
    ci = [f./crit;Inf];
end

% Store statistics in a structure
if nargout > 3
    stats.source = 'permuranova1';
    stats.method = arg.type;
    stats.sigmasq = mse;
    stats.colmeans = colmeans;
    stats.coln = n;
    stats.df = dfe;
end

% Create ANOVA table
if nargout > 4
    tbl = {
        'Source','SS','df','MS','F','Prob>F';
        'Columns',ssc,dfc,msc,f,p;
        'Error',sse,dfe,mse,[],[];
        'Total',sst,nobs-1,[],[],[]};
end