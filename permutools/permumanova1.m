function [lambda,p,ci,stats,tbl,dist] = permumanova1(x,group,varargin)
%PERMUMANOVA1  Permutation-based one-way multivariate ANOVA and rank MANOVA.
%   LAMBDA = PERMUMANOVA1(X,GROUP) performs a permutation-based one-way
%   multivariate analysis of variance (MANOVA) for comparing the  means of
%   two or more groups of data in matrix X, and returns the test statistic.
%   X is an N-by-P matrix where N is observations and P is dependent
%   variables. GROUP is an N-by-1 grouping vector. If X is a 3D matrix
%   (N-by-P-by-V), the test evaluates all V sets of variables
%   simultaneously, and a vector of results is returned.
%
%   For non-normally distributed data, the raw data may be transformed
%   to rank orders in order to compute a rank MANOVA test by setting the
%   'type' parameter to 'rank'.
%
%   PERMUMANOVA1 treats NaNs as missing values, and ignores them.
%
%   [LAMBDA,P] = PERMUMANOVA1(...) returns the probability (i.e. p-value)
%   of observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data, no assumption is made about the shape of the distributions that
%   the data come from, except that they have equal variances. P-values are
%   automatically adjusted for multiple comparisons using the max
%   correction method.
%
%   [LAMBDA,P,CI] = PERMUMANOVA1(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for Wilks' Lambda. CIs are also
%   adjusted for multiple comparisons using the max correction method.
%
%
%   [LAMBDA,P,CI,STATS] = PERMUMANOVA1(...) returns a structure with the
%   following fields:
%       'source'    -- the function used to compute the MANOVA
%       'method'    -- the statistical method used
%       'gnames'    -- the group names
%       'n'         -- the group sample sizes
%       'df'        -- the degrees of freedom (between and within)
%
%   [LAMBDA,P,CI,STATS,TBL] = PERMUMANOVA1(...) returns the MANOVA table
%   contents as a cell array.
%
%   [LAMBDA,P,CI,STATS,TBL,DIST] = PERMUMANOVA1(...) returns the permuted
%   sampling distribution of the test statistic.
%
%   [...] = PERMUMANOVA1(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'manova1'       one-way MANOVA (default)
%                       'rank'          rank-transformed MANOVA
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   control FWER using max correction (default=1).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator (default=shuffle).
%
%   See also MANOVA1 PERMUANOVA1 PERMUREGRESS.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(group)
    error('A grouping variable is required for MANOVA.')
end

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'manova1';
end

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Handle missing values
nanidx = any(isnan(x),2);
if ndims(x) == 3
    nanidx = any(nanidx,3);
end
if any(nanidx)
    x(nanidx,:,:) = [];
    group(nanidx) = [];
end

% Get data dimensions
[nobs,pnum,v] = size(x);

% Transform raw data to rank orders if specified
switch lower(arg.type)
    case 'manova1'
    case 'rank'
        x = reshape(tiedrank(reshape(x,nobs,[])),nobs,pnum,v);
    otherwise
        error('The TYPE parameter value must be MANOVA1, or RANK.')
end

% Parse grouping variable
[group,gnames] = grp2idx(group);
groups = unique(group)';
k = numel(groups);

% Dummy indicator matrix for groups
G = zeros(nobs,k);
for i = 1:k
    G(:,i) = (group == groups(i));
end
n = sum(G,1);
diagn = diag(n);

% Preallocate memory
lambda = zeros(1,v);
T = zeros(pnum,pnum,v);
detT = zeros(1,v);
gm = zeros(1,pnum,v);

% Compute true Wilks' Lambda statistics
for kidx = 1:v

    xk = x(:,:,kidx);
    gm(:,:,kidx) = sum(xk,1)/nobs;

    % Compute Total SSCP matrix
    xctot = xk-gm(:,:,kidx);
    Tk = xctot'*xctot;
    T(:,:,kidx) = Tk;
    detT(kidx) = det(Tk);

    % Compute Between-group SSCP matrix (Vectorized)
    S = G'*xk;
    M = S./n';
    Mc = M-gm(:,:,kidx);
    Bk = Mc'*diagn*Mc;

    % Compute Within-group SSCP matrix
    Wk = Tk-Bk;
    lambda(kidx) = det(Wk)/detT(kidx);

end

if nargout > 1

    % Generate random permutations
    rng(arg.seed);

    % Estimate sampling distribution
    dist = zeros(arg.nperm,v);
    for i = 1:arg.nperm
        randidx = randperm(nobs);
        Gp = G(randidx,:);
        for kidx = 1:v
            xk = x(:,:,kidx);
            Sp = Gp'*xk;
            Mp = Sp./n';
            Mc = Mp-gm(:,:,kidx);
            Bp = Mc'*diagn*Mc;
            Wp = T(:,:,kidx)-Bp;
            dist(i,kidx) = det(Wp)/detT(kidx);
        end
    end

    % Apply max-correction if specified
    if arg.correct
        distmin = min(dist,[],2);
    else
        distmin = dist;
    end

    % Compute p-values
    p = (sum(lambda>=distmin,1)+1)/(arg.nperm+1);

end

if nargout > 2
    if arg.correct
        crit = prctile(distmin,100*arg.alpha);
        crit = crit*ones(1,v);
    else
        crit = prctile(dist,100*arg.alpha,1);
    end
    ci = [zeros(1,v);lambda./crit];
end

% Store statistics in a structure
if nargout > 3
    stats.source = 'permumanova1';
    stats.method = arg.type;
    stats.gnames = gnames;
    stats.n = n;
    stats.df = [k-1,nobs-k];
end

% Create MANOVA table
if nargout > 4
    if v == 1
        tbl = {
            'Statistic','Value','df','Prob>Lambda';
            'Wilks Lambda',lambda(1),k-1,p(1)};
    else
        tbl = cell(v+1,4);
        tbl(1,:) = {'Statistic','Value','df','Prob>Lambda'};
        for kidx = 1:v
            tbl(kidx+1,:) = {['Variable ',num2str(kidx)],...
                lambda(kidx),k-1,p(kidx)};
        end
    end
end