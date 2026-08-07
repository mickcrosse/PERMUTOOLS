function [f,p,ci,stats,tbl,dist] = permuanovan(x,group,varargin)
%PERMUANOVAN  Permutation-based N-way ANOVA and rank ANOVA.
%   F = PERMUANOVAN(X,GROUP) performs a permutation-based N-way analysis
%   of variance (ANOVA) for comparing the means of data in matrix X, and
%   returns the test statistics. X is an N-by-V matrix where N is
%   observations and V is dependent variables. GROUP is an N-by-G matrix
%   or cell array containing G categorical factors.
%
%   For non-normally distributed data, the raw data may be transformed
%   to rank orders in order to compute a rank ANOVA test by setting the
%   'type' parameter to 'rank'.
%
%   PERMUANOVAN treats NaNs as missing values, and ignores them.
%
%   [F,P] = PERMUANOVAN(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data using the Freedman-Lane algorithm, no assumption is made about the
%   shape of the distributions that the data come from. P-values are
%   automatically adjusted for multiple comparisons using the max
%   correction method.
%
%   [F,P,CI] = PERMUANOVAN(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for the F-statistic. CIs are also adjusted for multiple
%   comparisons using the max correction method.
%
%   [F,P,CI,STATS] = PERMUANOVAN(...) returns a structure with the
%   following fields:
%       'source'    -- the function used to compute the ANOVA
%       'gnames'    -- the group names for all factors
%       'df'        -- the degrees of freedom [Factors 1 to G, Error]
%       'method'    -- the statistical method used
%
%   [F,P,CI,STATS,TBL] = PERMUANOVAN(...) returns the ANOVA table
%   contents as a cell array.
%
%   [F,P,CI,STATS,TBL,DIST] = PERMUANOVAN(...) returns the permuted
%   sampling distribution of the test statistics.
%
%   [...] = PERMUANOVAN(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'anovan'        N-way ANOVA (default)
%                       'rank'          rank-transformed ANOVA
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   control FWER using max correction (default=1).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator (default=shuffle).
%
%   See also ANOVAN PERMUANOVA1 PERMUANOVA2 PERMUREGRESS.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(group)
    error('A grouping matrix or cell array is required for ANOVAN.')
end

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'anovan';
end

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Standardize grouping variable format
if iscell(group)
    gmat = zeros(size(x,1),numel(group));
    for i = 1:numel(group)
        gmat(:,i) = grp2idx(group{i});
    end
    group = gmat;
end

% Parse grouping variables and handle missing values
nfactors = size(group,2);
g = zeros(size(group));
gnames = cell(1,nfactors);
for i = 1:nfactors
    [g(:,i),gnames{i}] = grp2idx(group(:,i));
end
nanidx = any(isnan(x),2) | any(isnan(g),2);
if any(nanidx)
    x(nanidx,:) = [];
    g(nanidx,:) = [];
end

% Get data dimensions
[nobs,v] = size(x);

% Transform raw data to rank orders if specified
switch lower(arg.type)
    case 'anovan'
    case 'rank'
        x = reshape(tiedrank(x),nobs,v);
    otherwise
        error('The TYPE parameter value must be ANOVAN, or RANK.')
end

% Construct Type III effect-coded design matrices
Xfactors = cell(1,nfactors);
for i = 1:nfactors
    k = max(g(:,i));
    Xg = zeros(nobs,k-1);
    for j = 1:k-1
        Xg(g(:,i)==j,j) = 1;
    end
    Xg(g(:,i)==k,:) = -1;
    Xfactors{i} = Xg;
end

% Define full and reduced models
XF = ones(nobs,1);
for i = 1:nfactors
    XF = [XF,Xfactors{i}];
end
Xr = cell(1,nfactors);
for i = 1:nfactors
    Xred = ones(nobs,1);
    for j = 1:nfactors
        if i ~= j
            Xred = [Xred,Xfactors{j}];
        end
    end
    Xr{i} = Xred;
end

% Precompute exact projection matrices for Freedman-Lane
I = eye(nobs);
MF = I-XF*pinv(XF);
Mr = cell(1,nfactors);
for i = 1:nfactors
    Mr{i} = I-Xr{i}*pinv(Xr{i});
end

% Compute degrees of freedom
rankF = rank(XF);
dfF = zeros(1,nfactors);
for i = 1:nfactors
    dfF(i) = rankF-rank(Xr{i});
end
dfE = nobs-rankF;

% Preallocate memory
f = zeros(nfactors,v);
R = cell(1,nfactors);
SSfactors = zeros(nfactors,v);

% Compute observed statistics
SSE = sum((MF*x).^2,1);
MSE = SSE/dfE;
for i = 1:nfactors
    SSred = sum((Mr{i}*x).^2,1);
    SSfactors(i,:) = SSred-SSE;
    MSfactor = SSfactors(i,:)/dfF(i);
    f(i,:) = MSfactor./MSE;
    R{i} = Mr{i}*x;
end

if nargout > 1

    rng(arg.seed);

    % Estimate sampling distribution
    dist = zeros(arg.nperm,nfactors,v);
    for i = 1:arg.nperm
        randidx = randperm(nobs);
        for j = 1:nfactors
            Rp = R{j}(randidx,:);
            SSEp = sum((MF*Rp).^2,1);
            SSredp = sum((Mr{j}*Rp).^2,1);
            SSeffectp = SSredp-SSEp;
            dist(i,j,:) = (SSeffectp/dfF(j))./(SSEp/dfE);
        end
    end

    % Apply max correction if specified
    if arg.correct
        distmax = max(dist,[],3);
    else
        distmax = dist;
    end

    % Compute p-values
    p = zeros(nfactors,v);
    for j = 1:nfactors
        p(j,:) = (sum(f(j,:)<=distmax(:,j),1)+1)/(arg.nperm+1);
    end

end

% Compute confidence intervals
if nargout > 2
    ci = zeros(2,nfactors,v);
    crit = zeros(nfactors,v);
    for j = 1:nfactors
        if arg.correct
            cf = prctile(distmax(:,j),100*(1-arg.alpha));
            crit(j,:) = cf*ones(1,v);
        else
            crit(j,:) = prctile(dist(:,j,:),100*(1-arg.alpha),1);
        end
    end
    ci(1,:,:) = f./crit;
    ci(2,:,:) = Inf;
end

% Store statistics in a structure
if nargout > 3
    stats.source = 'permuanovan';
    stats.gnames = gnames;
    stats.df = [dfF,dfE];
    stats.method = arg.type;
end

% Create ANOVA table
if nargout > 4
    if v == 1
        tbl = cell(nfactors+2,6);
        tbl(1,:) = {'Source','SS','df','MS','F','Prob>F'};
        for i = 1:nfactors
            tbl(i+1,:) = {['Factor ',num2str(i)],SSfactors(i),...
                dfF(i),SSfactors(i)/dfF(i),f(i),p(i)};
        end
        tbl(nfactors+2,:) = {'Error',SSE,dfE,MSE,[],[]};
    else
        tbl = cell(v*nfactors+1,6);
        tbl(1,:) = {'Source','SS','df','MS','F','Prob>F'};
        row = 2;
        for kidx = 1:v
            for i = 1:nfactors
                tbl(row,:) = {['Var ',num2str(kidx),' - Factor ',...
                    num2str(i)],SSfactors(i,kidx),dfF(i),...
                    SSfactors(i,kidx)/dfF(i),f(i,kidx),p(i,kidx)};
                row = row+1;
            end
        end
    end
end