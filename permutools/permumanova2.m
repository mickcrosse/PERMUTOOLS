function [lambda,p,ci,stats,tbl,dist] = permumanova2(x,group,varargin)
%PERMUMANOVA2  Permutation-based two-way multivariate ANOVA and rank MANOVA.
%   LAMBDA = PERMUMANOVA2(X,GROUP) performs a permutation-based two-way
%   multivariate analysis of variance (MANOVA) for comparing the means of
%   data in matrix X, and returns the test statistics. X is an N-by-P
%   matrix where N is observations and P is dependent variables. GROUP is
%   an N-by-2 grouping matrix or cell array containing two categorical
%   factors. If X is a 3D matrix (N-by-P-by-V), the test evaluates all V
%   sets of variables simultaneously, and a matrix of results is returned.
%
%   For non-normally distributed data, the raw data may be transformed to
%   rank orders in order to compute a rank MANOVA test by setting the
%   'type' parameter to 'rank'.
%
%   PERMUMANOVA2 treats NaNs as missing values, and ignores them.
%
%   [LAMBDA,P] = PERMUMANOVA2(...) returns the probability (i.e. p-value)
%   of observing the given result by chance if the null hypothesis is true.
%   As the null distribution is generated empirically by permuting the
%   data, no assumption is made about the shape of the distributions that
%   the data come from. P-values are automatically adjusted for multiple
%   comparisons using the max correction method.
%
%   [LAMBDA,P,CI] = PERMUMANOVA2(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for Wilks' Lambda. CIs are also adjusted for multiple
%   comparisons using the max correction method.
%
%   [LAMBDA,P,CI,STATS] = PERMUMANOVA2(...) returns a structure with the
%   following fields:
%       'source'    -- the function used to compute the MANOVA
%       'method'    -- the statistical method used
%       'gnames'    -- the group names for both factors
%       'df'        -- the degrees of freedom [A, B, AxB, Error]
%
%   [LAMBDA,P,CI,STATS,TBL] = PERMUMANOVA2(...) returns the MANOVA table
%   contents as a cell array.
%
%   [LAMBDA,P,CI,STATS,TBL,DIST] = PERMUMANOVA2(...) returns the permuted
%   sampling distribution of the test statistics.
%
%   [...] = PERMUMANOVA2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'manova2'       two-way MANOVA (default)
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
%   See also MANOVA1 PERMUMANOVA1 PERMUANOVA2 PERMUREGRESS.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(group) || size(group,2) < 2
    error('A grouping matrix with at least 2 columns is required.')
end

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'manova2';
end

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Parse grouping variables and handle missing values
[g1,gnames1] = grp2idx(group(:,1));
[g2,gnames2] = grp2idx(group(:,2));
nanidx = any(isnan(x),2) | isnan(g1) | isnan(g2);
if ndims(x) == 3
    nanidx = any(nanidx,3);
end
if any(nanidx)
    x(nanidx,:,:) = [];
    g1(nanidx) = [];
    g2(nanidx) = [];
end

% Get data dimensions
[nobs,pnum,v] = size(x);

% Transform raw data to rank orders if specified
switch lower(arg.type)
    case 'manova2'
    case 'rank'
        x = reshape(tiedrank(reshape(x,nobs,[])),nobs,pnum,v);
    otherwise
        error('The TYPE parameter value must be MANOVA2, or RANK.')
end

% Construct effect-coded design matrices
k1 = max(g1);
k2 = max(g2);
X1 = zeros(nobs,k1-1);
X2 = zeros(nobs,k2-1);
for i = 1:k1-1
    X1(g1==i,i) = 1;
end
X1(g1==k1,:) = -1;
for i = 1:k2-1
    X2(g2==i,i) = 1;
end
X2(g2==k2,:) = -1;
X3 = zeros(nobs,(k1-1)*(k2-1));
col = 1;
for i = 1:k1-1
    for j = 1:k2-1
        X3(:,col) = X1(:,i).*X2(:,j);
        col = col+1;
    end
end

% Define full and reduced models (Type III SS)
XF = [ones(nobs,1),X1,X2,X3];
Xr1 = [ones(nobs,1),X2,X3];
Xr2 = [ones(nobs,1),X1,X3];
Xr3 = [ones(nobs,1),X1,X2];

% Precompute exact projection matrices for Freedman-Lane
I = eye(nobs);
MF = I-XF*pinv(XF);
Mr1 = I-Xr1*pinv(Xr1);
Mr2 = I-Xr2*pinv(Xr2);
Mr3 = I-Xr3*pinv(Xr3);

% Compute degrees of freedom
rankF = rank(XF);
dfA = rankF-rank(Xr1);
dfB = rankF-rank(Xr2);
dfAB = rankF-rank(Xr3);
dfE = nobs-rankF;

% Preallocate memory
lambda = zeros(3,v);
R1 = zeros(nobs,pnum,v);
R2 = zeros(nobs,pnum,v);
R3 = zeros(nobs,pnum,v);

% Compute true Wilks' Lambda and isolate residuals
for kidx = 1:v
    xk = x(:,:,kidx);
    % Total error SSCP (Full model)
    WT = xk'*MF*xk;
    detWT = det(WT);
    % Factor A
    WA = xk'*Mr1*xk;
    lambda(1,kidx) = detWT/det(WA);
    R1(:,:,kidx) = Mr1*xk;
    % Factor B
    WB = xk'*Mr2*xk;
    lambda(2,kidx) = detWT/det(WB);
    R2(:,:,kidx) = Mr2*xk;
    % Interaction AB
    WAB = xk'*Mr3*xk;
    lambda(3,kidx) = detWT/det(WAB);
    R3(:,:,kidx) = Mr3*xk;
end

if nargout > 1

    % Generate random permutations
    rng(arg.seed);

    % Estimate sampling distribution
    dist = zeros(arg.nperm,3,v);
    for i = 1:arg.nperm
        randidx = randperm(nobs);
        for kidx = 1:v
            % Factor A permutation
            R1p = R1(randidx,:,kidx);
            dist(i,1,kidx) = det(R1p'*MF*R1p)/det(R1p'*Mr1*R1p);
            % Factor B permutation
            R2p = R2(randidx,:,kidx);
            dist(i,2,kidx) = det(R2p'*MF*R2p)/det(R2p'*Mr2*R2p);
            % Interaction AB permutation
            R3p = R3(randidx,:,kidx);
            dist(i,3,kidx) = det(R3p'*MF*R3p)/det(R3p'*Mr3*R3p);
        end
    end

    % Apply max-correction if specified
    if arg.correct
        distmin = min(dist,[],3);
    else
        distmin = dist;
    end

    % Compute p-values
    p = zeros(3,v);
    for fidx = 1:3
        p(fidx,:) = (sum(lambda(fidx,:)>=distmin(:,fidx),1)+1)/(arg.nperm+1);
    end

end

% Compute confidence interval
if nargout > 2
    ci = zeros(2,3,v);
    crit = zeros(3,v);
    for fidx = 1:3
        if arg.correct
            cf = prctile(distmin(:,fidx),100*arg.alpha);
            crit(fidx,:) = cf*ones(1,v);
        else
            crit(fidx,:) = prctile(dist(:,fidx,:),100*arg.alpha,1);
        end
    end
    ci(2,:,:) = lambda./crit;
end

% Store statistics in a structure
if nargout > 3
    stats.source = 'permumanova2';
    stats.method = arg.type;
    stats.gnames = {gnames1,gnames2};
    stats.df = [dfA,dfB,dfAB,dfE];
end

% Create MANOVA table
if nargout > 4
    if v == 1
        tbl = {
            'Source','Wilks Lambda','df','Prob>Lambda';
            'Factor A',lambda(1),dfA,p(1);
            'Factor B',lambda(2),dfB,p(2);
            'A x B',lambda(3),dfAB,p(3)};
    else
        tbl = cell(v*3+1,4);
        tbl(1,:) = {'Source','Wilks Lambda','df','Prob>Lambda'};
        row = 2;
        for kidx = 1:v
            tbl(row,:) = {['Var ',num2str(kidx),' - Factor A'],...
                lambda(1,kidx),dfA,p(1,kidx)};
            tbl(row+1,:) = {['Var ',num2str(kidx),' - Factor B'],...
                lambda(2,kidx),dfB,p(2,kidx)};
            tbl(row+2,:) = {['Var ',num2str(kidx),' - A x B'],...
                lambda(3,kidx),dfAB,p(3,kidx)};
            row = row+3;
        end
    end
end