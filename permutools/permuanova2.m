function [f,p,ci,stats,tbl,dist] = permuanova2(x,reps,varargin)
%PERMUANOVA2  Two-way permutation-based analysis of variance (ANOVA).
%   F = PERMUANOVA2(X) performs a balanced two-way permutation-based ANOVA
%   for comparing the means of two or more columns and two or more rows of
%   data in matrix X, and returns the test statistics for the columns, rows
%   and interactions (if any), respectively.
%
%   PERMUANOVA2 treats NaNs as missing values, and ignores them.
%
%   F = PERMUANOVA2(X,REPS) groups the rows of X according to the number of
%   replicates REPS for each combination of factor groups. To test for an
%   interaction effect, REPS must be greater than 1. REPS must be a scalar.
%
%   [F,P] = PERMUANOVA2(...) returns the probabilities (i.e. p-values) of
%   observing the given results by chance if the null hypothesis (that the
%   means of the groups are equal) is true. As the null distribution is
%   generated empirically by permuting the data, no assumption is made
%   about the shape of the distributions that the data come from, except
%   that they have equal variances.
%
%   [F,P,CI] = PERMUANOVA2(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for the true difference of population means.
%
%   [F,P,CI,STATS] = PERMUANOVA2(...) returns a structure with the
%   following fields:
%       'source'    -- the function used to compute the ANOVA
%       'sigmasq'   -- the error mean square
%       'colmeans'  -- the column means
%       'coln'      -- the column sample sizes
%       'rowmeans'  -- the row means
%       'rown'      -- the row sample sizes
%       'inter'     -- the inclusion of an interaction term
%       'pval'      -- the interaction p-value
%       'df'        -- the error degrees of freedom
%
%   [F,P,CI,STATS,TBL] = PERMUANOVA2(...) returns the ANOVA table contents
%   as a cell array.
%
%   [F,P,CI,STATS,TBL,PDIST] = PERMUANOVA2(...) returns the permuted
%   sampling distributions of the test statistics.
%
%   [...] = PERMUANOVA2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator. By default, the
%                   generator is initialised based on the current time,
%                   resulting in a different permutation on each call.
%
%   See also ANOVA2 PERMUANOVA1 PERMUTTEST2 BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(reps)
    reps = 1;
end

% Parse input arguments
arg = ptparsevarargin(varargin);

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Orient data column-wise
if arg.dim==2
    x = x';
end

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Get data dimensions, ignoring NaNs
shapex = size(x);
nobs = sum(~isnan(x),'all');
[rows,cols] = size(x);
if reps > 1
    grows = rows/reps;
    xm = zeros(grows,cols);
    for i = 1:grows
        idx = reps*(i-1);
        xm(i,:) = mean(x(idx+1:idx+reps,:));
    end
else
    grows = rows;
    xm = x;
end

% Compute group sample sizes
coln = reps*grows;
rown = reps*cols;

% Compute grand mean
gm = sum(xm,'all',nanflag)/nobs;

% Compute column sum of squares
colmeans = mean(xm,1,nanflag);
css = coln*sum((colmeans-gm).^2,'all',nanflag);

% Compute row sum of squares
rowmeans = mean(xm,2,nanflag);
rss = rown*sum((rowmeans-gm).^2,'all',nanflag);

% Compute interaction sum of squares
factor = reps*cols*grows*gm^2;
iss = reps*sum(xm.^2,'all',nanflag)-css-rss-factor;

% Compute total sum of squares
tss = sum(x.^2,'all',nanflag)-factor;

% Compute error sum of squares
if reps > 1
    ess = tss-css-rss-iss;
else
    ess = iss;
end

% Compute degrees of freedom
dft = nobs-1;
dfc = cols-1;
dfr = grows-1;
if reps > 1
    dfe = (reps-1)*cols*grows;
else
    dfe = dfc*(rows-1);
end
dfi = dfc*dfr;

% Compute mean squares
msc = css/dfc;
msr = rss/dfr;
msi = iss/dfi;
mse = ess/dfe;

% Compute F-statistics
fc = msc/mse;
fr = msr/mse;
fi = msi/mse;
if reps > 1
    f = [fc,fr,fi];
else
    f = [fc,fr];
end

if nargout > 1

    % Concatenate groups
    x = x(:);

    % Generate random permutations
    rng(arg.seed);
    maxnobs = numel(x);
    [~,idx] = sort(rand(maxnobs,arg.nperm));

    % Estimate sampling distributions
    distc = zeros(arg.nperm,1);
    distr = zeros(arg.nperm,1);
    if reps > 1
        disti = zeros(arg.nperm,1);
    end
    for i = 1:arg.nperm
        xp = x(idx(:,i));
        xp = reshape(xp,shapex);
        xm = zeros(grows,cols);
        for j = 1:grows
            idxp = reps*(j-1);
            xm(j,:) = mean(xp(idxp+1:idxp+reps,:));
        end
        colmeansp = mean(xm,1,nanflag);
        cssp = coln*sum((colmeansp-gm).^2,'all',nanflag);
        rowmeansp = mean(xm,2,nanflag);
        rssp = rown*sum((rowmeansp-gm).^2,'all',nanflag);
        issp = reps*sum(xm.^2,'all',nanflag)-cssp-rssp-factor;
        if reps > 1
            tssp = sum(xp.^2,'all',nanflag)-factor;
            essp = tssp-cssp-rssp-issp;
        else
            essp = issp;
        end
        msep = essp/dfe;
        distc(i) = cssp/dfc/msep;
        distr(i) = rssp/dfr/msep;
        if reps > 1
            disti(i) = issp/dfi/msep;
        end
    end

    % Compute p-values
    pc = (sum(fc<=distc)+1)/(arg.nperm+1);
    pr = (sum(fr<=distr)+1)/(arg.nperm+1);
    if reps > 1
        pint = (sum(fi<=disti)+1)/(arg.nperm+1);
        p = [pc,pr,pint];
    else
        pint = NaN;
        p = [pc,pr];
    end

    % Compute CIs
    if nargout > 2
        crit = prctile(distc,100*(1-arg.alpha));
        cic = [fc./crit;Inf];
        crit = prctile(distr,100*(1-arg.alpha));
        cir = [fr./crit;Inf];
        if reps > 1
            crit = prctile(disti,100*(1-arg.alpha));
            cii = [fi./crit;Inf];
            ci = [cic,cir,cii];
        else
            ci = [cic,cir];
        end
    end

end

% Store statistics in a structure
if nargout > 3
    stats.source = 'permuanova2';
    stats.sigmasq = mse;
    stats.colmeans = colmeans;
    stats.coln = coln;
    stats.rowmeans = rowmeans';
    stats.rown = rown;
    stats.inter = reps>1;
    stats.pval = pint;
    stats.df = dfe;
end

% Create ANOVA table
if nargout > 4
    tbl = {
        'Source','SS','df','MS','F','Prob>F';
        'Columns',css,dfc,msc,fc,pc;
        'Rows',rss,dfr,msr,fr,pr;
        'Interaction',iss,dfi,msi,fi,pint;
        'Error',ess,dfe,mse,[],[];
        'Total',tss,dft,[],[],[]};
    if reps == 1
        tbl(4,:) = [];
    end
end

% Arrange permutation distributions
if nargout > 5
    if reps > 1
        dist = [distc,distr,disti];
    else
        dist = [distc,distr];
    end
end