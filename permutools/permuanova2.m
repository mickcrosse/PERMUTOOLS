function [f,p,ci,stats,tbl,dist] = permuanova2(x,reps,varargin)
%PERMUANOVA2  Permutation-based two-way ANOVA and aligned rank transform test.
%   F = PERMUANOVA2(X) performs a balanced permutation-based two-way
%   analysis of variance (ANOVA) for comparing the means of two or more
%   columns and two or more rows of data in matrix X, and returns the test
%   statistics for the columns, rows and interactions, respectively.
%
%   For non-normally distributed data, the raw data may be transformed
%   to rank orders in order to compute an aligned rank transform (ART) test
%   by setting the 'type' parameter to 'alignrank' or 'rank'.
%
%   PERMUANOVA2 requires a fully balanced design and does not support
%   missing values (NaNs).
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
%   intervals (CIs) for the true difference of population means.
%
%   [F,P,CI,STATS] = PERMUANOVA2(...) returns a structure with the
%   following fields:
%       'source'    -- the function used to compute the ANOVA
%       'method'    -- the statistical method used
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
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'anova2'    two-way ANOVA (default)
%                       'alignrank' aligned rank transform (ART) test
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

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(reps)
    reps = 1;
end

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'anova2';
end

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Orient data column-wise
if arg.dim==2
    x = x';
end

% Balanced design check
if any(isnan(x(:)))
    error('Two-way ANOVA requires a perfectly balanced design. Data must not contain NaNs.')
end

% Get data dimensions
shapex = size(x);
nobs = numel(x);
[rows,cols] = size(x);
if reps > 1
    grows = rows/reps;
else
    grows = rows;
end

% Aligned rank transform (ART) data prep
switch arg.type
    case 'anova2'
        % Use same raw data for all three effects
        xc = x; xr = x; xi = x;
    case {'alignrank', 'rank'}
        % Calculate cell, column, and row means
        gm = mean(x,'all');
        MC = repmat(mean(x,1),rows,1);
        if reps > 1
            MR = zeros(rows,cols);
            MCR = zeros(rows,cols);
            for i = 1:grows
                idx = reps*(i-1)+(1:reps);
                MR(idx,:) = mean(x(idx,:),'all');
                MCR(idx,:) = repmat(mean(x(idx,:),1),reps,1);
            end
        else
            MR = repmat(mean(x,2),1,cols);
            MCR = x;
        end
        % Align and rank data
        YC = x-MCR+MC;
        YR = x-MCR+MR;
        YCR = x-MC-MR+2*gm;
        xc = reshape(tiedrank(YC(:)),rows,cols);
        xr = reshape(tiedrank(YR(:)),rows,cols);
        xi = reshape(tiedrank(YCR(:)),rows,cols);
    otherwise
        error('The TYPE parameter value must be ANOVA2, or ALIGNRANK.')
end

% Stack datasets into 3D array
X3 = cat(3,xc,xr,xi);

% Compute block means for 3D array if reps > 1
if reps > 1
    X3m = zeros(grows,cols,3);
    for i = 1:grows
        idx = reps*(i-1)+(1:reps);
        X3m(i,:,:) = mean(X3(idx,:,:),1);
    end
else
    X3m = X3;
end

% Compute group sample sizes
coln = reps*grows;
rown = reps*cols;

% Compute global 3D means and factors
gm3 = sum(X3m,[1,2])/nobs;
factor3 = nobs*gm3.^2;

% Compute sum of squares across 3rd dimension
colmeans3 = mean(X3m,1);
css3 = coln*sum((colmeans3-gm3).^2,[1,2]);
rowmeans3 = mean(X3m,2);
rss3 = rown*sum((rowmeans3-gm3).^2,[1,2]);
iss3 = reps*sum(X3m.^2,[1,2])-css3-rss3-factor3;
tss3 = sum(X3.^2,[1,2])-factor3;

% Compute error SS
if reps > 1
    ess3 = tss3-css3-rss3-iss3;
else
    ess3 = iss3;
end

% Extract relevant ART statistics
css = css3(1);
rss = rss3(2);
iss = iss3(3);
ess = ess3(3);

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

% Compute mean squares using effect-specific error SS
msc = css/dfc;
msr = rss/dfr;
msi = iss/dfi;
msec = ess3(1)/dfe;
mser = ess3(2)/dfe;
msei = ess3(3)/dfe;

% Compute F-statistics
fc = msc/msec;
fr = msr/mser;
fi = msi/msei;
if reps > 1
    f = [fc,fr,fi];
else
    f = [fc,fr];
end

if nargout > 1

    % Concatenate groups
    xc = xc(:);
    xr = xr(:);
    xi = xi(:);

    % Generate random permutations
    rng(arg.seed);
    [~,idx] = sort(rand(nobs,arg.nperm));

    % Estimate sampling distributions
    distc = zeros(arg.nperm,1);
    distr = zeros(arg.nperm,1);
    if reps > 1
        disti = zeros(arg.nperm,1);
    end

    for i = 1:arg.nperm

        % Apply permutation index to all 3 datasets
        xcp = xc(idx(:,i));
        xrp = xr(idx(:,i));
        xip = xi(idx(:,i));

        % Reshape and stack into 3D array
        X3p = cat(3,reshape(xcp,shapex),reshape(xrp,shapex),...
            reshape(xip,shapex));

        % Compute block means for 3D array
        if reps > 1
            X3mp = zeros(grows,cols,3);
            for j = 1:grows
                idxp = reps*(j-1)+(1:reps);
                X3mp(j,:,:) = mean(X3p(idxp,:,:),1);
            end
        else
            X3mp = X3p;
        end

        % Compute global 3D means and correction factor
        gm3p = sum(X3mp,[1,2])/nobs;
        factor3p = nobs*gm3p.^2;

        % Compute sum of squares across 3rd dimension
        colmeans3p = mean(X3mp,1);
        css3p = coln*sum((colmeans3p-gm3p).^2,[1,2]);
        rowmeans3p = mean(X3mp,2);
        rss3p = rown*sum((rowmeans3p-gm3p).^2,[1,2]);
        iss3p = reps*sum(X3mp.^2,[1,2])-css3p-rss3p-factor3p;
        if reps > 1
            tss3p = sum(X3p.^2,[1,2])-factor3p;
            ess3p = tss3p-css3p-rss3p-iss3p;
        else
            ess3p = iss3p;
        end

        % Compute permuted mean squares using effect-specific error SS
        msecp = ess3p(1)/dfe;
        mserp = ess3p(2)/dfe;

        % Store permuted F-statistics
        distc(i) = (css3p(1)/dfc)/msecp;
        distr(i) = (rss3p(2)/dfr)/mserp;
        if reps > 1
            mse_ip = ess3p(3)/dfe;
            disti(i) = (iss3p(3)/dfi)/mse_ip;
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
end

% Compute confidence intervals
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

% Store statistics in a structure
if nargout > 3
    stats.source = 'permuanova2';
    stats.method = arg.type;
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
    mse = ess/dfe;
    tss = css+rss+iss+ess;
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