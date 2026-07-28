function [d,ci,stats,dist] = booteffectsize(x,m,varargin)
%BOOTEFFECTSIZE  Effect size with bootstrapped confidence intervals.
%   D = BOOTEFFECTSIZE(X) returns the effect size measure for a single
%   sample X based on Cohen's d. By default, Cohen's d is bias-corrected
%   for sample size, also known as Hedges' g. For ordinal data, Cliff's
%   delta can be computed by setting the 'effect' parameter to 'cliff'. If
%   X is a matrix, separate effect sizes are measured along each column of
%   X, and a vector of results is returned. If the 'compare' parameter is
%   set to 'pairwise', the effect sizes between every pair of columns in X
%   are measured, and a matrix of results is returned.
%
%   BOOTEFFECTSIZE treats NaNs as missing values, and ignores them.
%
%   D = BOOTEFFECTSIZE(X,M) returns the effect size measure for a single
%   sample X with a known mean M. M must be a scalar.
%
%   D = BOOTEFFECTSIZE(X,Y) returns the effect size between two dependent
%   samples X and Y using the pooled standard deviation. X and Y can be
%   treated as independent samples by setting the 'paired' parameter to 0.
%   If X and Y are independent samples with significantly different
%   variances, an estimate based on the control sample's variance (Glass'
%   delta) can be computed by setting the 'effect' parameter to 'glass'.
%   For this, the control sample should be entered as X, and the test
%   sample as Y.
%
%   [D,CI] = BOOTEFFECTSIZE(...) returns the bootstrapped, bias-corrected
%   confidence intervals using the percentile method.
%
%   [D,CI,STATS] = BOOTEFFECTSIZE(...) returns a structure with the
%   following fields:
%       'method'    -- the effect size metric used
%       'df'        -- the degrees of freedom of each measure
%       'sd'        -- the pooled standard deviation, or of X for a one-
%                      sample or Glass' delta measure
%
%   [D,CI,STATS,DIST] = BOOTEFFECTSIZE(...) returns the bootstrapped
%   sampling distribution of the effect size measure.
%
%   [...] = BOOTEFFECTSIZE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'alpha'     A scalar between 0 and 1 specifying the confidence
%                   level as 100*(1-ALPHA)% (default=0.05 for 95%
%                   confidence).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'paired'    A numeric scalar (0,1) or logical indicating whether
%                   the data in X and Y are paired: pass in 1 for paired
%                   samples (default), or 0 for unpaired samples.
%       'effect'    A string specifying the effect size metric to compute:
%                       'cohen'      standardised mean difference based on
%                                    Cohen's d (default)
%                       'glass'      standardised mean difference based on
%                                    Glass' delta for comparing independent
%                                    samples with significantly different
%                                    variances
%                       'cliff'      unstandardised mean difference based
%                                    on Cliff's delta for comparing ordinal
%                                    data
%                       'meandiff'   unstandardised mean difference
%                       'mediandiff' unstandardised median difference
%       'vartype'   A string specifying the variance equivalence of X and Y
%                   to determine the SD and degrees of freedom:
%                       'equal'   	assume equal variances (default)
%                       'unequal' 	assume unequal variances
%       'compare'   A string specifying what to compare each variable to
%                   when only X is entered:
%                       'zero'      compare each column of X to zero and
%                                   return a vector of results (default)
%                       'pairwise'  compare every pair of columns in X to
%                                   each other and return a matrix of
%                                   results
%       'nboot'     An integer scalar specifying the number of bootstraps
%                   used to estimate the confidence intervals (default=
%                   10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   bias-correct the effect size and confidence intervals
%                   according to sample size (default=true). Note, this
%                   only applies to standardised effect size measures.
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'seed'      A scalar integer specifying the seed used to initialise
%                   the bootstrap generator. By default, the generator is
%                   initialised based on the current time, resulting in a
%                   different bootstrap on each call.
%
%   See also MEANEFFECTSIZE BOOTCI PERMUTTEST PERMUTTEST2 PERMUVARTEST2.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%       [2] Hentschke H, Stuttgen MC (2011) Computation of measures of
%           effect size for neuroscience data sets. Eur J Neurosci,
%           34:1887–1894.
%       [3] Cohen J (1969) Statistical power for the behavioural sciences.
%           London: Academic Press.
%       [4] Hedges LV, Olkin I (1985) Statistical methods for meta-
%           analysis. San Diego, CA: Academic Press.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(m)
    y = [];
elseif isscalar(m)
    y = [];
    x = x-m;
elseif ismatrix(m)
    y = m;
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

% Set up comparison
if isempty(y)
    switch arg.compare
        case 'zero'
            y = zeros(size(x));
        case 'pairwise'
            warning('Comparing all columns of X...')
            [x,y] = ptpaircols(x);
            arg.mat = true;
    end
else
    switch arg.compare
        case 'pairwise'
            error('The PAIRWISE option can only be used with one sample.')
    end
end
if size(x,2)~=size(y,2)
    error('X and Y must have the same number of variables.')
end

% Get data dimensions, ignoring NaNs
nvar = size(x,2);
nobsx = sum(~isnan(x));
nobsy = sum(~isnan(y));

% Compute degrees of freedom
dfx = nobsx-1;
dfy = nobsy-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:))) || any(isnan(y(:)))
    nanflag = 'omitnan';
else
    nanflag = 'includenan';
end

% Compute sample variance using fast algo
sumx = sum(x,nanflag);
sumy = sum(y,nanflag);
varx = (sum(x.^2,nanflag)-(sumx.^2)./nobsx)./dfx;
vary = (sum(y.^2,nanflag)-(sumy.^2)./nobsy)./dfy;

if arg.paired

    % Check input parameters
    switch arg.effect
        case 'glass'
            error('GLASS can only be used for independent samples.')
    end
    switch arg.vartype
        case 'unequal'
            error('Cannot assume unequal variance for dependent samples.')
    end

    % Compute difference between samples
    switch arg.effect
        case 'cliff'
            xdy = zeros(max(nobsx)^2,nvar);
            for i = 1:nvar
                signxdy = sign(x(:,i)-y(:,i)');
                xdy(:,i) = signxdy(:);
            end
            xdy(1:max(nobsx)+1:end,:) = 0;
        otherwise
            xdy = x-y;
    end

    % Use only rows with no NaN values if specified
    switch arg.rows
        case 'complete'
            xdy = xdy(~any(isnan(xdy),2),:);
    end

    % Get data dimensions, ignoring NaNs
    switch arg.effect
        case 'cliff'
            nobs = nobsx.*dfx;
        otherwise
            nobs = sum(~isnan(xdy));
    end

    % Compute mean difference
    switch arg.effect
        case 'mediandiff'
            mu = median(xdy,nanflag);
        otherwise
            mu = sum(xdy,nanflag)./nobs;
    end

    % Compute standard deviation
    df = nobs-1;
    sd = sqrt((varx+vary)/2);

else

    % Use only rows with no NaN values if specified
    switch arg.rows
        case 'complete'
            x = x(~any(isnan(x),2),:);
            y = y(~any(isnan(y),2),:);
    end

    % Compute difference between samples
    switch arg.effect
        case 'cliff'
            xdy = zeros(max(nobsx)*max(nobsy),nvar);
            for i = 1:nvar
                signxdy = sign(x(:,i)-y(:,i)');
                xdy(:,i) = signxdy(:);
            end
        otherwise
            xdy = [x;y];
    end

    % Get data dimensions, ignoring NaNs
    nobs = sum(~isnan(xdy));

    % Compute mean difference
    switch arg.effect
        case 'cliff'
            mu = sum(xdy,nanflag)./nobs;
        case 'mediandiff'
            mu = median(x,nanflag)-median(y,nanflag);
        otherwise
            mu = sumx./nobsx-sumy./nobsy;
    end

    % Compute standard deviation
    switch arg.effect
        case 'glass'
            switch arg.vartype
                case 'equal'
                    warning(['GLASS option should only be used for '...
                        'samples with unequal variances.'])
            end
            df = dfx;
            sd = sqrt(varx);
        otherwise
            switch arg.vartype
                case 'equal'
                    df = nobs-2;
                    sd = sqrt((dfx.*varx+dfy.*vary)./df);
                case 'unequal'
                    df = (dfx.*dfy.*(varx+vary).^2)./(dfy.*varx.^2+...
                        dfx.*vary.^2);
                    sd = sqrt((varx+vary)/2);
            end
    end

end

% Compute effect size
switch arg.effect
    case {'cliff','meandiff','mediandiff'}
        d = mu;
    otherwise
        d = mu./sd;
end

if nargout > 1

    rng(arg.seed);
    ci = zeros(2,nvar);
    dist = zeros(arg.nboot,nvar);

    for i = 1:nvar

        % Isolate valid data to eliminate NaNs
        if arg.paired
            valid = ~isnan(x(:,i))&~isnan(y(:,i));
            xi = x(valid,i);
            yi = y(valid,i);
            ni = sum(valid);
            idx = randi(ni,[ni,arg.nboot],'uint32');
            xb = xi(idx);
            yb = yi(idx);
            nxb = ni; nyb = ni;
            dfxb = ni-1; dfyb = ni-1; dfb = ni-1;
        else
            validx = ~isnan(x(:,i));
            xi = x(validx,i);
            nxb = sum(validx);
            idxxi = randi(nxb,[nxb,arg.nboot],'uint32');
            xb = xi(idxxi);
            dfxb = nxb-1;
            validy = ~isnan(y(:,i));
            yi = y(validy,i);
            nyb = sum(validy);
            idxyi = randi(nyb,[nyb,arg.nboot],'uint32');
            yb = yi(idxyi);
            dfyb = nyb-1;
            dfb = nxb+nyb-2;
        end

        % Estimate sampling distribution
        sumxb = sum(xb,1);
        sumyb = sum(yb,1);
        varxb = (sum(xb.^2,1)-(sumxb.^2)./nxb)./dfxb;
        varyb = (sum(yb.^2,1)-(sumyb.^2)./nyb)./dfyb;

        if arg.paired
            switch arg.effect
                case 'cliff'
                    try  % fast algo
                        xb3d = reshape(xb,nxb,1,arg.nboot);
                        yb3d = reshape(yb,1,nxb,arg.nboot);
                        signxdy = sign(xb3d-yb3d);
                        clear xb3d yb3d;
                        xdyb = reshape(signxdy,nxb^2,arg.nboot);
                        clear signxdy;
                        xdyb(1:nxb+1:end,:) = 0;
                    catch  % memory-efficient algo
                        xdyb = zeros(nxb^2,arg.nboot);
                        for j = 1:arg.nboot
                            signxdy = sign(xb(:,j)-yb(:,j)');
                            xdyb(:,j) = signxdy(:);
                        end
                        xdyb(1:nxb+1:end,:) = 0;
                    end
                otherwise
                    xdyb = xb-yb;
            end
            switch arg.effect
                case 'mediandiff'
                    mub = median(xdyb,1);
                otherwise
                    mub = sum(xdyb,1)./nxb;
            end
            sdb = sqrt((varxb+varyb)/2);
        else
            switch arg.effect
                case 'cliff'
                    try  % fast algo
                        xb3d = reshape(xb,nxb,1,arg.nboot);
                        yb3d = reshape(yb,1,nyb,arg.nboot);
                        signxdy = sign(xb3d-yb3d);
                        clear xb3d yb3d;
                        xdyb = reshape(signxdy,nxb*nyb,arg.nboot);
                        clear signxdy;
                    catch  % memory-efficient algo
                        xdyb = zeros(nxb*nyb,arg.nboot);
                        for j = 1:arg.nboot
                            signxdy = sign(xb(:,j)-yb(:,j)');
                            xdyb(:,j) = signxdy(:);
                        end
                    end
            end
            switch arg.effect
                case 'cliff'
                    mub = sum(xdyb,1)./(nxb*nyb);
                case 'mediandiff'
                    mub = median(xb,1)-median(yb,1);
                otherwise
                    mub = sumxb./nxb-sumyb./nyb;
            end
            switch arg.effect
                case 'glass'
                    sdb = sqrt(varxb);
                case 'cohen'
                    switch arg.vartype
                        case 'equal'
                            sdb = sqrt((dfxb.*varxb+dfyb.*varyb)./dfb);
                        case 'unequal'
                            sdb = sqrt((varxb+varyb)/2);
                    end
            end
        end
        switch arg.effect
            case {'cliff','meandiff','mediandiff'}
                dist(:,i) = mub;
            otherwise
                dist(:,i) = mub./sdb;
        end

        % Compute confidence interval
        ci(:,i) = prctile(dist(:,i),100*[arg.alpha/2;1-arg.alpha/2]);

    end

end

% Store statistics in a structure
if nargout > 2
    stats.method = arg.effect;
    stats.df = df;
    stats.sd = sd;
end

% Bias-correct standardised measures
if arg.correct
    switch arg.effect
        case {'cohen','glass'}
            factor = exp(gammaln(df/2)-log(sqrt(df/2))-gammaln((df-1)/2));
            d = factor.*d;
            if nargout > 1
                ci = [factor;factor].*ci;
            end
            if nargout > 3
                dist = factor.*dist;
            end
    end
end

% Arrange results in a matrix if specified
if arg.mat
    d = ptvec2mat(d);
    if nargout > 1
        ci = ptvec2mat(ci);
    end
    if nargout > 2
        stats.df = ptvec2mat(df);
        stats.sd = ptvec2mat(sd);
    end
end