function [d,ci,stats,dist] = booteffectsize(x,m,varargin)
%BOOTEFFECTSIZE  Effect size with bootstrapped confidence intervals.
%   D = BOOTEFFECTSIZE(X) returns the effect size measure for a single
%   sample X based on Cohen's d. By default, Cohen's d is bias-corrected
%   for sample size, also known as Hedges' g. For ordinal data, Cliff's
%   delta can be computed by setting the 'effect' parameter to 'cliff'.
%
%   If X is a matrix, separate effect sizes are measured along each column
%   of X, and a vector of results is returned. If the 'compare' parameter
%   is set to 'pairwise', the effect sizes between every pair of columns in
%   X are measured, and a matrix of results is returned.
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
%       'effect'    A string specifying the effect size to measure:
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
%       'nboot'     A scalar specifying the number of bootstraps used to
%                   estimate the confidence intervals (default=10,000).
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
%                   different bootstrap each time.
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

%   © 2018-2023 Mick Crosse <crossemj@tcd.ie>
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
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute sample variance using fast algo
smx = sum(x,nanflag);
smy = sum(y,nanflag);
varx = (sum(x.^2,nanflag)-(smx.^2)./nobsx)./dfx;
vary = (sum(y.^2,nanflag)-(smy.^2)./nobsy)./dfy;

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
            diffxy = zeros(max(nobsx)^2,nvar);
            for i = 1:nvar
                diffi = sign(x(:,i)-y(:,i)');
                diffxy(:,i) = diffi(:);
            end
            diffxy(1:max(nobsx)+1:end,:) = 0;
        otherwise
            diffxy = x-y;
    end

    % Use only rows with no NaN values if specified
    switch arg.rows
        case 'complete'
            diffxy = diffxy(~any(isnan(diffxy),2),:);
    end

    % Get data dimensions, ignoring NaNs
    switch arg.effect
        case 'cliff'
            nobs = nobsx.*dfx;
        otherwise
            nobs = sum(~isnan(diffxy));
    end

    % Compute mean difference
    switch arg.effect
        case 'mediandiff'
            mu = median(diffxy,nanflag);
        otherwise
            mu = sum(diffxy,nanflag)./nobs;
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
            diffxy = zeros(max(nobsx)*max(nobsy),nvar);
            for i = 1:nvar
                diffi = sign(x(:,i)-y(:,i)');
                diffxy(:,i) = diffi(:);
            end
        otherwise
            diffxy = [x;y];
    end

    % Get data dimensions, ignoring NaNs
    nobs = sum(~isnan(diffxy));

    % Compute mean difference
    switch arg.effect
        case 'cliff'
            mu = sum(diffxy,nanflag)./nobs;
        case 'mediandiff'
            mu = median(x,nanflag)-median(y,nanflag);
        otherwise
            mu = smx./nobsx-smy./nobsy;
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

    for i = 1:nvar

        % Generate random bootstraps
        xb = x(:,i);
        yb = y(:,i);
        idx = ceil(nobsx(i)*rand([nobsx(i),arg.nboot]));
        xb = xb(idx);
        if ~arg.paired
            idx = ceil(nobsy(i)*rand([nobsy(i),arg.nboot]));
        end
        yb = yb(idx);

        % Estimate sampling distribution
        smxb = sum(xb,nanflag);
        smyb = sum(yb,nanflag);
        varxb = (sum(xb.^2,nanflag)-(smxb.^2)/nobsx(i))/dfx(i);
        varyb = (sum(yb.^2,nanflag)-(smyb.^2)/nobsy(i))/dfy(i);
        if arg.paired
            switch arg.effect
                case 'cliff'
                    diffxyb = zeros(max(nobsx)^2,arg.nboot);
                    for j = 1:arg.nboot
                        diffi = sign(xb(:,j)-yb(:,j)');
                        diffxyb(:,j) = diffi(:);
                    end
                    diffxyb(1:max(nobsx)+1:end,:) = 0;
                otherwise
                    diffxyb = xb-yb;
            end
            switch arg.effect
                case 'mediandiff'
                    mub = median(diffxyb,nanflag);
                otherwise
                    mub = sum(diffxyb,nanflag)/nobs(i);
            end
            sdb = sqrt((varxb+varyb)/2);
        else
            switch arg.effect
                case 'cliff'
                    diffxyb = zeros(max(nobsx)*max(nobsy),arg.nboot);
                    for j = 1:arg.nboot
                        diffi = sign(xb(:,j)-yb(:,j)');
                        diffxyb(:,j) = diffi(:);
                    end
            end
            switch arg.effect
                case 'cliff'
                    mub = sum(diffxyb,nanflag)/nobs(i);
                case 'mediandiff'
                    mub = median(xb,nanflag)-median(yb,nanflag);
                otherwise
                    mub = smxb/nobsx(i)-smyb/nobsy(i);
            end
            switch arg.effect
                case 'glass'
                    sdb = sqrt(varxb);
                case 'cohen'
                    switch arg.vartype
                        case 'equal'
                            sdb = sqrt((dfx(i)*varxb+dfy(i)*varyb)/df(i));
                        case 'unequal'
                            sdb = sqrt((varxb+varyb)/2);
                    end
            end
        end
        switch arg.effect
            case {'cliff','meandiff','mediandiff'}
                dist = mub;
            otherwise
                dist = mub./sdb;
        end

        % Compute CI
        ci(:,i) = prctile(dist,100*[arg.alpha/2;1-arg.alpha/2]);

    end

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
    end
end

% Arrange results in a matrix if specified
if arg.mat
    d = ptvec2mat(d);
    if nargout > 1
        ci = ptvec2mat(ci);
    end
    if nargout > 2
        df = ptvec2mat(df);
        sd = ptvec2mat(sd);
    end
end

% Store statistics in a structure
if nargout > 2
    stats = struct('df',df,'sd',sd);
end