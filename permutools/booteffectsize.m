function [d,ci,stats] = booteffectsize(x,varargin)
%BOOTEFFECTSIZE  Effect size with bootstrapped confidence intervals.
%   D = BOOTEFFECTSIZE(X) returns the effect size measure for a single
%   sample X based on Cohen's d. By default, Cohen's d is bias-corrected
%   for sample size, also known as Hedges' g.
%
%   If X is a matrix, separate effect sizes are measured along each column
%   of X, and a vector of results is returned. If the parameter TEST is set
%   to 'PAIRWISE', the effect sizes between every pair of columns in X are
%   measured, and a matrix of results is returned.
%
%   BOOTEFFECTSIZE treats NaNs as missing values, and ignores them.
%
%   D = BOOTEFFECTSIZE(X,Y) returns the effect size between two dependent
%   samples X and Y using the pooled standard deviation. X and Y can be
%   treated as independent samples by setting the PAIRED parameter to 0. If
%   X and Y are independent samples with significantly different variances,
%   an estimate based on the control sample's variance (Glass' Δ) can be
%   computed by setting the EFFECT parameter to 'GLASS'. For this measure,
%   the control sample should be entered as X, and the test sample as Y.
%
%   [D,CI] = BOOTEFFECTSIZE(...) returns the bootstrapped, bias-corrected
%   confidence intervals using the percentile method.
%
%   [D,CI,STATS] = BOOTEFFECTSIZE(...) returns a structure with the
%   following fields:
%       'df'        -- the degrees of freedom of each measure
%       'sd'    	-- the pooled standard deviation, or of X for Glass' Δ
%
%   [...] = BOOTEFFECTSIZE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'nboot'     A scalar specifying the number of bootstraps used to
%                   estimate the confidence intervals (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   bias-correct the effect size and confidence intervals
%                   according to sample size (default=true).
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'paired'    A numeric scalar (0,1) or logical indicating whether
%                   the data in X and Y are paired: pass in 1 for paired
%                   samples (default), or 0 for unpaired samples.
%       'vartype'   A string specifying the variance equivalence of X and Y
%                   to determine the SD and degrees of freedom:
%                       'equal'   	assume equal variances (default)
%                       'unequal' 	assume unequal variances
%       'effect'    A string specifying the effect size to compute:
%                       'Cohen'     compute Cohen's d (default)
%                       'Glass'     compute Glass' Δ when comparing
%                                   independent samples with significantly
%                                   different variances
%       'test'      A string specifying whether to compute a one-sample
%                   measure or a pairwise measure when only X is entered:
%                       'one'       compare each column of X to zero and
%                                   return a vector of results (default)
%                       'pairwise'  compare each pair of columns in X and
%                                   return a matrix of results
%       'seed'      A scalar integer specifying the seed used to initialise
%                   the bootstrap generator. By default, the generator is
%                   initialised based on the current time, resulting in a
%                   different bootstrap each time.
%
%   See also MEANEFFECTSIZE BOOTCI PERMUTTEST PERMUTTEST2 PERMUVARTEST2.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Hentschke H, Stuttgen MC (2011) Computation of measures of
%           effect size for neuroscience data sets. Eur J Neurosci,
%           34:1887–1894.
%       [2] Cohen J (1969) Statistical power for the behavioural sciences.
%           London: Academic Press.
%       [3] Hedges LV, Olkin I (1985) Statistical methods for meta-
%           analysis. San Diego, CA: Academic Press.

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

if nargin>=2 && ~ischar(varargin{1})
    y = varargin{1};
    varargin = varargin(2:end);
else
    y = [];
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

% Set up permutation test
if isempty(y)
    switch arg.test
        case 'one'
            y = zeros(size(x));
        case 'pairwise'
            warning('Comparing all columns of X...')
            [x,y] = ptpaircols(x);
            arg.mat = true;
    end
else
    switch arg.test
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
        case 'Glass'
            error('GLASS can only be used for independent samples.')
    end
    switch arg.vartype
        case 'unequal'
            error('Cannot assume unequal variance for dependent samples.')
    end

    % Compute difference between samples
    diffxy = x-y;

    % Use only rows with no NaN values if specified
    switch arg.rows
        case 'complete'
            diffxy = diffxy(~any(isnan(diffxy),2),:);
    end

    % Get data dimensions, ignoring NaNs
    nobs = sum(~isnan(diffxy));

    % Compute mean difference
    mu = sum(diffxy,nanflag)./nobs;

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

    % Get data dimensions, ignoring NaNs
    nobs = sum(~isnan([x;y]));

    % Compute mean difference
    mu = smx./nobsx-smy./nobsy;

    % Compute standard deviation
    switch arg.effect
        case 'Cohen'
            switch arg.vartype
                case 'equal'
                    df = nobs-2;
                    sd = sqrt((dfx.*varx+dfy.*vary)./df);
                case 'unequal'
                    df = (dfx.*dfy.*(varx+vary).^2)./(dfy.*varx.^2+...
                        dfx.*vary.^2);
                    sd = sqrt((varx+vary)/2);
            end
        case 'Glass'
            switch arg.vartype
                case 'equal'
                    warning(['GLASS option should only be used for '...
                        'samples with unequal variances.'])
            end
            df = dfx;
            sd = sqrt(varx);
    end

end

% Compute effect size
d = mu./sd;

if nargout > 1

    % Estimate sampling distribution
    rng(arg.seed);
    ci = zeros(2,nvar);
    for i = 1:nvar

        % Bootstrap samples
        xb = x(:,i);
        yb = y(:,i);
        idx = ceil(nobsx(i)*rand([nobsx(i),arg.nboot]));
        xb = xb(idx);
        if ~arg.paired
            idx = ceil(nobsy(i)*rand([nobsy(i),arg.nboot]));
        end
        yb = yb(idx);

        % Compute sample variance using fast algo
        smxb = sum(xb,nanflag);
        smyb = sum(yb,nanflag);
        varxb = (sum(xb.^2,nanflag)-(smxb.^2)/nobsx(i))/dfx(i);
        varyb = (sum(yb.^2,nanflag)-(smyb.^2)/nobsy(i))/dfy(i);

        if arg.paired

            % Compute mean difference
            mub = sum(xb-yb,nanflag)/nobs(i);

            % Compute standard deviation
            sdb = sqrt((varxb+varyb)/2);

        else

            % Compute mean difference
            mub = smxb/nobsx(i)-smyb/nobsy(i);

            % Compute standard deviation
            switch arg.effect
                case 'Cohen'
                    switch arg.vartype
                        case 'equal'
                            sdb = sqrt((dfx(i)*varxb+dfy(i)*varyb)/df(i));
                        case 'unequal'
                            sdb = sqrt((varxb+varyb)/2);
                    end
                case 'Glass'
                    sdb = sqrt(varxb);
            end

        end

        % Compute effect size
        db = mub./sdb;

        % Compute confidence intervals
        ci(:,i) = prctile(db,100*[arg.alpha/2;1-arg.alpha/2]);

    end

end

% Bias-correct results for sample size
if arg.correct
    factor = exp(gammaln(df/2)-log(sqrt(df/2))-gammaln((df-1)/2));
    d = factor.*d;
    if nargout > 1
        ci = [factor;factor].*ci;
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