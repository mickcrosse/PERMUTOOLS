function [b,p,ci,stats,dist] = permuglm(x,y,varargin)
%PERMUGLM  Permutation-based generalized linear models (GLM).
%   B = PERMUGLM(X,Y) performs a permutation-based generalized linear
%   regression of the responses in Y on the predictors in X, and returns
%   the estimated coefficients. X is an N-by-P design matrix, Y is an
%   N-by-V matrix of dependent variables, and B is a P-by-V matrix of
%   coefficients.
%
%   PERMUGLM leverages the Freedman-Lane algorithm by default using a one-
%   step Newton-Raphson approximation to isolate the unique variance of
%   each predictor, properly controlling for collinearity among covariates.
%   Alternatively, the Manly method may be selected for simple GLMs,
%   evaluating overall model fit, or fitting Binomial or Poisson response
%   distributions by setting the 'type' parameter to 'manly'.
%
%   For severe outliers or heavy-tailed data, raw observations may be
%   transformed to rank orders by setting the 'type' parameter to 'rank'.
%
%   [B,P] = PERMUGLM(...) returns the probabilities (i.e. p-values) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by permuting the data,
%   no assumption is made about the shape of the distribution that the data
%   come from. P-values are automatically adjusted for multiple comparisons
%   using the max correction method.
%
%   [B,P,CI] = PERMUGLM(...) returns the 100*(1-ALPHA)% confidence
%   intervals (CI) for the regression coefficients. CI is a 2-by-V-by-P
%   array. CIs are adjusted for multiple comparisons using max correction.
%
%   [B,P,CI,STATS] = PERMUGLM(...) returns a structure with the
%   following fields:
%       'method'       -- the statistical method used
%       'distribution' -- the GLM distribution family
%       'df'           -- the degrees of freedom (error)
%       'se'           -- the standard errors of the coefficients
%       'tstat'        -- the observed pseudo t-statistics
%       'r2'           -- the pseudo R-squared (squared correlation)
%
%   [B,P,CI,STATS,DIST] = PERMUGLM(...) returns the permuted sampling
%   distributions of the t-statistics. DIST is an NPERM-by-V-by-P array.
%
%   [...] = PERMUGLM(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter       Value
%       'distribution'  A string specifying the distribution family:
%                           'normal'        identity link for continuous
%                                           data (default)
%                           'binomial'      logit link for proportion or
%                                           binary data (logistic regress.)
%                           'poisson'       log link for non-negative count
%                                           data (Poisson regression)
%       'type'          A string specifying the type of permutation to use:
%                           'freedmanlane'  permutes reduced residuals for
%                                           GLMs with covariates (default)
%                           'manly'         permutes raw response data for
%                                           simple or Binomial/Poisson GLMs
%                           'rank'          permutes rank-transformed data
%                                           for outliers or heavy tails
%       'alpha'         A scalar between 0 and 1 specifying the
%                       significance level as 100*ALPHA% (default=0.05).
%       'tail'          A string specifying the alternative hypothesis:
%                           'both'      coefficient is not 0 (default)
%                           'right'     coefficient is greater than 0
%                           'left'      coefficient is less than 0
%       'intercept'     A numeric scalar (0,1) or logical indicating
%                       whether to automatically include a constant
%                       intercept (default=1).
%       'nperm'         An integer scalar specifying the number of
%                       permutations (default=10,000).
%       'correct'       A numeric scalar (0,1) or logical indicating
%                       whether to control FWER using max correction
%                       (default=1).
%       'seed'          An integer scalar specifying the seed value used to
%                       initialise the permutation generator
%                       (default=shuffle).
%
%   See also FITGLM PERMUREGRESS PERMUCORR PERMUANOVA1.
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
    arg.type = 'freedmanlane';
end
if isempty(arg.distribution)
    arg.distribution = 'normal';
end

% Handle missing values
nanidx = any(isnan(x),2) | any(isnan(y),2);
if any(nanidx(:))
    x(nanidx,:) = [];
    y(nanidx,:) = [];
end

% Transform raw data to rank orders if specified
switch arg.type
    case {'freedmanlane', 'manly'}
    case 'rank'
        for k = 1:size(y,2)
            y(:,k) = tiedrank(y(:,k));
        end
    otherwise
        error('The TYPE parameter value must be FREEDMANLANE, MANLY, or RANK.')
end

% Append intercept if requested and missing
if arg.intercept && ~any(all(x==1,1))
    x = [ones(size(x,1),1),x];
end

% Get data dimensions
[nobs,npred] = size(x);
nvar = size(y,2);

% Compute observed statistics
b = zeros(npred,nvar);
se = zeros(npred,nvar);
mu = zeros(nobs,nvar);
dispersion = zeros(1,nvar);
df = nobs-npred;
for k = 1:nvar
    [b(:,k),Wfull,eta] = ptlocalfitglm(x,y(:,k),arg.distribution);
    XW = x.*Wfull;
    invfisher = pinv(XW'*x);
    switch arg.distribution
        case 'normal'
            res = y(:,k)-eta;
            dispersion(k) = sum(res.^2)/df;
        otherwise
            dispersion(k) = 1;
    end
    se(:,k) = sqrt(diag(invfisher).*dispersion(k));
    switch arg.distribution
        case 'normal'
            mu(:,k) = eta;
        case 'binomial'
            mu(:,k) = 1./(1+exp(-eta));
        case 'poisson'
            mu(:,k) = exp(eta);
    end
end
tstat = b./se;

if nargout > 1

    rng(arg.seed);

    % Preallocate output arrays
    p = zeros(npred,nvar);
    if nargout > 2
        ci = zeros(2,nvar,npred);
    end
    if nargout > 4
        if arg.correct
            dist = zeros(arg.nperm,npred);
        else
            dist = zeros(arg.nperm,nvar,npred);
        end
    end

    % Iterate across predictors
    for j = 1:npred

        % Define permutation target
        xred = x;
        xred(:,j) = [];
        Mj = zeros(nvar,nobs);
        sestar = zeros(nvar,1);
        ered = zeros(nobs,nvar);
        etared = zeros(nobs,nvar);
        for k = 1:nvar
            switch arg.type
                case {'freedmanlane','rank'}
                    [~,Wred,etared(:,k),zred] = ptlocalfitglm(xred,...
                        y(:,k),arg.distribution);
                    ered(:,k) = zred-etared(:,k);
                case 'manly'
                    [~,Wred,~,zred] = ptlocalfitglm(ones(nobs,1),...
                        y(:,k),arg.distribution);
                    ered(:,k) = zred;
            end

            % Compute exact WLS projection matrix for fast permutation
            XW = x.*Wred;
            invfisher = pinv(XW'*x);
            Mfull = invfisher*XW';
            Mj(k,:) = Mfull(j,:);
            sestar(k) = sqrt(invfisher(j,j).*dispersion(k));

        end

        % Estimate sampling distribution
        distj = zeros(arg.nperm,nvar);
        for i = 1:arg.nperm
            randidx = randperm(nobs);
            eperm = ered(randidx,:);
            zperm = etared+eperm;
            bstarj = sum(Mj'.*zperm,1);
            distj(i,:) = bstarj./sestar';
        end

        % Apply max-correction if specified
        if arg.correct
            switch arg.tail
                case 'both'
                    refdist = max(abs(distj),[],2);
                case 'right'
                    refdist = max(distj,[],2);
                case 'left'
                    refdist = min(distj,[],2);
            end
        else
            refdist = distj;
        end

        % Compute p-values
        switch arg.tail
            case 'both'
                p(j,:) = (sum(abs(tstat(j,:))<=abs(refdist),1)+1)./...
                    (arg.nperm+1);
            case 'right'
                p(j,:) = (sum(tstat(j,:)<=refdist,1)+1)./(arg.nperm+1);
            case 'left'
                p(j,:) = (sum(tstat(j,:)>=refdist,1)+1)./(arg.nperm+1);
        end

        % Compute confidence intervals
        if nargout > 2
            switch arg.tail
                case 'both'
                    crit = prctile(abs(refdist),100*(1-arg.alpha),1);
                    ci(:,:,j) = [b(j,:)-crit.*se(j,:);b(j,:)+crit.*se(j,:)];
                case 'right'
                    crit = prctile(refdist,100*(1-arg.alpha),1);
                    ci(:,:,j) = [b(j,:)-crit.*se(j,:);Inf(1,nvar)];
                case 'left'
                    crit = prctile(refdist,100*arg.alpha,1);
                    ci(:,:,j) = [-Inf(1,nvar);b(j,:)-crit.*se(j,:)];
            end
        end

        % Store null distribution
        if nargout > 4
            if arg.correct
                dist(:,j) = refdist;
            else
                dist(:,:,j) = distj;
            end
        end

    end

end

% Store statistics in a structure
if nargout > 3
    r = permucorr(y,mu);
    stats.distribution = arg.distribution;
    stats.df = df;
    stats.r2 = r.^2;
    stats.tstat = tstat;
    stats.res = y-mu;
    stats.se = se;
    stats.method = arg.type;
end

end

function [b,W,eta,z] = ptlocalfitglm(X,y,dist)
% Local function to perform iteratively reweighted least squares (IRLS)
[N,P] = size(X);
if P == 0
    b = [];
    eta = zeros(N,1);
    if strcmpi(dist,'binomial')
        mu = 0.5*ones(N,1);
        W = 0.25*ones(N,1);
    elseif strcmpi(dist,'poisson')
        mu = mean(y)*ones(N,1);
        W = mu; W(W<1e-5) = 1e-5;
    else
        mu = zeros(N,1);
        W = ones(N,1);
    end
    if strcmpi(dist,'normal')
        z = y;
    else
        z = eta+(y-mu)./W;
    end
    return;
end
if strcmpi(dist,'normal')
    b = pinv(X)*y;
    eta = X*b;
    W = ones(N,1);
    z = y;
    return;
end

% IRLS Initialization
if strcmpi(dist,'binomial')
    mu = (y+0.5)./2;
    eta = log(mu./(1-mu));
else
    mu = y+0.1;
    eta = log(mu);
end
b = zeros(P,1);
for iter = 1:50
    if strcmpi(dist,'binomial')
        W = mu.*(1-mu);
    else
        W = mu;
    end
    W(W < 1e-5) = 1e-5;
    z = eta+(y-mu)./W;
    XW = X.*W;
    bnew = pinv(XW'*X)*(XW'*z);
    if max(abs(bnew-b)) < 1e-6
        b = bnew;
        break;
    end
    b = bnew;
    eta = X*b;
    if strcmpi(dist,'binomial')
        mu = 1./(1+exp(-eta));
    else
        mu = exp(eta);
    end
end
end