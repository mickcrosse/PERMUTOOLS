function [b,p,ci,stats,dist] = permuglm(x,y,varargin)
%PERMUGLM  Permutation-based generalized linear models (GLM).
%   B = PERMUGLM(X,Y) performs a permutation-based generalized linear
%   regression of the responses in Y on the predictors in X, and returns
%   the estimated coefficients. X is an N-by-P design matrix, Y is an
%   N-by-V matrix of dependent variables, and B is a P-by-V matrix of
%   coefficients.
%
%   PERMUGLM leverages the Freedman-Lane algorithm using a one-step
%   Newton-Raphson approximation to isolate the unique variance of each
%   predictor, properly controlling for collinearity among covariates.
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
%                           'normal'    identity link (default)
%                           'binomial'  logit link (logistic regression)
%                           'poisson'   log link (count regression)
%       'type'          A string specifying the type of permutation to use:
%                           'freedmanlane'  permutes working residuals (default)
%                           'manly'         permutes raw data unrestricted
%                           'rank'          forces normal distribution & ranks Y
%       'alpha'         A scalar between 0 and 1 specifying the significance
%                       level as 100*ALPHA% (default=0.05).
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
    case 'freedmanlane'
    case 'manly'
    case 'rank'
        for k = 1:size(y,2)
            y(:,k) = tiedrank(y(:,k));
        end
    otherwise
        error('The TYPE parameter value must be FREEDMANLANE, or RANK.')
end

% Append intercept if requested and missing
if arg.intercept && ~any(all(x==1,1))
    x = [ones(size(x,1),1),x];
end

% Get data dimensions
[nobs,pnum] = size(x);
v = size(y,2);
df = nobs-pnum;

% Compute observed statistics
b = zeros(pnum,v);
se = zeros(pnum,v);
mu_full = zeros(nobs,v);
for k = 1:v
    [b(:,k),W_full,eta,~] = ptlocalfitglm(x,y(:,k),arg.distribution);
    XW = x.*W_full;
    inv_fisher = pinv(XW'*x);
    se(:,k) = sqrt(diag(inv_fisher));
    if strcmpi(arg.distribution,'binomial')
        mu_full(:,k) = 1./(1+exp(-eta));
    elseif strcmpi(arg.distribution,'poisson')
        mu_full(:,k) = exp(eta);
    else
        mu_full(:,k) = eta;
    end
end
tstat = b./se;

if nargout > 1

    % Initialize random number generator
    rng(arg.seed);

    % Preallocate output arrays
    p = zeros(pnum,v);
    if nargout > 2
        ci = zeros(2,v,pnum);
    end
    if nargout > 4
        dist = zeros(arg.nperm,v,pnum);
    end

    % Iterate across predictors
    for j = 1:pnum
        % Define permutation target
        xred = x;
        xred(:,j) = [];

        Mj = zeros(v,nobs);
        sestar = zeros(v,1);
        ered_all = zeros(nobs,v);
        etared_all = zeros(nobs,v);

        for k = 1:v
            switch arg.type
                case {'freedmanlane', 'rank'}
                    [~,Wred,etared_all(:,k),zred] = ptlocalfitglm(xred,...
                        y(:,k),arg.distribution);
                    ered_all(:,k) = zred-etared_all(:,k);
                case 'manly'
                    [~,Wred,~,zred] = ptlocalfitglm(ones(nobs,1),y(:,k),...
                        arg.distribution);
                    etared_all(:,k) = 0;
                    ered_all(:,k) = zred;
            end

            % Compute exact WLS projection matrix for fast permutation
            XW = x.*Wred;
            inv_fisher = pinv(XW'*x);
            Mfull = inv_fisher*XW';
            Mj(k,:) = Mfull(j,:);
            sestar(k) = sqrt(inv_fisher(j,j));
        end

        % Estimate sampling distribution
        distj = zeros(arg.nperm,v);
        for i = 1:arg.nperm
            randidx = randperm(nobs);
            eperm = ered_all(randidx,:);
            zperm = etared_all+eperm;
            bstar_j = sum(Mj'.*zperm,1);
            distj(i,:) = bstar_j./sestar';
        end

        % Apply max-correction if specified
        if arg.correct
            switch arg.tail
                case {'both','two'}
                    distj = max(abs(distj),[],2);
                case 'right'
                    distj = max(distj,[],2);
                case 'left'
                    distj = min(distj,[],2);
            end
        else
            switch arg.tail
                case {'both','two'}
                    distj = abs(distj);
            end
        end

        % Compute p-values
        switch arg.tail
            case {'both','two'}
                p(j,:) = (sum(abs(tstat(j,:))<=distj,1)+1)./(arg.nperm+1);
            case 'right'
                p(j,:) = (sum(tstat(j,:)<=distj,1)+1)./(arg.nperm+1);
            case 'left'
                p(j,:) = (sum(tstat(j,:)>=distj,1)+1)./(arg.nperm+1);
        end

        % Compute confidence intervals
        if nargout > 2
            switch arg.tail
                case {'both','two'}
                    crit=prctile(distj,100*(1-arg.alpha),1);
                    ci(:,:,j) = [b(j,:)-crit.*se(j,:);b(j,:)+crit.*se(j,:)];
                case 'right'
                    crit=prctile(distj,100*(1-arg.alpha),1);
                    ci(:,:,j) = [b(j,:)-crit.*se(j,:);Inf(1,v)];
                case 'left'
                    crit=prctile(distj,100*arg.alpha,1);
                    ci(:,:,j) = [-Inf(1,v);b(j,:)-crit.*se(j,:)];
            end
        end

        % Store null distribution
        if nargout > 4
            dist(:,:,j) = distj;
        end

    end

end

% Store statistics in a structure
if nargout > 3
    r2 = zeros(1,v);
    for k = 1:v
        % Calculate squared correlation as pseudo-R2
        c = corrcoef(y(:,k),mu_full(:,k));
        r2(k) = c(1,2)^2;
    end
    stats.method = arg.type;
    stats.distribution = arg.distribution;
    stats.df = df;
    stats.se = se;
    stats.tstat = tstat;
    stats.r2 = r2;
end

end

function [b, W, eta, z] = ptlocalfitglm(X, y, dist)
% Local function to perform iteratively reweighted least squares (IRLS)
[N, P] = size(X);
if P == 0
    b = [];
    eta = zeros(N,1);
    if strcmpi(dist, 'binomial')
        mu = 0.5 * ones(N,1);
        W = 0.25 * ones(N,1);
    elseif strcmpi(dist, 'poisson')
        mu = mean(y) * ones(N,1);
        W = mu; W(W < 1e-5) = 1e-5;
    else
        mu = zeros(N,1);
        W = ones(N,1);
    end
    if strcmpi(dist, 'normal')
        z = y;
    else
        z = eta + (y - mu) ./ W;
    end
    return;
end
if strcmpi(dist, 'normal')
    b = pinv(X) * y;
    eta = X * b;
    W = ones(N, 1);
    z = y;
    return;
end

% IRLS Initialization
if strcmpi(dist, 'binomial')
    mu = (y + 0.5) ./ 2;
    eta = log(mu ./ (1 - mu));
else
    mu = y + 0.1;
    eta = log(mu);
end
b = zeros(P, 1);
for iter = 1:50
    if strcmpi(dist, 'binomial')
        W = mu .* (1 - mu);
    else
        W = mu;
    end
    W(W < 1e-5) = 1e-5; % regularize
    z = eta + (y - mu) ./ W;
    XW = X .* W;
    b_new = pinv(XW' * X) * (XW' * z);

    if max(abs(b_new - b)) < 1e-6
        b = b_new;
        break;
    end
    b = b_new;
    eta = X * b;
    if strcmpi(dist, 'binomial')
        mu = 1 ./ (1 + exp(-eta));
    else
        mu = exp(eta);
    end
end
end