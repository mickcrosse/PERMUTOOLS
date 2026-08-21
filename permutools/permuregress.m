function [b,p,ci,stats,dist] = permuregress(x,y,varargin)
%PERMUREGRESS  Permutation-based multiple linear regression.
%   B = PERMUREGRESS(X,Y) performs a permutation-based multiple linear
%   regression of the responses in Y on the predictors in X, and returns
%   the estimated coefficients. X is an N-by-P design matrix, Y is an
%   N-by-V matrix of dependent variables, and B is a P-by-V matrix of
%   coefficients.
%
%   PERMUREGRESS leverages the Freedman-Lane algorithm by default to
%   isolate the unique variance of each predictor, properly controlling
%   for collinearity among covariates during permutation. Alternatively,
%   the Manly method may be selected for simple linear regression or when
%   evaluating overall model fit by setting the 'type' parameter to
%   'manly'.
%
%   For severe outliers or heavy-tailed data, raw observations may be
%   transformed to rank orders by setting the 'type' parameter to 'rank'.
%
%   [B,P] = PERMUREGRESS(...) returns the probabilities (i.e. p-values) of
%   observing the given result by chance if the null hypothesis is true. As
%   the null distribution is generated empirically by permuting the data,
%   no assumption is made about the shape of the distribution that the data
%   come from. P-values are automatically adjusted for multiple comparisons
%   using the max correction method.
%
%   [B,P,CI] = PERMUREGRESS(...) returns the 100*(1-ALPHA)% confidence
%   interval (CI) for each coefficient. CIs are also adjusted for multiple
%   comparisons using the max correction method.
%
%   [B,P,CI,STATS] = PERMUREGRESS(...) returns a structure with the
%   following fields:
%       'df'        -- the degrees of freedom (error)
%       'se'        -- the standard errors of the coefficients
%       'tstat'     -- the observed t-statistics
%       'res'       -- the residual values
%       'r2'        -- the R-square statistic (coefficient of determination)
%       'fstat'     -- the F-statistic for the overall model
%       'mse'       -- the estimated error variance
%       'method'    -- the statistical method used
%
%   [B,P,CI,STATS,DIST] = PERMUREGRESS(...) returns the permuted sampling
%   distributions of the t-statistics. DIST is an NPERM-by-V-by-P array.
%
%   [...] = PERMUREGRESS(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation to use:
%                       'freedmanlane'  permutes reduced model residuals
%                                       to control for covariates (default)
%                       'manly'         permutes raw response data for
%                                       simple regression or overall fit
%                       'rank'          permutes rank-transformed data for
%                                       severe outliers or heavy tails
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      coefficient is not 0 (default)
%                       'right'     coefficient is greater than 0
%                       'left'      coefficient is less than 0
%       'intercept' A numeric scalar (0,1) or logical specifying whether to
%                   automatically include a constant intercept (default=1).
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether to
%                   control FWER using max correction (default=1).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator (default=shuffle).
%
%   See also REGRESS FITLM PERMUCORR PERMUANOVA1.
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

% Handle missing values
nanidx = any(isnan(x),2) | any(isnan(y),2);
if any(nanidx(:))
    x(nanidx,:) = [];
    y(nanidx,:) = [];
end

% Transform raw data to rank orders if specified
switch arg.type
    case {'freedmanlane','manly'}
        % use raw data
    case 'rank'
        y = tiedrank(y);
        for j = 1:size(x,2)
            if length(unique(x(:,j))) > 2
                x(:,j) = tiedrank(x(:,j));
            end
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

% Compute projection matrices
invxtx = pinv(x'*x);
H = invxtx*x';
D = diag(invxtx);

% Compute observed statistics
df = nobs-npred;
b = H*y;
yhat = x*b;
res = y-yhat;
mse = sum(res.^2,1)./df;
se = sqrt(D.*mse);
tstat = b./se;

if nargout > 1

    rng(arg.seed);

    % Generate random permutations
    [~,idx] = sort(rand(nobs,arg.nperm),1);

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
        switch arg.type
            case {'freedmanlane','rank'}
                xred = x;
                xred(:,j) = [];
                if isempty(xred)
                    resred = y;
                else
                    betared = xred\y;
                    resred = y-xred*betared;
                end
            case 'manly'
                resred = y;
        end

        % Estimate sampling distribution
        distj = zeros(arg.nperm,nvar);
        for i = 1:arg.nperm
            resperm = resred(idx(:,i),:);
            bp = H*resperm;
            bstarj = bp(j,:);
            resstar = resperm-x*bp;
            msestar = sum(resstar.^2,1)./df;
            sestar = sqrt(D(j).*msestar);
            distj(i,:) = bstarj./sestar;
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

        % Store sampling distribution
        if nargout > 4
            if arg.correct
                dist(:,j) = refdist;
            else
                dist(:,:,j) = distj;
            end
        end

    end

end

% Compute & format descriptive statistics
if nargout > 3
    dfmodel = npred-any(all(x==1,1));
    sst = sum((y-mean(y,1)).^2,1);
    sse = sum(res.^2,1);
    r2 = 1-(sse./sst);
    fstat = (r2./dfmodel)./((1-r2)./df);
    stats.df = df;
    stats.se = se;
    stats.tstat = tstat;
    stats.res = res;
    stats.r2 = r2;
    stats.fstat = fstat;
    stats.mse = mse;
    stats.method = arg.type;
end