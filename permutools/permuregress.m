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
%   for collinearity among covariates during permutation.
% 
%   [B,P] = PERMUREGRESS(...) returns the probabilities (i.e. p-values) of 
%   observing the given result by chance if the null hypothesis is true. As 
%   the null distribution is generated empirically by permuting the data, 
%   no assumption is made about the shape of the distribution that the data 
%   come from. P-values are automatically adjusted for multiple comparisons 
%   using the max correction method.
%
%   [B,P,CI] = PERMUREGRESS(...) returns the 100*(1-ALPHA)% confidence
%   intervals (CI) for the regression coefficients. CI is a 2-by-V-by-P
%   array. CIs are adjusted for multiple comparisons using max-correction.
%
%   [B,P,CI,STATS] = PERMUREGRESS(...) returns a structure with the
%   following fields:
%       'method'    -- the statistical method used
%       'df'        -- the degrees of freedom (error)
%       'se'        -- the standard errors of the coefficients
%       'tstat'     -- the observed t-statistics
%       'r2'        -- the R-squared (coefficient of determination)
%       'fstat'     -- the F-statistic for the overall model
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
%                       'freedmanlane'  permutes reduced residuals (default)
%                       'manly'         permutes raw data unrestricted
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      two-tailed test (default)
%                       'right'     right-tailed test
%                       'left'      left-tailed test
%       'intercept' A numeric scalar (0,1) or logical specifying whether to
%                   automatically include a constant intercept (default=1).
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'correct'   A numeric scalar (0,1) or logical indicating whether
%                   to control FWER using max correction (default=1).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator.
%
%   See also FITLM REGRESS PERMUCORR PERMUANOVA1.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

% Parse input arguments
arg = ptparsevarargin(varargin);
if ~any(strcmpi(varargin,'type'))
    arg.type = 'freedmanlane';
end
if ~any(strcmpi(varargin,'intercept'))
    arg.intercept = 1;
end

% Handle missing values (Listwise deletion)
nanidx = any(isnan(x),2) | any(isnan(y),2);
if any(nanidx(:))
    x(nanidx,:) = [];
    y(nanidx,:) = [];
end

% Append intercept if requested and missing
if arg.intercept && ~any(all(x==1,1))
    x = [ones(size(x,1),1),x];
end

% Get data dimensions
[nobs,pnum] = size(x);
v = size(y,2);
df = nobs-pnum;

% Compute projection matrices
H = pinv(x'*x)*x';
D = diag(pinv(x'*x));

% Compute observed statistics
b = H*y;
yhat = x*b;
res = y-yhat;
mse = sum(res.^2,1)./df;
se = sqrt(D.*mse);
tstat = b./se;

if nargout > 1

    % Initialize random number generator
    rng(arg.seed);

    % Preallocate output arrays
    p = zeros(pnum,v);
    if nargout>2
        ci = zeros(2,v,pnum);
    end
    if nargout > 4
        dist = zeros(arg.nperm,v,pnum);
    end

    % Precompute full model identity projection
    IdH = eye(nobs)-x*H;

    % Iterate across predictors
    for j = 1:pnum

        % Define permutation target based on method
        switch arg.tail
            case 'freedmanlane'
                xred = x;
                xred(:,j) = [];
                if isempty(xred)
                    resred = y;
                else
                    Hred = pinv(xred'*xred)*xred';
                    betared = Hred*y;
                    resred = y-xred*betared;
            end
            case 'manly'
                resred = y;
        end

        Hj = H(j,:);
        distj = zeros(arg.nperm,v);

        % Permutation loop
        for i = 1:arg.nperm
            randidx = randperm(nobs);
            resperm = resred(randidx,:);
            bstar_j = Hj*resperm;
            resstar = IdH*resperm;
            msestar = sum(resstar.^2,1)./df;
            sestar = sqrt(D(j).*msestar);
            distj(i,:) = bstar_j./sestar;
        end

        % Apply max-correction if specified
        if arg.correct
            switch arg.tail
                case 'both'
                    distj = max(abs(distj),[],2);
                case 'right'
                    distj = max(distj,[],2);
                case 'left'
                    distj = min(distj,[],2);
            end
        else
            switch arg.tail
                case 'both'
                    distj = abs(distj);
            end
        end

        % Compute p-values
        switch arg.tail
            case 'both'
                p(j,:) = (sum(abs(tstat(j,:))<=distj,1)+1)./(arg.nperm+1);
            case 'right'
                p(j,:) = (sum(tstat(j,:)<=distj,1)+1)./(arg.nperm+1);
            case 'left'
                p(j,:) = (sum(tstat(j,:)>=distj,1)+1)./(arg.nperm+1);
        end

        % Compute confidence intervals
        if nargout > 2
            switch arg.tail
                case 'both'
                    crit=prctile(distj,100*(1-arg.alpha),1);
                    ci(1,:,j)=b(j,:)-crit.*se(j,:);
                    ci(2,:,j)=b(j,:)+crit.*se(j,:);
                case 'right'
                    crit=prctile(distj,100*(1-arg.alpha),1);
                    ci(1,:,j)=b(j,:)-crit.*se(j,:);
                    ci(2,:,j)=Inf;
                case 'left'
                    crit=prctile(distj,100*arg.alpha,1);
                    ci(1,:,j)=-Inf;
                    ci(2,:,j)=b(j,:)-crit.*se(j,:);
            end
        end

        % Store null distribution
        if nargout > 4
            dist(:,:,j) = distj;
        end

    end
end

% Compute overall model fit statistics
if nargout > 3
    sst = sum((y-mean(y,1)).^2,1);
    sse = sum(res.^2,1);
    r2 = 1-(sse./sst);
    fstat = (r2./(pnum-1))./((1-r2)./df);

    % Store statistics in a structure
    stats.method = arg.type;
    stats.df = df;
    stats.se = se;
    stats.tstat = tstat;
    stats.r2 = r2;
    stats.fstat = fstat;
end