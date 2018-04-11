function [rval,corx,orig,stats] = st_rmaxperm(x,y,nperm,tail,alpha,type)
%st_rmaxperm StatsTools correlation permutation test with Rmax correction
%   RVAL = ST_RMAXPERM(X,Y,NPERM,TAIL,ALPHA,TYPE) returns the pairwise
%   correlation coefficient of a nonparametric permutation test based on
%   the Pearson correlation coefficient or Spearman's rank correlation
%   coefficient. If X and Y are matrices, multiple permutation tests are
%   performed simultaneously between each pair of columns in X and Y and
%   family-wise error rate is controlled using the max statistic correction
%   method (Blair et al., 1994; Westfall & Young, 1993; Groppe et al., 
%   2011). 
%
%   [...,CORX] = ST_RMAXPERM(...) returns a structure containing the
%   corrected test statistics of the permutation tests.
%
%   [...,ORIG] = ST_RMAXPERM(...) returns a structure containing the
%   original, uncorrected test statistics of the permutation tests.
%
%   [...,STATS] = ST_RMAXPERM(...) returns a structure containing some
%   general data statistics including the effect size.
%
%   Inputs:
%   x     - column vector or matrix of data (observations by variables)
%   y     - column vector or matrix of data (observations by variables)
%   nperm - scalar specifying the number of permutations (default=10,000 or 
%           all possible permutations are computed for less than 8 obs.)
%   tail  - string specifying the alternative hypothesis
%           'both'  - correlation is not zero (default)
%           'right' - correlation is greater than zero
%           'left'  - correlation is less than zero
%   alpha - scalar between 0 and 1 specifying the significance level as 
%           100*ALPHA% (default=0.05)
%   type  - string specifying the type of correlation measure
%           'Pearson'  - Pearson correlation coefficient (default)
%           'Spearman' - Spearman's rank correlation coefficient
%
%   Outputs:
%   rval  - scalar or vector containing the correlation coefficients
%   corx  - structure of corrected test statistics containing the following
%           fields:
%           p     - probability of observing the given result by chance
%           rcrit - lower and upper critical r-values for given alpha level
%           estal - estimated alpha level of each test
%   orig  - structure of original, uncorrected test statistics containing
%           the same fields as CORX
%   stats - structure of data statistics containing the following fields:
%           df    - degrees of freedom of each test
%           sdx   - estimated population standard deviation of X
%           sdy   - estimated population standard deviation of Y
%
%   See README for examples of use.
%
%   See also ST_TMAXPERM ST_TMAXPERM2 ST_BOXDOTPLOT ST_PLOTPOLYCI.
%
%   StatsTools https://github.com/mickcrosse/StatsTools

%   References:
%      [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%          Multivariate Permutation Tests Which May Replace Hotelling's T2
%          Test in Prescribed Circumstances. Mult Behav Res, 29(2):141-163.
%      [2] Westfall PH, Young SS (1993) Resampling-Based Multiple Testing:
%           Examples and Methods for p-Value Adjustment. Wiley, New York.
%      [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis of
%          event-related brain potentials/fields I: A critical tutorial
%          review. Psychophysiology, 48(12):1711-1725.
%      [4] Cohen J (1988) Statistical Power Analysis for the Behavioral 
%          Sciences, 2nd Ed. Lawrence Earlbaum Associates, Hilsdale, NJ.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 20-Mar-2018

if ~exist('y','var') || isempty(y)
    error('st_rmaxperm: Requires at least two input arguments')
end
if ~exist('nperm','var') || isempty(nperm)
    nperm = 1e4;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
if ~exist('tail','var') || isempty(tail)
    tail = 'both';
end
if ~exist('type','var') || isempty(type)
    type = 'Pearson'; %#ok<NASGU>
elseif strcmpi(type,'Spearman')
    x = tiedrank(x);
    y = tiedrank(y);
end

% Compute some constants
[nobs,nvar] = size(x);
muxy = sum(x).*sum(y)/nobs;
sdxy = sqrt((sum(y.^2)-(sum(y).^2)/nobs).*(sum(x.^2)-(sum(x).^2)/nobs));
if size(y,1)~=nobs
    error('X and Y must have the same number of observations')
elseif size(y,2)~=nvar
    error('X and Y must have the same number of variables')
end

% Compute correlation coefficient
rval = (sum(x.*y)-muxy)./sdxy;

% Permute data and generate distribution of r-values
if nobs < 8
    nperm = factorial(nobs);
    idx = perms(1:nobs)';
    fprintf('Calculating all possible permutations because of small N.')
else
    [~,idx] = sort(rand(nobs,nperm));
end
rp = zeros(nperm,nvar);
for i = 1:nperm
    rp(i,:) = (sum(x(idx(:,i),:).*y)-muxy)./sdxy;
end

% Compute adjusted test statistics using Tmax correction
if strcmpi(tail,'both')
    [~,idx] = max(abs(rp),[],2);
    csvar = [0;cumsum(ones(nperm-1,1)*nvar)];
    rpT = rp'; rmax = rpT(idx+csvar);
    p = zeros(1,nvar);
    p(rval>0) = mean(abs(rval)<rmax)*2;
    p(rval<=0) = mean(abs(rval)>rmax)*2;
    rcrit(1) = prctile(rmax,100*alpha/2);
    rcrit(2) = prctile(rmax,100-100*alpha/2);
    estal = mean(rcrit(2)<rmax)+mean(rcrit(1)>rmax);
elseif strcmpi(tail,'right')
    rmax = max(rp,[],2);
    p = mean(rval<rmax);
    rcrit = prctile(rmax,100-100*alpha);
    estal = mean(rcrit<rmax);
elseif strcmpi(tail,'left')
    rmax = min(rp,[],2);
    p = mean(rval>rmax);
    rcrit = prctile(rmax,100*alpha);
    estal = mean(rcrit>rmax);
end
p(isnan(rcrit)) = NaN;

% Store values in structure
corx = struct('p',p,'rcrit',rcrit,'estal',estal);

% Execute if user specifies uncorrected test
if nargout > 2
    
    % Clear variables
    clear p rcrit estal
    
    % Compute original test statistics without correction
    if strcmpi(tail,'both')
        p = zeros(1,nvar);
        p(rval>0) = mean(abs(rval)<rp)*2;
        p(rval<=0) = mean(abs(rval)>rp)*2;
        rcrit(1,:) = prctile(rp,100*alpha/2);
        rcrit(2,:) = prctile(rp,100-100*alpha/2);
        estal = mean(rcrit(2,:)<rp)+mean(rcrit(1,:)>rp);
    elseif strcmpi(tail,'right')
        p = mean(rval<rp);
        rcrit = prctile(rp,100-100*alpha);
        estal = mean(rcrit<rp);
    elseif strcmpi(tail,'left')
        p = mean(rval>rp);
        rcrit = prctile(rp,100*alpha);
        estal = mean(rcrit>rp);
    end
    p(isnan(rcrit(1,:))) = NaN;
    
    % Store values in structure
    orig = struct('p',p,'rcrit',rcrit,'estal',estal);
    
end

% Execute if user specifies general data statistics
if nargout > 3
    df = nobs-2;
    if nvar>1
        df = repmat(df,1,nvar);
    end
    stats = struct('df',df,'sdx',std(x),'sdy',std(y));
end