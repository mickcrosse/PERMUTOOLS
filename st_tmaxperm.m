function [tstat,corx,orig,stats] = st_tmaxperm(x1,x2,nperm,tail,alpha,m)
%st_tmaxperm StatsTools one-sample permutation test with Tmax correction
%   TSTAT = ST_TMAXPERM(X1,X2,NPERM,TAIL,ALPHA,M) returns the t-statistic
%   of a one-sample or paired-sample nonparametric permutation test. If X1
%   and X2 are matrices, multiple permutation tests are performed
%   simultaneously between each pair of columns in X1 and X2 and family-
%   wise error rate is controlled using the Tmax correction method (Blair
%   et al., 1994; Westfall & Young, 1993). This method is suitable for
%   multivariate or multiple permutation tests in psychophysics (Gondan,
%   2010) and physiology (Blair & Karniski, 1993; Groppe et al., 2011). For
%   one-sample tests enter the data in X1 and leave X2 empty. For paired-
%   sample tests enter the matched data in X1 and X2.
%
%   [...,CORX] = ST_TMAXPERM(...) returns a structure containing the
%   corrected test statistics of the permutation tests.
%
%   [...,ORIG] = ST_TMAXPERM(...) returns a structure containing the
%   original, uncorrected test statAVistics of the permutation tests.
%
%   [...,STATS] = ST_TMAXPERM(...) returns a structure containing some
%   general data statistics including the effect size.
%
%   Inputs:
%   x1    - column vector or matrix of data (observations by variables)
%   x2    - column vector or matrix of data (observations by variables)
% 
%   Optional Inputs:
%   nperm - scalar specifying the number of permutations (default=10,000) 
%   tail  - string specifying the alternative hypothesis
%           'both'  - mean is not M (two-tailed test, default)
%           'right' - mean is greater than M (right-tailed test)
%           'left'  - mean is less than M (left-tailed test)
%   alpha - scalar between 0 and 1 specifying the significance level as 
%           100*ALPHA% (default=0.05)
%   m     - scalar or row vector specifying the mean of the null hypothesis 
%           for each variable (default=0)
%
%   Outputs:
%   tstat - scalar or vector containing t-statistic of each permutation test
%   corx  - structure of corrected test statistics containing the following
%           fields:
%           h     - test results (H=0: cannot reject, H=1: can reject)
%           p     - probability of observing the given result by chance
%           tcrit - lower and upper critical t-values for given alpha level
%           ci    - 100*(1-ALPHA)% confidence interval for the true mean of
%                   X1 or of X1-X2 for a paired test
%           estal - estimated alpha level of each test
%   orig  - structure of original, uncorrected test statistics containing
%           the same fields as CORX
%   stats - structure of data statistics containing the following fields:
%           df    - degrees of freedom of each test
%           sd    - estimated population standard deviation of X1 or of
%                   X1-X2 for a paired test
%
%   See README for examples of use.
%
%   See also ST_TMAXPERM2 ST_RMAXPERM ST_BOXDOTPLOT ST_PLOTPOLYCI.
%
%   StatsTools https://github.com/mickcrosse/StatsTools

%   References:
%      [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%          Multivariate Permutation Tests Which May Replace Hotelling's T2
%          Test in Prescribed Circumstances. Mult Behav Res, 29(2):141-163.
%      [2] Westfall PH, Young SS (1993) Resampling-Based Multiple Testing:
%           Examples and Methods for p-Value Adjustment. Wiley, New York.
%      [3] Gondan M (2010) A permutation test for the race model inequality
%          Behav Res Methods, 42(1):23-28.
%      [4] Blair RC, Karniski W (1993) An alternative method for
%          significance testing of waveform difference potentials.
%          Psychophysiology, 30:518-524.
%      [5] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis of
%          event-related brain potentials/fields I: A critical tutorial
%          review. Psychophysiology, 48(12):1711-1725.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 20-Mar-2018

if ~exist('x2','var') || isempty(x2)
    x = x1;
elseif size(x1,1)~=size(x2,1)
    error('X1 and X2 must have the same number of observations')
elseif size(x1,2)~=size(x2,2)
    error('X1 and X2 must have the same number of variables')
else
    x = x1-x2;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
elseif alpha<=0 || alpha>=1
    error('Value of ALPHA must be between 0 and 1')
end
if ~exist('nperm','var') || isempty(nperm)
    nperm = 1e4;
elseif nperm<1e3 && alpha<=0.05
    warning('Number of permutations may be too low for chosen alpha level')
elseif nperm<5e3 && alpha<=0.01
    warning('Number of permutations may be too low for chosen alpha level')
end
if ~exist('tail','var') || isempty(tail)
    tail = 'both';
end
if ~exist('m','var') || isempty(m)
    m = 0;
end

% Compute some constants
[nobs,nvar] = size(x); 
df = nobs-1;
mx = sum(x)/nobs;
se = std(x)/sqrt(nobs);
dfp = sqrt(nobs*df);

% Remove mean of null hypothesis from data
if isscalar(m)
    x = x-m;
else
    x = x-repmat(m,nobs,1);
end

% Compute t-statistic
tstat = sum(x)/nobs./se;

% Permute data and generate distribution of t-scores
signx = sign(rand(nobs,nperm)-0.5);
tp = zeros(nperm,nvar);
for i = 1:nperm
    xp = x.*repmat(signx(:,i),1,nvar); sm = sum(xp);
    tp(i,:) = sm/nobs./(sqrt(sum(xp.^2)-(sm.^2)/nobs)/dfp);
end

% Compute Tmax without sign and add negative values
tmax = max(abs(tp),[],2);
tmax(nperm+1:2*nperm) = -tmax;

% Compute corrected test statistics using Tmax correction
if strcmpi(tail,'both')
    p = mean(abs(tstat)<tmax)*2;
    tcrit(1) = prctile(tmax,100*alpha/2);
    tcrit(2) = -tcrit(1);
    ci = [mx+tcrit(1).*se;mx+tcrit(2).*se];
    estal = mean(tcrit(2)<tmax)*2;
elseif strcmpi(tail,'right')
    p = mean(tstat<tmax);
    tcrit = prctile(tmax,100-100*alpha);
    ci = [mx+tcrit.*se;Inf(1,nvar)];
    estal = mean(tcrit<tmax);
elseif strcmpi(tail,'left')
    p = mean(tstat>tmax);
    tcrit = prctile(tmax,100*alpha);
    ci = [-Inf(1,nvar);mx+tcrit.*se];
    estal = mean(tcrit>tmax);
end

% Determine if adjusted p-values exceed desired alpha level
h = cast(p<alpha,'like',p);
h(isnan(tcrit)) = NaN;
p(isnan(tcrit)) = NaN;

% Store values in structure
corx = struct('h',h,'p',p,'tcrit',tcrit,'ci',ci,'estal',estal);

% Execute if user specifies uncorrected test statistics
if nargout > 2
    
    % Clear variables
    clear h p tcrit ci estal
    
    % Add negative values
    tp(nperm+1:2*nperm,:) = -tp;
    
    % Compute original test statistics without correction
    if strcmpi(tail,'both')
        p = mean(abs(tstat)<tp)*2;
        tcrit(1,:) = prctile(tp,100*alpha/2);
        tcrit(2,:) = -tcrit(1,:);
        ci = [mx+tcrit(1).*se;mx+tcrit(2).*se];
        estal = mean(tcrit(2,:)<tp)*2;
    elseif strcmpi(tail,'right')
        p = mean(tstat<tp);
        tcrit = prctile(tp,100-100*alpha);
        ci = [mx+tcrit.*se;Inf(1,nvar)];
        estal = mean(tcrit<tp);
    elseif strcmpi(tail,'left')
        p = mean(tstat>tp);
        tcrit = prctile(tp,100*alpha);
        ci = [-Inf(1,nvar);mx+tcrit.*se];
        estal = mean(tcrit>tp);
    end
    
    % Determine if original p-values exceed desired alpha level
    h = cast(p<alpha,'like',p);
    h(isnan(tcrit(1,:))) = NaN;
    p(isnan(tcrit(1,:))) = NaN;
    
    % Store values in structure
    orig = struct('h',h,'p',p,'tcrit',tcrit,'ci',ci,'estal',estal);
    
end

% Execute if user specifies general data statistics
if nargout > 3
    if nvar>1
        df = repmat(df,1,nvar);
    end
    stats = struct('df',df,'sd',std(x));
end