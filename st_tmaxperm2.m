function [tstat,corx,orig,stats] = st_tmaxperm2(x1,x2,nperm,tail,alpha,varx)
%st_tmaxperm2 StatsTools two-sample permutation test with Tmax correction
%   TSTAT = ST_TMAXPERM2(X1,X2,NPERM,TAIL,ALPHA,VARX) returns the t-statistic
%   of an unpaired two-sample nonparametric permutation test. If X1 and X2
%   are matrices, multiple permutation tests are performed simultaneously
%   between each pair of columns in X1 and X2 and family-wise error rate is
%   controlled using the Tmax correction method (Blair et al., 1994;
%   Westfall & Young, 1993). This method is suitable for multivariate or
%   multiple permutation tests in psychophysics (Gondan, 2010) and
%   psychophysiology (Blair & Karniski, 1993; Groppe et al., 2011).
%
%   [...,CORX] = ST_TMAXPERM(...) returns a structure containing the
%   corrected test statistics of the permutation tests.
%
%   [...,ORIG] = ST_TMAXPERM2(...) returns a structure containing the
%   original, uncorrected test statistics of the permutation tests.
%
%   [...,STATS] = ST_TMAXPERM2(...) returns a structure containing some
%   general data statistics including the effect size.
%
%   Inputs:
%   x1    - column vector or matrix of data (observations by variables)
%   x2    - column vector or matrix of data (observations by variables)
% 
%   Optional Inputs:
%   nperm - scalar specifying the number of permutations (default=10,000) 
%   tail  - string specifying the alternative hypothesis
%           'both'  - means are not equal (two-tailed test, default)
%           'right' - mean of X1 greater than mean of X2 (right-tailed test)
%           'left'  - mean of X1 less than mean of X2 (left-tailed test)
%   alpha - scalar between 0 and 1 specifying the significance level as 
%           100*ALPHA% (default=0.05)
%   varx  - string specifying the variance equivalence of X1 and X2
%           'equal'   - assume samples have equal variances (default)
%           'unequal' - assume samples have unequal variances
%
%   Outputs:
%   tstat - scalar or vector containing t-statistic of each permutation test
%   corx  - structure of corrected test statistics containing the following
%           fields:
%           h     - test results (H=0: cannot reject, H=1: can reject)
%           p     - probability of observing the given result by chance
%           tcrit - lower and upper critical t-values for given alpha level
%           ci    - 100*(1-ALPHA)% confidence interval for the true
%                   difference of population means
%           estal - estimated alpha level of each test
%   orig  - structure of original, uncorrected test statistics containing
%           the same fields as CORX
%   stats - structure of data statistics containing the following fields:
% 
%           df    - degrees of freedom of each test
%           sd    - pooled estimate (equal variances) or unpooled estimates
%                   (unequal variances) of population standard deviation
%
%   See README for examples of use.
%
%   See also ST_TMAXPERM ST_RMAXPERM ST_BOXDOTPLOT ST_PLOTPOLYCI.
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
    error('st_tmaxperm2: Requires at least two input arguments')
elseif size(x1,2)~=size(x2,2)
    error('X1 and X2 must have the same number of variables')
else
    x = [x1;x2];
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
if ~exist('varx','var') || isempty(varx)
    varx = 'equal';
end

% Compute some constants
[nobs,nvar] = size(x);
nobs1 = size(x1,1); df1 = nobs1-1;
nobs2 = size(x2,1); df2 = nobs2-1;
dfp = nobs/(nobs1*nobs2);

% Compute t-statistic
sm1 = sum(x1); sm2 = sum(x2);
diffmx = sm1/nobs1-sm2/nobs2;
var1 = (sum(x1.^2)-(sm1.^2)/nobs1)/df1;
var2 = (sum(x2.^2)-(sm2.^2)/nobs2)/df2;
if strcmpi(varx,'equal')
    df = nobs-2;
    sd = sqrt((df1.*var1+df2.*var2)./df);
    se = sd.*sqrt(dfp);
elseif strcmpi(varx,'unequal')
    se2x1 = var1/nobs1;
    se2x2 = var2/nobs2;
    df = (se2x1+se2x2).^2./(se2x1.^2/df1+se2x2.^2/df2);
    se = sqrt(se2x1+se2x2);
end
tstat = diffmx./se;

% Permute data and generate distribution of t-scores
[~,idx] = sort(rand(nobs,nperm));
tp = zeros(nperm,nvar);
for i = 1:nperm
    xp1 = x(idx(1:nobs1,i),:); sm1 = sum(xp1);
    xp2 = x(idx(nobs1+1:nobs,i),:); sm2 = sum(xp2);
    s2x1 = (sum(xp1.^2)-(sm1.^2)/nobs1)/df1;
    s2x2 = (sum(xp2.^2)-(sm2.^2)/nobs2)/df2;
    if strcmpi(varx,'equal')
        se = sqrt((df1.*s2x1+df2.*s2x2)./df).*sqrt(dfp);
    elseif strcmpi(varx,'unequal')
        se = sqrt(s2x1/nobs1+s2x2/nobs2);
    end
    tp(i,:) = (sm1/nobs1-sm2/nobs2)./se;
end

% Compute Tmax with sign
[~,idx] = max(abs(tp),[],2);
csvar = [0;cumsum(ones(nperm-1,1)*nvar)];
tpT = tp';
tmax = tpT(idx+csvar);

% Compute corrected test statistics using Tmax correction
if strcmpi(tail,'both')
    tmax = abs(tmax);
    p = mean(abs(tstat)<tmax);
    tcrit(2) = prctile(tmax,100-100*alpha);
    tcrit(1) = -tcrit(2);
    ci = [diffmx+tcrit(1).*se;diffmx+tcrit(2).*se];
    estal = mean(tcrit(2)<tmax);
elseif strcmpi(tail,'right')
    p = mean(tstat<tmax);
    tcrit = prctile(tmax,100-100*alpha);
    ci = [diffmx+tcrit.*se;Inf(1,nvar)];
    estal = mean(tcrit<tmax);
elseif strcmpi(tail,'left')
    p = mean(tstat>tmax);
    tcrit = prctile(tmax,100*alpha);
    ci = [-Inf(1,nvar);diffmx+tcrit.*se];
    estal = mean(tcrit>tmax);
end

% Determine if adjusted p-values exceed desired alpha level
h = cast(p<alpha,'like',p);
h(isnan(tcrit)) = NaN;
p(isnan(tcrit)) = NaN;

% Store values in structure
corx = struct('h',h,'p',p,'tcrit',tcrit,'ci',ci,'estal',estal);

% Execute if user specifies uncorrected test
if nargout > 2
    
    % Clear variables
    clear h p tcrit estal
    
    % Compute original test statistics without correction
    if strcmpi(tail,'both')
        tpmag = abs(tp);
        p = mean(abs(tstat)<tpmag);
        tcrit(1,:) = prctile(tpmag,100-100*alpha);
        tcrit(2,:) = -tcrit(1,:);
        ci = [diffmx+tcrit(1).*se;diffmx+tcrit(2).*se];
        estal = mean(tcrit(2,:)<tpmag);
    elseif strcmpi(tail,'right')
        p = mean(tstat<tp);
        tcrit = prctile(tp,100-100*alpha);
        ci = [diffmx+tcrit.*se;Inf(1,nvar)];
        estal = mean(tcrit<tp);
    elseif strcmpi(tail,'left')
        p = mean(tstat>tp);
        tcrit = prctile(tp,100*alpha);
        ci = [-Inf(1,nvar);diffmx+tcrit.*se];
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
    if strcmpi(varx,'equal')
        stats = struct('df',df,'sd',sd);
    elseif strcmpi(varx,'unequal')
        stats = struct('df',df,'sd',sqrt([var1;var2]));
    end
end