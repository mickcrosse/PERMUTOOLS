function [corx,orig,stats] = cnl_tmaxperm2(x1,x2,nperm,tail,alpha,varx)
%cnl_tmaxperm2 CNL Toolbox two-sample permutation test with Tmax correction
%   CORX = CNL_TMAXPERM2(X1,X2,NPERM,TAIL,ALPHA,VARX) returns a structure
%   containing the test statistics of a two-sample nonparametric 
%   permutation test based on the t-statistic. If X1 and X2 are matrices, 
%   multiple permutation tests are performed simultaneously between each 
%   pair of columns in X1 and X2 and family-wise error rate is controlled
%   using the Tmax correction method (Blair & Karnisky, 1993). Calculations 
%   are based on Gondan (2010) and Groppe et al. (2011).
%
%   [...,ORIG] = CNL_TMAXPERM2(...) returns a structure containing the
%   original, uncorrected test statistics of the permutation tests.
%
%   [...,STATS] = CNL_TMAXPERM2(...) returns a structure containing some
%   general data statistics.
%
%   Inputs:
%   x1    - column vector or matrix of data (observations by measures)
%   x2    - column vector or matrix of data (observations by measures)
%   nperm - number of permutations (default=10,000)
%   tail  - string specifying the alternative hypothesis
%           'both'  - means are not equal (two-tailed test, default)
%           'right' - mean of X1 greater than mean of X2 (right-tailed test)
%           'left'  - mean of X1 less than mean of X2 (left-tailed test)
%   alpha - significance level between 0 and 1 (default=0.05)
%   varx  - variance equivalence of independent samples
%           'equal'   - assume samples have equal variances (default)
%           'unequal' - assume samples have unequal variances
%
%   Outputs:
%   corx  - structure of corrected test statistics containing the following
%           fields:
%           h     - test results (H=0: cannot reject, H=1: can reject)
%           p     - probability of observing the given result by chance
%           tcrit - lower and upper critical t-values for given alpha level
%           estal - estimated alpha level of each test
%   orig  - structure of original, uncorrected test statistics containing
%           the same fields as CORX
%   stats - structure of data statistics containing the following fields:
%           tstat - t-statistic of each test
%           df    - degrees of freedom of each test
%           sd    - estimated population standard deviation
%
%   See README for examples of use.
%
%   See also CNL_TMAXPERM CNL_TMAXCORR CNL_RACEMODEL CNL_CALCRMV
%   RMIPERM MULT_COMP_PERM_T1 MULT_COMP_PERM_T2 STATCOND.

%   References:
%      [1] Blair RC, Karniski W (1993) An alternative method for
%          significance testing of waveform difference potentials.
%          Psychophysiology, 30:518-524.
%      [2] Gondan M (2010) A permutation test for the race model inequality
%          Behav Res Methods, 42(1):23-28.
%      [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis of
%          event-related brain potentials/fields I: A critical tutorial
%          review. Psychophysiology, 48(12):1711-1725.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Feb 2018; Last Revision: 21-Feb-2018

if ~exist('x2','var') || isempty(x2)
    error('cnl_tmaxperm2: Requires at least two input arguments');
elseif size(x1,2)~=size(x2,2)
    error('X1 and X2 must have the same number of measures')
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
if ~exist('vartype','var') || isempty(varx)
    varx = 'equal';
end

% Compute some constants
[nobs,nmeas] = size(x); df = nobs-2;
nobs1 = size(x1,1); df1 = nobs1-1; 
nobs2 = size(x2,1); df2 = nobs2-1;
dfp = nobs/(nobs1*nobs2);

% Compute t-statistic
sm1 = sum(x1); ss1 = sum(x1.^2)-(sm1.^2)/nobs1;
sm2 = sum(x2); ss2 = sum(x2.^2)-(sm2.^2)/nobs2;
if strcmpi(varx,'equal')
    tstat = (sm1/nobs1-sm2/nobs2)./sqrt((ss1+ss2)/df*dfp);
elseif strcmpi(varx,'unequal')
    tstat = sm1/nobs1./sqrt(ss1/df1)-sm2/nobs2./sqrt(ss2/df2);
end

% Permute data and generate distribution of t-scores 
[~,idx] = sort(rand(nobs,nperm));
tp = zeros(nperm,nmeas);
for i = 1:nperm
    xp1 = x(idx(1:nobs1,i),:);
    xp2 = x(idx(nobs1+1:nobs,i),:);
    sm1 = sum(xp1); ss1 = sum(xp1.^2)-(sm1.^2)/nobs1;
    sm2 = sum(xp2); ss2 = sum(xp2.^2)-(sm2.^2)/nobs2;
    if strcmpi(varx,'equal')
        tp(i,:) = (sm1/nobs1-sm2/nobs2)./sqrt((ss1+ss2)/df*dfp);
    elseif strcmpi(varx,'unequal')
        tp(i,:) = sm1/nobs1./sqrt(ss1/df1)-sm2/nobs2./sqrt(ss2/df2);
    end
end

% Compute Tmax with sign and don't assume symmetric distribution
[~,idx] = max(abs(tp),[],2);
csvar = [0;cumsum(ones(nperm-1,1)*nmeas)];
tpT = tp';
tmax = tpT(idx+csvar);

% Compute adjusted test statistics using Tmax correction
if strcmpi(tail,'both')
    tmax = abs(tmax);
    p = mean(abs(tstat)<tmax);
    tcrit(2) = prctile(tmax,100-100*alpha);
    tcrit(1) = -tcrit(2);
    estal = mean(tcrit(2)<tmax);
elseif strcmpi(tail,'right')
    p = mean(tstat<tmax);
    tcrit = prctile(tmax,100-100*alpha);
    estal = mean(tcrit<tmax);
elseif strcmpi(tail,'left')
    p = mean(tstat>tmax);
    tcrit = prctile(tmax,100*alpha);
    estal = mean(tcrit>tmax);
end

% Determine if adjusted p-values exceed desired alpha level
h = cast(p<alpha,'like',p);
h(isnan(tcrit)) = NaN;
p(isnan(tcrit)) = NaN;

% Store values in structure
corx = struct('h',h,'p',p,'tcrit',tcrit,'estal',estal);

% Execute if user specifies uncorrected test
if nargout > 1
    
    % Clear variables
    clear h p tcrit estal
    
    % Compute original test statistics without correction
    if strcmpi(tail,'both')
        tpmag = abs(tp);
        p = mean(abs(tstat)<tpmag);
        tcrit(1,:) = prctile(tpmag,100-100*alpha);
        tcrit(2,:) = -tcrit(1,:);
        estal = mean(tcrit(2,:)<tpmag);
    elseif strcmpi(tail,'right')
        p = mean(tstat<tp);
        tcrit = prctile(tp,100-100*alpha);
        estal = mean(tcrit<tp);
    elseif strcmpi(tail,'left')
        p = mean(tstat>tp);
        tcrit = prctile(tp,100*alpha);
        estal = mean(tcrit>tp);
    end
    
    % Determine if original p-values exceed desired alpha level
    h = cast(p<alpha,'like',p);
    h(isnan(tcrit(1,:))) = NaN;
    p(isnan(tcrit(1,:))) = NaN;
    
    % Store values in structure
    orig = struct('h',h,'p',p,'tcrit',tcrit,'estal',estal);
    
end

% Execute if user specifies general data statistics
if nargout > 2
    if isscalar(df) && nmeas>1
        df = repmat(df,1,nmeas);
    end
    stats = struct('tstat',tstat,'df',df,'sd',std(x));
end