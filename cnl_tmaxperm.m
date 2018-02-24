function [corx,orig,stats] = cnl_tmaxperm(x1,x2,nperm,tail,alpha,m)
%cnl_tmaxperm CNL Toolbox one-sample permutation test with Tmax correction
%   CORX = CNL_TMAXPERM(X1,X2,NPERM,TAIL,ALPHA,M) returns a structure
%   containing the test statistics of a one-sample or paired-sample 
%   nonparametric permutation test based on the t-statistic. If X1 and X2 
%   are matrices, multiple permutation tests are performed simultaneously 
%   between each pair of columns in X1 and X2 and family-wise error rate is 
%   controlled using the Tmax correction method (Blair & Karnisky, 1993). 
%   Calculations based on Gondan (2010) and Groppe et al. (2011). For one-
%   sample tests, enter data in X1 and leave X2 empty. For paired-sample 
%   tests, enter matched data in X1 and X2.
%
%   [...,ORIG] = CNL_TMAXPERM(...) returns a structure containing the
%   original, uncorrected test statistics.
%
%   [...,STATS] = CNL_TMAXPERM(...) returns a structure containing some
%   general data statistics.
%
%   Inputs:
%   x1    - column vector or matrix of data (observations by measures)
%   x2    - column vector or matrix of data (observations by measures)
%   nperm - number of permutations (default=10,000)
%   tail  - string specifying the alternative hypothesis
%           'both'  - mean is not M (two-tailed test, default)
%           'right' - mean is greater than M (right-tailed test)
%           'left'  - mean is less than M (left-tailed test)
%   alpha - significance level between 0 and 1 (default=0.05)
%   m     - scalar or row vector of the mean of the null hypothesis for
%           each measure (default=0)
%
%   Outputs:
%   corx  - structure of corrected test statistics containing the following
%           fields:
%           h     - test results (H=0: cannot reject, H=1: can reject)
%           p     - probability of observing the given result by chance
%           tcrit - lower and upper critical t-values for given alpha level
%           ci    - 100*(1-ALPHA)% confidence interval for the true mean
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
%   See also CNL_TMAXPERM2 CNL_TMAXCORR CNL_RACEMODEL CNL_CALCRMV
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
%   Jan 2018; Last Revision: 02-Feb-2018

if ~exist('x2','var') || isempty(x2)
    x = x1;
elseif size(x1,2)~=size(x2,2)
    error('X1 and X2 must have the same number of measures')
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
[nobs,nmeas] = size(x); df = nobs-1;
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
tp = zeros(nperm,nmeas);
for i = 1:nperm
    xp = x.*repmat(signx(:,i),1,nmeas); sm = sum(xp);
    tp(i,:) = sm/nobs./(sqrt(sum(xp.^2)-(sm.^2)/nobs)/dfp);
end

% Compute Tmax without sign and add negative values
tmax = max(abs(tp),[],2);
tmax(nperm+1:2*nperm) = -tmax;

% Compute adjusted test statistics using Tmax correction
if strcmpi(tail,'both')
    p = mean(abs(tstat)<tmax)*2;
    tcrit(1) = prctile(tmax,100*alpha/2);
    tcrit(2) = -tcrit(1);
    ci = [mx+tcrit(1).*se;mx+tcrit(2).*se];
    estal = mean(tcrit(2)<tmax)*2;
elseif strcmpi(tail,'right')
    p = mean(tstat<tmax);
    tcrit = prctile(tmax,100-100*alpha);
    ci = [mx+tcrit.*se;Inf(1,nmeas)];
    estal = mean(tcrit<tmax);
elseif strcmpi(tail,'left')
    p = mean(tstat>tmax);
    tcrit = prctile(tmax,100*alpha);
    ci = [-Inf(1,nmeas);mx+tcrit.*se];
    estal = mean(tcrit>tmax);
end

% Determine if adjusted p-values exceed desired alpha level
h = cast(p<alpha,'like',p);
h(isnan(tcrit)) = NaN;
p(isnan(tcrit)) = NaN;

% Store values in structure
corx = struct('h',h,'p',p,'tcrit',tcrit,'ci',ci,'estal',estal);

% Execute if user specifies uncorrected test statistics
if nargout > 1
    
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
        ci = [mx+tcrit.*se;Inf(1,nmeas)];
        estal = mean(tcrit<tp);
    elseif strcmpi(tail,'left')
        p = mean(tstat>tp);
        tcrit = prctile(tp,100*alpha);
        ci = [-Inf(1,nmeas);mx+tcrit.*se];
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
if nargout > 2
    if isscalar(df) && nmeas>1
        df = repmat(df,1,nmeas);
    end
    stats = struct('tstat',tstat,'df',df,'sd',std(x));
end