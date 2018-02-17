function [corx,orig,stats] = cnl_tmaxperm(x1,x2,nperm,tail,alpha,m)
%cnl_tmaxperm CNL Toolbox paired permutation test with Tmax correction
%   CORX = CNL_TMAXPERM(X1,X2,NPERM,TAIL,ALPHA,M) returns a structure
%   containing the corrected test statistics CORX for one/paired sample
%   nonparametric permutation tests based on the t-statistic. If X1/X2 are
%   matrices, multiple permutation tests are performed simultaneously and 
%   corrected for multiple comparisons using the Tmax correction method 
%   (Blair & Karnisky, 1993). For one-sample tests, enter data column-wise 
%   in X1 and for paired-sample tests, enter matched data column-wise in X1 
%   and X2 or difference in X1. This function is based on work by Gondan 
%   (2010) and Groppe et al. (2011).
%
%   [...,ORIG] = CNL_TMAXPERM(...) returns a structure containing the
%   original, uncorrected test statistics.
%
%   [...,STATS] = CNL_TMAXPERM(...) returns a structure containing some
%   general data statistics.
%
%   Inputs:
%   x1    - sample 1 data or difference between paired samples 1 and 2 
%           (vector of size obs by 1 or matrix of size obs by var)
%   x2    - sample 2 data (same size as X1)
%   nperm - number of permutations used to estimate the distribution of
%           the null hypothesis (default==10,000)
%   tail  - string specifying alternative hypothesis ('both'==two-tailed
%           test, 'right'==right-tailed test, 'left'==left-tailed test)
%   alpha - significance level (default==0.05)
%   m     - mean of the null hypothesis for each variable (scaler or vector
%           of size 1 by nvar, default==0)
%
%   Outputs:
%   corx  - structure of corrected test statisitcs containing the following
%           fields: 
%           h     - results of permutation tests (can't reject==0, can==1)
%           p     - p-values for each test
%           tcrit - lower and upper critical t-values for given alpha level
%           ci    - 100*(1-ALPHA)% confidence interval for the true mean 
%           estal - estimated alpha level of each test
%   orig  - structure of original, uncorrected test statistics containing 
%           the same fields as CORX
%   stats - structure of data statistics containing the following fields:
%           tstat - t-statistic for each variable
%           df    - degrees of freedom of the test
%           sd    - estimated population standard deviation
%
%   See README for examples of use.
%
%   See also CNL_TMAXPERM2 CNL_RACEMODEL CNL_CALCRMV
%   RMIPERM MULT_COMP_PERM_T1 MULT_COMP_PERM_T2 STATCOND.

%   References:
%      [1] Blair RC, Karniski W (1993) An alternative method for
%          significance testing of waveform difference potentials.
%          Psychophysiology, 30:518-524.
%      [2] Gondan M (2010) A permutation test for the race model
%          inequality. Behav Res Methods, 42(1):23-28.
%      [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis of
%          event-related brain potentials/fields I: A critical tutorial
%          review. Psychophysiology, 48(12):1711-1725.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 02-Feb-2018

if ~exist('x2','var') || isempty(x2)
    x = x1; clear x1 x2
else
    x = x1-x2; clear x1 x2
end
if ~exist('nperm','var') || isempty(nperm)
    nperm = 1e4;
end
if ~exist('tail','var') || isempty(tail)
    tail = 'both';
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
if ~exist('m','var') || isempty(m)
    m = 0;
end

% Remove mean from data
[nobs,nvar] = size(x);
mx = sum(x)/nobs;
if isscalar(m)
    x = x-m;
else
    x = x-repmat(m,nobs,1);
end

% Compute general data statistics
df = nobs-1;
se = std(x)/sqrt(nobs);

% Compute t-statistic
tstat = sum(x)/nobs./se;

% Compute permutations
tp = zeros(nperm,nvar);
sqrtn2 = sqrt(nobs*df);
rsign = sign(rand(nobs,nperm)-0.5);
for i = 1:nperm
    xp = x.*repmat(rsign(:,i),1,nvar); sm = sum(xp);
    tp(i,:) = sm/nobs./(sqrt(sum(xp.^2)-(sm.^2)/nobs)/sqrtn2);
end

% Compute Tmax without sign and add negative values
tmax = max(abs(tp),[],2);
tmax(nperm+1:2*nperm) = -tmax;

% Compute adjusted test statistics using Tmax correction
if strcmp(tail,'both')
    p = mean(abs(tstat)<tmax)*2;
    tcrit(1) = prctile(tmax,100*alpha/2);
    tcrit(2) = -tcrit(1);
    ci = [mx+tcrit(1).*se;mx+tcrit(2).*se];
    estal = mean(tcrit(2)<tmax)*2;
elseif strcmp(tail,'right')
    p = mean(tstat<tmax);
    tcrit = prctile(tmax,100-100*alpha);
    ci = [mx+tcrit.*se;Inf(1,nvar)];
    estal = mean(tcrit<tmax);
elseif strcmp(tail,'left')
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

% Execute if user specifies uncorrected test
if nargout > 1
    
    % Clear variables
    clear h p tcrit ci estal

    % Add negative values
    tp(nperm+1:2*nperm,:) = -tp;

    % Compute original test statistics without correction
    if strcmp(tail,'both')
        p = mean(abs(tstat)<tp)*2;
        tcrit(1,:) = prctile(tp,100*alpha/2);
        tcrit(2,:) = -tcrit(1,:);
        ci = [mx+tcrit(1).*se;mx+tcrit(2).*se];
        estal = mean(tcrit(2,:)<tp)*2;
    elseif strcmp(tail,'right')
        p = mean(tstat<tp);
        tcrit = prctile(tp,100-100*alpha);
        ci = [mx+tcrit.*se;Inf(1,nvar)];
        estal = mean(tcrit<tp);
    elseif strcmp(tail,'left')
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
if nargout > 2
    stats = struct('tstat',tstat,'df',df,'sd',std(x));
    if isscalar(df) && ~isscalar(tstat)
        stats.df = repmat(stats.df,size(tstat));
    end
end