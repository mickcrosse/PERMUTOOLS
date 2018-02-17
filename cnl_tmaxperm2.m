function [orig,corx,stats] = cnl_tmaxperm2(x1,x2,nperm,tail,alpha,varx)
%cnl_tmaxperm CNL Toolbox unpaired permutation test with Tmax correction
%   ORIG = CNL_TMAXPERM(X1,X2,NPERM,TAIL,ALPHA) returns a structure
%   containing the original test statistics ORIG for unpaired sample
%   nonparametric permutation tests based on the t-statistic. Enter 
%   independent samples column-wise in X1 and X2. If X1/X2 are matrices, 
%   multiple permutation tests are performed simultaneously. This function 
%   is based on work by Gondan (2010) and Groppe et al. (2011).
%
%   [...,CORX] = CNL_TMAXPERM(...) returns a structure containing the
%   adjusted test statistics using the Tmax correction method (Blair &
%   Karnisky, 1993).
%
%   [...,STATS] = CNL_TMAXPERM(...) returns a structure containing the
%   general data statistics.
%
%   Inputs:
%   x1    - matrix of data for sample 1 (obs by var)
%   x2    - matrix of data for sample 2 (obs by var)
%   nperm - number of permutations used to estimate the distribution of
%           the null hypothesis (default==10,000)
%   tail  - string specifying alternative hypothesis ('both'==two-tailed
%           test, 'right'==right-tailed test, 'left'==left-tailed test).
%   alpha - significance level (default==0.05)
%   varx  - variance equivalance of X1 and X2
%
%   Outputs:
%   orig  - structure of original test statistics containing the following
%           fields:
%           p     - p-values for each test
%           tcrit - lower and upper critical t-values for given alpha level
%           estal - estimated alpha level of each test
%           h     - results of permutation tests (can't reject==0, can==1)
%   corx  - structure of adjusted test statisitcs using the tmax correction
%           method (contains the same fields as ORIG)
%   stats - structure of data statistics containing the following fields:
%           tstat - t-statistic for each variable
%           df    - degrees of freedom of the test
%           sd    - estimated population standard deviation
%
%   See README for examples of use.
%
%   See also CNL_TMAXPERM CNL_RACEMODEL CNL_CALCRMV
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
%   Jan 2018; Last Revision: 29-Jan-2018

if ~exist('x2','var') || isempty(x2)
    error(message('stats:ttest2:TooFewInputs'));
else
    x = [x1;x2];
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
if ~exist('vartype','var') || isempty(varx)
    varx = 'equal';
end

% Compute general data statistics
[nobs1,nvar] = size(x1);
[nobs2] = size(x2,1);
nobs = nobs1+nobs2;
smprod = nobs/(nobs1*nobs2);
df1 = nobs1-1; df2 = nobs2-1;
stats.df = nobs-2;
stats.sd = std(x);

% Compute t-statistic
sm1 = sum(x1); ss1 = sum(x1.^2)-(sm1.^2)/nobs1;
sm2 = sum(x2); ss2 = sum(x2.^2)-(sm2.^2)/nobs2;
stats.tstat = (sm1/nobs1-sm2/nobs2)./sqrt((ss1+ss2)/stats.df*smprod);

% Compute permutations
tp = zeros(nperm,nvar);
[~,ridx] = sort(rand(nobs,nperm));
for i = 1:nperm
    xp1 = x(ridx(1:nobs1,i),:);
    xp2 = x(ridx(nobs1+1:nobs,i),:);
    sm1 = sum(xp1); ss1 = sum(xp1.^2)-(sm1.^2)/nobs1;
    sm2 = sum(xp2); ss2 = sum(xp2.^2)-(sm2.^2)/nobs2;
    if strcmp(varx,'equal')
        tp(i,:) = (sm1/nobs1-sm2/nobs2)./sqrt((ss1+ss2)/stats.df*smprod);
    elseif srtcmp(varx,'unequal')
        tp(i,:) = sm1/nobs1./sqrt(ss1/df1)-sm2/nobs2./sqrt(ss2/df2);
    end
end

% Compute original test statistics without correction
if strcmp(tail,'both')
    tpmag = abs(tp);  
    orig.p = mean(abs(stats.tstat)<tpmag);
    orig.tcrit(1,:) = prctile(tpmag,100-100*alpha);
    orig.tcrit(2,:) = -orig.tcrit(1,:);
    orig.estal = mean(orig.tcrit(2,:)<tpmag);
%     orig.p = mean(abs(stats.tstat)<tp)*2;
%     orig.tcrit(1,:) = prctile(tp,100*alpha/2);
%     orig.tcrit(2,:) = -orig.tcrit(1,:);
%     orig.estal = mean(orig.tcrit(2,:)<tp)*2;
elseif strcmp(tail,'right')
    orig.p = mean(stats.tstat<tp);
    orig.tcrit = prctile(tp,100-100*alpha);
    orig.estal = mean(orig.tcrit<tp);
elseif strcmp(tail,'left')
    orig.p = mean(stats.tstat>tp);
    orig.tcrit = prctile(tp,100*alpha);
    orig.estal = mean(orig.tcrit>tp);
end

% Determine if original p-values exceed desired alpha level
orig.h = cast(orig.p<alpha,'like',orig.p);
orig.p(isnan(orig.tcrit(1,:))) = NaN;
orig.h(isnan(orig.p)) = NaN;

% Execute if user specifies Tmax corrected output
if nargout > 1
    
    % Compute Tmax with sign (don't assume symmetric distribution)
    [~,idx] = max(abs(tp),[],2);
    csvar = [0;cumsum(ones(nperm-1,1)*nvar)];
    tpT = tp';
    tmax = tpT(idx+csvar);

    % Compute adjusted test statistics using Tmax correction
    if strcmp(tail,'both')
        tmax = abs(tmax);   
        corx.p = mean(abs(stats.tstat)<tmax);
        corx.tcrit(2) = prctile(tmax,100-100*alpha);
        corx.tcrit(1) = -corx.tcrit(2);
        corx.estal = mean(corx.tcrit(2)<tmax);
%         corx.p = mean(abs(stats.tstat)<tmax)*2;
%         corx.tcrit(1) = prctile(tmax,100*alpha/2);
%         corx.tcrit(2) = -corx.tcrit(1);
%         corx.estal = mean(corx.tcrit(2)<tmax)*2;
    elseif strcmp(tail,'right')
        corx.p = mean(stats.tstat<tmax);
        corx.tcrit = prctile(tmax,100-100*alpha);
        corx.estal = mean(corx.tcrit<tmax);
    elseif strcmp(tail,'left')
        corx.p = mean(stats.tstat>tmax);
        corx.tcrit = prctile(tmax,100*alpha);
        corx.estal = mean(corx.tcrit>tmax);
    end

    % Determine if adjusted p-values exceed desired alpha level
    corx.h = cast(corx.p<alpha,'like',corx.p);
    corx.p(isnan(corx.tcrit)) = NaN;
    corx.h(isnan(corx.p)) = NaN;

end