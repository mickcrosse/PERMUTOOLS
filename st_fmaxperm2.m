function [fstat,corx,orig,stats] = st_fmaxperm2(x1,x2,nperm,tail,alpha)
%st_fmaxperm2 StatsTools two-sample F permutation test with Fmax correction
%   FSTAT = ST_FMAXPERM2(X1,X2,NPERM,TAIL,ALPHA) returns the F-statistic
%   of an unpaired two-sample nonparametric permutation test. If X1 and X2
%   are matrices, multiple permutation tests are performed simultaneously
%   between each pair of columns in X1 and X2 and family-wise error rate is
%   controlled using the Fmax correction method (Blair et al., 1994;
%   Westfall & Young, 1993). This method is suitable for multivariate or
%   multiple permutation tests in psychophysics (Gondan, 2010) and
%   psychophysiology (Blair & Karniski, 1993; Groppe et al., 2011).
%
%   [...,CORX] = ST_FMAXPERM2(...) returns a structure containing the
%   corrected test statistics of the permutation tests.
%
%   [...,ORIG] = ST_FMAXPERM2(...) returns a structure containing the
%   original, uncorrected test statistics of the permutation tests.
%
%   [...,STATS] = ST_FMAXPERM2(...) returns a structure containing some
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
%   fstat - scalar or vector containing t-statistic of each permutation test
%   corx  - structure of corrected test statistics containing the following
%           fields:
%           h     - test results (H=0: cannot reject, H=1: can reject)
%           p     - probability of observing the given result by chance
%           fcrit - lower and upper critical F-values for given alpha level
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
%   See also ST_TMAXPERM ST_TMAXPERM2 ST_RMAXPERM ST_BOXDOTPLOT.
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

% Compute some constants
[nobs,nvar] = size(x);
nobs1 = size(x1,1); df1 = nobs1-1;
nobs2 = size(x2,1); df2 = nobs2-1;

% Compute F-statistic
var1 = (sum(x1.^2)-(sum(x1).^2)/nobs1)/df1;
var2 = (sum(x2.^2)-(sum(x2).^2)/nobs2)/df2;
fstat = var1./var2;

% Permute data and generate distribution of F-scores
[~,idx] = sort(rand(nobs,nperm));
fp = zeros(nperm,nvar);
for i = 1:nperm
    xp1 = x(idx(1:nobs1,i),:);
    xp2 = x(idx(nobs1+1:nobs,i),:);
    var1 = (sum(xp1.^2)-(sum(xp1).^2)/nobs1)/df1;
    var2 = (sum(xp2.^2)-(sum(xp2).^2)/nobs2)/df2;
    fp(i,:) = var1./var2;
end

% Compute Fmax with sign
[~,idx] = max(abs(fp),[],2);
csvar = [0;cumsum(ones(nperm-1,1)*nvar)];
fpT = fp';
fmax = fpT(idx+csvar);

% Compute corrected test statistics using Tmax correction
if strcmpi(tail,'both')
    p = mean(fstat<fmax)*2;
    fcrit(1) = prctile(fmax,100*alpha/2);
    fcrit(2) = prctile(fmax,100-100*alpha/2);
    ci = [fstat.*finv(alpha/2,df2,df1);fstat./finv(alpha/2,df1,df2)];
    estal = mean(fcrit(2)<fmax)+mean(fcrit(1)>fmax);
elseif strcmpi(tail,'right')
    p = mean(fstat<fmax);
    fcrit = prctile(fmax,100-100*alpha);
    ci = [fstat.*finv(alpha/2,df2,df1);Inf(1,nvar)];
    estal = mean(fcrit<fmax);
elseif strcmpi(tail,'left')
    p = mean(fstat>fmax);
    fcrit = prctile(fmax,100*alpha);
    ci = [zeros(1,nvar);fstat./finv(alpha/2,df1,df2)];
    estal = mean(fcrit>fmax);
end

% Determine if adjusted p-values exceed desired alpha level
h = cast(p<alpha,'like',p);
h(isnan(fcrit)) = NaN;
p(isnan(fcrit)) = NaN;

% Store values in structure
corx = struct('h',h,'p',p,'fcrit',fcrit,'ci',ci,'estal',estal);

% Execute if user specifies uncorrected test
if nargout > 2
    
    % Clear variables
    clear h p fcrit estal
    
    % Compute original test statistics without correction
    if strcmpi(tail,'both')
        p = mean(fstat<fp)*2;
        fcrit(1,:) = prctile(fp,100*alpha/2);
        fcrit(2,:) = prctile(fp,100-100*alpha/2);
        ci = [fstat.*finv(alpha/2,df2,df1);fstat./finv(alpha/2,df1,df2)];
        estal = mean(fcrit(2)<fp)+mean(fcrit(1)>fp);
    elseif strcmpi(tail,'right')
        p = mean(fstat<fp);
        fcrit = prctile(fp,100-100*alpha);
        ci = [fstat.*finv(alpha/2,df2,df1);Inf(1,nvar)];
        estal = mean(fcrit<fp);
    elseif strcmpi(tail,'left')
        p = mean(fstat>fp);
        fcrit = prctile(fp,100*alpha);
        ci = [zeros(1,nvar);fstat./finv(alpha/2,df1,df2)];
        estal = mean(fcrit>fp);
    end
    
    % Determine if original p-values exceed desired alpha level
    h = cast(p<alpha,'like',p);
    h(isnan(fcrit(1,:))) = NaN;
    p(isnan(fcrit(1,:))) = NaN;
    
    % Store values in structure
    orig = struct('h',h,'p',p,'fcrit',fcrit,'ci',ci,'estal',estal);
    
end

% Execute if user specifies general data statistics
if nargout > 3
    if nvar>1
        df1 = repmat(df1,1,nvar);
        df2 = repmat(df2,1,nvar);
    end
    stats = struct('df1',df1,'df2',df2,'sd1',sqrt(var1),'sd2',sqrt(var2));
end