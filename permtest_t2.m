function [t,stats,orig,param] = permtest_t2(x1,x2,varargin)
%permtest_t2 unpaired two-sample permutation test with Tmax correction
%   T = PERMTEST_T2(X1,X2) returns the t-statistic of a two-sample
%   permutation test. If X1 and X2 are matrices, multiple permutation tests
%   are performed simultaneously between each corresponding pair of columns
%   in X1 and X2 and family-wise error rate (FWER) is controlled using the
%   Tmax correction method (Blair et al., 1994). This method provides
%   strong control of FWER, even for small sample sizes, and is much more
%   powerful than traditional correction methods (Groppe et al., 2011a). It
%   is also rather insensitive to differences in population variance when
%   samples of equal size are used (Groppe et al., 2011b). For samples of
%   unequal size or variance, Welch's t-statistic may be used as it is less
%   sensitive to differences in variance (but also less sensitive to
%   differences in means). If X2 is not entered, a permutation test between
%   each pair of columns in X1 is performed and output as a matrix. Samples
%   of different sizes may be used by replacing missing values with NaNs.
%   This function treats NaNs as missing values, and ignores them.
%
%   [...,STATS] = permtest_t2(...) returns the adjusted test statistics in
%   a structure containing the following fields:
%       h       test results: H=0 indicates the null hypothesis cannot be
%               rejected, H=1 indicates the null hypothesis can be rejected
%       p       the probability of observing the given result by chance
%       ci    	100*(1-ALPHA)% confidence interval for the true difference
%               of sample means
%       tcrit 	critical t-value for the given alpha level (for two-tailed
%               tests, the lower t-value is equal to -1*TCRIT)
%       estal 	the estimated alpha level of each test
%
%   [...,ORIG] = permtest_t2(...) returns the original, unadjusted test
%   statistics in a structure containing the same fields as STATS.
%
%   [...,PARAM] = permtest_t2(...) returns some statistical parameters in a
%   structure containing the following fields:
%       df      the degrees of freedom of each test
%       sd    	pooled estimate (equal variances) or unpooled estimates
%               (unequal variances) of population standard deviation
%
%   [...] = permtest_t2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'alpha'     a scalar between 0 and 1 specifying the significance level
%               as 100*ALPHA% (default=0.05)
%   'nperm'     a scalar specifying the number of permutations (default=
%               10,000)
%   'tail'      a string specifying the alternative hypothesis
%                   'both'      means are not equal (default)
%                   'right'     mean of X1 is greater than mean of X2
%                   'left'      mean of X1 is less than mean of X2
%   'rows'      a string specifying the rows to use in the case of any
%               missing values (NaNs)
%                   'all'       use all rows, and ignore any NaNs (default)
%                   'complete'  use only rows with no NaNs
%   'varx'      a string specifying the variance equivalence of X1 and X2
%               to determine the SD and test statistic estimation method
%                   'equal'   	assume equal variances (default)
%                   'unequal' 	assume unequal variances
%
%   Example 1: generate multivariate data for 2 samples, each with 20
%   variables and 30 observations and perform unpaired permutation tests
%   between the corresponding variables of each sample.
%       x1 = randn(30,20);
%       x2 = randn(30,20);
%       x2(:,1:8) = x2(:,1:8)-1;
%       [t,stats] = permtest_t2(x1,x2)
%
%   Example 2: generate univariate data for 5 samples, each with 30
%   observations and perform unpaired permutation tests between every
%   sample (5 samples = 10 comparisons). Note that each column of X
%   represents an independent sample and may contain NaNs for samples with
%   smaller number of observations.
%       x = randn(30,5);
%       x(:,3:5) = x(:,3:5)-1;
%       [t,stats] = permtest_t2(x)
%
%   See also PERMTEST_T PERMTEST_F PERMTEST_CORR EFFECTSIZE_D.
%
%   StatsTools https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Groppe DM, Urbach TP, Kutas M (2011a) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.
%       [3] Groppe DM, Urbach TP, Kutas M (2011b) Mass univariate analysis
%           of event-related brain potentials/fields II: Simulation studies
%           Psychophysiology, 48(12):1726-1737.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 18-Jan-2019

% Decode input variable arguments
[alpha,nperm,tail,rows,varx] = decode_varargin(varargin);

% Set up permutation test
mat = false;
if nargin<2 || isempty(x2)
    mat = true;
    warning('Comparing all columns of X1 using two-tailed test...')
    [x1,x2] = paircols(x1); tail = 'both';
elseif ~isnumeric(x2)
    error('X2 must be numeric or empty.')
elseif size(x1,2)~=size(x2,2)
    error('X1 and X2 must have the same number of variables.')
end

% Use only rows with no NaN values if specified
if strcmpi(rows,'complete')
    x1 = x1(~any(isnan(x1),2),:);
    x2 = x2(~any(isnan(x2),2),:);
end

% Concatenate data
x = [x1;x2];

% Get data dimensions, ignoring NaNs
[maxnobs,nvar] = size(x);
maxnobs1 = size(x1,1);
nobs = sum(~isnan(x));
nobs1 = sum(~isnan(x1));
nobs2 = sum(~isnan(x2));

% Compute degrees of freedom
df1 = nobs1-1; df2 = nobs2-1;
dfp = nobs./(nobs1.*nobs2);

% Compute t-statistic
sm1 = nansum(x1); sm2 = nansum(x2);
var1 = (nansum(x1.^2)-(sm1.^2)./nobs1)./df1;
var2 = (nansum(x2.^2)-(sm2.^2)./nobs2)./df2;
if strcmpi(varx,'equal')
    df = nobs-2;
    sd = sqrt((df1.*var1+df2.*var2)./df);
    se = sd.*sqrt(dfp);
elseif strcmpi(varx,'unequal')
    se2x1 = var1./nobs1;
    se2x2 = var2./nobs2;
    df = (se2x1+se2x2).^2./(se2x1.^2./df1+se2x2.^2./df2);
    sd = sqrt([var1;var2]);
    se = sqrt(se2x1+se2x2);
end
mx = sm1./nobs1-sm2./nobs2;
t = mx./se;

% Permute data and generate distribution of t-values
[~,idx] = sort(rand(maxnobs,nperm));
tp = zeros(nperm,nvar);
if any(isnan(x(:))) % for efficiency, only use NANSUM if necessary
    for i = 1:nperm
        xp1 = x(idx(1:maxnobs1,i),:); sm1 = nansum(xp1);
        xp2 = x(idx(maxnobs1+1:maxnobs,i),:); sm2 = nansum(xp2);
        s2x1 = (nansum(xp1.^2)-(sm1.^2)./nobs1)./df1;
        s2x2 = (nansum(xp2.^2)-(sm2.^2)./nobs2)./df2;
        if strcmpi(varx,'equal')
            se = sqrt((df1.*s2x1+df2.*s2x2)./df).*sqrt(dfp);
        elseif strcmpi(varx,'unequal')
            se = sqrt(s2x1./nobs1+s2x2./nobs2);
        end
        tp(i,:) = (sm1./nobs1-sm2./nobs2)./se;
    end
else
    for i = 1:nperm
        xp1 = x(idx(1:maxnobs1,i),:); sm1 = sum(xp1);
        xp2 = x(idx(maxnobs1+1:maxnobs,i),:); sm2 = sum(xp2);
        s2x1 = (sum(xp1.^2)-(sm1.^2)./nobs1)./df1;
        s2x2 = (sum(xp2.^2)-(sm2.^2)./nobs2)./df2;
        if strcmpi(varx,'equal')
            se = sqrt((df1.*s2x1+df2.*s2x2)./df).*sqrt(dfp);
        elseif strcmpi(varx,'unequal')
            se = sqrt(s2x1./nobs1+s2x2./nobs2);
        end
        tp(i,:) = (sm1./nobs1-sm2./nobs2)./se;
    end
end

% Compute Tmax with sign
[~,idx] = max(abs(tp),[],2);
csvar = [0;cumsum(ones(nperm-1,1)*nvar)];
tmax = tp'; tmax = tmax(idx+csvar);

% Compute adjusted test statistics using Tmax correction
if strcmpi(tail,'both')
    p = mean(abs(t)<tmax)*2;
    tcrit = prctile(tmax,100*(1-alpha/2));
    ci = [mx-tcrit.*se;mx+tcrit.*se];
    estal = mean(tcrit<tmax)*2;
elseif strcmpi(tail,'right')
    p = mean(t<tmax);
    tcrit = prctile(tmax,100*(1-alpha));
    ci = [mx+tcrit.*se;Inf(1,nvar)];
    estal = mean(tcrit<tmax);
elseif strcmpi(tail,'left')
    p = mean(t>tmax);
    tcrit = prctile(tmax,100*alpha);
    ci = [-Inf(1,nvar);mx+tcrit.*se];
    estal = mean(tcrit>tmax);
end

% Determine if adjusted p-values exceed desired alpha level
h = cast(p<alpha,'like',p);
h(isnan(tcrit)) = NaN;
p(isnan(tcrit)) = NaN;

% Arrange test results in a matrix if specified
if mat==true
    h = vec2mat(h);
    p = vec2mat(p);
    ciLwr = vec2mat(ci(1,:));
    ciUpr = vec2mat(ci(2,:));
    ci = cat(3,ciLwr,ciUpr);
    ci = permute(ci,[3,1,2]);
end

% Store values in structure
stats = struct('h',h,'p',p,'ci',ci,'tcrit',tcrit,'estal',estal);

% Execute if user requests unadjusted test statistics
if nargout > 2
    
    % Clear variables
    clear h p tcrit ciLwr ciUpr ci estal
    
    % Compute unadjusted test statistics
    if strcmpi(tail,'both')
        p = mean(abs(t)<tp)*2;
        tcrit = prctile(tp,100*(1-alpha/2));
        ci = [mx-tcrit.*se;mx+tcrit.*se];
        estal = mean(tcrit<tp)*2;
    elseif strcmpi(tail,'right')
        p = mean(t<tp);
        tcrit = prctile(tp,100*(1-alpha));
        ci = [mx+tcrit.*se;Inf(1,nvar)];
        estal = mean(tcrit<tp);
    elseif strcmpi(tail,'left')
        p = mean(t>tp);
        tcrit = prctile(tp,100*alpha);
        ci = [-Inf(1,nvar);mx+tcrit.*se];
        estal = mean(tcrit>tp);
    end
    
    % Determine if unadjusted p-values exceed desired alpha level
    h = cast(p<alpha,'like',p);
    h(isnan(tcrit(1,:))) = NaN;
    p(isnan(tcrit(1,:))) = NaN;
    
    % Arrange test results in a matrix if specified
    if mat==true
        h = vec2mat(h);
        p = vec2mat(p);
        tcrit = vec2mat(tcrit);
        ciLwr = vec2mat(ci(1,:));
        ciUpr = vec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
        estal = vec2mat(estal);
    end
    
    % Store values in structure
    orig = struct('h',h,'p',p,'ci',ci,'tcrit',tcrit,'estal',estal);
    
end

% Execute if user requests additional statistical parameters
if nargout > 3
    if mat==true
        df = vec2mat(df);
        sd = vec2mat(sd);
    end
    param = struct('df',df,'sd',sd);
end

% Arrange t-values in a matrix if specified
if mat==true
    t = vec2mat(t);
end

function [y1,y2] = paircols(x)
%paircols pair matrix columns and output as two separate matrices
%   [Y1,Y2] = PAIRCOLS(X) returns matrices Y1 and Y2 whose paired columns
%   correspond to every combination of column pairs in X. For efficiency,
%   repeated column pairs are skipped.

% Get matrix dimensions
[nobs,nvar] = size(x);

% Preallocate memory
y1 = zeros(nobs,(nvar^2-nvar)/2);
y2 = zeros(nobs,(nvar^2-nvar)/2);

% Initialize counters
ctr = 1;
jctr = 2;

% Generate paired matrices
for i = 1:nvar
    j = jctr;
    while j <= nvar
        y1(:,ctr) = x(:,i);
        y2(:,ctr) = x(:,j);
        j = j+1;
        ctr = ctr+1;
    end
    jctr = jctr+1;
end

function [y] = vec2mat(x)
%vec2mat convert vector output to matrix format
%   [Y] = VEC2MAT(X) returns a matrix Y by rearranging the values in vector
%   X according to their position as determined by PAIRCOLS. The values in
%   X may represent the output of some statistical test between every pair
%   of rows and columns in Y.

% Compute matrix dimensions
nvar = ceil(sqrt(length(x)*2));

% Preallocate memory
y = NaN(nvar,nvar);

% Initialize counters
ctr = 1;
jctr = 2;

% Generate matrix
for i = 1:nvar
    j = jctr;
    while j <= nvar
        y(i,j) = x(ctr);
        y(j,i) = x(ctr);
        j = j+1;
        ctr = ctr+1;
    end
    jctr = jctr+1;
end

function [alpha,nperm,tail,rows,varx] = decode_varargin(varargin)
%decode_varargin decode input variable arguments
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of the main function.

varargin = varargin{1,1};
if any(strcmpi(varargin,'alpha')) && ~isempty(varargin{find(strcmpi(varargin,'alpha'))+1})
    alpha = varargin{find(strcmpi(varargin,'alpha'))+1};
    if ~isscalar(alpha) || ~isnumeric(alpha) || isnan(alpha) || alpha<=0 || alpha>=1
        error('ALPHA must be a scalar between 0 and 1.')
    end
else
    alpha = 0.05;
end
if any(strcmpi(varargin,'nperm')) && ~isempty(varargin{find(strcmpi(varargin,'nperm'))+1})
    nperm = varargin{find(strcmpi(varargin,'nperm'))+1};
    if ~isscalar(nperm) || ~isnumeric(nperm) || isnan(nperm) || isinf(nperm) || floor(nperm)~=nperm || nperm<=0
        error('NPERM must be a positive integer.')
    elseif (nperm<1e3 && alpha<=0.05) || (nperm<5e3 && alpha<=0.01)
        warning('Number of permutations may be too low for chosen ALPHA.')
    end
else
    nperm = 1e4;
end
if any(strcmpi(varargin,'tail')) && ~isempty(varargin{find(strcmpi(varargin,'tail'))+1})
    tail = varargin{find(strcmpi(varargin,'tail'))+1};
    if ~any(strcmpi(tail,{'left','both','right'}))
        error('Invalid value for argument TAIL. Valid values are: ''left'', ''both'', ''right''.')
    end
else
    tail = 'both';
end
if any(strcmpi(varargin,'rows')) && ~isempty(varargin{find(strcmpi(varargin,'rows'))+1})
    rows = varargin{find(strcmpi(varargin,'rows'))+1};
    if ~any(strcmpi(rows,{'all','complete'}))
        error('Invalid value for argument ROWS. Valid values are: ''all'', ''complete''.')
    end
else
    rows = 'all';
end
if any(strcmpi(varargin,'varx')) && ~isempty(varargin{find(strcmpi(varargin,'varx'))+1})
    varx = varargin{find(strcmpi(varargin,'varx'))+1};
    if ~any(strcmpi(varx,{'equal','unequal'}))
        error('Invalid value for argument VARX. Valid values are: ''equal'', ''unequal''.')
    end
else
    varx = 'equal';
end