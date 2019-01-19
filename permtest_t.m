function [t,stats,orig,param] = permtest_t(x1,x2,varargin)
%permtest_t one-sample/paired-sample permutation test with Tmax correction
%   T = PERMTEST_T(X1,X2) returns the t-statistic of a paired-sample
%   permutation test. If X1 and X2 are matrices, multiple permutation tests
%   are performed simultaneously between each corresponding pair of columns
%   in X1 and X2 and family-wise error rate (FWER) is controlled using the
%   Tmax correction method (Blair et al., 1994). This method provides
%   strong control of FWER, even for small sample sizes, and is much more
%   powerful than traditional correction methods (Gondan, 2010; Groppe et
%   al., 2011). For one-sample tests, enter the data column-wise in X1 and
%   leave X2 empty. This function treats NaNs as missing values, and
%   ignores them.
%
%   [...,STATS] = permtest_t(...) returns the adjusted test statistics in
%   a structure containing the following fields:
%       h       test results: H=0 indicates the null hypothesis cannot be
%               rejected, H=1 indicates the null hypothesis can be rejected
%       p       the probability of observing the given result by chance
%       ci    	100*(1-ALPHA)% confidence interval for the true mean of X1,
%               or of X1-X2 for a paired test
%       tcrit 	critical t-value for the given alpha level (for two-tailed
%               tests, the lower critical t-value is equal to -1*TCRIT)
%       estal 	the estimated alpha level of each test
%
%   [...,ORIG] = permtest_t(...) returns the original, unadjusted test
%   statistics in a structure containing the same fields as STATS.
%
%   [...,PARAM] = permtest_t(...) returns some statistical parameters in a
%   structure containing the following fields:
%       df      the degrees of freedom of each test
%       sd    	the estimated population standard deviation of X1, or of
%               X1-X2 for a paired test
%
%   [...] = permtest_t(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%   Parameter   Value
%   'm'         a scalar or row vector specifying the mean of the null
%               hypothesis for each variable (default=0)
%   'alpha'     a scalar between 0 and 1 specifying the significance level
%               as 100*ALPHA% (default=0.05)
%   'nperm'     a scalar specifying the number of permutations (default=
%               10,000 or all possible permutations for less than 14 obs.)
%   'tail'      a string specifying the alternative hypothesis
%                   'both'      mean is not M (two-tailed test, default)
%                   'right'     mean is greater than M (right-tailed test)
%                   'left'      mean is less than M (left-tailed test)
%   'rows'      a string specifying the rows to use in the case of any
%               missing values (NaNs)
%                   'all'       use all rows, and ignore any NaNs (default)
%                   'complete'  use only rows with no NaNs
%   'sample'    a string specifying whether to perform a one-sample test or
%               a paired-sample test when only X1 is entered
%                   'one'       compare each column of X1 to zero and
%                               output the results as a vector (default)
%                   'paired'    compare each pair of columns in X1 to each
%                               other and output the results as a matrix
%
%   Example 1: generate multivariate data for 2 conditions, each with 20
%   variables and 30 observations and perform paired permutation tests
%   between the corresponding variables of each condition.
%       x1 = randn(30,20);
%       x2 = randn(30,20);
%       x2(:,1:8) = x2(:,1:8)-1;
%       [t,stats] = permtest_t(x1,x2) % paired-sample test
%       or
%       [t,stats] = permtest_t(x1-x2) % one-sample test
%
%   Example 2: generate univariate data for 5 conditions, each with 30
%   observations and perform paired permutation tests between every pair of
%   conditions (5 conditions = 10 comparisons).
%       x = randn(30,5);
%       x(:,3:5) = x(:,3:5)-1;
%       [t,stats] = permtest_t(x,[],'sample','paired') % paired-sample test
%
%   See also PERMTEST_T2 PERMTEST_F PERMTEST_CORR EFFECTSIZE_D.
%
%   StatsTools https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Gondan M (2010) A permutation test for the race model
%           inequality. Behav Res Methods, 42(1):23-28.
%       [3] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 18-Jan-2019

% Decode input variable arguments
[alpha,nperm,tail,rows,sample,m] = decode_varargin(varargin);

% Set up permutation test
if nargin<2 || (isempty(x2) && strcmpi(sample,'one'))
    x2 = 0;
elseif isempty(x2) && strcmpi(sample,'paired')
    warning('Comparing all columns of X1 using two-tailed test...')
    [x1,x2] = paircols(x1); tail = 'both';
elseif ~isempty(x2) && strcmpi(sample,'paired')
    error('PAIRED option only applies to X1.')
elseif ~isnumeric(x2)
    error('X2 must be numeric or empty.')
elseif size(x1)~=size(x2)
    error('The data in a paired test must be the same size.')
end

% Compute difference between samples
x = x1-x2;

% Use only rows with no NaN values if specified
if strcmpi(rows,'complete')
    x = x(~any(isnan(x),2),:);
end

% Get data dimensions, ignoring NaNs
[maxnobs,nvar] = size(x);
nobs = sum(~isnan(x));

% Use all possible permutations if less than 14 observations
if min(nobs)<14
    nperm = 2^min(nobs);
end

% Compute degrees of freedom
df = nobs-1;
dfp = sqrt(nobs.*df);

% Remove mean of null hypothesis from data
if isscalar(m)
    x = x-m;
else
    x = x-repmat(m,maxnobs,1);
end

% Compute t-statistic
sd = nanstd(x);
mx = nansum(x)./nobs;
se = sd./sqrt(nobs);
t = mx./se;

% Permute data and generate distribution of t-values
signx = sign(rand(maxnobs,nperm)-0.5);
tp = zeros(nperm,nvar);
if any(isnan(x(:))) % for efficiency, only use NANSUM if necessary
    for i = 1:nperm
        xp = x.*repmat(signx(:,i),1,nvar); sm = nansum(xp);
        tp(i,:) = sm./nobs./(sqrt(sum(xp.^2)-(sm.^2)./nobs)./dfp);
    end
else
    for i = 1:nperm
        xp = x.*repmat(signx(:,i),1,nvar); sm = sum(xp);
        tp(i,:) = sm./nobs./(sqrt(sum(xp.^2)-(sm.^2)./nobs)./dfp);
    end
end

% Compute Tmax without sign and add negative values
tmax = max(abs(tp),[],2);
tmax(nperm+1:2*nperm) = -tmax;

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
if strcmpi(sample,'paired')
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
    clear h p ci tcrit estal ciLwr ciUpr
    
    % Add negative values
    tp(nperm+1:2*nperm,:) = -tp;
    
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
    if strcmpi(sample,'paired')
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
    if strcmpi(sample,'paired')
        df = vec2mat(df);
        sd = vec2mat(sd);
    end
    param = struct('df',df,'sd',sd);
end

% Arrange t-values in a matrix if specified
if strcmpi(sample,'paired')
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

function [alpha,nperm,tail,rows,sample,m] = decode_varargin(varargin)
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
if any(strcmpi(varargin,'sample')) && ~isempty(varargin{find(strcmpi(varargin,'sample'))+1})
    sample = varargin{find(strcmpi(varargin,'sample'))+1};
    if ~any(strcmpi(sample,{'one','paired'}))
        error('Invalid value for argument SAMPLE. Valid values are: ''one'', ''paired''.')
    end
else
    sample = 'one';
end
if any(strcmpi(varargin,'m')) && ~isempty(varargin{find(strcmpi(varargin,'m'))+1})
    m = varargin{find(strcmpi(varargin,'m'))+1};
    if ~isnumeric(m) || any(isnan(m)) || any(isinf(m))
        error('M must be a scalar or vector of numeric values.')
    end
else
    m = 0;
end