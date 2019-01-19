function [f,stats,orig,param] = permtest_f2(x1,x2,varargin)
%permtest_f2 permutation-based F-test with max statistic correction
%   F = PERMTEST_F2(X1,X2) returns the F-statistic of a two-sample
%   permutation test. If X1 and X2 are matrices, multiple permutation tests
%   are performed simultaneously between each corresponding pair of columns
%   in X1 and X2 and family-wise error rate (FWER) is controlled using the 
%   Tmax correction method (Blair et al., 1994). This method provides
%   strong control of FWER, even for small sample sizes, and is much more
%   powerful than traditional correction methods (Groppe et al., 2011).
%   If X2 is not entered, a permutation test between each pair of columns 
%   in X1 is performed and output as a matrix. Samples of different sizes 
%   may be used by replacing missing values with NaNs. This function treats 
%   NaNs as missing values, and ignores them.
%
%   [...,STATS] = permtest_f2(...) returns the adjusted test statistics in
%   a structure containing the following fields:
%       h       test results: H=0 indicates the null hypothesis cannot be
%               rejected, H=1 indicates the null hypothesis can be rejected
%       p       the probability of observing the given result by chance
%       ci    	100*(1-ALPHA)% confidence interval for the true difference
%               of sample variances
%       fcrit 	critical F-value for the given alpha level (for two-tailed
%               tests, the lower t-value is equal to -1*fcrit)
%       estal 	the estimated alpha level of each test
%
%   [...,ORIG] = permtest_f2(...) returns the original, unadjusted test
%   statistics in a structure containing the same fields as STATS.
%
%   [...,PARAM] = permtest_f2(...) returns some statistical parameters in a
%   structure containing the following fields:
%       df1     the numerator degrees of freedom of each test
%       df1     the denominator degrees of freedom of each test
%
%   [...] = permtest_f2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%
%   Example 1: generate multivariate data for 2 samples, each with 20
%   variables and 30 observations and perform unpaired permutation tests
%   between the corresponding variables of each sample.
%       x1 = randn(30,20);
%       x2 = randn(30,20);
%       x2(:,1:8) = x2(:,1:8)-1;
%       [t,stats] = permtest_f2(x1,x2)
%
%   Example 2: generate univariate data for 5 samples, each with 30
%   observations and perform unpaired permutation tests between every
%   sample (5 samples = 10 comparisons). Note that each column of X
%   represents an independent sample and may contain NaNs for samples with
%   smaller number of observations.
%       x = randn(30,5);
%       x(:,3:5) = x(:,3:5)-1;
%       [t,stats] = permtest_f2(x)
%
%   See also PERMTEST_T PERMTEST_T2 PERMTEST_CORR EFFECTSIZE_D.
%
%   StatsTools https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 18-Jan-2019

% Decode input variable arguments
[alpha,nperm,tail,rows] = decode_varargin(varargin);

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
nobs1 = sum(~isnan(x1));
nobs2 = sum(~isnan(x2));

% Compute degrees of freedom
df1 = nobs1-1; df2 = nobs2-1;

% Compute F-statistic
var1 = (nansum(x1.^2)-(nansum(x1).^2)./nobs1)./df1;
var2 = (nansum(x2.^2)-(nansum(x2).^2)./nobs2)./df2;
f = var1./var2;

% Permute data and generate distribution of F-values
[~,idx] = sort(rand(maxnobs,nperm));
fp = zeros(nperm,nvar);
if any(isnan(x(:))) % for efficiency, only use NANSUM if necessary
    for i = 1:nperm
        xp1 = x(idx(1:maxnobs1,i),:);
        xp2 = x(idx(maxnobs1+1:maxnobs,i),:);
        var1 = (nansum(xp1.^2)-(nansum(xp1).^2)./nobs1)./df1;
        var2 = (nansum(xp2.^2)-(nansum(xp2).^2)./nobs2)./df2;
        fp(i,:) = var1./var2;
    end
else
    for i = 1:nperm
        xp1 = x(idx(1:maxnobs1,i),:);
        xp2 = x(idx(maxnobs1+1:maxnobs,i),:);
        var1 = (sum(xp1.^2)-(sum(xp1).^2)./nobs1)./df1;
        var2 = (sum(xp2.^2)-(sum(xp2).^2)./nobs2)./df2;
        fp(i,:) = var1./var2;
    end
end

% Compute Fmax with sign
[~,idx] = max(abs(fp),[],2);
csvar = [0;cumsum(ones(nperm-1,1)*nvar)];
fmax = fp'; fmax = fmax(idx+csvar);

% Compute corrected test statistics using max statistic correction
if strcmpi(tail,'both')
    p = mean(abs(f)<fmax)*2;
    fcrit = prctile(fmax,100*(1-alpha/2));
    ci = [f.*finv(alpha/2,df2,df1);f./finv(alpha/2,df1,df2)];
    estal = mean(fcrit<fmax)*2;
elseif strcmpi(tail,'right')
    p = mean(f<fmax);
    fcrit = prctile(fmax,100*(1-alpha));
    ci = [f.*finv(alpha/2,df2,df1);Inf(1,nvar)];
    estal = mean(fcrit<fmax);
elseif strcmpi(tail,'left')
    p = mean(f>fmax);
    fcrit = prctile(fmax,100*alpha);
    ci = [zeros(1,nvar);f./finv(alpha/2,df1,df2)];
    estal = mean(fcrit>fmax);
end

% Determine if adjusted p-values exceed desired alpha level
h = cast(p<alpha,'like',p);
h(isnan(fcrit)) = NaN;
p(isnan(fcrit)) = NaN;

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
stats = struct('h',h,'p',p,'ci',ci,'fcrit',fcrit,'estal',estal);

% Execute if user requests unadjusted test statistics
if nargout > 2
    
    % Clear variables
    clear h p fcrit ciLwr ciUpr ci estal
    
    % Compute unadjusted test statistics
    if strcmpi(tail,'both')
        p = mean(abs(f)<fp)*2;
        fcrit = prctile(fp,100*(1-alpha/2));
        ci = [f.*finv(alpha/2,df2,df1);f./finv(alpha/2,df1,df2)];
        estal = mean(fcrit<fp)*2;
    elseif strcmpi(tail,'right')
        p = mean(f<fp);
        fcrit = prctile(fp,100*(1-alpha));
        ci = [f.*finv(alpha/2,df2,df1);Inf(1,nvar)];
        estal = mean(fcrit<fp);
    elseif strcmpi(tail,'left')
        p = mean(f>fp);
        fcrit = prctile(fp,100*alpha);
        ci = [zeros(1,nvar);f./finv(alpha/2,df1,df2)];
        estal = mean(fcrit>fp);
    end
    
    % Determine if unadjusted p-values exceed desired alpha level
    h = cast(p<alpha,'like',p);
    h(isnan(fcrit(1,:))) = NaN;
    p(isnan(fcrit(1,:))) = NaN;
    
    % Arrange test results in a matrix if specified
    if mat==true
        h = vec2mat(h);
        p = vec2mat(p);
        fcrit = vec2mat(fcrit);
        ciLwr = vec2mat(ci(1,:));
        ciUpr = vec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
        estal = vec2mat(estal);
    end
    
    % Store values in structure
    orig = struct('h',h,'p',p,'ci',ci,'fcrit',fcrit,'estal',estal);
    
end

% Execute if user requests additional statistical parameters
if nargout > 3
    if mat==true
        df1 = vec2mat(df1);
        df2 = vec2mat(df2);
    end
    param = struct('df1',df1,'df2',df2);
end

% Arrange F-values in a matrix if specified
if mat==true
    f = vec2mat(f);
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

function [alpha,nperm,tail,rows] = decode_varargin(varargin)
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