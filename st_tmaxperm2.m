function [tstat,corx,orig,stats] = st_tmaxperm2(x1,x2,varargin)
%st_tmaxperm2 independent/two-sample permutation test with Tmax correction
%   TSTAT = ST_TMAXPERM2(X1,X2) returns the t-statistic of a two-sample
%   permutation test. If X1 and X2 are matrices, multiple permutation tests
%   are performed simultaneously between each corresponding pair of columns
%   in X1 and X2 and family-wise error rate is controlled using the Tmax
%   correction method (Blair et al., 1994). See Groppe et al. (2011) for
%   comparisons with alternative correction methods and for simulations of
%   its sensitivity to differences in population variance. This function
%   treats NaNs as missing values, and ignores them.
%
%   This function performs an unpaired two-sample permutation test. For a
%   paired test, use the ST_TMAXPERM function.
%
%   [...,CORX] = ST_TMAXPERM2(...) returns the corrected test statistics in
%   a structure containing the following fields:
%       h       test results: H=0 indicates the null hypothesis cannot be
%               rejected, H=1 indicates the null hypothesis can be rejected
%       p       the probability of observing the given result by chance
%       tcrit 	critical t-value for the given alpha level (for two-tailed
%               tests, the lower t-value is equal to -1*TCRIT)
%       ci    	100*(1-ALPHA)% confidence interval for the true difference
%               of population means
%       estal 	the estimated alpha level of each test
%
%   [...,ORIG] = ST_TMAXPERM2(...) returns the original, uncorrected test
%   statistics in a structure containing the same fields as CORX.
%
%   [...,STATS] = ST_TMAXPERM2(...) returns some general data statistics in
%   a structure containing the following fields:
%       df      the degrees of freedom of each test
%       sd    	pooled estimate (equal variances) or unpooled estimates
%               (unequal variances) of population standard deviation
%
%   [...] = ST_TMAXPERM2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
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
%   'rows'      a string specifying the rows to use in the permutation test
%                   'all'       use all rows, and ignore any NaN values
%                   'complete'  use only rows with no NaN values
%   'sample'    a string specifying whether to perform a two-sample test
%               between X1 and X2 or an unpaired test between every pair of
%               columns in X1 when X2 is not entered
%                   'two'       compare each column of X1 to each column of
%                               X2 and output results as a vector (default)
%                   'unpaired'  compare each column of X1 to every other
%                               column of X1 and output results as a matrix
%   'varx'      a string specifying the variance equivalence of X1 and X2
%                   'equal'   	assume equal variances (default)
%                   'unequal' 	assume unequal variances
%
%   Example 1: generate multivariate data for 2 populations, each with 20
%   variables and 30 observations and perform unpaired permutation tests
%   between the corresponding variables of each population.
%       x1 = randn(30,20);
%       x2 = randn(30,20);
%       x2(:,1:8) = x2(:,1:8)-1;
%       [tstat,corx,orig] = st_tmaxperm2(x1,x2)
%
%   Example 2: generate univariate data for 5 populations, each with 1
%   variable and 30 observations and perform unpaired permutation tests
%   between every population (5 populations = 10 comparisons). Note that
%   each column of X represents an independent sample.
%       x = randn(30,5);
%       x(:,3:5) = x(:,3:5)-1;
%       [tstat,corx,orig] = st_tmaxperm2(x,[],'sample','unpaired')
%
%   See also ST_TMAXPERM ST_FMAXPERM2 ST_RMAXPERM.
%
%   StatsTools https://github.com/mickcrosse/StatsTools

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields II: Simulation studies
%           Psychophysiology, 48(12):1726-1737.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 26-Sep-2018

% Decode input variable arguments
[alpha,nperm,tail,rows,sample,varx] = st_decode_varargin(varargin);

% Check number of variables is same
if nargin<2 || (isempty(x2) && strcmpi(sample,'one'))
    error('Requires at least two input arguments.')
elseif isempty(x2) && strcmpi(sample,'unpaired')
    warning('Comparing all columns of X1 using two-tailed test...')
    [x1,x2] = st_paircols(x1);
    tail = 'both';
elseif ~isempty(x2) && strcmpi(sample,'unpaired')
    error('UNPAIRED option only applies to X1.')
elseif size(x1,2)~=size(x2,2)
    error('X1 and X2 must have the same number of variables.')
end

% Use only rows with no NaN values if specified
if strcmpi(rows,'complete')
    x1 = x1(any(~isnan(x1),2),:);
    x2 = x2(any(~isnan(x2),2),:);
end

% Concatenate data
x = [x1;x2];

% Get data dimensions, ignoring NaNs
nvar = size(x,2);
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
tstat = mx./se;

% Permute data and generate distribution of t-scores
[~,idx] = sort(rand(max(nobs),nperm));
tp = zeros(nperm,nvar);
if any(isnan(x(:))) % for efficiency, only use NANSUM if necessary
    for i = 1:nperm
        xp1 = x(idx(1:nobs1,i),:); sm1 = nansum(xp1);
        xp2 = x(idx(nobs1+1:nobs,i),:); sm2 = nansum(xp2);
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
        xp1 = x(idx(1:nobs1,i),:); sm1 = sum(xp1);
        xp2 = x(idx(nobs1+1:nobs,i),:); sm2 = sum(xp2);
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

% Compute corrected test statistics using Tmax correction
if strcmpi(tail,'both')
    p = mean(abs(tstat)<tmax)*2;
    tcrit = prctile(tmax,100*(1-alpha/2));
    ci = [mx-tcrit.*se;mx+tcrit.*se];
    estal = mean(tcrit<tmax)*2;
elseif strcmpi(tail,'right')
    p = mean(tstat<tmax);
    tcrit = prctile(tmax,100*(1-alpha));
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

% Arrange test results in a matrix if specified
if strcmpi(sample,'unpaired')
    h = st_ttestmat(h);
    p = st_ttestmat(p);
    ciLwr = st_ttestmat(ci(1,:));
    ciUpr = st_ttestmat(ci(2,:));
    ci = cat(3,ciLwr,ciUpr);
    ci = permute(ci,[3,1,2]);
end

% Store values in structure
corx = struct('h',h,'p',p,'tcrit',tcrit,'ci',ci,'estal',estal);

% Execute if user specifies uncorrected test statistics
if nargout > 2
    
    % Clear variables
    clear h p tcrit ciLwr ciUpr ci estal
    
    % Compute original test statistics without correction
    if strcmpi(tail,'both')
        p = mean(abs(tstat)<tp)*2;
        tcrit = prctile(tp,100*(1-alpha/2));
        ci = [mx-tcrit.*se;mx+tcrit.*se];
        estal = mean(tcrit<tp)*2;
    elseif strcmpi(tail,'right')
        p = mean(tstat<tp);
        tcrit = prctile(tp,100*(1-alpha));
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
    
    % Arrange test results in a matrix if specified
    if strcmpi(sample,'unpaired')
        h = st_ttestmat(h);
        p = st_ttestmat(p);
        tcrit = st_ttestmat(tcrit);
        ciLwr = st_ttestmat(ci(1,:));
        ciUpr = st_ttestmat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
        estal = st_ttestmat(estal);
    end
    
    % Store values in structure
    orig = struct('h',h,'p',p,'tcrit',tcrit,'ci',ci,'estal',estal);
    
end

% Execute if user specifies general data statistics
if nargout > 3
    if strcmpi(sample,'unpaired')
        df = st_ttestmat(df);
        sd = st_ttestmat(sd);
    end
    stats = struct('df',df,'sd',sd);
end

% Arrange t-values in a matrix if specified
if strcmpi(sample,'unpaired')
    tstat = st_ttestmat(tstat);
end

function [alpha,nperm,tail,rows,sample,varx] = st_decode_varargin(varargin)
% st_decode_varargin decode input variable arguments
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
    if ~isscalar(nperm) || ~isnumeric(nperm) || isnan(nperm) || isinf(nperm) || floor(nperm)~=nperm
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
    if ~any(strcmpi(sample,{'two','unpaired'}))
        error('Invalid value for argument SAMPLE. Valid values are: ''two'', ''unpaired''.')
    end
else
    sample = 'two';
end
if any(strcmpi(varargin,'varx')) && ~isempty(varargin{find(strcmpi(varargin,'varx'))+1})
    varx = varargin{find(strcmpi(varargin,'varx'))+1};
    if ~any(strcmpi(varx,{'equal','unequal'}))
        error('Invalid value for argument VARX. Valid values are: ''equal'', ''unequal''.')
    end
else
    varx = 'equal';
end

function [y1,y2] = st_paircols(x)
%st_paircols pair matrix columns and output as two separate matrices

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

function [y] = st_ttestmat(x)
%st_ttestmat generate a t-test matrix

% Compute matrix dimensions
nvar = ceil(sqrt(length(x)*2));

% Preallocate memory
y = NaN(nvar,nvar);

% Initialize counters
ctr = 1;
jctr = 2;

% Generate t-test matrix
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