function [t,stats,orig,params] = permuttest(x,y,varargin)
%PERMUTTEST  One/paired-sample permutation test with tmax correction.
%   T = PERMUTTEST(X,Y) returns the t-statistic of a paired-sample
%   permutation test. If X and Y are matrices, multiple permutation tests
%   are performed simultaneously between each corresponding pair of columns
%   in X and Y and family-wise error rate (FWER) is controlled using the
%   tmax correction method (Blair et al., 1994). This method provides
%   strong control of FWER, even for small sample sizes, and is much more
%   powerful than traditional correction methods (Gondan, 2010; Groppe et
%   al., 2011). For one-sample tests, enter the data column-wise in X and
%   leave Y empty. This function treats NaNs as missing values, and ignores
%   them.
%
%   [T,STATS] = PERMUTTEST(...) returns the adjusted test statistics in a
%   structure with the following fields:
%       'h'         -- test results. H=0 indicates the null hypothesis
%                      cannot be rejected, H=1 indicates the null
%                      hypothesis can be rejected.
%       'p'         -- the probability of observing the result by chance
%       'ci'    	-- 100*(1-ALPHA)% confidence interval for the true mean
%                      of X, or of X-Y for a paired test
%       'tcrit'     -- critical t-value for the given alpha level. For two-
%                      tailed tests, the lower value is equal to -1*TCRIT.
%       'estal' 	-- the estimated alpha level of each test
%
%   [T,STATS,ORIG] = PERMUTTEST(...) returns the original, unadjusted test
%   statistics in a structure with the same fields as STATS.
%
%   [T,STATS,ORIG,PARAMS] = PERMUTTEST(...) returns other parameters in a
%   structure with the following fields:
%       'df'        -- the degrees of freedom of each test
%       'sd'    	-- the estimated population standard deviation of X, or
%                      of X-Y for a paired test
%
%   [...] = PERMUTTEST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'nperm'     A scalar integer specifying the number of permutations
%                   (default=10,000, or all possible permutations for less
%                   than 14 observations).
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      mean is not M (two-tailed, default)
%                       'right'     mean is greater than M (right-tailed)
%                       'left'      mean is less than M (left-tailed)
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'sample'    A string specifying whether to perform a one-sample
%                   test or a paired-sample test when only X is entered:
%                       'one'       compare each column of X to zero and
%                                   store the results in a vector (default)
%                       'paired'    compare each pair of columns in X and
%                                   store the results in a matrix
%       'm'         A scalar or row vector specifying the mean of the null
%                   hypothesis for each variable (default=0).
%       'seed'      A scalar integer specifying the seed used to initialise
%                   the permutation generator. By default, the generator is
%                   initialised based on the current time, resulting in a
%                   different permutation each time.
%
%   Example 1: generate multivariate data for 2 conditions, each with 20
%   variables and 30 observations and perform paired permutation tests
%   between the corresponding variables of each condition.
%       rng(42);
%       x = randn(30,20);
%       y = randn(30,20);
%       y(:,1:8) = y(:,1:8)-1;
%       [t,stats] = permuttest(x,y) % paired-sample test
%       or
%       [t,stats] = permuttest(x) % one-sample test
%
%   Example 2: generate univariate data for 5 conditions, each with 30
%   observations and perform paired permutation tests between every pair of
%   conditions (5 conditions = 10 comparisons).
%       rng(42);
%       x = randn(30,5);
%       x(:,3:5) = x(:,3:5)-1;
%       [t,stats] = permuttest(x,[],'sample','paired') % paired-sample test
%
%   See also PERMUTTEST2 PERMUVARTEST2 PERMUCORR BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

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

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

if nargin<2
    y = [];
end

% Parse input arguments
arg = parsevarargin(varargin);

% Validate input parameters
validateparamin(x,y,arg)

% Orient data column-wise
if arg.dim==2 || isrow(x)
    x = x';
end
if ~isempty(y) && (arg.dim==2 || isrow(y))
    y = y';
end

% Set up permutation test
switch arg.sample
    case 'paired'
        if isempty(y)
            warning('Comparing all columns of X using a two-tailed test...')
            [x,y] = paircols(x);
            arg.tail = 'both';
            arg.mat = true;
        else
            error('The paired-sample test option only applies to X.')
        end
    otherwise
        if isempty(y)
            y = zeros(size(x));
        end
end
if size(x)~=size(y)
    error('X and Y must be the same size.')
end

% Compute difference between samples
x = x-y;

% Use only rows with no NaN values if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(x),2),:);
end

% Get data dimensions, ignoring NaNs
[maxnobs,nvar] = size(x);
nobs = sum(~isnan(x));

% Compute degrees of freedom
df = nobs-1;
dfp = sqrt(nobs.*df);

% Remove mean of null hypothesis from data
if isscalar(arg.m)
    x = x-arg.m;
else
    x = x-repmat(arg.m,maxnobs,1);
end

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute t-statistic
sd = std(x,nanflag);
mx = sum(x,nanflag)./nobs;
se = sd./sqrt(nobs);
t = mx./se;

% Execute if user requests adjusted test statistics
if nargout > 1

    % Use all possible permutations if less than 14 observations
    if min(nobs) < 14
        warning('Computing all possible permutations due to small N.')
        arg.nperm = 2^min(nobs);
    end

    % Permute data and generate distribution of t-values
    rng(arg.seed);
    signx = sign(rand(maxnobs,arg.nperm)-0.5);
    tp = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        xp = x.*repmat(signx(:,i),1,nvar);
        sm = sum(xp,nanflag);
        tp(i,:) = sm./nobs./(sqrt(sum(xp.^2,nanflag)-(sm.^2)./nobs)./dfp);
    end

    % Compute tmax without sign and add negative values
    tmax = max(abs(tp),[],2);
    tmax(arg.nperm+1:2*arg.nperm) = -tmax;

    % Compute adjusted test statistics using tmax correction
    switch arg.tail
        case 'both'
            p = 2*(sum(abs(t)<=tmax)+1)/(arg.nperm+1);
            tcrit = prctile(tmax,100*(1-arg.alpha/2));
            ci = [mx-tcrit.*se;mx+tcrit.*se];
            estal = mean(tcrit<tmax)*2;
        case 'right'
            p = (sum(t<=tmax)+1)/(arg.nperm+1);
            tcrit = prctile(tmax,100*(1-arg.alpha));
            ci = [mx+tcrit.*se;Inf(1,nvar)];
            estal = mean(tcrit<tmax);
        case 'left'
            p = (sum(t>=tmax)+1)/(arg.nperm+1);
            tcrit = prctile(tmax,100*arg.alpha);
            ci = [-Inf(1,nvar);mx+tcrit.*se];
            estal = mean(tcrit>tmax);
    end

    % Determine if adjusted p-values exceed desired alpha level
    h = cast(p<arg.alpha,'like',p);
    h(isnan(tcrit)) = NaN;
    p(isnan(tcrit)) = NaN;
    p(p>1) = 1;

    % Arrange test results in a matrix if specified
    if arg.mat
        h = vec2mat(h);
        p = vec2mat(p);
        ciLwr = vec2mat(ci(1,:));
        ciUpr = vec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
    end

    % Store values in a structure
    stats = struct('h',h,'p',p,'ci',ci,'tcrit',tcrit,'estal',estal);

end

% Execute if user requests unadjusted test statistics
if nargout > 2

    clear h p ci tcrit estal ciLwr ciUpr

    % Add negative values
    tp(arg.nperm+1:2*arg.nperm,:) = -tp;

    % Compute unadjusted test statistics
    switch arg.tail
        case 'both'
            p = 2*(sum(abs(t)<=tp)+1)/(arg.nperm+1);
            tcrit = prctile(tp,100*(1-arg.alpha/2));
            ci = [mx-tcrit.*se;mx+tcrit.*se];
            estal = mean(tcrit<tp)*2;
        case 'right'
            p = (sum(t<=tp)+1)/(arg.nperm+1);
            tcrit = prctile(tp,100*(1-arg.alpha));
            ci = [mx+tcrit.*se;Inf(1,nvar)];
            estal = mean(tcrit<tp);
        case 'left'
            p = (sum(t>=tp)+1)/(arg.nperm+1);
            tcrit = prctile(tp,100*arg.alpha);
            ci = [-Inf(1,nvar);mx+tcrit.*se];
            estal = mean(tcrit>tp);
    end

    % Determine if unadjusted p-values exceed desired alpha level
    h = cast(p<arg.alpha,'like',p);
    h(isnan(tcrit(1,:))) = NaN;
    p(isnan(tcrit(1,:))) = NaN;
    p(p>1) = 1;

    % Arrange test results in a matrix if specified
    if arg.mat
        h = vec2mat(h);
        p = vec2mat(p);
        tcrit = vec2mat(tcrit);
        ciLwr = vec2mat(ci(1,:));
        ciUpr = vec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
        estal = vec2mat(estal);
    end

    % Store values in a structure
    orig = struct('h',h,'p',p,'ci',ci,'tcrit',tcrit,'estal',estal);

end

% Execute if user requests parameters
if nargout > 3
    if arg.mat
        df = vec2mat(df);
        sd = vec2mat(sd);
    end
    params = struct('df',df,'sd',sd);
end

% Arrange t-values in a matrix if specified
if arg.mat
    t = vec2mat(t);
end

function [y1,y2] = paircols(x)
%PAIRCOLS  Pair matrix columns and output as two separate matrices.
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
%VEC2MAT  Convert vector output to matrix format.
%   Y = VEC2MAT(X) returns a matrix Y by rearranging the values in vector
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

function validateparamin(x,y,arg)
%VALIDATEPARAMIN  Validate input parameters.
%   VALIDATEPARAMIN(X,Y,ARG) validates the input parameters of the main
%   function.

if ~isnumeric(x)
    error('X must be numeric.')
elseif ~isnumeric(y)
    error('Y must be numeric or empty.')
end
if (arg.nperm<1e3 && arg.alpha<=0.05) || (arg.nperm<5e3 && arg.alpha<=0.01)
    warning('Number of permutations may be too low for chosen ALPHA.')
end

function arg = parsevarargin(varargin)
%PARSEVARARGIN  Parse input arguments.
%   [PARAM1,PARAM2,...] = PARSEVARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   parses the input arguments of the main function.

% Create parser object
p = inputParser;

% Dimension to work along
errorMsg = 'It must be a positive integer scalar within indexing range.';
validFcn = @(x) assert(x==1||x==2,errorMsg);
addParameter(p,'dim',1,validFcn);

% Alpha level
errorMsg = 'It must be a scalar between 0 and 1.';
validFcn = @(x) assert(x>0&&x<1,errorMsg);
addParameter(p,'alpha',0.05,validFcn);

% Number of permutations
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x)&&x>0,errorMsg);
addParameter(p,'nperm',1e4,validFcn);

% Alternative hypothesis
tailOptions = {'left','both','right'};
validFcn = @(x) any(validatestring(x,tailOptions));
addParameter(p,'tail','both',validFcn);

% Rows to use if NaNs
rowsOptions = {'all','complete'};
validFcn = @(x) any(validatestring(x,rowsOptions));
addParameter(p,'rows','all',validFcn);

% Sample type
sampleOptions = {'one','paired'};
validFcn = @(x) any(validatestring(x,sampleOptions));
addParameter(p,'sample','one',validFcn);

% Null hypothesis mean
errorMsg = 'It must be a scalar or row vector.';
validFcn = @(x) assert(isnumeric(x),errorMsg);
addParameter(p,'m',0,validFcn);

% Permutation generator seed
errorMsg = 'It must be an integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'seed','shuffle',validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'mat',false,validFcn); % matrix flag

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially matched strings
arg.tail = validatestring(arg.tail,tailOptions);
arg.rows = validatestring(arg.rows,rowsOptions);
arg.sample = validatestring(arg.sample,sampleOptions);