function [t,stats,orig,params] = permuttest2(x,y,varargin)
%PERMUTTEST2  Unpaired two-sample permutation test with tmax correction.
%   T = PERMUTTEST2(X,Y) returns the t-statistic of a two-sample
%   permutation test. If X and Y are matrices, multiple permutation tests
%   are performed simultaneously between each corresponding pair of columns
%   in X and Y and family-wise error rate (FWER) is controlled using the
%   tmax correction method (Blair et al., 1994). This method provides
%   strong control of FWER, even for small sample sizes, and is much more
%   powerful than traditional correction methods (Groppe et al., 2011a). It
%   is also rather insensitive to differences in population variance when
%   samples of equal size are used (Groppe et al., 2011b). For samples of
%   unequal size or variance, Welch's t-statistic may be used as it is less
%   sensitive to differences in variance (but also less sensitive to
%   differences in means). If Y is not entered, a permutation test between
%   each pair of columns in X is performed and output as a matrix. Samples
%   of different sizes may be used by replacing missing values with NaNs.
%   This function treats NaNs as missing values, and ignores them.
%
%   [T,STATS] = PERMUTTEST2(...) returns the adjusted test statistics in a
%   structure with the following fields:
%       'h'         -- test results. H=0 indicates the null hypothesis
%                      cannot be rejected, H=1 indicates the null
%                      hypothesis can be rejected.
%       'p'         -- the probability of observing the result by chance
%       'ci'    	-- 100*(1-(1/NPERM+ALPHA))% confidence interval for the
%                      true difference of population means (Groppe, 2016)
%       'tcrit'     -- critical t-value for the given alpha level. For two-
%                      tailed tests, the lower value is equal to -1*TCRIT.
%       'estal' 	-- the estimated alpha level of each test
%
%   [T,STATS,ORIG] = PERMUTTEST2(...) returns the original, unadjusted test
%   statistics in a structure with the same fields as STATS.
%
%   [T,STATS,ORIG,PARAMS] = PERMUTTEST2(...) returns other parameters in a
%   structure with the following fields:
%       'df'        -- the degrees of freedom of each test
%       'sd'    	-- the pooled estimate (equal variances) or unpooled
%                      estimates (unequal variances) of population standard
%                      deviation
%
%   [...] = PERMUTTEST2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'nperm'     A scalar specifying the number of permutations
%                   (default=10,000).
%       'tail'      A string specifying the alternative hypothesis:
%                       'both'      means are not equal (default)
%                       'right'     mean of X is greater than mean of Y
%                       'left'      mean of X is less than mean of Y
%       'rows'      A string specifying the rows to use in the case of any
%                   missing values (NaNs):
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'vartype'   A string specifying the variance equivalence of X and Y
%                   to determine the SD and t-statistic estimation method:
%                       'equal'   	assume equal variances (default)
%                       'unequal' 	assume unequal variances
%
%   Example 1: generate multivariate data for 2 samples, each with 20
%   variables and 30 observations and perform unpaired permutation tests
%   between the corresponding variables of each sample.
%       x = randn(30,20);
%       y = randn(30,20);
%       y(:,1:8) = y(:,1:8)-1;
%       [tstat,stats] = permuttest2(x,y)
%
%   Example 2: generate univariate data for 5 samples, each with 30
%   observations and perform unpaired permutation tests between every
%   sample (5 samples = 10 comparisons). Note that each column of X
%   represents an independent sample and may contain NaNs for samples with
%   smaller number of observations.
%       x = randn(30,5);
%       x(:,3:5) = x(:,3:5)-1;
%       [tstat,stats] = permuttest2(x)
%
%   See also PERMUTTEST PERMUVARTEST2 PERMUCORR DEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

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
%       [4] Groppe DM (2016) Combating the scientific decline effect with
%           confidence (intervals). Psychophysiology, 54(1):139-145.

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
if isempty(y)
    warning('Comparing all columns of X using two-tailed test...')
    [x,y] = paircols(x);
    arg.tail = 'both';
    arg.mat = true;
end
if size(x,2)~=size(y,2)
    error('X and Y must have the same number of variables.')
end

% Use only rows with no NaNs if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(x),2),:);
        y = y(~any(isnan(y),2),:);
end

% Get data dimensions, ignoring NaNs
maxnobsx = size(x,1);
nobsx = sum(~isnan(x));
nobsy = sum(~isnan(y));

% Compute degrees of freedom
dfx = nobsx-1;
dfy = nobsy-1;

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:))) || any(isnan(y(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Compute sample variance
smx = sum(x,nanflag);
smy = sum(y,nanflag);
varx = (sum(x.^2,nanflag)-(smx.^2)./nobsx)./dfx;
vary = (sum(y.^2,nanflag)-(smy.^2)./nobsy)./dfy;

% Concatenate samples
x = [x;y];
[maxnobs,nvar] = size(x);
nobs = sum(~isnan(x));
dfp = nobs./(nobsx.*nobsy);

% Compute t-statistic
switch arg.vartype
    case 'equal'
        df = nobs-2;
        sd = sqrt((dfx.*varx+dfy.*vary)./df);
        se = sd.*sqrt(dfp);
    case 'unequal'
        se2X = varx./nobsx;
        se2Y = vary./nobsy;
        df = (se2X+se2Y).^2./(se2X.^2./dfx+se2Y.^2./dfy);
        sd = sqrt([varx;vary]);
        se = sqrt(se2X+se2Y);
end
mx = smx./nobsx-smy./nobsy;
t = mx./se;

% Execute if user requests adjusted test statistics
if nargout > 1

    % Permute data and generate distribution of t-values
    [~,idx] = sort(rand(maxnobs,arg.nperm));
    tp = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        xp = x(idx(1:maxnobsx,i),:);
        yp = x(idx(maxnobsx+1:maxnobs,i),:);
        smx = sum(xp,nanflag);
        smy = sum(yp,nanflag);
        s2x = (sum(xp.^2,nanflag)-(smx.^2)./nobsx)./dfx;
        s2y = (sum(yp.^2,nanflag)-(smy.^2)./nobsy)./dfy;
        switch arg.vartype
            case 'equal'
                sep = sqrt((dfx.*s2x+dfy.*s2y)./df).*sqrt(dfp);
            case 'unequal'
                sep = sqrt(s2x./nobsx+s2y./nobsy);
        end
        tp(i,:) = (smx./nobsx-smy./nobsy)./sep;
    end

    % Compute tmax with sign
    [~,idx] = max(abs(tp),[],2);
    csvar = [0;cumsum(ones(arg.nperm-1,1)*nvar)];
    tmax = tp';
    tmax = tmax(idx+csvar);

    % Compute adjusted test statistics using tmax correction
    switch arg.tail
        case 'both'
            tmax = abs(tmax);
            p = (sum(abs(t)<=tmax)+1)/(arg.nperm+1);
            tcrit = prctile(tmax,100*(1-(1/arg.nperm+arg.alpha)));
            ci = [mx-tcrit.*se;mx+tcrit.*se];
            estal = mean(tmax>=tcrit);
        case 'right'
            p = (sum(t<=tmax)+1)/(arg.nperm+1);
            tcrit = prctile(tmax,100*(1-(1/arg.nperm+arg.alpha)));
            ci = [mx+tcrit.*se;Inf(1,nvar)];
            estal = mean(tmax>=tcrit);
        case 'left'
            p = (sum(t>=tmax)+1)/(arg.nperm+1);
            tcrit = prctile(tmax,100*(1/arg.nperm+arg.alpha));
            ci = [-Inf(1,nvar);mx+tcrit.*se];
            estal = mean(tmax<=tcrit);
    end

    % Determine if adjusted p-values exceed desired alpha level
    h = cast(p<arg.alpha,'like',p);
    h(isnan(tcrit)) = NaN;
    p(isnan(tcrit)) = NaN;

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

    % Compute unadjusted test statistics
    switch arg.tail
        case 'both'
            tp = abs(tp);
            p = (sum(abs(t)<=tp)+1)/(arg.nperm+1);
            tcrit = prctile(tp,100*(1-(1/arg.nperm+arg.alpha)));
            ci = [mx-tcrit.*se;mx+tcrit.*se];
            estal = mean(tp>=tcrit);
        case 'right'
            p = (sum(t<=tp)+1)/(arg.nperm+1);
            tcrit = prctile(tp,100*(1-(1/arg.nperm+arg.alpha)));
            ci = [mx+tcrit.*se;Inf(1,nvar)];
            estal = mean(tp>=tcrit);
        case 'left'
            p = (sum(t>=tp)+1)/(arg.nperm+1);
            tcrit = prctile(tp,100*(1/arg.nperm+arg.alpha));
            ci = [-Inf(1,nvar);mx+tcrit.*se];
            estal = mean(tp<=tcrit);
    end

    % Determine if unadjusted p-values exceed desired alpha level
    h = cast(p<arg.alpha,'like',p);
    h(isnan(tcrit(1,:))) = NaN;
    p(isnan(tcrit(1,:))) = NaN;

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

% XY variance equivalence
vartypeOptions = {'equal','unequal'};
validFcn = @(x) any(validatestring(x,vartypeOptions));
addParameter(p,'vartype','equal',validFcn);

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
arg.vartype = validatestring(arg.vartype,vartypeOptions);