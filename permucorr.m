function [r,stats,orig] = permucorr(x,y,varargin)
%PERMUCORR  Permutation-based correlation with max statistic correction.
%   R = PERMUCORR(X,Y) returns the pairwise linear correlation
%   coefficient of a permutation test based on Pearson's r. If X and Y are
%   matrices, multiple permutation tests are performed simultaneously
%   between each corresponding pair of columns in X and Y and family-wise
%   error rate (FWER) is controlled using the max statistic correction
%   method (Blair et al., 1994). This method provides strong control of
%   FWER, even for small sample sizes, and is much more powerful than
%   traditional correction methods (Groppe et al., 2011). For nonlinear
%   correlations, the raw data may be transformed to rank orders using a
%   Spearman's or a rankit transformation (Bliss, 1967; Bishara & Hittner,
%   2012). If Y is not entered, the correlation between each pair of
%   columns in X is computed and output as a correlation matrix.
%
%   [R,STATS] = PERMUCORR(...) returns the adjusted test statistics in a
%   structure with the following fields:
%       'p'         -- the probability of observing the result by chance
%       'ci'    	-- 100*(1-ALPHA)% confidence interval. For Pearson and
%                      rankit correlations, the Fisher z' method is used
%                      (Fisher, 1958) and for Spearman rank correlations,
%                      the Fieller et al. (1957) estimate is used (Bishara
%                      & Hittner, 2017)
%       'rcrit'     -- critical r-value for the given alpha level (for two-
%                      tailed tests, the lower t-value is equal to
%                      -1*RCRIT)
%       'estal' 	-- the estimated alpha level of each test
%
%   [R,STATS,ORIG] = PERMUCORR(...) returns the original, unadjusted test
%   statistics in a structure with the same fields as STATS.
%
%   [...] = PERMUCORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'alpha'     a scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05)
%       'nperm'     a scalar specifying the number of permutations
%                   (default=10,000 or all possible permutations for less
%                   than 14 observations)
%       'tail'      string specifying the alternative hypothesis
%                       'both'      correlation is not zero (default)
%                       'right'     correlation is greater than zero
%                       'left'      correlation is less than zero
%       'rows'      a string specifying the rows to use in the case of any
%                   missing values (NaNs)
%                       'all'       use all rows, even with NaNs (default)
%                       'complete'  use only rows with no NaNs
%       'type'      string specifying the type of correlation measure
%                       'Pearson'   Pearson's correlation coefficient (def)
%                       'Spearman'  Spearman's rank correlation coefficient
%                       'Rankit'    Bliss's rankit correlation coefficient
%
%   Example 1: generate multivariate data for 2 conditions, each with 20
%   variables and 30 observations and calculate the correlation between the
%   corresponding variables of each condition.
%       x = randn(30,20);
%       y = randn(30,20);
%       y(:,1:5) = y(:,1:5)+0.5*x(:,1:5);
%       y(:,15:20) = y(:,15:20)-x(:,15:20);
%       [r,stats] = permucorr(x,y)
%
%   Example 2: generate univariate data for 5 conditions, each with 1
%   variable and 30 observations and calculate the correlation between
%   every pair of conditions (5 conditions = 10 correlations).
%       x = randn(30,5);
%       x(:,1:2) = x(:,1:2)-0.5*x(:,4:5);
%       [r,stats] = permucorr(x)
%
%   See also PERMUTTEST PERMUTTEST2 PERMUVARTEST2 DEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Blair RC, Higgins JJ, Karniski W, Kromrey JD (1994) A Study of
%           Multivariate Permutation Tests Which May Replace Hotelling's T2
%           Test in Prescribed Circumstances. Multivariate Behav Res,
%           29(2):141-163.
%       [2] Groppe DM, Urbach TP, Kutas M (2011) Mass univariate analysis
%           of event-related brain potentials/fields I: A critical tutorial
%           review. Psychophysiology, 48(12):1711-1725.
%       [3] Bliss CI (1967) Statistics in biology. New York: McGraw-Hill.
%       [4] Bishara AJ, Hittner JB, (2012) Testing the Significance of a
%           Correlation With Nonnormal Data: Comparison of Pearson,
%           Spearman, Transformation, and Resampling Approaches. Psychol
%           Methods, 17(3):399-417.
%       [5] Fisher RA (1921) On the "probable error" of a coefficient of
%           correlation deduced from a small sample. Metron, 1:3-32.
%       [6] Fieller EC, Hartley HO, Pearson ES (1957) Tests for Rank
%           Correlation Coefficients. I. Biometrika, 44(3/4):470-481.
%       [7] Bishara AJ, Hittner JB, (2017) Confidence intervals for
%           correlations when data are not normal. Behav Res, 49:294-309.

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
if arg.dim == 2 || isrow(x)
    x = x';
end
if ~isempty(y) && (arg.dim == 2 || isrow(y))
    y = y';
end

% Set up permutation test
if nargin<2 || isempty(y)
    warning('Comparing all columns of X in a correlation matrix...')
    [x,y] = paircols(x);
    arg.tail = 'both';
    arg.corrmat = true;
end
if size(x)~=size(y)
    error('X and Y must be the same size.')
end

% Use only rows with no NaN values if specified
switch arg.rows
    case 'complete'
        x = x(~any(isnan(y),2),:);
        y = y(~any(isnan(y),2),:);
        y = y(~any(isnan(x),2),:);
        x = x(~any(isnan(x),2),:);
    case 'all'
        if any(any(isnan(x))) || any(any(isnan(y)))
            error('X or Y is missing values. Set ROWS to ''complete''.')
        end
end

% Get data dimensions
[nobs,nvar] = size(x);

% Transform raw data to rank-orders if specified
switch arg.type
    case 'Rankit'
        x = norminv((tiedrank(x)-0.5)/nobs);
        y = norminv((tiedrank(y)-0.5)/nobs);
    case 'Spearman'
        x = tiedrank(x);
        y = tiedrank(y);
end

% Compute correlation coefficient
muxy = sum(x).*sum(y)/nobs;
sdxy = sqrt((sum(x.^2)-(sum(x).^2)/nobs).*(sum(y.^2)-(sum(y).^2)/nobs));
r = (sum(x.*y)-muxy)./sdxy;

% Execute if user requests adjusted test statistics
if nargout > 1

    % Permute data and generate distribution of r-values
    if nobs < 8
        warning('Computing all possible permutations due to small N.')
        arg.nperm = factorial(nobs);
        idx = perms(1:nobs)';
    else
        [~,idx] = sort(rand(nobs,arg.nperm));
    end
    rp = zeros(arg.nperm,nvar);
    for i = 1:arg.nperm
        rp(i,:) = (sum(x(idx(:,i),:).*y)-muxy)./sdxy;
    end

    % Compute adjusted test statistics using max statistic correction
    switch arg.tail
        case 'both'
            [~,idx] = max(abs(rp),[],2);
            csvar = [0;cumsum(ones(arg.nperm-1,1)*nvar)];
            rmax = rp'; rmax = rmax(idx+csvar);
            p = zeros(1,nvar);
            rcrit = zeros(2,1);
            if any(r>0)
                p(r>0) = 2*(sum(r(r>0)<=rmax)+1)/(arg.nperm+1);
            end
            if any(r<=0)
                p(r<=0) = 2*(sum(r(r<=0)>=rmax)+1)/(arg.nperm+1);
            end
            rcrit(1) = prctile(rmax,100*arg.alpha/2);
            rcrit(2) = prctile(rmax,100*(1-arg.alpha/2));
            zci = [norminv(arg.alpha/2);norminv(1-arg.alpha/2)];
            estal = (sum(rcrit(1)>=rmax)+sum(rcrit(2)<=rmax)+2)...
                /(arg.nperm+1);
        case 'left'
            rmax = min(rp,[],2);
            p = (sum(r>=rmax)+1)/(arg.nperm+1);
            rcrit = prctile(rmax,100*arg.alpha);
            zci = [-Inf(1,nvar);norminv(1-arg.alpha)];
            estal = (sum(rcrit>=rmax)+1)/(arg.nperm+1);
        case 'right'
            rmax = max(rp,[],2);
            p = (sum(r<=rmax)+1)/(arg.nperm+1);
            rcrit = prctile(rmax,100*(1-arg.alpha));
            zci = [norminv(arg.alpha);Inf(1,nvar)];
            estal = (sum(rcrit<=rmax)+1)/(arg.nperm+1);
    end
    p(isnan(rcrit)) = NaN;

    % Compute confidence intervals
    zr = log((1+r)./(1-r))/2;
    switch arg.type
        case 'Spearman'
            zci = zr+zci*sqrt(1.06/(nobs-3));
        otherwise
            zci = zr+zci/sqrt(nobs-3);
    end
    ci = (exp(2*zci)-1)./(exp(2*zci)+1);

    % Arrange test results in a matrix if specified
    if arg.corrmat
        p = vec2mat(p);
        ciLwr = vec2mat(ci(1,:));
        ciUpr = vec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
    end

    % Store values in structure
    stats = struct('p',p,'ci',ci,'rcrit',rcrit,'estal',estal);

end

% Execute if user requests unadjusted test statistics
if nargout > 2

    clear p rcrit ciLwr ciUpr ci estal

    % Compute unadjusted test statistics
    switch arg.tail
        case 'both'
            p = zeros(1,nvar);
            rcrit = zeros(2,nvar);
            if any(r>0)
                p(r>0) = 2*(sum(r(r>0)<=rp(r>0))+1)/(arg.nperm+1);
            end
            if any(r<=0)
                p(r<=0) = 2*(sum(r(r<=0)>=rp(r<=0))+1)/(arg.nperm+1);
            end
            rcrit(1,:) = prctile(rp,100*arg.alpha/2);
            rcrit(2,:) = prctile(rp,100-100*arg.alpha/2);
            zci = [norminv(arg.alpha/2);norminv(1-arg.alpha/2)];
            estal = (sum(rcrit(1,:)>=rp)+sum(rcrit(2,:)<=rp)+2)...
                /(arg.nperm+1);
        case 'left'
            p = (sum(r>=rp)+1)/(arg.nperm+1);
            rcrit = prctile(rp,100*arg.alpha);
            zci = [-Inf(1,nvar);norminv(1-arg.alpha)];
            estal = (sum(rcrit>=rp)+1)/(arg.nperm+1);
        case 'right'
            p = (sum(r<=rp)+1)/(arg.nperm+1);
            rcrit = prctile(rp,100*(1-arg.alpha));
            zci = [norminv(arg.alpha);Inf(1,nvar)];
            estal = (sum(rcrit<=rp)+1)/(arg.nperm+1);
    end
    p(isnan(rcrit(1,:))) = NaN;

    % Compute confidence intervals
    zr = log((1+r)./(1-r))/2;
    switch arg.type
        case 'Spearman'
            zci = zr+zci*sqrt(1.06/(nobs-3));
        otherwise
            zci = zr+zci/sqrt(nobs-3);
    end
    ci = (exp(2*zci)-1)./(exp(2*zci)+1);

    % Arrange test results in a matrix if specified
    if arg.corrmat
        p = vec2mat(p);
        rcrit = vec2mat(rcrit);
        ciLwr = vec2mat(ci(1,:));
        ciUpr = vec2mat(ci(2,:));
        ci = cat(3,ciLwr,ciUpr);
        ci = permute(ci,[3,1,2]);
        estal = vec2mat(estal);
    end

    % Store values in structure
    orig = struct('p',p,'ci',ci,'rcrit',rcrit,'estal',estal);

end

% Arrange r-values in a matrix if specified
if arg.corrmat
    r = vec2mat(r);
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

% Correlation type
typeOptions = {'Pearson','Spearman','Rankit'};
validFcn = @(x) any(validatestring(x,typeOptions));
addParameter(p,'type','Pearson',validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'corrmat',false,validFcn); % correlation matrix

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially matched strings
arg.tail = validatestring(arg.tail,tailOptions);
arg.rows = validatestring(arg.rows,rowsOptions);
arg.type = validatestring(arg.type,typeOptions);