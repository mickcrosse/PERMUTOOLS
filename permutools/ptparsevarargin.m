function arg = ptparsevarargin(varargin)
%PTPARSEVARARGIN  Parse input arguments of PERMUTOOLS functions.
%   [PARAM1,PARAM2,...] = PTPARSEVARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   parses the input arguments of the main PERMUTOOLS function.
%
%   See also PTVALIDATEPARAMIN.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

% Create parser object
p = inputParser;

% Alpha level
errorMsg = 'It must be a scalar between 0 and 1.';
validFcn = @(x) assert(x>0&&x<1,errorMsg);
addParameter(p,'alpha',0.05,validFcn);

% Dimension to work along
errorMsg = 'It must be a positive integer scalar within indexing range.';
validFcn = @(x) assert(x==1||x==2,errorMsg);
addParameter(p,'dim',1,validFcn);

% Alternative hypothesis
tailOptions = {'left','both','right'};
validFcn = @(x) any(validatestring(x,tailOptions));
addParameter(p,'tail','both',validFcn);

% Number of permutations
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x)&&x>0,errorMsg);
addParameter(p,'nperm',1e4,validFcn);

% Number of bootstraps
errorMsg = 'It must be a positive integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x)&&x>0,errorMsg);
addParameter(p,'nboot',1e4,validFcn);

% Rows to use if NaNs
rowsOptions = {'all','complete'};
validFcn = @(x) any(validatestring(x,rowsOptions));
addParameter(p,'rows','all',validFcn);

% Permutation generator seed
errorMsg = 'It must be an integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'seed','shuffle',validFcn);

% Null hypothesis mean
errorMsg = 'It must be a numeric scalar or row vector.';
validFcn = @(x) assert(isnumeric(x),errorMsg);
addParameter(p,'m',0,validFcn);

% Test type
testOptions = {'one','pairwise'};
validFcn = @(x) any(validatestring(x,testOptions));
addParameter(p,'test','one',validFcn);

% Variance equivalence
vartypeOptions = {'equal','unequal'};
validFcn = @(x) any(validatestring(x,vartypeOptions));
addParameter(p,'vartype','equal',validFcn);

% Correlation type
typeOptions = {'Pearson','Spearman','Rankit'};
validFcn = @(x) any(validatestring(x,typeOptions));
addParameter(p,'type','Pearson',validFcn);

% Effect size measure
effectOptions = {'Cohen','Glass','Cliff','meandiff','mediandiff'};
validFcn = @(x) any(validatestring(x,effectOptions));
addParameter(p,'effect','Cohen',validFcn);

% Boolean arguments
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'mat',false,validFcn); % matrix conversion
addParameter(p,'correct',true,validFcn); % correction
addParameter(p,'paired',true,validFcn); % paired samples
addParameter(p,'verbose',true,validFcn); % verbose mode

% Parse input arguments
parse(p,varargin{1,1}{:});
arg = p.Results;

% Redefine partially matched strings
arg.tail = validatestring(arg.tail,tailOptions);
arg.rows = validatestring(arg.rows,rowsOptions);
arg.test = validatestring(arg.test,testOptions);
arg.vartype = validatestring(arg.vartype,vartypeOptions);
arg.type = validatestring(arg.type,typeOptions);
arg.effect = validatestring(arg.effect,effectOptions);