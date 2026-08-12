function arg = ptparsevarargin(varargin)
%PTPARSEVARARGIN  Parse input arguments of PERMUTOOLS functions.
%   [PARAM1,PARAM2,...] = PTPARSEVARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   parses the input arguments of the main PERMUTOOLS functions.
%
%   See also PTVALIDATEPARAMIN.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Define valid string options
strOpts.tail = {'both','left','right','two'};
strOpts.type = {'','rank','anova1','kruskalwallis','anova2','alignrank',...
    'anovan','ranova1','ranova2','friedman','manova1','manova2',...
    'ttest','signrank','ttest2','ranksum','ftest','squarerank',...
    'pearson','spearman','rankit','freedmanlane','manly'};
strOpts.distribution = {'normal','binomial','poisson'};
strOpts.rows = {'all','complete'};
strOpts.compare = {'zero','pairwise'};
strOpts.vartype = {'equal','unequal'};
strOpts.effect = {'cohen','glass','cliff','meandiff','mediandiff'};

% Create parser object
p = inputParser;
p.CaseSensitive = false;

% Alpha level
errorMsg = 'It must be a scalar between 0 and 1.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x)&&x>0&&x<1,errorMsg);
addParameter(p,'alpha',0.05,validFcn);

% Dimension to work along
errorMsg = 'It must be a positive integer scalar within indexing range.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x)&&(x==1||x==2),errorMsg);
addParameter(p,'dim',1,validFcn);

% Resampling counts
errorMsg = 'It must be a positive integer scalar.';
validCount = @(x) assert(isnumeric(x)&&isscalar(x)&&x>0&&mod(x,1)==0,...
    errorMsg);
addParameter(p,'nperm',1e4,validCount);
addParameter(p,'nboot',1e4,validCount);

% Permutation generator seed
errorMsg = 'It must be an integer scalar.';
validFcn = @(x) assert(isnumeric(x)&&isscalar(x),errorMsg);
addParameter(p,'seed','shuffle',validFcn);

% Null hypothesis mean
errorMsg = 'It must be a numeric scalar or row vector.';
validFcn = @(x) assert(isnumeric(x),errorMsg);
addParameter(p,'m',0,validFcn);

% String options
addParameter(p,'type','',@(x) any(validatestring(x,strOpts.type)));
addParameter(p,'distribution','normal',@(x) any(validatestring(x,strOpts.distribution)));
addParameter(p,'tail','both',@(x) any(validatestring(x,strOpts.tail)));
addParameter(p,'rows','all',@(x) any(validatestring(x,strOpts.rows)));
addParameter(p,'compare','zero',@(x) any(validatestring(x,strOpts.compare)));
addParameter(p,'vartype','equal',@(x) any(validatestring(x,strOpts.vartype)));
addParameter(p,'effect','cohen',@(x) any(validatestring(x,strOpts.effect)));

% Boolean flags
errorMsg = 'It must be a numeric scalar (0,1) or logical.';
validFcn = @(x) assert(x==0||x==1||islogical(x),errorMsg);
addParameter(p,'correct',true,validFcn); % control FWER via max correction
addParameter(p,'intercept',true,validFcn); % include regression intercept
addParameter(p,'paired',true,validFcn); % compute effect size for paired samples
addParameter(p,'matrix',false,validFcn); % return results in a matrix
addParameter(p,'verbose',true,validFcn); % run in verbose mode

% Parse input
parse(p,varargin{1,1}{:});
arg = p.Results;

% Normalize all string options to canonical lowercase
fields = fieldnames(strOpts);
for i = 1:numel(fields)
    fn = fields{i};
    arg.(fn) = validatestring(arg.(fn), strOpts.(fn));
end