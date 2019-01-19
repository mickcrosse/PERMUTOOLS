function [alpha,nperm,tail,rows,sample,varx,type,m] = decode_varargin(varargin)
%decode_varargin decode input variable arguments
%   [PARAM1,PARAM2,...] = DECODE_VARARGIN('PARAM1',VAL1,'PARAM2',VAL2,...)
%   decodes the input variable arguments of various functions of the
%   PERMUTOOLS statistical toolbox. To define certain parameters, specify
%   both the parameter name and its value. Certain parameters only apply to
%   specific functions. See function documentation for help. Valid
%   parameters are the following:
%
%   Parameter   Value
%   'alpha'     a scalar between 0 and 1 specifying the significance level
%               as 100*ALPHA% (default=0.05)
%   'nperm'     a scalar specifying the number of permutations (default=
%               10,000 or all possible permutations for less than 14 obs.)
%   'tail'      a string specifying the alternative hypothesis
%                   'both'      two-tailed test (default)
%                   'right'     right-tailed test
%                   'left'      left-tailed test
%   'rows'      a string specifying the rows to use in the case of any
%               missing values (NaNs)
%                   'all'       use all rows, regardless of NaNs (default)
%                   'complete'  use only rows with no NaNs
%   'sample'    a string specifying whether to perform a one-sample test or
%               a paired-sample test when only X1 is entered
%                   'one'       compare each column of X1 to zero and
%                               output the results as a vector (default)
%                   'paired'    compare every pair of columns of X1 to each
%                               other and output the results as a matrix
%   'varx'      a string specifying the variance equivalence of X1 and X2
%                   'equal'   	assume equal variances (default)
%                   'unequal' 	assume unequal variances
%   'type'      string specifying the type of correlation measure
%                   'Pearson'   Pearson's correlation coefficient (default)
%                   'Spearman'  Spearman's rank correlation coefficient
%                   'Rankit'    Bliss's rankit correlation coefficient
%   'm'         a scalar or row vector specifying the mean of the null
%               hypothesis for each variable (default=0)
%
%   See also PERMTEST_T PERMTEST_T2 PERMTEST_F PERMTEST_CORR EFFECTSIZE_D.
%
%   StatsTools https://github.com/mickcrosse/PERMUTOOLS

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Sep 2018; Last Revision: 18-Jan-2019

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
if any(strcmpi(varargin,'varx')) && ~isempty(varargin{find(strcmpi(varargin,'varx'))+1})
    varx = varargin{find(strcmpi(varargin,'varx'))+1};
    if ~any(strcmpi(varx,{'equal','unequal'}))
        error('Invalid value for argument VARX. Valid values are: ''equal'', ''unequal''.')
    end
else
    varx = 'equal';
end
if any(strcmpi(varargin,'type')) && ~isempty(varargin{find(strcmpi(varargin,'type'))+1})
    type = varargin{find(strcmpi(varargin,'type'))+1};
    if ~any(strcmpi(type,{'Pearson','Spearman','Rankit'}))
        error('Invalid value for argument TYPE. Valid values are: ''Pearson'', ''Spearman''.')
    end
else
    type = 'Pearson';
end
if any(strcmpi(varargin,'m')) && ~isempty(varargin{find(strcmpi(varargin,'m'))+1})
    m = varargin{find(strcmpi(varargin,'m'))+1};
    if ~isnumeric(m) || any(isnan(m)) || any(isinf(m))
        error('M must be a scalar or vector of numeric values.')
    end
else
    m = 0;
end