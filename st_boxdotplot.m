function st_boxdotplot(x,y,dotCol,ciMeth,nBoot,bxCol,bxFA,bxWth,dotSz,jitSz,lnWth)

%st_boxdotplot StatsTools box dot plot
%
%   Inputs:
%   x      - vector of x-axis values (1 value per variable)
%   y      - cell array of data (1 cell per variable)
%   dotCol - matrix containing colour of datapoints (variables by RGB)
%
%   Optional Inputs:
%   ciMeth - method for estimating the confidence interval (CI)
%           'iqr'       - interquartile range (default)
%           'bootstrap' - bootstrap of the median using the bias corrected
%                         and accelerated percentile method
%   nBoot  - number of bootstraps (default=1e4)
%   bxCol  - CI box colour (default=[0.6,0.6,0.6])
%   bxFA   - CI box face alpha level (default=0.5)
%   bxWth  - CI box width (default=0.5)
%   dotSz  - dot size of datapoints (default=2)
%   jitSz  - value between 0 and 1 specifying size of jitter (default=0.2)
%   lnWth  - line width of median marker (default=2)
%
%   See README for examples of use.
%
%   See also ST_TMAXPERM ST_TMAXPERM2 ST_RMAXPERM ST_PLOTPOLYCI.
%
%   StatsTools https://github.com/mickcrosse/StatsTools

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Feb 2018; Last Revision: 22-Mar-2018

if strcmpi(ciMeth,'bootstrap') && (~exist('nBoot','var') || isempty(nBoot))
    nBoot = 1e4;
elseif exist('nBoot','var') && (~exist('ciMeth','var') || isempty(ciMeth))
    ciMeth = 'bootstrap';
elseif ~exist('ciMeth','var') || isempty(ciMeth)
    ciMeth = 'iqr';
end
if ~exist('bxCol','var') || isempty(bxCol)
    bxCol = [0.5,0.5,0.5];
end
if ~exist('bxFA','var') || isempty(bxFA)
    bxFA = 0.5;
end
if ~exist('bxWth','var') || isempty(bxWth)
    bxWth = 0.5;
end
if ~exist('dotSz','var') || isempty(dotSz)
    dotSz = 2;
end
if ~exist('jitSz','var') || isempty(jitSz)
    jitSz = 0.2;
end
if ~exist('lnWth','var') || isempty(lnWth)
    lnWth = 2;
end

% Loop through variables
for i = 1:length(y)
    
    % Generate jitter of datapoints
    jitter = randn(size(y{i}));
    jitter = jitter/max(abs(jitter))*bxWth*jitSz;
    
    % Plot dots
    scatter(x(i)+jitter,y{i},dotSz,dotCol(i,:),'fill')
    
    % Compute CI
    if strcmpi(ciMeth,'iqr')
        ci = prctile(y{i},[25,75]);
    elseif strcmpi(ciMeth,'bootstrap')
        ci = bootci(nBoot,@median,y{i});
    end
    
    % Create CI box
    X = [x(i)-bxWth/2,x(i)-bxWth/2,x(i)+bxWth/2,x(i)+bxWth/2];
    Y = [ci(1),ci(2),ci(2),ci(1)];
    
    % Plot CI box
    fill(X,Y,bxCol,'edgecolor','none','facealpha',bxFA)
    
    % Plot median line
    plot([x(i)-bxWth/1.8,x(i)+bxWth/1.8],[median(y{i}),median(y{i})],'k','linewidth',lnWth)
    
end