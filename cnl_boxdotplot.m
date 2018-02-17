function [med,ymin,ymax] = cnl_boxdotplot(x,y,bxTyp,nBoot,bxWth,lnWth,dotSz,color)
%cnl_boxdotplot CNL Toolbox box-dotplot
% 
%   Inputs:
%   x     - x-value for datapoints (scaler) 
%   y     - y-values (column vector)
%   bxTyp - CI box type {'bootstrap'==BCa boostrapping (default), 'iqr'==
%           inter-quartile range}
%   nBoot - number of bootstraps
%   bxWth - width of CI box
%   lnWth - width of horizontal median line
%   dotSz - size of y-value dots
%   color - colour of y-value dots
%   
%   Outputs:
%   med    - median value 
%   ymin   - median value 
%   ymin   - median value 

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Feb 2018; Last Revision: 16-Feb-2018

if ~isscaler(x)
    x = x(1);
end
if ~exist('bxTyp','var') || isempty(bxTyp)
    bxTyp = 'bootstrap';
end
if ~exist('nPerm','var') || isempty(nBoot)
    nBoot = 1e4;
end
if ~exist('bxWth','var') || isempty(bxWth)
    bxWth = 0.5;
end
if ~exist('lnWth','var') || isempty(lnWth)
    lnWth = 2;
end
if ~exist('dotSz','var') || isempty(dotSz)
    dotSz = 2;
end
if ~exist('color','var') || isempty(color)
    color = [0,0.447,0.741];
end

% Compute median vlaue
med = median(y);

% Define box type
if strcmp(bxTyp,'iqr')
    ymin = prctile(y,25);
    ymax = prctile(y,75);
elseif strcmp(bxTyp,'bootstrap')
    ci = bootci(nBoot,@median,y);
    ymin = ci(1);
    ymax = ci(2);
end

% Create grey box
X = [x-bxWth/2,x-bxWth/2,x+bxWth/2,x+bxWth/2];
Y = [ymin,ymax,ymax,ymin];
C = [0.6,0.6,0.6];

% Plot grey box
fill(X,Y,C,'edgecolor','none')

% Plot y-value dots
scatter(ones(size(y))*x,y,dotSz,color,'fill')

% Plot median horizontal line
plot([x-bxWth/1.8,x+bxWth/1.8],[med,med],'k','linewidth',lnWth)