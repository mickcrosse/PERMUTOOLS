function st_plotpolyci(x,y,n,alpha,lnWth,lnCol,ciCol,ciFA)

%st_plotpolyci StatsTools plot confidence interval of N-degree polynomial
%
%   Inputs:
%   x     - column vector of data (observations by 1)
%   y     - column vector of data (observations by 1)
%
%   Optional Inputs:
%   n     - degree of polynomial fit (default=1)
%   alpha - significance level between 0 and 1 (default=0.05)
%   lnWth - line width of polynomial (default=1)
%   lnCol - line colour of polynomial (default=[0,0,0])
%   ciCol - colour of shaded confidence interval (default=[0.5,0.5,0.5])
%   ciFA  - face alpha level of shaded confidence interval (default=0.5)
%
%   See README for examples of use.
%
%   See also ST_TMAXPERM ST_TMAXPERM2 ST_RMAXPERM ST_BOXDOTPLOT.
%
%   StatsTools https://github.com/mickcrosse/StatsTools

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Feb 2018; Last Revision: 22-Mar-2018

if ~exist('n','var') || isempty(n)
    n = 1;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
if ~exist('lnWth','var') || isempty(lnWth)
    lnWth = 1;
end
if ~exist('lnCol','var') || isempty(lnCol)
    lnCol = [0,0,0];
end
if ~exist('ciCol','var') || isempty(ciCol)
    ciCol = [0.5,0.5,0.5];
end
if ~exist('ciFA','var') || isempty(ciFA)
    ciFA = 0.5;
end

% Orientate data
if size(x,1)==1 && size(x,2)>1
    x = x';
end
if size(y,1)==1 && size(y,2)>1
    y = y';
end

% Fit polynomial of degree N
[b,s] = polyfit(x,y,n);

% Compute confidence interval of polynomial fit
[y,d] = polyconf(b,x,s,'predopt','curve','alpha',alpha);

% Sort data X
[x,ix] = sort(x);

% Plot confidence interval
h = fill([x;flipud(x)],[y(ix)-d(ix);flipud(y(ix)+d(ix))],ciCol,'edgecolor','none');
set(h,'facealpha',ciFA)

% Plot polynomial
line(x,y(ix),'color',lnCol,'linewidth',lnWth)