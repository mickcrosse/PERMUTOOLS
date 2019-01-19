function [y1,y2] = paircols(x)
%paircols pair matrix columns and output as two separate matrices
%   [Y1,Y2] = PAIRCOLS(X) returns matrices Y1 and Y2 whose paired columns
%   correspond to every combination of column pairs in X. For efficiency,
%   repeated column pairs are skipped.
%
%   See also VEC2MAT PERMTEST_T PERMTEST_T2 PERMTEST_F PERMTEST_CORR.
%
%   StatsTools https://github.com/mickcrosse/PERMUTOOLS

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Sep 2018; Last Revision: 14-Sep-2018

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