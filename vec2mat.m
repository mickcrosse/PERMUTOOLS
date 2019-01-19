function [y] = vec2mat(x)
%vec2mat convert vector output to matrix format
%   [Y] = VEC2MAT(X) returns a matrix Y by rearranging the values in vector
%   X according to their position as determined by PAIRCOLS. The values in
%   X may represent the output of some statistical test between every pair
%   of rows and columns in Y.
%
%   See also PAIRCOLS PERMTEST_T PERMTEST_T2 PERMTEST_F PERMTEST_CORR.
%
%   StatsTools https://github.com/mickcrosse/PERMUTOOLS

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Sep 2018; Last Revision: 14-Sep-2018

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