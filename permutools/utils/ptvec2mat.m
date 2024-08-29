function [y] = ptvec2mat(x)
%PTVEC2MAT  Convert vector output to matrix format.
%   Y = PTVEC2MAT(X) returns a matrix Y by rearranging the values in vector
%   X according to their position as determined by PTPAIRCOLS. The values
%   in X may represent the output of some statistical test between every
%   pair of rows and columns in Y.
%
%   See also PTPAIRCOLS.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

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