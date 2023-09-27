function ptvalidateparamin(x,y,arg)
%PTVALIDATEPARAMIN  Validate input parameters of PERMUTOOLS functions.
%   PTVALIDATEPARAMIN(X,Y,ARG) validates the input parameters of the main
%   PERMUTOOLS function.
%
%   See also PTPARSEVARARGIN.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

if ~isnumeric(x)
    error('X must be numeric.')
elseif ~isnumeric(y)
    error('Y must be numeric or empty.')
end
if (arg.nperm<1e3 && arg.alpha<=0.05) || (arg.nperm<5e3 && arg.alpha<=0.01)
    warning('Number of permutations may be too low for chosen ALPHA.')
end