function [f,p,ci,stats,tbl,dist] = permuanova1(x,group,varargin)
%PERMUANOVA1  One-way permutation-based analysis of variance (ANOVA).
%   F = PERMUANOVA1(X) performs a one-way permutation-based ANOVA for
%   comparing the means of two or more groups of data in matrix X, and
%   returns the test statistic. The columns of X can have different lengths
%   by in including NaN values.
%
%   PERMUANOVA1 treats NaNs as missing values, and ignores them.
%
%   F = PERMUANOVA1(X,GROUP) groups the columns of X according to the
%   values in vector GROUP. The values in GROUP can be categorical,
%   numeric, logical or strings, with each value corresponding to a column
%   of X. Columns with the same corresponding group values are placed in
%   the same group.
%
%   [F,P] = PERMUANOVA1(...) returns the probability (i.e. p-value) of
%   observing the given result by chance if the null hypothesis (that the
%   means of the groups are equal) is true. As the null distribution is
%   generated empirically by permuting the data, no assumption is made
%   about the shape of the distributions that the data come from, except
%   that they have equal variances.
%
%   [F,P,CI] = PERMUANOVA1(...) returns a 100*(1-ALPHA)% confidence
%   interval (CI) for the true difference of population means.
%
%   [F,P,CI,STATS] = PERMUANOVA1(...) returns a structure with the
%   following fields:
%       'gnames'    -- the group names
%       'n'         -- the group sample sizes
%       'source'    -- the function used to compute the ANOVA
%       'means'     -- the group means
%       'df'        -- the error degrees of freedom
%       's'         -- the root mean square
%
%   [F,P,CI,STATS,TBL] = PERMUANOVA1(...) returns the ANOVA table contents
%   as a cell array.
%
%   [F,P,CI,STATS,TBL,DIST] = PERMUANOVA1(...) returns the permuted
%   sampling distribution of the test statistic.
%
%   [...] = PERMUANOVA1(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'dim'       A scalar specifying the dimension to work along: pass
%                   in 1 to work along the columns (default), or 2 to work
%                   along the rows. Applies to both X and Y.
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator. By default, the
%                   generator is initialised based on the current time,
%                   resulting in a different permutation on each call.
%
%   See also ANOVA1 PERMUANOVA2 PERMUTTEST2 BOOTEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

if nargin < 2 || isempty(group)
    group = 1:size(x,2);
end
[group,gnames] = grp2idx(group);

% Parse input arguments
arg = ptparsevarargin(varargin);

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Orient data column-wise
if arg.dim==2 || iscolumn(x)
    x = x';
end

% For efficiency, only omit NaNs if necessary
if any(isnan(x(:)))
    nanflag = 'omitmissing';
else
    nanflag = 'includemissing';
end

% Get data dimensions, ignoring NaNs
shapex = size(x);
nobs = sum(~isnan(x),'all');

% Compute grand mean
gm = sum(x,'all',nanflag)/nobs;

% Compute within and between-group sum of squares
groups = unique(group)';
n = zeros(1,numel(groups));
means = zeros(1,numel(groups));
ess = 0; rss = 0;
for i = groups
    xg = x(:,group==i);
    n(i) = sum(~isnan(xg),'all',nanflag);
    means(i) = sum(xg,'all',nanflag)/n(i);
    ess = ess + sum((xg-means(i)).^2,'all',nanflag);
    rss = rss + n(i)*sum((means(i)-gm).^2,'all',nanflag);
end

% Compute total sum of squares
tss = ess+rss;

% Compute degrees of freedom
dft = nobs-1;
dfr = numel(groups)-1;
dfe = dft-dfr;

% Compute mean squares
msr = rss/dfr;
mse = ess/dfe;

% Compute F-statistic
f = msr/mse;

if nargout > 1

    % Concatenate groups
    x = x(:);

    % Generate random permutations
    rng(arg.seed);
    maxnobs = numel(x);
    [~,idx] = sort(rand(maxnobs,arg.nperm));

    % Estimate sampling distribution
    dist = zeros(arg.nperm,1);
    for i = 1:arg.nperm
        xp = x(idx(:,i));
        xp = reshape(xp,shapex);
        essp = 0; rssp = 0;
        for j = groups
            xg = xp(:,group==j);
            np = sum(~isnan(xg),'all',nanflag);
            mp = sum(xg,'all',nanflag)/np;
            essp = essp + sum((xg-mp).^2,'all',nanflag);
            rssp = rssp + np*sum((mp-gm).^2,'all',nanflag);
        end
        dist(i) = (rssp/dfr)/(essp/dfe);
    end

    % Compute p-value & CI
    p = (sum(f<=dist)+1)/(arg.nperm+1);
    if nargout > 2
        crit = prctile(dist,100*(1-arg.alpha));
        ci = [f./crit;Inf];
    end

end

% Store statistics in a structure
if nargout > 3
    stats.gnames = gnames;
    stats.n = n;
    stats.source = 'permuanova1';
    stats.means = means;
    stats.df = dfe;
    stats.s = sqrt(mse);
end

% Create ANOVA table
if nargout > 4
    tbl = {
        'Source','SS','df','MS','F','Prob>F';
        'Groups',rss,dfr,msr,f,p;
        'Error',ess,dfe,mse,[],[];
        'Total',tss,dft,[],[],[]};
end