function [f,p,ci,stats,tbl,dist] = permuranova2(x,varargin)
%PERMURANOVA2  Two-way permutation repeated-measures ANOVA and Friedman test.
%   F = PERMURANOVA2(X) performs a two-way permutation-based repeated-
%   measures analysis of variance (RM ANOVA) for comparing the means of
%   data in a 3D matrix X (where rows represent subjects, columns represent
%   Factor A, and pages/depth represent Factor B), and returns the test
%   statistics F for Factor A, Factor B, and the A x B interaction.
%
%   For non-normally distributed data, the raw data may be transformed
%   to rank orders in order to compute a two-way Friedman test by setting
%   the 'type' parameter to 'friedman' or 'rank'.
%
%   PERMURANOVA2 requires a fully balanced design and does not support
%   missing values (NaNs).
%
%   [F,P] = PERMURANOVA2(...) returns the probabilities (i.e. p-values) of
%   observing the given results by chance if the null hypothesis is true.
%   As the null distributions are generated empirically by permuting the
%   data within each subject, no assumption is made about distributions.
%
%   [F,P,CI] = PERMURANOVA2(...) returns 100*(1-ALPHA)% confidence
%   intervals (CIs) for the true difference of population means.
%
%   [F,P,CI,STATS] = PERMURANOVA2(...) returns a structure with the
%   following fields:
%       'source'    -- the function used to compute the ANOVA
%       'method'    -- the type of test performed
%       'sigmasq'   -- the error mean squares [ErrorA,ErrorB,ErrorAxB]
%       'df'        -- the error degrees of freedom [dfEA,dfEB,dfEAB]
%
%   [F,P,CI,STATS,TBL] = PERMURANOVA2(...) returns the ANOVA table contents
%   as a cell array.
%
%   [F,P,CI,STATS,TBL,DIST] = PERMURANOVA2(...) returns the permuted
%   sampling distributions of the test statistics.
%
%   [...] = PERMURANOVA2(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values. Valid parameters are the
%   following:
%
%       Parameter   Value
%       'type'      A string specifying the type of permutation test to
%                   perform:
%                       'ranova2'   two-way repeated-measures ANOVA (default)
%                       'friedman'  two-way Friedman rank test
%       'alpha'     A scalar between 0 and 1 specifying the significance
%                   level as 100*ALPHA% (default=0.05).
%       'nperm'     An integer scalar specifying the number of permutations
%                   (default=10,000).
%       'seed'      An integer scalar specifying the seed value used to
%                   initialise the permutation generator.
%
%   See also RANOVA FRIEDMAN FITRM PERMURANOVA1 PERMUANOVA2.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Parse input arguments
arg = ptparsevarargin(varargin);
if isempty(arg.type)
    arg.type = 'ranova2';
end

% Validate input parameters
ptvalidateparamin(x,[],arg)

% Check matrix dimensions
if ndims(x) ~= 3
    error('Two-way RM ANOVA requires a 3D matrix (Subjects x Factor A x Factor B).')
end

% Balanced design check
if any(isnan(x(:)))
    error('Repeated-measures ANOVA requires a perfectly balanced design. Data must not contain NaNs.')
end

% Get data dimensions
[n,a,b] = size(x);
nobs = n*a*b;

% Transform raw data to rank orders if specified
switch arg.type
    case 'ranova2'
    case {'friedman','rank'}
        x2d = reshape(x,n,a*b);
        x2d = tiedrank(x2d')';
        x = reshape(x2d,n,a,b);
    otherwise
        error('The TYPE parameter value must be RANOVA2, or FRIEDMAN.')
end

% Compute grand mean
gm = sum(x,'all')/nobs;

% Compute marginal means
subjmeans = mean(x,[2,3]);
Ameans = mean(x,[1,3]);
Bmeans = mean(x,[1,2]);
ABmeans = mean(x,1);
subjAmeans = mean(x,3);
subjBmeans = mean(x,2);

% Compute sums of squares
sst = sum((x-gm).^2,'all');
sssubj = a*b*sum((subjmeans-gm).^2,'all');
ssa = n*b*sum((Ameans-gm).^2,'all');
ssb = n*a*sum((Bmeans-gm).^2,'all');
ssab = n*sum((ABmeans-Ameans-Bmeans+gm).^2,'all');

% Compute specific error sums of squares
ssea = b*sum((subjAmeans-subjmeans-Ameans+gm).^2,'all');
sseb = a*sum((subjBmeans-subjmeans-Bmeans+gm).^2,'all');
sseab = sst-sssubj-ssa-ssb-ssab-ssea-sseb;

% Compute degrees of freedom
dfa = a-1;
dfb = b-1;
dfab = (a-1)*(b-1);
dfea = (n-1)*(a-1);
dfeb = (n-1)*(b-1);
dfeab = (n-1)*(a-1)*(b-1);

% Compute mean squares
msa = ssa/dfa;
msb = ssb/dfb;
msab = ssab/dfab;
msea = ssea/dfea;
mseb = sseb/dfeb;
mseab = sseab/dfeab;

% Compute F-statistics [Factor A, Factor B, Interaction]
fa = msa/msea;
fb = msb/mseb;
fab = msab/mseab;
f = [fa, fb, fab];

if nargout > 1

    % Generate random permutations
    rng(arg.seed);
    x2d = reshape(x,n,a*b);
    rowidx = repmat((1:n)',1,a*b);

    % Estimate sampling distributions
    dista = zeros(arg.nperm,1);
    distb = zeros(arg.nperm,1);
    distab = zeros(arg.nperm,1);
    for i = 1:arg.nperm
        [~,idx] = sort(rand(n,a*b),2);
        linidx = rowidx+(idx-1)*n;
        xp2d = x2d(linidx);
        xp = reshape(xp2d,n,a,b);
        Ap = mean(xp,[1,3]);
        Bp = mean(xp,[1,2]);
        ABp = mean(xp,1);
        subjAp = mean(xp,3);
        subjBp = mean(xp,2);
        ssap = n*b*sum((Ap-gm).^2,'all');
        ssbp = n*a*sum((Bp-gm).^2,'all');
        ssabp = n*sum((ABp-Ap-Bp+gm).^2,'all');
        sseap = b*sum((subjAp-subjmeans-Ap+gm).^2,'all');
        ssebp = a*sum((subjBp-subjmeans-Bp+gm).^2,'all');
        sseabp = sst-sssubj-ssap-ssbp-ssabp-sseap-ssebp;
        dista(i) = (ssap/dfa)/(sseap/dfea);
        distb(i) = (ssbp/dfb)/(ssebp/dfeb);
        distab(i) = (ssabp/dfab)/(sseabp/dfeab);
    end

    % Compute p-values
    pa = (sum(fa<=dista)+1)/(arg.nperm+1);
    pb = (sum(fb<=distb)+1)/(arg.nperm+1);
    pab = (sum(fab<=distab)+1)/(arg.nperm+1);
    p = [pa,pb,pab];
    dist = [dista,distb,distab];

end

% Compute confidence intervals
if nargout > 2
    crita = prctile(dista,100*(1-arg.alpha));
    critb = prctile(distb,100*(1-arg.alpha));
    critab = prctile(distab,100*(1-arg.alpha));
    ci = [fa/crita,fb/critb,fab/critab;Inf,Inf,Inf];
end

% Store statistics in a structure
if nargout > 3
    stats.source = 'permuranova2';
    stats.method = arg.type;
    stats.sigmasq = [msea,mseb,mseab];
    stats.df = [dfea,dfeb,dfeab];
end

% Create ANOVA table
if nargout > 4
    tbl = {
        'Source','SS','df','MS','F','Prob>F';
        'Factor A',ssa,dfa,msa,fa,pa;
        'Factor B',ssb,dfb,msb,fb,pb;
        'A x B',ssab,dfab,msab,fab,pab;
        'Error(A)',ssea,dfea,msea,[],[];
        'Error(B)',sseb,dfeb,mseb,[],[];
        'Error(AxB)',sseab,dfeab,mseab,[],[];
        'Total',sst,nobs-1,[],[],[]};
end