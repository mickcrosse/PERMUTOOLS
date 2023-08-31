function [d,ci,stats] = deffectsize(x1,x2,dep,varx,bias,nboot,alpha)
%DEFFECTSIZE  Effect size measure with bootsrapped confidence intervals
%   D = DEFFECTSIZE(X1,X2) returns the effect size measure based on
%   Cohen's d, bias-corrected according to sample size (Hedges and Olkin, 
%   1985).

% 	If X1 and X2 are matrices, multiple effect size 
%   measures are computed simultaneously between each corresponding pair of 
%   columns in X1 and X2.
%   This function treats NaNs as missing values, and ignores them.

%   [...,CI] = DEFFECTSIZE(...) returns the confidence intervals using a bias-corrected and accelerated bootstrap approach (Davison and Hinkley, 1997). 
% 
%   [...,PARAM] = DEFFECTSIZE(...) returns some statistical parameters in a
%   structure containing the following fields:
%       df      the degrees of freedom of each test
%       sd    	the estimated population standard deviation of X1, or of
%               X1-X2 for a paired test
% 
%   See also PERMTEST_T PERMTEST_T2 PERMTEST_F2 PERMTEST_CORR.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Hentschke H, Stuttgen MC (2011) Computation of measures of 
%           effect size for neuroscience data sets. Eur J Neurosci, 
%           34:1887–1894.

%   Author: Mick Crosse
%   Email: mickcrosse@gmail.com
%   Cognitive Neurophysiology Laboratory,
%   Albert Einstein College of Medicine, NY
%   Jan 2018; Last Revision: 18-Jan-2019

% Decode input variable arguments
[alpha,nboot,bias,varx] = decode_varargin(varargin);

if dep==1 && strcmpi(varx,'unequal')
    warning('Assuming unequal variance for dependent samples')
end

var1 = (sum(x1.^2)-(sm1.^2)/nobs1)/df1;
var2 = (sum(x2.^2)-(sm2.^2)/nobs2)/df2;

% Compute effect size
if dep==0

    % Use only rows with no NaN values if specified
    if strcmpi(rows,'complete')
        x1 = x1(~any(isnan(x1),2),:);
        x2 = x2(~any(isnan(x2),2),:);
    end

    % Concatenate samples
    x = [x1;x2];
    
    % Get data dimensions, ignoring NaNs
    nobs = sum(~isnan(x));    
    nobs1 = sum(~isnan(x1));
    nobs2 = sum(~isnan(x2));

    % Compute degrees of freedom
    df1 = nobs1-1; 
	df2 = nobs2-1;
    dfp = nobs./(nobs1.*nobs2);
    
    % Compute Cohen's d
    sm1 = nansum(x1); sm2 = nansum(x2);
    mx = sm1/nobs1-sm2/nobs2;
    if strcmpi(varx,'equal')
        df = nobs-2;
        sd = sqrt((df1.*var1+df2.*var2)./df);
        se = sd.*sqrt(dfp);
        d = mx./sd;
    elseif strcmpi(varx,'unequal')
        se2x1 = var1/nobs1;
        se2x2 = var2/nobs2;
        df = (se2x1+se2x2).^2./(se2x1.^2/df1+se2x2.^2/df2);
        se = sqrt(se2x1+se2x2);
        d = mx./sqrt(var1);
    end
    
    % Compute t-statistic
    tstat = mx./se;
    
elseif dep==1
    
    % Compute difference between samples
    x = x1-x2;
    
    % Use only rows with no NaN values if specified
    if strcmpi(rows,'complete')
        x = x(~any(isnan(x),2),:);
    end
    
    % Get data dimensions, ignoring NaNs
    nobs = sum(~isnan(x));
    
    % Compute degrees of freedom
    df = nobs-1;
        
    % Compute t-statistic
    sd = nanstd(x);
    mx = nansum(x)./nobs;
    se = sd./sqrt(nobs);
    tstat = mx./se;
    
    % Compute Cohen's d
    d = tstat.*sqrt(2*sd.^2./(nobs.*(var1+var2)));
    
end

% Compute bias-corrected values
if bias==1 && strcmpi(varx,'equal')
    factor = 1-3./(4.*(nobs1+nobs2)-9);
    d = d*factor;
end

% Compute approximate CIs
if ~exist('nboot','var') || isempty(nboot)
    if dep==0
        if strcmpi(varx,'equal')
            se = sqrt((nobs1+nobs2)/(nobs1*nobs2)+(d.^2/(2*nobs1+2*nobs2-4)));
        elseif strcmpi(varx,'unequal')
            se = sqrt(d.^2/(2*nobs2-2)+(nobs1+nobs2)/(nobs1*nobs2));
        end
    elseif dep==1
        se = sqrt((2-2*corr(x1,x2))/nobs1+d.^2/(2*df1));
    end
    invt = -tinv(alpha/2,df);
    ci = cat(1,d-invt.*se,d+invt.*se);
    if bias==1 && strcmpi(varx,'equal') && dep==0
        ci = ci*factor;
    end
end

% Compute bootstrapped CIs
if ~isempty(nboot)
    px1 = ceil(nobs1*rand([nobs1,nboot]));
    x1 = x1(px1);
    if dep==0
        px2 = ceil(nobs2*rand([nobs2,nboot]));
        x2 = x2(px2);
        nobs = nobs1+nobs2;
        sm1 = sum(x1); sm2 = sum(x2);
        var1 = (sum(x1.^2)-(sm1.^2)/nobs1)/df1;
        var2 = (sum(x2.^2)-(sm2.^2)/nobs2)/df2;
        mx = sm1/nobs1-sm2/nobs2;
        if strcmpi(varx,'equal')
            df = nobs-2;
            sd = sqrt((df1.*var1+df2.*var2)./df);
            dp = mx./sd;
        elseif strcmpi(varx,'unequal')
            dp = mx./sqrt(var1);
        end
    elseif dep==1
        x2 = x2(px1);
        x = x1-x2;
        nobs = size(x,1);
        mx = sum(x)/nobs;
        se = std(x)/sqrt(nobs);
        dp = mx./se.*sqrt(2*std(x).^2./(nobs.*(var1+var2)));
    end
    ci = prctile(dp,[alpha/2,1-alpha/2]'*100);
end

% Execute if user specifies general data statistics
if nargout > 2
    if strcmpi(varx,'equal')
        stats = struct('tstat',tstat,'df',df,'sd',sd);
    elseif strcmpi(varx,'unequal')
        stats = struct('tstat',tstat,'df',df,'sd',sqrt([var1;var2]));
    end
end
