function run_permuchi2_examples
%RUN_PERMUCHI2_EXAMPLES  Run permutation Chi-square test examples.
%   Generates random multivariate categorical data for a sample X,
%   containing 20 variables. The first 10 variables contain a weak
%   association, while the remaining 10 variables are completely
%   independent. Each variable has 100 observations. Permutation tests of
%   independence based on the Chi-square statistic are performed on each
%   variable for two-tailed, right-tailed and left-tailed tests. The
%   results are compared to those of the equivalent parametric statistical
%   test (Chi-square test) using MATLAB's chi2cdf.m.
%
%   See also PERMUCHI2 CROSSTAB CHI2GOF.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS
%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

close all; clc;

% Set up experiment
insert_nan = false;
nobs = 100; nvar = 20; df = 1;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = 'permuchi2';

% Generate random categorical data with a weak effect in the first half
rng(42);
x = zeros(2,2,nvar);
varA = [ones(nobs/2,1);2*ones(nobs/2,1)];
for j = 1:nvar
    if j <= round(nvar/2)
        % Weak association: ~60% match
        varB = varA;
        flip = rand(nobs,1)>0.60;
        varB(flip) = 3-varB(flip);
    else
        % No association: ~50% match (chance)
        varB = randi([1,2],nobs,1);
    end
    x(:,:,j) = accumarray([varA,varB],1,[2,2]);
end
if insert_nan
    x(1,1) = NaN;
end

toc1 = zeros(numel(tail),1);
toc2 = zeros(numel(tail),1);
toc3 = zeros(numel(tail),1);

f1 = figure('Name',['One-sample ',type,': test statistic'],...
    'NumberTitle','off');
set(gcf,'color','w')
f2 = figure('Name',['One-sample ',type,': p-values'],'NumberTitle','off');
set(gcf,'color','w')

for i = 1:numel(tail)

    % Permutation test (uncorrected)
    tic
    [chi2,p2] = permuchi2(x,'tail',tail{i},'correct',0);
    toc2(i) = toc;

    % Parametric test
    tic
    [nrow,ncol,~] = size(x);
    [I,J] = ndgrid(1:nrow,1:ncol);
    p1 = zeros(1,nvar);
    chi1 = zeros(1,nvar);
    for j = 1:nvar
        xslice = x(:,:,j);
        var1 = repelem(I(:),xslice(:));
        var2 = repelem(J(:),xslice(:));
        switch tail{i}
            case 'right'
                [~,chi1(j),p1(j)] = crosstab(var1,var2);
            case 'left'
                [~,chi1(j)] = crosstab(var1,var2);
                p1(j) = chi2cdf(chi2(j),df);
            case 'both'
                [~,chi1(j)] = crosstab(var1,var2);
                pl = chi2cdf(chi2(j),df);
                pr = 1-chi2cdf(chi2(j),df);
                p1(j) = min(1,2*min(pl,pr));
        end
    end
    toc1(i) = toc;

    % Permutation test (max-corrected)
    tic
    [chi3,p3] = permuchi2(x,'tail',tail{i},'correct',1);
    toc3(i) = toc;

    % Plot statistic & CIs
    figure(f1)
    subplot(3,2,i+i-1), hold on
    plot(xaxis,chi1,xaxis,chi2,'--','LineWidth',2)
    plot(xaxis(p1<=alpha),chi2(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),chi2(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,25]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 1
        legend('Chi2 statistic','Location','best')
    end
    subplot(3,2,i+i), hold on
    plot(xaxis,chi1,xaxis,chi3,'--','LineWidth',2)
    plot(xaxis(p1<=alpha),chi3(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p3<=alpha),chi3(p3<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,25]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end

    % Plot p-values
    figure(f2)
    subplot(3,2,i+i-1), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 1
        legend('{\itp}-value (param.)','{\itp}-value (perm.)',...
            'Location','best')
    end
    subplot(3,2,i+i), hold on
    plot(xaxis,p1,'k',xaxis,p3,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end

end

fprintf('Parametric (uncorrect): %.1f ms\n',mean(toc1)*1e3)
fprintf('Permutation (uncorrect): %.1f ms\n',mean(toc2)*1e3)
fprintf('Permutation (max-corr.): %.1f ms\n',mean(toc3)*1e3)