function run_permuvartest_examples
%RUN_PERMUVARTEST_EXAMPLES  Run one-sample bootstrap chi-squared test examples.
%   Generates random multivariate data for samples X. Each sample has 20
%   variables, each with a standard deviation of 1, except for the first 10
%   variables of X which have a standard deviation of 2. Each variable has
%   30 observations. One-sample bootstrap tests based on the chi-squared
%   statistic are performed on each variables of X for two-tailed, right-
%   tailed and left-tailed tests. The results are compared to those of the
%   equivalent parametric statistical test (one-sample test of variance)
%   using MATLAB's vartest.m.

%   See also PERMUVARTEST VARTEST.
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
nobs = 30; nvar = 20; v = 2;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = 'chisqtest';
test_metric = 'chisq-value';

% Generate random data
rng(42);
x = randn(nobs,nvar);
x(:,1:round(nvar/2)) = x(:,1:round(nvar/2))*sqrt(v);
if insert_nan
    x(1,1) = NaN;
end

f1 = figure('Name',['One-sample ',type,': statistic & CIs'],...
    'NumberTitle','off');
set(gcf,'color','w')
f2 = figure('Name',['One-sample ',type,': p-values'],'NumberTitle','off');
set(gcf,'color','w')

disp(type)

toc1 = zeros(numel(tail),1);
toc2 = zeros(numel(tail),1);
toc3 = zeros(numel(tail),1);

for i = 1:numel(tail)

    % Parametric test
    tic
    [~,p1,ci1,stats1] = vartest(x,v,'tail',tail{i});
    toc1(i) = toc;

    % Permutation test (uncorrected)
    tic
    [~,p2,ci2,stats2] = permuvartest(x,v,'tail',tail{i},'correct',0,...
        'verbose',0);
    toc2(i) = toc;

    % Permutation test (max-corrected)
    tic
    [~,p3,ci3,stats3] = permuvartest(x,v,'tail',tail{i},'correct',1,...
        'verbose',0);
    toc3(i) = toc;

    % Plot statistic & CIs
    figure(f1)
    subplot(3,2,i+i-1), hold on
    plot(xaxis,stats2.var,'LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),stats2.var(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),stats2.var(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,10]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 1
        legend('variance','95% CI (param.)','','95% CI (boot.)',...
            'Location','best')
    end
    subplot(3,2,i+i), hold on
    plot(xaxis,stats3.var,'LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci3,'--r')
    plot(xaxis(p1<=alpha),stats2.var(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p3<=alpha),stats2.var(p3<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,10]), box on, grid on
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
        legend('{\itp}-value (param.)','{\itp}-value (boot.)',...
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

% Plot descriptive statistics
figure('Name',['Paired ',type,': descriptive statistics'],...
    'NumberTitle','off');
set(gcf,'color','w')
plot(xaxis,stats1.chisqstat,xaxis,stats3.chisqstat,'--','LineWidth',2)
xlim([0,nvar+1]), box on, grid on
title('Test Statistic')
xlabel('variable')
ylabel(test_metric)
legend('param.','perm.','Location','best')

fprintf('Parametric (uncorrect): %.1f ms\n',mean(toc1)*1e3)
fprintf('Permutation (uncorrect): %.1f ms\n',mean(toc2)*1e3)
fprintf('Permutation (max-corr.): %.1f ms\n',mean(toc3)*1e3)