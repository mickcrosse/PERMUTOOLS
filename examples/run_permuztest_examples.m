function run_permuztest_examples
%RUN_PERMUZTEST_EXAMPLES  Run one-sample permutation Z-test examples.
%   Generates random multivariate data for a sample X, containing 20
%   variables, each with a mean value of 0, except for the first 10
%   variables which have a mean value of -1. Each variable has 30
%   observations. One-sample permutation tests based on the Z-statistic are
%   performed on each variable for two-tailed, right-tailed and left-tailed
%   tests. The results are compared to those of the equivalent parametric
%   statistical tests (i.e. paired Z-tests) using ztest.m.
%
%   See also PERMUZTEST ZTEST.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

close all; clc;

% Check version
info = ver;
isoctave = any(ismember({info.Name},'Octave'));

% Set up experiment
insert_nan = false;
nobs = 30; nvar = 20;
m = 0; sigma = 1;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = 'ztest';

% Generate random data
rng(42);
x = randn(nobs,nvar);
x(:,1:round(nvar/2)) = x(:,1:round(nvar/2))-1;
if insert_nan
    x(1,1) = NaN;
end

disp(type)

toc1 = zeros(numel(tail),1);
toc2 = zeros(numel(tail),1);
toc3 = zeros(numel(tail),1);

f1 = figure('Name',['One-sample ',type,': mean & CIs'],'NumberTitle','off');
set(gcf,'color','w')
f2 = figure('Name',['One-sample ',type,': p-values'],'NumberTitle','off');
set(gcf,'color','w')

for i = 1:numel(tail)

    % Parametric test
    if isoctave
        tic
        p1 = zeros(1,nvar);
        ci1 = zeros(2,nvar);
        for j = 1:nvar
            [~,p1(j),ci1(:,j)] = ztest(x(:,j),m,sigma,'tail',tail{i});
        end
        toc1(i) = toc;
    else
        tic
        [~,p1,ci1,zval1] = ztest(x,m,sigma,'tail',tail{i});
        sd1 = std(x);
        toc1(i) = toc;
    end

    % Permutation test (uncorrected)
    tic
    [~,p2,ci2,stats2] = permuztest(x,m,sigma,'tail',tail{i},'correct',0,...
        'verbose',0);
    toc2(i) = toc;

    % Permutation test (max-corrected)
    tic
    [~,p3,ci3,stats3] = permuztest(x,m,sigma,'tail',tail{i},'correct',1,...
        'verbose',0);
    toc3(i) = toc;

    % Plot central tendancy & CIs
    figure(f1)
    subplot(3,2,i+i-1), hold on
    plot(xaxis,stats2.mean,'LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),stats2.mean(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),stats2.mean(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([-2,2]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 1
        legend('mean','95% CI (param.)','','95% CI (perm.)',...
            'Location','best')
    end
    subplot(3,2,i+i), hold on
    plot(xaxis,stats3.mean,'LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci3,'--r')
    plot(xaxis(p1<=alpha),stats3.mean(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p3<=alpha),stats3.mean(p3<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([-2,2]), box on, grid on
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

% Plot descriptive statistics
figure('Name',['One-sample ',type,': descriptive statistics'],...
    'NumberTitle','off');
set(gcf,'color','w')
subplot(2,2,1), hold on
plot(xaxis,zval1,xaxis,stats2.zstat,'--','LineWidth',2)
xlim([0,nvar+1]), ylim([-7,7]), box on, grid on
title('Test Statistic')
xlabel('variable')
ylabel('z-value')
legend('param.','perm.','Location','best')
subplot(2,2,2), hold on
plot(xaxis,sd1,xaxis,stats2.sd,'--','LineWidth',2)
xlim([0,nvar+1]), ylim([0,2]), box on, grid on
title('Variance')
xlabel('variable')
ylabel('SD')

fprintf('Parametric (uncorrect): %.1f ms\n',mean(toc1)*1e3)
fprintf('Permutation (uncorrect): %.1f ms\n',mean(toc2)*1e3)
fprintf('Permutation (max-corr.): %.1f ms\n',mean(toc3)*1e3)