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

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate random data
rng(42);
x = randn(30,20);
x(:,1:10) = x(:,1:10)-1;
m = 0; sigma = 1;
xaxis = 1:20; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};

% Plot parametric & permutation CIs
figure('Name','One-sample Test: mean value & CIs','NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [~,p1,ci1] = ztest(x,m,sigma,'tail',tail{i});
    [~,p2,ci2,stats2] = permuztest(x,m,sigma,'tail',tail{i},'correct',0,...
        'verbose',0);
    subplot(3,2,i+i-1), hold on
    plot(xaxis,stats2.mu,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),stats2.mu(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),stats2.mu(p2<=alpha),'xr','LineWidth',2)
    xlim([0,21]), ylim([-2,2]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 2
        legend('mean value','95% CI (param.)','','95% CI (perm.)')
    end
    [~,p2,ci2,stats2] = permuztest(x,m,sigma,'tail',tail{i},'correct',1,...
        'verbose',0);
    subplot(3,2,i+i), hold on
    plot(xaxis,stats2.mu,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),stats2.mu(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),stats2.mu(p2<=alpha),'xr','LineWidth',2)
    xlim([0,21]), ylim([-2,2]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end
end

% Plot parametric & permutation p-values
figure('Name','One-sample Test: p-values','NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [~,p1] = ztest(x,m,sigma,'tail',tail{i});
    [~,p2] = permuztest(x,m,sigma,'tail',tail{i},'correct',0,'verbose',0);
    subplot(3,2,i+i-1), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,21]), ylim([0,1]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 2
        legend('{\itp}-value (param.)','{\itp}-value (perm.)')
    end
    [~,p2] = permuztest(x,m,sigma,'tail',tail{i},'correct',1,'verbose',0);
    subplot(3,2,i+i), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,21]), ylim([0,1]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end
end