function run_bootvartest_examples
%RUN_BOOTVARTEST_EXAMPLES  Run one-sample bootstrapped chi-squared test examples.
%   Generates random multivariate data for samples X. Each sample has 20
%   variables, each with a standard deviation of 1, except for the first 10
%   variables of X which have a standard deviation of 2. Each variable has
%   30 observations. One-sample bootstrapped tests based on the chi-squared
%   statistic are performed on each variables of X for two-tailed, right-
%   tailed and left-tailed tests. The results are compared to those of the
%   equivalent parametric statistical tests (i.e. one-sample chi-squared 
%   tests) using vartest.m.

%   See also BOOTVARTEST VARTEST.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate random data
rng(42);
x = randn(30,20);
v = 2;
x(:,1:10) = x(:,1:10)*sqrt(v);
xaxis = 1:20; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};

% Plot parametric & bootstrapped CIs
figure('Name','One-sample Test: variance & CIs','NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [~,p1,ci1] = vartest(x,v,'tail',tail{i});
    [~,p2,ci2,stats2] = bootvartest(x,v,'tail',tail{i},'correct',0,...
        'verbose',0);
    subplot(3,2,i+i-1), hold on
    plot(xaxis,stats2.varx,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),stats2.varx(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),stats2.varx(p2<=alpha),'xr','LineWidth',2)
    xlim([0,21]), ylim([0,10]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 2
        legend('variance','95% CI (param.)','','95% CI (boot.)')
    end
    [~,p2,ci2,stats2] = bootvartest(x,v,'tail',tail{i},'correct',1,...
        'verbose',0);
    subplot(3,2,i+i), hold on
    plot(xaxis,stats2.varx,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),stats2.varx(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),stats2.varx(p2<=alpha),'xr','LineWidth',2)
    xlim([0,21]), ylim([0,10]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end
end

% Plot parametric & bootstrapped p-values
figure('Name','One-sample Test: p-values','NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [~,p1] = vartest(x,v,'tail',tail{i});
    [~,p2] = bootvartest(x,v,'tail',tail{i},'correct',0,'verbose',0);
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
        legend('{\itp}-value (param.)','{\itp}-value (boot.)')
    end
    [~,p2] = bootvartest(x,v,'tail',tail{i},'correct',1,'verbose',0);
    subplot(3,2,i+i), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,21]), ylim([0,1]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end
end