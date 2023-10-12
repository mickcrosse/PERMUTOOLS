function run_permuttest_examples
%RUN_PERMUTTEST_EXAMPLES  Run paired permutation t-test examples.
%   Generates random multivariate data for 2 "dependent" samples X and Y.
%   Each sample has 20 variables, each with a mean value of 0, except for
%   the first 10 variables of Y which have a mean value of -1. Each
%   variable has 30 observations. Paired-sample permutation tests based on
%   the t-statistic are performed between the corresponding variables of
%   each sample for two-tailed, right-tailed and left-tailed tests. The
%   results are compared to those of the equivalent parametric statistical
%   tests (i.e. paired t-tests) using ttest.m.
%
%   See also PERMUTTEST TTEST.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

% Generate random data
rng(42);
x = randn(30,20);
y = randn(30,20);
y(:,1:10) = y(:,1:10)-1;
md = mean(x-y);
xaxis = 1:20;
tail = {'both','right','left'};
label = {'two','right','left'};

% Plot parametric & permutation CIs
figure('Name','Paired Test: mean difference & CIs','NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [h1,~,ci1] = ttest(x,y,'tail',tail{i});
    [h2,~,ci2] = permuttest(x,y,'tail',tail{i},'correct',0,'verbose',0);
    subplot(3,2,i+i-1), hold on
    plot(xaxis,md,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(logical(h1)),md(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),md(logical(h2)),'xr','LineWidth',2)
    xlim([0,21]), ylim([-3,3]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 2
        legend('mean difference','parametric CI','','permutation CI')
    end
    [h2,~,ci2] = permuttest(x,y,'tail',tail{i},'correct',1,'verbose',0);
    subplot(3,2,i+i), hold on
    plot(xaxis,md,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(logical(h1)),md(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),md(logical(h2)),'xr','LineWidth',2)
    xlim([0,21]), ylim([-3,3]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end
end

% Plot parametric & permutation p-values
figure('Name','Paired Test: p-values','NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [h1,p1] = ttest(x,y,'tail',tail{i});
    [h2,p2] = permuttest(x,y,'tail',tail{i},'correct',0,'verbose',0);
    subplot(3,2,i+i-1), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    plot(xaxis(logical(h1)),p1(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),p2(logical(h2)),'xr','LineWidth',2)
    xlim([0,21]), ylim([0,1]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 2
        legend('parametric {\itp}','permutation {\itp}')
    end
    [h2,p2] = permuttest(x,y,'tail',tail{i},'correct',1,'verbose',0);
    subplot(3,2,i+i), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    plot(xaxis(logical(h1)),p1(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),p2(logical(h2)),'xr','LineWidth',2)
    xlim([0,21]), ylim([0,1]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end
end