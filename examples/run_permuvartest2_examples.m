function run_permuvartest2_examples
%RUN_PERMUVARTEST2_EXAMPLES  Run two-sample permutation F-test examples.
%   Generates random multivariate data for 2 "independent" samples X and Y.
%   Each sample has 20 variables, each with a standard deviation of 1,
%   except for the first 10 variables of Y which have a standard deviation
%   of 2. Each variable has 30 observations. Two-sample permutation tests
%   based on the F-statistic are performed between the corresponding
%   variables of each sample for two-tailed, right-tailed and left-tailed
%   tests. The results are compared to those of the equivalent parametric
%   statistical tests (i.e. two-sample F-tests) using vartest2.m.
%
%   See also PERMUVARTEST2 VARTEST2.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

% Generate random data
rng(42);
x = randn(30,20);
y = randn(30,20);
y(:,1:10) = y(:,1:10)*2;
xaxis = 1:20;
tail = {'both','right','left'};

% Plot parametric & permutation CIs
figure
for i = 1:numel(tail)
    [h1,~,ci1,stats1] = vartest2(x,y,'tail',tail{i});
    [h2,~,ci2,stats2] = permuvartest2(x,y,'tail',tail{i},'correct',false);
    subplot(3,2,i+i-1)
    hold on
    plot(xaxis,stats2.fstat,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(logical(h1)),stats1.fstat(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),stats2.fstat(logical(h2)),'xr','LineWidth',2)
    hold off
    ylim([0,6])
    xlim([0,21])
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([tail{i},'-tailed'])
    if i == 2
        legend('F-statistic','parametric CIs','','permutation CIs')
    end
    [h2,~,ci2,stats2] = permuvartest2(x,y,'tail',tail{i},'correct',true);
    subplot(3,2,i+i)
    hold on
    plot(xaxis,stats2.fstat,'LineWidth',3)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(logical(h1)),stats1.fstat(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),stats2.fstat(logical(h2)),'xr','LineWidth',2)
    hold off
    ylim([0,6])
    xlim([0,21])
    if i == 1
        title('Corrected')
    elseif i == 3
        xlabel('variable')
    end
end

% Plot parametric & permutation p-values
figure
for i = 1:numel(tail)
    [h1,p1] = vartest2(x,y,'tail',tail{i});
    [h2,p2] = permuvartest2(x,y,'tail',tail{i},'correct',false);
    subplot(3,2,i+i-1)
    hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    plot(xaxis(logical(h1)),p1(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),p2(logical(h2)),'xr','LineWidth',2)
    hold off
    ylim([0,1])
    xlim([0,21])
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([tail{i},'-tailed'])
    if i == 2
        legend('parametric p-value','permutation p-value')
    end
    [h2,p2] = permuvartest2(x,y,'tail',tail{i},'correct',true);
    subplot(3,2,i+i)
    hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    plot(xaxis(logical(h1)),p1(logical(h1)),'ok','LineWidth',2)
    plot(xaxis(logical(h2)),p2(logical(h2)),'xr','LineWidth',2)
    hold off
    ylim([0,1])
    xlim([0,21])
    if i == 1
        title('Corrected')
    elseif i == 3
        xlabel('variable')
    end
end