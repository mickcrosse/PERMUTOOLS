function run_permuttest2_examples
%RUN_PERMUTTEST2_EXAMPLES  Run two-sample permutation t-test examples.
%   Generates random multivariate data for 2 "independent" samples X and Y.
%   Each sample has 20 variables, each with a mean value of 0, except for
%   the first 10 variables of Y which have a mean value of -1. Each
%   variable has 30 observations. Two-sample permutation tests based on the
%   t-statistic are performed between the corresponding variables of each
%   sample for two-tailed, right-tailed and left-tailed tests, as well as
%   samples of equal and unequal variances. The results are compared to
%   those of the equivalent parametric statistical tests (i.e. two-sample
%   t-tests) using ttest2.m.
%
%   See also PERMUTTEST2 TTEST2.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

% Generate random data
rng(42);
x = randn(30,20);
y = randn(30,20);
y(:,1:10) = y(:,1:10)-1;
diffxy = mean(x-y);
xaxis = 1:20;
tail = {'both','right','left'};
vartype = {'equal','unequal'};

for n = 1:numel(vartype)

    % Make variance unequal
    if n == 2
        y = y*1.25;
        diffxy = mean(x-y);
    end

    % Plot parametric & permutation CIs
    figure
    for i = 1:numel(tail)
        [h1,~,ci1] = ttest2(x,y,'tail',tail{i},'vartype',vartype{n});
        [h2,~,ci2] = permuttest2(x,y,'tail',tail{i},'vartype',vartype{n},...
            'correct',false);
        subplot(3,2,i+i-1)
        hold on
        plot(xaxis,diffxy,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(logical(h1)),diffxy(logical(h1)),'ok','LineWidth',2)
        plot(xaxis(logical(h2)),diffxy(logical(h2)),'xr','LineWidth',2)
        hold off
        ylim([-3,3])
        xlim([0,21])
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([tail{i},'-tailed'])
        if i == 2
            legend('mean(X−Y)','parametric CIs','','permutation CIs')
        end
        [h2,~,ci2] = permuttest2(x,y,'tail',tail{i},'vartype',vartype{n},...
            'correct',true);
        subplot(3,2,i+i)
        hold on
        plot(xaxis,diffxy,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(logical(h1)),diffxy(logical(h1)),'ok','LineWidth',2)
        plot(xaxis(logical(h2)),diffxy(logical(h2)),'xr','LineWidth',2)
        hold off
        ylim([-3,3])
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
        [h1,p1] = ttest2(x,y,'tail',tail{i},'vartype',vartype{n});
        [h2,p2] = permuttest2(x,y,'tail',tail{i},'vartype',vartype{n},...
            'correct',false);
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
        [h2,p2] = permuttest2(x,y,'tail',tail{i},'vartype',vartype{n},...
            'correct',true);
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

end