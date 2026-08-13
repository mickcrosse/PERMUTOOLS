function run_permuttest_examples
%RUN_PERMUTTEST_EXAMPLES  Run paired permutation t-test examples.
%   Generates random multivariate data for 2 "dependent" samples X and Y.
%   Each sample has 20 variables, each with a mean value of 0, except for
%   the first 10 variables of Y which have a mean value of -1. Each
%   variable has 30 observations. Paired-sample permutation tests based on
%   the t-statistic are performed between the corresponding variables of
%   each sample for two-tailed, right-tailed and left-tailed tests. The
%   results are compared to those of the equivalent parametric statistical
%   tests (i.e. paired t-tests) using ttest.m, and non-parametric
%   statistical tests (i.e. Wilcoxon signed-rank tests) using signrank.m.
%
%   See also PERMUTTEST TTEST SIGNRANK.
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
nobs = 30; nvar = 20;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = {'ttest','signrank'};
test_metric = {'t-value','z-value'};
ct_metric = {'mean','median'};
var_metric = {'SD','IQR'};

% Generate random data
rng(42);
x = randn(nobs,nvar);
y = randn(nobs,nvar);
y(:,1:round(nvar/2)) = y(:,1:round(nvar/2))-1;
if insert_nan
    x(1,1) = NaN;
end

disp('Mean elapsed time:')

for t = 1:numel(type)

    disp(type{t})

    toc1 = zeros(numel(tail),1);
    toc2 = zeros(numel(tail),1);
    toc3 = zeros(numel(tail),1);

    f1 = figure('Name',['Paired ',type{t},': ',ct_metric{t},' & CIs'],...
        'NumberTitle','off');
    set(gcf,'color','w')
    f2 = figure('Name',['Paired ',type{t},': p-values'],...
        'NumberTitle','off');
    set(gcf,'color','w')

    for i = 1:numel(tail)

        % Parametric test
        tic
        switch type{t}
            case 'ttest'
                [~,p1,ci1,stats1] = ttest(x,y,'tail',tail{i});
            case 'signrank'
                p1 = zeros(nvar,1);
                zval1 = zeros(nvar,1);
                iqr1 = zeros(nvar,1);
                for j = 1:nvar
                    [p1(j),~,stats1] = signrank(x(:,j),y(:,j),...
                        'tail',tail{i});
                    zval1(j) = stats1.zval;
                    iqr1(j) = iqr(x(:,j)-y(:,j));
                end
                ci1 = nan(2,nvar);
        end
        toc1(i) = toc;

        % Permutation test (uncorrected)
        tic
        [~,p2,ci2,stats2] = permuttest(x,y,'tail',tail{i},'correct',0,...
            'verbose',0,'type',type{t});
        toc2(i) = toc;

        % Permutation test (max-corrected)
        tic
        [~,p3,ci3,stats3] = permuttest(x,y,'tail',tail{i},'correct',1,...
            'verbose',0,'type',type{t});
        toc3(i) = toc;

        % Plot central tendancy & CIs
        figure(f1)
        switch type{t}
            case 'ttest'
                ct = stats2.mean;
            case 'signrank'
                ct = stats2.median;
        end
        subplot(3,2,i+i-1), hold on
        plot(xaxis,ct,'LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),ct(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2<=alpha),ct(p2<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim([-3,3]), box on, grid on
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([label{i},'-tailed'])
        if i == 1
            switch type{t}
                case 'ttest'
                    legend(ct_metric{t},'95% CI (param.)','','95% CI (perm.)',...
                        'Location','best')
                case 'signrank'
                    legend(ct_metric{t},'Location','best')
            end
        end
        subplot(3,2,i+i), hold on
        plot(xaxis,ct,'LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci3,'--r')
        plot(xaxis(p1<=alpha),ct(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p3<=alpha),ct(p3<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim([-3,3]), box on, grid on
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
    figure('Name',['Paired ',type{t},': descriptive statistics'],...
        'NumberTitle','off');
    set(gcf,'color','w')
    switch type{t}
        case 'ttest'
            stat1 = stats1.tstat;
            stat2 = stats2.tstat;
            var1 = stats1.sd;
            var2 = stats2.sd;
        case 'signrank'
            stat1 = zval1;
            stat2 = stats2.zstat;
            var1 = iqr1;
            var2 = stats2.iqr;
    end
    subplot(2,2,1), hold on
    plot(xaxis,stat1,xaxis,stat2,'--','LineWidth',2)
    xlim([0,nvar+1]), ylim([-7,7]), box on, grid on
    title('Test Statistic')
    xlabel('variable')
    ylabel(test_metric{t})
    legend('param.','perm.','Location','best')
    subplot(2,2,2), hold on
    plot(xaxis,var1,xaxis,var2,'--','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,3]), box on, grid on
    title('Variance')
    xlabel('variable')
    ylabel(var_metric{t})

    fprintf('Parametric (uncorrect): %.1f ms\n',mean(toc1)*1e3)
    fprintf('Permutation (uncorrect): %.1f ms\n',mean(toc2)*1e3)
    fprintf('Permutation (max-corr.): %.1f ms\n',mean(toc3)*1e3)

end