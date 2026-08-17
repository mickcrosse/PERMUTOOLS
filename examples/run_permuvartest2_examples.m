function run_permuvartest2_examples
%RUN_PERMUVARTEST2_EXAMPLES  Run two-sample permutation F-test examples.
%   Generates random multivariate data for 2 "independent" samples X and Y.
%   Each sample has 20 variables, each with a standard deviation of 1,
%   except for the first 10 variables of Y which have a standard deviation
%   of 2. Each variable has 30 observations. Two-sample permutation tests
%   based on the F-statistic are performed between the corresponding
%   variables of each sample for two-tailed, right-tailed and left-tailed
%   tests. The results are compared to those of the equivalent parametric
%   statistical test (two-sample F-test) using MATLAB's vartest2.m.
%
%   See also PERMUVARTEST2 VARTEST2.
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
nobs = 30; nvar = 20;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = {'ftest','squarerank'};

% Generate random data
rng(42);
x = randn(nobs,nvar);
y = randn(nobs,nvar);
y(:,1:round(nvar/2)) = y(:,1:round(nvar/2))*2;

for t = 1:numel(type)

    disp(type{t})

    toc1 = zeros(numel(tail),1);
    toc2 = zeros(numel(tail),1);
    toc3 = zeros(numel(tail),1);

    f1 = figure('Name',['Two-sample ',type{t},': statistic & CIs'],...
        'NumberTitle','off');
    set(gcf,'color','w')
    f2 = figure('Name',['Two-sample ',type{t},': p-values'],...
        'NumberTitle','off');
    set(gcf,'color','w')

    for i = 1:numel(tail)

        % Parametric test
        tic
        switch type{t}
            case 'ftest'
                [~,p1,ci1,stats1] = vartest2(x,y,'tail',tail{i});
            case 'squarerank'
                p1 = nan(nvar,1);
                ci1 = nan(nvar,2);
        end
        toc1(i) = toc;

        % Permutation test (uncorrected)
        tic
        [~,p2,ci2,stats2] = permuvartest2(x,y,'type',type{t},...
            'tail',tail{i},'correct',0);
        toc2(i) = toc;

        % Permutation test (max-corrected)
        tic
        [~,p3,ci3,stats3] = permuvartest2(x,y,'type',type{t},...
            'tail',tail{i},'correct',1);
        toc3(i) = toc;

        switch type{t}
            case 'ftest'
                stat1 = stats1.fstat;
                stat2 = stats2.fstat;
                stat3 = stats3.fstat;
                ylims = [0,6];
            case 'squarerank'
                stat1 = nan(nvar,1);
                stat2 = stats2.tstat;
                stat3 = stats3.tstat;
                ylims = [-6,6];
        end

        % Plot statistic & CIs
        figure(f1)
        subplot(3,2,i+i-1), hold on
        plot(xaxis,stat1,xaxis,stat2,'--','LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),stat1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2<=alpha),stat2(p2<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim(ylims), box on, grid on
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([label{i},'-tailed'])
        if i == 1
            legend('test stat. (param.)','test stat. (perm.)',...
                '95% CI (param.)','','95% CI (perm.)','Location','best')
        end
        subplot(3,2,i+i), hold on
        plot(xaxis,stat1,xaxis,stat3,'--','LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci3,'--r')
        plot(xaxis(p1<=alpha),stat1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p3<=alpha),stat3(p3<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim(ylims), box on, grid on
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

end