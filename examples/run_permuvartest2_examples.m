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

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate random data
rng(42);
nobs = 30; nvar = 20;
x = randn(nobs,nvar);
y = randn(nobs,nvar);
y(:,1:round(nvar/2)) = y(:,1:round(nvar/2))*2;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = {'ftest','squarerank'};

for t = 1:numel(type)

    % Plot parametric & permutation CIs
    figure('Name',[type{t},' test: F-statistic & CIs'],...
        'NumberTitle','off')
    set(gcf,'color','w')
    for i = 1:numel(tail)
        switch type{t}
            case 'ftest'
                [~,p1,ci1] = vartest2(x,y,'tail',tail{i});
            case 'squarerank'
                p1 = nan(nvar,1);
                ci1 = nan(nvar,2);
        end
        [f2,p2,ci2] = permuvartest2(x,y,'tail',tail{i},'correct',0,...
            'type',type{t});
        subplot(3,2,i+i-1), hold on
        plot(xaxis,f2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),f2(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2<=alpha),f2(p2<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim([0,6]), box on, grid on
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([label{i},'-tailed'])
        if i == 2
            legend('{\itF}-statistic','95% CI (param.)','','95% CI (perm.)')
        end
        [f2,p2,ci2] = permuvartest2(x,y,'tail',tail{i},'type',type{t},...
            'correct',1);
        subplot(3,2,i+i), hold on
        plot(xaxis,f2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),f2(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2<=alpha),f2(p2<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim([0,6]), box on, grid on
        if i == 1
            title('Max-corrected')
        elseif i == 3
            xlabel('variable')
        end
    end

    % Plot parametric & permutation p-values
    figure('Name',[type{t},' test: p-values'],'NumberTitle','off')
    set(gcf,'color','w')
    for i = 1:numel(tail)
        switch type{t}
            case 'ftest'
                [~,p1] = vartest2(x,y,'tail',tail{i});
            case 'squarerank'
                p1 = nan(nvar,1);
        end
        [~,p2] = permuvartest2(x,y,'tail',tail{i},'correct',0,...
            'type',type{t});
        subplot(3,2,i+i-1), hold on
        plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
        xlim([0,nvar+1]), ylim([0,1]), box on, grid on
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([label{i},'-tailed'])
        if i == 2
            legend('{\itp}-value (param.)','{\itp}-value (perm.)')
        end
        [~,p2] = permuvartest2(x,y,'tail',tail{i},'type',type{t},...
            'correct',1);
        subplot(3,2,i+i), hold on
        plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
        xlim([0,nvar+1]), ylim([0,1]), box on, grid on
        if i == 1
            title('Max-corrected')
        elseif i == 3
            xlabel('variable')
        end
    end

end