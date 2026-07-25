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
%   t-tests) using ttest2.m, and non-parametric statistical tests (i.e.
%   Mann-Whitney U / Wilcoxon rank-sum tests) using ranksum.m.
%
%   See also PERMUTTEST2 TTEST2 RANKSUM.
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
y(:,1:round(nvar/2)) = y(:,1:round(nvar/2))-1;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = {'mean','rank'};
vartype = {'equal','unequal'};

for t = 1:numel(type)

    for v = 1:numel(vartype)

        % Make variance unequal
        if v == 2
            y = y*1.25;
        end

        % Plot parametric & permutation CIs
        figure('Name',['Unpaired ',type{t},'-based test (',vartype{v},...
            ' SD): mean difference & CIs'],'NumberTitle','off')
        set(gcf,'color','w')
        for i = 1:numel(tail)
            switch type{t}
                case 'mean'
                    [~,p1,ci1] = ttest2(x,y,'tail',tail{i},'vartype',vartype{v});
                    ylims = [-3,3];
                case 'rank'
                    p1 = zeros(nvar,1);
                    for j = 1:nvar
                        p1(j) = ranksum(x(:,j),y(:,j),'tail',tail{i});
                    end
                    ci1 = nan(nvar,2);
                    ylims = [-40,40];
            end
            [~,p2,ci2,stats2] = permuttest2(x,y,'tail',tail{i},...
                'type',type{t},'vartype',vartype{v},'correct',0);
            subplot(3,2,i+i-1), hold on
            plot(xaxis,stats2.mu,'LineWidth',3)
            plot(xaxis,ci1,'k',xaxis,ci2,'--r')
            plot(xaxis(p1<=alpha),stats2.mu(p1<=alpha),'ok','LineWidth',2)
            plot(xaxis(p2<=alpha),stats2.mu(p2<=alpha),'xr','LineWidth',2)
            xlim([0,nvar+1]), ylim(ylims), box on, grid on
            if i == 1
                title('Uncorrected')
            elseif i == 3
                xlabel('variable')
            end
            ylabel([label{i},'-tailed'])
            if i == 2
                legend('mean difference','95% CI (param.)','','95% CI (perm.)')
            end
            [~,p2,ci2,stats2] = permuttest2(x,y,'tail',tail{i},...
                'type',type{t},'vartype',vartype{v},'correct',1);
            subplot(3,2,i+i), hold on
            plot(xaxis,stats2.mu,'LineWidth',3)
            plot(xaxis,ci1,'k',xaxis,ci2,'--r')
            plot(xaxis(p1<=alpha),stats2.mu(p1<=alpha),'ok','LineWidth',2)
            plot(xaxis(p2<=alpha),stats2.mu(p2<=alpha),'xr','LineWidth',2)
            xlim([0,nvar+1]), ylim(ylims), box on, grid on
            if i == 1
                title('Max-corrected')
            elseif i == 3
                xlabel('variable')
            end
        end

        % Plot parametric & permutation p-values
        figure('Name',['Unpaired ',type{t},'-based test (',vartype{v},...
            ' SD): p-values'],'NumberTitle','off')
        set(gcf,'color','w')
        for i = 1:numel(tail)
            switch type{t}
                case 'mean'
                    [~,p1] = ttest2(x,y,'tail',tail{i},'vartype',vartype{v});
                case 'rank'
                    p1 = zeros(nvar,1);
                    for j = 1:nvar
                        p1(j) = ranksum(x(:,j),y(:,j),'tail',tail{i});
                    end
            end
            [~,p2] = permuttest2(x,y,'tail',tail{i},'type',type{t},...
                'vartype',vartype{v},'correct',0);
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
            [~,p2] = permuttest2(x,y,'tail',tail{i},'type',type{t},...
                'vartype',vartype{v},'correct',1);
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

end