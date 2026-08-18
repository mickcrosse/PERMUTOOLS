function run_permuregress_examples
%RUN_PERMUREGRESS_EXAMPLES  Run permutation multiple linear regression examples.
%   Generates random multivariate data for a predictor sample X and
%   dependent variables Y. Y contains 20 variables, each with a correlation
%   of 0 to X, except for the first 5 variables which have a positive
%   relationship, and the second 5 variables which have a negative
%   relationship. Each variable has 30 observations. Permutation-based
%   multiple linear regression is performed using the Freedman-Lane, Manly,
%   and Rank methods for two-tailed, right-tailed, and left-tailed tests.
%   The results are compared to those of the equivalent parametric
%   statistical test (multiple linear regression via least squares) using 
%   MATLAB's regress.m.
%
%   See also PERMUREGRESS REGRESS.
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
type = {'freedmanlane','manly','rank'};

% Generate random data
rng(42);
x = randn(nobs,1);
y = randn(nobs,nvar);
y(:,1:5) = y(:,1:5)+x*0.8;
y(:,6:10) = y(:,6:10)-x*0.8;

for t = 1:numel(type)

    disp(type{t})

    toc1 = zeros(numel(tail),1);
    toc2 = zeros(numel(tail),1);
    toc3 = zeros(numel(tail),1);

    f1 = figure('Name',['Regression (',type{t},'): coefficients & CIs'],...
        'NumberTitle','off');
    set(gcf,'color','w')
    f2 = figure('Name',['Regression (',type{t},'): p-values'],...
        'NumberTitle','off');
    set(gcf,'color','w')

    for i = 1:numel(tail)

        % Parametric test
        tic
        switch type{t}
            case 'rank'
                yparam = tiedrank(y);
                xparam = [ones(nobs,1),tiedrank(x)];
            otherwise
                yparam = y;
                xparam = [ones(nobs,1),x];
        end
        b1 = zeros(1,nvar);
        p1 = zeros(1,nvar);
        ci1 = zeros(2,nvar);
        switch tail{i}
            case 'both'
                cialpha = alpha;
            case {'right','left'}
                cialpha = alpha*2;
        end
        for v = 1:nvar
            [b,bint,~,~,stats] = regress(yparam(:,v),xparam,cialpha);
            b1(v) = b(2);
            switch tail{i}
                case 'both'
                    p1(v) = stats(3);
                    ci1(:,v) = bint(2,:);
                case 'right'
                    if b(2) > 0
                        p1(v) = stats(3)/2;
                    else
                        p1(v) = 1-(stats(3)/2);
                    end
                    ci1(:,v) = [bint(2,1);Inf];
                case 'left'
                    if b(2) < 0
                        p1(v) = stats(3)/2;
                    else
                        p1(v) = 1-(stats(3)/2);
                    end
                    ci1(:,v) = [-Inf;bint(2,2)];
            end
        end
        toc1(i) = toc;

        % Permutation test (uncorrected)
        tic
        [b2,p2,ci2] = permuregress(x,y,'tail',tail{i},...
            'type',type{t},'correct',0);
        b2 = b2(2,:); p2 = p2(2,:); ci2 = squeeze(ci2(:,:,2));
        toc2(i) = toc;

        % Permutation test (max-corrected)
        tic
        [b3,p3,ci3] = permuregress(x,y,'tail',tail{i},'type',type{t},...
            'correct',1);
        b3 = b3(2,:); p3 = p3(2,:); ci3 = squeeze(ci3(:,:,2));
        toc3(i) = toc;

        % Plot statistic & CIs
        figure(f1)
        subplot(3,2,i+i-1), hold on
        plot(xaxis,b1,xaxis,b2,'--','LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2<=alpha),b2(p2<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim([-2,2]), box on, grid on
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([label{i},'-tailed'])
        if i == 1
            legend('\beta coef. (param.)','\beta coef. (perm.)',...
                '95% CI (param.)','','95% CI (perm.)','Location','best')
        end
        subplot(3,2,i+i), hold on
        plot(xaxis,b1,xaxis,b3,'--','LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci3,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p3<=alpha),b3(p3<=alpha),'xr','LineWidth',2)
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

    fprintf('Parametric (uncorrect): %.1f ms\n',mean(toc1)*1e3)
    fprintf('Permutation (uncorrect): %.1f ms\n',mean(toc2)*1e3)
    fprintf('Permutation (max-corr.): %.1f ms\n',mean(toc3)*1e3)

end