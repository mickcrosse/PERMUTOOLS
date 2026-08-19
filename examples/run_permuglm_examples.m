function run_permuglm_examples
%RUN_PERMUGLM_EXAMPLES  Run permutation generalized linear model examples.
%   Generates random multivariate binary data for a predictor sample X and
%   dependent variables Y. Y contains 20 variables, each with a
%   relationship of 0 to X, except for the first 5 variables which have a
%   positive logistic relationship, and the second 5 variables which have a
%   negative logistic relationship. Each variable has 60 observations.
%
%   Permutation-based logistic regression (binomial GLM) is performed using
%   the Freedman-Lane, Manly, and Rank methods for two-tailed,
%   right-tailed, and left-tailed tests. The results are compared to those
%   of the equivalent parametric statistical test (generalized linear
%   model) using MATLAB's fitglm.m.
%
%   See also PERMUGLM FITGLM FITLM.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS
%
%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%
%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

close all; clc;

% Set up experiment
nobs = 60; nvar = 20;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = {'freedmanlane','manly','rank'};

% Generate random binary data using a logistic link
rng(42);
x = randn(nobs,1);
y = zeros(nobs,nvar);
% True underlying betas for the logistic function
beta_true = [repmat(1.5,1,5), repmat(-1.5,1,5), zeros(1,10)];
for v = 1:nvar
    % Calculate probabilities using logistic function
    prob = 1 ./ (1 + exp(-(0.5 + beta_true(v)*x)));
    % Generate binary outcomes based on probabilities
    y(:,v) = rand(nobs,1) < prob;
end

for t = 1:numel(type)

    disp(type{t})

    toc1 = zeros(numel(tail),1);
    toc2 = zeros(numel(tail),1);
    toc3 = zeros(numel(tail),1);

    % Set up dynamic y-limits (Rank coefficients scale much higher)
    if strcmp(type{t}, 'rank')
        y_lims = [-25, 25];
    else
        y_lims = [-4, 4];
    end

    f1 = figure('Name',['GLM (',type{t},'): coefficients & CIs'],...
        'NumberTitle','off');
    set(gcf,'color','w')
    f2 = figure('Name',['GLM (',type{t},'): p-values'],...
        'NumberTitle','off');
    set(gcf,'color','w')

    for i = 1:numel(tail)

        % Parametric test
        tic
        switch type{t}
            case 'rank'
                yparam = zeros(nobs,nvar);
                for v = 1:nvar
                    yparam(:,v) = tiedrank(y(:,v));
                end
                dist = 'normal';
            otherwise
                yparam = y;
                dist = 'binomial';
        end
        b1 = zeros(1,nvar);
        p1 = zeros(1,nvar);
        ci1 = zeros(2,nvar);
        tstat = zeros(1,nvar);
        res = zeros(1,nvar);
        r2 = zeros(1,nvar);
        se = zeros(1,nvar);
        df = nobs-2;
        for v = 1:nvar
            if strcmp(dist,'binomial')
                mdl = fitglm(x, yparam(:,v),'Distribution',dist,...
                    'DispersionFlag',false);
            else
                mdl = fitglm(x, yparam(:,v),'Distribution',dist);
            end
            b1(v) = mdl.Coefficients.Estimate(2);
            p1(v) = mdl.Coefficients.pValue(2);
            r2(v) = mdl.Rsquared.Ordinary;
            tstat(v) = mdl.Coefficients.tStat(2);
            res(v) = mdl.Residuals.Raw(2);
            se(v) = mdl.Coefficients.SE(2);
            switch tail{i}
                case 'right'
                    if b1(v) > 0
                        p1(v) = p1(v)/2;
                    else
                        p1(v) = 1-(p1(v)/2);
                    end
                case 'left'
                    if b1(v) < 0
                        p1(v) = p1(v)/2;
                    else
                        p1(v) = 1-(p1(v)/2);
                    end
            end
        end
        switch tail{i}
            case 'both'
                crit = tinv(1-alpha/2,df);
                ci1 = [b1-crit.*se;b1+crit.*se];
            case 'right'
                crit = tinv(1-alpha,df);
                ci1 = [b1-crit.*se;inf(1,nvar)];
            case 'left'
                crit = tinv(1-alpha,df);
                ci1 = [-inf(1,nvar);b1+crit.*se];
        end
        toc1(i) = toc;

        % Permutation test (uncorrected)
        tic
        [b2,p2,ci2,stats2] = permuglm(x,y,'distribution',dist,'tail',tail{i},...
            'type',type{t},'correct',0);
        b2 = b2(2,:); p2 = p2(2,:); ci2 = squeeze(ci2(:,:,2));
        toc2(i) = toc;

        % Permutation test (max-corrected)
        tic
        [b3,p3,ci3,stats3] = permuglm(x,y,'distribution',dist,'tail',tail{i},...
            'type',type{t},'correct',1);
        b3 = b3(2,:); p3 = p3(2,:); ci3 = squeeze(ci3(:,:,2));
        toc3(i) = toc;

        % Plot statistic & CIs
        figure(f1)
        subplot(3,2,i+i-1), hold on
        plot(xaxis,b1,xaxis,b2,'--','LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2<=alpha),b2(p2<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim(y_lims), box on, grid on
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([label{i},'-tailed'])
        if i == 1
            legend('\beta coeff. (param.)','\beta coeff. (perm.)',...
                '95% CI (param.)','','95% CI (perm.)','Location','best')
        end
        subplot(3,2,i+i), hold on
        plot(xaxis,b1,xaxis,b3,'--','LineWidth',2)
        plot(xaxis,ci1,'k',xaxis,ci3,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p3<=alpha),b3(p3<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim(y_lims), box on, grid on
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
    subplot(2,2,1), hold on
    plot(xaxis,r2,xaxis,stats3.r2,'--','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    title('R^2 Statistic')
    ylabel('\itR^2')
    legend('param.','perm.','Location','best')
    subplot(2,2,2), hold on
    plot(xaxis,tstat,xaxis,stats3.tstat(2,:),'--','LineWidth',2)
    xlim([0,nvar+1]), box on, grid on
    title('t Statistic')
    ylabel('\itt')
    subplot(2,2,3), hold on
    plot(xaxis,res,xaxis,stats3.res(2,:),'--','LineWidth',2)
    xlim([0,nvar+1]), box on, grid on
    title('Residuals')
    xlabel('variable')
    ylabel('\ite')
    subplot(2,2,4), hold on
    plot(xaxis,se,xaxis,stats3.se(2,:),'--','LineWidth',2)
    xlim([0,nvar+1]), box on, grid on
    title('Error Variance')
    xlabel('variable')
    ylabel('SE')

    fprintf('Parametric (uncorrect): %.1f ms\n',mean(toc1)*1e3)
    fprintf('Permutation (uncorrect): %.1f ms\n',mean(toc2)*1e3)
    fprintf('Permutation (max-corr.): %.1f ms\n',mean(toc3)*1e3)

end