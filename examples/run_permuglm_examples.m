function run_permuglm_examples
%RUN_PERMUGLM_EXAMPLES  Run permutation generalized linear model examples.
%   Generates random multivariate binary data for a predictor sample X and 
%   dependent variables Y. Y contains 20 variables, each with a relationship 
%   of 0 to X, except for the first 5 variables which have a positive 
%   logistic relationship, and the second 5 variables which have a negative 
%   logistic relationship. Each variable has 60 observations. 
%
%   Permutation-based logistic regression (binomial GLM) is performed using 
%   the Freedman-Lane, Manly, and Rank methods for two-tailed, right-tailed, 
%   and left-tailed tests. The results are compared to those of the 
%   equivalent parametric statistical tests using MATLAB's native fitglm.
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

% Set up experiment
nobs = 60; nvar = 20;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = {'freedmanlane','manly','rank'}; % Added rank type
dist = 'binomial'; % Base distribution

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
    
    % Set up dynamic y-limits (Rank coefficients scale much higher)
    if strcmp(type{t}, 'rank')
        y_lims = [-25, 25];
    else
        y_lims = [-4, 4];
    end
    
    % Prepare parametric data based on permutation type
    switch type{t}
        case 'rank'
            yparam = zeros(nobs,nvar);
            for v = 1:nvar
                yparam(:,v) = tiedrank(y(:,v));
            end
            distparam = 'normal'; % Rank forces normal distribution
        otherwise
            yparam = y;
            distparam = dist;
    end
    
    % Parametric generalized linear regression
    b1 = zeros(1,nvar);
    se1 = zeros(1,nvar);
    t1 = zeros(1,nvar);
    df = nobs-2;
    for v = 1:nvar
        if strcmp(distparam, 'binomial')
            mdl = fitglm(x, yparam(:,v), 'Distribution', distparam, 'DispersionFlag', false);
        else
            mdl = fitglm(x, yparam(:,v), 'Distribution', distparam);
        end
        b1(v) = mdl.Coefficients.Estimate(2);
        se1(v) = mdl.Coefficients.SE(2);
        t1(v) = mdl.Coefficients.tStat(2);
    end
    
    % Plot parametric & permutation CIs
    figure('Name',['GLM (',type{t},'): coefficients & CIs'],...
        'NumberTitle','off')
    set(gcf,'color','w')
    for i = 1:numel(tail)
        ci1 = zeros(2,nvar);
        switch tail{i}
            case 'both'
                p1 = 2*tcdf(-abs(t1),df);
                crit = tinv(1-alpha/2,df);
                ci1 = [b1-crit.*se1;b1+crit.*se1];
            case 'right'
                p1 = tcdf(-t1,df);
                crit = tinv(1-alpha,df);
                ci1 = [b1-crit.*se1;inf(1,nvar)];
            case 'left'
                p1 = tcdf(t1,df);
                crit = tinv(1-alpha,df);
                ci1 = [-inf(1,nvar);b1+crit.*se1];
        end
        
        [b2,p2,ci_perm] = permuglm(x,y,'distribution',dist,'tail',tail{i},...
            'type',type{t},'correct',0);
        b2 = b2(2,:); p2_uncorr = p2(2,:); ci2 = squeeze(ci_perm(:,:,2));
        
        subplot(3,2,i+i-1), hold on
        plot(xaxis,b2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2_uncorr<=alpha),b2(p2_uncorr<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim(y_lims), box on, grid on
        if i == 1
            title('Uncorrected')
        elseif i == 3
            xlabel('variable')
        end
        ylabel([label{i},'-tailed'])
        if i == 1
            legend('\beta coefficient','95% CI (param.)','',...
                '95% CI (perm.)','Location','best')
        end
        
        [~,p2,ci_perm] = permuglm(x,y,'distribution',dist,'tail',tail{i},...
            'type',type{t},'correct',1);
        p2_corr = p2(2,:); ci2 = squeeze(ci_perm(:,:,2));
        
        subplot(3,2,i+i), hold on
        plot(xaxis,b2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2_corr<=alpha),b2(p2_corr<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim(y_lims), box on, grid on
        if i == 1
            title('Max-corrected')
        elseif i == 3
            xlabel('variable')
        end
    end
    
    % Plot parametric & permutation p-values
    figure('Name',['GLM (',type{t},'): p-values'],...
        'NumberTitle','off')
    set(gcf,'color','w')
    for i = 1:numel(tail)
        switch tail{i}
            case 'both'
                p1 = 2*tcdf(-abs(t1),df);
            case 'right'
                p1 = tcdf(-t1,df);
            case 'left'
                p1 = tcdf(t1,df);
        end
        
        [~,p2] = permuglm(x,y,'distribution',dist,'tail',tail{i},...
            'type',type{t},'correct',0);
        p2_uncorr = p2(2,:);
        
        subplot(3,2,i+i-1), hold on
        plot(xaxis,p1,'k',xaxis,p2_uncorr,'--r','LineWidth',2)
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
        
        [~,p2] = permuglm(x,y,'distribution',dist,'tail',tail{i},...
            'type',type{t},'correct',1);
        p2_corr = p2(2,:);
        
        subplot(3,2,i+i), hold on
        plot(xaxis,p1,'k',xaxis,p2_corr,'--r','LineWidth',2)
        xlim([0,nvar+1]), ylim([0,1]), box on, grid on
        if i == 1
            title('Max-corrected')
        elseif i == 3
            xlabel('variable')
        end
    end
end