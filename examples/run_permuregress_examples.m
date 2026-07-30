function run_permuregress_examples
%RUN_PERMUREGRESS_EXAMPLES  Run permutation multiple linear regression examples.
%   Generates random multivariate data for a predictor sample X and 
%   dependent variables Y. Y contains 20 variables, each with a correlation 
%   of 0 to X, except for the first 5 variables which have a positive 
%   relationship, and the second 5 variables which have a negative 
%   relationship. Each variable has 30 observations. 
%
%   Permutation-based multiple linear regression is performed using both 
%   the Freedman-Lane and Manly methods for two-tailed, right-tailed, and 
%   left-tailed tests. The results are compared to those of the equivalent 
%   parametric statistical tests (i.e. Ordinary Least Squares regression 
%   and Student's t-distribution).
%
%   See also PERMUREGRESS FITLM REGRESS.
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
nobs = 30; nvar = 20;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = {'freedmanlane','manly'};

% Generate random data
rng(42);
x = randn(nobs,1);
y = randn(nobs,nvar);
y(:,1:5) = y(:,1:5)+x*0.8;
y(:,6:10) = y(:,6:10)-x*0.8;

% Parametric multiple linear regression (Ordinary Least Squares)
xparam = [ones(nobs,1),x];
bparam = xparam\y;
df = nobs-2;
mse = sum((y-xparam*bparam).^2)/df;
separam = sqrt(diag(inv(xparam'*xparam)).*mse);
b1 = bparam(2,:);
se1 = separam(2,:);
t1 = b1./se1;

for t = 1:numel(type)
    
    % Plot parametric & permutation CIs
    figure('Name',['Regression (',type{t},'): coefficients & CIs'],...
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
        [b2,p2,ci_perm] = permuregress(x,y,'tail',tail{i},...
            'type',type{t},'correct',0);
        b2 = b2(2,:); p2_uncorr = p2(2,:); ci2 = squeeze(ci_perm(:,:,2));
        subplot(3,2,i+i-1), hold on
        plot(xaxis,b2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2_uncorr<=alpha),b2(p2_uncorr<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim([-2,2]), box on, grid on
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
        [~,p2,ci_perm] = permuregress(x,y,'tail',tail{i},'type',type{t},...
            'correct',1);
        p2_corr = p2(2,:); ci2 = squeeze(ci_perm(:,:,2));
        subplot(3,2,i+i), hold on
        plot(xaxis,b2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        plot(xaxis(p1<=alpha),b1(p1<=alpha),'ok','LineWidth',2)
        plot(xaxis(p2_corr<=alpha),b2(p2_corr<=alpha),'xr','LineWidth',2)
        xlim([0,nvar+1]), ylim([-2,2]), box on, grid on
        if i == 1
            title('Max-corrected')
        elseif i == 3
            xlabel('variable')
        end
    end
    
    % Plot parametric & permutation p-values
    figure('Name',['Regression (',type{t},'): p-values'],...
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
        [~,p2] = permuregress(x,y,'tail',tail{i},'type',type{t},...
            'correct',0);
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
        [~,p2] = permuregress(x,y,'tail',tail{i},'type',type{t},...
            'correct',1);
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