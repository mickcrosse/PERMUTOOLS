function run_booteffectsize_examples
%RUN_BOOTEFFECTSIZE_EXAMPLES  Run bootstrapped effect size examples.
%   Generates random multivariate data for 2 samples X and Y. Each sample
%   has 20 variables, each with a mean value of 0, except for the first 10
%   variables of Y which have a mean value of -1. Each variable has 30
%   observations. Effect sizes based on Cohen's d and Glass' Δ with
%   bootstrapped confidence intervals (CIs) are measured between the
%   corresponding variables of each sample for independent and dependent
%   samples, as well as samples with equal and unequal variances. The
%   results are compared to those of the equivalent parametric measures
%   (i.e. CIs using the Student's t-distribution) using meanEffectSize.m.
%
%   Note: The parametric effect sizes and CIs generated by meanEffectSize.m
%   for the uncorrected plots are actually bias-corrected (as there is no
%   option to specify bias-corrected in meanEffectSize.m) and are thus
%   slightly overestimated.
%
%   See also BOOTEFFECTSIZE MEANEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate random data
seed = 42;
rng(seed);
x = randn(30,20);
y = randn(30,20);
y(:,1:10) = y(:,1:10)-1;
xaxis = 1:20;
paired = [false,true];
samples = {'indep.','dep.'};
vartype = {'equal','unequal'};

figure('Name','Effect Size Analysis: effect sizes & CIs','NumberTitle','off')
set(gcf,'color','w')
k = 1;
for i = 1:numel(paired)

    for n = 1:numel(vartype)

        % Make variance unequal
        if n == 2
            yp = y*1.25;
        else
            yp = y;
        end

        % Set effect size measure
        if paired(i)==true && strcmp(vartype{n},'unequal')
            paired(i) = false;
            samples{i} = 'indep.';
            effect = 'Glass';
        else
            effect = 'Cohen';
        end

        % Plot parametric & uncorrected bootstrapped CIs
        d1 = zeros(1,20);
        ci1 = zeros(2,20);
        for j = 1:20
            stats1 = meanEffectSize(x(:,j),yp(:,j),'Effect',effect,...
                'Paired',paired(i),'VarianceType',vartype{n},...
                'ConfidenceIntervalType','exact');
            d1(j) = stats1.Effect;
            ci1(:,j) = stats1.ConfidenceIntervals';
        end
        [d2,ci2] = booteffectsize(x,yp,'paired',paired(i),...
            'vartype',vartype{n},'effect',effect,'correct',0);
        subplot(4,2,k), hold on
        plot(xaxis,d2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        xlim([0,21]), ylim([-2,5]), box on, grid on
        if i == 1 && n == 1
            title('Uncorrected')
        end
        if i == 2 && n == 2
            xlabel('variable')
        end
        ylabel({[samples{i},' samples'];['w/ ',vartype{n},' SD']})
        switch effect
            case 'Cohen'
                legend('Cohen''s {\itd}','95% CI (param.)','',...
                    '95% CI (boot.)')
            case 'Glass'
                legend('Glass'' {\itΔ}','95% CI (param.)','',...
                    '95% CI (boot.)')
        end
        [d2,ci2] = booteffectsize(x,yp,'paired',paired(i),...
            'vartype',vartype{n},'effect',effect,'correct',1);
        subplot(4,2,k+1), hold on
        plot(xaxis,d2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        xlim([0,21]), ylim([-2,5]), box on, grid on
        if i == 1 && n == 1
            title('Corrected')
        end
        if i == 2 && n == 2
            xlabel('variable')
        end
        switch effect
            case 'Cohen'
                legend('Hedges'' {\itg}','95% CI (param.)','',...
                    '95% CI (boot.)')
            case 'Glass'
                legend('Glass'' {\itΔ}','95% CI (param.)','',...
                    '95% CI (boot.)')
        end

        k = k+2;

    end

end