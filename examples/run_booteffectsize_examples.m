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
%   See also BOOTEFFECTSIZE MEANEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018 Mick Crosse <mickcrosse@gmail.com>
%   CNL, Albert Einstein College of Medicine, NY.

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

figure
k = 1;
for i = 1:numel(paired)

    for n = 1:numel(vartype)

        % Make variance unequal
        if n == 2
            yp = y*1.25;
        else
            yp = y;
        end

        % Set Cohen's d or Glass' Δ params
        if paired(i)==true && strcmp(vartype{n},'unequal')
            pairedi = false;
            samplesi = 'indep.';
            effect = 'Glass';
            symbol = 'Δ';
        else
            pairedi = paired(i);
            samplesi = samples{i};
            effect = 'Cohen';
            symbol = 'd';
        end

        % Plot parametric & bootstrapped CIs
        d1 = zeros(1,20);
        ci1 = zeros(2,20);
        for j = 1:20
            stats1 = meanEffectSize(x(:,j),yp(:,j),'Effect',effect,...
                'Paired',pairedi,'VarianceType',vartype{n},...
                'ConfidenceIntervalType','exact');
            d1(j) = stats1.Effect;
            ci1(:,j) = stats1.ConfidenceIntervals';
        end
        [d2,ci2] = booteffectsize(x,yp,'paired',pairedi,...
            'vartype',vartype{n},'effect',effect,'correct',false);
        subplot(4,2,k)
        hold on
        plot(xaxis,d2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        hold off
        ylim([-2,5])
        xlim([0,21])
        if i == 1 && n == 1
            title('Uncorrected')
        end
        if i == 2 && n == 2
            xlabel('variable')
        end
        ylabel({[samplesi,' samples'];['w/ ',vartype{n},' SD']})
        legend([effect,'''s ',symbol],'parametric CIs','','boostrapped CIs')
        [d2,ci2] = booteffectsize(x,yp,'paired',pairedi,...
            'vartype',vartype{n},'effect',effect,'correct',true);
        subplot(4,2,k+1)
        hold on
        plot(xaxis,d2,'LineWidth',3)
        plot(xaxis,ci1,'k',xaxis,ci2,'--r')
        hold off
        ylim([-2,5])
        xlim([0,21])
        if i == 1 && n == 1
            title('Corrected')
        end
        if i == 2 && n == 2
            xlabel('variable')
        end

        k = k+2;

    end

end