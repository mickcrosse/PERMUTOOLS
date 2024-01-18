function run_paper_examples
%RUN_PAPER_EXAMPLES  Run arXiv paper examples.
%   Generates random multivariate data for 2 "independent" samples X and Y.
%   Each sample has 20 variables, each with a mean value of 0, except for
%   the first 10 variables of Y which have a mean value of -1. Each
%   variable has 30 observations. Two-sample permutation tests based on the
%   t-statistic are performed between the corresponding variables of each
%   sample for a two-tailed test, assuming the samples of equal variances.
%   The results are compared to those of the equivalent parametric
%   statistical tests (i.e. two-sample t-tests) using ttest2.m. Effect
%   sizes based on Hedges's g with bootstrapped confidence intervals (CIs)
%   are measured between the corresponding variables of each independent
%   sample, assuming equal variances. The results are compared to those of
%   the equivalent parametric measures (i.e. CIs using the Student's
%   t-distribution) using meanEffectSize.m.
%
%   See also PERMUTTEST2 BOOTEFFECTSIZE TTEST2 MEANEFFECTSIZE.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate random data
rng(42);
x = randn(30,20);
y = randn(30,20);
y(:,1:10) = y(:,1:10)-1;

% Run MATLAB's two-sample parametric t-test
[h2,p1,ci1,stats1] = ttest2(x,y);

% Run PERMUTOOLS' two-sample permutation t-test
[t2,p2,ci2,stats2] = permuttest2(x,y);

% Run MATLAB's parametric effect size analysis
d3 = zeros(1,20);
ci3 = zeros(2,20);
for j = 1:20
    stats3 = meanEffectSize(x(:,j),y(:,j),'effect','cohen','paired',0);
    d3(j) = stats3.Effect;
    ci3(:,j) = stats3.ConfidenceIntervals';
end

% Run PERMUTOOLS' bootstrapped effect size analysis
[d4,ci4,stats4] = booteffectsize(x,y,'effect','cohen','paired',0);

% Set up figure
figure('Name','Permutation Tests & Effect Size Analysis','NumberTitle','off')
set(gcf,'color','w')
xaxis = 1:20;

% Plot test statistic & CI
subplot(1,3,1), hold on
plot(xaxis,stats2.mu,'LineWidth',2.5)
plot(xaxis,ci1,'k',xaxis,ci2,'--r','LineWidth',1)
xlim([0,21]), ylim([-3,3]), box on, grid on
title('Test Statistic'), xlabel('variable'), ylabel('X−Y')
legend('mean difference','95% CI (param.)','','95% CI (perm.)','Location','southwest')

% Plot p-value
subplot(1,3,2), hold on
plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',1.5)
ylim([0,1]), xlim([0,21]), box on, grid on
title('{\itP}-value'), xlabel('variable'), ylabel('probability')
legend('{\itp}-value (param.)','{\itp}-value (perm.)','Location','northwest')

% Plot effect size & CI
subplot(1,3,3), hold on
plot(xaxis,d4,'LineWidth',2.5)
plot(xaxis,ci3,'k',xaxis,ci4,'--r','LineWidth',1)
ylim([-1,4]), xlim([0,21]), box on, grid on
title('Effect Size'), xlabel('variable'), ylabel('standardised effect size')
legend('Hedges''{\itg}','95% CI (param.)','','95% CI (boot.)')