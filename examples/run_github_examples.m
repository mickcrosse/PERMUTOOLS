function run_github_examples

close all; clear; clc;

% Generate random data
rng(42);
x = randn(30,20);
y = randn(30,20);

% Make the first 10 variables of Y have a mean of -1
y(:,1:10) = y(:,1:10)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run MATLAB's two-sample variance test (F-test)
[h,p,ci,stats] = vartest2(x,y);

% Run PERMUTOOLS' two-sample permutation variance test (uncorrected)
[hu,pu,ciu,statsu] = permuvartest2(x,y,'correct',0);

% Run PERMUTOOLS' two-sample permutation variance test (max-corrected)
[hc,pc,cic,statsc] = permuvartest2(x,y,'correct',1);

% Get F-statistic
f = statsc.fstat;
xaxis = 1:size(f,2);

% Plot parametric & uncorrected permutation CIs
figure('Name','Unpaired Tests: F-statistic & CIs','NumberTitle','off')
set(gcf,'color','w')
subplot(2,2,1), hold on
plot(xaxis,f,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,ciu,'--r','LineWidth',1)
plot(xaxis(logical(h)),f(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hu)),f(logical(hu)),'xr','LineWidth',2)
xlim([0,21]), ylim([0,6]), box on, grid on
title('Uncorrected'), ylabel('{\itF}-value')
legend('{\itF}-statistic','parametric CI','','permutation CI','Location','northeast')

% Plot parametric & corrected permutation CIs
subplot(2,2,2), hold on
plot(xaxis,f,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,cic,'--r','LineWidth',1)
plot(xaxis(logical(h)),f(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hc)),f(logical(hc)),'xr','LineWidth',2)
xlim([0,21]), ylim([0,6]), box on, grid on
title('Max-corrected')

% Plot parametric & uncorrected permutation p-values
subplot(2,2,3), hold on
plot(xaxis,p,'k',xaxis,pu,'--r','LineWidth',2)
plot(xaxis(logical(h)),p(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hu)),pu(logical(hu)),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable'), ylabel('probability')
legend('parametric {\itp}','permutation {\itp}','Location','northeast')

% Plot parametric & corrected permutation p-values
subplot(2,2,4), hold on
plot(xaxis,p,'k',xaxis,pc,'--r','LineWidth',2)
plot(xaxis(logical(h)),p(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hc)),pc(logical(hc)),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run MATLAB's two-sample t-test
[h,p,ci,stats] = ttest2(x,y);

% Run PERMUTOOLS' two-sample permutation test (uncorrected)
[hu,pu,ciu,statsu] = permuttest2(x,y,'correct',0);

% Run PERMUTOOLS' two-sample permutation test (max-corrected)
[hc,pc,cic,statsc] = permuttest2(x,y,'correct',1);

% Compute the mean difference
md = mean(x)-mean(y);
xaxis = 1:size(md,2);

% Plot parametric & uncorrected permutation CIs
figure, set(gcf,'color','w')
subplot(2,2,1), hold on
plot(xaxis,md,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,ciu,'--r','LineWidth',1)
plot(xaxis(logical(h)),md(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hu)),md(logical(hu)),'xr','LineWidth',2)
xlim([0,21]), ylim([-3,3]), box on, grid on
title('Uncorrected'), ylabel('X−Y')
legend('mean difference','parametric CI','','permutation CI','Location','southwest')

% Plot parametric & corrected permutation CIs
subplot(2,2,2), hold on
plot(xaxis,md,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,cic,'--r','LineWidth',1)
plot(xaxis(logical(h)),md(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hc)),md(logical(hc)),'xr','LineWidth',2)
xlim([0,21]), ylim([-3,3]), box on, grid on
title('Max-corrected')

% Plot parametric & uncorrected permutation p-values
subplot(2,2,3), hold on
plot(xaxis,p,'k',xaxis,pu,'--r','LineWidth',2)
plot(xaxis(logical(h)),p(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hu)),pu(logical(hu)),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable'), ylabel('probability')
legend('parametric {\itp}','permutation {\itp}','Location','southwest')

% Plot parametric & corrected permutation p-values
subplot(2,2,4), hold on
plot(xaxis,p,'k',xaxis,pc,'--r','LineWidth',2)
plot(xaxis(logical(h)),p(logical(h)),'ok','LineWidth',2)
plot(xaxis(logical(hc)),pc(logical(hc)),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run MATLAB's effect size measure
d = zeros(1,20); ci = zeros(2,20);
for j = 1:20
    stats1 = meanEffectSize(x(:,j),y(:,j),'effect','cohen','paired',0);
    d(j) = stats1.Effect;
    ci(:,j) = stats1.ConfidenceIntervals';
end

% Run PERMUTOOLS' effect size measure (uncorrected)
[du,ciu] = booteffectsize(x,y,'effect','cohen','paired',0,'correct',0);

% Run PERMUTOOLS' effect size measure (bias-corrected)
[dc,cic] = booteffectsize(x,y,'effect','cohen','paired',0,'correct',1);

% Plot parametric & uncorrected bootstrapped measures
figure, set(gcf,'color','w')
subplot(2,2,1), hold on
plot(xaxis,du,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,ciu,'--r','LineWidth',1)
xlim([0,21]), ylim([-2,6]), box on, grid on
title('Uncorrected'), xlabel('variable'), ylabel('effect size')
legend('Cohen''s {\itd}','parametric CI','','boostrapped CI')

% Plot parametric & corrected bootstrapped measures
subplot(2,2,2), hold on
plot(xaxis,dc,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,cic,'--r','LineWidth',1)
xlim([0,21]), ylim([-2,6]), box on, grid on
title('Bias-corrected'), xlabel('variable')
legend('Hedges'' {\itg}','parametric CI','','boostrapped CI')
