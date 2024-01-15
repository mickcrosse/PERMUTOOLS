function run_paper_examples

% Generate random data
rng(42);
x = randn(30,20);
y = randn(30,20);
y(:,1:10) = y(:,1:10)-1;

% Run MATLAB's two-sample parametric t-test
[h2,p1,ci1,stats1] = ttest2(x,y);

% Run PERMUTOOLS' two-sample permutation t-test
[t2,p2,ci2,stats2] = permuttest2(x,y);

% Run MATLAB's parametric effect size measure
d3 = zeros(1,20); 
ci3 = zeros(2,20);
for j = 1:20
    stats3 = meanEffectSize(x(:,j),y(:,j),'effect','cohen','paired',0);
    d3(j) = stats3.Effect;
    ci3(:,j) = stats3.ConfidenceIntervals';
end

% Run PERMUTOOLS' bootstrapped effect size measure
[d4,ci4,stats4] = booteffectsize(x,y,'effect','cohen','paired',0);

% Set up figure
figure('Name','Two-sample Permutation Test','NumberTitle','off')
set(gcf,'color','w')
xaxis = 1:20; 
alpha = 0.05;

% Plot test statistics & CIs
subplot(1,3,1), hold on
plot(xaxis,stats2.mu,'LineWidth',2.5)
plot(xaxis,ci1,'k',xaxis,ci2,'--r','LineWidth',1)
plot(xaxis(p1<=alpha),stats2.mu(p1<=alpha),'ok','LineWidth',2)
plot(xaxis(p2<=alpha),stats2.mu(p2<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([-3,3]), box on, grid on
title('Test Statistic'), xlabel('variable'), ylabel('X−Y')
legend('mean difference','95% CI (param.)','','95% CI (perm.)','Location','southwest')

% Plot p-values
subplot(1,3,2), hold on
plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',1.5)
ylim([0,1]), xlim([0,21]), box on, grid on
title('{\itP}-value'), xlabel('variable'), ylabel('probability')
legend('{\itp}-value (param.)','{\itp}-value (perm.)','Location','northwest')

% Plot effect sizes & CIs
subplot(1,3,3), hold on
plot(xaxis,d4,'LineWidth',2.5)
plot(xaxis,ci3,'k',xaxis,ci4,'--r','LineWidth',1)
ylim([-1,4]), xlim([0,21]), box on, grid on
title('Effect Size'), xlabel('variable'), ylabel('standardised effect size')
legend('Hedges''{\itg}','95% CI (param.)','','95% CI (boot.)')