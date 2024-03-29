function run_github_examples
%RUN_GITHUB_EXAMPLES  Run GitHub examples.
%   Generates random multivariate data for 2 samples X and Y and runs the
%   following examples for GitHub README:
%       1. Plots the permutation distribution and test statistic
%       2. Tests for differences in multivariate sample variances (F-test)
%       3. Tests for differences in multivariate sample means (t-test)
%       4. Computes the effect sizes between multivariate samples
%       5. Measures the correlations between multivariate samples
%
%   See also PERMUVARTEST2 PERMUTTEST2 BOOTEFFECTSIZE PERMUCORR.
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

% Make the first 10 variables of Y have a mean of -1
y(:,1:10) = y(:,1:10)-1;

%% 1. Permutation Distribution Example

% Run PERMUTOOLS' permutation correlation measure (uncorrected)
[~,pu,~,~,distu] = permuttest2(x,y,'correct',0,'seed',7);

% Run PERMUTOOLS' permutation correlation measure (max-corrected)
[tc,pc,~,~,distc] = permuttest2(x,y,'correct',1,'seed',7);

% Get test statistic for third variable
t = tc(3);

% Plot parametric & uncorrected permutation CIs
figure('Name','Unpaired tests based on the t-statistic','NumberTitle','off')
set(gcf,'renderer','Painters')
% set(gca,'FontName','Helvetica')
set(gcf,'color','w'), hold on
histogram(distu(:,3),100,'FaceAlpha',0.5)
h = histogram(distc,100,'FaceAlpha',0.5,'FaceAlpha', 0.5);
uistack(h, 'top')
plot([t,t],[0,500],'--k','LineWidth',2)
xlim([-5,5]), ylim([0,500]), box on, grid on
title('Permutation Distribution'), xlabel('{\itt}-statistic'), ylabel('frequency')
legend(['uncorrected ({\itp} = ',num2str(round(pu(3),3)),')'],...
    ['max-corrected ({\itp} = ',num2str(round(pc(3),2)),')'],...
    ['test statistic ({\itt} = ',num2str(round(t,2)),')'],...
    'Location','northwest')

%% 2. F-test Example

% Run MATLAB's two-sample parametric variance test (F-test)
[~,p,ci] = vartest2(x,y);

% Run PERMUTOOLS' two-sample permutation variance test (uncorrected)
[f,pu,ciu] = permuvartest2(x,y,'correct',0);

% Run PERMUTOOLS' two-sample permutation variance test (max-corrected)
[~,pc,cic] = permuvartest2(x,y,'correct',1);

% Set up figure
figure('Name','Two-sample permutation F-test','NumberTitle','off')
set(gcf,'color','w')
xaxis = 1:size(f,2);
alpha = 0.05;

% Plot parametric & uncorrected permutation CIs
subplot(2,2,1), hold on
plot(xaxis,f,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,ciu,'--r','LineWidth',1)
plot(xaxis(p<=alpha),f(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pu<=alpha),f(pu<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,6]), box on, grid on
title('Uncorrected'), ylabel('{\itF}-value')
legend('{\itF}-statistic','parametric CI','','permutation CI','Location','northeast')

% Plot parametric & corrected permutation CIs
subplot(2,2,2), hold on
plot(xaxis,f,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,cic,'--r','LineWidth',1)
plot(xaxis(p<=alpha),f(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pc<=alpha),f(pc<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,6]), box on, grid on
title('Max-corrected')

% Plot parametric & uncorrected permutation p-values
subplot(2,2,3), hold on
plot(xaxis,p,'k',xaxis,pu,'--r','LineWidth',2)
plot(xaxis(p<=alpha),p(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pu<=alpha),pu(pu<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable'), ylabel('probability')
legend('parametric {\itp}','permutation {\itp}','Location','northeast')

% Plot parametric & corrected permutation p-values
subplot(2,2,4), hold on
plot(xaxis,p,'k',xaxis,pc,'--r','LineWidth',2)
plot(xaxis(p<=alpha),p(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pc<=alpha),pc(pc<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable')

%% 3. t-test Example

% Run MATLAB's two-sample parametric t-test
[~,p,ci] = ttest2(x,y);

% Run PERMUTOOLS' two-sample permutation test (uncorrected)
[~,pu,ciu,stats] = permuttest2(x,y,'correct',0);

% Run PERMUTOOLS' two-sample permutation test (max-corrected)
[~,pc,cic] = permuttest2(x,y,'correct',1);

% Set up figure
figure('Name','Two-sample permutation t-test','NumberTitle','off')
set(gcf,'color','w')
xaxis = 1:size(f,2);
alpha = 0.05;

% Plot parametric & uncorrected permutation CIs
subplot(2,2,1), hold on
plot(xaxis,stats.mu,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,ciu,'--r','LineWidth',1)
plot(xaxis(p<=alpha),stats.mu(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pu<=alpha),stats.mu(pu<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([-3,3]), box on, grid on
title('Uncorrected'), ylabel('X−Y')
legend('mean difference','parametric CI','','permutation CI','Location','southwest')

% Plot parametric & corrected permutation CIs
subplot(2,2,2), hold on
plot(xaxis,stats.mu,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,cic,'--r','LineWidth',1)
plot(xaxis(p<=alpha),stats.mu(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pc<=alpha),stats.mu(pc<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([-3,3]), box on, grid on
title('Max-corrected')

% Plot parametric & uncorrected permutation p-values
subplot(2,2,3), hold on
plot(xaxis,p,'k',xaxis,pu,'--r','LineWidth',2)
plot(xaxis(p<=alpha),p(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pu<=alpha),pu(pu<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable'), ylabel('probability')
legend('parametric {\itp}','permutation {\itp}','Location','southwest')

% Plot parametric & corrected permutation p-values
subplot(2,2,4), hold on
plot(xaxis,p,'k',xaxis,pc,'--r','LineWidth',2)
plot(xaxis(p<=alpha),p(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pc<=alpha),pc(pc<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable')

%% 4. Effect Size Example

% Run MATLAB's parametric effect size measure
d = zeros(1,20); ci = zeros(2,20);
for j = 1:20
    stats1 = meanEffectSize(x(:,j),y(:,j),'effect','cohen','paired',0);
    d(j) = stats1.Effect;
    ci(:,j) = stats1.ConfidenceIntervals';
end

% Run PERMUTOOLS' bootstrapped effect size measure (uncorrected)
[du,ciu] = booteffectsize(x,y,'effect','cohen','paired',0,'correct',0);

% Run PERMUTOOLS' bootstrapped effect size measure (bias-corrected)
[dc,cic] = booteffectsize(x,y,'effect','cohen','paired',0,'correct',1);

% Plot parametric & uncorrected bootstrapped measures
figure('Name','Effect size measures based on Cohen''s d','NumberTitle','off')
set(gcf,'color','w')
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

%% 5. Correlation Example

% Generate random data
rng(42);
x = randn(30,20);
y = randn(30,20);

% Make the some variables positively and negatively correlated
y(:,1:5) = y(:,1:5)+x(:,1:5)/2;
y(:,6:10) = y(:,6:10)-x(:,6:10);
xaxis = 1:20; alpha = 0.05;

% Run MATLAB's parametric correlation measure
[r,p] = corr(x,y);
r = diag(r);
p = diag(p);
ci = zeros(2,20);
for j = 1:20
    [~,~,clwr,cupr] = corrcoef(x(:,j),y(:,j));
    ci(:,j) = [clwr(2);cupr(2)];
end

% Run PERMUTOOLS' permutation correlation measure (uncorrected)
[ru,pu,ciu] = permucorr(x,y,'correct',0);

% Run PERMUTOOLS' permutation correlation measure (max-corrected)
[rc,pc,cic] = permucorr(x,y,'correct',1);

% Plot parametric & uncorrected permutation CIs
figure('Name','Correlation measures based on Pearson''s r','NumberTitle','off')
set(gcf,'color','w')
subplot(2,2,1), hold on
plot(xaxis,ru,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,ciu,'--r','LineWidth',1)
plot(xaxis(p<=alpha),r(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pu<=alpha),ru(pu<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([-1,1]), box on, grid on
title('Uncorrected'), ylabel('correlation')
legend('Pearson''s {\itr}','parametric CI','','permutation CI')

% Plot parametric & corrected permutation CIs
subplot(2,2,2), hold on
plot(xaxis,rc,'LineWidth',3)
plot(xaxis,ci,'k',xaxis,cic,'--r','LineWidth',1)
plot(xaxis(p<=alpha),r(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pc<=alpha),rc(pc<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([-1,1]), box on, grid on
title('Max-corrected')

% Plot parametric & uncorrected permutation p-values
subplot(2,2,3), hold on
plot(xaxis,p,'k',xaxis,pu,'--r','LineWidth',2)
plot(xaxis(p<=alpha),p(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pu<=alpha),pu(pu<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable'), ylabel('probability')
legend('parametric {\itp}','permutation {\itp}')

% Plot parametric & corrected permutation p-values
subplot(2,2,4), hold on
plot(xaxis,p,'k',xaxis,pc,'--r','LineWidth',2)
plot(xaxis(p<=alpha),p(p<=alpha),'ok','LineWidth',2)
plot(xaxis(pc<=alpha),pc(pc<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
xlabel('variable')
