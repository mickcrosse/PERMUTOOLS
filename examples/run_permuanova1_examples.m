function run_permuanova1_examples
%RUN_PERMUANOVA1_EXAMPLES  Run permutation one-way ANOVA examples.
%   Generates random data for 5 "independent" groups along the columns of
%   X. Each group has a mean value of 0, except for the first group which
%   has a mean value of 1. Each group has 30 observations. The group labels
%   are randomly sampled with replacement 20 times to generate 20 different
%   scenarios. One-way permutation-based ANOVAs are performed on each of
%   the datasets. The results are compared to those of the equivalent
%   parametric statistical test (i.e. one-way ANOVA) using anova1.m.
%
%   See also PERMUANOVA1 ANOVA1.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Generate random data
rng(42);
x = randn(30,5);
x(:,1) = x(:,1)+1;
xaxis = 1:20; alpha = 0.05;
group_labels = 1:5;

% Compute ANOVA
p1 = zeros(1,20);
f2 = zeros(1,20);
p2 = zeros(1,20);
ci2 = zeros(2,20);
s = RandStream('mlfg6331_64');
for i = 1:20
    group = datasample(s,group_labels,numel(group_labels),'Replace',true);
    [p1(i)] = anova1(x,group,'off');
    [f2(i),p2(i),ci2(:,i)] = permuanova1(x,group);
end

% Set up figure
figure('Name','One-way permutation-based ANOVA','NumberTitle','off')
set(gcf,'color','w')

% Plot F-statistic & CIs
subplot(2,2,1), hold on
plot(xaxis,f2,'LineWidth',3)
plot(xaxis,ci2,'k')
plot(xaxis(p1<=alpha),f2(p1<=alpha),'ok','LineWidth',2)
plot(xaxis(p2<=alpha),f2(p2<=alpha),'xr','LineWidth',2)
xlim([0,21]), ylim([0,15]), box on, grid on
title('Test Statistic'), xlabel('permutation'), ylabel('{\itF}-value')
legend('{\itF}-statistic','95% CI (perm.)')

% Plot p-values
subplot(2,2,2), hold on
plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
xlim([0,21]), ylim([0,1]), box on, grid on
title('{\itP}-values'), xlabel('variable'), ylabel('probability')
legend('{\itp}-value (param.)','{\itp}-value (perm.)')