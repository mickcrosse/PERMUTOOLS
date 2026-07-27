function run_permuanova1_examples
%RUN_PERMUANOVA1_EXAMPLES  Run one-way permutation ANOVA examples.
%   Generates random data for 5 groups along the columns of X. Each group
%   has a mean value of 0, except for the first group which has a mean
%   value of 1. Each group has 30 observations. The group labels are
%   randomly sampled with replacement 20 times to generate 20 different
%   scenarios. One-way permutation-based ANOVAs are performed on each of
%   the datasets. The results are compared to those of the equivalent
%   parametric statistical test (i.e. one-way ANOVA) using anova1.m, and
%   non-parametric statistical test (i.e. Kruskal-Wallis test) using
%   kruskalwallis.m.
%
%   See also PERMUANOVA1 ANOVA1 KRUSKALWALLIS.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Check version
info = ver;
isoctave = any(ismember({info.Name},'Octave'));

% Generate random data
rng(42);
nobs = 30; nvar = 5; nperm = 10;
x = randn(nobs,nvar);
x(:,1) = x(:,1)+1;
xaxis = 1:nperm; alpha = 0.05;
group_labels = 1:5;
type = {'anova1','kruskalwallis'};

% Compute ANOVA
p1 = zeros(1,nperm);
f2 = zeros(1,nperm);
p2 = zeros(1,nperm);
ci2 = zeros(2,nperm);
if ~isoctave
    s = RandStream('mlfg6331_64');
end

for t = 1:numel(type)

    for i = 1:nperm
        if isoctave
            group = datasample(group_labels,numel(group_labels),...
                'Replace',true);
        else
            group = datasample(s,group_labels,numel(group_labels),...
                'Replace',true);
        end
        switch type{t}
            case 'anova1'
                p1(i) = anova1(x,group,'off');
            case 'kruskalwallis'
                p1(i) = kruskalwallis(x,group,'off');
        end
        [f2(i),p2(i),ci2(:,i)] = permuanova1(x,group,'type',type{t});
    end

    % Set up figure
    figure('Name',['One-way ',type{t},' test'],'NumberTitle','off')
    set(gcf,'color','w')

    % Plot F-statistic & CIs
    subplot(2,2,1), hold on
    plot(xaxis,f2,'LineWidth',3)
    plot(xaxis,ci2,'k')
    plot(xaxis(p1<=alpha),f2(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),f2(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nperm+1]), ylim([0,15]), box on, grid on
    title('Test Statistic'), xlabel('permutation'), ylabel('{\itF}-value')
    legend('{\itF}-statistic','95% CI (perm.)')

    % Plot p-values
    subplot(2,2,2), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,nperm+1]), ylim([0,1]), box on, grid on
    title('{\itP}-values'), xlabel('permutation'), ylabel('probability')
    legend('{\itp}-value (param.)','{\itp}-value (perm.)')

end
