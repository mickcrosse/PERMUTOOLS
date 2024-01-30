function run_permuanova2_examples
%RUN_PERMUANOVA2_EXAMPLES  Run permutation two-way ANOVA examples.
%   Generates random data for 5 groups along the columns of X and 3 groups
%   along the rows of X, each with 2 replicates. Each group has a mean
%   value of 0, except for the first two columns which have a mean value of
%   1. The row labels are randomly sampled without replacement 20 times to
%   generate 20 different scenarios. Two-way permutation-based ANOVAs are
%   performed on each of the datasets. The results are compared to those of
%   the equivalent parametric statistical test (i.e. two-way ANOVAs) using
%   anova2.m.
%
%   See also PERMUANOVA2 ANOVA2.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2024 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Check version
info = ver;
isoctave = any(ismember({info.Name},'Octave'));

% Generate random data
rng(42);
x = randn(6,5);
x(:,1:2) = x(:,1:2)+1;
xaxis = 1:20;
alpha = 0.05;
reps = 2;
if reps > 1
    dim = 3;
else
    dim = 2;
end
ylabels = {'columns','rows','interaction'};

% Compute ANOVA
p1 = zeros(dim,20);
f2 = zeros(dim,20);
p2 = zeros(dim,20);
ci2 = zeros(2,dim,20);
if ~isoctave
    s = RandStream('mlfg6331_64');
end
for i = 1:20
    if isoctave
        idx = datasample(1:6,6,'Replace',false);
    else
        idx = datasample(s,1:6,6,'Replace',false);
    end
    [p1(:,i)] = anova2(x(idx,:),reps,'off');
    [f2(:,i),p2(:,i),ci2(:,:,i)] = permuanova2(x(idx,:),reps);
end

% Set up figure
figure('Name','Two-way permutation-based ANOVA','NumberTitle','off')
set(gcf,'color','w')

k = 1;
for i = 1:dim

    % Plot F-statistic & CIs
    subplot(3,2,k), hold on
    plot(xaxis,f2(i,:),'LineWidth',3)
    plot(xaxis,squeeze(ci2(:,i,:)),'k')
    plot(xaxis(p1(i,:)<=alpha),f2(i,p1(i,:)<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2(i,:)<=alpha),f2(i,p2(i,:)<=alpha),'xr','LineWidth',2)
    xlim([0,21]), ylim([0,6]), box on, grid on
    if i == 1
        title('Test Statistic')
    end
    ylabel(ylabels{i})
    if i ==2
        legend('{\itF}-statistic','95% CI (perm.)')
    end
    if i == dim
        xlabel('permutation')
    end

    % Plot p-values
    subplot(3,2,k+1), hold on
    plot(xaxis,p1(i,:),'k',xaxis,p2(i,:),'--r','LineWidth',2)
    xlim([0,21]), ylim([0,1]), box on, grid on
    if i == 1
        title('{\itP}-values')
        legend('{\itp}-value (param.)','{\itp}-value (perm.)')
    end
    if i == dim
        xlabel('permutation')
    end

    k = k+2;

end