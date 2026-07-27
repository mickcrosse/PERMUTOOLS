function run_permuranova1_examples
%RUN_PERMURANOVA1_EXAMPLES  Run one-way permutation RM ANOVA examples.
%   Generates random repeated-measures data for 5 conditions along the
%   columns of X. Each condition has a mean value of 0, except for the
%   first condition which has a mean value of 1. A random subject baseline
%   effect is added to introduce within-subject correlation. Each condition
%   has 30 observations (subjects). The subjects are randomly sampled with
%   replacement 20 times to generate 20 different scenarios. One-way
%   permutation-based repeated-measures ANOVAs are performed on each of
%   the datasets. The results are compared to those of the equivalent
%   parametric statistical test (i.e. RM ANOVA) using fitrm.m and ranova.m,
%   and the non-parametric statistical test (i.e. Friedman test) using
%   friedman.m.
%
%   See also PERMURANOVA1 RANOVA FRIEDMAN.
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

% Generate random repeated-measures data
rng(42);
nobs = 30; nvar = 5; nperm = 20;

% Add subject baseline effects to create correlated within-subject data
sub_eff = randn(nobs,1)*2;
x = randn(nobs,nvar) + sub_eff;
x(:,1) = x(:,1)+1;

xaxis = 1:nperm; alpha = 0.05;
subjects = 1:nobs;
type = {'ranova1','friedman'};

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
            idx = datasample(subjects,nobs,'Replace',true);
        else
            idx = datasample(s,subjects,nobs,'Replace',true);
        end
        x_curr = x(idx,:);

        switch type{t}
            case 'ranova1'
                if isoctave
                    p1(i) = NaN;
                else
                    varNames = {'Y1','Y2','Y3','Y4','Y5'};
                    tbl = array2table(x_curr,'VariableNames',varNames);
                    Meas = table((1:nvar)','VariableNames',{'Condition'});
                    rm = fitrm(tbl,'Y1-Y5~1','WithinDesign',Meas);
                    ranovatbl = ranova(rm);
                    p1(i) = ranovatbl.pValue(1);
                end
            case 'friedman'
                p1(i) = friedman(x_curr,1,'off');
        end
        [f2(i),p2(i),ci2(:,i)] = permuranova1(x_curr,'type',type{t});

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

    if strcmp(type{t}, 'ranova') && isoctave
        legend('{\itp}-value (param. unavailable)','{\itp}-value (perm.)')
    else
        legend('{\itp}-value (param.)','{\itp}-value (perm.)')
    end

end