function run_permuranova2_examples
%RUN_PERMURANOVA2_EXAMPLES  Run two-way permutation RM ANOVA examples.
%   Generates random 3D repeated-measures data (Subjects x Factor A x
%   Factor B). Factor A has 2 levels and Factor B has 3 levels. Main
%   effects and an interaction effect are injected along with a subject
%   baseline effect to introduce within-subject correlation. 30 subjects
%   are randomly sampled with replacement 10 times to generate 10 scenarios.
%   Two-way permutation-based repeated-measures ANOVAs are performed.
%   The results are compared to the equivalent parametric statistical test
%   (RM ANOVA) using ranova.m and fitrm.m.
%
%   See also PERMURANOVA2 RANOVA FITRM.
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

% Generate random 3D repeated-measures data (30 subjects, 2x3 design)
rng(42);
nobs = 30; a = 2; b = 3; nperm = 10;

% Add subject baseline effects, main effects, and interaction
sub_eff = randn(nobs,1)*2;
x = randn(nobs,a,b) + sub_eff;
x(:,1,:) = x(:,1,:) + 0.4;
x(:,:,1) = x(:,:,1) + 0.3;
x(:,1,1) = x(:,1,1) + 0.5;
xaxis = 1:nperm; alpha = 0.05;
subjects = 1:nobs;
type = {'ranova2','friedman'};

% Setup MATLAB native fitrm variables
if ~isoctave
    varNames = {'Y1','Y2','Y3','Y4','Y5','Y6'};
    FactorA = categorical([1;2;1;2;1;2]);
    FactorB = categorical([1;1;2;2;3;3]);
    Meas = table(FactorA,FactorB,'VariableNames',{'FactorA','FactorB'});
end

% Initialize matrices for 3 effects (A, B, AxB)
p1 = zeros(nperm,3);
f2 = zeros(nperm,3);
p2 = zeros(nperm,3);
ci2 = zeros(2,3,nperm);

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
        x_curr = x(idx,:,:);

        if isoctave
            p1(i,:) = [NaN, NaN, NaN];
        else
            x_param = reshape(x_curr,nobs,a*b);
            tbl = array2table(x_param,'VariableNames',varNames);
            rm = fitrm(tbl,'Y1-Y6~1','WithinDesign',Meas);
            ranovatbl = ranova(rm,'WithinModel','FactorA*FactorB');
            rows = ranovatbl.Properties.RowNames;
            idxA  = contains(rows, 'FactorA') & ~contains(rows, 'FactorB') & ~contains(rows, 'Error');
            idxB  = contains(rows, 'FactorB') & ~contains(rows, 'FactorA') & ~contains(rows, 'Error');
            idxAB = contains(rows, 'FactorA') &  contains(rows, 'FactorB') & ~contains(rows, 'Error');
            p1(i,1) = ranovatbl.pValue(idxA);
            p1(i,2) = ranovatbl.pValue(idxB);
            p1(i,3) = ranovatbl.pValue(idxAB);
        end

        [f2(i,:),p2(i,:),ci_tmp] = permuranova2(x_curr,'type',type{t});
        ci2(:,:,i) = ci_tmp;

    end

    % Set up figure
    figure('Name',['Two-way ',type{t},' test'],'NumberTitle','off',...
        'Position',[100, 100, 800, 900])
    set(gcf,'color','w')
    effects = {'Factor A', 'Factor B', 'A x B Interaction'};

    for e = 1:3

        % Plot F-statistic & CIs
        subplot(3,2,e*2-1), hold on
        plot(xaxis,f2(:,e),'LineWidth',3)
        plot(xaxis,squeeze(ci2(:,e,:)),'k')
        plot(xaxis(p1(:,e)<=alpha),f2(p1(:,e)<=alpha,e),'ok','LineWidth',2)
        plot(xaxis(p2(:,e)<=alpha),f2(p2(:,e)<=alpha,e),'xr','LineWidth',2)
        xlim([0,nperm+1]), box on, grid on
        title([effects{e},' Test Statistic']), xlabel('permutation'), ylabel('{\itF}-value')
        if e == 1
            legend('{\itF}-statistic','95% CI (perm.)','Location','best')
        end

        % Plot p-values
        subplot(3,2,e*2), hold on
        plot(xaxis,p1(:,e),'k',xaxis,p2(:,e),'--r','LineWidth',2)
        xlim([0,nperm+1]), ylim([0,1]), box on, grid on
        title([effects{e},' {\itP}-values']), xlabel('permutation'), ylabel('probability')
        if e == 1
            if isoctave
                legend('{\itp}-value (param. unavailable)','{\itp}-value (perm.)','Location','best')
            elseif strcmp(type{t},'rank')
                legend('{\itp}-value (param. raw data)','{\itp}-value (perm. rank)','Location','best')
            else
                legend('{\itp}-value (param.)','{\itp}-value (perm.)','Location','best')
            end
        end

    end

end