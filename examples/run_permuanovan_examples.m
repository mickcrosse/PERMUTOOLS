function run_permuanovan_examples
%RUN_PERMUANOVAN_EXAMPLES  Run permutation-based N-way ANOVA examples.
%   Generates random data for a 3-way design (2x2x2) along the columns
%   of X. X contains 20 dependent variables. The groups have identical
%   means of 0 for most variables,except for the first 5 variables where
%   Factor 1 induces a shift. Each cell in the 3-way design has 5
%   observations (40 total).
%
%   Permutation-based N-way ANOVAs are performed on the datasets
%   using standard and rank-transformed methods. The results are compared
%   to those of the equivalent parametric statistical test using MATLAB's
%   native anovan.m function.
%
%   See also PERMUANOVAN ANOVAN.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ,Foxe JJ,Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL,Albert Einstein College of Medicine,NY.
%   TCBE,Trinity College Dublin,Ireland.

% Set up experiment
nobs = 40;
nvar = 20;
xaxis = 1:nvar;
alpha = 0.05;

% 3-way design: 2x2x2 (Factor 1,Factor 2,Factor 3)
group1 = [ones(20,1); 2*ones(20,1)];
group2 = repmat([ones(10,1); 2*ones(10,1)],2,1);
group3 = repmat([ones(5,1); 2*ones(5,1)],4,1);
group = [group1,group2,group3];
type = {'anovan','rank'};

% Generate random 2D data
rng(42);
x = randn(nobs,nvar);
% Add main effect to Factor 1 for the first 5 variables
x(group(:,1)==1,1:5) = x(group(:,1)==1,1:5) + 1.0;

for t = 1:numel(type)

    % Set up figure
    figure('Name',['N-way (',type{t},') - Factor 1'],'NumberTitle','off')
    set(gcf,'color','w')

    % Run parametric ANOVAN equivalent
    tic
    switch type{t}
        case 'anovan'
            xparam = x;
        case 'rank'
            xparam = reshape(tiedrank(x),nobs,nvar);
    end
    f1 = zeros(1,nvar);
    p1 = zeros(1,nvar);
    fcrit1 = zeros(1,nvar);
    for v = 1:nvar
        [panova,tblanova] = anovan(xparam(:,v),group,'model','linear','display','off');
        p1(v) = panova(1);
        colF = find(strcmp(tblanova(1,:),'F'));
        colDF = find(strcmp(tblanova(1,:),'d.f.'));
        f1(v) = tblanova{2,colF};
        dfA = tblanova{2,colDF};
        dfE = tblanova{end-1,colDF};
        fcrit1(v) = finv(1-alpha,dfA,dfE);
    end
    ci1 = f1 ./ fcrit1;
    fprintf('Elapsed time (Param.): %.3f s\n',toc)

    % Run permutation N-way ANOVA (uncorrected)
    tic
    [f2,p2,ci2] = permuanovan(x,group,'type',type{t},'correct',0);
    fprintf('Elapsed time (Perm. uncorr.): %.3f s\n',toc)

    % Isolate Factor 1 for plotting
    f2a = f2(1,:);
    p2a = p2(1,:);
    ci2a = squeeze(ci2(1,1,:))';

    % Plot uncorrected test statistic & CIs
    subplot(2,2,1),hold on
    plot(xaxis,f1,'LineWidth',2)
    plot(xaxis,f2a,'--','LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci2a,'--r')
    plot(xaxis(p1<=alpha),f1(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2a<=alpha),f2a(p2a<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,20]), box on, grid on
    title('Uncorrected'), ylabel('{\itF}-value')
    legend('{\itF}-statistic (param.)','{\itF}-statistic (perm.)',...
        '95% CI (param.)','95% CI (perm.)','Location','best')

    % Plot uncorrected p-values
    subplot(2,2,3),hold on
    plot(xaxis,p1,'k',xaxis,p2a,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    xlabel('variable'), ylabel('probability')
    legend('{\itp}-value (param.)','{\itp}-value (perm.)','Location','best')

    % Run permutation N-way ANOVA (max-corrected)
    tic
    [f2,p2,ci2] = permuanovan(x,group,'type',type{t},'correct',1);
    fprintf('Elapsed time (Perm. max-corr.): %.3f s\n',toc)

    % Isolate Factor 1 for plotting
    f2a = f2(1,:);
    p2a = p2(1,:);
    ci2a = squeeze(ci2(1,1,:))';

    % Plot max-corrected test statistic & CIs
    subplot(2,2,2),hold on
    plot(xaxis,f1,'LineWidth',2)
    plot(xaxis,f2a,'--','LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci2a,'--r')
    plot(xaxis(p1<=alpha),f1(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2a<=alpha),f2a(p2a<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,20]), box on, grid on
    title('Max-corrected')

    % Plot max-corrected p-values
    subplot(2,2,4),hold on
    plot(xaxis,p1,'k',xaxis,p2a,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    xlabel('variable')

end