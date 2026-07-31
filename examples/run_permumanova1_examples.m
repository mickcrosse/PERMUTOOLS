function run_permumanova1_examples
%RUN_PERMUMANOVA1_EXAMPLES  Run permutation-based one-way MANOVA examples.
%   Generates random multivariate data for 3 groups along the rows of X.
%   X contains 20 sets of dependent variables (2 dimensions per set).
%   The groups have identical multivariate means of 0 for most variables,
%   except for the first 5 variables where Group 1 and Group 2 are shifted
%   in different dimensions. Each group has 10 observations (30 total).
%   Permutation-based one-way MANOVAs are performed on the datasets
%   using standard and rank-transformed methods. The results are compared
%   to those of the equivalent parametric statistical test (one-way MANOVA)
%   using manova1.m.
%
%   See also PERMUMANOVA1 MANOVA1.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS

%   References:
%       [1] Crosse MJ,Foxe JJ,Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.

%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL,Albert Einstein College of Medicine,NY.
%   TCBE,Trinity College Dublin,Ireland.

% Set up experiment
nobs = 30; nvar = 20; ndep = 2;
xaxis = 1:nvar; alpha = 0.05;
group = repmat(1:3,1,10)';
type = {'manova1','rank'};

% Generate random 3D data
rng(42);
x = randn(nobs,ndep,nvar);
x(group==1,1,1:5) = x(group==1,1,1:5) + 1.5;
x(group==2,2,1:5) = x(group==2,2,1:5) + 1.5;

for t = 1:numel(type)

    % Set up figure
    figure('Name',['One-way (',type{t},')'],'NumberTitle','off')
    set(gcf,'color','w')

    % Run native parametric MANOVA
    tic
    switch type{t}
        case 'manova1'
            xparam = x;
        case 'rank'
            xparam = reshape(tiedrank(reshape(x,nobs,[])),nobs,ndep,nvar);
    end
    p1 = zeros(1,nvar);
    lambda1 = zeros(1,nvar);
    for v = 1:nvar
        [~,pmanova,stats] = manova1(xparam(:,:,v),group);
        p1(v) = pmanova(1);
        lambda1(v) = stats.lambda(1);
    end
    kgroups = numel(unique(group));
    df = ndep*(kgroups-1);
    bartlettC = nobs-1-(ndep+kgroups)/2;
    chi2crit = chi2inv(1-alpha,df);
    crit = exp(-chi2crit/bartlettC);
    ci1 = [zeros(1,nvar);lambda1./crit];
    fprintf('Elapsed time (Param.): %.3f s\n',toc)

    % Run Permutation MANOVA (Uncorrected)
    tic
    [lambda2,p2,ci2] = permumanova1(x,group,'type',type{t},'correct',0);
    fprintf('Elapsed time (Perm. uncorr.): %.3f s\n',toc)
    ci2 = squeeze(ci2);

    % Plot uncorrected test statistic & CIs
    subplot(2,2,1),hold on
    plot(xaxis,lambda1,'LineWidth',2)
    plot(xaxis,lambda2,'--','LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),lambda1(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),lambda2(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,2]), box on, grid on
    title('Uncorrected'),ylabel('Wilks'' \Lambda')
    legend('Wilks'' \Lambda (param.)','Wilks'' \Lambda (perm.)',...
        '95% CI (param.)','95% CI (perm.)','Location','best')

    % Plot uncorrected p-values
    subplot(2,2,3),hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    xlabel('variable'), ylabel('probability')
    legend('{\itp}-value (param.)','{\itp}-value (perm.)','Location','best')

    % Run permutation MANOVA (max-corrected)
    tic
    [lambda2,p2,ci2] = permumanova1(x,group,'type',type{t},'correct',1);
    fprintf('Elapsed time (Perm. max-corr.): %.3f s\n',toc)
    ci2 = squeeze(ci2);

    % Plot max-corrected test statistic & CIs
    subplot(2,2,2),hold on
    plot(xaxis,lambda1,'LineWidth',2)
    plot(xaxis,lambda2,'--','LineWidth',2)
    plot(xaxis,ci1,'k',xaxis,ci2,'--r')
    plot(xaxis(p1<=alpha),lambda1(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),lambda2(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,2]), box on, grid on
    title('Max-corrected')

    % Plot max-corrected p-values
    subplot(2,2,4),hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    xlabel('variable')

end