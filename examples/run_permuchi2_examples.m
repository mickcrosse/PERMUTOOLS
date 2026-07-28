function run_permuchi2_examples
%RUN_PERMUCHI2_EXAMPLES  Run permutation Chi-square test examples.
%   Generates random multivariate categorical data for a sample X,
%   containing 20 variables. The first 10 variables contain a weak
%   association, while the remaining 10 variables are completely
%   independent. Each variable has 100 observations. Permutation tests of
%   independence based on the Chi-square statistic are performed on each
%   variable for two-tailed, right-tailed and left-tailed tests. The
%   results are compared to those of the equivalent parametric statistical
%   tests using chi2cdf.
%
%   See also PERMUCHI2 CROSSTAB CHI2GOF.
%
%   PERMUTOOLS https://github.com/mickcrosse/PERMUTOOLS
%   References:
%       [1] Crosse MJ, Foxe JJ, Molholm S (2024) PERMUTOOLS: A MATLAB
%           Package for Multivariate Permutation Testing. arXiv 2401.09401.
%   © 2018-2026 Mick Crosse <crossemj@tcd.ie>
%   CNL, Albert Einstein College of Medicine, NY.
%   TCBE, Trinity College Dublin, Ireland.

% Set up experiment
nobs = 100; nvar = 20; df = 1;
xaxis = 1:nvar; alpha = 0.05;
tail = {'both','right','left'};
label = {'two','right','left'};
type = 'permuchi2';

% Generate random categorical data with a weak effect in the first half
rng(42);
x = zeros(2,2,nvar);
varA = [ones(nobs/2,1);2*ones(nobs/2,1)];
for j = 1:nvar
    if j <= round(nvar/2)
        % Weak association: ~60% match
        varB = varA;
        flip = rand(nobs,1)>0.60;
        varB(flip) = 3-varB(flip);
    else
        % No association: ~50% match (chance)
        varB = randi([1,2],nobs,1);
    end
    x(:,:,j) = accumarray([varA,varB],1,[2,2]);
end

% Plot parametric & permutation CIs and statistics
figure('Name',['One-sample ',type,': statistic & CIs'],'NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [chi2,p2] = permuchi2(x,'tail',tail{i},'correct',0);
    p1 = zeros(1,nvar);
    for j = 1:nvar
        switch tail{i}
            case 'right'
                p1(j) = 1-chi2cdf(chi2(j),df);
            case 'left'
                p1(j) = chi2cdf(chi2(j),df);
            case 'both'
                pl = chi2cdf(chi2(j),df);
                pr = 1-chi2cdf(chi2(j),df);
                p1(j) = min(1,2*min(pl,pr));
        end
    end
    subplot(3,2,i+i-1), hold on
    plot(xaxis,chi2,'LineWidth',3)
    plot(xaxis(p1<=alpha),chi2(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),chi2(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,25]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 1
        legend('Chi2 statistic','Location','best')
    end
    [~,p2] = permuchi2(x,'tail',tail{i},'correct',1);
    subplot(3,2,i+i), hold on
    plot(xaxis,chi2,'LineWidth',3)
    plot(xaxis(p1<=alpha),chi2(p1<=alpha),'ok','LineWidth',2)
    plot(xaxis(p2<=alpha),chi2(p2<=alpha),'xr','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,25]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end

end

% Plot parametric & permutation p-values
figure('Name',['One-sample ',type,': p-values'],'NumberTitle','off')
set(gcf,'color','w')
for i = 1:numel(tail)
    [chi2,p2] = permuchi2(x,'tail',tail{i},'correct',0);
    p1 = zeros(1,nvar);
    for j = 1:nvar
        switch tail{i}
            case 'right'
                p1(j) = 1 - chi2cdf(chi2(j), df);
            case 'left'
                p1(j) = chi2cdf(chi2(j), df);
            case 'both'
                pl = chi2cdf(chi2(j), df);
                pr = 1 - chi2cdf(chi2(j), df);
                p1(j) = min(1, 2*min(pl, pr));
        end
    end
    subplot(3,2,i+i-1), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    if i == 1
        title('Uncorrected')
    elseif i == 3
        xlabel('variable')
    end
    ylabel([label{i},'-tailed'])
    if i == 1
        legend('{\itp}-value (param.)','{\itp}-value (perm.)',...
            'Location','best')
    end
    [~,p2] = permuchi2(x,'tail',tail{i},'correct',1);
    subplot(3,2,i+i), hold on
    plot(xaxis,p1,'k',xaxis,p2,'--r','LineWidth',2)
    xlim([0,nvar+1]), ylim([0,1]), box on, grid on
    if i == 1
        title('Max-corrected')
    elseif i == 3
        xlabel('variable')
    end
end