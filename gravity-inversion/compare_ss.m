% Master script to compare traditional MCMC to ABC for a gravity inversion

clear;clc

addpath('../dram')
addpath('../dram/utils')
addpath('../abc/abcutils')

% delete('abc_misfit','analytical_misfit.csv');

unix('rm abc_misfit.csv analytical_misfit.csv abc_chain.csv analytical_chain.csv');


%% Set data - common to both ABC and Analytical inversion

true_model = reshape(linspace(-1500,1500,32),[4,8]);
data = forward_grav(true_model);
data = data +(randn(1,length(data))*sqrt(2));


%% Only consider MCMC sampling

nsimu = 10000;
adaptint = 0;
drscale = 0;
seed = round(rand*100);
rng(2)
start = zeros(1,32)+(rand(1,32)*3000)-1500;
rng(seed)
bounds = zeros(2,32);
bounds(1,:) = -1500;
bounds(2,:) = 1500;

qcov = eye(32)*500;

ss_analytical = @(parameters,data,toprint) gravfun(parameters,data,toprint);

ss_prior = @(parameters) smoothness(parameters);


%% For analytical call

clear model params options

model.ssfun = ss_analytical;
model.priorfun = ss_prior;

first_data = forward_grav(reshape(start,[4,8]));
first_misfit = sum((data-first_data).^2);
dlmwrite('analytical_misfit.csv',first_misfit);

params.par0 = start;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = qcov;

% method = analytical (1)
method = 1;

[anaresults,anachain] = gravtestblockeddram(model,data,params,options,method);


%% For ABC call

summary = {'The lot','Residual and marginals','Original','Residual and line','No residual'};

for ii = 1:5
    
    clear range tolerance ss_abc
    
    unix('rm abc_misfit.csv');
    % delete('abc_misfit');
    
    % set tolerance
    range = mcrange2(true_model,ii);
    tolerance = ones(1,length(range))*10^(-2);
    tolerance = tolerance.*range;
    
    % set forward
    ss_abc = @(parameters,data,toprint) abcgravfun2(parameters,data,tolerance,toprint,ii);
    
    clear model params options

    model.ssfun = ss_abc;
    model.priorfun = ss_prior;

    first_data = forward_grav(reshape(start,[4,8]));
    first_misfit = sum((data-first_data).^2);
    dlmwrite('abc_misfit.csv',first_misfit);

    params.par0 = start;
    params.bounds = bounds;

    options.nsimu = nsimu;
    options.adaptint = adaptint;
    options.drscale = drscale;
    options.qcov = qcov;

    % method = abc (2)
    method = 2;

    [abcresults,abcchain] = gravtestblockeddram(model,data,params,options,method);
    
    
    ana = csvread('analytical_misfit.csv');
    abc = csvread('abc_misfit.csv');
    
    figure(1)
    subplot(2,3,ii)
    colors = brewermap(2,'set1');
    plot(ana(1:nsimu),'color',colors(1,:),'linewidth',2)
    hold on
    plot(abc(1:nsimu),'color',colors(2,:),'linewidth',2)
    legend({'Analytical MCMC','ABC-MCMC'})
    legend boxoff
    text(0.45,0.65,['Analytical \tau = ',num2str(round(mean(iact(anachain))))],'Units','normalized')
    text(0.45,0.6,['Acceptance rate = ',num2str(anaresults.accepted*100),'%'],'Units','normalized')
    text(0.45,0.8,['ABC \tau = ',num2str(round(mean(iact(abcchain))))],'Units','normalized')
    text(0.45,0.75,['Acceptance rate = ',num2str(abcresults.accepted*100),'%'],'Units','normalized')
    ylabel('misfit')
    set(gca,'Ytick',[])
    xlabel('time step')
    title(summary{ii})
    
    figure(3+ii)
    subplot(5,8,1:8)
    plot(1:8,data,'k-')
    hold on
    plot(1:8,data,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
    abc_solution = forward_grav(reshape(mean(abcchain),[4,8]));
    abc_smoothness = smoothness(reshape(mean(abcchain),[4,8]));
    abc_misfit = sum((abc_solution-data).^2);
    plot(1:8,abc_solution,'-','color',colors(2,:))
    plot(1:8,abc_solution,'o','MarkerEdgeColor', colors(2,:),'MarkerFaceColor',colors(2,:))
    text(0.05,0.8,['Smoothness = ',num2str(round(abc_smoothness))],'Units','normalized')
    text(0.05,0.5,['Misfit = ',num2str(round(abc_misfit))],'Units','normalized')
    title('ABC-MCMC')
    subplot(5,8,9:40)
    colormap(brewermap([],'PiYG'));
    imagesc(reshape(mean(abcchain),[4,8]))
    set(gca,'Ytick',[],'Xtick',[])
    title(summary{ii})
    
end


%% Plotting and comparison

% Plot chains at their mean positions
figure(2)
subplot(5,8,1:8)
plot(1:8,data,'k-')
hold on
plot(1:8,data,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
analytical_solution = forward_grav(reshape(mean(anachain),[4,8]));
analytical_smoothness = smoothness(reshape(mean(anachain),[4,8]));
analytical_misfit = sum((analytical_solution-data).^2);
plot(1:8,analytical_solution,'-','color',colors(1,:))
plot(1:8,analytical_solution,'o','MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
text(0.05,0.8,['Smoothness = ',num2str(round(analytical_smoothness))],'Units','normalized')
text(0.05,0.5,['Misfit = ',num2str(round(analytical_misfit))],'Units','normalized')
title('Analytical MCMC')
subplot(5,8,9:40)
colormap(brewermap([],'PiYG'));
imagesc(reshape(mean(anachain),[4,8]))
set(gca,'Ytick',[],'Xtick',[])

rand