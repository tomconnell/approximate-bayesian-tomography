% Master script to compare traditional MCMC to ABC for a gravity inversion

clear;clc

addpath('../dram')
addpath('../dram/utils')
addpath('../abc/abcutils')

delete('abc_misfit.csv','analytical_misfit.csv');

% unix('rm abc_misfit.csv analytical_misfit.csv');


%% Set data - common to both ABC and Analytical inversion

true_model = reshape(linspace(-1500,1500,32),[4,8]);

data = forward_grav(true_model);
data = data +(randn(1,length(data))*sqrt(2));


%% Only consider MCMC sampling

nsimu = 5000;
adaptint = 0;
drscale = 0;

start = zeros(1,32)+(rand(1,32)*3000)-1500;

bounds = zeros(2,32);
bounds(1,:) = -1500;
bounds(2,:) = 1500;

qcov = eye(32)*500;

ss_analytical = @(parameters,data,toprint) gravfun(parameters,data,toprint);

% set tolerance
range = mcrange(true_model);
tolerance = ones(1,length(range))*10^-(1.0);
tolerance = tolerance.*range;
    
ss_abc = @(parameters,data,toprint) abcgravfun(parameters,data,tolerance,toprint);

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

update = {'traditional','shift','skew','trunc'};

for ii = 1:4
    
    % unix('rm abc_misfit.csv');
    delete('abc_misfit.csv');
    
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

    [abcresults,abcchain] = dynamicabcdram(model,data,params,options,update{ii});
    
    comparison = figure(1);
    subplot(2,2,ii)
    abc = csvread('abc_misfit.csv');
    ana = csvread('analytical_misfit.csv');
    colors = brewermap(2,'set1');
    plot(ana,'color',colors(1,:),'linewidth',2)
    hold on
    plot(abc,'color',colors(2,:),'linewidth',2)
    legend({'Analytical MCMC','ABC-MCMC'})
    legend boxoff
    text(0.55,0.65,['Analytical \tau = ',num2str(round(mean(iact(anachain))))],'Units','normalized')
    text(0.55,0.6,['Acceptance rate = ',num2str(anaresults.accepted*100),'%'],'Units','normalized')
    text(0.55,0.8,['ABC \tau = ',num2str(round(mean(iact(abcchain))))],'Units','normalized')
    text(0.55,0.75,['Acceptance rate = ',num2str(abcresults.accepted*100),'%'],'Units','normalized')
    ylabel('misfit')
    set(gca,'Ytick',[])
    xlabel('time step')
    title(update{ii})
    
    
    figure(30+ii)
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
    title(update{ii})
    set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
    print(gcf,'-dpdf','-painters',['abc_im_',update{ii},'.pdf'])
    
    figure(20+ii)
    for jj = 1:32
        subplot(8,4,jj)
        par2plot = jj;
        plot(1:length(abcchain),abcchain(:,par2plot),'color',colors(2,:))
        hold on
        plot(linspace(1,length(abcchain)),...
            linspace(true_model(par2plot),true_model(par2plot)),'k-')
        ylabel(['Parameter ',num2str(par2plot)])
        xlabel(['Time',update{ii}])
    end
    set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
    print(gcf,'-dpdf','-painters',['abc_trace_',update{ii},'.pdf'])
    
    figure(10+ii);
    for jj = 1:32
        subplot(8,4,jj)

        par2plot = jj;

        h = histogram(abcchain(:,par2plot),'facecolor',colors(2,:),'facealpha',0.5,'Normalization','pdf');
        topoff = max(h.Values)+0.2*max(h.Values);
        hold on
        plot(linspace(true_model(par2plot),true_model(par2plot)),linspace(0,topoff),'k--','linewidth',2)
        xlabel(['Parameter ',num2str(par2plot)])
        ylim([0,topoff])
        set(gca,'Ytick',[])
        xlim([-1500,1500])
    end
    set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
    print(gcf,'-dpdf','-painters',['abc_hist_',update{ii},'.pdf'])
    
    
end


%% Plotting and comparison

set(comparison,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print(comparison,'-dpdf','-painters','comparison.pdf')

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
set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print(gcf,'-dpdf','-painters','ana_im.pdf')

ana_hist = figure(4);
for ii = 1:32
    subplot(8,4,ii)
    h = histogram(anachain(:,ii),'facecolor',colors(1,:),'facealpha',0.5,'Normalization','pdf');
    topoff = max(h.Values)+0.2*max(h.Values);
    hold on
    plot(linspace(true_model(ii),true_model(ii)),linspace(0,topoff),'k--','linewidth',2)
    xlabel(['Parameter ',num2str(ii)])
    ylim([0,topoff])
    set(gca,'Ytick',[])
    xlim([-1500,1500])
end
set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print(ana_hist,'-dpdf','-painters','ana_hist.pdf')

% Trace plots
ana_trace = figure(6);
for ii = 1:32
    subplot(8,4,ii)
    
    plot(anachain(:,ii),'color',colors(1,:))
    hold on
    plot(linspace(1,length(anachain)),linspace(true_model(ii),true_model(ii)),'k-')
    xlabel(['Parameter ',num2str(ii)])
    ylim([-1500,1500])
    
end
set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print(ana_trace,'-dpdf','-painters','ana_trace.pdf')