% Master script to compare traditional MCMC to ABC for a gravity inversion

clear;clc

addpath('../dram')
addpath('../dram/utils')
addpath('../abc/abcutils')

delete('abc_misfit.csv','analytical_misfit.csv','abc_chain.csv','analytical_chain.csv');
%unix('rm abc_misfit.csv analytical_misfit.csv abc_chain.csv analytical_chain.csv');


%% Set data - common to both ABC and Analytical inversion

true_model = reshape(linspace(2000,3500,32),[4,8]);
% true_model = reshape(linspace(-750,750,32),[4,8]);

data = forward_grav(true_model);
no_edge_effect = data;
%data = data +(randn(1,length(data))*sqrt(2));
edge_effect = edge_effect_forward_grav([2750,0,2750]);
data = data+edge_effect;

% plot true model
figure(8)
subplot(5,8,1:7)
plot(1:8,data,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
plot(1:8,data,'k-','HandleVisibility','off')
set(gca,'xtick',[])
ylabel('\Delta g')
xlim([0.5,8.25])
true_model_smoothness = smoothness(true_model);
text(0.025,0.8,['Smoothness = ',num2str(round(true_model_smoothness))]...
    ,'Units','normalized')
subplot(5,8,9:40)
colormap(brewermap([],'purples'))
imagesc(true_model/1000)
colorbar()
set(gca,'ytick',[],'xtick',[])
set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print('-dpdf','-painters','true_model.pdf')


%% Only consider MCMC sampling

nsimu = 100000;
adaptint = 0;
drscale = 0;
seed = round(rand*100);
rng(56)
start = zeros(1,32)+(rand(1,32)*1500)+2000;
rng(seed)
bounds = zeros(2,32);
bounds(1,:) = 2000;
bounds(2,:) = 3500;

qcov = eye(32)*250;

ss_analytical = @(parameters,data,toprint) gravfun(parameters,data,toprint);

% set tolerance
range = mcrange(data);
% factor = 10^-(2);
factor = 10^(-1.5);
tolerance = ones(1,length(range))*factor;
tolerance = tolerance.*range;

ss_abc = @(parameters,data,toprint) abcgravfun(parameters,data,tolerance,toprint);
ss_prior = @(parameters) smoothness(parameters)*(10^(-1));
ss_prior_abc = @(parameters) smoothness(parameters);


%% For analytical call

clear model params options

model.ssfun = ss_analytical;
model.priorfun = ss_prior;

first_data = forward_grav(reshape(start,[4,8]))+edge_effect;
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
csvwrite('analytical_chain.csv',anachain)

%% For ABC call

clear model params options

model.ssfun = ss_abc;
model.priorfun = ss_prior_abc;

% first_data = forward_grav(reshape(start,[4,8]));
% first_misfit = sum((data-first_data).^2);
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
csvwrite('abc_chain.csv',abcchain)


%% Plotting and comparison

ana = csvread('analytical_misfit.csv');
abc = csvread('abc_misfit.csv');

figure(1)
colors = brewermap(2,'set1');
plot(ana(1:10000),'color',colors(1,:),'linewidth',2)
hold on
plot(abc(1:10000),'color',colors(2,:),'linewidth',2)
legend({'Analytical MCMC','ABC-MCMC'})
legend boxoff
text(0.55,0.8,['Analytical \tau = ',num2str(round(mean(iact(anachain))))],'Units','normalized')
text(0.55,0.75,['Acceptance rate = ',num2str(anaresults.accepted*100),'%'],'Units','normalized')
text(0.55,0.65,['ABC \tau = ',num2str(round(mean(iact(abcchain))))],'Units','normalized')
text(0.55,0.6,['Acceptance rate = ',num2str(abcresults.accepted*100),'%'],'Units','normalized')
text(0.55,0.55,['tolerance = ',num2str(factor)],'Units','normalized')
ylabel('misfit')
set(gca,'Ytick',[])
xlabel('time step')

set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print('-dpdf','-painters','comparison.pdf')

% Plot chains at their mean positions
ana_im = figure(2);
subplot(5,8,1:8)
plot(1:8,data,'k-')
hold on
plot(1:8,data,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
xlim([0.5,8.5])
analytical_solution = forward_grav(reshape(mean(anachain),[4,8]))+edge_effect;
analytical_smoothness = smoothness(reshape(mean(anachain),[4,8]));
analytical_misfit = sum((analytical_solution-data).^2);
plot(1:8,analytical_solution,'-','color',colors(1,:))
plot(1:8,analytical_solution,'o','MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
text(0.05,0.8,['Smoothness = ',num2str(round(analytical_smoothness))],'Units','normalized')
text(0.05,0.5,['Misfit = ',num2str(round(analytical_misfit))],'Units','normalized')
title('Analytical MCMC')
subplot(5,8,9:40)
colormap(brewermap([],'purples'));
imagesc(reshape(mean(anachain),[4,8]))
set(gca,'Ytick',[],'Xtick',[])

set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print(ana_im,'-dpdf','-painters','ana_im.pdf')

abc_im = figure(3);
subplot(5,8,1:8)
plot(1:8,data,'k-')
hold on
plot(1:8,data,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k')
xlim([0.5,8.5])
abc_solution = forward_grav(reshape(mean(abcchain),[4,8]))+edge_effect;
abc_smoothness = smoothness(reshape(mean(abcchain),[4,8]));
abc_misfit = sum((abc_solution-data).^2);
plot(1:8,abc_solution,'-','color',colors(2,:))
plot(1:8,abc_solution,'o','MarkerEdgeColor', colors(2,:),'MarkerFaceColor',colors(2,:))
text(0.05,0.8,['Smoothness = ',num2str(round(abc_smoothness))],'Units','normalized')
text(0.05,0.5,['Misfit = ',num2str(round(abc_misfit))],'Units','normalized')
title('ABC-MCMC')
subplot(5,8,9:40)
colormap(brewermap([],'purples'));
imagesc(reshape(mean(abcchain),[4,8]))
set(gca,'Ytick',[],'Xtick',[])

set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print(abc_im,'-dpdf','-painters','abc_im.pdf')

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
    xlim([2000,3500])
end
set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print(ana_hist,'-dpdf','-painters','ana_hist.pdf')

abc_hist = figure(5);
for ii = 1:32
    subplot(8,4,ii)
    h = histogram(abcchain(:,ii),'facecolor',colors(2,:),'facealpha',0.5,'Normalization','pdf');
    topoff = max(h.Values)+0.2*max(h.Values);
    hold on
    plot(linspace(true_model(ii),true_model(ii)),linspace(0,topoff),'k--','linewidth',2)
    xlabel(['Parameter ',num2str(ii)])
    ylim([0,topoff])
    set(gca,'Ytick',[])
    xlim([2000,3500])
end
set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print(abc_hist,'-dpdf','-painters','abc_hist.pdf')

ana_trace = figure(6);
for ii = 1:32
    subplot(8,4,ii)
    
    plot(anachain(:,ii),'color',colors(1,:))
    hold on
    plot(linspace(1,length(anachain)),linspace(true_model(ii),true_model(ii)),'k-')
    xlabel(['Parameter ',num2str(ii)])
    ylim([2000,3500])
    
end
set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print(ana_trace,'-dpdf','-painters','ana_trace.pdf')

abc_trace = figure(7);
for ii = 1:32
    subplot(8,4,ii)
    
    plot(abcchain(:,ii),'color',colors(2,:))
    hold on
    plot(linspace(1,length(anachain)),linspace(true_model(ii),true_model(ii)),'k-')
    xlabel(['Parameter ',num2str(ii)])
    ylim([2000,3500])
    
end
set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print(abc_trace,'-dpdf','-painters','abc_trace.pdf')

% New plot for thesis
marginals = figure(8);
pars = 2:4:30;
for ii = 1:8
    
    par2plot = pars(ii);
    
    subplot(16,8,[(((ii*2)-1)*8)-7,(((ii*2)-1)*8)-6,((ii*2)*8)-7,((ii*2)*8)-6]);
    
    h = histogram(anachain(:,par2plot),20,'normalization','pdf',...
        'facecolor',colors(1,:),'facealpha',0.3);
    hold on
    histogram(abcchain(:,par2plot),20,'normalization','pdf',...
        'facecolor',colors(2,:),'facealpha',0.3);
    
    topoff = max(h.Values)+0.2*max(h.Values);
    ylim([0,topoff])
    
    xlabel(['Parameter ',num2str(par2plot)])
    
    set(gca,'view',[-90,90],'xlim',[2000,3500]...
        ,'xtick',[],'ytick',[])
    
    subplot(16,8,[(((ii*2)-1)*8)-5:((ii*2)-1)*8,((ii*2)*8)-5:(ii*2)*8])
    
    plot(anachain(:,par2plot),'color',colors(1,:))
    hold on
    plot(abcchain(:,par2plot),'color',colors(2,:))
    plot(linspace(1,length(anachain)),linspace(true_model(par2plot),...
        true_model(par2plot)),'k-')
    ylim([2000,3500])
    set(gca,'Xtick',[],'Ytick',[])
    
end
set(gcf,'units','centimeters','position',[0,0,25,30],'papersize',[25,30])
print(marginals,'-dpdf','-painters','marginals.pdf')