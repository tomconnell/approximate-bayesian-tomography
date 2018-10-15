% Master script for joint inversion

clear;clc

addpath('../dram')
addpath('../dram/utils')
addpath('../abc/abcutils')

delete('abc_misfit.csv','ana_misfit.csv','abc_tom_misfit.csv',...
    'abc_gravity_misfit.csv','gravity_misfit.csv','tom_misfit.csv',...
    'gravity_likelihood.csv','tom_likelihood.csv',...
    'abc_gravity_likelihood.csv','abc_tom_likelihood.csv');


%% Set data - common to both ABC and Analytical inversion

true_model_density = reshape(linspace(2000,3500,120),[6,20]); % units kg/m3

edge_effect = edge_effect_forward_grav([2750,0,2750]);

observed_data_gravity = forward_grav(true_model_density) + edge_effect;

% brocher function includes conversion of density from kg/m3 to g/cm3
true_model_vp = brocherizegrid(true_model_density); % units km/s

n_receivers = 10;
n_sources = 4;
observed_data_tomography = tom(true_model_vp,n_receivers,n_sources);

%% Plot true density model with gravity data

figure(8)
subplot(5,8,1:7)
plot(1:20,observed_data_gravity,'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
plot(1:20,observed_data_gravity,'k-','HandleVisibility','off')
set(gca,'xtick',[])
ylabel('\Delta g')
xlim([0.5,20])
true_model_smoothness = smoothness(true_model_density,20,6);
text(0.025,0.8,['Smoothness = ',num2str(round(true_model_smoothness))]...
    ,'Units','normalized')
subplot(5,8,9:40)
colormap(brewermap([],'purples'))
imagesc(true_model_density/1000)
colorbar()
set(gca,'ytick',[],'xtick',[])
set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print('-dpdf','-painters','true_model_gravity.pdf')


%% Sampling considerations

% nsimu = 100000;
nsimu = 60000;
adaptint = 0;
drscale = 0;

seed = round(rand*100);
rng(15)
% uniform random in the range 2000-3500 kg/m3
oldchain = csvread('abchalfchain.csv');
start = oldchain(40000,:);
rng(seed)
bounds = zeros(2,120);
bounds(1,:) = 2000;
bounds(2,:) = 3500;

qcov = eye(120)*250; % 1/6 of the total parameter space

ss_analytical = @(parameters,data,toprint) ...
    anafun(parameters,data,edge_effect,toprint);

% range for gravity
grange = mcrange(true_model_density,edge_effect);
% range for tomography
trange = tomrange(true_model_density);

% common tolerance factor
factor = 0.05;

% set tolerance values
gtolerance = grange*factor;
ttolerance = trange*factor;

ss_abc = @(parameters,data,toprint) ...
    abcfun(parameters,data,edge_effect,gtolerance,ttolerance,toprint);

% prior_ss = @(parameters) smoothness(reshape(parameters,[6,20]),6,20);


%% Analytical call

clear model params options

model.ssfun = ss_analytical;
% model.priorfun = prior_ss;

data = [observed_data_gravity,observed_data_tomography];

first_gravity = forward_grav(reshape(start,[6,20])) + edge_effect;

first_tomography = tom(brocherizegrid(reshape(start,[6,20]))...
    ,n_receivers,n_sources);

first_misfit = sum((data-[first_gravity,first_tomography]).^2);

dlmwrite('ana_misfit.csv',first_misfit);
dlmwrite('abc_gravity_misfit.csv',sum((observed_data_gravity-first_gravity).^2));
dlmwrite('abc_tom_misfit.csv',sum((observed_data_tomography-first_tomography).^2));
dlmwrite('gravity_misfit.csv',sum((observed_data_gravity-first_gravity).^2));
dlmwrite('tom_misfit.csv',sum((observed_data_tomography-first_tomography).^2));

dlmwrite('gravity_likelihood.csv',0);
dlmwrite('tom_likelihood.csv',0);

params.par0 = start;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = qcov;

% method = analytical (1)
method = 1;

[anaresults,anachain] = gravtestblockeddram(model,data,params,options,method);
csvwrite('anachain.csv',anachain);

%% ABC call

clear model params options

model.ssfun = ss_abc;
% model.priorfun = prior_ss;

dlmwrite('abc_misfit.csv',first_misfit);
dlmwrite('abc_gravity_likelihood.csv',0);
dlmwrite('abc_tom_likelihood.csv',0);

params.par0 = start;
params.bounds = bounds;

options.nsimu = nsimu;
options.adaptint = adaptint;
options.drscale = drscale;
options.qcov = qcov;

% method = abc (2)
% method = 2;

% [abcresults,abcchain] = gravtestblockeddram(model,data,params,options,method);

[abcresults,abcchain] = jointdynamicdram(model,data,params,options,'trunc');
csvwrite('abcchain.csv',abcchain);

% [abcresults,abcchain] = skewdram(model,data,params,options);


%% Plotting

anachain = csvread('anachain.csv');

ana = csvread('ana_misfit1.csv');

% abc = csvread('abc_misfit.csv');

tom_abc = csvread('abc_tom_misfit.csv');

grav_abc = csvread('abc_gravity_misfit.csv');

tom_ana = csvread('tom_misfit.csv');

grav_ana = csvread('gravity_misfit.csv');

%% Stitch together ABC chains

abc1 = csvread('abcchain1.csv');
abc2 = csvread('abcchain2.csv');
abcchain = vertcat(abc1,abc2);
csvwrite('abcchain.csv',abcchain)

% figure(1)

colors = brewermap(2,'set1');

% plot(ana(1:10000),'color',colors(1,:),'linewidth',2)

% hold on

% plot(abc(1:10000),'color',colors(2,:),'linewidth',2)

% legend({'Analytical MCMC','ABC-MCMC'})

% legend boxoff

% text(0.7,0.8,['Analytical \tau = ',num2str(round(mean(iact(anachain))))],'Units','normalized')

% text(0.7,0.75,['Acceptance rate = ',num2str(anaresults.accepted*100),'%'],'Units','normalized')

% text(0.7,0.65,['ABC \tau = ',num2str(round(mean(iact(abcchain))))],'Units','normalized')

% text(0.7,0.6,['Acceptance rate = ',num2str(abcresults.accepted*100),'%'],'Units','normalized')

% text(0.7,0.55,['tolerance = ',num2str(factor)],'Units','normalized')

% ylabel('misfit')

% set(gca,'Ytick',[])

% xlabel('time step')




%% Comparison figure

% figure(1)
% 
% colors2 = brewermap(6,'paired');
% 
% plot(normalize(grav_ana(2:10000)),'color',colors2(6,:),'linewidth',2)
% 
% hold on
% 
% plot(normalize(tom_ana(2:10000)),'color',colors2(5,:),'linewidth',2)
% 
% plot(normalize(grav_abc(2:10000)),'color',colors2(2,:),'linewidth',2)
% 
% plot(normalize(tom_abc(2:10000)),'color',colors2(1,:),'linewidth',2)
% 
% legend({'Gravity - Analytical MCMC','Tomography - Analytical MCMC',...
%     'Gravity - ABC-MCMC','Tomography - ABC-MCMC'})
% 
% legend boxoff
% 
% text(0.7,0.7,['Analytical \tau = ',num2str(round(mean(iact(anachain))))],'Units','normalized')
% 
% text(0.7,0.65,['Acceptance rate = ',num2str(anaresults.accepted*100),'%'],'Units','normalized')
% 
% text(0.7,0.55,['ABC \tau = ',num2str(round(mean(iact(abcchain))))],'Units','normalized')
% 
% text(0.7,0.5,['Acceptance rate = ',num2str(abcresults.accepted*100),'%'],'Units','normalized')
% 
% text(0.7,0.45,['tolerance = ',num2str(factor)],'Units','normalized')
% 
% ylabel('normalized misfit')
% 
% set(gca,'Ytick',[])
% 
% xlabel('time step')
% 
% set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
% 
% print('-dpdf','-painters','comparison.pdf')

%% ANA IM

% Plot chains at their mean positions
ana_im = figure(2);
subplot(7,20,1:18)
plot(1:20,observed_data_gravity,'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
plot(1:20,observed_data_gravity,'k-')
analytical_solution = forward_grav(reshape(mean(anachain),[6,20]))+edge_effect;
analytical_smoothness = smoothness(reshape(mean(anachain),[6,20]),6,20);
analytical_misfit = sum((analytical_solution-observed_data_gravity).^2);
plot(1:20,analytical_solution,'-','color',colors(1,:))
plot(1:20,analytical_solution,'x','MarkerEdgeColor', colors(1,:),'MarkerFaceColor',colors(1,:))
text(0.05,0.8,['Smoothness = ',num2str(round(analytical_smoothness))],'Units','normalized')
text(0.05,0.5,['Misfit = ',num2str(round(analytical_misfit))],'Units','normalized')
title('Analytical MCMC')
xlim([0.75,20.25])
subplot(7,20,21:140)
colormap(brewermap([],'purples'));
imagesc(reshape(mean(anachain)/1000,[6,20]))
colorbar()
set(gca,'Ytick',[],'Xtick',[])

set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print(ana_im,'-dpdf','-painters','ana_im.pdf')

%% ABC IM

abc_im = figure(3);
subplot(7,20,1:18)
plot(1:20,observed_data_gravity,'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
hold on
plot(1:20,observed_data_gravity,'k-')
abc_solution = forward_grav(reshape(mean(abcchain),[6,20]))+edge_effect;
abc_smoothness = smoothness(reshape(mean(abcchain),[6,20]),6,20);
abc_misfit = sum((abc_solution-observed_data_gravity).^2);
plot(1:20,abc_solution,'-','color',colors(2,:))
plot(1:20,abc_solution,'x','MarkerEdgeColor', colors(2,:),'MarkerFaceColor',colors(2,:))
text(0.05,0.8,['Smoothness = ',num2str(round(abc_smoothness))],'Units','normalized')
text(0.05,0.5,['Misfit = ',num2str(round(abc_misfit))],'Units','normalized')
title('ABC-MCMC')
xlim([0.75,20.25])
subplot(7,20,21:140)
colormap(brewermap([],'purples'));
imagesc(reshape(mean(abcchain)/1000,[6,20]))
colorbar()
set(gca,'Ytick',[],'Xtick',[])

set(gcf,'units','centimeters','position',[0,0,20,10],'papersize',[20,10])
print(abc_im,'-dpdf','-painters','abc_im.pdf')

%% Tomography plots

tomographyplot(brocherizegrid(mean(anachain)),observed_data_tomography,'ana','ana.pdf')
tomographyplot(brocherizegrid(mean(abcchain)),observed_data_tomography,'abc','abc.pdf')
