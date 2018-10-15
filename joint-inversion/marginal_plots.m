% marginal plots

abcchain = csvread('abcchain.csv');
anachain = csvread('anachain.csv');

par2plot = round(linspace(1,120,32));

colors = brewermap(2,'set1');

true_model_density = reshape(linspace(2000,3500,120),[6,20]); % units kg/m3
true_model = true_model_density/1000;

figure
for ii = 1:32
    
    
    subplot(8,4,ii) 
    
    h = histogram(anachain(:,par2plot(ii))/1000,30,'normalization','pdf',...
        'facecolor',colors(1,:),'facealpha',0.3);
    hold on
    histogram(abcchain(:,par2plot(ii))/1000,30,'normalization','pdf',...
        'facecolor',colors(2,:),'facealpha',0.3)
    
    topoff = max(h.Values)+0.3*max(h.Values);
    
    plot(linspace(true_model(par2plot(ii)),true_model(par2plot(ii))),...
        linspace(0,topoff),'k--','linewidth',2)
    
    ylim([0,topoff])
    xlim([2,3.5])
    set(gca,'Ytick',[])
    xlabel(['Parameter ',num2str(par2plot(ii))])
    
end

set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print(gcf,'-dpdf','-painters','joint_hist.pdf')

figure
pars = round(linspace(2,118,8));
for ii = 1:8
    
    par2plot = pars(ii);
    
    subplot(16,8,[(((ii*2)-1)*8)-7,(((ii*2)-1)*8)-6,((ii*2)*8)-7,((ii*2)*8)-6]);
    
    h = histogram(anachain(:,par2plot),30,'normalization','pdf',...
        'facecolor',colors(1,:),'facealpha',0.3);
    hold on
    histogram(abcchain(:,par2plot),30,'normalization','pdf',...
        'facecolor',colors(2,:),'facealpha',0.3);
    
    topoff = max(h.Values)+0.3*max(h.Values);
    
    plot(linspace(true_model_density(par2plot),...
        true_model_density(par2plot)),linspace(0,topoff),'k')
    
    ylim([0,topoff])
    
    xlabel(['Parameter ',num2str(par2plot)])
    
    set(gca,'view',[-90,90],'xlim',[2000,3500]...
        ,'xtick',[],'ytick',[])
    
    subplot(16,8,[(((ii*2)-1)*8)-5:((ii*2)-1)*8,((ii*2)*8)-5:(ii*2)*8])
    
    plot(anachain(:,par2plot),'color',colors(1,:))
    hold on
    plot(abcchain(:,par2plot),'color',colors(2,:))
    plot(linspace(1,length(anachain)),linspace(true_model_density(par2plot),...
        true_model_density(par2plot)),'k-')
    ylim([2000,3500])
    set(gca,'Xtick',[],'Ytick',[])
    
end
set(gcf,'units','centimeters','position',[0,0,25,30],'papersize',[25,30])
print(gcf,'-dpdf','-painters','joint_marginals.pdf')




