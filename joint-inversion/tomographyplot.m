function [] = tomographyplot(parameters,observed_data,method,file_name)
clf
tom_image = reshape(parameters,[6,20]);

subplot(6,3,1:6)
colormap(brewermap([],'YlOrRd'));
imagesc(tom_image)
colorbar()
set(gca,'ytick',[],'xtick',[])

colors = brewermap(2,'set1');
if strcmp(method,'ana')
    c_to_use = colors(1,:);
end
if strcmp(method,'abc')
    c_to_use = colors(2,:);
end

simulated_data = tom(tom_image,10,4);

for ii = 1:12
    
    subplot(6,3,6+ii)
   
    
    plot(1:10,observed_data((ii*10)-9:ii*10),...
        'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
    hold on
    plot(1:10,simulated_data((ii*10)-9:ii*10)...
        ,'x','MarkerEdgeColor',c_to_use,'MarkerFaceColor',c_to_use)
    
    
    plot(simulated_data((ii*10)-9:ii*10),'color',c_to_use)
   
    plot(observed_data((ii*10)-9:ii*10),'color','k')

    
    xlabel(['source ',num2str(ii)])
end

set(gcf,'units','centimeters','position',[0,0,20,30],'papersize',[20,30])
print('-dpdf','-painters',file_name)
clf
