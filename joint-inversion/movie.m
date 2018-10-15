% movie of moves

abcchain = csvread('abcchain.csv');
abc_tom_misfit = csvread('abc_tom_misfit.csv');

for ii = 2:5:5000
    
   subplot(1,2,1)
   
   colormap(brewermap([],'YlOrRd'));
   imagesc(brocherizegrid(reshape(abcchain(ii,:),[6,20])));
   
   subplot(1,2,2)
   
   plot(2:ii,abc_tom_misfit(2:ii))
   xlim([2,5000])
   ylim([0,6*10^5])
   
   set(gcf,'units','centimeters','position',[0,0,30,10],'papersize',[30,10])
   print('-dpng',['moves',filesep,'time',num2str(ii)]);
   clf
    
end
    
    