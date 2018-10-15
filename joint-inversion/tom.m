function [data] = tom(slowness,n_receivers,n_sources)
% A simple tomography example
% The rays are straight lines


%% Establish grid dimensions (0 to width of length width)

[l,w] = size(slowness);

% x-direction
x_min = 0;
x_max =  l;
dx = (x_max - x_min)/(l - 1);
Dx = x_max-x_min;

% y-direction
ymin = 0;
ymax = w;
dy=(ymax-ymin)/(w-1);
Dy = ymax-ymin;


%% Rearrange the model into a vector
% By rows not matlab default columns

holder = slowness';
mtrue = holder(:);


%% Define receiver locations

% Container for combined receivers and sources
xs = zeros(n_receivers + 3*n_sources,2);

% Set receiver locations at the surface
xs(1:n_receivers,1) = x_min; 
xs(1:n_receivers,2) = linspace(1,w-1,n_receivers);

% Set source locations at the base
xs(n_receivers+1:n_receivers+n_sources,1) = 3*x_max;
xs(n_receivers+1:n_receivers+n_sources,2) = ...
    ymin + Dy*(1:n_sources)'/(n_sources+1);

% Set source locations on the left side
xs(n_receivers+n_sources+1:n_receivers+(2*n_sources),1) = ...
    linspace(l*0.5,x_max*3,n_sources);
xs(n_receivers+n_sources+1:n_receivers+(2*n_sources),2) = -5;

% Set source locations on the right side
xs(n_receivers+(2*n_sources)+1:n_receivers+(3*n_sources),1) = ...
    linspace(l*0.5,x_max*3,n_sources);
xs(n_receivers+(2*n_sources)+1:n_receivers+(3*n_sources),2) = w+5;


%% Source reciever pairs

% sets paths between receivers for every receiver
% ray [py,px,qy,qx] not sure if p or q is source/receiver

% Number of rays
N = 3*n_sources*n_receivers;
ray = zeros(N,4);

% Ray counter
k = 1;
% Loop to refine all source/receiver pairs possible
for ii = 1:n_sources*3
    for jj = 1:n_receivers
        
        ray(k,1) = xs(jj,1);
        ray(k,2) = xs(jj,2);
        ray(k,3) = xs(ii+n_receivers,1);
        ray(k,4) = xs(ii+n_receivers,2);
        
        k = k + 1;
        
    end
end

% colormap(brewermap([],'RdBu'))
% imagesc(slowness)
% hold on
% for k = 1:N
%     plot( [ray(k,2) ray(k,4)], [ray(k,1), ray(k,3)], 'k-', 'LineWidth', 1 );
% end

%% Build data kernel

% order or model parameters in vector: index = (ix-1)*Ny + iy
M = w*l;
G=spalloc(N,M,2*N*w);
Nr = 2000;
bins= 1:w*l;
for k = 1:N
    x1 = ray(k,1);
    y1 = ray(k,2);
    x2 = ray(k,3);
    y2 = ray(k,4);
    r = sqrt( (x2-x1)^2 + (y2-y1)^2 );
    dr = r/Nr;
% I use a sloppy way of computing the length of the ray
% in each pixel.  I subdivide the ray into Nr pieces, and
% assign each piece to exactly one pixel, the one its
% closest to

if 0 % slow but sure way
    for i = 1:Nr
        x = x1 + (x2-x1)*i/Nr;
        y = y1 + (y2-y1)*i/Nr;
        ix = 1+floor( (Nx-1)*(x-xmin)/Dx );
        iy = 1+floor( (Ny-1)*(y-ymin)/Dy );
        q = (ix-1)*Ny + iy;
        G(k,q) = G(k,q) + dr;
    end
else % faster way, or so we hope
    % calculate the array indices of all the ray pieces
    xv = x1 + (x2-x1)*[1:Nr]'/Nr;
    yv = y1 + (y2-y1)*[1:Nr]'/Nr;
    ixv = 1+floor( (w-1)*(xv-x_min)/Dx );
    iyv = 1+floor( (l-1)*(yv-ymin)/Dy );
    qv = (ixv-1)*l + iyv;
    % now count of the ray segments in each pixel of the
    % image, and use the count to increment the appropriate
    % element of G.  The use of the hist() function to do
    % the counting is a bit weird, but it seems to work
    count=hist(qv,bins); 
    icount = find( count~=0 );
    G(k,icount) = G(k,icount) + count(icount)*dr;
end
end



%% compute true travel times
d = G*mtrue;

data = d';


% figure
% 
% subplot(4,4,1)
% plot(1:10,d(11:20),...
%         'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
% hold on
% plot(1:10,d(11:20),'k-')
% title('Source 2')
% ylabel('\Delta t')
% xlabel('receiver')
% set(gca,'ytick',[])
% 
% subplot(4,4,2)
% plot(1:10,d(41:50),...
%         'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
% hold on
% plot(1:10,d(41:50),'k-')
% title('Source 5')
% ylabel('\Delta t')
% xlabel('receiver')
% set(gca,'ytick',[])
% 
% subplot(4,4,3)
% plot(1:10,d(71:80),...
%         'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
% hold on    
% plot(1:10,d(71:80),'k-')
% title('Source 8')
% ylabel('\Delta t')
% xlabel('receiver')
% set(gca,'ytick',[])
% 
% subplot(4,4,4)
% plot(1:10,d(101:110),...
%         'marker','^','MarkerEdgeColor','k','MarkerFaceColor','k')
% hold on
% plot(1:10,d(101:110),'k-')
% title('Source 11')
% ylabel('\Delta t')
% xlabel('receiver')
% set(gca,'ytick',[])
% 
% 
% subplot(4,4,5:16)
% colormap(brewermap([],'YlOrRd'))
% imagesc(slowness)
% hold on
% for k = 1:N
%     plot( [ray(k,2) ray(k,4)], [ray(k,1), ray(k,3)], 'k-', 'LineWidth', 1 );
% end
% ylim([0.5,6])
% xlim([0.5,19.5])
% colorbar()
% set(gca,'ytick',[],'xtick',[])
% 
% set(gcf,'units','centimeters','position',[0,0,20,15],'papersize',[20,15])
% print('-dpdf','-painters','observed_tomography.pdf')
