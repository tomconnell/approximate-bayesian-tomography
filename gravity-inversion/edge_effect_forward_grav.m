function [simulated_data] = edge_effect_forward_grav(simulated_parameters)

% simulated_data = forward_grav(simulated_parameters)
% This function solves the forward gravity problem, vertical gravity
% anomaly for a 2D discretized subsurface for a set model domain
% it uses the density assigned in the simulated_parameters argument
% the output is a set of equidistant measurements across the fixed domain

% This code uses Talwani's equations (via gpoly.m) to compute the
% contribution from each individual rectangle

% Written by Tom Connell
% 2017
% Masters of Research (year 1)
%% define mesh and observation points
% MESH

height=4e4;    % rectangles vertical dimension (m)
width=(1.6e5);    % rectangles horizontal dimension (m)
nx = 3;    % number of rectangles in the x-direction
ny = 1;    % number of rectangles in the y-direction
ntot = nx*ny;   % the total number of rectangles
domain_width = nx*width;   % total width of the model domain
% domain_height = ny*height;  % total height of the model domain

% OBSERVATION PROFILE
number_of_observations = 8;
obs_points = linspace((1.6e5)+10000,(2*(1.6e5))-10000,number_of_observations);
topography = zeros(1,number_of_observations);

% initialise matrix containing the total effect of each pixel at each
% observation point
gtotal = zeros(1,number_of_observations);

%% Main loop
for obs_location = 1:number_of_observations     % start loop over observations
    current_gtotal = 0;     % initialise holder for gravity @ current obsevation 
    
    for current_y_location = 1:ny       % loop over pixels by depth
        for current_x_location = 1:nx       % loop over pixels by width
            
            % Set co-ordinates of current pixel
            % [(j-1)*Dh  j*Dh  j*Dh  (j-1)*Dh (j-1)*Dh];
            xcorn = [(current_x_location-1)*width, current_x_location*width, current_x_location*width, (current_x_location-1)*width, (current_x_location-1)*width];
            zcorn = [(current_y_location-1)*height, (current_y_location-1)*height, (current_y_location)*height, (current_y_location)*height, (current_y_location-1)*height];
            
            g = gpoly(obs_points(obs_location),topography(obs_location),xcorn, zcorn, 4, simulated_parameters(current_y_location,current_x_location));
            
            current_gtotal = current_gtotal+g;
            
        end         % end loop over pixels by width        
    end         % end loop over pixels by depth
    
    gtotal(obs_location) = current_gtotal;
    
end     % end loop over observation

%csvwrite('simulated_observations.dat', gtotal);
simulated_data = gtotal;