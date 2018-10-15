function [direction] = worstdatadirection(index,pars,data)
% returns -1 if data is too large and 1 if data is too low in the worst
% location
% indicates the direction the data needs to go

grav_data = data(1:20);

parameters = reshape(pars,[6,20]);

sim_data = forward_grav(parameters)+edge_effect_forward_grav([2750,0,2750]);

error = grav_data(index)-sim_data(index);

direction = error*(abs(1/error));