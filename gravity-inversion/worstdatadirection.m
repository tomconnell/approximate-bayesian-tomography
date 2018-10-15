function [direction] = worstdatadirection(index,pars,data)
% returns -1 if data is too large and 1 if data is too low in the worst
% location
% indicates the direction the data needs to go

parameters = reshape(pars,[4,8]);

sim_data = forward_grav(parameters);

error = data(index)-sim_data(index);

direction = error*(abs(1/error));