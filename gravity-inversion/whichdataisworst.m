function [index] = whichdataisworst(pars,data)

parameters = reshape(pars,[4,8]);

sim_data = forward_grav(parameters);

error = abs(sim_data-data);

[~,index] = max(error);