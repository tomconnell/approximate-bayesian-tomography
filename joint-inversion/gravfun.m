function [loglikelihood,e] = gravfun(parameters,data,edge_effect)

simulated_data = forward_grav(reshape(parameters,[6,20]))+edge_effect;

misfit = sum((simulated_data-data).^2);

% sigma = error term
sigma = 2;
e = misfit;
loglikelihood = misfit/sigma;