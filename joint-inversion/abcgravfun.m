function [loglikelihood,e] = abcgravfun(parameters,data,edge_effect,tolerance)

reshaped_pars = reshape(parameters,[6,20]);

simulated_data = forward_grav(reshaped_pars) + edge_effect;

simulated_data = simulated_data + mean(randn(100,20)*sqrt(2));

summaries_sim = gravsummaries(simulated_data,data);
summaries_data = gravsummaries(data,data);

distance = abs(summaries_data-summaries_sim);

error = abs(data-simulated_data);

e = sum(error.^2);
loglikelihood = gaussiankernel(distance,tolerance);