function [loglikelihood] = abcgravfun(parameters,data,tolerance,toprint)

reshaped_pars = reshape(parameters,[4,8]);

simulated_data = forward_grav(reshaped_pars);
edge_effect = edge_effect_forward_grav([2750,0,2750]);
simulated_data = simulated_data + edge_effect;
simulated_data = simulated_data+mean(randn(100,8)*sqrt(2));

summaries_sim = gravsummaries(simulated_data,data);
summaries_data = gravsummaries(data,data);

distance = abs(summaries_data-summaries_sim);

error = abs(data-simulated_data);

if toprint
    dlmwrite('abc_misfit.csv',sum(error.^2),'-append')
end

loglikelihood = gaussiankernel(distance,tolerance);
