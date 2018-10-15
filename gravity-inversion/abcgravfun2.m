function [loglikelihood] = abcgravfun2(parameters,data,tolerance,toprint,option)

reshaped_pars = reshape(parameters,[4,8]);

simulated_data = forward_grav(reshaped_pars);

simulated_data = simulated_data+(randn(1,length(simulated_data))*sqrt(2));

summaries_sim = variablesummaries(simulated_data,data,option);
summaries_data = variablesummaries(data,data,option);

distance = abs(summaries_data-summaries_sim);

error = abs(data-simulated_data);

if toprint
    dlmwrite('abc_misfit.csv',sum(error.^2),'-append')
end

loglikelihood = gaussiankernel(distance,tolerance);