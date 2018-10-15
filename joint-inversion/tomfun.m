function [loglikelihood,e] = tomfun(parameters,data)

simulated_data = tom(reshape(parameters,[6,20]),10,4);

misfit = sum((simulated_data-data).^2);

% sigma = error term
sigma = 2;
e = misfit;
loglikelihood = misfit/sigma;