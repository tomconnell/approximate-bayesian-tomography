function [loglikelihood] = gravfun(parameters,data,toprint)

simulated_data = forward_grav(reshape(parameters,[4,8]));

edge_effect = edge_effect_forward_grav([2750,0,2750]);
simulated_data = simulated_data + edge_effect;

misfit = sum((simulated_data-data).^2);

if toprint
    % dlmwrite('analytical_misfit1.csv',misfit,'-append')
    dlmwrite('analytical_misfit.csv',misfit,'-append')
end

% sigma = error term
sigma = 2;

loglikelihood = misfit/sigma;