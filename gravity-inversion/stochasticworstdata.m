function [index] = stochasticworstdata(pars,data)

parameters = reshape(pars,[4,8]);

sim_data = forward_grav(parameters);

edge_effect = edge_effect_forward_grav([2750,0,2750]);
sim_data = sim_data + edge_effect;

l2error = (sim_data-data).^2;

max_error = max(l2error);

indicator = 1;

index_holder = 0;

while indicator
    
    trial = round((rand*length(data))+0.5);
    
    if rand*max_error < l2error(trial)
        
        indicator = 0;
        
        index_holder = trial;
        
    end
    
end

index = index_holder;