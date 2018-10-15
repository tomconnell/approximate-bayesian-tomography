function [range] = mcrange(true_model,edge_effect)

% Monte Carlo approximation for the range between summary statistics

[l,w] = size(true_model);

nSims = 10000;

% data
data = forward_grav(true_model)+edge_effect;
observed_summaries = gravsummaries(data,data);

summary_container = zeros(nSims,length(observed_summaries));
distance_container = zeros(nSims,length(observed_summaries));

for ii = 1:nSims
   
   model = (rand(120,1)*1500)+2000;
   model = reshape(model,[l,w]);
   
   % simualtion step   
   sim = forward_grav(model)+edge_effect;
   
   % compute summaries
   sim_summaries = gravsummaries(sim,data);
   
   summary_container(ii,:) = sim_summaries;
   
   distance_container(ii,:) = abs(sim_summaries-observed_summaries);
   
end

normalization_term = zeros(1,length(observed_summaries));

for ii = 1:length(normalization_term)
    % correct
    normalization_term(ii) = std(distance_container(:,ii));
    % incorrect
    % normalization_term(ii) = std(distance_container(ii,:));
end

range = normalization_term;