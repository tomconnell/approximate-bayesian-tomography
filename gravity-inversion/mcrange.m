function [range] = mcrange(data)

% Monte Carlo approximation for the range between summary statistics

nSims = 10000;

% edge_effect = edge_effect_forward_grav([2750,0,2750]);

% data
% data = forward_grav(true_model)+edge_effect;
observed_summaries = gravsummaries(data,data);

summary_container = zeros(nSims,length(observed_summaries));
distance_container = zeros(nSims,length(observed_summaries));

for ii = 1:nSims
   
   model = (rand(32,1)*3500)+2000;
   model = reshape(model,[4,8]);
   
   edge_effect = edge_effect_forward_grav([2750,0,2750]);
   
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