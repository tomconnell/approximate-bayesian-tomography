function [range] = mcrange2(true_model,option)

% Monte Carlo approximation for the range between summary statistics

nSims = 10000;

% data
data = forward_grav(true_model);
observed_summaries = variablesummaries(data,data,option);

summary_container = zeros(nSims,length(observed_summaries));
distance_container = zeros(nSims,length(observed_summaries));

for ii = 1:nSims
   
   model = (rand(32,1)*3000)-1500;
   model = reshape(model,[4,8]);
   
   % simualtion step   
   sim = forward_grav(model);
   
   % compute summaries
   sim_summaries = variablesummaries(sim,data,option);
   
   summary_container(ii,:) = sim_summaries;
   
   distance_container(ii,:) = abs(sim_summaries-observed_summaries);
   
end

normalization_term = zeros(1,length(observed_summaries));

for ii = 1:length(normalization_term)
    normalization_term(ii) = std(distance_container(:,ii));
end

range = normalization_term;