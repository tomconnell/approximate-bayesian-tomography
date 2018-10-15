function [range] = tomrange(true_model)

% Monte Carlo approximation for the range between summary statistics

[l,w] = size(true_model);

nSims = 1000;

% data
data = tom(brocherizegrid(true_model),10,4);

observed_summaries = zeros(12,3);

for ii = 1:12
    
    observed_summaries(ii,:) = ...
        tomsummaries(data((ii*10)-9:ii*10),data((ii*10)-9:ii*10));
    
end

summary_container = zeros(nSims*12,3);
distance_container = zeros(nSims*12,3);

for ii = 1:nSims
   
   model = (rand(120,1)*6.6)+2;
   model = reshape(model,[l,w]);
   
   % simualtion step   
   sim = tom(model,10,4);
   
   for jj = 1:12
       
       % compute summaries
       sim_summaries = tomsummaries(sim((jj*10)-9:jj*10),data((jj*10)-9:jj*10));
   
       summary_container((ii*12)-11+jj,:) = sim_summaries;
   
       distance_container((ii*12)-11+jj,:) = ...
           abs(sim_summaries-observed_summaries(jj,:));
       
   end
   
end

normalization_term = zeros(1,3);

for ii = 1:length(normalization_term)
    normalization_term(ii) = std(distance_container(:,ii));
end

range = normalization_term;