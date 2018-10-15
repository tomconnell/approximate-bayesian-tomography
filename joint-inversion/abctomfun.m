function [loglikelihood,e] = abctomfun(parameters,data,tolerance)

simulated_data = tom(reshape(parameters,[6,20]),10,4)...
    + mean(randn(100,120)*sqrt(2));

total_loglikelihood = 0;

for ii = 1:12
    
    summaries_sim = tomsummaries(simulated_data((ii*10)-9:ii*10),data((ii*10)-9:ii*10));
    summaries_data = tomsummaries(data((ii*10)-9:ii*10),data((ii*10)-9:ii*10));

    distance = abs(summaries_data-summaries_sim);

    ll = gaussiankernel(distance,tolerance);
    
    total_loglikelihood = total_loglikelihood + ll;
    
end

error = abs(data-simulated_data);

e = sum(error.^2);
loglikelihood = total_loglikelihood;