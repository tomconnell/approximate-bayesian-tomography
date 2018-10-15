function [summaries] = tomsummaries(simulation,data)

% x = [((1:length(simulation)).^2)',(1:length(simulation))',ones(length(simulation),1)]\simulation';

 x = [(1:length(simulation))',ones(length(simulation),1)]\simulation';

residual = sum(abs(simulation-data));

summaries = [x(1),x(2),residual];