function [summaries] = gravsummaries(simulation,data)

% x = [((1:length(simulation)).^2)',(1:length(simulation))',ones(length(simulation),1)]\simulation';
 x = [(1:length(simulation))',ones(length(simulation),1)]\simulation';
residual = sum(abs(simulation-data));

sigma = std(simulation);

m = mean(simulation);

% summaries = [x(1),x(2),x(3),residual,sigma];
summaries = [x(1),x(2),residual,sigma,m];