function [summaries] = variablesummaries(simulation,data,option)

% x = [((1:length(simulation)).^2)',(1:length(simulation))',ones(length(simulation),1)]\simulation';
 x = [(1:length(simulation))',ones(length(simulation),1)]\simulation';
residual = sum(abs(simulation-data));

sigma = std(simulation);

m = mean(simulation);

% The lot
if option == 1
    summaries = [x(1),x(2),residual,sigma,m];
end

% Residual and marginals
if option == 2
    summaries = [residual,sigma,m];
end

% Orginal
if option == 3
    summaries = [x(1),x(2),residual,sigma];
end

% Residual and straight line
if option == 4
        summaries = [x(1),x(2),residual];
end

% No residual
if option == 5
        summaries = [x(1),x(2),sigma,m];
end

