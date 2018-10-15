function [indicies] = blockselection()

% This is a script to randomly select blocks in the 2-dimensional 4x8 
% subsurface to update in an MCMC algorithm
% Each individual block is 2x2

n = round((rand*10)+0.5);

n = (n*12)-11;

sim = rand;

if sim <= 0.33
    n = n + 2;
end

if sim > 0.66
    n = n + 4;
end

i = zeros(1,120);
i(n) = 1;
i(n+1) = 1;
i(n+6) = 1;
i(n+7) = 1;
indicies = i;