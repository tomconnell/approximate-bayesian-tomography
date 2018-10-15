function [indicies] = blockselection()

% This is a script to randomly select blocks in the 2-dimensional 4x8 
% subsurface to update in an MCMC algorithm
% Each individual block is 2x2

n = round((rand*4)+0.5);

n = (n*8)-7;

if rand < 0.5
    n = n +2;
end

i = zeros(1,32);
i(n) = 1;
i(n+1) = 1;
i(n+4) = 1;
i(n+5) = 1;
indicies = i;