function [par] = jcwhichpartochange(highest_error)

highest_error = (highest_error*4)-3;

% which depth ? 1:4
% n = round((rand*4)+0.5);
% 
% par = 4*(highest_error-1)+n;

to_change = zeros(1,32);
to_change(highest_error:highest_error+3) = 1;

par = to_change;