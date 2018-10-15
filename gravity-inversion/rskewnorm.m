function [sample] = rskewnorm(location,scale,shape)

% work out what the approx highest probability will be
x = linspace(-scale,scale);
y = zeros(100,1);
for ii = 1:100
    y(ii) = skewnormal(x(ii),location,scale,shape);
end

% highest probability the skew normal pdf defined by skewnromal.m will
% reach
max_prob = max(y);

accept = 1;

while accept

    % monte carlo generate sample in the range
    x = (rand*10*scale)-(5*scale);

    % accept that mc sample with probability = to skew normal dist at
    % location x
    sim = skewnormal(x,location,scale,shape);

    if sim > rand*max_prob
        
        accept = 0;
        
    end
    
end

sample = x;
    
    