function [probability] = skewnormal(x,location,scale,shape)

xi = location;
omega = scale;
alpha = shape;

% anonymous function names based on what appears on the wikipedia page for
% "skew normal distribution"
upphi = @(x) 0.5*(1+erf(x/sqrt(2)));
lowphi = @(x) (1/sqrt(2*pi))*exp(-(x^2)/2);

probability = (2/omega)*lowphi((x-xi)/omega)*upphi(alpha*((x-xi)/omega));