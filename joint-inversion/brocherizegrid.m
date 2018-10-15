function [vp_grid] = brocherizegrid(density_grid)
% function to take in a density grid and covert it to a Vp (compressional
% wave speed) grid based on the empirical relationship defined by a
% 5th order polynomial by Brocher (2005)
% This relationship is valid for densities between 2 and 3.5 g/cm^3

density_grid = density_grid/1000;

% From Brocher (2005)
a = 0.8228;
b = -9.1819;
c = 37.083;
d = -63.064;
e = 39.128;

D = density_grid;

% Compute compressional wavespeed grid 
vp_grid = a*(D.^5) + b*(D.^4) + c*(D.^3) + d*(D.^2) + e*D;