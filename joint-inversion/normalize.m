function [normalized] = normalize(input)
 
normalized = (input-min(input))/(max(input)-min(input));
