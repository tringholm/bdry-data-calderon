function out = c(x,y)
%out = 1*( 2 + cos(10*(x - y))); % Oscillating case
out = 3./(1+exp(2*(x + y))); % Non-oscillating case
end