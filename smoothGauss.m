function yfilt = smoothGauss(y)
sigma = 0.5;
sz = 13;    % length of gaussFilter vector
yext = [y(end-(sz-1)/2+1:end); y; y(1:(sz-1)/2)];
x = linspace(-(sz-1)/2,(sz-1)/2,sz);
gaussFilter = exp(-x.^2/(2*sigma^2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
yfilt = conv(yext,gaussFilter,'same');
yfilt = yfilt((sz-1)/2+1:end-(sz-1)/2);
end