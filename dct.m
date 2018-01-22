function [numcoef,rmse,y2,featfrac]=dct(y,perc)
% apply dct
set(0,'RecursionLimit',10000);
X = dct(y);

[XX, ind] = sort(X.^2,'descend'); %squaring for energy measure

totalE = sum(XX);
partialE = cumsum(XX);

%find the index where the remaining energy fraction is greater than perc
numcoef = find(partialE/totalE >= perc,1); 

%set lowest coefficients to zero and reconstruct the signal
X(ind(numcoef+1:end)) = 0;
y2 = idct(X);
rmse = sqrt(mean((y-y2).^2));
featfrac = round(numcoef/length(y)*1000)/10; %in new matlab versions, you can use 100*round(numcoef/length(y),3)

end