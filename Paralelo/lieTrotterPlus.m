%%tengo el vector k con todos los 
%%2*pi/N/
function xnew = lieTrotterPlus(h, X, k)
    
    %%parte lineal
    X = X.*exp(1i*k.^3*h);
    
    %%parte no lineal
    xnew = X - (3i*k*h).*fft(real(ifft(X)).^2);

end