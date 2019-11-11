%%tengo el vector k con todos los 
%%2*pi/N/
function xnew = lieTrotterMinus(h, X, k)
    
        
    %%parte no lineal
    X = X - (3i*k*h).*fft(real(ifft(X)).^2);
    %%parte lineal
    xnew = X.*exp(1i*k.^3*h);

end