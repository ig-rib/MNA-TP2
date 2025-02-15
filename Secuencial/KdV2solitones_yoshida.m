%Yoshida 4
%Kdv Dos_solitones
function result = KdV2solitones_yoshida(func, xmin, xmax, N, tmax, Dt, wGraph)


x = linspace(xmin,xmax,N);
Dx = x(2) - x(1);
delta_k = 2*pi/(N*Dx);

k = [0:delta_k:(N/2-1)*delta_k,0,-(N/2-1)*delta_k:delta_k:-delta_k];
c_1=13;
c_2 =3;

%los coeficientes de Yoshida
%https://www.asc.tuwien.ac.at/~winfried/splitting/index.php?rc=0&ab=y-44-ab&name=Y%204-4
w1 = 1/(2-2^(1/3));
w0 = 1-2*w1;


u = func(x);

%u = ;

t=0;
plot(x,u,'LineWidth',2)
axis([-10 10 0 10])
xlabel('x')
ylabel('u')
text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',14)
drawnow

nplt = floor((tmax/100)/Dt); nmax = round(tmax/Dt);
udata = u'; tdata = 0;
U = fft(u);

for n = 1:nmax
    t = n*Dt;
    
    % lineal
    U = U.*exp(1i*k.^3*Dt*w1/2);
    
    % no lineal
    
    U = U  - (3i*k*Dt*w1).*fft((real(ifft(U))).^2);
    
    % lineal
    U = U.*exp(1i*k.^3*Dt*(w1+w0)/2);
    
    % no lineal
    
    U = U  - (3i*k*Dt*w0).*fft((real(ifft(U))).^2);
    
    % lineal
    U = U.*exp(1i*k.^3*Dt*(w1+w0)/2);
    
    % no lineal
    
    U = U  - (3i*k*Dt*w1).*fft((real(ifft(U))).^2);
    
    % lineal
    U = U.*exp(1i*k.^3*Dt*w1/2);
    
    % no lineal
    
    %U = U  - (3i*k*delta_t*w1/2).*fft((real(ifft(U))).^2);
    
    
    if mod(n,round(nplt/8)) == 0
        u = real(ifft(U));
        udata = [udata u']; tdata = [tdata t];
        if wGraph
            if mod(n,round(nplt/8)) == 0
                plot(x,u,'LineWidth',2)
                axis([xmin xmax 0 10])
                xlabel('x')
                ylabel('u')
                text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
                drawnow
            end
        end
    end
end

if wGraph
    figure

    waterfall(x,tdata(1:4:end),udata(:,1:4:end)')
    xlabel x, ylabel t, axis([xmin xmax 0 tmax 0 10]), grid off
    zlabel u
end
end
