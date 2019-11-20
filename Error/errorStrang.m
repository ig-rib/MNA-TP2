function result = errorStrang(func, exactF, xmin, xmax, N, tmax, Dt, wGraph)

x = linspace(xmin,xmax,N);
Dx = x(2) - x(1);
Dk = 2*pi/(N*Dx);

k = [0:Dk:(N/2-1)*Dk,0,-(N/2-1)*Dk:Dk:-Dk];

u = func(x);
%Dt = 0.0001;

U = fft(u);

%tmax = 1.25; 
nplt = floor((tmax/100)/Dt); nmax = round(tmax/Dt/5);
UData = u'; TData = 0;
%%gammas = [2/3 2/3 -1/6 -1/6];
Error = [];
for i = 1:nmax
    t = i*Dt;
   
    uorig = exactF(x, t);
    
    U = U.*exp(1i*k.^3*Dt/2);
    U = U  - (3i*k*Dt).*fft((real(ifft(U))).^2);
    U = U.*exp(1i*k.^3*Dt/2);  

    if mod(i,round(nplt)) == 0
        u = real(ifft(U));
        if mod(i,round(nplt*4))==0
            UData = [UData u']; TData = [TData t];
        end
        
        Error = [Error abs(u-uorig)'];
        
        if wGraph
            plot(x,u, x, uorig,x, Error(:, round(i/nplt))'),
                axis([xmin xmax 0 10])
                xlabel('x')
                ylabel('u')
                text(xmin + (xmax-xmin)/2,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
                text(xmin + (xmax-xmin)/2, 8.5, ['Dt = ', num2str(Dt, '%1.5g')], 'FontSize', 10);
                text(xmin + (xmax-xmin)/2, 8, ['Dx = ', num2str(Dx, '%1.5g')], 'FontSize', 10);
                text(xmin + (xmax-xmin)/2, 7.5, ['Dt/Dx = ', num2str(Dt/Dx, '%1.5g')], 'FontSize', 10);
                drawnow
        end
    end

end


if wGraph
figure
    plot(max(Error)), hold on
    plot(mean(Error)), 
    legend('Error Global (Maximo Valor Absoluto)', 'Media del Error', 'Location', 'southoutside'),
    xlabel(['Pasos Temporales (Dt = ', num2str(Dt, '%1.5g'), ')']);
end
    result = UData;
end