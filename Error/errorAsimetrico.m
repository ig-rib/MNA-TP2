function result = errorAsimetrico(q, func, exactF, xmin, xmax, N, tmax, Dt, wGraph)

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
gammas = gammasAsimetrico(q);
y = zeros(q, length(u));
Error = [];
for i = 1:nmax
    t = i*Dt;
   
    for j = 1:q
         y(j,:)=U;
    end
    
    uorig = exactF(x, t);
    
    for j = 1:q
         for s = 1:j
                    %Lie Trotter con FFT 
                    %parte lineal
                   
                    y(j,:) = y(j,:).*exp(1i*k.^3*Dt/j);
    
                    %%parte no lineal
                    y(j,:) = y(j,:) - (3i*k*Dt/j).*fft(real(ifft(y(j,:))).^2);
                    

%                     %Strang con FFT
%                         % lineal
%                     y(j,:) = y(j,:).*exp(1i*k.^3*Dt/j/2);
% 
%                     % no lineal
% 
%                     y(j,:) = y(j,:)  - (3i*k*Dt/j).*fft((real(ifft(y(j,:)))).^2);
% 
%                     % lineal 
% 
%                     y(j,:) = y(j,:).*exp(1i*k.^3*Dt/j/2);
         end
        y(j,:) = gammas(j)*y(j,:);
        
    end
    for s = 2:q
      y(1, :) = y(1, :) + y(s, :);
    end

    U = y(1, :);

    if mod(i,round(nplt)) == 0
        u = real(ifft(U));
        if mod(i,round(nplt*4))==0
            UData = [UData u']; TData = [TData t];
        end
        
        Error = [Error abs(u-uorig)'];
        
        if wGraph
            plot(x,u, x, uorig,x, Error(:, round(i/nplt))'),
                legend('Solucion Numerica', 'Solucion Analitica', 'Error', 'Location', 'southoutside');
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