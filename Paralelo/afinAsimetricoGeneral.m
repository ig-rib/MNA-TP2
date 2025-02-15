function result = afinAsimetricoGeneral(q, func, xmin, xmax, N, tmax, Dt)
c = parcluster;
c.NumWorkers = q;
parpool('local', q);
tic
x = linspace(xmin,xmax,N);
Dx = x(2) - x(1);
Dk = 2*pi/(N*Dx);

k = [0:Dk:(N/2-1)*Dk,0,-(N/2-1)*Dk:Dk:-Dk];

u = func(x);

%Dt = 0.0001;
t=0;

U = fft(u);

%tmax = 1.25; 
nplt = floor((tmax/100)/Dt); nmax = round(tmax/Dt);
UData = u'; TData = 0;
%%gammas = [2/3 2/3 -1/6 -1/6];
gammas = gammasAsimetrico(q);
y = zeros(q, length(u));
for i = 1:nmax
    t = i*Dt;
   
    for j = 1:q
         y(j,:)=U;
    end
    
    parfor ( j = 1:q)
         for s = 1:j
                    %y = lieTrotterPlus(Dt/labindex, y, k);
                    %parte lineal
                   
                    y(j,:) = y(j,:).*exp(1i*k.^3*Dt/j);
    
                    %%parte no lineal
                    y(j,:) = y(j,:) - (3i*k*Dt/j).*fft(real(ifft(y(j,:))).^2);
         end
        y(j,:) = gammas(j)*y(j,:);
        
    end
    for s = 2:q
      y(1, :) = y(1, :) + y(s, :);
    end

    U = y(1, :);

    if mod(i,round(nplt*4)) == 0
        u = real(ifft(U));
        UData(:, i+1) = u'; TData(i+1) = t;
        UData = [UData u']; TData = [TData t];
        plot(x,u,'LineWidth',2)
            axis([-10 10 0 10])
            xlabel('x')
            ylabel('u')
            text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
            drawnow
    end

end
toc
delete(gcp('nocreate'))

figure

waterfall(x,TData(1:4:end),UData(:,1:4:end)')
xlabel x, ylabel t, axis([-10 10 0 tmax 0 10]), grid off
zlabel u
    result = [];
end