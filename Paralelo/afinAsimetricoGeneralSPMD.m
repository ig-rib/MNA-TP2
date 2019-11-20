function result = afinAsimetricoGeneralSPMD(q, func, xmin, xmax, N, tmax, Dt, wGraph)

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
y = [];
for i = 1:nmax
    t = i*Dt;
    y = U;
    spmd(q)
        for j = 1:q
             if labindex == j
                 for s=1:labindex
                        %Lie Trotter con FFT 
                        y = y.*exp(1i*k.^3*Dt/labindex);
                        y = y - (3i*k*Dt/labindex).*fft(real(ifft(y)).^2);
                 end
             end
            y = gammas(labindex)*y;
        end
    end
    for s = 2:q
      y{1} = y{1} + y{s};
    end
    U = y{1};
    
    %plot    
    plot(real(ifft(U)));
    if mod(i,round(nplt)) == 0
        u = real(ifft(U));
        if mod(i,round(nplt*4))==0
            UData = [UData u']; TData = [TData t];
        end
        if wGraph
            plot(x,u,'LineWidth',2)
                axis([xmin xmax 0 10])
                xlabel('x')
                ylabel('u')
                text(xmin + (xmax-xmin)/2,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
                text(xmin + (xmax-xmin)/2, 8, ['dx = ', num2str(Dx, '%1.5g')], 'FontSize', 10);
                text(xmin + (xmax-xmin)/2, 7, ['Dt/Dx = ', num2str(Dt/Dx, '%1.5g')], 'FontSize', 10);
                drawnow
        end
    end

end

toc
delete(gcp('nocreate'))

if wGraph
figure
    waterfall(x,TData(1:4:end),UData(:,1:4:end)')
    xlabel x, ylabel t, axis([xmin xmax 0 tmax 0 10]), grid off
    zlabel u
end
    result = UData;
end