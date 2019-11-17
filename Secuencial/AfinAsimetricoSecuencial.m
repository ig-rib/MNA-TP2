%%Afin asimetrico de orden q
q=4;
N = 256;
x = linspace(-10,10,N);
Dx = x(2) - x(1);
Dk = 2*pi/(N*Dx);

k = [0:Dk:(N/2-1)*Dk,0,-(N/2-1)*Dk:Dk:-Dk];
c_1=13;
c_2 =3;

u = 1/2*c_1*(sech(sqrt(c_1)*(x+5)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x+3)/2)).^2;

Dt = 0.0001;
t=0;




U = fft(u);

tmax = 1.25; nplt = floor((tmax/100)/Dt); nmax = round(tmax/Dt);

UData = u'; TData = 0;
gammas4 = [2/3 2/3 -1/6 -1/6];

y = zeros(4, length(u));


for i = 1:nmax
    t = i*Dt;
    
    y(1,:) = U;
    y(2,:) = U;
    y(3,:) = U;
    y(4,:) = U;
        for j=1:q
                for s = 1:j
                    %y = lieTrotterPlus(Dt/labindex, y, k);
                    %parte lineal
                    y(j,:) = y(j,:).*exp(1i*k.^3*Dt/j);
    
                    %%parte no lineal
                    y(j,:) = y(j,:) - (3i*k*Dt/j).*fft(real(ifft(y(j,:))).^2);
                end
        y(j,:) = gammas4(j)*y(j,:);
        end
    for s = 2:q
       y(1,:) = y(1,:) + y(s,:); 
    end
    U = y(1,:);
    if mod(i,nplt) == 0
        u = real(ifft(U));
        UData = [UData u']; TData = [TData t];
        if mod(i,4*nplt) == 0
            plot(x,u,'LineWidth',2)
            axis([-10 10 0 10])
            xlabel('x')
            ylabel('u')
            text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
            drawnow
        end
    end
end

figure

waterfall(x,TData(1:4:end),UData(:,1:4:end)')
xlabel x, ylabel t, axis([-10 10 0 tmax 0 10]), grid off
zlabel u


