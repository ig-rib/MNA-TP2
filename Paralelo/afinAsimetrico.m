%%Afin asimetrico de orden q
o=4;
q=floor(o/2);
c = parcluster;
c.NumWorkers = 4;
parpool('local', o);
N = 256;
x = linspace(-10,10,N);
Dx = x(2) - x(1);
Dk = 2*pi/(N*Dx);

k = [0:Dk:(N/2-1)*Dk,0,-(N/2-1)*Dk:Dk:-Dk];
c_1=13;
c_2 =3;

u = 1/2*c_1*(sech(sqrt(c_1)*(x+8)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x+1)/2)).^2;

Dt = 0.0001;
t=0;

U = fft(u);



tmax = 1.25; nplt = floor((tmax/100)/Dt); nmax = round(tmax/Dt);
%UData = zeros(length(u), ceil((nmax-40000)/nplt)+1);
%TData = zeros(ceil((nmax-40000)/nplt)+1);
%UData(:,1) = u';
UData = u'; TData = 0;
gammas4 = [2/3 2/3 -1/6 -1/6];

%yaux = zeros(1, length(u));
for i = 1:nmax
    t = i*Dt;
    
     y(1,:) = U;
    y(2,:) = U;
    y(3,:) = U;
    y(4,:) = U;
    
   %y=U;
%     spmd(o)
%         for j=1:o
%             if labindex==j
%                 for s = 1:labindex
%                     %y = lieTrotterPlus(Dt/labindex, y, k);
%                     %parte lineal
%                     y = y.*exp(1i*k.^3*Dt/labindex);
%     
%                     %%parte no lineal
%                     y = y - (3i*k*Dt/labindex).*fft(real(ifft(y)).^2);
%                 end
%             end
%         y = gammas4(labindex)*y;
%         end
%     end

    parfor ( j = 1:o)
         for s = 1:j
                    %y = lieTrotterPlus(Dt/labindex, y, k);
                    %parte lineal
                    y(j,:) = y(j,:).*exp(1i*k.^3*Dt/j);
    
                    %%parte no lineal
                    y(j,:) = y(j,:) - (3i*k*Dt/j).*fft(real(ifft(y(j,:))).^2);
         end
        y(j,:) = gammas4(j)*y(j,:);
        
    end
    for s = 2:o
      y(1, :) = y(1, :) + y(s, :);
    end

%     for z = 1:q
%         yaux = yaux + y{z};
%     end

%     y1 = y{1};
%     y2 = y{2};
%     y3 = y{3};
%     y4 = y{4};
%     U = y1+y2+y3+y4;
    U = y(1, :);

    if mod(i,round(nplt*4)) == 0
        u = real(ifft(U));
        %UData(:, i+1) = u'; TData(i+1) = t;
        %UData = [UData u']; TData = [TData t];
        plot(x,u,'LineWidth',2)
            axis([-10 10 0 10])
            xlabel('x')
            ylabel('u')
            text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
            drawnow
    end
%     if mod(i,nplt/8) == 0
%             plot(x,u,'LineWidth',2)
%             axis([-10 10 0 10])
%             xlabel('x')
%             ylabel('u')
%             text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
%             drawnow
%     end
end
delete(gcp('nocreate'))

figure

waterfall(x,TData(1:4:end),UData(:,1:4:end)')
xlabel x, ylabel t, axis([-10 10 0 tmax 0 10]), grid off
zlabel u



