%%tests

c_1=13;
c_2 =3;
func1 = @(x)(1/2*c_1*(sech(sqrt(c_1)*(x+3)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x+1)/2)).^2);

%afinAsimetricoGeneralSecuencial(2, func1, -5, 5, 256, 5, 0.00001, 1);



% for i=4:7
% UData1 = afinAsimetricoGeneralSecuencial(i, func1, -5, 5, 128, 1, 0.0001, 1);
% UData2 =  afinAsimetricoGeneralSecuencial(i+1, func1, -5, 5, 128, 1, 0.0001, 1);
% %UData2 = UData2(1:2:2^(7+i), :);
% Error = max(abs(UData1-UData2));
% figure
% plot(Error),
% %%loglog(Error);
% end

Dt = 0.0001;
tmax = 1;
xmin = -8;
xmax = 8;
N = 128;
x = linspace(xmin,xmax,N);
Dx = x(2)-x(1);
nmax = round(tmax/Dt);
nplt = floor((tmax/100)/Dt);

func1 = @(x)(1/2*c_1*(sech(sqrt(c_1)*(mod(x+3, xmax-xmin)+xmin)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(mod(x+1,xmax-xmin)+xmin)/2)).^2);

% UOrig = [];
% figure
% for i = 0:nplt:nmax
%     u=(1/2*c_1*(sech(sqrt(c_1)*(mod(x+3-c_1*(i*Dt), xmax-xmin)+xmin)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(mod(x+1-c_2*(i*Dt), xmax-xmin)+xmin)/2)).^2);
%      plot(x,u,'LineWidth',2)
%                 axis([xmin xmax 0 10])
%                 xlabel('x')
%                 ylabel('u')
%                 text(xmin + (xmax-xmin)/2,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
%                 text(xmin + (xmax-xmin)/2, 8, ['dx = ', num2str(Dx, '%1.5g')], 'FontSize', 10);
%                 text(xmin + (xmax-xmin)/2, 7, ['Dt/Dx = ', num2str(Dt/Dx, '%1.5g')], 'FontSize', 10);
%                 drawnow
%     if mod(i, nplt*4)==0
%         UOrig = [UOrig u'];
%     end
% end
%     UData1 = afinAsimetricoGeneralSecuencial(4, func1, xmin, xmax, N, tmax, Dt, 1);
%     plot(max(abs(UOrig-UData1)));

%Error
func1 = @(x)(1/2*c_1*(sech(sqrt(c_1)*(mod(x-4, xmax-xmin)+xmin)/2)).^2);
errorAsimetrico(8, func1, xmin, xmax, N, tmax, Dt, 1);
