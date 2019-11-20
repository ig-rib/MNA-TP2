%%tests

c_1=13;
c_2 =3;
func1 = @(x)(1/2*c_1*(sech(sqrt(c_1)*(x+3)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x+1)/2)).^2);


%Secuencial vs Paralelo
 %afinAsimetricoGeneralSecuencial(4, func1, -5, 5, 128, 1, 0.0001, 1);
 %afinAsimetricoGeneral(4, func1, -5, 5, 128, 1, 0.0001);
%------------------------------------------------------------------------%
%Prueba de error entre distintos ordenes
% for i=4:7
% UData1 = afinAsimetricoGeneralSecuencial(i, func1, -5, 5, 128, 1, 0.0001, 1);
% UData2 =  afinAsimetricoGeneralSecuencial(i+1, func1, -5, 5, 128, 1, 0.0001, 1);
% %UData2 = UData2(1:2:2^(7+i), :);
% Error = max(abs(UData1-UData2));
% figure
% plot(Error),
% %%loglog(Error);
% end
%------------------------------------------------------------------------%
Dt = 0.0001;
tmax = 1;
xmin = -4;
xmax = 4;
N = 64;
x = linspace(xmin,xmax,N);
Dx = x(2)-x(1);
nmax = round(tmax/Dt);
nplt = floor((tmax/100)/Dt);

func1 = @(x)(1/2*c_1*(sech(sqrt(c_1)*(mod(x+3, xmax-xmin)+xmin)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(mod(x+1,xmax-xmin)+xmin)/2)).^2);

solnExacta = @(x, t)(1/2*c_1*(sech(sqrt(c_1)*(mod(x-2-c_1*(t), xmax-xmin)+xmin)/2)).^2);
%Error del afin asimetrico vs solucion analitica
 func1 = @(x)(1/2*c_1*(sech(sqrt(c_1)*(mod(x-2, xmax-xmin)+xmin)/2)).^2);
 %errorAsimetrico(2, func1, solnExacta, xmin, xmax, N, tmax, Dt, 1);
 %errorYoshida(func1, solnExacta, xmin, xmax, N, tmax, Dt, 1);
 errorStrang(func1, solnExacta, xmin, xmax, N, tmax, Dt, 1);
%------------------------------------------------------------------------%

%KdV2solitones_yoshida(func1, xmin, xmax, N, tmax, Dt, 1);

%------------------------------------------------------------------------%

%aagSPMD(2, func1, xmin, xmax, N, tmax, Dt, 1);