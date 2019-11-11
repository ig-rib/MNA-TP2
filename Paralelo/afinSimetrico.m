%%Afin simetrico de orden q

N = 256;
x = linspace(-10,10,N);
Dx = x(2) - x(1);
Dk = 2*pi/(N*Dx);

k = [0:Dk:(N/2-1)*Dk,0,-(N/2-1)*Dk:Dk:-Dk];
c_1=13;
c_2 =3;

u = 1/2*c_1*(sech(sqrt(c_1)*(x/16+8)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x/2+1)/2)).^2;

Dt = 0.4/N^2;
t=0;

UData = u';
TData = 0;

U = fft(u);

tmax = 1.5; nplt = floor((tmax/100)/Dt); nmax = round(tmax/Dt);

for n = 1:nmax-40000
    t = n*Dt;
    
    U = symmetricIntegrator(4, Dt, U, [-1/6 2/3]);
end