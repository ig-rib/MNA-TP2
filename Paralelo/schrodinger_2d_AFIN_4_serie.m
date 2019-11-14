
 % Resolución de 2D  Schr \"{o} dinger NL usando
 % un metodo de splitting - Afin de 4º orden -

 clear all; 
 format compact ; 
 format short ;
 
 tic %inicio de toma de tiempos 
 
 % la grilla
  
 Lx = 20; % periodo 2* pi*L
 Ly = 20; % periodo 2* pi*L
 Nx = 2*256; % numero de armonicos
 Ny = 2*256; % numero de armonicos
 Nt = 100; % pasos
 %dt2 = 0.5/Nt; % 1/2paso de tiempo
 dt = 0.25/Nt;%1/4 paso
 dt2= 0.5/Nt; %1/2`paso
 q=10; %1/varianza
 Es = 1.0;
 
 % inicializazion de las variables
 
 x = (2* pi/Nx)*(- Nx /2: Nx /2 -1)'*Lx;   % x 
 kx = 1i *[0: Nx /2 -2  -Nx /2+2: -1]'/ Lx; % vector de onda
 y = (2* pi/Ny)*(- Ny /2: Ny /2 -1)'*Ly; % y 
 ky = 1i *[0: Ny /2 -2  -Ny /2+2: -1]'/ Ly; % vector de onda
 [xx ,yy ]= meshgrid(x,y);
 [k2xm , k2ym ]= meshgrid(kx.^2, ky.^2);
 
 % condicion inicial
 
 u = exp (-q*( xx.^2+ yy.^2));%gaussiana de inicio co varianza q
 figure (1); clf; mesh(xx ,yy ,u); drawnow ;%ploteo de la condicion de inicio
 t=0; tdata(1)=t;

 % masa
 ma0 = sum(abs (u).^2) ;
 %ma0 = ma(1 ,1);
 
 % coeficientes afines
 
 gamma1=2/3;
 gamma2=2/3;
 gamma3=-1/6;
 gamma4=-1/6;
 
 % calculo con LT y ploteo
 
 for n =2: Nt +1 
     
% 4LT LNLN
  v1 = fft2(u);%transformada de fourier
  vna1=exp (0.5 *1i*dt *( k2xm + k2ym )).*v1;%transformada a la primera aplicacion lineal
  una1= ifft2(vna1);% antitransformada de la 1º aplicación lineal
  pota1=Es *(( abs(una1)).^2);%1º composición con lo no lineal (el modulo)
  unb1=exp (-1i*dt*pota1).*una1;%1º composición con lo no lineal abs(u)^2.u
  vnb1= fft2 (unb1);%transformada a la 2º aplicación
  v11=exp(0.5 *1i*dt *( k2xm + k2ym )).*vnb1;%transformada a la 2º aplicación
  unc1= ifft2(v11);% antitransformada de la 2º aplicación lineal
  potc1=Es *(( abs(unc1)).^2);%2º composición con lo no lineal (el modulo)
  u1=exp (-1i*dt*potc1).*unc1;%2º composición con lo no lineal abs(u)^2.u

% 4LT NLNL
  pot2=Es *(( abs(u)).^2);%composición con lo no lineal (el modulo)
  und2=exp (-1i*dt*pot2).*u;%compocición con lo no lineal abs(u)^2.u
  vnd2= fft2 (und2);%transformada a la 1º aplicación
  vna2=exp (0.5 *1i*dt *( k2xm + k2ym )).*vnd2;%transformada a la primera aplicacion lineal
  una2= ifft2(vna2);% antitransformada de la 1º aplicación lineal
  pota2=Es *(( abs(una2)).^2);%1º composición con lo no lineal (el modulo)
  unb2=exp (-1i*dt*pota2).*una2;%1º composición con lo no lineal abs(u)^2.u
  vnb2= fft2 (unb2);%transformada a la 2º aplicación
  v2=exp(0.5 *1i*dt *( k2xm + k2ym )).*vnb2;%transformada a la 2º aplicación
  unc2= ifft2(v2);% antitransformada de la 2º aplicación lineal
  potc2=Es *(( abs(unc2)).^2);%2º composición con lo no lineal (el modulo)
  u2=exp (-1i*dt*potc2).*unc2;%2º composición con lo no lineal abs(u)^2.u
 
% 2 Lie Trotter LN
  
  v3 = fft2(u);%transformada de fourier
  vna3=exp (0.5 *1i*dt2 *( k2xm + k2ym )).*v3;%transformada a la primera aplicacion lineal
  una3= ifft2(vna3);% antitransformada de la 1º aplicación lineal
  pot3=Es *(( abs(una3)).^2);%composición con lo no lineal (el modulo)
  u3=exp (-1i*dt*pot3).*una3;%compocición con lo no lineal abs(u)^2.u

% 2  Lie Trotter NL
  
  pot4=Es *(( abs(u)).^2);%composición con lo no lineal (el modulo)
  unb4=exp (-1i*dt2*pot4).*u;%compocición con lo no lineal abs(u)^2.u
  vnb4= fft2 (unb4);%transformada a la 2º aplicación
  v4=exp(0.5 *1i*dt2 *( k2xm + k2ym )).*vnb4;%transformada a la 2º aplicación
  u4= ifft2(v4);% antitransformada de la 2º aplicación lineal

% combinación lineal de las aplicaciones parciales

  u = gamma1*u1 + gamma2*u2 + gamma3*u3 + gamma4*u4;
  
  t=(n -1)*dt;% tiempos
  tdata(n)=t;%creación de una sucesión de slides
  if ( mod(n ,10) ==0)%cada 10
  figure(2); clf ; mesh(xx ,yy ,abs(u).^2); title( num2str(t));
  drawnow ;
  
  ma =sum(abs(u).^2);
  %ma = ma(1,1);
  test = log10(abs(1- ma/ma0 ))%conservación por linea de comandos
  
 end
 end
 figure(4); clf; mesh(xx ,yy ,abs(u).^2) ;%figura final
 toc %medición de tiempo de corrida