%%tests

c_1=13;
c_2 =3;
func1 = @(x)(1/2*c_1*(sech(sqrt(c_1)*(x+8)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x+1)/2)).^2);

%%afinAsimetricoGeneral(4, func1, -10, 10, 256);
afinSimetricoGeneral(4, func1, -10, 10, 256);