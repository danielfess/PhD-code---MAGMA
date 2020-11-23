//Calculating the equation of the trace cubed = 0 variety inside \tilde{S}.

Q := RationalField();

P<f0,f1,f2,f3,f4,f5,y1,y2,y3,y4,y5> := PolynomialRing(Q,11);

invt := [ ];
for i := 1 to 6 do
	Append(~invt, [ ]);
	for j := 1 to 6 do
		Append(~invt[i], [ P | 0, 0, 0, 0, 0, 0 ]);
	end for;
end for;

invt[1, 1, 1] := -4*(100*f0*f3*f5 - 40*f0*f4^2 - 20*f1*f2*f5 + 4*f1*f3*f4 + 2*f2^2*f4 - f2*f3^2)/75;
invt[1, 1, 2] := (2500*f0^2*f5^2 - 600*f0*f1*f4*f5 + 500*f0*f2*f3*f5 - 320*f0*f2*f4^2 + 80*f0*f3^2*f4 - 320*f1^2*f3*f5 + 228*f1^2*f4^2 + 80*f1*f2^2*f5 - 92*f1*f2*f3*f4 + 16*f1*f3^3 + 16*f2^3*f4 - 3*f2^2*f3^2)/250;
invt[1, 1, 3] := (1000*f0*f1*f5^2 - 1200*f0*f2*f4*f5 + 950*f0*f3^2*f5 - 160*f0*f3*f4^2 + 200*f1^2*f4*f5 - 380*f1*f2*f3*f5 + 112*f1*f2*f4^2 - 2*f1*f3^2*f4 + 120*f2^3*f5 - 32*f2^2*f3*f4 + 7*f2*f3^3)/250;
invt[1, 1, 4] := -(4000*f0*f2*f5^2 - 2400*f0*f3*f4*f5 + 640*f0*f4^3 - 2000*f1^2*f5^2 + 800*f1*f2*f4*f5 + 200*f1*f3^2*f5 - 160*f1*f3*f4^2 - 160*f2^2*f3*f5 - 96*f2^2*f4^2 + 112*f2*f3^2*f4 - 21*f3^4)/250;
invt[1, 1, 5] := -4*(20*f1*f3*f5 - 8*f1*f4^2 - 10*f2^2*f5 + 4*f2*f3*f4 - f3^3)/25;
invt[1, 2, 1] := -(20*f1*f5 - 8*f2*f4 + 3*f3^2)/15;
invt[1, 2, 2] := (325*f0*f3*f5 - 160*f0*f4^2 - 50*f1*f2*f5 + 25*f1*f3*f4 - 4*f2^2*f4 - f2*f3^2)/75;
invt[1, 2, 3] := -(100*f0*f4*f5 - 100*f1*f3*f5 + 20*f1*f4^2 + 40*f2^2*f5 - 6*f2*f3*f4 + f3^3)/50;
invt[1, 2, 4] := -(250*f0*f5^2 - 50*f1*f4*f5 + 5*f2*f3*f5 + 8*f2*f4^2 - 3*f3^2*f4)/25;
invt[1, 2, 5] := 0;
invt[1, 3, 1] := -(50*f0*f5 - 6*f1*f4 + f2*f3)/15;
invt[1, 3, 2] := (50*f0*f2*f5 - 40*f1^2*f5 + 10*f1*f2*f4 - 6*f1*f3^2 + f2^2*f3)/25;
invt[1, 3, 3] := 2*(275*f0*f3*f5 - 80*f0*f4^2 - 70*f1*f2*f5 - f1*f3*f4 + 16*f2^2*f4 - 5*f2*f3^2)/75;
invt[1, 3, 4] := (100*f0*f4*f5 + 140*f1*f3*f5 - 76*f1*f4^2 - 80*f2^2*f5 + 42*f2*f3*f4 - 11*f3^3)/25;
invt[1, 3, 5] := -(20*f1*f5 - 8*f2*f4 + 3*f3^2)/5;
invt[1, 4, 1] := 2*(20*f0*f4 - 8*f1*f3 + 3*f2^2)/15;
invt[1, 4, 2] := (50*f0*f1*f5 - 120*f0*f2*f4 + 55*f0*f3^2 + 38*f1^2*f4 - 17*f1*f2*f3 + 4*f2^3)/25;
invt[1, 4, 3] := -(450*f0*f2*f5 - 80*f0*f3*f4 - 200*f1^2*f5 + 26*f1*f2*f4 + 2*f1*f3^2 - 3*f2^2*f3)/50;
invt[1, 4, 4] := -(875*f0*f3*f5 - 320*f0*f4^2 - 190*f1*f2*f5 + 23*f1*f3*f4 + 28*f2^2*f4 - 11*f2*f3^2)/75;
invt[1, 4, 5] := (50*f0*f5 - 6*f1*f4 + f2*f3)/5;
invt[1, 5, 1] := 4*(100*f0*f2*f5 - 20*f0*f3*f4 - 40*f1^2*f5 + 4*f1*f2*f4 + 2*f1*f3^2 - f2^2*f3)/75;
invt[1, 5, 2] := -(1000*f0^2*f4*f5 - 1200*f0*f1*f3*f5 + 200*f0*f1*f4^2 + 950*f0*f2^2*f5 - 380*f0*f2*f3*f4 + 120*f0*f3^3 - 160*f1^2*f2*f5 + 112*f1^2*f3*f4 - 2*f1*f2^2*f4 - 32*f1*f2*f3^2 + 7*f2^3*f3)/250;
invt[1, 5, 3] := -(1250*f0^2*f5^2 - 100*f0*f1*f4*f5 + 450*f0*f2*f3*f5 - 400*f0*f2*f4^2 + 110*f0*f3^2*f4 - 400*f1^2*f3*f5 + 210*f1^2*f4^2 + 110*f1*f2^2*f5 - 54*f1*f2*f3*f4 + 4*f1*f3^3 + 4*f2^3*f4 + f2^2*f3^2)/250;
invt[1, 5, 4] := -(1000*f0*f1*f5^2 - 1200*f0*f2*f4*f5 + 950*f0*f3^2*f5 - 160*f0*f3*f4^2 + 200*f1^2*f4*f5 - 380*f1*f2*f3*f5 + 112*f1*f2*f4^2 - 2*f1*f3^2*f4 + 120*f2^3*f5 - 32*f2^2*f3*f4 + 7*f2*f3^3)/250;
invt[1, 5, 5] := 4*(100*f0*f3*f5 - 40*f0*f4^2 - 20*f1*f2*f5 + 4*f1*f3*f4 + 2*f2^2*f4 - f2*f3^2)/75;
invt[2, 1, 1] := -(20*f1*f5 - 8*f2*f4 + 3*f3^2)/15;
invt[2, 1, 2] := (325*f0*f3*f5 - 160*f0*f4^2 - 50*f1*f2*f5 + 25*f1*f3*f4 - 4*f2^2*f4 - f2*f3^2)/75;
invt[2, 1, 3] := -(100*f0*f4*f5 - 100*f1*f3*f5 + 20*f1*f4^2 + 40*f2^2*f5 - 6*f2*f3*f4 + f3^3)/50;
invt[2, 1, 4] := -(250*f0*f5^2 - 50*f1*f4*f5 + 5*f2*f3*f5 + 8*f2*f4^2 - 3*f3^2*f4)/25;
invt[2, 1, 5] := 0;
invt[2, 2, 1] := -2*f4;
invt[2, 2, 2] := 2*(10*f1*f5 - f2*f4)/15;
invt[2, 2, 3] := (5*f2*f5 - f3*f4)/5;
invt[2, 2, 4] := 2*(5*f3*f5 - 2*f4^2)/5;
invt[2, 2, 5] := -10*f5;
invt[2, 3, 1] := 2*f3;
invt[2, 3, 2] := -2*(25*f0*f5 + 3*f1*f4 - f2*f3)/15;
invt[2, 3, 3] := -(20*f1*f5 + 4*f2*f4 - 3*f3^2)/15;
invt[2, 3, 4] := -2*(5*f2*f5 - f3*f4)/5;
invt[2, 3, 5] := 4*f4;
invt[2, 4, 1] := -f2;
invt[2, 4, 2] := -2*(10*f0*f4 - 7*f1*f3 + 3*f2^2)/15;
invt[2, 4, 3] := -(25*f0*f5 - 5*f1*f4 + f2*f3)/5;
invt[2, 4, 4] := -2*(10*f1*f5 - 7*f2*f4 + 3*f3^2)/15;
invt[2, 4, 5] := -f3;
invt[2, 5, 1] := (50*f0*f5 - 6*f1*f4 + f2*f3)/5;
invt[2, 5, 2] := -(875*f0*f2*f5 - 190*f0*f3*f4 - 320*f1^2*f5 + 23*f1*f2*f4 + 28*f1*f3^2 - 11*f2^2*f3)/75;
invt[2, 5, 3] := -(450*f0*f3*f5 - 200*f0*f4^2 - 80*f1*f2*f5 + 26*f1*f3*f4 + 2*f2^2*f4 - 3*f2*f3^2)/50;
invt[2, 5, 4] := (50*f0*f4*f5 - 120*f1*f3*f5 + 38*f1*f4^2 + 55*f2^2*f5 - 17*f2*f3*f4 + 4*f3^3)/25;
invt[2, 5, 5] := 2*(20*f1*f5 - 8*f2*f4 + 3*f3^2)/15;
invt[3, 1, 1] := -(50*f0*f5 - 6*f1*f4 + f2*f3)/15;
invt[3, 1, 2] := (50*f0*f2*f5 - 40*f1^2*f5 + 10*f1*f2*f4 - 6*f1*f3^2 + f2^2*f3)/25;
invt[3, 1, 3] := 2*(275*f0*f3*f5 - 80*f0*f4^2 - 70*f1*f2*f5 - f1*f3*f4 + 16*f2^2*f4 - 5*f2*f3^2)/75;
invt[3, 1, 4] := (100*f0*f4*f5 + 140*f1*f3*f5 - 76*f1*f4^2 - 80*f2^2*f5 + 42*f2*f3*f4 - 11*f3^3)/25;
invt[3, 1, 5] := -(20*f1*f5 - 8*f2*f4 + 3*f3^2)/5;
invt[3, 2, 1] := 2*f3;
invt[3, 2, 2] := -2*(25*f0*f5 + 3*f1*f4 - f2*f3)/15;
invt[3, 2, 3] := -(20*f1*f5 + 4*f2*f4 - 3*f3^2)/15;
invt[3, 2, 4] := -2*(5*f2*f5 - f3*f4)/5;
invt[3, 2, 5] := 4*f4;
invt[3, 3, 1] := -4*f2;
invt[3, 3, 2] := 2*(20*f0*f4 - 4*f1*f3 + f2^2)/5;
invt[3, 3, 3] := 2*(100*f0*f5 - f2*f3)/15;
invt[3, 3, 4] := 2*(20*f1*f5 - 4*f2*f4 + f3^2)/5;
invt[3, 3, 5] := -4*f3;
invt[3, 4, 1] := 4*f1;
invt[3, 4, 2] := -2*(5*f0*f3 - f1*f2)/5;
invt[3, 4, 3] := -(20*f0*f4 + 4*f1*f3 - 3*f2^2)/15;
invt[3, 4, 4] := -2*(25*f0*f5 + 3*f1*f4 - f2*f3)/15;
invt[3, 4, 5] := 2*f2;
invt[3, 5, 1] := -(20*f0*f4 - 8*f1*f3 + 3*f2^2)/5;
invt[3, 5, 2] := (100*f0*f1*f5 + 140*f0*f2*f4 - 80*f0*f3^2 - 76*f1^2*f4 + 42*f1*f2*f3 - 11*f2^3)/25;
invt[3, 5, 3] := 2*(275*f0*f2*f5 - 70*f0*f3*f4 - 80*f1^2*f5 - f1*f2*f4 + 16*f1*f3^2 - 5*f2^2*f3)/75;
invt[3, 5, 4] := (50*f0*f3*f5 - 40*f0*f4^2 + 10*f1*f3*f4 - 6*f2^2*f4 + f2*f3^2)/25;
invt[3, 5, 5] := -(50*f0*f5 - 6*f1*f4 + f2*f3)/15;
invt[4, 1, 1] := 2*(20*f0*f4 - 8*f1*f3 + 3*f2^2)/15;
invt[4, 1, 2] := (50*f0*f1*f5 - 120*f0*f2*f4 + 55*f0*f3^2 + 38*f1^2*f4 - 17*f1*f2*f3 + 4*f2^3)/25;
invt[4, 1, 3] := -(450*f0*f2*f5 - 80*f0*f3*f4 - 200*f1^2*f5 + 26*f1*f2*f4 + 2*f1*f3^2 - 3*f2^2*f3)/50;
invt[4, 1, 4] := -(875*f0*f3*f5 - 320*f0*f4^2 - 190*f1*f2*f5 + 23*f1*f3*f4 + 28*f2^2*f4 - 11*f2*f3^2)/75;
invt[4, 1, 5] := (50*f0*f5 - 6*f1*f4 + f2*f3)/5;
invt[4, 2, 1] := -f2;
invt[4, 2, 2] := -2*(10*f0*f4 - 7*f1*f3 + 3*f2^2)/15;
invt[4, 2, 3] := -(25*f0*f5 - 5*f1*f4 + f2*f3)/5;
invt[4, 2, 4] := -2*(10*f1*f5 - 7*f2*f4 + 3*f3^2)/15;
invt[4, 2, 5] := -f3;
invt[4, 3, 1] := 4*f1;
invt[4, 3, 2] := -2*(5*f0*f3 - f1*f2)/5;
invt[4, 3, 3] := -(20*f0*f4 + 4*f1*f3 - 3*f2^2)/15;
invt[4, 3, 4] := -2*(25*f0*f5 + 3*f1*f4 - f2*f3)/15;
invt[4, 3, 5] := 2*f2;
invt[4, 4, 1] := -10*f0;
invt[4, 4, 2] := 2*(5*f0*f2 - 2*f1^2)/5;
invt[4, 4, 3] := (5*f0*f3 - f1*f2)/5;
invt[4, 4, 4] := 2*(10*f0*f4 - f1*f3)/15;
invt[4, 4, 5] := -2*f1;
invt[4, 5, 1] := 0;
invt[4, 5, 2] := -(250*f0^2*f5 - 50*f0*f1*f4 + 5*f0*f2*f3 + 8*f1^2*f3 - 3*f1*f2^2)/25;
invt[4, 5, 3] := -(100*f0*f1*f5 - 100*f0*f2*f4 + 40*f0*f3^2 + 20*f1^2*f4 - 6*f1*f2*f3 + f2^3)/50;
invt[4, 5, 4] := (325*f0*f2*f5 - 50*f0*f3*f4 - 160*f1^2*f5 + 25*f1*f2*f4 - 4*f1*f3^2 - f2^2*f3)/75;
invt[4, 5, 5] := -(20*f0*f4 - 8*f1*f3 + 3*f2^2)/15;
invt[5, 1, 1] := 4*(100*f0*f2*f5 - 20*f0*f3*f4 - 40*f1^2*f5 + 4*f1*f2*f4 + 2*f1*f3^2 - f2^2*f3)/75;
invt[5, 1, 2] := -(1000*f0^2*f4*f5 - 1200*f0*f1*f3*f5 + 200*f0*f1*f4^2 + 950*f0*f2^2*f5 - 380*f0*f2*f3*f4 + 120*f0*f3^3 - 160*f1^2*f2*f5 + 112*f1^2*f3*f4 - 2*f1*f2^2*f4 - 32*f1*f2*f3^2 + 7*f2^3*f3)/250;
invt[5, 1, 3] := -(1250*f0^2*f5^2 - 100*f0*f1*f4*f5 + 450*f0*f2*f3*f5 - 400*f0*f2*f4^2 + 110*f0*f3^2*f4 - 400*f1^2*f3*f5 + 210*f1^2*f4^2 + 110*f1*f2^2*f5 - 54*f1*f2*f3*f4 + 4*f1*f3^3 + 4*f2^3*f4 + f2^2*f3^2)/250;
invt[5, 1, 4] := -(1000*f0*f1*f5^2 - 1200*f0*f2*f4*f5 + 950*f0*f3^2*f5 - 160*f0*f3*f4^2 + 200*f1^2*f4*f5 - 380*f1*f2*f3*f5 + 112*f1*f2*f4^2 - 2*f1*f3^2*f4 + 120*f2^3*f5 - 32*f2^2*f3*f4 + 7*f2*f3^3)/250;
invt[5, 1, 5] := 4*(100*f0*f3*f5 - 40*f0*f4^2 - 20*f1*f2*f5 + 4*f1*f3*f4 + 2*f2^2*f4 - f2*f3^2)/75;
invt[5, 2, 1] := (50*f0*f5 - 6*f1*f4 + f2*f3)/5;
invt[5, 2, 2] := -(875*f0*f2*f5 - 190*f0*f3*f4 - 320*f1^2*f5 + 23*f1*f2*f4 + 28*f1*f3^2 - 11*f2^2*f3)/75;
invt[5, 2, 3] := -(450*f0*f3*f5 - 200*f0*f4^2 - 80*f1*f2*f5 + 26*f1*f3*f4 + 2*f2^2*f4 - 3*f2*f3^2)/50;
invt[5, 2, 4] := (50*f0*f4*f5 - 120*f1*f3*f5 + 38*f1*f4^2 + 55*f2^2*f5 - 17*f2*f3*f4 + 4*f3^3)/25;
invt[5, 2, 5] := 2*(20*f1*f5 - 8*f2*f4 + 3*f3^2)/15;
invt[5, 3, 1] := -(20*f0*f4 - 8*f1*f3 + 3*f2^2)/5;
invt[5, 3, 2] := (100*f0*f1*f5 + 140*f0*f2*f4 - 80*f0*f3^2 - 76*f1^2*f4 + 42*f1*f2*f3 - 11*f2^3)/25;
invt[5, 3, 3] := 2*(275*f0*f2*f5 - 70*f0*f3*f4 - 80*f1^2*f5 - f1*f2*f4 + 16*f1*f3^2 - 5*f2^2*f3)/75;
invt[5, 3, 4] := (50*f0*f3*f5 - 40*f0*f4^2 + 10*f1*f3*f4 - 6*f2^2*f4 + f2*f3^2)/25;
invt[5, 3, 5] := -(50*f0*f5 - 6*f1*f4 + f2*f3)/15;
invt[5, 4, 1] := 0;
invt[5, 4, 2] := -(250*f0^2*f5 - 50*f0*f1*f4 + 5*f0*f2*f3 + 8*f1^2*f3 - 3*f1*f2^2)/25;
invt[5, 4, 3] := -(100*f0*f1*f5 - 100*f0*f2*f4 + 40*f0*f3^2 + 20*f1^2*f4 - 6*f1*f2*f3 + f2^3)/50;
invt[5, 4, 4] := (325*f0*f2*f5 - 50*f0*f3*f4 - 160*f1^2*f5 + 25*f1*f2*f4 - 4*f1*f3^2 - f2^2*f3)/75;
invt[5, 4, 5] := -(20*f0*f4 - 8*f1*f3 + 3*f2^2)/15;
invt[5, 5, 1] := -4*(20*f0*f2*f4 - 10*f0*f3^2 - 8*f1^2*f4 + 4*f1*f2*f3 - f2^3)/25;
invt[5, 5, 2] := -(4000*f0^2*f3*f5 - 2000*f0^2*f4^2 - 2400*f0*f1*f2*f5 + 800*f0*f1*f3*f4 + 200*f0*f2^2*f4 - 160*f0*f2*f3^2 + 640*f1^3*f5 - 160*f1^2*f2*f4 - 96*f1^2*f3^2 + 112*f1*f2^2*f3 - 21*f2^4)/250;
invt[5, 5, 3] := (1000*f0^2*f4*f5 - 1200*f0*f1*f3*f5 + 200*f0*f1*f4^2 + 950*f0*f2^2*f5 - 380*f0*f2*f3*f4 + 120*f0*f3^3 - 160*f1^2*f2*f5 + 112*f1^2*f3*f4 - 2*f1*f2^2*f4 - 32*f1*f2*f3^2 + 7*f2^3*f3)/250;
invt[5, 5, 4] := (2500*f0^2*f5^2 - 600*f0*f1*f4*f5 + 500*f0*f2*f3*f5 - 320*f0*f2*f4^2 + 80*f0*f3^2*f4 - 320*f1^2*f3*f5 + 228*f1^2*f4^2 + 80*f1*f2^2*f5 - 92*f1*f2*f3*f4 + 16*f1*f3^3 + 16*f2^3*f4 - 3*f2^2*f3^2)/250;
invt[5, 5, 5] := -4*(100*f0*f2*f5 - 20*f0*f3*f4 - 40*f1^2*f5 + 4*f1*f2*f4 + 2*f1*f3^2 - f2^2*f3)/75;

key_invts := [invt[1,2,2],invt[2,4,4],invt[3,3,3],invt[4,2,2],invt[5,4,4]];

mult_coeffs := invt;

mult_coeffs[3,5,5] := invt[3,3,4]/4;
mult_coeffs[5,3,3] := invt[3,3,2]/4;
mult_coeffs[4,4,4] := 4*invt[2,4,3];
mult_coeffs[2,3,3] := 0;
mult_coeffs[6,5,5] := 0;

key_mult := [mult_coeffs[2,3,3], mult_coeffs[3,5,5], mult_coeffs[4,4,4], mult_coeffs[5,3,3], mult_coeffs[6,5,5]];


for i := 2 to 6 do
	for j := 2 to 6 do
		for k := 2 to 6 do
			if i eq j then
				if i eq k then
					if i ne 4 then
						mult_coeffs[i,j,k] := invt[i-1,j-1,k-1] - 2*key_invts[i-1] + 2*key_mult[i-1];
					end if;
				else
					mult_coeffs[i,j,k] := invt[i-1,j-1,k-1];
				end if;
			else
				if i eq k then
					if j ne 4 then
						mult_coeffs[i,j,i] := invt[i-1,j-1,i-1] - key_invts[j-1] + key_mult[j-1];
					else
						mult_coeffs[i,4,i] := mult_coeffs[4,4,4]/2 - invt[3,3,3]/2 + invt[i-1,3,i-1];
					end if;
				elif j eq k then
					if i ne 4 then
						mult_coeffs[i,j,j] := invt[i-1,j-1,j-1] - key_invts[i-1] + key_mult[i-1];
					else
						mult_coeffs[4,j,j] := mult_coeffs[4,4,4]/2 - invt[3,3,3]/2 + invt[3,j-1,j-1];
					end if;
				else
					mult_coeffs[i,j,k] := invt[i-1,j-1,k-1];
				end if;
			end if;
		end for;
	end for;
end for;

for i := 1 to 6 do
	for k := 1 to 6 do
		if i eq k then
			mult_coeffs[i,1,k] := 1;
			mult_coeffs[1,i,k] := 1;
		else
			mult_coeffs[i,1,k] := 0;
			mult_coeffs[1,i,k] := 0;
		end if;
	end for;
end for;

not_i := [3,4,5,6,2];

for i := 2 to 6 do
	for j := 2 to 6 do
		k := not_i[i-1];
		constant_term := 0;
		for r := 2 to 6 do
			constant_term +:= mult_coeffs[j,k,r]*mult_coeffs[r,i,k] - mult_coeffs[i,j,r]*mult_coeffs[r,k,k];
		end for;
		mult_coeffs[i,j,1] := constant_term;
	end for;
end for;

S := AssociativeAlgebra< P, 6 | mult_coeffs >;

//Trace of beta_i beta_j, indexed from 1 to 6
trace := function(i,j)
	tr := 0;
	for k := 1 to 6 do
		for l := 1 to 6 do
			tr +:= mult_coeffs[i,j,k]*mult_coeffs[k,l,l];
		end for;
	end for;
	return tr;
end function;

matrix_mult := function(A,B)
	product := [];
	for i := 1 to #A do
		Append(~product,[]);
		for j := 1 to #B[1] do
			entry := 0;
			for k := 1 to #A[1] do
				entry +:= A[i,k]*B[k,j];
			end for;
			Append(~product[i],entry);
		end for;
	end for;
	return product;
end function;

submatrix := function(A,i,j)
	Asub := A[1..i-1] cat A[i+1..#A];
	for k := 1 to #A-1 do
        Asub[k] := Asub[k][1..j-1] cat Asub[k][j+1..#A[1]];
    end for;
    return Asub;
end function;

determinant := function(A) 
    total := 0;

    if #A eq 1 then
        return A[1,1];
    end if;

    for col := 1 to #A do
        Asub := submatrix(A,1,col);
        subdet := $$(Asub);
        sign := (-1) ^ (col+1 mod 2);
        total +:= sign * A[1,col] * subdet;
    end for;
    return total;
end function;

cofactor_matrix := function(A)
	cofactor := A;
	if #A eq 1 then
		return A;
	end if;
	for i := 1 to #A do
		for j := 1 to #A do
			sub := submatrix(A,i,j);
			sign := (-1) ^ (i+j mod 2);
			cofactor[j,i] := sign*determinant(sub);
		end for;
	end for;
	return cofactor;
end function;

T := [];
for i := 1 to 6 do
	Append(~T,[]);
	for j := 1 to 6 do
		Append(~T[i],trace(i,j));
	end for;
end for;
print T;

T_inv := cofactor_matrix(T);
det := determinant(T);
print( T[1,4]^2 - 4*T[1,3]*T[1,5]);
G := 4*T[2,6]^2 - 4*T[2,2]*T[6,6];
print G;
Factorisation(G);

disc := Discriminant(f0*y1^5 + f1*y1^4 + f2*y1^3 + f3*y1^2 + f4*y1 + f5,y1);

//Divide entries of cofactor matrix by disc^2:
T_inv2 := T_inv;
for i := 1 to 6 do
	for j := 1 to 6 do
		T_inv2[i,j] := T_inv[i,j]/(disc^2);
	end for;
end for;

//Vector of generic element of \tilde{S}: (y_1 \beta_1^* + ... + y_5 \beta_5^*)*disc*16^3.
//Adjusting yi's so that we are working with \Phi(f), not \Phi'(f).
//This amounts to multiplying by a specific element of GL5.
y := [y1 - 2/5*f2*y2 - f3/10*y3 - 2/5*f4*y4,y2,y3,y4,y5 - 2/5*f1*y2 - f2/10*y3 - 2/5*f3*y4];
elt := [ P | 0,0,0,0,0,0];
elt2 := [ P | 0,0,0,0,0,0];
elt3 := [ P | 0,0,0,0,0,0];
for k := 1 to 6 do
	entry := 0;
	for l := 1 to 5 do
		entry +:= y[l]*T_inv2[k,l+1];
	end for;
	print entry;
	elt[k] := entry;
end for;

sextic_prod := function(x,y)
	prod := [P | 0,0,0,0,0,0];
	for k := 1 to 6 do
		prod_coeff := 0;
		for i := 1 to 6 do
			for j := 1 to 6 do
				prod_coeff +:= mult_coeffs[i,j,k]*x[i]*y[j];
			end for;
		end for;
	prod[k] := prod_coeff;
	end for;
	return prod;
end function;

elt2 := sextic_prod(elt,elt);
elt3 := sextic_prod(elt,elt2);

//Calculate trace cubed = trace(elt3), then divide by disc^2
trace_cubed := 0;
for k := 1 to 6 do
	trace_cubed +:= elt3[k]*trace(k,1);
end for;

print trace_cubed;

trace_cubed2 := Factorisation(trace_cubed)[1,1];
print trace_cubed2;

segre := f0*f2*y2^3 + f0*f3*y2^2*y3 + f0*f4*y2*y3^2 - f0*f5*y2*y3*y4 + f0*f5*y3^3 - f0*y1*y2^2 + f1*f3*y2^2*y4 + f1*f4*y2*y3*y4 + f1*f5*y3^2*y4 - f1*y2^2*y5 + f2*f4*y2*y4^2 + f2*f5*y3*y4^2 - f2*y1*y2*y4 + f3*f5*y4^3 - f3*y2*y4*y5 - f4*y1*y4^2 - f5*y4^2*y5 + y1^2*y4 - y1*y3*y5 + y2*y5^2;
I := Ideal([trace_cubed2,segre]);


