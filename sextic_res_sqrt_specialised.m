//Magma code to find a square root of beta_2 or beta_4 (up to scalars) for a specific f
//Use other code to check output.

//Strategy
//For each k, let pk = sum_i,j d_ij^k ai aj.  These are the coordinates of a square in our sextic resolvent.
//We want to see if beta_2 is a square (up to a scalar), so we look for solutions to all pk = 0 for k ne 2.
//The pk generate an ideal I.  We want to find solutions to I.
//Approach 1: Find Groebner basis of I......this is too slow.
//Approach 2: Primary decomposition of I.  This works for specialised f, but I haven't yet had success for general f.


//Note: beta_2 is the 3rd element of the basis, and because of how Magma indexing works, there will be indexing adjustments below.
//beta_2, beta_3, beta_4 normalised so that 4 beta_2 beta_4 = beta_3^2.

Q := RationalField();
P<a0,a1,a2,a3,a4,a5> := PolynomialRing(Q,6);

f0 := 3;
f1 := 10;
f2 := -20;
f3 := -70;
f4 := 60;
f5 := 11;

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

a := [a0, a1, a2, a3, a4, a5];
r := 2;
//r can be 2 or 4, depending on whether we want to compute a square root of beta_2 or beta_4.
square_polys := [P | ];
for k := 1 to 6 do
	square_coeff := 0;
	for i := 1 to 6 do
		for j := 1 to 6 do
			square_coeff +:= mult_coeffs[i,j,k]*a[i]*a[j];
		end for;
	end for;
	if k ne r+1 then
		Append(~square_polys,square_coeff);
	else
		coeff_r := square_coeff;
	end if;
end for;

print(P);
print(square_polys);
print(coeff_r);

I := Ideal(square_polys);
PrimaryDecomposition(I);
print(I);

