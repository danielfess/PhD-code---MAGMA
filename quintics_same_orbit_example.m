//Equations for quintics landing in same orbit

//Structure of code:
//We work in the sextic resolvent ring S_f.
//Say x = \sum x_i \beta_i^*, similarly for y in terms of yi.
//Then we look for c,d,e in \mathbb{Z}^5 such that the matrix (x,c,d,e,y) is in GL5(\mathbb{Z})
//and c,d,e are orthonormal wrt 8 \Delta(f) (x^2,2xy,y^2).

//These are necessary conditions for being of the form \Phi(g)

//Note: Trivial case is x = e_1, y = e_5, c = e_2, d = e_3, e = e_4 (e_i are unit vectors).
//Less non-trivial but still general case comes from GL2 action.

Q := RationalField();
//P<f0,f1,f2,f3,f4,f5,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5> := PolynomialRing(Q,16);
P<x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,t1,t2,t3,t4> := PolynomialRing(Q,14);

f0 := 1;
f1 := 0;
f2 := 0;
f3 := 0;
f4 := 0;
f5 := 1;

Phi_f := [[0,t3,-t2,t1,0],[-t3,0,-f0*t1-f1*t2,-f2*t2-f3*t3,-t4],[t2,f0*t1+f1*t2,0,-f4*t3-f5*t4,t3],[-t1,f2*t2+f3*t3,f4*t3+f5*t4,0,-t2],[0,t4,-t3,t2,0]];

//Matrix functions

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

matrix_scalar := function(A,k)
	product := [];
	for i := 1 to #A do
		Append(~product,[]);
		for j := 1 to #A[1] do
			Append(~product[i],A[i,j]*k);
		end for;
	end for;
	return product;
end function;

transpose := function(A)
	A_t := [];
	for j := 1 to #(A[1]) do
		Append(~A_t,[]);
		for i := 1 to #A do
			A_t[j,i] := A[i,j];
		end for;
	end for;
	return A_t;
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

//Invariants and mult coeffs.  Need to change basis, new b1 = - old b5, new b5 = old b1.
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

//Haven't yet changed basis, here we go:
mult_coeffs2 := mult_coeffs;
swap := [1,6,3,4,5,2];
for i := 1 to 6 do
	for j := 1 to 6 do
		for k := 1 to 6 do
			I := swap[i];
			J := swap[j];
			K := swap[k];
			sign := 1;
			if i eq 2 then
				sign := - sign;
			end if;
			if j eq 2 then
				sign := - sign;
			end if;
			if k eq 2 then
				sign := - sign;
			end if;
			mult_coeffs2[i,j,k] := mult_coeffs[I,J,K]*sign;
		end for;
	end for;
end for;

mult_coeffs2[2,3,4] + mult_coeffs[6,3,4];
mult_coeffs[3,2,4] + mult_coeffs[6,3,4];
mult_coeffs2[4,3,6] - mult_coeffs[4,3,2];
mult_coeffs2[3,4,5] - mult_coeffs[3,4,5];
mult_coeffs2[2,6,4] + mult_coeffs[2,6,4];
mult_coeffs2[2,6,6] + mult_coeffs[2,6,2];
mult_coeffs2[2,6,2] - mult_coeffs[2,6,6];

//Trace of beta_i beta_j
trace := function(i,j)
	tr := 0;
	for k := 1 to 6 do
		for l := 1 to 6 do
			tr +:= mult_coeffs2[i,j,k]*mult_coeffs2[k,l,l];
		end for;
	end for;
	return tr;
end function;

//Matrix of trace pairing
T := [];
for i := 1 to 6 do
	Append(~T,[]);
	for j := 1 to 6 do
		Append(~T[i],trace(i,j));
	end for;
end for;
print T;

determinant(T);
// Equals (16*5^5)^3, as should be.  Good.

//Squaring map from dual lattice to ring, for original basis of form b1* = y, b5* = x, (b2,b3,b4) = 8 Delta (x^2,-2xy,y^2)
F1 := [[0,0,0,1,0],[0,-f0,0,-f2/2,0],[0,0,0,0,-1/2],[1,-f2/2,0,-f4,0],[0,0,-1/2,0,0]];
F2 := [[0,-f0,0,-f2/2,0],[-f0,3*f0*f2,f0*f3,f1*f3,-f1],[0,f0*f3,f0*f4,(f1*f4-f0*f5)/2,0],[-f2/2,f1*f3,(f1*f4-f0*f5)/2,f2*f4,-f3/2],[0,-f1,0,-f3/2,1]];
F3 := [[0,0,0,0,-1/2],[0,f0*f3,f0*f4,(f1*f4-f0*f5)/2,0],[0,f0*f4,3*f0*f5,f1*f5,0],[0,(f1*f4-f0*f5)/2,f1*f5,f2*f5,0],[-1/2,0,0,0,0]];
F4 := [[1,-f2/2,0,-f4,0],[-f2/2,f1*f3,(f1*f4-f0*f5)/2,f2*f4,-f3/2],[0,(f1*f4-f0*f5)/2,f1*f5,f2*f5,0],[-f4,f2*f4,f2*f5,3*f3*f5,-f5],[0,-f3/2,0,-f5,0]];
F5 := [[0,0,-1/2,0,0],[0,-f1,0,-f3/2,1],[-1/2,0,0,0,0],[0,-f3/2,0,-f5,0],[0,1,0,0,0]];

F := [F1,F2,F3,F4,F5];

//New squaring maps, after changing basis by b1 = -b5, b5 = b1, and correspondingly b1* = -b5*, b5* = b1*:
C := [[0,0,0,0,1],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[-1,0,0,0,0]];
G1 := matrix_mult(transpose(C),matrix_mult(matrix_scalar(F5,-1),C));
G2 := matrix_mult(transpose(C),matrix_mult(F2,C));
G3 := matrix_mult(transpose(C),matrix_mult(F3,C));
G4 := matrix_mult(transpose(C),matrix_mult(F4,C));
G5 := matrix_mult(transpose(C),matrix_mult(F1,C));

G := [G1,G2,G3,G4,G5];

dual_mult := function(a,b)
	product := [[P | 0,0,0,0,0]];
	for k := 1 to 5 do
		product[1,k] := matrix_mult(matrix_mult(a,G[k]),transpose(b))[1,1];
	end for;
	return product;
end function;

x := [[x1,x2,x3,x4,x5]];
y := [[y1,y2,y3,y4,y5]];
//x := [[1,0,x3,x4,x5]];
//y := [[0,1,y3,y4,y5]];

print transpose(x);

print dual_mult(x,x);

print matrix_mult(dual_mult(x,x),transpose(x));

segre := function(a)
	return matrix_mult(dual_mult(a,a),transpose(a))[1,1];
end function;

basis_vector := function(i)
	v := [0,0,0,0,0];
	for k := 1 to 5 do
		if k eq i then
			v[k] := 1;
		end if;
	end for;
	return v;
end function;

quad_subdet := function(i,j)
	B := [basis_vector(i),dual_mult(x,x)[1],dual_mult(x,y)[1],dual_mult(y,y)[1],basis_vector(j)];
	det := 2*determinant(B);
	return det;
end function;

print dual_mult(x,y);

pluecker := [ P | ];
for i := 1 to 5 do
	for j := i+1 to 5 do
		Append(~pluecker, x[1,i]*y[1,j] - x[1,j]*y[1,i]);
	end for;
end for;

J := ideal< P | pluecker >;

//Copying D_6 linear equations from phi_line_0.m (with pij replaced by qij):

eqns := [ P | ];

q12 := x1*y2 - x2*y1;
q13 := x1*y3 - x3*y1;
q14 := x1*y4 - x4*y1;
q15 := x1*y5 - x5*y1;
q23 := x2*y3 - x3*y2;
q24 := x2*y4 - x4*y2;
q25 := x2*y5 - x5*y2;
q34 := x3*y4 - x4*y3;
q35 := x3*y5 - x5*y3;
q45 := x4*y5 - x5*y4;

Append(~eqns,q14 - f0*q23);
Append(~eqns,q45 + q13 + f1*q23 + f2*q24);
Append(~eqns,q12 + q35 - f3*q24 - f4*q34);
Append(~eqns,q25 + f5*q34);

J := Ideal(eqns);

//Append(~eqns, q15^2 - 11*q15*q24 - q24^2 - 1);

I := Ideal(eqns);
IsPrime(I);

quad_subdet_polys := [];
for i := 1 to 5 do
	for j := i+1 to 5 do
		print i,j;
		g := quad_subdet(i,j);
		pij := x[1,i]*y[1,j] - x[1,j]*y[1,i];
		f := g - pij;
		Append(~eqns,f);
		Append(~quad_subdet_polys,f);
		if not IsIrreducible(f) then
			print "reducible";
		end if;
		if f in I then
			print "in I";
		end if;
//		Append(~eqns, x[1,i]*y[1,j] - x[1,j]*y[1,i] + quad_subdet(i,j));
		h := quad_subdet(i,j) - pij*(q15^2 - 11*q15*q24 - q24^2);
		if h in J then
			print "quad(i,j) = pij*unit eqn mod pluecker";
		end if;
	end for;
end for;

K := Ideal(eqns);
IsPrime(K);

//Note: For f = x^5 + y^5, the polynomials p_ij - quad_subdet(i,j) lie in the ideal generated by
//the four linear equations and the unit equation!
//Some subtlety in the signs of +-1 in unit equation and +-1 in front of quad_subdet(i,j).
//They need to be the same.  This must reflect the sign of det(M).

//Idea: Compute the quad_subdet's on the open set p12 \neq 0, map to the polynomial ring
//generated by the other pij, and rehomogenise using p12 to obtain a candidate for quad_subdet.

S<g0,g1,g2,g3,g4,g5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45> := PolynomialRing(Q,16);
h := hom< P -> S | g0,g1,g2,g3,g4,g5,1,0,-p23,-p24,-p25,0,1,p13,p14,p15,0>;
g_hom := hom< S -> P | f0,f1,f2,f3,f4,f5,q12,q13,q14,q15,q23,q24,q25,q34,q35,q45>;

//The form G = 4*<x^2,y^2> - <2xy,2xy> is a polynomial in pij, as follows:
G := -4*q12^2 - 4*q12*q13 - q13^2 - 4*q14^2 - 8*q14*q15 - q15^2 - 8*q15*q25 + 4*q14*q34 + 8*q12*q24 + 12*q13*q35 + 8*q12*q24 - q23^2 + 4*q23*q25 - 8*q12*q24 - q24^2 - 4*q25^2 + 12*q23*q34 + 8*q24*q45 - q34^2 - q35^2 - 4*q35*q45 - 4*q45^2 + 8*q13*q24 + 4*q13*q25 - 2*q14*q25 - 4*q12*q34 - 8*q12*q35 + 4*q14*q35 - 2*q12*q45 - 8*q13*q45 + 8*q24*q35 - 4*q23*q45;

R<p15,p24,p13,p23,p34,p35> := PolynomialRing(Q,6);
p12 := - p35;
p14 := p23;
p25 := - p34;
p45 := - p13;

GG := -4*p12^2 - 4*p12*p13 - p13^2 - 4*p14^2 - 8*p14*p15 - p15^2 - 8*p15*p25 + 4*p14*p34 + 8*p12*p24 + 12*p13*p35 + 8*p12*p24 - p23^2 + 4*p23*p25 - 8*p12*p24 - p24^2 - 4*p25^2 + 12*p23*p34 + 8*p24*p45 - p34^2 - p35^2 - 4*p35*p45 - 4*p45^2 + 8*p13*p24 + 4*p13*p25 - 2*p14*p25 - 4*p12*p34 - 8*p12*p35 + 4*p14*p35 - 2*p12*p45 - 8*p13*p45 + 8*p24*p35 - 4*p23*p45;

//Idea: If the trace pairing is non-degenerate on the span of beta_2^0, beta_3^0, beta_4^0, we can find beta_2^*, beta_3^*, beta_4^*
//in this span, modulo beta_1^*, beta_5^*.  Then, we find the coordinates of these with respect to a fixed basis of \tilde{S}
//and see if - together with beta_1^*, beta_5^* - they span \tilde{S}.

//Trace-zero translation, given coordinates in basis of S
tr_general := function(s)
	tr := 0;
	for i := 1 to 6 do
		tr +:= s[i]*trace(i,1);
	end for;
	return tr;
end function;

tr0 := function(s)
	return [s[1] - tr_general(s)/6,s[2],s[3],s[4],s[5],s[6]];
end function;

tr_general(tr0([1,2,3,4,5,6]));
xx := dual_mult(x,x);
xy := dual_mult(x,y);
yy := dual_mult(y,y);

xx0 := tr0([0,xx[1,1],xx[1,2],xx[1,3],xx[1,4],xx[1,5]]);
xy0 := tr0([0,2*xy[1,1],2*xy[1,2],2*xy[1,3],2*xy[1,4],2*xy[1,5]]);
yy0 := tr0([0,yy[1,1],yy[1,2],yy[1,3],yy[1,4],yy[1,5]]);

C_234 := [xx0,xy0,yy0];
T_234 := matrix_mult(C_234,matrix_mult(T,transpose(C_234)));
T_234_inv := cofactor_matrix(T_234);

//Want beta_2^0 etc in terms of dual basis of \tilde{S}
//This function takes coords in S and outputs dual coords
dual_coords := function(s);
	v := [];
	for j := 1 to 6 do
		tr := 0;
		for i := 1 to 6 do
			tr +:= s[i]*trace(i,j);
		end for;
		Append(~v,tr);
	end for;
	return v;
end function;

v2 := dual_coords(xx0);
v3 := dual_coords(xy0);
v4 := dual_coords(yy0);

V := [v2[2..6],v3[2..6],v4[2..6]];
dual_234 := matrix_mult(T_234_inv,V);
dual := [x[1],dual_234[1],dual_234[2],dual_234[3],y[1]];

//det := determinant(dual);
//Too big

quad_subdet2 := function(i,j)
	B := [basis_vector(i),dual_234[1],dual_234[2],dual_234[3],basis_vector(j)];
	det := determinant(B);
	return det;
end function;

M := [x[1],V[1],V[2],V[3],y[1]];
d := determinant(M);

//Evaluate \Phi(f)(x,beta_2^*) etc.

g11 := matrix_mult(x,matrix_mult(Phi_f,transpose([dual_234[1]])))[1,1];
g21 := matrix_mult(x,matrix_mult(Phi_f,transpose([dual_234[2]])))[1,1];
g31 := matrix_mult(x,matrix_mult(Phi_f,transpose([dual_234[3]])))[1,1];
g12 := matrix_mult(y,matrix_mult(Phi_f,transpose([dual_234[1]])))[1,1];
g22 := matrix_mult(y,matrix_mult(Phi_f,transpose([dual_234[2]])))[1,1];
g32 := matrix_mult(y,matrix_mult(Phi_f,transpose([dual_234[3]])))[1,1];

//Look for equations mod g(x,y) = 0.
//Need to separate g(x,y) into t1,..,t4 components.

eqns2 := [];
Append(~eqns2,q14 - f0*q23);
Append(~eqns2,q45 + q13 + f1*q23 + f2*q24);
Append(~eqns2,q12 + q35 - f3*q24 - f4*q34);
Append(~eqns2,q25 + f5*q34);
I2 := Ideal(eqns2);

//Note: segre(x) not in I2, but is in the non-trivial component of its
//primary decomposition, i.e. the component where some pluecker coord is non-zero.

K2 := PrimaryDecomposition(I2)[1];
Append(~eqns2,segre(x));
Append(~eqns2,matrix_mult(dual_mult(x,x),transpose(y))[1,1]);
Append(~eqns2,matrix_mult(dual_mult(y,y),transpose(x))[1,1]);
Append(~eqns2,segre(y));
J2 := Ideal(eqns2);
J2 eq K2;
//true: J2 = K2.

g11 + g22 in J2;
g21 + g32 in J2;

//Computing det of matrix x,beta_2_0,beta_3_0,beta_4_0,y in \tilde{S} coords
//....already computed - see matrix M above

//d;

//Computing det of this msame matrix but in S coords
//Want x,y in terms of basis of S
//First function computes trace(beta_i^*, beta_j^*) * Disc(S)^2.
//Second function takes coords in S^* and outputs S coords * Disc(S)^2 / constant
T_inv := cofactor_matrix(T);

trace_dual := function(i,j)
	b_i_dual := [T_inv[i]];
	b_j_dual := [T_inv[j]];
	return matrix_mult(b_i_dual,matrix_mult(T,transpose(b_j_dual)))[1,1];
end function;

print trace_dual(1,1);
print trace_dual(1,4);
print trace_dual(3,5);

ring_coords := function(s);
	v := [];
	for j := 1 to 6 do
		tr := 0;
		for i := 1 to 6 do
			tr +:= s[i]*trace_dual(i,j);
		end for;
		Append(~v,tr/(10^21*5^6));
	end for;
	return v;
end function;

//Add in zero entry for beta_0^* in x,y
x_plus0 := [0,x1,x2,x3,x4,x5];
y_plus0 := [0,y1,y2,y3,y4,y5];

N := [ring_coords(x_plus0)[2..6],xx0[2..6],xy0[2..6],yy0[2..6],ring_coords(y_plus0)[2..6]];
n := determinant(N);

tr_x2 := 0;
tr_xy := 0;
tr_y2 := 0;
for i := 1 to 5 do
	for j := 1 to 5 do
		tr_x2 +:= x[1,i]*x[1,j]*trace_dual(i+1,j+1)/(10^21*5^6);
		tr_xy +:= x[1,i]*y[1,j]*trace_dual(i+1,j+1)/(10^21*5^6);
		tr_y2 +:= y[1,i]*y[1,j]*trace_dual(i+1,j+1)/(10^21*5^6);
	end for;
end for;

I4 := tr_x2*tr_y2 - tr_xy^2;

//Taking into account our scalings here and there,
//the correct function to consider is n - I.

FF := n - I4;
FF in J2;
//false

//FF in Ideal(quad_subdet_polys);
//too long to compute

quad_subdet_polys_plus := quad_subdet_polys;
Append(~quad_subdet_polys_plus,q14 - f0*q23);
Append(~quad_subdet_polys_plus,q45 + q13 + f1*q23 + f2*q24);
Append(~quad_subdet_polys_plus,q12 + q35 - f3*q24 - f4*q34);
Append(~quad_subdet_polys_plus,q25 + f5*q34);
FF in Ideal(quad_subdet_polys_plus);
//true
//so (quad_subdet(i,j) - pij) equations are more fundamental than the equation FF

quad := q15^2 - 11*q15*q24 - q24^2 - 1;
Append(~eqns2,quad);
L2 := Ideal(eqns2);

FF in L2;
//true

FF2 := n - I4*(q15^2 - 11*q15*q24 - q24^2);
FF2 in J2;
//true

eqns3 := eqns2[1..8];
Append(~eqns3,FF);
O := Ideal(eqns3);
I4*quad in O;
//true - because we know n - I4*(quad + 1) is in J2

//Is the condition g(beta_1^*, beta_5^*) = 0 redundant?
//i.e. does this lie in the ideal generated by the quad_subdet equations?
quad_subdet_polys2 := Append(quad_subdet_polys,quad);
basis_ideal1 := Ideal(quad_subdet_polys2);
basis_ideal1 subset L2;

//Does basis_ideal1 = L2?
//L2 is ideal with matrix vanishing eqns, plus quad (and redundant segre eqns).
//basis_ideal also needs segre equations: it only says that beta_1^*, beta_5^*
//extend to a basis with the orthonormality property when we already know
//that beta_1^*, beta_5^* are orthogonal to beta_2, beta_3, beta_4.
//So, we need to add segre equations for this to be the correct ideal of basis equations.
basis_polys := quad_subdet_polys2;
for i := 1 to 4 do
	Append(~basis_polys,eqns2[4+i]);
end for;
basis_ideal2 := Ideal(basis_polys);
basis_ideal2 subset L2;

//Is quad captured by segre and other basis equations?
//Note: This is slow!!!!
segre_polys := eqns2[5..8];
segre_ideal := Ideal(segre_polys);

for i := 1 to 5 do
	for j := i+1 to 5 do
		print i,j;
		pij := x[1,i]*y[1,j] - x[1,j]*y[1,i];
		h := quad_subdet(i,j) - pij*(q15^2 - 11*q15*q24 - q24^2);
		if h in segre_ideal then
			print "quad(i,j) = pij*unit eqn mod segre_ideal";
		end if;
	end for;
end for;
//Result is that no h are in segre_ideal.
//Suggests that del pezzo equations are stronger


