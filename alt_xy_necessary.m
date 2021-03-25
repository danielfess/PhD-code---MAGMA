//Equations for quintics landing in same orbit

//Structure of code:
//We work in the sextic resolvent ring S_f.
//Say x = \sum x_i \beta_i^*, similarly for y in terms of yi.
//Then we look for c,d,e in \mathbb{Z}^5 such that the matrix (x,c,d,e,y) is in GL5(\mathbb{Z})
//and c,d,e are orthonormal wrt 8 \Delta(f) (x^2,2xy,y^2).

//These are necessary conditions for being of the form \Phi(g)

//Note: Trivial case is x = e_1, y = e_5, c = e_2, d = e_3, e = e_4 (e_i are unit vectors).
//Less non-trivial but still general case comes from GL2 action.

//In this script:::::
//We investigate whether the D6 equations follow from the basis equations.

Q := RationalField();
//P<x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,t1,t2,t3,t4> := PolynomialRing(Q,14);
P<x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,f0,f1,f2,f3,f4,f5> := PolynomialRing(Q,16);
R<x1,x2,x3,x4,x5,y1,y2,y3,y4,y5> := PolynomialRing(Q,10);

//f0 := 1;
//f1 := 0;
//f2 := 0;
//f3 := 0;
//f4 := 0;
//f5 := 1;

//Phi_f := [[0,t3,-t2,t1,0],[-t3,0,-f0*t1-f1*t2,-f2*t2-f3*t3,-t4],[t2,f0*t1+f1*t2,0,-f4*t3-f5*t4,t3],[-t1,f2*t2+f3*t3,f4*t3+f5*t4,0,-t2],[0,t4,-t3,t2,0]];

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
//Change of basis comes later
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

//Trace of beta_i beta_j, in new basis.
trace := function(i,j)
	tr := 0;
	for k := 1 to 6 do
		for l := 1 to 6 do
			tr +:= mult_coeffs2[i,j,k]*mult_coeffs2[k,l,l];
		end for;
	end for;
	return tr;
end function;

//Matrix of trace pairing, in new basis.
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

quad_subdet := function(x,y,i,j)
	B := [basis_vector(i),dual_mult(x,x)[1],dual_mult(x,y)[1],dual_mult(y,y)[1],basis_vector(j)];
	det := 2*determinant(B);
	return det;
end function;

print dual_mult(x,y);

//Function evaluating D_6 equations and returning true if they all vanish:

D6_polys := [];
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
Append(~D6_polys,q14 - f0*q23);
Append(~D6_polys,q45 + q13 + f1*q23 + f2*q24);
Append(~D6_polys,q12 + q35 - f3*q24 - f4*q34);
Append(~D6_polys,q25 + f5*q34);

D6_hom := function(s1,s2,s3,s4,s5,t1,t2,t3,t4,t5)
	evals := [];
	h := hom< P -> Q | s1,s2,s3,s4,s5,t1,t2,t3,t4,t5 >;
	for p in D6_polys do
		Append(~evals,h(p));
	end for;
	if evals eq [0,0,0,0] then
		return true;
	else
		return false;
	end if;
end function;

D6_basic := function(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,map)
	eqns := [];
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
	Append(~eqns,map(q14 - f0*q23));
	Append(~eqns,map(q45 + q13 + f1*q23 + f2*q24));
	Append(~eqns,map(q12 + q35 - f3*q24 - f4*q34));
	Append(~eqns,map(q25 + f5*q34));
	if eqns eq [0,0,0,0] then
		return true;
	else
		return false;
	end if;
end function;

//From comparing D6_basic and D6_hom below, we see that D6_basic is faster:

stress_test_D6 := function(lower,upper,boolean)
	for s1 := lower to upper do
		for s2 := lower to upper do
			for s3 := lower to upper do
				for s4 := lower to upper do
					for s5 := lower to upper do
						s := [[s1,s2,s3,s4,s5]];
						for t1 := lower to upper do
							for t2 := lower to upper do
								for t3 := lower to upper do
									for t4 := lower to upper do
										for t5 := lower to upper do
											t := [[t1,t2,t3,t4,t5]];
//											if D6_basic(s1,s2,s3,s4,s5,t1,t2,t3,t4,t5) ne D6_hom(s1,s2,s3,s4,s5,t1,t2,t3,t4,t5) then
//												print(s);
//												print(t);
//												print("ERROR");
//											end if;
											if boolean then
												D6_basic(s1,s2,s3,s4,s5,t1,t2,t3,t4,t5);
											else
												D6_hom(s1,s2,s3,s4,s5,t1,t2,t3,t4,t5);
											end if;
										end for;
									end for;
								end for;
							end for;
						end for;
					end for;
				end for;
			end for;
		end for;
	end for;
	return "done";
end function;



//Function evaluating basis equations and returning true if they all vanish

basis_polys := function(x,y,map)
	for i := 1 to 5 do
		for j := i+1 to 5 do
			g := map(quad_subdet(x,y,i,j));
			pij := x[1,i]*y[1,j] - x[1,j]*y[1,i];
			f := g - pij;
			if f eq 0 then
				continue;
			else
				return false;
			end if;
		end for;
	end for;
	segre_eval := [];
	Append(~segre_eval,map(segre(x)));
	Append(~segre_eval,map(matrix_mult(dual_mult(x,x),transpose(y))[1,1]));	
	Append(~segre_eval,map(matrix_mult(dual_mult(y,y),transpose(x))[1,1]));
	Append(~segre_eval,map(segre(y)));
	if segre_eval eq [0,0,0,0] then
		return true;
	else
		return false;
	end if;
end function;

find_exception := function(lower,upper)
	for g0 := lower to upper do
		for g1 := lower to upper do
			for g2 := lower to upper do
				for g3 := lower to upper do
					for g4 := lower to upper do
						for g5 := lower to upper do
							h := hom< P -> R | x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,g0,g1,g2,g3,g4,g5 >;
							for s2 := lower to upper do
								for s3 := lower to upper do
									for s4 := lower to upper do
										for s5 := lower to upper do
											s := [[0,s2,s3,s4,s5]];
											if h(segre(s)) ne 0 or s eq [[0,0,0,0,0]] then
												continue;
											end if;
											for t1 := lower to upper do
												for t2 := lower to upper do
													for t3 := lower to upper do
														for t4 := lower to upper do
															for t5 := lower to upper do
																t := [[t1,t2,t3,t4,t5]];
																if h(segre(t)) ne 0 or t eq [[0,0,0,0,0]] or t eq s then
																	continue;
																end if;
																if not D6_basic(0,s2,s3,s4,s5,t1,t2,t3,t4,t5,h) then
																	if basis_polys(s,t,h) then
																		print(s);
																		print(t);
																		print([g0,g1,g2,g3,g4,g5]);
																		print("not D6");
																	end if;
																end if;
															end for;
														end for;
													end for;
												end for;
											end for;
										end for;
									end for;
								end for;
							end for;
						end for;
					end for;
				end for;
			end for;
		end for;
	end for;
	return "done";
end function;




