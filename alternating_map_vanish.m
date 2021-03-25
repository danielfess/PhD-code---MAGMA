//If g(x,y) = 0, what does this tell us about associated determinants?
//In particular, is the 3x3 determinant, formed in a similar manner but with x^2, xy, y^2, also zero?

Q := RationalField();
R<D> := FunctionField(Q);
P<x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c2,c3,c4,c5,c6> := PolynomialRing(R,30);

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

M := [[x1^2 + x2^2, x3^2 + x6^2, x4^2 + x5^2],[x1*y1 + x2*y2, x3*y3 + x6*y6, x4*y4 + x5*y5],[y1^2 + y2^2, y3^2 + y6^2, y4^2 + y5^2]];
det := determinant(M);
Factorization(det);
//So, this determinant is irreducible, hence in general non-zero even when alt(x,y) = 0.

alt := function(x,y)
	N := [[1,1,1],[x[1] + x[2], x[3] + x[6], x[4] + x[5]],[y[1] + y[2], y[3] + y[6], y[4] + y[5]]];
	return determinant(N);
end function;

//Playing with alt(x,y), alt(x,a) and associated determinants:
x := [x1,x2,x3,x4,x5,x6];
y := [y1,y2,y3,y4,y5,y6];
IsIrreducible(alt(x,y));

G := function(a,b)
	return a[1]*b[2] + a[2]*b[1] + a[3]*b[6] + a[6]*b[3] + a[4]*b[5] + a[5]*b[4];
end function;

xx := [x1^2,x2^2,x3^2,x4^2,x5^2,x6^2];
xy := [x1*y1,x2*y2,x3*y3,x4*y4,x5*y5,x6*y6];
yy := [y1^2,y2^2,y3^2,y4^2,y5^2,y6^2];
A := [A1,A2,A3,A4,A5,A6];

L := [[x1^2 + x2^2 + x3^2 + x6^2 + x4^2 + x5^2, x1*y1 + x2*y2 + x3*y3 + x6*y6 + x4*y4 + x5*y5, y1^2 + y2^2 + y3^2 + y6^2 + y4^2 + y5^2],[G(x,xx),G(x,xy),G(x,yy)],[G(A,xx),G(A,xy),G(A,yy)]];
h := determinant(L);
IsIrreducible(h);

//Is alt(x,a) = alt(y,b) etc. an algebraic consequence of the other conditions on x,y,a,b,c?
//i.e. alt(x,y) = 0, and a,b,c 'orthonormal' to x^2,xy,y^2.

polys := [];
Append(~polys,alt(x,y));

a := [a1,a2,a3,a4,a5,a6];
b := [b1,b2,b3,b4,b5,b6];
c := [c1,c2,c3,c4,c5,c6];

//trace of a*b (sum over conjugates)
tr := function(a,b)
	sum := 0;
	for i := 1 to 6 do
		sum +:= a[i]*b[i];
	end for;
	return sum;
end function;

print(tr(a,b));

Append(~polys,8*D*tr(a,xx) - 1);
Append(~polys,-16*D*tr(a,xy));
Append(~polys,8*D*tr(a,yy));
Append(~polys,8*D*tr(b,xx));
Append(~polys,-16*D*tr(b,xy) - 1);
Append(~polys,8*D*tr(b,yy));
//Append(~polys,8*D*tr(c,xx));
//Append(~polys,-16*D*tr(c,xy));
//Append(~polys,8*D*tr(c,yy) - 1);

I := Ideal(polys);

alt(x,a) - alt(y,b) in I;
