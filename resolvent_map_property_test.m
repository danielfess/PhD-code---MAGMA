//If g(x,y) = 0, what does this tell us about associated determinants?
//In particular, is the 3x3 determinant, formed in a similar manner but with x^2, xy, y^2, also zero?

Q := RationalField();
P<x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6> := PolynomialRing(Q,12);

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

alt := function(x,y)
	N := [[1,1,1],[x[1] + x[2], x[3] + x[6], x[4] + x[5]],[y[1] + y[2], y[3] + y[6], y[4] + y[5]]];
	return determinant(N);
end function;

product := function(x,y)
	prod := x;
	for i := 1 to 6 do
		prod[i] := x[i]*y[i];
	end for;
	return prod;
end function;

tr := function(x)
	sum := 0;
	for i := 1 to 6 do
		sum +:= x[i];
	end for;
	return sum;
end function;

tr0 := function(x)
	y := x;
	trace := tr(x);
	for i := 1 to 6 do
		y[i] := x[i] - trace/6;
	end for;
	return y;
end function;

//Playing with alt(x,y), alt(x,a) and associated determinants:
x := [x1,x2,x3,x4,x5,x6];
y := [y1,y2,y3,y4,y5,y6];
IsIrreducible(alt(x,y));

xx := product(x,x);
xy := product(product(x,y),[2,2,2,2,2,2]);
yy := product(y,y);

xx0 := tr0(xx);
xy0 := tr0(xy);
yy0 := tr0(yy);

v := [xx0,xy0,yy0];

T := [];
for i := 1 to 3 do
	Append(~T,[]);
	for j := 1 to 3 do
		Append(~T[i], tr(product(v[i],v[j])));
	end for;
end for;

cof := cofactor_matrix(T);

g := [[alt(x,xx0),alt(y,xx0)],[alt(x,xy0),alt(y,xy0)],[alt(x,yy0),alt(y,yy0)]];

h := matrix_mult(cof,g);
p1 := h[1,1] + h[2,2];
p2 := h[2,1] + h[3,2];

polys := [];
Append(~polys,tr(x));
Append(~polys,tr(y));
Append(~polys,alt(x,y));

J := Ideal(polys);
IsPrime(J);

Append(~polys,tr(product(x,xx)));
Append(~polys,tr(product(y,xx)));
Append(~polys,tr(product(x,yy)));
Append(~polys,tr(product(y,yy)));

I := Ideal(polys);
p1 in I;
p2 in I;

F := hom< P -> P | x1,x3,x4,x5,x6,x2,y1,y3,y4,y5,y6,y2 >;

tr_R := function(r)
	conjugate := r;
	sum := r;
	for i := 1 to 4 do
		conjugate := F(conjugate);
		sum +:= conjugate;
	end for;
	return sum;
end function;

r := [h[1,1],h[1,2],h[3,1],h[3,2]];

//T_R := [];
//for i := 1 to 4 do
//	Append(~T_R,[]);
//	for j := 1 to 4 do
//		print i,j;
//		Append(~T_R[i],tr_R(r[i]*r[j]));
//	end for;
//end for;

//D := determinant(T_R);

polys2 := polys[1..3];
alt_conjugate := F(alt(x,y));
for i := 1 to 4 do
	Append(~polys2,alt_conjugate);
	alt_conjugate := F(alt_conjugate);
end for;

alt(x,y) - alt_conjugate;

polys2;

K := Ideal(polys2);
L := PrimaryDecomposition(K)[1];
//Note: Second PrimaryDecomp factor is generated by pij = 0 and tr(x) = tr(y) = 0
//Can check that Segre polynomials are in L, but L does not equal I.
//This means that the additional alternating equations cut out a smaller subvariety than the trace cubed equations and alt(x,y).

Append(~polys2,tr(product(x,xx)));
Append(~polys2,tr(product(y,xx)));
Append(~polys2,tr(product(x,yy)));
Append(~polys2,tr(product(y,yy)));
O := Ideal(polys2);

O eq L;
//true
//so, the non-trivial solutions x,y to the alternating equations are those with all trace cubed equations = 0.

//Checking that alt(x,y) is stronger than segre:
segre_polys := [tr(x),tr(y),tr(product(x,xx)),tr(product(y,xx)),tr(product(x,yy)),tr(product(y,yy))];
I_segre := Ideal(segre_polys);
alt(x,y) in I_segre;
//False

