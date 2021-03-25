//Take a basis with the x,y property
//We know that the matrix A is of the form A_15 = 0, the top row
//looks like (0,a,b,c,0) and the bottom row like (0,-d,-a,-b,0)
//We also know the derivatives of Segre evaluated at y1 e1 + y5 e5,
//because Segre relates to the squaring map from S^* to S, and
//we have the x,y property.

//Do these derivatives tell us that a,b,c,d form a basis of \mathbb{Z}^4?


//Matrix functions:
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

//Segre setup and calculations:

Q := RationalField();
P<a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,c1,c2,c3,c4,c5,d1,d2,d3,d4,d5,e1,e2,e3,e4,e5,f1,f2,f3,f4,f5,z1,z2,z3,z4,z5,y1,y5> := PolynomialRing(Q,37);

basis_vector := function(i)
	v := [0,0,0,0,0];
	for k := 1 to 5 do
		if k eq i then
			v[k] := 1;
		end if;
	end for;
	return v;
end function;

A1 := [[0,a1,b1,c1,0],[-a1,0,0,0,-d1],[-b1,0,0,0,-e1],[-c1,0,0,0,-f1],[0,d1,e1,f1,0]];
A2 := [[0,a2,b2,c2,0],[-a2,0,0,0,-d2],[-b2,0,0,0,-e2],[-c2,0,0,0,-f2],[0,d2,e2,f2,0]];
A3 := [[0,a3,b3,c3,0],[-a3,0,0,0,-d3],[-b3,0,0,0,-e3],[-c3,0,0,0,-f3],[0,d3,e3,f3,0]];
A4 := [[0,a4,b4,c4,0],[-a4,0,0,0,-d4],[-b4,0,0,0,-e4],[-c4,0,0,0,-f4],[0,d4,e4,f4,0]];
A := [A1,A2,A3,A4];

y := [[y1,0,0,0,y5]];
z := [[z1,z2,z3,z4,z5]];

//segre_matrix_ij computes the j-th matrix
//in the computation of segre_dyi below
segre_matrix_ij := function(i,j)
	M := [];
	for k := 1 to 4 do
		if k eq j then
			Append(~M,matrix_mult([basis_vector(i)],A[k])[1]);
		else
			Append(~M,matrix_mult(y,A[k])[1]);
		end if;
	end for;
	Append(~M,z[1]);
	return M;
end function;

//segre_di calculates i-th derivative times linear product <y,z>
//evaluated at y1 e1 + y5 e5
segre_di := function(i)
	sum := 0;
	for j := 1 to 4 do
		sum +:= determinant(segre_matrix_ij(i,j));
	end for;
	return sum;
end function;

segre_di_factor := function(i)
	if segre_di(i) eq 0 then
		return 0;
	end if;
	return Factorisation(segre_di(i))[2,1];
end function;

matrix_Z4 := [[a1,a2,a3,a4],[b1,b2,b3,b4],[c1,c2,c3,c4],[d1,d2,d3,d4]];
d := determinant(matrix_Z4);

//segre_di(2) - y1^2*(y1*z1 + y5*z5)*d;
//segre_di(3) - y1*y5*(y1*z1 + y5*z5)*d;
//segre_di(4) - y5^2*(y1*z1 + y5*z5)*d;

a := [a1,a2,a3,a4];
b := [b1,b2,b3,b4];
c := [c1,c2,c3,c4];
d := [d1,d2,d3,d4];
e := [e1,e2,e3,e4];
f := [f1,f2,f3,f4];

//These coefficients of Segre are zero.
//Calculating them as determinant of a,b,c,d,e,f
Coefficients(segre_di_factor(2),y1)[1] - y5^2*determinant([a,d,e,f]);
Coefficients(segre_di_factor(2),y1)[2] + y5*determinant([a,b,d,f]) - y5*determinant([a,c,d,e]);
Coefficients(segre_di_factor(3),y1)[1] - y5^2*determinant([b,d,e,f]);
Coefficients(segre_di_factor(3),y1)[3] - determinant([a,b,c,e]);
Coefficients(segre_di_factor(4),y1)[2] + y5*determinant([a,c,e,f]) - y5*determinant([b,c,d,f]);
Coefficients(segre_di_factor(4),y1)[3] - determinant([a,b,c,f]);

//These are the non-zero coefficients.
//Calculating them as determinant of a,b,c,d,e,f
Coefficients(segre_di_factor(2),y1)[3] - determinant([a,b,c,d]);
Coefficients(segre_di_factor(3),y1)[2] + y5*determinant([a,b,e,f]) - y5*determinant([b,c,d,e]);
Coefficients(segre_di_factor(4),y1)[1] - y5^2*determinant([c,d,e,f]);

I := Ideal([determinant([a,d,e,f]), - determinant([a,b,d,f]) + determinant([a,c,d,e]), determinant([b,d,e,f]), determinant([a,b,c,e]), - determinant([a,c,e,f]) + determinant([b,c,d,f]), determinant([a,b,c,f]), determinant([a,b,c,d]) - 1, - determinant([a,b,e,f]) + determinant([b,c,d,e]) - 1, determinant([c,d,e,f]) - 1]);

segre_di(2) - (y1*z1+y5*z5)*( - y1^2*determinant([a,b,c,d]) + y1*y5*(determinant([a,b,d,f]) - determinant([a,c,d,e])) - y5^2*determinant([a,d,e,f]));
segre_di(3) - (y1*z1+y5*z5)*( - y1^2*determinant([a,b,c,e]) + y1*y5*(determinant([a,b,e,f]) - determinant([b,c,d,e])) - y5^2*determinant([b,d,e,f]));
segre_di(4) - (y1*z1+y5*z5)*( - y1^2*determinant([a,b,c,f]) + y1*y5*(determinant([a,c,e,f]) - determinant([b,c,d,f])) - y5^2*determinant([c,d,e,f]));

