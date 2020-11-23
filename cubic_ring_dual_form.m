//Given a cubic ring R = R(f), compute the form Tr( (x \omega^* + y \theta^*)^3 ).
//This should be a covariant of f.

Q := RationalField();
P<a,b,c,d,x,y> := PolynomialRing(Q,6);

m1 := [[1,0,0],[0,1,0],[0,0,1]];
m2 := [[0,-a*c,-a*d],[1,b,0],[0,-a,0]];
m3 := [[0,-a*d,-b*d],[0,0,d],[1,0,-c]];
mult_coeffs := [m1,m2,m3];

//Trace of alpha_i alpha_j, indexed from 1 to 3
trace := function(i,j)
	tr := 0;
	for k := 1 to 3 do
		for l := 1 to 3 do
			tr +:= mult_coeffs[i,k,j]*mult_coeffs[k,l,l];
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

matrix_add := function(A,B)
	sum := [];
	for i := 1 to #A do
		Append(~sum,[]);
		for j := 1 to #A[1] do
			Append(~sum[i], A[i,j] + B[i,j]);
		end for;
	end for;
	return sum;
end function;

matrix_scalar := function(A,k)
	product := [];
	for i := 1 to #A do
		Append(~product,[]);
		for j := 1 to #A[1] do
			Append(~product[i], k*A[i,j]);
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
for i := 1 to 3 do
	Append(~T,[]);
	for j := 1 to 3 do
		Append(~T[i],trace(i,j));
	end for;
end for;
print T;

disc := determinant(T);

T_inv := cofactor_matrix(T);

b2 := [T_inv[2]];
b3 := [T_inv[3]];
elt := matrix_add(matrix_scalar(b2,x),matrix_scalar(b3,y));

m := matrix_scalar(m1,elt[1,1]);
m := matrix_add(m,matrix_scalar(m2,elt[1,2]));
m := matrix_add(m,matrix_scalar(m3,elt[1,3]));

m_cubed := matrix_mult(m,matrix_mult(m,m));

g := m_cubed[1,1] + m_cubed[2,2] + m_cubed[3,3];
g;

