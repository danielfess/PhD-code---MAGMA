//Equations for quintics landing in same orbit

//Structure of code:
//We work in the sextic resolvent ring S_f.
//Say x = \sum x_i \beta_i^*, similarly for y in terms of yi.
//Then we look for c,d,e in \mathbb{Z}^5 such that the matrix (x,c,d,e,y) is in GL5(\mathbb{Z})
//and c,d,e are orthonormal wrt 8 \Delta(f) (x^2,-2xy,y^2).

//These are necessary conditions for being of the form \Phi(g)

//Note: Trivial case is x = e_5, y = e_1, c = e_2, d = e_3, e = e_4 (e_i are unit vectors).
//Less non-trivial but still general case comes from GL2 action.

Q := RationalField();
P<f0,f1,f2,f3,f4,f5,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5> := PolynomialRing(Q,16);

F1 := [[0,0,0,1,0],[0,-f0,0,-f2/2,0],[0,0,0,0,-1/2],[1,-f2/2,0,-f4,0],[0,0,-1/2,0,0]];
F2 := [[0,-f0,0,-f2/2,0],[-f0,3*f0*f2,f0*f3,f1*f3,-f1],[0,f0*f3,f0*f4,(f1*f4-f0*f5)/2,0],[-f2/2,f1*f3,(f1*f4-f0*f5)/2,f2*f4,-f3/2],[0,-f1,0,-f3/2,1]];
F3 := [[0,0,0,0,-1/2],[0,f0*f3,f0*f4,(f1*f4-f0*f5)/2,0],[0,f0*f4,3*f0*f5,f1*f5,0],[0,(f1*f4-f0*f5)/2,f1*f5,f2*f5,0],[-1/2,0,0,0,0]];
F4 := [[1,-f2/2,0,-f4,0],[-f2/2,f1*f3,(f1*f4-f0*f5)/2,f2*f4,-f3/2],[0,(f1*f4-f0*f5)/2,f1*f5,f2*f5,0],[-f4,f2*f4,f2*f5,3*f3*f5,-f5],[0,-f3/2,0,-f5,0]];
F5 := [[0,0,-1/2,0,0],[0,-f1,0,-f3/2,1],[-1/2,0,0,0,0],[0,-f3/2,0,-f5,0],[0,1,0,0,0]];

F := [F1,F2,F3,F4,F5];

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

dual_mult := function(a,b)
	product := [[P | 0,0,0,0,0]];
	for k := 1 to 5 do
		product[1,k] := matrix_mult(matrix_mult(a,F[k]),transpose(b))[1,1];
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
	det := -2*determinant(B);
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

Append(~eqns,q45 - f0*q23);
Append(~eqns,q14 - q35 - f1*q23 - f2*q24);
Append(~eqns,q25 - q13 - f3*q24 - f4*q34);
Append(~eqns,q12 - f5*q34);

J := Ideal(eqns);

//Append(~eqns,q24^2 + 11*q15*q24 - q15^2 - 1);

I := Ideal(eqns);
IsPrime(I);

quad_subdet(1,5) in I;

for i := 1 to 5 do
	for j := i+1 to 5 do
		print i,j;
		g := quad_subdet(i,j);
		pij := x[1,i]*y[1,j] - x[1,j]*y[1,i];
		f := g - pij;
		if not IsIrreducible(f) then
			print "reducible";
		end if;
		if f in I then
			print "in I";
		end if;
		Append(~eqns, f);
//		h := quad_subdet(i,j) - pij*(q24^2 + 11*q15*q24 - q15^2);
//		if h in J then
//			print "quad(i,j) = pij*unit eqn mod pluecker";
//		end if;
	end for;
end for;

K := Ideal(eqns);
IsPrime(K);
RadicalDecomposition(K);

//Note: For f = x^5 + y^5, the polynomials p_ij - quad_subdet(i,j) lie in the ideal generated by
//the four linear equations and the unit equation!
//Some subtlety in the signs of +-1 in unit equation and +-1 in front of quad_subdet(i,j).
//They need to be the same.  This must reflect the sign of det(M).

//Idea: Compute the quad_subdet's on the open set p12 \neq 0, map to the polynomial ring
//generated by the other pij, and rehomogenise using p12 to obtain a candidate for quad_subdet.

S<g0,g1,g2,g3,g4,g5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45> := PolynomialRing(Q,16);
h := hom< P -> S | g0,g1,g2,g3,g4,g5,1,0,-p23,-p24,-p25,0,1,p13,p14,p15,0>;
g_hom := hom< S -> P | f0,f1,f2,f3,f4,f5,q12,q13,q14,q15,q23,q24,q25,q34,q35,q45>;





