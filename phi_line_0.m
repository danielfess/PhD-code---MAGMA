//Given a line l, \Phi vanishing on l gives 4 linear equations in pluecker coordinates.
//In \mathbb{P}^9, these 4 linear equations plus the pluecker equations cut out a
//2-dimensional subvariety of \mathbb{P}^9.

//Such lines automatically lie on the Segre cubic.  Could this be the del Pezzo D6?

//To find forms lying in the same orbit after applying \Phi, we then need to introduce
//the determinant 1 equations.  What variety does this result in?

Z := Integers();
Q := RationalField();
//P<p15,p13,p35,p23,p24,p34,f0,f1,f2,f3,f4,f5> := PolynomialRing(Q,12);
P<p15,p13,p35,p23,p24,p34> := PolynomialRing(Q,6);
//A<p15,p13,p35,p23,p24,p34> := AffineSpace(Z,6);
//A<p13,p35,p23,p34> := AffineSpace(Q,4);
//P := CoordinateRing(A);

//p15 := 11;
//p24 := 1;

f0 := 1;
f1 := 0;
f2 := 0;
f3 := 0;
f4 := 0;
f5 := 1;

p14 := f0*p23;
p45 := - p13 - f1*p23 - f2*p24;
p12 := - p35 + f3*p24 + f4*p34;
p25 := - f5*p34;

p := [[0,p12,p13,p14,p15],[-p12,0,p23,p24,p25],[-p13,-p23,0,p34,p35],[-p14,-p24,-p34,0,p45],[-p15,-p25,-p35,-p45,0]];

p;

pluecker := [];

for i := 1 to 2 do
	for j := i+1 to 3 do
		for k := j+1 to 4 do
			for l := k+1 to 5 do
				poly := p[i,j]*p[k,l] + p[j,k]*p[i,l] - p[i,k]*p[j,l];
				Append(~pluecker,poly);
			end for;
		end for;
	end for;
end for;

print(pluecker);

I := Ideal(pluecker);
//IsRadical(I);
IsPrime(I);

eqns := pluecker;
Append(~eqns, p15^2 - 11*p15*p24 - p24^2 - 1);

C := Scheme(A,eqns);

