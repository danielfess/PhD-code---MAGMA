//Playing with equations for f = x^5 + y^5

Q := RationalField();
P<p23,p34,p13,p35,p15,p24> := PolynomialRing(Q,6);

eqns := [];
Append(~eqns,-p13*p24+p35*p23+p34^2);
Append(~eqns,p15*p23+p35*p34-p13^2);
Append(~eqns,p15*p24-p35*p13+p23*p34);
Append(~eqns,p15*p34+p13*p23-p35^2);
Append(~eqns,p13*p34-p35*p24+p23^2);
eqns2 := eqns;
Append(~eqns,p15^2 - 11*p15*p24 - p24^2 - 1);
Append(~eqns2,p15^2 - 11*p15*p24 - p24^2);

I := Ideal(eqns);
I;
IsPrime(I);
Groebner(I);

J := Ideal(eqns2);
J;
IsPrime(J);
