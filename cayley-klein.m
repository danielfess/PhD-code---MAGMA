//Min poly of Cayley-Klein

Q := RationalField();
R<x1,x2,x3,x4,x5> := PolynomialRing(Q,5);
y1 := x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1 - x1*x3 - x3*x5 - x5*x2 - x2*x4 - x4*x1;
y2 := x1*x3 + x3*x2 + x2*x5 + x5*x4 + x4*x1 - x1*x2 - x2*x4 - x4*x3 - x3*x5 - x5*x1;
print y2;
h := hom< R -> R | x2,x3,x4,x5,x1 >;
print y2;
y3 := h(y2);
print y2;
y4 := h(y3);
y5 := h(y4);
y6 := h(y5);

Disc := (y1-y2)^2*(y1-y3)^2*(y1-y4)^2*(y1-y5)^2*(y1-y6)^2*(y2-y3)^2*(y2-y4)^2*(y2-y5)^2*(y2-y6)^2*(y3-y4)^2*(y3-y5)^2*(y3-y6)^2*(y4-y5)^2*(y4-y6)^2*(y5-y6)^2;
Disc2 := (x1-x2)^2*(x1-x3)^2*(x1-x4)^2*(x1-x5)^2*(x2-x3)^2*(x2-x4)^2*(x2-x5)^2*(x3-x4)^2*(x3-x5)^2*(x4-x5)^2;

#(Disc^Sym(5));

S<e1,e2,e3,e4,e5> := PolynomialRing(Q,5);
l1, p1 := IsSymmetric(y1+y2+y3+y4+y5+y6,S);
l2, p2 := IsSymmetric(y1*y2+y1*y3+y1*y4+y1*y5+y1*y6+y2*y3+y2*y4+y2*y5+y2*y6+y3*y4+y3*y5+y3*y6+y4*y5+y4*y6+y5*y6,S);
l3, p3 := IsSymmetric(y1^3+y2^3+y3^3+y4^3+y5^3+y6^3,S);
l4, p4 := IsSymmetric(y1^4+y2^4+y3^4+y4^4+y5^4+y6^4,S);
l5, p5 := IsSymmetric(y1^5+y2^5+y3^5+y4^5+y5^5+y6^5,S);
l6, p6 := IsSymmetric(y1*y2*y3*y4*y5*y6,S);

//Product of images of x and x^{-1}

cayley := function(x1,x2,x3,x4,x5)
	return x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1 - x1*x3 - x3*x5 - x5*x2 - x2*x4 - x4*x1;
end function;

//cayley_inverse applies cayley to x^-1 but having cancelled denominator.
cayley_inverse := function(x1,x2,x3,x4,x5)
	return x3*x4*x5 + x4*x5*x1 + x5*x1*x2 + x1*x2*x3 + x2*x3*x4 - x2*x4*x5 - x4*x1*x2 - x1*x3*x4 - x3*x5*x1 - x5*x2*x3;
end function;