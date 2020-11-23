Q := RationalField();
//P<f0,f1,f2,f3,f4,f5,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45> := PolynomialRing(Q,26);
P<x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,p12,p13,p14,p15,p23,p24,p25,p34,p35,p45> := PolynomialRing(Q,20);

f0 := 1;
f1 := 0;
f2 := 0;
f3 := 0;
f4 := 0;
f5 := 1;

F1 := f0*f2*y2^3 + f0*f3*y2^2*y3 + f0*f4*y2*y3^2 - f0*f5*y2*y3*y4 + f0*f5*y3^3 - f0*y1*y2^2 + f1*f3*y2^2*y4 + f1*f4*y2*y3*y4 + f1*f5*y3^2*y4 - f1*y2^2*y5 + f2*f4*y2*y4^2 + f2*f5*y3*y4^2 - f2*y1*y2*y4 + f3*f5*y4^3 - f3*y2*y4*y5 - f4*y1*y4^2 - f5*y4^2*y5 + y1^2*y4 - y1*y3*y5 + y2*y5^2;
F2 := 3*f0*f2*x2*y2^2 + 2*f0*f3*x2*y2*y3 + f0*f3*x3*y2^2 + f0*f4*x2*y3^2 + 2*f0*f4*x3*y2*y3 - f0*f5*x2*y3*y4 - f0*f5*x3*y2*y4 + 3*f0*f5*x3*y3^2 - f0*f5*x4*y2*y3 - f0*x1*y2^2 - 2*f0*x2*y1*y2 + 2*f1*f3*x2*y2*y4 + f1*f3*x4*y2^2 + f1*f4*x2*y3*y4 + f1*f4*x3*y2*y4 + f1*f4*x4*y2*y3 + 2*f1*f5*x3*y3*y4 + f1*f5*x4*y3^2 - 2*f1*x2*y2*y5 - f1*x5*y2^2 + f2*f4*x2*y4^2 + 2*f2*f4*x4*y2*y4 + f2*f5*x3*y4^2 + 2*f2*f5*x4*y3*y4 - f2*x1*y2*y4 - f2*x2*y1*y4 - f2*x4*y1*y2 + 3*f3*f5*x4*y4^2 - f3*x2*y4*y5 - f3*x4*y2*y5 - f3*x5*y2*y4 - f4*x1*y4^2 - 2*f4*x4*y1*y4 - 2*f5*x4*y4*y5 - f5*x5*y4^2 + 2*x1*y1*y4 - x1*y3*y5 + x2*y5^2 - x3*y1*y5 + x4*y1^2 - x5*y1*y3 + 2*x5*y2*y5;
F3 := 3*f0*f2*x2^2*y2 + f0*f3*x2^2*y3 + 2*f0*f3*x2*x3*y2 + 2*f0*f4*x2*x3*y3 + f0*f4*x3^2*y2 - f0*f5*x2*x3*y4 - f0*f5*x2*x4*y3 + 3*f0*f5*x3^2*y3 - f0*f5*x3*x4*y2 - 2*f0*x1*x2*y2 - f0*x2^2*y1 + f1*f3*x2^2*y4 + 2*f1*f3*x2*x4*y2 + f1*f4*x2*x3*y4 + f1*f4*x2*x4*y3 + f1*f4*x3*x4*y2 + f1*f5*x3^2*y4 + 2*f1*f5*x3*x4*y3 - f1*x2^2*y5 - 2*f1*x2*x5*y2 + 2*f2*f4*x2*x4*y4 + f2*f4*x4^2*y2 + 2*f2*f5*x3*x4*y4 + f2*f5*x4^2*y3 - f2*x1*x2*y4 - f2*x1*x4*y2 - f2*x2*x4*y1 + 3*f3*f5*x4^2*y4 - f3*x2*x4*y5 - f3*x2*x5*y4 - f3*x4*x5*y2 - 2*f4*x1*x4*y4 - f4*x4^2*y1 - f5*x4^2*y5 - 2*f5*x4*x5*y4 + x1^2*y4 - x1*x3*y5 + 2*x1*x4*y1 - x1*x5*y3 + 2*x2*x5*y5 - x3*x5*y1 + x5^2*y2;
F4 := f0*f2*x2^3 + f0*f3*x2^2*x3 + f0*f4*x2*x3^2 - f0*f5*x2*x3*x4 + f0*f5*x3^3 - f0*x1*x2^2 + f1*f3*x2^2*x4 + f1*f4*x2*x3*x4 + f1*f5*x3^2*x4 - f1*x2^2*x5 + f2*f4*x2*x4^2 + f2*f5*x3*x4^2 - f2*x1*x2*x4 + f3*f5*x4^3 - f3*x2*x4*x5 - f4*x1*x4^2 - f5*x4^2*x5 + x1^2*x4 - x1*x3*x5 + x2*x5^2;
G12 := p12 - x1*y2 + x2*y1;
G13 := p13 - x1*y3 + x3*y1;
G14 := p14 - x1*y4 + x4*y1;
G15 := p15 - x1*y5 + x5*y1;
G23 := p23 - x2*y3 + x3*y2;
G24 := p24 - x2*y4 + x4*y2;
G25 := p25 - x2*y5 + x5*y2;
G34 := p34 - x3*y4 + x4*y3;
G35 := p35 - x3*y5 + x5*y3;
G45 := p45 - x4*y5 + x5*y4;

line_polys := [P | ];
Append(~line_polys,F1);
Append(~line_polys,F2);
Append(~line_polys,F3);
Append(~line_polys,F4);
Append(~line_polys,G12);
Append(~line_polys,G13);
Append(~line_polys,G14);
Append(~line_polys,G15);
Append(~line_polys,G23);
Append(~line_polys,G24);
Append(~line_polys,G25);
Append(~line_polys,G34);
Append(~line_polys,G35);
Append(~line_polys,G45);


I := Ideal(line_polys);
print(I);
PrimaryDecomposition(I);



