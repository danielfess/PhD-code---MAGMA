//Calculating the equation of the trace cubed = 0 variety inside \tilde{S}.

Q := RationalField();

P<t1,t2,t3,t4,f0,f1,f2,f3,f4,f5> := PolynomialRing(Q,10);

//Cayley map for \Phi'(f).
//According to HCL IV, should be unsigned sub-pfaffians, but signed sub-pfaffians are what work. 
cayley := function(s1,s2,s3,s4)
	Q1 := -f0*s1^2 - f1*s1*s2 - f2*s2^2 - f3*s2*s3 - f4*s3^2 - f5*s3*s4;
	Q2 := s2^2 - s1*s3;
	Q3 := s1*s4 - s2*s3;
	Q4 := s3^2 - s2*s4;
	Q5 := - f0*s1*s2 - f1*s2^2 - f2*s2*s3 - f3*s3^2 - f4*s3*s4 - f5*s4^2;
	Q := [Q1 + 2/5*f2*Q2 - f3/10*Q3 + 2/5*f4*Q4,-Q2,Q3,-Q4,Q5 + 2/5*f1*Q2 - f2/10*Q3 + 2/5*f3*Q4];
	return Q;
end function;

Q := cayley(t1,t2,t3,t4);
y1 := Q[1];
y2 := Q[2];
y3 := Q[3];
y4 := Q[4];
y5 := Q[5];

trace_cubed2 := 60*f0*f2*y2^3 + 90*f0*f3*y2^2*y3 - 40*f0*f4*y2^2*y4 + 100*f0*f4*y2*y3^2 - 100*f0*f5*y2*y3*y4 + 100*f0*f5*y3^3 - 100*f0*y1*y2^2 - 24*f1^2*y2^3 - 18*f1*f2*y2^2*y3 + 52*f1*f3*y2^2*y4 - 4*f1*f3*y2*y3^2 + 84*f1*f4*y2*y3*y4 - 40*f1*f5*y2*y4^2 + 100*f1*f5*y3^2*y4 - 40*f1*y1*y2*y3 - 20*f1*y2^2*y5 - 24*f2^2*y2^2*y4 - 3*f2^2*y2*y3^2 - 20*f2*f3*y2*y3*y4 - f2*f3*y3^3 + 52*f2*f4*y2*y4^2 - 4*f2*f4*y3^2*y4 + 90*f2*f5*y3*y4^2 - 20*f2*y1*y2*y4 - 10*f2*y1*y3^2 - 20*f2*y2*y3*y5 - 24*f3^2*y2*y4^2 - 3*f3^2*y3^2*y4 - 18*f3*f4*y3*y4^2 + 60*f3*f5*y4^3 - 20*f3*y1*y3*y4 - 20*f3*y2*y4*y5 - 10*f3*y3^2*y5 - 24*f4^2*y4^3 - 20*f4*y1*y4^2 - 40*f4*y3*y4*y5 - 100*f5*y4^2*y5 + 100*y1^2*y4 - 100*y1*y3*y5 + 100*y2*y5^2;
print trace_cubed2;

segre := f0*f2*y2^3 + f0*f3*y2^2*y3 + f0*f4*y2*y3^2 - f0*f5*y2*y3*y4 + f0*f5*y3^3 - f0*y1*y2^2 + f1*f3*y2^2*y4 + f1*f4*y2*y3*y4 + f1*f5*y3^2*y4 - f1*y2^2*y5 + f2*f4*y2*y4^2 + f2*f5*y3*y4^2 - f2*y1*y2*y4 + f3*f5*y4^3 - f3*y2*y4*y5 - f4*y1*y4^2 - f5*y4^2*y5 + y1^2*y4 - y1*y3*y5 + y2*y5^2;

//I := Ideal([trace_cubed2,segre]);
