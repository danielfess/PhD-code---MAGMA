//Trying some calculations to see how Cayley-Klein and fundamental resolvent maps
//behave with each other.

Q := RationalField();
P<u1,u2,u3,u4,v1,v2,v3,v4,t1,t2,t3,t4,f0,f1,f2,f3,f4,f5> := PolynomialRing(Q,18);

cayley := function(s1,s2,s3,s4)
	Q1 := -f0*s1^2 - f1*s1*s2 - f2*s2^2 - f3*s2*s3 - f4*s3^2 - f5*s3*s4;
	Q2 := s2^2 - s1*s3;
	Q3 := s1*s4 - s2*s3;
	Q4 := s3^2 - s2*s4;
	Q5 := - f0*s1*s2 - f1*s2^2 - f2*s2*s3 - f3*s3^2 - f4*s3*s4 - f5*s4^2;
	Q := [Q1,Q2,Q3,Q4,Q5];
	return Q;
end function;

print cayley(u1,u2,u3,u4);

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

A := [cayley(u1,u2,u3,u4)];
B := [[0,t4,-t3,t2,0],]