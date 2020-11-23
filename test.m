//MAGMA test

x := 1;
f := function(x,y)
	return x+y;
end function;

print f(2,2);

//So MAGMA creates local variables for numbers in definitions.

x := [3,4];
g1 := function(x1)
	sum := 0;
	for num in x do
		sum +:= num;
	end for;
	return sum;
end function;

print g1([1,1,2]);

g2 := function(x)
	sum := 0;
	for num in x do
		sum +:= num;
	end for;
	return sum;
end function;

print g2([1,1,2]);

//So MAGMA creates local variables for lists in definitions,
//and will use global variables if the variable name does
//not appear in the definition.

x := [[3,4],[1]];
h := procedure(~x,y)
	x[1,2] := y;
end procedure;

h(~x,2);
print x;

//Above we see how to change global variables.
//Use a tilde in the function/procedure definition and
//in the function/procedure call.