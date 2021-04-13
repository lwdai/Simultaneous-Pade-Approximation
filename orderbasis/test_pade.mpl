read "pade.mpl";

printlevel := -1;
n := 4;
d := 7;
g := [ seq(modp1(ConvertIn(x^d,x),p), i=1..n) ];
p := 2;
#N := [5,3,4,5];
roll := rand(d-3..d);
N := [ seq(roll(), i=1..n+1) ];
#S := [x^4+x^2+1, x^4+1, x^4+x^3+1];

S:= [ seq(randpoly(x, degree = d-1), i=1..n) ];

S := map(y -> modp1(ConvertIn(y,x),p), S);

print("p = ", p);
print("S = ", S);
print("N = ", N);

lambda1, delta1 := RecursiveSimPade(S,x,n,g,p,N);

print("lambda0 = ", lambda1);
print("delta0 = ", delta1);

lambda1,delta1 := SimPade(S,x,n,d,p,N);

print("lambda1 = ", lambda1);
print("delta1 = ", delta1);

lambda2,delta2 := SimPade_2(S,x,n,d,p,N);

print("lambda2 = ", lambda2);
print("delta2 = ", delta2);


printlevel := 0;
