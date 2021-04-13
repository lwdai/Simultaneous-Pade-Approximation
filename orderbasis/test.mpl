#read `minapproxbasis.mpl`;
read "recur.mpl";
read "../maple_examples/modpops.mpl";

printlevel := -1;
p := 101;
n := 10;
#s := [10,9,8,7,6,5,4,3,2,1];
s := [0,0,0,0,0,0,0,0,0,0];
#s := [4,3,2,1];
m := 2;
d := 100;
deg := 99;

k := n;
l := n;

A := Matrix(n, m);
	
for i to n do 
  for j to m do
    A[i, j] := modp1(ConvertIn(randpoly(x, degree=deg), x), p) 
  od;
od;

print("A = ");
print(A);

M,news := MinApproxBasis(A, x, n, m, d, s, p):

print("M A = 0 mod x^d in Z_p[x], M = ");
#print(M);


deg_M := Matrix(n, n);



for i to n do
  for j to n do
    deg_M[i, j] := modp1(Degree(M[i,j]), p)
  od
od;

print("deg_M = ");
print(deg_M);

MA := Matrix(n, m);

modpmult(M, A, x, n, n, m, p, MA); 

for i to n do
  for j to m do
    MA[i,j] := modp1(Rem(MA[i,j], modp1(ConvertIn(x^d,x),p)),p)
  od
od;

# 
print("MA = 0 mod x^d in Z_p[x] =");
print(MA);
printlevel := 0;

#M;
s;
news;
