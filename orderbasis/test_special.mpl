read "special.mpl";
read "../maple_examples/modpops.mpl";

printlevel := -1;
p := 101;
n := 5;
s := [0, 0, 0, 0, 0];
m := 1;
d := 10;
deg := 9;

k := 1;
l := 1;

A := Matrix(n, m);
	
for i to n do 
  for j to m do
    A[i, j] := modp1(ConvertIn(randpoly(x, degree=deg), x), p) 
  od;
od;

print("A = ");
print(A);

M, N, delta := MinApproxBasis(A, x, n, m, n, n, d, s, p);
print("M A = 0 mod x^d in Z_p[x], M = ");
print(M);

MA := Matrix(n, m);
modpmult(M, A, x, n, n, m, p, MA);
print("MA = 0 mod x^d in Z_p[x] =");
print(MA);

deg_M := Matrix(n, n);

for i to n do
  for j to n do
    deg_M[i, j] := modp1(Degree(M[i,j]), p)
  od
od;

print("deg_M = ");
print(deg_M);

det_M := modp1(Det(M), p);

MN := Matrix(n, n);
modpmult(M, N, x, n, n, n, p, MN);
print("N = adj(M) =");
print(N);
print("MN = M N = det(M) I = ");
print(MN);
print("det(M) = ");
print(det_M);


M2, N2, delta2 := MinApproxBasis(A, x, n, m, k, l, d, s, p);

print("M2 = M[1..n, 1..k] = ");
print(M2);
if submatrix(M, M2, n, n, n, k) then
  print("M2 is the first k columns of M");
else print("ERROR: M2 is NOT the first k columns of M!"); fi;

print("N2 = N[1..l, 1..n] = ");
print(N2);
if submatrix(N, N2, n, n, l, n) then
  print("N2 is the first l rows of N");
else print("ERROR: N2 is NOT the first l columns of N!"); fi;
print("delta = ", delta2);
printlevel := 0;

