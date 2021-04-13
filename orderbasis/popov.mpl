macro(AddMultiple = LinearAlgebra:-Modular:-AddMultiple):
macro(Create = LinearAlgebra:-Modular:-Create):

#
# Input: A - an n x m matrix filled with modp1 polynomials modulo p
#        s - a list of integers of length m (a shift)
#
# Output: (P,delta) where 
#          P - the s-Popov form of A
#          delta - the s-degrees of the rows of P
#
# Note: The definition of s-Popov form we are using is such
#       that if A is square and nonsingular, then P will be
#       be s-reduced with (row) leading matrix lower triangular,
#       and off-diagonal entries in each column strictly less than
#       the diagonal entry in the same column.  In other words,
#       the pivots will be on the diagonal of the matrix.
#
popov := proc(A,x,n,m,s,p)
   local i,j,c,d,dmax,P,Q,B,R,t,k,delta,l,ws,simple,AA,nmin,shift;

   # Add c x^e times row i to row j.
   # - variables p, B, P, and m are lexically scoped
   simple := proc(c,e,i,j)
      AddMultiple(p,c,B,P[j][3]..P[j][3],e*m+1..(e+P[i][2]+1)*m,
                      B,P[i][3]..P[i][3],1..(P[i][2]+1)*m,
                      B,P[j][3]..P[j][3],e*m+1..(e+P[i][2]+1)*m);
   end:

   nmin := min(op(s));
   shift := [seq(s[i]-nmin,i=1..m)];

   AA := applyShift(p,A,shift)[[seq(i,i=1..n)],[seq(m+1-j,j=1..m)]];

   dmax := -infinity;
   for i to n do 
      P[i] := [-1,-infinity,i];
      for j to m do
         if modp1(IsZero(AA[i,j]),p) then d := -infinity else d := modp1(Degree(AA[i,j]),p) fi;
         if d>P[i][2] then P[i] := [j,d,i] fi;
      od;
      dmax := max(dmax,P[i][2]);
   od;

   P := table(sort([seq(P[i],i=1..n)],(x,y)->evalb(x[1]<y[1] or (x[1]=y[1] and x[2]<y[2]))));
   ws := kernelopts(wordsize);
   if p<2^25 then
      B := Create(p,n,m*(dmax+1),float[8]);
   elif (ws=32 and p<2^16) or (ws=64 and p<2^32) then
      B := Create(p,n,m*(dmax+1),integer[]);
   else
      B := Create(p,n,m*(dmax+1),integer);
   fi;
   for d from 0 to dmax do for i to n do for j to m do 
      B[i,d*m+j] := modp1(Coeff(AA[i,j],d),p) 
   od od od;

   # Transform to weak-Popov form.
   do
      for i to n-1 while P[i][2]=-infinity or P[i][1]<>P[i+1][1] do od;
      if i=n then break fi;
      c := modp(-round(B[P[i+1][3],m*P[i+1][2]+P[i+1][1]])/round(B[P[i][3],m*P[i][2]+P[i][1]]),p);
      simple(c,P[i+1][2]-P[i][2],i,i+1);
      for d from P[i+1][2] by -1 to 0 do
         for j to m while B[P[i+1][3],d*m+j]=0 do od;
         if j<=m then break fi;
      od; if j<=m then P[i+1] := [j,d,P[i+1][3]] else P[i+1] := [-1,-infinity,P[i+1][3]] fi;
      P := table(sort([seq(P[i],i=1..n)],(x,y)->evalb(x[1]<y[1] or (x[1]=y[1] and x[2]<y[2]))));
   od;

   # Transform to Popov form.
   P := table(sort([seq(P[i],i=1..n)],(x,y)->evalb(x[2]<y[2] or (x[2]=y[2] and y[1]<x[1]))));
   for t to n while P[t][1]=-1 do od;
   for k from t to n do
      do
         delta := -infinity;
         for i from t to k-1 do
            for d from P[k][2] by -1 to 0 while B[P[k][3],d*m+P[i][1]]=0 do od;
            if d>=0 and d-P[i][2]>delta then l,delta,c := i,d-P[i][2],B[P[k][3],d*m+P[i][1]] fi;
         od;
         if delta>=0 then simple(modp(-round(c),p),delta,l,k) else break fi;
      od;
      c := modp(1/round(B[P[k][3],m*P[k][2]+P[k][1]]),p);
      simple(c-1,0,k,k);
   od;
   P := table(sort([seq(P[i],i=1..n)],(x,y)->evalb(x[1]>y[1])));

   R := Matrix(n,m);
   for i to n do for j to m do 
         R[i,j] := modp1(ConvertIn([seq(round(B[P[i][3],d*m+j]),d=0..P[i][2])],x),p);
   od od;

   P := map2(subsop,3=NULL,P);
   
   R := R[[seq(i,i=1..n)],[seq(m+1-j,j=1..m)]];
   
   return applyShift(p,R,-shift),[seq(P[i][2]+nmin,i=1..n)];

end:

applyShift := proc(p,A,shift)
   local i,j,B,n,m,x,s;

   n := LinearAlgebra:-RowDimension(A);
   m := LinearAlgebra:-ColumnDimension(A);
   if n>0 and m>0 then x := modp1(Indeterminate(A[1,1],p)) fi;
   B := Matrix(n,m):
   for i to n do for j to m do
      if shift[j]<0 then
         B[i,j] := modp1(Quo(A[i,j],ConvertIn(x^(-shift[j]),x)),p)
      else
         s := modp1(ConvertIn(x^(shift[j]),x),p);
         B[i,j] := modp1('Multiply'(A[i,j],s),p)
      fi
   od od;

   return eval(B);

end:


