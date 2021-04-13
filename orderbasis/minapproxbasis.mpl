#
# Input: A - an n x m modp1 matrix over Z_p[x]
#        d - desired order of the approximant
#        s - an integer tuple of size n (a shift)
#      
# Output: (F, delta) 
#
#          with
#      
#         delta = deg_s F 
# 
#           where
#
#         F is minimal approximant an order d with shift s
#

MinApproxBasis := proc(A,x,n,m,d,s,p)
   local F,i,j,ip,jp,delta,B,ord, pivot, mult, c1, delta_index;
   # c1 - temporarily stores the constant term at the pivot (row) in B
   # delta_min - delta[delta_min][2] is the pivoting row number

   # Initialize F to be an n x n indentity matrix of type modp1.
   F := Matrix(n,n):
   for i to n do for j to n do F[i,j] := modp1(Zero(x),p) od; F[i,i] := modp1(One(x),p) od:
 
   # Initialize list delta = [ [s[1],1], [s[2],2],...,s[n],n] ] of shifts.
   delta := [ seq([s[i],i], i=1..n) ];

   # Initialize B to be a copy of A.
   B := LinearAlgebra:-Copy(A);

   ord := 0;
   while ord < d do 
       # Sort delta in nondecreasing order according to first component.
       delta := sort(delta, proc(a,b) evalb(a[1] <= b[1]) end);		
	
       # Check if all constant coefficients of B are zero
       #  - if so, divide B by x, increase ord and continue
       for i to n do 
           for j to m do
	       c1 := modp1(Coeff(B[delta[i][2],j],0),p);
               if c1 <>0 then break fi;
           od;
           if j <= m then break fi;
       od;

       
       if i = n+1 then 
          # All constant coefficients are zero.
          # Divide B by x (use modp1/Shift), increase ord and continue.
          # AS use modp1/Shift
          for i to n do for j to m do
              B[i, j] := modp1(Shift(B[i, j], -1), p);
          od od;
          ord := ord + 1;
          next
       fi;


       # Have found i,j with i minimal such that constant coeff of B[delta[i][2],j].
       # - use B[delta[i][2],j] to zero out constant coefficients in column j in other rows of B
       # - multiply row i of F and B by x and decrease delta[i][1] by one
       delta_index := i;
       ip, jp := delta[i][2], j;
       mult := modp1(ConvertIn(modp(-1/c1, p), x), p);
       # now mult * c1 = -1 mod p.
               
       # AS Multiply row ip of B by mult
       for j to m do B[ip, j] := modp1(Multiply(B[ip, j], mult), p) od;
       # AS Multipy row ip of F by mult
       for j to n do F[ip, j] := modp1(Multiply(F[ip, j], mult), p) od;

       for i to n do	 
           if i = ip then next fi;
           mult := modp1(ConvertIn(modp1(Coeff(B[i,jp],0),p),x),p);
           # Add mult times row ip to row i of B.
           for j to m do B[i, j] := modp1(Add(B[i, j], modp1(Multiply(B[ip, j], mult), p)), p) od;
           # Add mult times row ip to row i of F
           for j to n do F[i, j] := modp1(Add(F[i, j], modp1(Multiply(F[ip, j], mult), p)), p) od;
       od;
       # Multiply row ip of B and F by x. (Use modp/Shift)
       for j to m do B[ip, j] := modp1(Shift(B[ip, j], 1), p) od;
       for j to n do F[ip, j] := modp1(Shift(F[ip, j], 1), p) od;	
       
       delta[delta_index][1] := delta[delta_index][1] + 1;


      od;

   delta := sort(delta, proc(a,b) evalb(a[2]<b[2]) end);
   delta := map(proc(a) a[1] end,delta);

   return F, delta;
 
end:
