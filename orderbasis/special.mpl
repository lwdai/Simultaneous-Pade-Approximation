#
# Input: A - an n x m modp1 matrix over Z_p[x]
#        d - desired order of the approximant
#        s - an integer tuple of size n (a shift)
#        k - number of initial columns of M to compute
#        l - number of initial rows of adj(M) to compute
#      
# Output: (M, N, delta) 
#
#          with
#      
#         M - the first k columns of F
#         N - the first l rows of adjoint(F)
#         delta - deg_s F 
# 
#           where
#
#         F is minimal approximant an order d with shift s
#
MinApproxBasis_2 := proc(A,x,n,m,k,l,d,s,p)

 # Initialization

 # M - the first k columns of I_n
 # N - first l rows of I_n
 # delta = s
 # B a copy of A
local I_n, M, N,i,j, ip, jp, delta,B,ord, pivot, mult, c1, delta_index; #factor;
   

   # Initialize I_n to be an n x n indentity matrix of type modp1.
   I_n := Matrix(n,n):
   for i to n do for j to n do I_n[i,j] := modp1(Zero(x),p) od; I_n[i,i] := modp1(One(x),p) od:

   M := Matrix(n, k):
   for i to n do for j to k do M[i, j] := I_n[i, j] od od:

   N := Matrix(l, n):
   for i to l do for j to n do N[i, j] := I_n[i, j] od od:
   

   # Initialize list delta = [ [s[1],1], [s[2],2],...,s[n],n] ] of shifts.    #flage here "s"
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
       # AS Multipy row ip of M by mult
       for j to k do M[ip, j] := modp1(Multiply(M[ip, j], mult), p) od;
       # Multiply all columns of N, except col ip, by mult
       for i to l do 
         for j to n do
           if j <> ip then N[i, j] := modp1(Multiply(N[i, j], mult), p) fi;
         od;
       od;

       for i to n do	 
           if i = ip then next fi;
           mult := modp1(ConvertIn(modp1(Coeff(B[i,jp],0),p),x),p);
           # Add mult times row ip to row i of B.
           for j to m do B[i, j] := modp1(Add(B[i, j], modp1(Multiply(B[ip, j], mult), p)), p) od;
           # Add mult times row ip to row i of M
           for j to k do M[i, j] := modp1(Add(M[i, j], modp1(Multiply(M[ip, j], mult), p)), p) od;
	   # Subtract mult times col i from col ip of N
           for j to l do N[j,ip] := modp1(Subtract(N[j, ip], modp1(Multiply(N[j, i], mult), p)), p) od; 
       od;
       # Multiply row ip of B and M by x. (Use modp/Shift)
       for j to m do B[ip, j] := modp1(Shift(B[ip, j], 1), p) od;
       for j to k do M[ip, j] := modp1(Shift(M[ip, j], 1), p) od;
       # Multiply all columns, except col ip of N by x.
       for i to l do 
         for j to n do
           if j <> ip then N[i, j] := modp1(Shift(N[i, j], 1), p) fi;
         od;
       od;
       
       delta[delta_index][1] := delta[delta_index][1] + 1;
    od;

    delta := sort(delta, proc(a,b) evalb(a[2] <= b[2]) end);	
    
    return M, N, [seq(delta[i][1], i=1..n)];

end:

#Below are notes for record.

 # for ord from 0 to d-1 do 
   
   # while there exists an entry B[i,j] that has a nonzero constant coeff do

      # Check if all constant coefficients of B are zero
      # for i to n do
      #     for j to m do
      #         if modp1(Coeff(B[delta[i][2],j],0),p)<>0 then break fi;
      #     od;
      #     if j <= m then break fi;
      # od;
      # if i < n+1 then
            
      #      ord := ord + 1;
      #      next;

      # find entry B[i,j] such that
      #  - B[i,j] has nonzero constant coefficient
      #  - delta[i] is minimal among all such entries
      #
      # use B[i,j] to zero out constant coefficients of B[1..i-1,j]
      # and B[i+1..n,j] while
      #  - recording row operations in M
      #  - recording the inverse of the row operations in N
      #
      # multiply row i of M and B by x 
      # multiply all column of N by x except for column i
      # increase deta[i] by one
   
   # od

   # divide all entries in B by x

 # od

 # return  M, N, delta
 

