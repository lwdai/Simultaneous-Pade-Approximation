read "../maple_examples/modpops.mpl";

#
# Input: A - an m by n matrix
# Output: (A_rref, r, rp)
#         A_rref - the row reduced echelon form of A mod p
#         r - the rank of A
#         rp - the rank profile
#
RREF := proc(A, m, n, p)
    local B, U, Q, rp, r, i, temp, A_rref, de; 
    B := LinearAlgebra:-Copy(A); # m by n
    Q, rp, de := LinearAlgebra[Modular]:-RowEchelonTransform(p, B, true, true, true, false);
         
    r := nops(rp);

    if type(A[1,1], float[8]) then
      U := Matrix(m, m, datatype=float[8]);
    else
      U := Matrix(m, m);
    fi;

    # resize
    if n >= m then
       U[1..m,1..m] := B[1..m,1..m];
    else
       U[1..m,1..n] := B[1..m,1..n];
    fi;

    U[1..m,r+1..m] := 0;
 
    for i from r+1 to m do U[i,i] :=1 od;

    for i to ArrayTools:-Size(Q, 2) do
      if i <> Q[i] then
        temp := A[i];
        A[i] := A[Q[i]];
        A[Q[i]] := temp;
      fi;
    od;
       
    A_rref := LinearAlgebra[Modular]:-Multiply(p, U, A);

    return A_rref, r, rp;
end;

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
    local recur, base_case, delta, M_all, A_0;
    # Initialize list delta = [ [s[1],1], [s[2],2],...,s[n],n] ] of shifts.
    delta := [ seq([s[i],i], i=1..n) ];

    if n*(p-1)^2 < 2^53-1 then # not sure how this condition is built
       A_0 := Matrix(n, m, datatype=float[8]);
    else
       A_0 := Matrix(n, m);
    fi;


# returns M such that M A = 0 mod x^d
recur := proc(A, d)
    local d_1, d_2, R, R_1, R_2, M_1, M_2, i, j,
          row_sorted, A_t_rref, r, rp, T, M ;

    if d = 1 then # base case
     delta := sort(delta, (a,b)->(a[1] < b[1]));
  
     row_sorted := [ seq(delta[i][2], i=1..n) ];

     A_0[1..n,1..m] := LinearAlgebra:-Map(y -> modp1(Coeff(y, 0),p), A[row_sorted]);
     A_t_rref, r, rp:=  RREF(LinearAlgebra:-Transpose(A_0), m, n, p);

     T := LinearAlgebra:-IdentityMatrix(n, compact=false);
     for i to r do
       j := rp[i];
       T[j+1..n,j] := - A_t_rref[i, j+1...n];
     od;
     
     # applying sort operation as before to A_0, but converting to column operations
     M := Matrix(n, n);
     for j to n do   
       M[1..n, row_sorted[j]] := T[1..n, j];
     od;

     M := LinearAlgebra:-Map(y -> modp1(ConvertIn(round(y), x), p), M);
     for i to r do
       j := rp[i];
       # multiply pivot rows by x and increment the s-degree
       delta[j][1] := delta[j][1] + 1;  
       M[j] := LinearAlgebra:-Map(y -> modp1(Shift(y, 1), p), M[j]);
     od;

     for i to n do delta[i][2] := i; od; # update row number

     return M;		

    fi;

    # recursion
    d_1 := floor(d/2);
    d_2 := d - d_1;
   
    M_1 := recur(A, d_1);
    R := Matrix(n, m); 
    modpmult(M_1,A,x,n,n,m,p,R); # R = M_1 A = 0 mod x^d_1
    
    R := LinearAlgebra:-Map(y -> modp1(Shift(y,-d_1),p), R); # R = M_1 A / (x^d_1)	
    M_2 := recur(R, d_2); # M_2 M_1 A = 0 mod x^d

    M := Matrix(n, n);
    modpmult(M_2,M_1,x,n,n,n,p,M); # M_result = M_2 M_1


    return M;
      
end;

    M_all := recur(A, d);
    delta := sort(delta, proc(a,b) evalb(a[2]<b[2]) end);
    delta := map(a -> a[1],delta);

    return M_all, delta;
end;
    
