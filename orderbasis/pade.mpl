read "special.mpl";
read "popov.mpl";
read "recur.mpl";
read "../maple_examples/modpops.mpl";

#
# Input: S[1..n]   - modp1 polynomials in terms of x mod p
#        N[1..n+1] - degree bounds
# Output: (lambda, delta)
#        lambda - a solution basis
#        for any i,j, deg(lambda[i]) < N[1] and
#          deg(lambda[i] S[j] mod x^d) < N[j+1]
#
#        delta - (-N)-degrees of the rows of completion matrix
#
SimPade := proc(S,x,n,d,p,N)
    local B, M, M_adj, delta, delta2, N_sum, delta_sum, lambda, i, j,zipped;

    B := Matrix(n+1, 1);
    B[1,1] := modp1(ConvertIn(1, x), p);
    
    for i to n do
        B[i+1,1] := S[i];
    od;

    M,M_adj,delta2 := MinApproxBasis_2(B,x,n+1,1,1,1,d,N,p);

    N_sum := sum(N[j],j=1..n+1);
    delta_sum := sum(delta2[j],j=1..n+1);
 
    zipped := select(q -> q[2] < 0, 
                    zip((a,b)->[a,delta_sum - N_sum - b], convert(M_adj,list), delta2));

    return zipped[1..-1,1], zipped[1..-1,2];
end;


#
# The same input and output as SimPade
# Also does a testing for the clain in section 4.1
#
SimPade_2 := proc(S,x,n,d,p,N)
    local i,j,JB, NJ, M, delta1, delta2, G, N_sum, delta_sum, delta, G_t, J_adj_G_t_J, lambda,result,
       zipped;
    JB := Matrix(n+1,1);
    JB[n+1,1] := modp1(Constant(1,x),p);
    
    for i to n do
        JB[i,1] := S[n+1-i];
    od;
    
    NJ := [ seq(N[n+2-i], i=1..n+1) ];
   
    M, delta1 := MinApproxBasis(JB,x,n+1,1,d,NJ,p);
    G, delta2 := popov(M,x,n+1,n+1,NJ,p);

    N_sum := sum(N[j],j=1..n+1);
    delta_sum := sum(delta2[j],j=1..n+1);

    delta := [ seq(delta_sum - N_sum - delta2[n+1-i], i=0..n) ];

    G_t := map(y -> modp1(ConvertOut(y,x),p),LinearAlgebra:-Transpose(G));
    J_adj_G_t_J := ArrayTools:-FlipDimension(ArrayTools:-FlipDimension(
	map(y -> modp1(ConvertIn(y, x), p), LinearAlgebra:-Adjoint(G_t)),1),2);

    lambda := [ seq(J_adj_G_t_J[i,1], i=1..n+1) ];

    # test:
    # Dx^d + rem(lambda B^t, x^d) = J adj(G^t) J = result

    result := Matrix(n+1,n+1);
    for i to n+1 do
      for j to n+1 do
        result[i,j] := modp1(Rem(Multiply(lambda[i], JB[n+2-j,1]),ConvertIn(x^d,x)),p);
      od;
    od;
    
    #print("full delta in SimPade_2 = ", delta);
    for i to n+1 do
      if delta[i] = d - N[i] then result[i,i] := modp1(Add(result[i,i],ConvertIn(x^d,x)),p); fi;
    od;

    print("test: result = Dx^d + rem(lambda B^t, x^d) = ", result);
    print("J adj(G^t) J = ", J_adj_G_t_J);

    if submatrix(result, J_adj_G_t_J, n+1,n+1,n+1,n+1) 
      then print("test succeeds."); 
    else print("test fails.") 
    fi;
    
    zipped := select(q -> q[2] < 0, 
               zip((a,b)->[a,b], lambda, delta));

    return zipped[1..-1,1], zipped[1..-1,2];
end;

    

#
# Combine the solutions of 2 subproblems.
# N_0 is the degree bound for all lambdas
#
Intersect := proc(lambda_1, delta_1, lambda_2, delta_2, N_0, x, p)
    local k_1,k_2,k_0,s_0,d,H,M,M_adj,delta,zipped,i;

    k_1 := nops(delta_1);
    k_2 := nops(delta_2);
    k_0 := k_1 + k_2 + 1;

    s_0 := [-N_0, op(delta_1),op(delta_2)];
    d := N_0 - min(op(delta_1),op(delta_2));

    H := Matrix(k_0, 2);
    H[1..-1,1..-1] := modp1(ConvertIn(0, x),p);
    H[1,1..2] := modp1(ConvertIn(1,x),p);  
    
    for i to k_1 do
      H[i+1,1] := modp1(Multiply(Constant(-1,x), lambda_1[i]),p);
    od;

    for i to k_2 do
      H[k_1+i+1,2] := modp1(Multiply(Constant(-1,x),lambda_2[i]),p);
    od;

    M,M_adj,delta := MinApproxBasis_2(H,x,k_0,2,1,1,d,s_0,p);

    zipped := select(q->q[2]<0, zip((a,b)->[a,b],convert(M,list),delta));

    return zipped[1..-1,1], zipped[1..-1,2];

end;



#
# SimPade for general case.
# The difference is that g[1..n] is a sequence of modp1 polynomials instead of x^d. 
#  
# for any i,j, deg(lambda[i]) < N[1] and
#   deg(lambda[i] S[j] mod g[j]) < N[j+1]
#

RecursiveSimPade := proc(S,x,n,g,p,N)
    local k, lambda_1, lambda_2, delta_1, delta_2, B, M;

    if n = 1 then 
      B := Matrix(2,2);
      B[1,1] := modp1(ConvertIn(1,x),p);
      B[1,2] := S[1];
      B[2,1] := modp1(Zero(x), p);
      B[2,2] := g[1];
      # not sure how to do this
      M, delta_1 := popov(B,x,2,2,[-N[1], -N[2]],p);
      lambda_1 := [ M[1,1],M[2,1] ];
      return lambda_1, delta_1;
    fi;

    k := floor(n/2);

    lambda_1, delta_1 := RecursiveSimPade(S[1..k],x,k,g[1..k],p,N[1..k+1]);
    lambda_2, delta_2 := RecursiveSimPade(S[k+1..n],x,n-k,g[k+1..n],p,[N[1],op(N[k+2..n+1])]);

    return Intersect(lambda_1, delta_1, lambda_2, delta_2, N[1], x, p);

end;




