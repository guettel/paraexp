function y = lusolver(M,v,solnr)
% Solves system M*v = y, where the coefficient matrix
% is identified by the solnr.
% This function requires a global structure called sol,
% which is initialized as sol(solnr).init = 0.

    global sol
    if issparse(M),
       if ~sol(solnr).init,
           %disp(['sparse factorization nr ' num2str(solnr)])
           [sol(solnr).L,sol(solnr).U,sol(solnr).P,sol(solnr).Q,sol(solnr).R] = lu(M);
       end;
       y = sol(solnr).Q * (sol(solnr).U \ (sol(solnr).L \ (sol(solnr).P * (sol(solnr).R \ v))));      
    else
       if ~sol(solnr).init,
           disp(['sparse factorization nr ' num2str(solnr)])
           [sol(solnr).L,sol(solnr).U,sol(solnr).p] = lu(M,'vector');
       end;
       y = sol(solnr).U \ (sol(solnr).L \ v(sol(solnr).p));
    end;
    sol(solnr).init = sol(solnr).init + 1;
    
end
