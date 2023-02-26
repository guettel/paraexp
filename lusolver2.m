function y = lusolver2(M,v,res)
% Solves system M*v = y, where the coefficient matrix
% is identified by its (1,1) element.
% To reset call lusolver2 with empty argument list.

    persistent sol  indm
    
    verbose = 0;
    
    if nargin < 2,
        sol = []; indm = [];
        if verbose, disp('lusolver reset'); end;
        return
    end;
    
    indm(1) = NaN; %dummy element
    m = M(1,1);
    ind = find(indm == m,1,'first');
    
    if issparse(M),
       if isempty(ind),
           ind = length(indm) + 1;
           if verbose, disp(['sparse factorization nr ' num2str(ind)]); end;
           [sol(ind).L,sol(ind).U,sol(ind).P,sol(ind).Q,sol(ind).R] = lu(M);
           indm(ind) = m;
       end;
       y = sol(ind).Q * (sol(ind).U \ (sol(ind).L \ (sol(ind).P * (sol(ind).R \ v))));      
    else
       if isempty(ind),
           ind = length(indm) + 1;
           if verbose, disp(['full factorization nr ' num2str(ind)]); end;
           [sol(ind).L,sol(ind).U,sol(ind).p] = lu(M,'vector');
           indm(ind) = m;
       end;
       y = sol(ind).U \ (sol(ind).L \ v(sol(ind).p));
    end;
    
end

