% Convection--diffusion problem
% Requires Comsol!

clear all
close all

tab = '';
errors = [];
runs = 5;
p = 8;
global sol



for kdiffus = [100],
    for f = [3],

        ode = cd_ode(kdiffus,f);
        
        if 0,
            T = linspace(ode.t(1),ode.t(2),p+1);    % equispaced time decomposition
        else
            tol = 1e-1;
            opts = odeset('Mass',ode.M,'MStateDependence','none',...
                'Jacobian',ode.A,'AbsTol',tol,'RelTol',1e-13,'NormControl','on');
            serial = @(odefun,T,u0,ignore) ode15s(odefun,T,u0,opts);
            disp('get time partitioning')
            tic
            [t,y] = serial(@(t,u)ode.A*u+ode.g(t,u),ode.t,1*ode.u0,[]);
            toc
            %[t,y] = serial(@(t,u) ode.g(t,0*u),ode.t,0*ode.u0,[]);
            T = t(round(linspace(1,length(t),p+1))).';
        end;
        
        tol = 1e-3;
        opts = odeset('Mass',ode.M,'MStateDependence','none',...
            'Jacobian',ode.A,'AbsTol',tol,'RelTol',1e-13,'NormControl','on');
        serial = @(odefun,T,u0,ignore) ode15s(odefun,T,u0,opts);
        disp('serial integration over full interval')
        tic
            for k = 1:runs,
                [tout0,yout0] = serial(@(t,u)ode.A*u+ode.g(t,u),T,1*ode.u0,[]);
            end;
        t_ser = toc/runs;
        n_ser = length(tout0);
        
        
        % once more to get a higher accuracy "exact" solution
        opts = odeset(opts,'AbsTol',tol/10);
        serial = @(odefun,T,u0,ignore) ode15s(odefun,T,u0,opts);
        [toute,youte] = serial(@(t,u)ode.A*u+ode.g(t,u),T,1*ode.u0,[]);
        
       
        disp('serial integration over partitions')  % with initial value zero
        opts = odeset(opts,'AbsTol',tol/sqrt(p));
        serial = @(odefun,T,u0,ignore) ode15s(odefun,T,u0,opts);
        yout1 = zeros(p+1,length(ode.u0));
        for j = 1:p, 
            tic
            for k = 1:runs,
                [tout1,y] = serial(@(t,u)ode.A*u+ode.g(t,u),[T(j),T(j+1)],0*ode.u0,[]);
            end;
            yout1(j+1,:) = y(end,:);
            t_type1(j) = toc/runs;
            n_type1(j) = length(tout1);
        end;



        disp('exponential propagation')
        t_type2 = zeros(1,p); n_type2 = zeros(1,p);
        for k = 1:runs,
            yout1(1,:) = ode.u0.';
            yout2 = yout1;
            for j = 1:p, 
                tic
                y = yout1(j,:).';
                for s = j:p,
                    sol(1).init = 0;
                    %y = expm((T(s+1)-T(s))*ode.A)*y; m = 0;
                    
                    [y,m] = exparnoldi(ode.A,ode.M,y,T(s+1)-T(s));
                    %[y,m] = polycheby((T(s+1)-T(s))*ode.A,y,-1i*spec,1i*spec,1e-5);
                    %[y,m] = rcexpmv((T(s+1)-T(s))*ode.A,y,1e-4,@(M,v)lusolver(M,v,1));
                    %[y,m] = siexpmv((T(s+1)-T(s))*ode.A,y,1e-5,@(M,v)lusolver(M,v,1));
                    yout2(s+1,:) = yout2(s+1,:) + y.';
                end;
                t_type2(j) = t_type2(j) + toc;
                n_type2(j) = n_type2(j) + m;
            end;
        end;
        t_type2 = t_type2/runs;
        n_type2 = n_type2/runs;


        
        %errors = [ errors ; max(max(abs(yout0 - youte))) , max(max(abs(yout2 - youte))) ]
        err = [ max(max(abs(yout0 - youte))) , max(max(abs(yout2 - youte))) ];
        
        speedup = t_ser / max(t_type1 + sum(t_type2)/p);
        efficiency = speedup/p*100;
        v = [ kdiffus , f , t_ser, err(1) , max(t_type1) ,  max(t_type2) , err(2) , efficiency];
        tab = [ tab sprintf('%5.4g & %5.4g & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e & %8.0f \\%% \\\\ \n',v) ]
        
        [ t_ser , err(1) , min(t_type1) , max(t_type1) , sum(t_type2)/p , err(2) , efficiency ]
    end;
end;

