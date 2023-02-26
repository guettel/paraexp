% wave equation with oscillating sound source

clear all
close all

tab = '';
errors = [];
runs = 10;
p = 8;
q = 4;
global sol

serial = @rk4;

for kdiffus = [.1 , 1 , 10 ],
    for f = [1 , 5 , 25],

        ode = myode2(kdiffus,f);
        T = linspace(ode.t(1),ode.t(2),p+1);    % time decomposition
        h_ser = min(0.0005/sqrt(kdiffus),.0015/f);
        disp('serial integration over full interval')
        tic
            for k = 1:runs,
                [tout0,yout0,n_ser] = serial(@(t,u)ode.A*u+ode.g(t,u),T,1*ode.u0,h_ser);
            end;
        t_ser = toc/runs;
        
        
        % once more to get a higher accuracy "exact" solution
        [toute,youte] = serial(@(t,u)ode.A*u+ode.g(t,u),T,1*ode.u0,h_ser/5);
        
       
        disp('serial integration over partitions')  % with initial value zero
        h_ser2 = h_ser/sqrt(p)^(1/q);
        yout1 = zeros(p+1,length(ode.u0));
        for j = 1:p, 
            tic
            for k = 1:runs,
                [tout1,y,n] = serial(@(t,u)ode.A*u+ode.g(t,u),[T(j),T(j+1)],0*ode.u0,h_ser2);
            end;
            yout1(j+1,:) = y(end,:);
            t_type1(j) = toc/runs;
            n_type1(j) = n;
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
                    %y = expm((T(j+1)-T(j))*ode.A)*y; m = 0;
                    
                    spec = (T(j+1)-T(j))*sqrt(kdiffus)*202;
                    %[spec , max(imag(eig(full((T(j+1)-T(j))*ode.A))))]
                    
                    [y,m] = polycheby((T(j+1)-T(j))*ode.A,y,-1i*spec,1i*spec,1e-5);
                    %[y,m] = rcexpmv((T(j+1)-T(j))*ode.A,y,1e-4,@(M,v)lusolver(M,v,1));
                    %[y,m] = siexpmv((T(j+1)-T(j))*ode.A,y,1e-5,@(M,v)lusolver(M,v,1));
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
        
        speedup = t_ser/max([t_type1(p)+t_type2(1) , t_type1(1:p-1)+t_type2(2:p)]);
        efficiency = speedup/p*100;
        v = [ kdiffus , f , t_ser, err(1) , max(t_type1) ,  max(t_type2) , err(2) , efficiency];
        tab = [ tab sprintf('%5.4g & %5.4g & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e & %8.0f \\%% \\\\ \n',v) ]

        plot(youte(:,101:end).')
        drawnow
        shg
    end;
end;

