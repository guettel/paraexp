% compare various time stepping schemes; 
% creates the left of Figure 3.1 in the Paraexp paper

clear all
close all

mydefaults
N = 199;
D1 = -.5*(N+1)*spdiags(ones(N,1)*[-1 0 1],-1:1,N,N);
D2 = -(N+1)^2*gallery('tridiag',N);
k = 0.005;
A = D1 + k*D2;
A = 0.1*A;
v0 = randn(N,1); v0 = v0/norm(v0);
exact = expm(full(A))*v0;
I = speye(size(A));


%%
% explicit Euler
for n = 10:10:200,
    h = 1/n;
    v = v0;
    for j = 1:n,
        v = v + h*(A*v);
    end;
    err_ex1(n) = norm(v-exact);
end;
semilogy(find(err_ex1),err_ex1(err_ex1>0),'c-*','Color',[0 .5 1]);
hold on

%%
% RK4
for n = 10:2:50,
    h = 1/n;
    v = v0;
    for j = 1:n,
        
        k1 = A*v;
        k2 = A*(v+.5*h*k1);
        k3 = A*(v+.5*h*k2);
        k4 = A*(v+h*k3);
        v = v + 1/6*h*(k1+2*k2+2*k3+k4);
        
    end;
    err_ex4(n) = norm(v-exact);
end;
semilogy(4*find(err_ex4),err_ex4(err_ex4>0),'c+-','Color',[0 .5 1]);


%% Arnoldi

m = 120;
[V,H] = arnoldi(A,v0,m,1);
for n = 1:2:m,
    v = V(:,1:n)*(expm(H(1:n,1:n))*eye(n,1));
    err_pa(n) = norm(v-exact);
end;
semilogy(find(err_pa),err_pa(err_pa>0),'co-','Color',[0 .5 1]);


%% projection
f = (V'*exact);
v = 0*v0;
for n = 1:m,
    v = v + V(:,n)*f(n);
    err_pp(n) = norm(v-exact);
end;
semilogy(find(err_pp),err_pp(err_pp>0),'k-');



%%
% implicit Euler
for n = [1 , 10:10:200 ],
    h = 1/n;
    v = v0;
    for j = 1:n,
        v = (I - h*A)\v;
    end;
    err_im1(n) = norm(v-exact);
end;
semilogy(find(err_im1),err_im1(err_im1>0),'r-*');



%% rat Arnoldi

m = 60;
xi = 40;
AA = @(v) (A - xi*I)\v;
[V,H] = arnoldi(AA,v0,m,1);
for n = 1:2:m,
    v = V(:,1:n)*(expm(inv(H(1:n,1:n))+xi*eye(n))*eye(n,1));
    err_ra(n) = norm(v-exact);
end;
semilogy(find(err_ra),err_ra(err_ra>0),'ro-');
hold on

%% projection
f = (V'*exact);
v = 0*v0;
for n = 1:m,
    v = v + V(:,n)*f(n);
    err_rp(n) = norm(v-exact);
end;
semilogy(find(err_rp),err_rp(err_rp>0),'k-');




axis([0,200,eps/10,1])
xlabel('number of operations with A')
ylabel('2-norm error')
legend('expl Euler','expl RK4','poly Arnoldi (\sigma=\infty)   ','projection','impl Euler','RD Arnoldi (\sigma=40)','projection','Location','SouthEast');
grid on



figure(1)
myeps('../pics/compts199')
!epstopdf ../pics/compts199.eps



