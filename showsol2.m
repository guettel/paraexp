% creates the plots in Fig. 7.1 in the Paraexp paper

close all
clear all
mydefaults

ode = myode1(.01,1);

opts = odeset('AbsTol',1e-6,'RelTol',3e-6);
T = linspace(0,1,5);
[t,y] = ode15s(@(t,y) ode.A*y+ode.g(t,y),T,ode.u0,opts);
c = {'k-','b--','r-.','g:','m:'};
for j = 1:size(y,1),
    x = linspace(-1,1,102);
    u = [ 0 , y(j,:) , 0 ];
   
    plot(x,u,c{j}); axis([-1,1,0,1.5]); hold on; drawnow; shg; pause(0.05)
end;
legend('t = 0','t = 0.25','t = 0.5','t = 0.75','t = 1','Location','Best')
title('\alpha = 0.01, f = 1')
myeps('../pics/showsol1',1.4)
!epstopdf ../pics/showsol1.eps


%% 
figure
ode = myode1(.1,10);
opts = odeset('AbsTol',1e-6,'RelTol',3e-6);
T = linspace(0,1,5);
[t,y] = ode15s(@(t,y) ode.A*y+ode.g(t,y),T,ode.u0,opts);
for j = 1:size(y,1),
    x = linspace(-1,1,102);
    u = [ 0 , y(j,:) , 0 ];   
    plot(x,u,c{j}); axis([-1,1,0,1.5]); hold on; drawnow; shg; pause(0.05)
end
title('\alpha = 0.1, f = 10')
legend('t = 0','t = 0.25','t = 0.5','t = 0.75','t = 1','Location','Best')
myeps('../pics/showsol2',1.4)
!epstopdf ../pics/showsol2.eps


return

%%
ode = myode2(1,3);

opts = odeset('AbsTol',1e-4,'RelTol',3e-4);
T = linspace(0,1,5);
[t,y] = ode45(@(t,y) ode.A*y+ode.g(t,y),T,ode.u0,opts);
for j = 1:size(y,1),
    plot(y(j,101:end)); axis([0,100,-10,10]); hold on; drawnow; shg; pause(0.05)
end;
