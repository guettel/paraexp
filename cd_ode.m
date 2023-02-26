function ode = cd_ode(Pe,f)
% 2D linear convection--diffusion with mass matrix
% and source term


flreport('off');

% COMSOL Multiphysics Model M-file
% Generated by COMSOL 3.2 (COMSOL 3.2.0.222, $Date: 2005/09/01 18:02:30 $)

flclear fem
% COMSOL version
clear vrsn
vrsn.name = 'COMSOL 3.2';
vrsn.ext = '';
vrsn.major = 0;
vrsn.build = 222;
vrsn.rcs = '$Name:  $';
vrsn.date = '$Date: 2005/09/01 18:02:30 $';
fem.version = vrsn;

% Geometry
g1=rect2(1,1,'base','corner','pos',[-1,0]);
g2=rect2(1,1,'base','corner','pos',[0,0]);
g3=geomcomp({g1,g2},'ns',{'R1','R2'},'sf','R1+R2','edge','none');
clear s
s.objs={g3};
s.name={'CO1'};
s.tags={'g3'};

fem.draw=struct('s',s);
fem.geom=geomcsg(fem);

% Initialize mesh
fem.mesh=meshinit(fem,'hmax',[0.05]);

% Application mode 1
clear appl
appl.mode.class = 'FlConvDiff';
appl.assignsuffix = '_cd';
clear prop
prop.elemdefault='Lag2';
appl.prop = prop;
clear bnd
bnd.c0 = {0,0,'1-tanh(100)','1+tanh((2*x+1)*100)'};
bnd.c0 = {0,0,['1-tanh(' num2str(Pe) ')'],['max(0,cos(2*pi*t*' num2str(f) '))^.3333 * ( 1+tanh((2*x+1.5)*' num2str(Pe) ') )']};
bnd.type = {'Nc','cont','C','C'};
bnd.ind = [3,4,3,2,1,3,3];
appl.bnd = bnd;
clear equ
equ.D = 1/Pe;
equ.v = '-2*x*(1-y^2)';
equ.u = '2*y*(1-x^2)';
equ.ind = [1,1];
appl.equ = equ;
fem.appl{1} = appl;
fem.border = 1;
fem.outform = 'general';
clear units;
units.basesystem = 'SI';
fem.units = units;

% Solution form
fem.solform = 'general';

% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem);

%%
[K,L,M,N,D] = assemble(fem,'T',0);
[Null,Compl,Range] = flnull(N); 
Ke = Null'*K*Null; 
De = Null'*D*Null;
u0 = asseminit(fem); u0 = u0.u;

% in stationary case we have
% ud = Compl*((Range'*N*Compl)\(Range'*M));
% Le = Null'*(L-K*ud); 
% but here M is time-dependent

ud = @(t) Compl*((Range'*N*Compl)\(Range'*assemble(fem,'Out',{'M'},'T',t)));
Le = @(t) Null'*(L-K*ud(t));
f = @(t,u) De\(-Ke*u + Le(t));


ode.A = -Ke;
ode.M = De;
ode.u0 = Null'*u0;
ode.g = @(t,u) Le(t);
ode.t = [0,2];
ode.fem = fem;
ode.Null = Null;
ode.ud = @(t) ud(t);
ode.plot = @(t,y) postplot(fem, 'u', ode.Null*y(:) + ode.ud(t),...
                 'tridata',{'c','cont','internal'}, ...
                 'trimap','jet(1024)');
ode.descr = 'solve Mu'' = Au + g(t,u), plot solution with plot(t,u)';