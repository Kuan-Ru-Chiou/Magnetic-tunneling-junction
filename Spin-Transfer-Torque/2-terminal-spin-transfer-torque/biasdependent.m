clc;
clear all;
close all;
%non---equibrium  2017/11/19   check Torque-parallel V -; T-perpendicular not V^2 ;
%%%check ok
%%voltage asymmetry of spin transfer        contact1 is Z up and
%down.   contact 2 is x up and down

%constant and input
hbar = 1.06e-34; %plank constant
q = 1.6e-19;  %free electron charge
m0 = 9.1e-31; %free electon mass
mf =  0.73*m0        %FM effective mass
mO =  0.32*m0        %Oxside effective mass
a = 2e-10;       %lattice constant
zplus = 1i*1e-12;
kT = 1d-9; %ev in room temp.
IE = q*q/(2*pi*hbar); % A/ev
tO = (hbar^2)/(2*mO*(a^2)*q); %Oxside hopping energy (eV)  (oxide a and Fe a using the same value because i dont know real value)
tf = (hbar^2)/(2*mf*(a^2)*q); %FM hopping energy (eV)
exchange =2.15;  %exchange splitting energy 2.15
Ef = 2.25;  %fermi-energy 2.25
Ec = 0;  %conduction band minimum
NL =2;  NOx =4;  NR =2;
Np = NL+NOx+NR;
Ub = Ef+0.93        %oxide barrier ev
%pauli matrix
sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];


%rotation matrix
theta =pi/6; phi=0;
R = [cos(0.5*theta) -sin(0.5*theta); sin(0.5*theta) cos(0.5*theta)];



%construct hamiltonian

alphaL = [2*tf 0;0 2*tf] + (0.5*eye(2)*exchange-0.5*sz*exchange);
alphaL = kron(diag([ones(1,NL) zeros(1,NOx+NR)]),alphaL);



alphaox = [2*tO 0;0 2*tO];
alphaox = kron(diag([zeros(1,NL) ones(1,NOx) zeros(1,NR)]),alphaox);


alphaR = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*(sx*sin(theta)+sz*cos(theta))*exchange;
alphaR = kron(diag([zeros(1,NOx+NL) ones(1,NR)]),alphaR);



% beta = [ones(1,NL)*(-tf) ones(1,NOx-1)*(-tO) ones(1,NR)*(-tf)];
beta = [ones(1,NL-1)*(-tf) (-tf)  ones(1,NOx-1)*(-tO) (-tf) ones(1,NR-1)*(-tf)];

beta = kron(diag(beta,1),eye(2));

%potential barrier
UBB = [zeros(1,NL) Ub*ones(1,NOx) zeros(1,NR)];
UB = kron(diag(UBB),eye(2));

H = zeros(2*Np,2*Np);
H = H+alphaL+alphaR+alphaox+beta+beta';
% H = H+alphaL+alphaR+alphaox+beta+beta1;
 
%bias

ev =0.5;
V = linspace(-ev,ev,11);

% V(6) =0;
Nv = length(V);



sigL = zeros(2*Np,2*Np);
sigR = sigL;
  
%loop over every energy

Jxy1vec=zeros(Nv,1);
Jyy1vec=zeros(Nv,1);
Jzy1vec=zeros(Nv,1);
Jcharge =zeros(Nv,1);
for i =1:Nv
    
   mu1 = Ef+0.5*V(i);  mu2 = Ef-0.5*V(i);

   M = [0 0 1];
   m = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
%     m = [sin(theta)*cos(phi) sin(theta)*sin(phi) 0];
%    m = [0.5*cos(phi) 0.5*sin(phi) sqrt(3)/2];
   basein = cross(m,cross(M,m));
   baseout = cross(M,m);
   
   basisin(i,:) = basein;
   basisout(i,:) = baseout;
   
 
%   a=min(mu1,mu2);
   a = -40;
  
   b=max(mu1,mu2);
  %    b =Ef
% b=min(mu1,mu2);     
myfun=@(E) myintd(E,V(i),mu1,mu2,UBB,tf,R,H,NL,NOx,NR,exchange,sx,sy,sz,Np,zplus,theta);



 myQ =quadv(myfun,a,b,1.d-9);

 
 Jxy1vec(i,1)=myQ(1,1);
 Jyy1vec(i,1)=myQ(2,1);
 Jzy1vec(i,1)=myQ(3,1);
 Jcharge(i,1)=myQ(4,1);
   
end
 


     JJ = [Jxy1vec Jyy1vec Jzy1vec]; 
   
   
%    
   for i =1:Nv
       Jp(i) = dot(JJ(i,:),basisin(i,:));
%        Jp(i) = sum(JJ(i,:).*basisin(i,:));
       JO(i) = dot(JJ(i,:),basisout(i,:));
%         JO(i) = sum(JJ(i,:).*basisout(i,:));
   end
%      
%     Jp =Jp/(1-(0.75)^2)/0.5;
%     JO = JO/(1-(0.75)^2)/0.5;
%    
   


% 
% JOO =JO;




figure(1)
plot(V,Jcharge,'O-r','Linewidth',2);
legend('charge current');
xlabel('bias (volt)');
ylabel('charge current(q^2/h ampere)');

figure(3)
plot(V,JO,'O-r','Linewidth',2)
xlabel('bias (volt)')
ylabel('field like torque (q^2/h ampere)')

figure(2)
plot(V,Jp,'O-r','Linewidth',2)
xlabel('bias (volt)')
ylabel('anti-damping torque (q^2/h ampere)')






% Jdiff = JO-JOO;

  










