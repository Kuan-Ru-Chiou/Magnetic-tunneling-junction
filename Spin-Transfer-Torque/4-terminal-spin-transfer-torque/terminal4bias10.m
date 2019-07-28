clc;
clear all;
%%%non-equilibrium 2017/09/27
%%spin torque phi 0~2pi    theta = pi/6  non-equilibrium

%constant and input
hbar = 1.06e-34; %plank constant
q = 1.6e-19;  %free electron charge
m0 = 9.1e-31; %free electon mass
mf =  0.73*m0        %FM effective mass
mO =  0.32*m0        %Oxside effective mass
a = 2e-10;       %MgO unit cell length ( I don't know real value , This is come from some article. Please check)
zplus = 1i*1e-12;
kT = 1e-9; %ev in room temp.
IE = q*q/(2*pi*hbar); % A/ev
tO = (hbar^2)/(2*mO*(a^2)*q); %Oxside hopping energy (eV)  (oxide a and Fe a using the same value because i dont know real value)
tf = (hbar^2)/(2*mf*(a^2)*q); %FM hopping energy (eV)
exchange =2.15;  %exchange splitting energy
Ef = 2.25;  %fermi-energy
Ec = 0;  %conduction band minimum
N1 =2;  NOx =4;  N2 =2; N3 =2; N4 = 2;
Np = N1+NOx+N2+N3+N4;
Ub = Ef+0.93        %oxide barrier ev

%pauli matrix
sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];

%rotation matrix angle
theta = pi/6; 
phi = 0;

%%%construct hamiltonian

alpha3 = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*sz*exchange;
alpha3 = kron(diag([zeros(1,N1+NOx+N2) ones(1,N3) zeros(1,N4)]),alpha3);

alphaox = [2*tO 0;0 2*tO];
alphaox = kron(diag([zeros(1,N1) ones(1,NOx) zeros(1,N2) zeros(1,N3+N4) ]),alphaox);

alpha1 = [2*tO 0;0 2*tO];
alpha1 = kron(diag([ones(1,N1) zeros(1,NOx+N2+N3+N4)]),alpha1);

alpha2 = [2*tO 0;0 2*tO];
alpha2 = kron(diag([zeros(1,N1+NOx) ones(1,N2) zeros(1,N3+N4)]),alpha2);

% beta = [ones(1,N1)*(-tO) ones(1,NOx-1)*(-tO) ones(1,N2)*(-tO) 0 ones(1,N3-1)*(-tf) 0 ones(1,N4-1)*(-tf)];
beta = [ones(1,N1-1)*(-tO) -tO ones(1,NOx-1)*(-tO) -tO ones(1,N2-1)*(-tO) 0 ones(1,N3-1)*(-tf) 0 ones(1,N4-1)*(-tf)];

beta = kron(diag(beta,1),eye(2)) + kron(diag(beta,-1),eye(2));

beta(2*N1+1:2*N1+2,2*(Np-3)-1:2*(Np-3)) = [-tf 0;0 -tf];
beta(2*(Np-3)-1:2*(Np-3),2*N1+1:2*N1+2) = [-tf 0;0 -tf];

beta(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(Np-1)-1:2*(Np-1)) = [-tf 0;0 -tf];
beta(2*(Np-1)-1:2*(Np-1),2*(N1+NOx-1)+1:2*(N1+NOx-1)+2) = [-tf 0;0 -tf];

%%%% potential barrier
UBB = [Ub*ones(1,N1+NOx+N2) zeros(1,N3+N4)];
UB = kron(diag(UBB),eye(2));


%initialize hamiltonian
H = zeros(2*Np,2*Np);


%bias
ev = 0.5;
V =linspace(-ev,ev,11);

Nv = length(V);


  %hamitonian depend on phi in right contact
    alpha4 = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-...
        0.5*(sx*sin(theta)*cos(phi)+sy*sin(theta)*sin(phi)+sz*cos(theta))*exchange;
    alpha4 = kron(diag([zeros(1,N1+NOx+N2+N3) ones(1,N4)]),alpha4);

     H = alpha1+alphaox+alpha2+alpha3+alpha4+beta; %hamitonian (no bias barrier)

 
   %rotation matrix
   R = [cos(theta/2)  -sin(theta/2);...
        sin(theta/2) cos(theta/2)];



for i =1:Nv%loop bias
    
 

    %magmetic moment and in out plane basis
   M = [0 0 1];
   m = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
   basein = cross(m,cross(M,m));
   baseout = cross(M,m);
   
   basisin(i,:) = basein;
   basisout(i,:) = baseout;
   
   
    
  mu1 = Ef;  mu2 = Ef; 
  mu3 = Ef+0.5*V(i) ; mu4 = Ef-0.5*V(i);
  
%    a= min(mu3,mu4);
   a =-40;
   
   b =max(mu3,mu4);
%   b =min(mu3,mu4);
 myfun=@(E) myintd1(E,V(i),mu1,mu2,mu3,mu4,UBB,tf,R,H,N1,NOx,N2,N3,N4,exchange,sx,sy,sz,Np,zplus,tO);



 myQ=quadv(myfun,a,b,1.d-10);

 
 Jxy1vec(i,1)=myQ(1,1);
 Jyy1vec(i,1)=myQ(2,1);
 Jzy1vec(i,1)=myQ(3,1);
 Jcharge(i,1)=myQ(4,1);
 

   
  
     
   end
  
      

      
    

  
     JJ = [Jxy1vec Jyy1vec Jzy1vec]; 
  
   
   
   for i = 1:Nv
       Jin(i) = dot(JJ(i,:),basisin(i,:))
       Jout(i) = dot(JJ(i,:),basisout(i,:));
   end

figure(1)
plot(V,Jcharge,'O-r','Linewidth',2);
legend('charge current');
xlabel('bias (volt)');
ylabel('charge current(q/h ampere)');

figure(3)
plot(V,Jout,'O-r','Linewidth',2)
xlabel('bias (volt)')
ylabel('field like torque (q/h ampere)')

figure(2)
plot(V,Jin,'O-r','Linewidth',2)
xlabel('bias (volt)')
ylabel('anti-damping torque (q/h ampere)')    
     
   
   

