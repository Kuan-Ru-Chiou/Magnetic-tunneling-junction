%  modified double barrier  TMR resanant double barrier   2017/01/08
clc;
clear all;
close all;
%%%check ok


%double barrier FM1/oxide1/FM/oxide2/FM2 

%one FM layer


%constant and input
hbar = 1.06e-34; %plank constant
q = 1.6e-19;  %free electron charge
m0 = 9.1e-31; %free electon mass
mf =  0.73*m0        %FM effective mass
mO =  0.32*m0        %Oxside effective mass
a = 2e-10;       %MgO unit cell length ( I don't know real value , This is come from some article. Please check)
zplus = 1i*1e-15;
kT = 0.025; %ev in room temp.
% kT = 1d-9; %ev in room temp.
IE = q*q/(2*pi*hbar); % A/ev
tO = (hbar^2)/(2*mO*(a^2)*q); %Oxside hopping energy (eV)  (oxide a and Fe a using the same value because i dont know real value)
tf = (hbar^2)/(2*mf*(a^2)*q); %FM hopping energy (eV)
exchange =2.15;  %exchange splitting energy
Ef = 2.25;  %fermi-energy
Ec = 0;  %conduction band minimum
NL =2;  NOx1 =2;  NOx2 =2;  NF = 4; NR =2;
Np = NL+NOx1+NOx2+NR+NF;
Ub = Ef+0.93;    %oxide barrier ev
% Ub = Ef+5 % it will see the resanance peak in conductance function
DD =0;    %dephasing strength
%pauli matrix
sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];


%rotation matrix
% theta = pi/2; 
%theta = pi/4;
% R = [cos(0.5*theta) -sin(0.5*theta); sin(0.5*theta) cos(0.5*theta)];



%construct hamiltonian

alphaL = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*sz*exchange;
alphaL = kron(diag([ones(1,NL) zeros(1,NOx1+NR+NOx2+NF)]),alphaL);


alphaox = [2*tO 0;0 2*tO];
alphaox = kron(diag([zeros(1,NL) ones(1,NOx1) zeros(1,NF) ones(1,NOx2) zeros(1,NR)]),alphaox);


alphaF = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*sz*exchange;
alphaF = kron(diag([zeros(1,NL+NOx1) ones(1,NF) zeros(1,NR+NOx2)]),alphaF);


alphaR = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*sz*exchange;
alphaR = kron(diag([zeros(1,NOx1+NOx2+NL+NF) ones(1,NR)]),alphaR);

alphaAR = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange+0.5*sz*exchange;
alphaAR = kron(diag([zeros(1,NOx1+NOx2+NL+NF) ones(1,NR)]),alphaAR);


% beta = [ones(1,NL)*(-tf) ones(1,NOx1-1)*(-tO) -tf ones(1,NF)*(-tf) ones(1,NOx2-1)*(-tO)  ones(1,NR)*(-tf)];
beta = [ones(1,NL-1)*(-tf) (-tf) ones(1,NOx1-1)*(-tO) -tf ones(1,NF-1)*(-tf) (-tf) ones(1,NOx2-1)*(-tO) (-tf) ones(1,NR-1)*(-tf)];
beta = kron(diag(beta,1),eye(2));

% potential barrier
UBB = [zeros(1,NL) Ub*ones(1,NOx1) zeros(1,NF) Ub*ones(1,NOx2) zeros(1,NR)];
UB = kron(diag(UBB),eye(2));


H = zeros(2*Np,2*Np);
HA =H;

H = H+alphaL+alphaR+alphaox+alphaF+beta+beta';
HA = HA + alphaL+alphaAR+alphaox+alphaF+beta+beta';

 %energy grid
%   E = linspace(Ef-14*kT,Ef+5*kT,100);
% 
%   dE = E(2)-E(1);
%   NE = length(E); 
%bias


%gate voltage adding in the NF layer
% Ug = -0.1;
Ug=0;
Ugg = [zeros(1,NL) Ug*ones(1,NOx1) Ug*ones(1,NF) Ug*ones(1,NOx2)  zeros(1,NR)]; 





V = linspace(-0.5,0.5,8);
V(5:8) = -1*V(4:-1:1);
Nv = length(V);



% sigL = zeros(2*Np,2*Np);
% sigR = sigL;
% sigAR = sigL;
%loop over every energy
ii =1;  
for i = 1:Nv%bias
    
   mu1 = Ef+0.5*V(i);  mu2 = Ef-0.5*V(i);
   
   I = 0; 
   Ipup =0; 
   Ipdn=0;
   Ia = 0; 
   Iapup =0; 
   Iapdn =0;
   
   TT=0;
   TTA = 0;
   
   TTpup =0;
   TTpdn = 0;
   TTapup = 0;
   TTapdn =0;
       
       
       
  
   
   IL = 0; %check charge conservation for me
   IR = 0; %check charge conservation for me
   
   E = linspace(min(mu1,mu2)-20*kT,max(mu1,mu2)+20*kT,300); %400 point converge
   dE = E(2)-E(1);
   NE = length(E);
    
    
    
   
    
    U = [0.5*V(i)*ones(1,NL) V(i)*linspace(0.5,-0.5,NOx1+NOx2+NF) -0.5*V(i)*ones(1,NR)]+UBB; %no gate voltage
    
    
    U = U-[zeros(1,NL+NOx1) U(NL+NOx1+1:NL+NOx1+NF) zeros(1,NOx2+NR)];    %no gate voltage
    U = U-[zeros(1,NL+NOx1) U(NL+NOx1+1:NL+NOx1+NF) zeros(1,NOx2+NR)]+Ugg;    %with gate voltage
      
    UUU =U; %for device potential plot
    
    UUUU(i,:) =UUU;
      
      
      
      
      
      U = kron(diag(U),eye(2));
    
   
   
   A0 = 0;
  
   Np0 = zeros(2*Np,1); Na0 = zeros(2*Np,1); Ap0 = zeros(2*Np,1); Aa0 = zeros(2*Np,1);
   for j = 1:NE
     
       sigL = zeros(2*Np,2*Np);   sigLup = sigL;  sigLdn = sigL;
       sigR = sigL;               sigRup = sigL;  sigRdn = sigL;
       sigAR = sigL;              sigARup = sigL;   sigARdn = sigL;
       
       
      sigB =zeros(2*Np); sigBin = sigB; %initialize guess value for dephasing  
      sigBA = sigB;     sigBAin = sigB;
       
       
% fL = (2*mf*kT*q)/(2*pi*(hbar^2))*(log(1+exp((mu1-E(j))/kT))); %2D-transverse degenacy (unit 1/(meter*meter))
% fR = (2*mf*kT*q)/(2*pi*(hbar^2))*(log(1+exp((mu2-E(j))/kT))); %2D-transverse degenacy (unit 1/(meter*meter))
   
    fL = 1/(1+exp((E(j)-mu1)/kT));    %fermi contact 1  
    fR = 1/(1+exp((E(j)-mu2)/kT));    %fermi contact 2
    
   
    
    
      FFL(i,j)=fL;  %every bias fermi function
      FFR(i,j)= fR; %every bias fermi function
     
%      U = [0.5*V(i)*ones(1,NL) V(i)*linspace(-0.5,0.5,NOx) -0.5*V(i)*ones(1,NR)] + UBB;
%       U = [0.5*V(i)*ones(1,NL) V(i)*linspace(0.5,-0.5,NOx1+NOx2+NF) -0.5*V(i)*ones(1,NR)]+UBB;
%        % U = U-[zeros(1,NL+NOx1) U(NL+NOx1+NF) zeros(1,NOx2+NR)]; 
%       U = U-[zeros(1,NL+NOx1) U(NL+NOx1+1:NL+NOx1+NF) zeros(1,NOx2+NR)];
%       UUU =U; %for device potential plot
%       U = kron(diag(U),eye(2));
%     
     
     
     
  %self -energy resolved up and dn channel
  
  sigLup(1,1) = selfdatta1D(E(j),0.5*V(i),tf);
  sigLdn(2,2) = selfdatta1D(E(j),0.5*V(i)+exchange,tf);
  sigRup(2*Np-1,2*Np-1) = selfdatta1D(E(j),-0.5*V(i),tf);
  sigRdn(2*Np,2*Np) = selfdatta1D(E(j),-0.5*V(i)+exchange,tf);
  sigARdn(2*Np-1,2*Np-1) = selfdatta1D(E(j),-0.5*V(i)+exchange,tf);
  sigARup(2*Np,2*Np) = selfdatta1D(E(j),-0.5*V(i),tf);
  %%%%%%%%%%%%%%%%%%%%%
  
  sigL(1,1) = selfdatta1D(E(j),0.5*V(i),tf);
  sigL(2,2) = selfdatta1D(E(j),0.5*V(i)+exchange,tf);
  sigR(2*Np-1:2*Np,2*Np-1:2*Np) = [selfdatta1D(E(j),-0.5*V(i),tf) 0;0 selfdatta1D(E(j),-0.5*V(i)+exchange,tf)];
  sigAR(2*Np-1:2*Np,2*Np-1:2*Np) = [selfdatta1D(E(j),-0.5*V(i)+exchange,tf) 0; 0 selfdatta1D(E(j),-0.5*V(i),tf)];
  
  
  gamLup =  1i*(sigLup - sigLup');  gamLdn =  1i*(sigLdn - sigLdn');
  gamRup = 1i*(sigRup - sigRup');   gamRdn = 1i*(sigRdn - sigRdn');
  gamARdn = 1i*(sigARdn - sigARdn'); gamARup = 1i*(sigARup - sigARup');
  
  
  
  gamL = 1i*(sigL - sigL');
  gamR = 1i*(sigR - sigR');
  gamAR = 1i*(sigAR - sigAR');
  gamLin = fL*gamL;
  gamRin = fR*gamR;
  gamARin = fR*gamAR;
  
    
     
            change = 100;
     while change > 1e-8
         g = inv((E(j))*eye(2*Np)-H-U-sigL-sigR-sigB);
        sigBnew = diag(diag(DD*g));
         % sigBnew = DD*g%phase relaxation
        change = sum(sum(abs(sigBnew-sigB)))/sum(sum(abs(sigBnew)+abs(sigB)));
%          change = sum(sum(abs(sigBnew-sigB)))
        sigB = sigB+0.5*(sigBnew-sigB);
     end
     
     
       change = 100;
     while change > 1e-8
         gn = g*(gamLin+gamRin+sigBin)*g';
      sigBinnew = diag(diag(DD*gn));
       % sigBinnew = DD*gn; %phase relaxation  
        change = sum(sum(abs(sigBinnew-sigBin)))/sum(sum(abs(sigBinnew)+abs(sigBin)));
%         change = sum(sum(abs(sigBinnew-sigBin)));
        sigBin = sigBin+0.5*(sigBinnew-sigBin);
     end
     
     
     
     
            change = 100;
     while change > 1e-8
         gA = inv((E(j))*eye(2*Np)-HA-U-sigL-sigAR-sigBA);
        sigBAnew = diag(diag(DD*g));
         %sigBAnew = DD*gA%phase relaxation
        change = sum(sum(abs(sigBAnew-sigBA)))/sum(sum(abs(sigBAnew)+abs(sigBA)));
        sigBA = sigBA+0.5*(sigBAnew-sigBA);
     end
     
     
       change = 100;
     while change > 1e-8
         gnA = gA*(gamLin+gamARin+sigBAin)*gA';
      sigBAinnew = diag(diag(DD*gnA));
     % sigBAinnew = DD*gnA; %phase relaxation  
        change = sum(sum(abs(sigBAinnew-sigBAin)))/sum(sum(abs(sigBAinnew)+abs(sigBAin)))
        sigBAin = sigBAin+0.5*(sigBAinnew-sigBAin);
     end
    
      Gless =  1i*gn;              %%  antihermitian
      GAless =  1i*gnA;            %%  antihermitian
      
      

      
     %gn = g*(gamLin+gamRin)*g';
     %A = diag(1i*(g-g'));
      A = 1i*(g-g');
     AA = 1i*(gA - gA');
     
     
     %density(j)=real(trace(gn)); %electron density
     Tp(j)= real(trace(gamL*g*gamR*g'));  %parallel magnet
     Tap(j)= real(trace(gamL*gA*gamAR*gA'));  %parallel magnet
     
     
     Tpup(j) = real(trace(gamLup*g*gamRup*g'));
     Tpdn(j) = real(trace(gamLdn*g*gamRdn*g'));
     Tapup(j) = real(trace(gamLup*gA*gamARdn*gA'));
     Tapdn(j) = real(trace(gamLdn*gA*gamARup*gA'));
     
%      I = I + Tp(j)*(fL-fR)*dE;  %paralaell total current
%      Ia = Ia +Tap(j)*(fL-fR)*dE;  %paralaell total current
Ipup = Ipup + Tpup(j)*(fL-fR)*dE;
Ipdn = Ipdn +  Tpdn(j)*(fL-fR)*dE;
Iapup = Iapup + Tapup(j)*(fL-fR)*dE;
Iapdn = Iapdn + Tapdn(j)*(fL-fR)*dE;


%energy resolved current 
      IIpup(i,j)=Tpup(j)*(fL-fR);
      IIpdn(i,j)=Tpdn(j)*(fL-fR);
      IIapup(i,j)=Tapup(j)*(fL-fR);
      IIapdn(i,j)=Tapdn(j)*(fL-fR);




  %transmition intergrand
     TT = TT +Tp(j);
     TTA = TTA +Tap(j);
     
    TTpup = TTpup + Tpup(j);
    TTpdn = TTpdn + Tpdn(j);
    TTapup = TTapup + Tpup(j);
    TTapdn = TTapdn + Tpdn(j);
    
     
   %%%%%%%%%%%  
     
   
   
   %%%for calculate  charge current
   Tt(i,j) = Tp(j);
   TtA(i,j) = Tap(j);
   
   Ttup(i,j) = Tpup(j);
   Ttdn(i,j) = Tpdn(j);
   TtAup(i,j) = Tapup(j);
   TtAdn(i,j) = Tapdn(j);
   %%%%%%%%
     
     
     
     conductanceTT(i,j) = Tp(j);  %different bias energy resolved conductance
     conductanceTp(i,j) = Tap(j);
     
     
     %local electron density

Np0 = Np0 + (-1i)*diag(Gless);
    
Na0 = Na0 + (-1i)*diag(GAless);

%local  density of state
Ap0 = Ap0 + diag(A);
Aa0 = Aa0 + diag(AA);

     
NNNp0(j,:) = (-1i)*diag(Gless);
NNNa0(j,:) = (-1i)*diag(GAless); 

AAAp0(j,:) = diag(A);
AAAa0(j,:) = diag(AA);



HHHH = H+U;
HHHA = HA+U;
     
    %energy resolved charged current
    for kkk  = 1:Np-1
        JJJp0(j,kkk) = -1*trace(HHHH(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)*Gless(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2) -HHHH(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2)*Gless(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk));
        JJJa0(j,kkk) = -1*trace(HHHA(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)*GAless(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2) -HHHA(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2)*GAless(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk));
    end



     
   end
%       Ipp(i) = I;
%       Iap(i) = Ia;
%       
      Ippup(i) = Ipup;
      Ippdn(i) = Ipdn;
      Iaapup(i) = Iapup;
      Iaapdn(i) = Iapdn;
      
       TTT(i) = TT; 
       TTTA(i) = TTA;
       
       
       TTTup(i) =  TTpup;
       TTTdn(i) = TTpdn;
       TTTAup(i) = TTapup;
       TTTAdn(i) = TTapdn;
       
       
    
end



IIpp = Tt.*(FFL-FFR);
IIap = TtA.*(FFL-FFR);



for i =1:Nv
Ipp(i,:) = sum(IIpp(i,:),2);
Iap(i,:) = sum(IIap(i,:),2);
end
Ipp = Ipp*dE;
Iap = Iap*dE;















%check terminal current
IIIIp1 = trace(gamLin*A + gamL*1i*Gless);
IIIIp2 = trace(gamRin*A + gamR*1i*Gless);

IIIIa1 =  trace(gamLin*AA + gamL*1i*GAless);
IIIIa2 =  trace(gamARin*AA + gamAR*1i*GAless);




%Local electron density for some bias
for i = 1:Np
    NNp0(i) =Np0(2*i-1) + Np0(2*i);
    NNa0(i) =Na0(2*i-1) + Na0(2*i);
end

figure
plot(a*[1:1:Np],NNp0);
hold on 
plot(a*[1:1:Np],NNa0)
hold off
xlabel('transport position (nm)');
ylabel('Local electron density (1/eV)');
set(gca,'Fontsize',15)

%Local density of state for some bias
for i = 1:Np
    AAp0(i) =Ap0(2*i-1) + Ap0(2*i);
    AAa0(i) =Aa0(2*i-1) + Aa0(2*i);
end
figure
plot(a*[1:1:Np],AAp0);
hold on 
plot(a*[1:1:Np],AAa0);
hold off
xlabel('transport position (m)');
ylabel('Local density of state (1/eV)');
set(gca,'Fontsize',15)


%fermi-occupancy
fp =  NNp0./ AAp0 ;
fa =  NNa0./ AAa0;

figure
plot(a*[1:1:Np],fp);
hold on 
plot(a*[1:1:Np],fa)
hold off
xlabel(' transport position (m)');
ylabel('quasi-fermilevel');
set(gca,'Fontsize',15)


%%%%%%%%%%%%%%




%%%local electron density for contourf rasonace state
for i = 1:Np
    NNNNp0(:,i) =NNNp0(:,2*i-1) + NNNp0(:,2*i);
    NNNNa0(:,i) =NNNa0(:,2*i-1) + NNNa0(:,2*i);
end



figure
 [xxx,yyy] = meshgrid(a*[1:1:Np],E);
contourf(xxx,yyy,real(NNNNp0),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)');
ylabel('energy (eV)');
title('energy resolved local electron density (1/eV)');
set(gca,'Fontsize',15)






%%%local  density of state for contourf
for i = 1:Np
    AAAAp0(:,i) =AAAp0(:,2*i-1) + AAAp0(:,2*i);
    AAAAa0(:,i) =AAAa0(:,2*i-1) + AAAa0(:,2*i);
end


figure
 [xxx,yyy] = meshgrid(a*[1:1:Np],E);
contourf(xxx,yyy,real(AAAAp0),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)');
ylabel('energy (eV)');
title('energy resolved local electron density of state (1/eV)');
set(gca,'Fontsize',15)

%total density of states (E)
 Dosp0 = sum(AAAp0,2);
 Dosa0 = sum(AAAa0,2);
 figure
 plot(E,Dosp0,'b*');
 hold on 
 plot(E,Dosa0,'r');
 hold off
 xlabel('Energy (eV)');
 ylabel('total density of states (1/eV)');
 set(gca,'Fontsize',15)
 

%total electron density (E)
LDOEp0 = sum(NNNp0,2);
LDOEa0 = sum(NNNa0,2);
 figure
 plot(E,LDOEp0,'b*');
 hold on 
 plot(E,LDOEa0,'r');
 hold off
xlabel('Energy (eV)');
 ylabel('total electron density (1/eV)');
 set(gca,'Fontsize',15)
 
 
 
 %electronic-states-occupancy
 fpp  =  LDOEp0./Dosp0;
 fap =   LDOEa0./Dosa0;
 figure
 plot(E,fpp ,'b*');
 hold on 
 plot(E,fap,'r');
 hold off
 xlabel('Energy (eV)');
 ylabel('electronic-states-occupancy');
 set(gca,'Fontsize',15)












% 
% figure
% [xx,yy]= meshgrid([1:1:Np],E);
% contourf(xx,yy,real(conductanceTT));



figure(1)
plot(V,Ipp,'-rs',V,Iap,'-bs');
xlabel('bias (eV)');
ylabel(' current (q/h)');
legend('Ipp','Iap');
set(gca,'Fontsize',15)




% figure(2)
%  rrr =diff(V)./diff(Ipp);
%  plot(V(1:Nv-1),rrr);
  

%   xlabel('voltage');
%  ylabel('dV/dI');

  

% 
% figure(4)
%  kkk =diff(V)./diff(Iap);
%  plot(V(1:Nv-1),rrr);
%  
% 
%   xlabel('voltage');
%  ylabel('dV/dI'); 

%TMR calculation
%method 1

%  TTT = sum(Tt,2);
%  TTTA = sum(TtA,2);
%  TMR = ((TTT-TTTA)./TTTA)*100;
% 


%method 2
TMR = ((Ipp-Iap)./Ipp)*100;
figure(5)
plot(V,TMR,'b--'); 
xlabel('voltage (V)');
 ylabel('TMR %'); 
title('TMR')
set(gca,'Fontsize',15)


% plot(V,TMR1,'-bo',V,TMR2,'-ro',V,TMR3,'-ko');
% xlabel('voltage (V)');
%  ylabel('TMR %'); 
%  legend('D = 0','D = 0.0005','D = 0.001');
% title('TMR')

% figure(6)
% plot(V,TTT); %TTT-TTTA/TTT is also similar
% xlabel('bias');
% ylabel('parallel total conductance'); 
% 
% 
% figure(7)
% plot(V,TTTA); %TTT-TTTA/TTT is also similar
% xlabel('bias');
% ylabel('antiparallel total conductance'); 
% 

% 
% figure(8)
% plot(V, conductanceTT);
% xlabel('bias');
% ylabel('parallel total conductance'); 
%   
% figure(9)
% plot(V, conductanceTp); 
% xlabel('bias');
% ylabel('antiparallel total conductance'); 



%transmition funtion for different bias voltage
kk = input('bias number ');
figure(10)
plot(E,Ttup(kk,:),'-bo',E,Ttdn(kk,:),'-ro',E,conductanceTT(kk,:),'-ko');
% plot(E,log10(Ttup(kk,:)),'-bo',E,log10(Ttdn(kk,:)),'-ro',E,log10(conductanceTT(kk,:)),'-ko','Linewidth',2);
legend('Tup','Tdn','Tup+Tdn');
xlabel('energy (eV)');
ylabel('Transmition function')
title('Energy-resolved parallel Transmition function'); 
set(gca,'Fontsize',15)


figure(11)
plot(E,TtAup(kk,:),'-bo',E,TtAdn(kk,:),'-ro',E,conductanceTp(kk,:),'-ko');
legend('Tup','Tdn','Tup+Tdn');
xlabel('energy (eV)');
ylabel('Transmition function')
title('Energy-resolved antiparallel Transmition function'); 
set(gca,'Fontsize',15)


%%%%%%%%%%%%%%%%%%%%%%


%spin resolved current

figure(13)
plot(V,Ippup,'-ro',V,Ippdn,'-bo',V,Ipp,'-go');
xlabel('bias(eV)');
ylabel('current (q/h)')
legend('Ippup','Ippdn','Ipp');
title('parallel current');
set(gca,'Fontsize',15)





figure(14)
plot(V,Iaapup,'-ro',V,Iaapdn,'-bo',V,Iap,'-go');
xlabel('bias(eV)');
ylabel('current (q/h)')
legend('Iaapup','Iaapdn','Iap');
title('anti-parallel current');
set(gca,'Fontsize',15)



figure(19)
plot(V,Ippup,'-rs',V,Ippdn,'-bx',V,Ipp,'-gd',V,Iaapup,'-cs',V,Iaapdn,'-mx',V,Iap,'-kd');
xlabel('bias(eV)');
ylabel('current (q/h)')
legend('Ippup','Ippdn','Ipp','Iaapup','Iaapdn','Iap');
set(gca,'Fontsize',15)


%band diagram for up dn channel
xxx =a*[1:1:Np]; xxL =xxx([1:NL]);xxR =xxx(NL+NOx1+NOx2+NF+1:Np);

figure(15)
plot(xxx,UUUU(kk,:)+[0*ones(1,NL) zeros(1,NOx1) 0*ones(1,NF) zeros(1,NOx2) 0*ones(1,NR)],'-ro',xxx,UUUU(kk,:)+[exchange*ones(1,NL) zeros(1,NOx1) exchange*ones(1,NF) zeros(1,NOx2) exchange*ones(1,NR)],'-bo');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'k--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'k--');
plot(xxx,(Ef)*ones(1,Np),'k--');
hold off
legend('up-channel','dn-channel');
xlabel('position (m)');
ylabel('Conduction band edge potential (eV)');
title('parallel-spin-dependent energy band');
text(0,Ef+0.5*V(kk),'mu1','Fontsize',15)
text(Np*a,Ef-0.5*V(kk),'mu2','Fontsize',15)
set(gca,'Fontsize',15)



figure(16)
plot(xxx,UUUU(kk,:)+[0*ones(1,NL) zeros(1,NOx1) 0*ones(1,NF) zeros(1,NOx2) exchange*ones(1,NR)],'-ro',xxx,UUUU(kk,:)+[exchange*ones(1,NL) zeros(1,NOx1) exchange*ones(1,NF) zeros(1,NOx2) 0*ones(1,NR)],'-bo');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,2),'k--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'k--');
plot(xxx,(Ef)*ones(1,Np),'k--');
hold off
legend('up-channel','dn-channel');
xlabel('position (m)');
ylabel('Conduction band edge potential (eV)');
title('antiparallel-spin-dependent energy band');
text(0,Ef+0.5*V(kk),'mu1','Fontsize',15)
text(Np*a,Ef-0.5*V(kk),'mu2','Fontsize',15)
set(gca,'Fontsize',15)



%device potentital 
figure(12)
plot([1:1:Np],UUUU(kk,:),'-ro');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'k--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'k--');
plot(xxx,(Ef)*ones(1,Np),'k--');
hold off
xlabel('position (m)');
ylabel('Conduction band edge potential (eV)');
title('Oxide barrier energy band');
text(0,Ef+0.5*V(kk),'mu1','Fontsize',15)
text(Np,Ef-0.5*V(kk),'mu2','Fontsize',15)
set(gca,'Fontsize',15)


%energy resolved current 
figure(17)
plot(E,IIpup(kk,:),'-bd',E,IIpdn(kk,:),'-rs',E,IIpup(kk,:)+IIpdn(kk,:),'-kx');
legend('IIpup','IIpdn','IIpup+IIpdn');
xlabel('energy (eV)');
ylabel('current (1/ev)(q/h)')
title('Energy-resolved parallel current'); 
set(gca,'Fontsize',15)

figure(18)
plot(E,IIapup(kk,:),'-bd',E,IIapdn(kk,:),'-rs',E,IIapup(kk,:)+IIapdn(kk,:),'-kx');
legend('IIapup','IIapdn','IIapup+IIapdn');
xlabel('energy (eV)');
ylabel('current (1/ev)(q/h)')
title('Energy-resolved antiparallel current'); 
set(gca,'Fontsize',15)

figure(20)
plot(E,IIpup(kk,:),'-bd',E,IIpdn(kk,:),'-rs',E,IIpup(kk,:)+IIpdn(kk,:),'-kx',E,IIapup(kk,:),'-md',E,IIapdn(kk,:),'-cs',E,IIapup(kk,:)+IIapdn(kk,:),'-gx');
legend('IIpup','IIpdn','IIpup+IIpdn','IIapup','IIapdn','IIapup+IIapdn');
xlabel('energy (eV)');
ylabel('current (1/ev)(q/h)')
title('Energy-resolved  current'); 
set(gca,'Fontsize',15)



hhh = a*[1:1:Np-1];
[iiii,jjjj] = meshgrid(hhh,E);


figure(21)
contourf(iiii,jjjj,real(JJJp0),50,'linestyle','none');
ylim([2.5 2.6])
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (ev)')
title('P total energy-position resolved charge current ');
set(gca,'Fontsize',15)


figure(22)
contourf(iiii,jjjj,real(JJJa0),50,'linestyle','none');
ylim([2.5 2.6])
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (ev)')
title('AP total energy-position resolved charge current ');
set(gca,'Fontsize',15)




