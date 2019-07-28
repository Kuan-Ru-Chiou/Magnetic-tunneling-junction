clc;
clear all;
close all;
%  2018/1/3   TMR using Transmition function converge succesful
%%%check ok
%  contact1 is Z up and contact2 is down

%constant and input
hbar = 1.06e-34; %plank constant
q = 1.6e-19;  %free electron charge
m0 = 9.1e-31; %free electon mass
mf =  0.73*m0        %FM effective mass
mO =  0.32*m0        %Oxside effective mass
a = 2e-10;       %MgO unit cell length ( I don't know real value , This is come from some article. Please check)
zplus = 1i*1e-15;
kT = 1d-9; %ev in room temp.
% IE = q*q/(2*pi*hbar); % A/ev
IE = 1; % A/ev set 1 for check error
tO = (hbar^2)/(2*mO*(a^2)*q); %Oxside hopping energy (eV)  (oxide a and Fe a using the same value because i dont know real value)
tf = (hbar^2)/(2*mf*(a^2)*q); %FM hopping energy (eV)
exchange =2.15;  %exchange splitting energy
Ef = 2.25;  %fermi-energy
Ec = 0;  %conduction band minimum
NL =2;  NOx =2;  NR =2;
Np = NL+NOx+NR;
Ub = Ef+0.93        %oxide barrier ev
DD =0;

% J2NI =0.009%0~6 ev2-nm  
%pauli matrix
sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];


%rotation matrix
% theta = pi/2; 
%theta = pi/4;
% R = [cos(0.5*theta) -sin(0.5*theta); sin(0.5*theta) cos(0.5*theta)];



%construct hamiltonian

alphaL = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*sz*exchange;
alphaL = kron(diag([ones(1,NL) zeros(1,NOx+NR)]),alphaL);


alphaox = [2*tO 0;0 2*tO];
alphaox = kron(diag([zeros(1,NL) ones(1,NOx) zeros(1,NR)]),alphaox);


alphaR = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange-0.5*sz*exchange;
alphaR = kron(diag([zeros(1,NOx+NL) ones(1,NR)]),alphaR);

alphaAR = [2*tf 0;0 2*tf] + 0.5*eye(2)*exchange+0.5*sz*exchange;
alphaAR = kron(diag([zeros(1,NOx+NL) ones(1,NR)]),alphaAR);


beta = [ones(1,NL-1)*(-tf) (-tf) ones(1,NOx-1)*(-tO) (-tf) ones(1,NR-1)*(-tf)];
beta = kron(diag(beta,1),eye(2));

%potential barrier
UBB = [zeros(1,NL) Ub*ones(1,NOx) zeros(1,NR)];
UB = kron(diag(UBB),eye(2));


H = zeros(2*Np,2*Np);
HA =H;

H = H+alphaL+alphaR+alphaox+beta+beta';
HA = HA + alphaL+alphaAR+alphaox+beta+beta';

 %energy grid
%   E = linspace(Ef-10*kT,Ef+10*kT,100);
%   dE = E(2)-E(1);
%   NE = length(E); 
%bias

V = linspace(-0.4,0.4,8);
Nv = length(V);


sigL = zeros(2*Np,2*Np);
sigR = sigL;
sigAR = sigL;



%loop over every energy
ii =1;  
for i = 1:Nv%bias
    
   mu1 = Ef+0.5*V(i);  mu2 = Ef-0.5*V(i);
   
       I = 0;
       Ia = 0;
       Is1 =0; Is2 =0; IIS = 0; Ipup =0; Ipdn=0;
       Isa1 =0; Isa2 =0; IISa =0;  Iapup =0; Iapdn =0;
       TT=0;
       TTA = 0;
       TTs = 0;
       TTAs =0;
       
       TTpup =0;
       TTpdn = 0;
       TTapup = 0;
       TTapdn =0;
       
       
   
   
   E = linspace(min(mu1,mu2)-20*kT,max(mu1,mu2)+20*kT,1000);
   %dE = E(2)-E(1);
   dE = 1;  %set 1 for check error
   NE = length(E);
   
   
   IL = 0; %check charge conservation for me
   IR = 0; %check charge conservation for me
 
   LDOSP = 0;
   LDOSAP = 0;
   GGP =0;
   GGAP =0;
  
   for j = 1:NE
     

 
%Initialize self-energy
sigL = zeros(2*Np,2*Np);  
sigR = sigL;               
sigAR = sigL;             
sigLup = sigL;  
sigLdn = sigL;
sigRup = sigL;  
sigRdn = sigL;
sigARup = sigL;   
sigARdn = sigL;




 %spin-flip dephasing
 sigs = zeros(2*Np);
 sigAs = sigs;
 sigins = sigs;
 siginAs = sigs;
      
      
      
      
%fL = (2*mf*kT*q)/(2*pi*(hbar^2))*(log(1+exp((mu1-E(j))/kT))); %2D-transverse degenacy (unit 1/(meter*meter))
% fR = (2*mf*kT*q)/(2*pi*(hbar^2))*(log(1+exp((mu2-E(j))/kT))); %2D-transverse degenacy (unit 1/(meter*meter))
   
    fL = 1/(1+exp((E(j)-mu1)/kT));    %fermi contact 1  
    fR = 1/(1+exp((E(j)-mu2)/kT));    %fermi contact 2
    
   
    
    
      FFL(i,j)=fL;  %every bias fermi function
      FFR(i,j)= fR; %every bias fermi function
     
%      U = [0.5*V(i)*ones(1,NL) V(i)*linspace(0.5,-0.5,NOx) -0.5*V(i)*ones(1,NR)] + UBB;
     U = [0.5*V(i)*ones(1,NL) V(i)*(0.5-[1:1:NOx]/(NOx+1)) -0.5*V(i)*ones(1,NR)] + UBB;
     UUUU(i,:) =U;
     U = kron(diag(U),eye(2));
    
     
     
    
   %%%Self -energy resolved up and dn channel
   sigLup(1,1) = selfdatta1D(E(j),0.5*V(i),tf);
   sigLdn(2,2) = selfdatta1D(E(j),0.5*V(i)+exchange,tf);
   sigRup(2*Np-1,2*Np-1) = selfdatta1D(E(j),-0.5*V(i),tf);
   sigRdn(2*Np,2*Np) = selfdatta1D(E(j),-0.5*V(i)+exchange,tf);
   sigARdn(2*Np-1,2*Np-1) = selfdatta1D(E(j),-0.5*V(i)+exchange,tf);
   sigARup(2*Np,2*Np) = selfdatta1D(E(j),-0.5*V(i),tf);
   
   gamLup =  1i*(sigLup - sigLup');  gamLdn =  1i*(sigLdn - sigLdn');
   gamRup = 1i*(sigRup - sigRup');   gamRdn = 1i*(sigRdn - sigRdn');
   gamARdn = 1i*(sigARdn - sigARdn'); gamARup = 1i*(sigARup - sigARup');
   
   gamLinup = fL*gamLup; gamLindn = fL*gamLdn;
   gamRinup = fR*gamRup; gamRindn = fR*gamRdn;
   gamARinup = fR*gamARup; gamARindn = fR*gamARdn;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
   %Self-energy for L and R contact
   sigL(1,1) = selfdatta1D(E(j),0.5*V(i),tf);
   sigL(2,2) = selfdatta1D(E(j),0.5*V(i)+exchange,tf);
   sigR(2*Np-1:2*Np,2*Np-1:2*Np) = [selfdatta1D(E(j),-0.5*V(i),tf) 0;0 selfdatta1D(E(j),-0.5*V(i)+exchange,tf)];
   sigAR(2*Np-1:2*Np,2*Np-1:2*Np) = [selfdatta1D(E(j),-0.5*V(i)+exchange,tf) 0; 0 selfdatta1D(E(j),-0.5*V(i),tf)];
   
   gamL = 1i*(sigL - sigL');  
   gamR = 1i*(sigR - sigR'); 
   gamAR = 1i*(sigAR - sigAR');
   gamLin = fL*gamL;
   gamRin = fR*gamR;    
   gamARin = fR*gamAR;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     

     
     change = 100;
     while change > 1e-13
         gs = inv((E(j)+zplus)*eye(2*Np)-H-U-sigL-sigR-sigs);
         ggs = diag(gs);
         gggs = zeros(2*Np,2*Np);
         for kk =1:2:2*Np-1
             gggs(kk,kk) = 2*ggs(kk+1) + ggs(kk);
             gggs(kk+1,kk+1) = 2*ggs(kk) + ggs(kk+1);
             gggs(kk,kk+1) = -1*gs(kk,kk+1);
             gggs(kk+1,kk) = -1*gs(kk+1,kk);
         end
         sigsnew = 0.25*DD*gggs;
%          change = sum(sum(abs(sigsnew-sigs)))/sum(sum(abs(sigsnew+sigs)));
         change = sum(sum(abs(sigsnew-sigs)))/sum(sum(abs(sigsnew)+abs(sigs)));
         sigs = sigs+0.5*(sigsnew-sigs);
     end
     As =1i*(gs-gs');
     LDOSP = LDOSP + As;
     energyAs(j,:)= diag(As);
     
            change = 100;
     while change > 1e-13
         gsn = gs*(gamLin+gamRin+sigins)*gs';
         ggsn = diag(gsn);
         gggsn = zeros(2*Np,2*Np);
            for kk =1:2:2*Np-1
             gggsn(kk,kk) = 2*ggsn(kk+1) + ggsn(kk);
             gggsn(kk+1,kk+1) = 2*ggsn(kk) + ggsn(kk+1);
             gggsn(kk,kk+1) = -1*gsn(kk,kk+1);
             gggsn(kk+1,kk) = -1*gsn(kk+1,kk);
            end 
         sigsinnew = 0.25*DD*gggsn;
%         change = sum(sum(abs(sigsinnew-sigins)))/sum(sum(abs(sigsinnew+sigins)))
     change = sum(sum(abs(sigsinnew-sigins)))/sum(sum(abs(sigsinnew)+abs(sigins)));
     sigins = sigins+0.5*(sigsinnew-sigins);
     end
     
     GGP = GGP + gsn;
     energygsn(j,:)= diag(gsn);
     
     
     change = 100;
     while change > 1e-13
         gA = inv((E(j)+zplus)*eye(2*Np)-HA-U-sigL-sigAR-sigAs);
         ggA = diag(gA);
         gggA = zeros(2*Np,2*Np);
         for kk =1:2:2*Np-1
             gggA(kk,kk) = 2*ggA(kk+1) + ggA(kk);
             gggA(kk+1,kk+1) = 2*ggA(kk) + ggA(kk+1);
             gggA(kk,kk+1) = -1*gA(kk,kk+1);
             gggA(kk+1,kk) = -1*gA(kk+1,kk);
         end
         
         sigAsnew = 0.25*DD*gggA;
%                  change = sum(sum(abs(sigAsnew-sigAs)))/sum(sum(abs(sigAsnew+sigAs)));
         change = sum(sum(abs(sigAsnew-sigAs)))/sum(sum(abs(sigAsnew)+abs(sigAs)));
         sigAs = sigAs+0.5*(sigAsnew-sigAs);
     end
        AAs =1i*(gA-gA');
        LDOSAP = LDOSAP + AAs;
        energyAAs(j,:)= diag(AAs);
        
        change = 100;
        while change > 1e-13
            gAn = gA*(gamLin+gamARin+siginAs)*gA';
            ggAn = diag(gAn);
            gggAn = zeros(2*Np,2*Np);
            for kk =1:2:2*Np-1
                gggAn(kk,kk) = 2*ggAn(kk+1) + ggAn(kk);
                gggAn(kk+1,kk+1) = 2*ggAn(kk) + ggAn(kk+1);
                gggAn(kk,kk+1) = -1*gAn(kk,kk+1);
                gggAn(kk+1,kk) = -1*gAn(kk+1,kk);
            end
            sigAsinnew = 0.25*DD*gggAn;
%                     change = sum(sum(abs(sigAsinnew-siginAs)))/sum(sum(abs(sigAsinnew+siginAs)));
            change = sum(sum(abs(sigAsinnew-siginAs)))/sum(sum(abs(sigAsinnew)+abs(siginAs)))
            siginAs = siginAs+0.5*(sigAsinnew-siginAs);
        end
        
        GGAP = GGAP + gAn;
     
       energygAn(j,:)= diag(gAn);
       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %use for  spin-flip dephasing
       Gless =  1i*gsn;              %%  antihermitian
       GAless =  1i*gAn;            %%  antihermitian

     
  
%For coherent transport
Tsp(j)= real(trace(gamL*gs*gamR*gs'));  %parallel magnet
Tsap(j)= real(trace(gamL*gA*gamAR*gA'));  %antiparallel magnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
   
     Tpupup(j) = real(trace(gamLup*gs*gamRup*gs'));
     Tpdndn(j) = real(trace(gamLdn*gs*gamRdn*gs'));
     Tpupdn(j) = real(trace(gamLup*gs*gamRdn*gs'));
     Tpdnup(j) = real(trace(gamLdn*gs*gamRup*gs'));
     
     Tapupup(j) = real(trace(gamLup*gA*gamARdn*gA'));
     Tapdndn(j) = real(trace(gamLdn*gA*gamARup*gA'));
     Tapupdn(j) = real(trace(gamLup*gA*gamARup*gA'));
     Tapdnup(j) = real(trace(gamLdn*gA*gamARdn*gA'));
      
     %      I = I + (IE*dE)*(Tp(j)*(fL-fR));  %paralaell total current
     %      Ia = Ia + (IE*dE)*(Tap(j)*(fL-fR));  %anti-paralaell total current
     
     
     
     %%% For non-coherent transport 
     Is1 = Is1 + (IE*dE)*(trace(gamLin*As-gamL*gsn));%paralaell total current L
     Is2 = Is2 + (IE*dE)*(trace(gamRin*As-gamR*gsn));%paralaell total current R
     Isa1 = Isa1 + (IE*dE)*(trace(gamLin*AAs-gamL*gAn));  %antiparalaell total current L
     Isa2 = Isa2 + (IE*dE)*(trace(gamARin*AAs-gamAR*gAn));%antiparalaell total current R
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
      

     %be careful this is an inelastic no energy exchange so use this
     %current formula
     Ipup = Ipup + (trace(gamLinup*As-gamLup*gsn))*(IE*dE);
     Ipdn = Ipdn + (trace(gamLindn*As-gamLdn*gsn))*(IE*dE);  
     Iapup = Iapup + (trace(gamLinup*AAs-gamLup*gAn))*(IE*dE);
     Iapdn = Iapdn + (trace(gamLindn*AAs-gamLdn*gAn))*(IE*dE);

     %energy resolved current at source
      IIpup(i,j)=trace(gamLinup*As-gamLup*gsn);
      IIpdn(i,j)=trace(gamLindn*As-gamLdn*gsn);
      IIapup(i,j)=trace(gamLinup*AAs-gamLup*gAn);
      IIapdn(i,j)=trace(gamLindn*AAs-gamLdn*gAn);
      

    HHHH = H+U;
    HHHA = HA+U;
     
%      TT = TT +Tp(j);
%      TTA = TTA +Tap(j);
     
TTs  = TTs + Tsp(j);
TTAs = TTAs + Tsap(j);

Tt(i,j) = Tsp(j);
TtA(i,j)= Tsap(j);


TTpup  = TTpup + Tpupup(j);
TTpdn  = TTpdn + Tpdndn(j);
TTapup = TTapup + Tpupup(j);
TTapdn = TTapdn + Tpdndn(j);


Ttup(i,j) = Tpupup(j);
Ttdn(i,j) = Tpdndn(j);
TtAup(i,j)= Tapupup(j);
TtAdn(i,j)= Tapdndn(j);
       
       
conductanceTT(i,j) = Tsp(j);  %different bias energy resolved conductance
conductanceTp(i,j) = Tsap(j);


 
    %energy resolved charged current
    for kkk  = 1:Np-1
        JJJp0(j,kkk) = (-1*trace(HHHH(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)*Gless(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2) -HHHH(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2)*Gless(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)))*(IE*dE);
        JJJa0(j,kkk) = -(1*trace(HHHA(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)*GAless(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2) -HHHA(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2)*GAless(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)))*(IE*dE);
    end


    
    
   end
   
   
%       Ipp(i) = I;
%       Iap(i) = Ia;

Ispp1(i) = Is1;
Ispp2(i) = Is2;


Isap1(i) = Isa1;
Isap2(i) = Isa2;



Ippup(i) = Ipup;
Ippdn(i) = Ipdn;
Iaapup(i)= Iapup;
Iaapdn(i)= Iapdn;
      
      
      
TTTs(i)  = TTs;
TTTAs(i) = TTAs;

TTTup(i) =  TTpup;
TTTdn(i) = TTpdn;
TTTAup(i)= TTapup;
TTTAdn(i)= TTapdn;

       
       
end

%check in and out terminal current conservation
fprintf('charge conservation for parallel %g\n',sum(Ispp1+Ispp2));
fprintf('charge conservation for anti-parallel %g\n',sum(Isap1+Isap2));

fprintf('charge conservation for spin resolved parallel %g\n',sum(sum(((IIpup+IIpdn)*(IE*dE)),2)-Ispp1'));
fprintf('charge conservation for spin resolved anti-parallel %g\n',sum(sum(((IIapup+IIapdn)*(IE*dE)),2)-Isap1'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



xxx = a*[1:1:Np];

%total density of states up and dn 
LDOSP = diag(LDOSP);
LDOSAP = diag(LDOSAP);
Appup  = LDOSP(1:2:end);
Appdn = LDOSP(2:2:end);
Aptot = Appup+Appdn ;

Aapup = LDOSAP(1:2:end);
Aapdn =  LDOSAP(2:2:end);
Aaptot = Aapup+Aapdn ;

plot(xxx,Appup,xxx,Appdn,xxx,Aptot);

plot(xxx,Aapup,xxx,Aapdn,xxx,Aaptot);

%total electron density up and down

GGP = diag(GGP);
GGAP = diag(GGAP);

Gpup  = GGP(1:2:end);
Gpdn = GGP(2:2:end);
Gptot = Gpup + Gpdn;


Gapup = GGAP (1:2:end);
Gapdn = GGAP (2:2:end);
Gaptot = Gapup + Gapdn;


plot(xxx,Gpup,xxx,Gpdn,xxx,Gptot);

plot(xxx,Gapup,xxx,Gapdn,xxx,Gaptot);





%energy resolved position up and dn spectrum

energyAsup = energyAs(:,1:2:end);
energyAsdn = energyAs(:,2:2:end);
energyAstot = energyAsup + energyAsdn;



energyAAsup = energyAAs(:,1:2:end);
energyAAsdn = energyAAs(:,2:2:end);
energyAAstot = energyAAsup + energyAAsdn;



energygsnup = energygsn(:,1:2:end);
energygsndn = energygsn(:,2:2:end);
energygsntot = energygsnup + energygsndn;




energygAnup = energygAn(:,1:2:end);
energygAndn = energygAn(:,2:2:end);
energygAntot = energygAnup + energygAndn;



[iii,jjj] = meshgrid(xxx,E);

figure(4)
contourf(iii,jjj,real(energyAsup),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('P LDOS up(1/eV)');

figure(5)
contourf(iii,jjj,real(energyAsdn),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('P LDOS dn(1/eV)');


figure(6)
contourf(iii,jjj,real(energyAstot),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('P  LDOS(1/eV)');


figure(7)
contourf(iii,jjj,real(energyAAsup),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('AP LDOS up(1/eV)');


figure(8)
contourf(iii,jjj,real(energyAAsdn),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('AP LDOS dn(1/eV)');


figure(9)
contourf(iii,jjj,real(energyAAstot),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('AP LDOS(1/eV)');


figure(10)
contourf(iii,jjj,real(energygsnup),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('P electron up');


figure(11)
contourf(iii,jjj,real(energygsndn),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('P electron dn(1/eV)');


figure(12)
contourf(iii,jjj,real(energygsntot),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('P electron density(1/eV)');


figure(13)
contourf(iii,jjj,real(energygAnup),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('AP electron up(1/eV)');


figure(14)
contourf(iii,jjj,real(energygAndn),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('AP electron dn(1/eV)');


figure(15)
contourf(iii,jjj,real(energygAntot),100,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('AP electron density(1/eV)');



hhh = a*[1:1:Np-1];
[iiii,jjjj] = meshgrid(hhh,E);


figure(16)
contourf(iiii,jjjj,real(JJJp0),50,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('P total energy-position charge current(1/eV)');


figure(17)
contourf(iiii,jjjj,real(JJJa0),50,'linestyle','none');
shading flat;
colormap('hot');
colorbar;
xlabel('position (m)')
ylabel('energy (eV)')
title('AP total energy-position charge current(1/ev)');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% IIpp = Tt.*(FFL-FFR);
% IIap = TtA.*(FFL-FFR);
% 
% 
% 
% for i =1:Nv
% Ispp(i,:) = sum(IIpp(i,:),2);
% Isap(i,:) = sum(IIap(i,:),2);
% end
% Ispp = Ispp.*ddE;
% Isap = Isap.*ddE;


%check terminal current
IIIIp1 =    real(trace(gamLin*As + gamL*1i*Gless));
IIIIp2 =    real(trace(gamRin*As + gamR*1i*Gless));
% IIIIss =    real(trace(sigins*As + gamR*1i*Gless));   %wrong


IIIIa1 =    real(trace(gamLin*AAs + gamL*1i*GAless));
IIIIa2 =    real(trace(gamARin*AAs + gamAR*1i*GAless));
% IIIIssa =   real(trace(siginAs*AAs + gamAR*1i*GAless)); %wrong



%check local current density for Ap and P
for i  = 1:Np-1
    Jp0(i) = -1*trace(HHHH(2*i+1:2*i+2,2*i-1:2*i)*Gless(2*i-1:2*i,2*i+1:2*i+2) -HHHH(2*i-1:2*i,2*i+1:2*i+2)*Gless(2*i+1:2*i+2,2*i-1:2*i));
    Ja0(i) = -1*trace(HHHA(2*i+1:2*i+2,2*i-1:2*i)*GAless(2*i-1:2*i,2*i+1:2*i+2) -HHHA(2*i-1:2*i,2*i+1:2*i+2)*GAless(2*i+1:2*i+2,2*i-1:2*i));
end



% figure(1)
% plot(V,Ipp)
% xlabel('voltage');
% ylabel('current');
% 
% 
% 
% figure(2)
%  rrr =diff(V)./diff(Ipp);
%  plot(V(1:Nv-1),rrr);
%   ylim([0 500]);
% 
%   xlabel('voltage');
%  ylabel('dV/dI');
% 
%    figure(3)
% plot(V,Iap)
% xlabel('voltage');
% ylabel('current');
% 
% 
% 
% figure(4)
%  kkk =diff(V)./diff(Iap);
%  plot(V(1:Nv-1),rrr);
%   ylim([0 500]);
% 
%   xlabel('voltage');
%  ylabel('dV/dI'); 
% 
% 

%  TMR =(( TTTs- TTTAs)./TTTs);



TMR =(( TTTs- TTTAs)./TTTAs);
figure(3)
plot(V,TMR*100);
xlabel('voltage');
ylabel('TMR %'); 
title('TMR')
set(gca,'Fontsize',15)








 
 figure(1)
 plot(V,Ispp1,'-rs',V,Isap1,'-bx')
 legend('Ipp1','Iap1')
 xlabel('bias (V)');
 ylabel('current (q/h Ampere)');
 title('current vs. bias')
 set(gca,'Fontsize',15)
% 
% 
% figure(2)
%  rrr =diff(V)./diff(Ispp);
%  plot(V(1:Nv-1),rrr);
%   xlabel('voltage');
%  ylabel('dV/dI');

figure(2)
plot(V,-Ispp2,'-rs',V,-Isap2,'-bx');
legend('Ipp2','Iap2')
xlabel('bias (V)');
ylabel('current (q/h Ampere)');
title('current vs. bias')
set(gca,'Fontsize',15)


% figure(4)
%  kkk =diff(V)./diff(Isap);
%  plot(V(1:Nv-1),rrr);
%   xlabel('voltage');
%  ylabel('dV/dI'); 

% figure(5)
% plot(V,((Ispp-Isap)./Ispp)*100);
% xlabel('voltage');
%  ylabel('TMR %'); 


  
%band diagram for up dn channel
kk = input('bias number ');

xxx =a*[1:1:Np]; xxL =xxx([1:NL]); 
xxR =xxx(NL+NOx+1:Np);




figure(18)
plot(xxx,UUUU(kk,:)+[0*ones(1,NL) zeros(1,NOx) 0*ones(1,NR)],'-rs',xxx,UUUU(kk,:)+[exchange*ones(1,NL) zeros(1,NOx) exchange*ones(1,NR)],'-bd');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'g--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'g--');
plot(xxx,(Ef)*ones(1,Np),'g--');
hold off
legend('up-channel','dn-channel');
xlabel('position (m)');
ylabel('Conduction band edge potential (eV)');
title('parallel-spin-dependent energy band');
text(0*a,Ef+0.5*V(kk),'mu1','Fontsize',15)
text(Np*a,Ef-0.5*V(kk),'mu2','Fontsize',15)
set(gca,'Fontsize',15)

figure(19)
plot(xxx,UUUU(kk,:)+[0*ones(1,NL) zeros(1,NOx) exchange*ones(1,NR)],'-rs',xxx,UUUU(kk,:)+[exchange*ones(1,NL) zeros(1,NOx) 0*ones(1,NR)],'-bd');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'g--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'g--');
plot(xxx,(Ef)*ones(1,Np),'g--');
hold off
legend('up-channel','dn-channel');
xlabel('position (m)');
ylabel('Conduction band edge potential (eV)');
title('antiparallel-spin-dependent energy band');
text(0*a,Ef+0.5*V(kk),'mu1','Fontsize',15)
text(Np*a,Ef-0.5*V(kk),'mu2','Fontsize',15)
set(gca,'Fontsize',15)

%device potentital 
figure(20)
plot([1:1:Np],UUUU(kk,:),'-rs');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'g--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'g--');
plot(xxx,(Ef)*ones(1,Np),'g--');
hold off
xlabel('position (m)');
ylabel('potential (ev)');
title('Oxide barrier energy band');
text(0*a,Ef+0.5*V(kk),'mu1','Fontsize',15)
text(Np*a,Ef-0.5*V(kk),'mu2','Fontsize',15)
set(gca,'Fontsize',15)

%spin resolved current

figure(21)
plot(V,Ippup,'-rs',V,Ippdn,'-bx',V,Ispp1,'-gd');
xlabel('bias (eV)')
ylabel('current (q/h Ampere)')
legend('Ippup','Ippdn','Ipp');
title('parallel current');
set(gca,'Fontsize',15)


figure(22)
plot(V,Iaapup,'-rs',V,Iaapdn,'-bx',V,Isap1,'-gd');
xlabel('bias (eV)')
ylabel('current (q/h Ampere)')
legend('Iaapup','Iaapdn','Iap');
title('anti-parallel current');
set(gca,'Fontsize',15)


%transmition funtion for different bias voltage

figure(23)
plot(E,Ttup(kk,:),'-bd',E,Ttdn(kk,:),'-rs',E,conductanceTT(kk,:),'-kx');
legend('Ttup','Ttdn','Ttup+Ttdn');
xlabel('energy (eV)');
ylabel('Transmition function')
title('Energy-resolved parallel Transmition function'); 
set(gca,'Fontsize',15)


figure(24)
plot(E,TtAup(kk,:),'-bd',E,TtAdn(kk,:),'-rs',E,conductanceTp(kk,:),'-kx');
legend('TtAup','TtAdn','TtAup+TtAdn');
xlabel('energy (eV)');
ylabel('Transmition function')
title('Energy-resolved antiparallel Transmition function');
set(gca,'Fontsize',15)


figure(28)
plot(E,Ttup(kk,:),'-bd',E,Ttdn(kk,:),'-rs',E,conductanceTT(kk,:),'-kx',E,TtAup(kk,:),'-gd',E,TtAdn(kk,:),'-ms',E,conductanceTp(kk,:),'-cx');
legend('Ttup','Ttdn','Ttup+Ttdn','TtAup','TtAdn','TtAup+TtAdn');
xlabel('energy (eV)');
ylabel('Transmition function')
title('Energy-resolved Transmition function');
set(gca,'Fontsize',15)

%energy resolved current at source 
figure(25)
plot(E,IIpup(kk,:),'-bd',E,IIpdn(kk,:),'-rs',E,IIpup(kk,:)+IIpdn(kk,:),'-kx');
legend('IIpup','IIpdn','IIpup+IIpdn');
xlabel('energy (eV)');
ylabel('current ((q/h)/ev)')
title('Energy-resolved parallel current'); 
set(gca,'Fontsize',15)

figure(26)
plot(E,IIapup(kk,:),'-bd',E,IIapdn(kk,:),'-rs',E,IIapup(kk,:)+IIapdn(kk,:),'-kx');
legend('IIapup','IIapdn','IIapup+IIapdn');
xlabel('energy (eV)');
ylabel('current ((q/h)/ev )')
title('Energy-resolved antiparallel current'); 
set(gca,'Fontsize',15)


figure(27)
plot(E,IIpup(kk,:),'-bd',E,IIpdn(kk,:),'-r*',E,IIpup(kk,:)+IIpdn(kk,:),'-kx',E,IIapup(kk,:),'-gd',E,IIapdn(kk,:),'-ms',E,IIapup(kk,:)+IIapdn(kk,:),'-cx');
legend('IIpup','IIpdn','IIpup+IIpdn','IIapup','IIapdn','IIapup+IIapdn');
xlabel('energy (eV)');
ylabel('current ((q/h)/ev)')
title('Energy-resolved current'); 
set(gca,'Fontsize',15)




