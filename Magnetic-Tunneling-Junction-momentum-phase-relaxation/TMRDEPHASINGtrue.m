clc;
clear all;
close all;
%  2018/1/9 TMR using Transmition function because numerical error
%%%check ok  low bias TMR very exact
%  contact1 is Z up and contact2 is down
%TMR low bias is very good
%constant and input
hbar = 1.06e-34; %plank constant
q = 1.6e-19;  %free electron charge
m0 = 9.1e-31; %free electon mass
mf =  0.73*m0        %FM effective mass
mO =  0.32*m0        %Oxside effective mass
a = 2e-10;       %MgO unit cell length ( I don't know real value , This is come from some article. Please check)
zplus = 1i*1e-12;
kT = 1d-9; %ev in room temp.
IE = q*q/(2*pi*hbar); % A/ev
tO = (hbar^2)/(2*mO*(a^2)*q); %Oxside hopping energy (eV)  (oxide a and Fe a using the same value because i dont know real value)
tf = (hbar^2)/(2*mf*(a^2)*q); %FM hopping energy (eV)
exchange =2.15;  %exchange splitting energy
Ef = 2.25;  %fermi-energy
Ec = 0;  %conduction band minimum
NL =2;  NOx =2;  NR =2;
Np = NL+NOx+NR;
Ub = Ef+0.93     %oxide barrier ev

DD =0;    %dephasing strength
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

% potential barrier
UBB = [zeros(1,NL) Ub*ones(1,NOx) zeros(1,NR)];
UB = kron(diag(UBB),eye(2));


H = zeros(2*Np,2*Np);
HA =H;

H = H+alphaL+alphaR+alphaox+beta+beta';
HA = HA + alphaL+alphaAR+alphaox+beta+beta';

 %energy grid
%   E = linspace(Ef-14*kT,Ef+5*kT,100);
% 
%   dE = E(2)-E(1);
%   NE = length(E); 
%bias
V = linspace(-0.1,0.1,8);
V(5:8) = -1.0*V(4:-1:1);
Nv = length(V);



% sigL = zeros(2*Np,2*Np);
% sigR = sigL;
% sigAR = sigL;
%loop over every energy
ii =1;  gg=1; ddE = zeros(Nv,1);
for i = 1:Nv  %bias
    
   mu1 = Ef+0.5*V(i);  mu2 = Ef-0.5*V(i);
   

       I = 0;  %coherent
       Ia = 0;  %coherent
       
       II1 = 0;  %noncoherent terminal L
       IIa1 = 0;  %noncoherent terminal L
        II2 = 0;  %noncoherent terminal R
       IIa2 = 0;  %noncoherent terminal R
       
       
%        TT=0;
%        TTA = 0;


E = linspace(min(mu1,mu2)-20*kT,max(mu1,mu2)+20*kT,1000);
dE = E(2)-E(1);
NE = length(E);

   

   
   
   
   
  
   for j = 1:NE
     
       sigL = zeros(2*Np,2*Np);
       sigR = sigL;
       sigAR = sigL;
       
       
       
      %Initialize guess value for dephasing model 
      sigB =zeros(2*Np); 
      sigBin = sigB; 
      sigBA = sigB;     
      sigBAin = sigB;
       
       
% fL = (2*mf*kT*q)/(2*pi*(hbar^2))*(log(1+exp((mu1-E(j))/kT))); %2D-transverse degenacy (unit 1/(meter*meter))
% fR = (2*mf*kT*q)/(2*pi*(hbar^2))*(log(1+exp((mu2-E(j))/kT))); %2D-transverse degenacy (unit 1/(meter*meter))
%    
    fL = 1/(1+exp((E(j)-mu1)/kT));    %fermi contact 1  
    fR = 1/(1+exp((E(j)-mu2)/kT));    %fermi contact 2
    
% if ((E(j)-mu1)>0)
%     fL = 0;
% else
%     fL = 1;
% end
% 
% if ((E(j)-mu2)>0)
%     fR = 0;
% else
%     fR = 1;
% end

    
    
      FFL(i,j)=fL;  %every bias fermi function
      FFR(i,j)= fR; %every bias fermi function
     
%      U = [0.5*V(i)*ones(1,NL) V(i)*linspace(0.5,-0.5,NOx) -0.5*V(i)*ones(1,NR)] + UBB;
     U = [0.5*V(i)*ones(1,NL) V(i)*(0.5-[1:1:NOx]/(NOx+1)) -0.5*V(i)*ones(1,NR)] + UBB;
     UUUU(i,:) =U;
     U = kron(diag(U),eye(2));
    
     
     
     
   
    
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
              
    
     
            change = 100;
     while change > 1e-8
         g = inv(((E(j)+zplus)*eye(2*Np))-H-U-sigL-sigR-sigB);
        sigBnew = diag(diag(DD*g));
%          sigBnew = DD*g;%phase relaxation
        change = sum(sum(abs(sigBnew-sigB)))/sum(sum(abs(sigBnew+sigB)));
%          change = sum(sum(abs(sigBnew-sigB)));
        sigB = sigB+0.5*(sigBnew-sigB);
     end
     
     Ap = 1i*(g-g');
     
     
       change = 100;
     while change > 1e-8
         gn = g*(gamLin+gamRin+sigBin)*g';
      sigBinnew = diag(diag(DD*gn));
%        sigBinnew = DD*gn; %phase relaxation  
         change = sum(sum(abs(sigBinnew-sigBin)))/sum(sum(abs(sigBinnew+sigBin)));
%         change = sum(sum(abs(sigBinnew-sigBin)));
        sigBin = sigBin+0.5*(sigBinnew-sigBin);
     end
     
     
     
     
            change = 100;
     while change > 1e-8
         gA = inv(((E(j)+zplus)*eye(2*Np))-HA-U-sigL-sigAR-sigBA);
        sigBAnew = diag(diag(DD*g));
%          sigBAnew = DD*gA;%phase relaxation
        change = sum(sum(abs(sigBAnew-sigBA)))/sum(sum(abs(sigBAnew+sigBA)));
%         change = sum(sum(abs(sigBAnew-sigBA)));
        sigBA = sigBA+0.5*(sigBAnew-sigBA);
     end
     
     Apa = 1i*(gA-gA');
     
     
       change = 100;
     while change > 1e-8
         gnA = gA*(gamLin+gamARin+sigBAin)*gA';
      sigBAinnew = diag(diag(DD*gnA));
%      sigBAinnew = DD*gnA; %phase relaxation  
         change = sum(sum(abs(sigBAinnew-sigBAin)))/sum(sum(abs(sigBAinnew+sigBAin)))
        
%         change = sum(sum(abs(sigBAinnew-sigBAin)))
        sigBAin = sigBAin+0.5*(sigBAinnew-sigBAin);
     end
    
      Gless =  1i*gn;              %%  antihermitian
      GAless =  1i*gnA;            %%  antihermitian
     %gn = g*(gamLin+gamRin)*g';
     % A = 1i*(g-g');
     %density(j)=real(trace(gn)); %electron density
     
     
     %%%%%%coherent transport
     Tp(j)= real(trace(gamL*g*gamR*g'));  %parallel magnet
     Tap(j)= real(trace(gamL*gA*gamAR*gA'));  %parallel magnet
     I = I + (IE*dE)*Tp(j)*(fL-fR);  %paralaell total current
     Ia = Ia + (IE*dE)*Tap(j)*(fL-fR);  %anti-paralaell total current
    %%%%%%%%
    
    %%%%non-coherent transport
    II1 = II1 + (IE*dE)*trace(gamLin*Ap-gamL*gn);%paralaell total current L
    II2 = II2 + (IE*dE)*trace(gamRin*Ap-gamR*gn);%paralaell total current R
    IIa1 = IIa1 + (IE*dE)*trace(gamLin*Apa-gamL*gnA);%anti-paralaell total current L
    IIa2 = IIa2 + (IE*dE)*trace(gamARin*Apa-gamAR*gnA);%anti-paralaell total current R
    %%%%%%%%%%%%%%%%%
    
    
    
    
%      TT = TT +Tp(j);
%      TTA = TTA +Tap(j);
     
       TT(i,j) = Tp(j);
       TTA(i,j) = Tap(j); 
     conductanceTT(i,j) = Tp(j);  %different bias energy resolved conductance
     conductanceTp(i,j) = Tap(j);
     
     
     
     HHHH = H+U;
     HHHA = HA+U;
     
     for kkk  = 1:Np-1
         JJJp0(j,kkk) = (-1*trace(HHHH(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)*Gless(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2) -HHHH(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2)*Gless(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)))*(IE*dE);
         JJJa0(j,kkk) = -(1*trace(HHHA(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)*GAless(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2) -HHHA(2*kkk-1:2*kkk,2*kkk+1:2*kkk+2)*GAless(2*kkk+1:2*kkk+2,2*kkk-1:2*kkk)))*(IE*dE);
     end

     
     
     
     
     
   end

   
      Ipp(i) = I;
      Iap(i) = Ia;
%        TTT(i) = TT; 
%        TTTA(i) = TTA;
end


%%%%energy current spectrum 1/eV
IIpp = TT.*(FFL-FFR);
IIap = TTA.*(FFL-FFR); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% IIpp = sum(IIpp,2).*(IE*dE);
% IIap = sum(IIap,2).*(IE*dE);



%check in and out terminal current conservation
fprintf('charge conservation for parallel %.f\n',sum(II1+II2));
fprintf('charge conservation for anti-parallel %.f\n',sum(IIa1+IIa2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













figure(1)
plot(V,Ipp,'--rs',V,Iap,'--bs');
xlabel('bias(eV)');
ylabel(' current (A)');
legend('Ipp','Iap');
title('current vs. bias')
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
% ylabel('antiparallel current');


% 
% figure(4)
%  kkk =diff(V)./diff(Iap);
%  plot(V(1:Nv-1),rrr);
%   ylim([0 500]);
% 
%   xlabel('voltage');
%  ylabel('dV/dI'); 

%TMR calculate
figure(5)
TTT = sum(TT,2);
TTTA = sum(TTA,2);
% TMR = (TTT-TTTA)./(TTTA);
TMR =(TTT./TTTA)-1;
plot(V,TMR*100,'b--'); 
xlabel('bias(eV)');
ylabel('TMR %'); 
title('TMR')



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


% figure(8)
% plot(V, conductanceTT,'b--',V,conductanceTp,'r--');
% xlabel('bias');
% ylabel('Total transmition')
% legend('parallel total transmition','antiparallel total transmition');



%band diagram for up dn channel
kk = input('bias number ');

xxx =[1:1:Np]; xxL =xxx([1:NL]); 
xxR =xxx(NL+NOx+1:Np);




figure(18)
plot(xxx,UUUU(kk,:)+[0*ones(1,NL) zeros(1,NOx) 0*ones(1,NR)],'--rs',xxx,UUUU(kk,:)+[exchange*ones(1,NL) zeros(1,NOx) exchange*ones(1,NR)],'--bd');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'g--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'g--');
plot(xxx,(Ef)*ones(1,Np),'g--');
hold off
legend('up-channel','dn-channel');
xlabel('position(nm)');
ylabel('Conduction band edge potential (eV)');
title('parallel-spin-dependent energy band');
text(0,Ef+0.5*V(kk),'mu1')
text(Np,Ef-0.5*V(kk),'mu2')


figure(19)
plot(xxx,UUUU(kk,:)+[0*ones(1,NL) zeros(1,NOx) exchange*ones(1,NR)],'--rs',xxx,UUUU(kk,:)+[exchange*ones(1,NL) zeros(1,NOx) 0*ones(1,NR)],'--bd');
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'g--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'g--');
plot(xxx,(Ef)*ones(1,Np),'g--');
hold off
legend('up-channel','dn-channel');
xlabel('position(nm)');
ylabel('Conduction band edge potential (eV)');
title('antiparallel-spin-dependent energy band');
text(0,Ef+0.5*V(kk),'mu1')
text(Np,Ef-0.5*V(kk),'mu2')


%device potentital 
figure(20)
plot([1:1:Np],UUUU(kk,:));
hold on
plot(xxL,(Ef+0.5*V(kk))*ones(1,NL),'g--');
plot(xxR,(Ef-0.5*V(kk))*ones(1,NR),'g--');
plot(xxx,(Ef)*ones(1,Np),'g--');
hold off
xlabel('position (nm)');
ylabel('potential (ev)');
title('Oxide barrier energy band');
text(0,Ef+0.5*V(kk),'mu1')
text(Np,Ef-0.5*V(kk),'mu2')



figure(10)
plot(E, conductanceTT(kk,:),'b--',E,conductanceTp(kk,:),'r--'); 
xlabel('energy');
ylabel('transmition')
legend('Energy-resolved parallel transmition','Energy-resolved antiparallel transmition'); 







