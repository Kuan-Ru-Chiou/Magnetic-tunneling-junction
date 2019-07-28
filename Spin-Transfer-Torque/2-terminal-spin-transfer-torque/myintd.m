function out= myintd(E,V,mu1,mu2,UBB,tf,R,H,NL,NOx,NR,exchange,sx,sy,sz,Np,zplus,theta)
% 
if ((E-mu1)>0)
    fL = 0;
else
    fL = 1;
end

if ((E-mu2)>0)
    fR = 0;
else
    fR = 1;
end


%    fL = 1/(1+exp((E-mu1)/kT));
%    fR = 1/(1+exp((E-mu2)/kT));

U = [0.5*V*ones(1,NL) V*linspace(0.5,-0.5,NOx) -0.5*V*ones(1,NR)] + UBB;
%   U = [0.5*V*ones(1,NL) V*(0.5-[1:1:NOx]/(NOx+1)) -0.5*V*ones(1,NR)] + UBB;
U = kron(diag(U),eye(2));



% %datta-book-1D-self method 1
sigL = zeros(2*Np,2*Np); sigR = zeros(2*Np,2*Np);
sigL(1,1) = selfdatta11D(E,U(1,1),tf);
sigL(2,2) = selfdatta11D(E,U(2,2)+exchange,tf);
sigR(2*Np-1:2*Np,2*Np-1:2*Np) = R*[selfdatta11D(E,U(2*Np-1,2*Np-1),tf) 0;0 selfdatta11D(E,U(2*Np,2*Np)+exchange,tf)]*R';
      

%surface green function method2 iteration
% sigL = surface1Dite(E,tf,R,U,exchange,zplus,theta,sz,sx,Np,1)
% sigR = surface1Dite(E,tf,R,U,exchange,zplus,theta,sz,sx,Np,2)

%chen method 3
%      sigL = zeros(2*Np,2*Np);  sigR = zeros(2*Np,2*Np);
%      sigL(1,1) =  self1Dchern(E,U(1),tf);
%      sigL(2,2) =  self1Dchern(E,U(1)+exchange,tf);
%      sigR(2*Np-1:2*Np,2*Np-1:2*Np) = R*[self1Dchern(E,U(2*Np),tf) 0;0 self1Dchern(E,U(2*Np)+exchange,tf)]*R';



gamL = 1i*(sigL - sigL');  gamR = 1i*(sigR - sigR');
gamLin = fL*gamL;          gamRin = fR*gamR;

%green function non-self consistent
g = inv((E+zplus)*eye(2*Np)-sigL-sigR-H-U);
% g = inv((E)*eye(2*Np)-sigL-sigR-H-U);
%%%%using keldish Gless formula to check datta formula (checking ok)
Sigless = 1i*(gamL*fL + gamR*fR);   %%  antihermitian
Gless =  g*Sigless*g';              %%  antihermitian
% Gn = -1i*Gless;
%gn = g*(gamLin+gamRin)*g';
A = 1i*(g-g');
%density(j)=real(trace(gn)); %electron density
Tp = trace(gamL*g*gamR*g')*(fL-fR);  %parallel magnet

HH = H+U;

%spin current at interface point for calculation torque
%calculate Np-NR-3 spot's spincurrent

Jz1 = -trace(0.5*(HH(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx))*sz +sz*HH(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx)))*Gless(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1)))+...
    trace(0.5*(HH(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1))*sz +sz*HH(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1)))*Gless(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx)));

Jzy1 = Jz1;


Jx1 = -trace(0.5*(HH(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx))*sx +sx*HH(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx)))*Gless(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1)))+...
    trace(0.5*(HH(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1))*sx +sx*HH(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1)))*Gless(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx)));

Jxy1= Jx1;

Jy1 = -trace(0.5*(HH(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx))*sy +sy*HH(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx)))*Gless(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1)))+...
    trace(0.5*(HH(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1))*sy +sy*HH(2*(NL+NOx)-1:2*(NL+NOx),2*(NL+NOx+1)-1:2*(NL+NOx+1)))*Gless(2*(NL+NOx+1)-1:2*(NL+NOx+1),2*(NL+NOx)-1:2*(NL+NOx)));

Jyy1 = Jy1;


out=[Jxy1;Jyy1;Jzy1;Tp];

end