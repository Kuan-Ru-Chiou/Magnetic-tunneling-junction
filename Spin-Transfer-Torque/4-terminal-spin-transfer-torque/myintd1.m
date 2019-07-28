function out= myintd1(E,V,mu1,mu2,mu3,mu4,UBB,tf,R,H,N1,NOx,N2,N3,N4,exchange,sx,sy,sz,Np,zplus,tO)
% 

   
 
        
  
 %initialize contact self-energy
sig1 = zeros(2*Np,2*Np);
sig2 = sig1;
sig3 = sig1;
sig4 = sig1;
    
       
       
%  f3 = 1/(1+exp((E(j)-mu3)/kT));
%  f4 = 1/(1+exp((E(j)-mu4)/kT));
%  f1 = 1/(1+exp((E(j)-mu1)/kT));
% f2 = 1/(1+exp((E(j)-mu2)/kT));
       
if ((E-mu1)>0)
    f1 = 0;
else
    f1 = 1;
end

if ((E-mu2)>0)
    f2 = 0;
else
    f2 = 1;
end

if ((E-mu3)>0)
    f3 = 0;
else
    f3 = 1;
end

if ((E-mu4)>0)
    f4 = 0;
else
    f4 = 1;
end



 %%Total potential (bias +barrier)
 U = [zeros(1,N1) V*linspace(0.5,-0.5,NOx) zeros(1,N2) 0.5*V*ones(1,N3) -0.5*V*ones(1,N4)] + UBB;
%      U = [zeros(1,N1) V*(0.5-[1:1:NOx]/(NOx+1)) zeros(1,N2) 0.5*V*ones(1,N3) -0.5*V*ones(1,N4)] + UBB;
     U = kron(diag(U),eye(2));
    
      %%self-energy-1D-exact
     sig1(1,1) = selfdatta1D(E,U(1),tO);
     sig1(2,2) = selfdatta1D(E,U(1),tO);
     
     sig2(2*(N1+NOx+N2)-1,2*(N1+NOx+N2)-1) = selfdatta1D(E,U(2*(N1+NOx+N2)-1,2*(N1+NOx+N2)-1),tO);
     sig2(2*(N1+NOx+N2),2*(N1+NOx+N2)) = selfdatta1D(E,U(2*(N1+NOx+N2),2*(N1+NOx+N2)),tO);
     
     sig3(2*(N1+NOx+N2+N3-1)+1:2*(N1+NOx+N2+N3-1)+2,2*(N1+NOx+N2+N3-1)+1:2*(N1+NOx+N2+N3-1)+2) = [selfdatta1D(E,0.5*V,tf) 0;0 selfdatta1D(E,0.5*V+exchange,tf)];
     sig4(2*(Np-1)+1:2*(Np-1)+2,2*(Np-1)+1:2*(Np-1)+2) = R*[selfdatta1D(E,-0.5*V,tf) 0;0 selfdatta1D(E,-0.5*V+exchange,tf)]*R';
     
     
     
     gam1 = 1i*(sig1 - sig1');  gam2 = 1i*(sig2 - sig2'); 
     gam3 = 1i*(sig3-sig3');    gam4 = 1i*(sig4-sig4');  
     gam1in = f1*gam1;          gam2in = f2*gam2;
     gam3in = f3*gam3;          gam4in = f4*gam4;
     
     
      %green function non-self consistent
     g = inv((E+zplus)*eye(2*Np)-sig1-sig2-sig3-sig4-H-U);
     
      %%%%using keldish Gless formula to check datta formula (checking ok)
      Sigless = 1i*(gam1*f1 + gam2*f2 + gam3*f3 + gam4*f4);   %%  antihermitian
      Gless =  g*Sigless*g';                                  %%  antihermitian
%       
%       if(max(max(abs(Gless+Gless')))>=1e-14)
%           error('wrong')
%       end
%       
      
       A = 1i*(g-g');
       
       Tp= real(trace(gam3*g*gam4*g'))*(f3-f4);  %transmittion function
       
      
       
       
        
     %spin current at interface point for calculation torque
     %calculate Np-NR-3 spot's spincurrent
     HH = H + U;
    
     
        %%%%%%%%%check spin current conserve in n,n-1,n+1 points  spin-conservation
%%%%%%%%%
      Jz1 = -trace(0.5*(HH(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2)*sz +sz*HH(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2))*Gless(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2)-...
           0.5*(HH(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2)*sz +sz*HH(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2))*Gless(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2));
      Jzy1 = Jz1; 
     
       Jx1 = -trace(0.5*(HH(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2)*sx +sx*HH(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2))*Gless(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2)-...
           0.5*(HH(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2)*sx +sx*HH(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2))*Gless(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2));
    
     Jxy1 = Jx1;
     
      Jy1 = -trace(0.5*(HH(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2)*sy +sy*HH(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2))*Gless(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2)-...
           0.5*(HH(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2)*sy +sy*HH(2*(N1+NOx-1)+1:2*(N1+NOx-1)+2,2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2))*Gless(2*(N1+NOx+N2+N3)+1:2*(N1+NOx+N2+N3)+2,2*(N1+NOx-1)+1:2*(N1+NOx-1)+2));
   
     Jyy1 = Jy1;


out=[Jxy1;Jyy1;Jzy1;Tp];
end