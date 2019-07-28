function  sigma = selenerggy1D(E,U,t0)

EE = E-U;


rad = 4*t0^2-EE^2

if(rad>0.0)
    sigma = (EE - 1i*sqrt(rad))/2;
else
    sigma = (EE - sign(EE)*sqrt(-rad))/2;
end

end




% % test plot for self energy
% t0 =1;
% EE = linspace(-5,5,100);
% rad = 4*t0^2-(EE+2).^2;
% for i =1:100
% if(rad(i)>0.0)
%     sigma(i) = (EE(i) - 1i*sqrt(rad(i)))/2;
% else
%     sigma(i) = (EE(i) - sign(EE(i))*sqrt(-rad(i)))/2;
% end
% 
% end
% 
% plot(EE,real(sigma));
% hold on
% plot(EE,imag(sigma));
