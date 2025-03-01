function value = phi(T,alpha)

 phi_air = 3.5*log(T) - (2.8*(10^-5))*T + (1.12*(10^-8))*(T^2) + (3090/(T*((exp(3090/T))-1))) - log(((exp(3090/T))-1)/exp(3090/T));

 phi_fuel = 4.47659*log(T) + 8.01994*(10^-3)*T + 9.18648*(10^-7)*(T^2);
 
 %Pressure ratio (phi)
 value = (phi_air + alpha*phi_fuel)/(1+alpha);
 
end