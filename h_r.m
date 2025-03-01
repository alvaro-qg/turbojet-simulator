function value = h_r(T,alpha)

 h_r_air = 3.5*T - (1.4*(10^-5))*(T^2) + (7.467*(10^-9))*(T^3) + (3090/((exp(3090/T))-1));

 h_r_fuel = -149.054 + 4.47659*T + 4.00997*(10^-3)*(T^2) - 6.12432*(10^-7)*(T^3);
 
 %Enthalpy (h/r)
 value = (h_r_air + alpha*h_r_fuel)/(1+alpha);
 
end