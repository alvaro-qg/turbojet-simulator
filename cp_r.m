function value = cp_r(T,alpha)

 cp_r_air = 3.5 - (2.8*(10^-5))*T + (2.24*(10^-8))*(T^2) + ((3090/T)^2)*(exp(3090/T)/((exp(3090/T)-1)^2));
 
 cp_r_fuel = 4.47659 + (8.01994*(10^-3))*T - (1.8373*(10^-6))*(T^2);

 %Constant pressure specific heat (cp/r)
 value= (cp_r_air + alpha*cp_r_fuel) / (1+ alpha);
 
end