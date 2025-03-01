function value = hf_r(T)
 
 r=287.15;

 hf0=4.3095*(10^7);
 
 hfe_r = -1607.2 + 4.47659*T + 4.00997*(10^-3)*(T^2) - 6.12432*(10^-7)*(T^3);
 
 % Fuel effective lower heating value
 value= (hf0/r) - hfe_r;
 
end