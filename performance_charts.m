clear all;
options = optimset('Display','off');

%Turbojet(black line) design parameters 
%Tt4=1373;
%M0=2;
%delta_Tc=389,7;

%Design parameters
Tt4=1373;
M0=[0,2,3];
delta_Tc=[100:10:800];

f1=figure; % psi vs delta_Tc (Fixed MO)
f2=figure; % C_TS vs delta_Tc (Fixed M0)
f3=figure; % psi vs M0 (Fixed delta_Tc)
f4=figure; % C_TS vs M0 (Fixed delta Tc)

r=287.15;
h=11000;
if(h<11000) 
 T0=273.15+(15.04 - 0.00649*h);
 p0=101.29*((T0/288.08)^5.256);
elseif(h<25000 && h>=11000) 
 T0=273.15-56.46;
 p0=22.65*exp(1.73-0.000157*h); 
elseif(h>=25000)
 T0=273.15+(-131.21+0.00299*h);
 p0=2.488*((T0/216.6)^-11.388);
end

cp0=r*cp_r(T0,0);
gamma0=cp0/(cp0-r);
a=sqrt(gamma0*r*T0);

e_intake=0.075;
e_burner=0.06;
e_nozzle=0.02;
n_compressor=0.88;
n_turbine=0.93;
phi_nozzle=0.98;
x=0.07;

psi=[];
C_TS=[];

t=1;
i=1;

while(t<=length(M0))
 i=1;
 ind(t)=0;
 while(i<= length(delta_Tc))
 %Free stream 0
 h0=r*h_r(T0,0);
 phi0=phi(T0,0);
 p_r0=exp(phi0);

 v0=(M0(t))*a;
 ht0 = h0 + (v0^2)/2;
 T_t0_guess=390.337;
 fun=@(T) (r*h_r(T,0))-ht0;
 Tt0=fsolve(fun,T_t0_guess,options);
 phit0=phi(Tt0,0);
 p_rt0=exp(phit0);
 p_t0=p0*(p_rt0/p_r0);

 %Intake 0-2
 ht2=ht0;
 Tt2=Tt0;
 phit2=phit0;
 p_rt2=p_rt0;
 p_t2=p_t0*(1-e_intake);

 %Compressor 2-3
 tau_c=1+(delta_Tc(i)/Tt2);
 Tt3=tau_c*Tt2;
 ht3=r*h_r(Tt3,0);
 phit3=phi(Tt3,0);
 p_rt3=exp(phit3);
 p_t3=p_t2*((p_rt3/p_rt2)^n_compressor);

%Burner 3-4
 p_t4=p_t3*(1-e_burner);
 
 ht4_air=r*h_r(Tt4,0);
 hf_Tt4=r*hf_r(Tt4); %Lower heating value
 alpha=(ht4_air-ht3)/hf_Tt4; %Real fuel air ratio at the burner. Energetic balance at burner
 
 ht4=r*h_r(Tt4,alpha); %Air+fuel mixture at the exit of the burner
 phit4=phi(Tt4,alpha);
 p_rt4=exp(phit4);
 
 %Turbine 4-5
 ht5=ht4-((ht3-ht2)/((1+alpha)*(1-x))); %Turbine-compressor matching before bleed massflow reinjection
 ht5_mix=(ht5*(1+alpha)*(1-x)+ht3*x)/((1+alpha)*(1-x)+x); %Energetic balance of the mixing at the exit of the turbine after bleed air re-incorporation
 
 alpha_eff=alpha*(1-x); %Real fuel air ratio after the bleed re-incorporation
 
 Tt5_guess=1047.097;
 fun=@(T) (r*h_r(T,alpha_eff))-ht5_mix;
 Tt5=fsolve(fun,Tt5_guess,options);
 f_sol_ht5=r*h_r(Tt5,alpha_eff);
 residual_2=f_sol_ht5-ht5_mix;
 
 phit5=phi(Tt5,alpha_eff);
 p_rt5=exp(phit5);
 
 p_t5=p_t4*((p_rt5/p_rt4)^(1/n_turbine)); % Adiabatic expansion
 
 %Nozzle 5-9
 p9=p0; 
 p_t9=p_t5*(1-e_nozzle);
 ht9=ht5;
 Tt9=Tt5;
 p_rt9=p_rt5;
 p_r9i = p_rt9*(p9/p_t9);
 T9i_guess=493.62;
 fun=@(T) (exp(phi(T,alpha_eff)))-p_r9i;
 T9i=fsolve(fun,T9i_guess,options);
 h9i=r*h_r(T9i,alpha_eff);
 v9i=sqrt(2*(ht9-h9i));
 v9=phi_nozzle*v9i;

 psi(t,i)=((1-alpha_eff).*v9) - v0;
 if(psi(t,i)>=0)
  ind(t)=ind(t)+1;
  C_TS(t,i)=alpha_eff/psi(t,i);
 end
 
 i=i+1;
 end
 t=t+1;
end

% 1 line for each fixed M0 -> 3 lines
psi_1=[];C_TS_1=[];delta_Tc_1=[];
psi_2=[];C_TS_2=[];delta_Tc_2=[];
psi_3=[];C_TS_3=[];delta_Tc_3=[];

for i=1:ind(1)
 psi_1(i)=psi(1,i);
 C_TS_1(i)=C_TS(1,i);
 delta_Tc_1(i)=delta_Tc(i);
end

for i=1:ind(2)
 psi_2(i)=psi(2,i);
 C_TS_2(i)=C_TS(2,i);
 delta_Tc_2(i)=delta_Tc(i);
end

for i=1:ind(3)
 psi_3(i)=psi(3,i);
 C_TS_3(i)=C_TS(3,i);
 delta_Tc_3(i)=delta_Tc(i);
end

figure(f1)
hold on
title('psi vs delta T_c (Fixed MO)')
xlabel('T_c')
ylabel('Specific thrust')
plot(delta_Tc_1,psi_1,': black');
plot(delta_Tc_2,psi_2,'black');
plot(delta_Tc_3,psi_3,'green');
hold off

figure(f2)
hold on
xlabel('T_c')
ylabel('Thrust specific fuel consumption') 
title('C_T_S vs delta T_c (Fixed M0)')
plot(delta_Tc_1,C_TS_1,': black');
plot(delta_Tc_2,C_TS_2,'black');
plot(delta_Tc_3,C_TS_3,'green');
hold off

%Repeat the design parameters for diferent Tt4 values
Tt4=[1000,1800];
M0=2;
delta_Tc=[100:10:800];

psi=[];
C_TS=[];

t=1;
i=1;

while(t<=length(Tt4))
 i=1;
 ind(t)=0;
 while(i<= length(delta_Tc))
 %Free stream 0
 h0=r*h_r(T0,0);
 phi0=phi(T0,0);
 p_r0=exp(phi0);

 v0=M0*a;
 ht0 = h0 + (v0^2)/2;
 T_t0_guess=390.337;
 fun=@(T) (r*h_r(T,0))-ht0;
 Tt0=fsolve(fun,T_t0_guess,options);
 phit0=phi(Tt0,0);
 p_rt0=exp(phit0);
 p_t0=p0*(p_rt0/p_r0);

 %Intake 0-2
 ht2=ht0;
 Tt2=Tt0;
 phit2=phit0;
 p_rt2=p_rt0;
 p_t2=p_t0*(1-e_intake);

 %Compressor 2-3
 tau_c=1+(delta_Tc(i)/Tt2);
 Tt3=tau_c*Tt2;
 ht3=r*h_r(Tt3,0);
 phit3=phi(Tt3,0);
 p_rt3=exp(phit3);
 p_t3=p_t2*((p_rt3/p_rt2)^n_compressor);

 %Burner 3-4
 p_t4=p_t3*(1-e_burner);
 
 ht4_air=r*h_r(Tt4(t),0);
 hf_Tt4=r*hf_r(Tt4(t)); %Lower heating value
 alpha=(ht4_air-ht3)/hf_Tt4; %Real fuel air ratio at the burner. Energetic balance at burner
 
 ht4=r*h_r(Tt4(t),alpha); %Air+fuel mixture at the exit of the burner
 phit4=phi(Tt4(t),alpha);
 p_rt4=exp(phit4);
 
 %Turbine 4-5
 ht5=ht4-((ht3-ht2)/((1+alpha)*(1-x))); %Turbine-compressor matching before bleed massflow reinjection
 ht5_mix=(ht5*(1+alpha)*(1-x)+ht3*x)/((1+alpha)*(1-x)+x); %Energetic balance of the mixing at the exit of the turbine after bleed air re-incorporation
 
 alpha_eff=alpha*(1-x); %Real fuel air ratio after the bleed re-incorporation
 
 Tt5_guess=1047.097;
 fun=@(T) (r*h_r(T,alpha_eff))-ht5_mix;
 Tt5=fsolve(fun,Tt5_guess,options);
 f_sol_ht5=r*h_r(Tt5,alpha_eff);
 residual_2=f_sol_ht5-ht5_mix;
 
 phit5=phi(Tt5,alpha_eff);
 p_rt5=exp(phit5);
 
 p_t5=p_t4*((p_rt5/p_rt4)^(1/n_turbine)); % Adiabatic expansion
 
 %Nozzle 5-9
 p9=p0; 
 p_t9=p_t5*(1-e_nozzle);
 ht9=ht5;
 Tt9=Tt5;
 p_rt9=p_rt5;
 p_r9i = p_rt9*(p9/p_t9);
 T9i_guess=493.62;
 fun=@(T) (exp(phi(T,alpha_eff)))-p_r9i;
 T9i=fsolve(fun,T9i_guess,options);
 h9i=r*h_r(T9i,alpha_eff);
 v9i=sqrt(2*(ht9-h9i));
 v9=phi_nozzle*v9i;

 psi(t,i)=((1-alpha_eff).*v9) - v0;
 if(psi(t,i)>=0)
  ind(t)=ind(t)+1;
  C_TS(t,i)=alpha_eff/psi(t,i);
 end
 i=i+1;
 end
 t=t+1;
end

% 1 line for each Tt4 -> 2 lines
psi_1=[];C_TS_1=[];delta_Tc_1=[];
psi_2=[];C_TS_2=[];delta_Tc_2=[];

for i=1:ind(1)
 psi_1(i)=psi(1,i);
 C_TS_1(i)=C_TS(1,i);
 delta_Tc_1(i)=delta_Tc(i);
end

for i=1:ind(2)
 psi_2(i)=psi(2,i);
 C_TS_2(i)=C_TS(2,i);
 delta_Tc_2(i)=delta_Tc(i);
end

figure(f1)
hold on
plot(delta_Tc_1,psi_1,'blue');
plot(delta_Tc_2,psi_2,'red');
hold off

figure(f2)
hold on
plot(delta_Tc_1,C_TS_1,'blue');
plot(delta_Tc_2,C_TS_2,'red');
hold off

figure(f1)
hold on
legend('ground(M0=0)','turbojet','M0=3','Tt4=1000 K','Tt4=1800 K')

figure(f2)
hold on
legend('ground(M0=0)','turbojet','M0=3','Tt4=1000 K','Tt4=1800 K')
xlim([delta_Tc(1) delta_Tc(length(delta_Tc))])
ylim([0 8*10^-5])

Tt4=1373;
M0=[0:0.05:3];
delta_Tc=[100,389.7,800];

psi=[];
C_TS=[];

t=1;
i=1;

while(t<=length(delta_Tc))
 i=1;
 ind(t)=0;
 while(i<= length(M0))
 %Free stream 0
 h0=r*h_r(T0,0);
 phi0=phi(T0,0);
 p_r0=exp(phi0);

 v0=(M0(i))*a;
 ht0 = h0 + (v0^2)/2;
 T_t0_guess=390.337;
 fun=@(T) (r*h_r(T,0))-ht0;
 Tt0=fsolve(fun,T_t0_guess,options);
 phit0=phi(Tt0,0);
 p_rt0=exp(phit0);
 p_t0=p0*(p_rt0/p_r0);

 %Intake 0-2
 ht2=ht0;
 Tt2=Tt0;
 phit2=phit0;
 p_rt2=p_rt0;
 p_t2=p_t0*(1-e_intake);

 %Compressor 2-3
 tau_c=1+(delta_Tc(t)/Tt2);
 Tt3=tau_c*Tt2;
 ht3=r*h_r(Tt3,0);
 phit3=phi(Tt3,0);
 p_rt3=exp(phit3);
 p_t3=p_t2*((p_rt3/p_rt2)^n_compressor);

 %Burner 3-4
 p_t4=p_t3*(1-e_burner);
 
 ht4_air=r*h_r(Tt4,0);
 hf_Tt4=r*hf_r(Tt4); %Lower heating value
 alpha=(ht4_air-ht3)/hf_Tt4; %Real fuel air ratio at the burner. Energetic balance at burner
 
 ht4=r*h_r(Tt4,alpha); %Air+fuel mixture at the exit of the burner
 phit4=phi(Tt4,alpha);
 p_rt4=exp(phit4);

 %Turbine 4-5
 ht5=ht4-((ht3-ht2)/((1+alpha)*(1-x))); %Turbine-compressor matching before bleed massflow reinjection
 ht5_mix=(ht5*(1+alpha)*(1-x)+ht3*x)/((1+alpha)*(1-x)+x); %Energetic balance of the mixing at the exit of the turbine after bleed air re-incorporation
 
 alpha_eff=alpha*(1-x); %Real fuel air ratio after the bleed re-incorporation
 
 Tt5_guess=1047.097;
 fun=@(T) (r*h_r(T,alpha_eff))-ht5_mix;
 Tt5=fsolve(fun,Tt5_guess,options);
 f_sol_ht5=r*h_r(Tt5,alpha_eff);
 residual_2=f_sol_ht5-ht5_mix;
 
 phit5=phi(Tt5,alpha_eff);
 p_rt5=exp(phit5);
 
 p_t5=p_t4*((p_rt5/p_rt4)^(1/n_turbine)); % Adiabatic expansion

 %Nozzle 5-9
 p9=p0; 
 p_t9=p_t5*(1-e_nozzle);
 ht9=ht5;
 Tt9=Tt5;
 p_rt9=p_rt5;
 p_r9i = p_rt9*(p9/p_t9);
 T9i_guess=493.62;
 fun=@(T) (exp(phi(T,alpha_eff)))-p_r9i;
 T9i=fsolve(fun,T9i_guess,options);
 h9i=r*h_r(T9i,alpha_eff);
 v9i=sqrt(2*(ht9-h9i));
 v9=phi_nozzle*v9i;

 psi(t,i)=((1-alpha_eff).*v9) - v0;
 if(psi(t,i)>=0)
  ind(t)=ind(t)+1;
  C_TS(t,i)=alpha_eff/psi(t,i);
 end
 i=i+1;
 end
 t=t+1;
end

% 1 line for each delta Tc -> 3 lines 
psi_1=[];C_TS_1=[];M0_1=[];
psi_2=[];C_TS_2=[];M0_2=[];
psi_3=[];C_TS_3=[];M0_3=[];

for i=1:ind(1)
 psi_1(i)=psi(1,i);
 C_TS_1(i)=C_TS(1,i);
 M0_1(i)=M0(i);
end

for i=1:ind(2)
 psi_2(i)=psi(2,i);
 C_TS_2(i)=C_TS(2,i);
 M0_2(i)=M0(i);
end

for i=1:ind(3)
 psi_3(i)=psi(3,i);
 C_TS_3(i)=C_TS(3,i);
 M0_3(i)=M0(i);
end

figure(f3)
hold on
title('psi vs M0 (Fixed delta T_c)')
xlabel('M_0')
ylabel('Specific thrust')
plot(M0_1,psi_1,'cyan');
plot(M0_2,psi_2,'black');
plot(M0_3,psi_3,'magenta');
hold off

figure(f4)
hold on
title('C_T_S vs M_0 (Fixed delta T_c)')
xlabel('M_0')
ylabel('Thrust specific fuel consumption')
plot(M0_1,C_TS_1,'cyan');
plot(M0_2,C_TS_2,'black');
plot(M0_3,C_TS_3,'magenta');
hold off

%Repeat the design parameters for diferent Tt4 values
Tt4=[1000,1800];
M0=[0:0.05:3];
delta_Tc=389.7;

psi=[];
C_TS=[];

t=1;
i=1;

while(t<=length(Tt4))
 i=1;
 ind(t)=0;
 while(i<= length(M0))
 %Free stream 0
 h0=r*h_r(T0,0);
 phi0=phi(T0,0);
 p_r0=exp(phi0);

 v0=M0((i))*a;
 ht0 = h0 + (v0^2)/2;
 T_t0_guess=390.337;
 fun=@(T) (r*h_r(T,0))-ht0;
 Tt0=fsolve(fun,T_t0_guess,options);
 phit0=phi(Tt0,0);
 p_rt0=exp(phit0);
 p_t0=p0*(p_rt0/p_r0);

 %Intake 0-2
 ht2=ht0;
 Tt2=Tt0;
 phit2=phit0;
 p_rt2=p_rt0;
 p_t2=p_t0*(1-e_intake);

 %Compressor 2-3
 tau_c=1+(delta_Tc/Tt2);
 Tt3=tau_c*Tt2;
 ht3=r*h_r(Tt3,0);
 phit3=phi(Tt3,0);
 p_rt3=exp(phit3);
 p_t3=p_t2*((p_rt3/p_rt2)^n_compressor);

 %Burner 3-4
 p_t4=p_t3*(1-e_burner);
 
 ht4_air=r*h_r(Tt4(t),0);
 hf_Tt4=r*hf_r(Tt4(t)); %Lower heating value
 alpha=(ht4_air-ht3)/hf_Tt4; %Real fuel air ratio at the burner. Energetic balance at burner
 
 ht4=r*h_r(Tt4(t),alpha); %Air+fuel mixture at the exit of the burner
 phit4=phi(Tt4(t),alpha);
 p_rt4=exp(phit4);
 
 %Turbine 4-5
 ht5=ht4-((ht3-ht2)/((1+alpha)*(1-x))); %Turbine-compressor matching before bleed massflow reinjection
 ht5_mix=(ht5*(1+alpha)*(1-x)+ht3*x)/((1+alpha)*(1-x)+x); %Energetic balance of the mixing at the exit of the turbine after bleed air re-incorporation
 
 alpha_eff=alpha*(1-x); %Real fuel air ratio after the bleed re-incorporation
 
 Tt5_guess=1047.097;
 fun=@(T) (r*h_r(T,alpha_eff))-ht5_mix;
 Tt5=fsolve(fun,Tt5_guess,options);
 f_sol_ht5=r*h_r(Tt5,alpha_eff);
 residual_2=f_sol_ht5-ht5_mix;
 
 phit5=phi(Tt5,alpha_eff);
 p_rt5=exp(phit5);
 
 p_t5=p_t4*((p_rt5/p_rt4)^(1/n_turbine)); % Adiabatic expansion

 %Nozzle 5-9
 p9=p0; 
 p_t9=p_t5*(1-e_nozzle);
 ht9=ht5;
 Tt9=Tt5;
 p_rt9=p_rt5;
 p_r9i = p_rt9*(p9/p_t9);
 T9i_guess=493.62;
 fun=@(T) (exp(phi(T,alpha_eff)))-p_r9i;
 T9i=fsolve(fun,T9i_guess,options);
 h9i=r*h_r(T9i,alpha_eff);
 v9i=sqrt(2*(ht9-h9i));
 v9=phi_nozzle*v9i;

 psi(t,i)=((1-alpha_eff).*v9) - v0;
 if(psi(t,i)>=0)
  ind(t)=ind(t)+1;
  C_TS(t,i)=alpha_eff/psi(t,i);
 end
 i=i+1;
 end
 t=t+1;
end

% 1 line for each Tt4 -> 2 lines
psi_1=[];C_TS_1=[];M0_1=[];
psi_2=[];C_TS_2=[];M0_2=[];

for i=1:ind(1)
 psi_1(i)=psi(1,i);
 C_TS_1(i)=C_TS(1,i);
 M0_1(i)=M0(i);
end

for i=1:ind(2)
 psi_2(i)=psi(2,i);
 C_TS_2(i)=C_TS(2,i);
 M0_2(i)=M0(i);
end

figure(f3)
hold on
plot(M0_1,psi_1,'blue');
plot(M0_2,psi_2,'red');
hold off

figure(f4)
hold on
plot(M0_1,C_TS_1,'blue');
plot(M0_2,C_TS_2,'red');
hold off

figure(f3)
hold on
legend('T_c=100','turbojet','T_c=800','Tt4=1000 K','Tt4=1800 K')

figure(f4)
hold on
legend('T_c=100','turbojet','T_c=800','Tt4=1000 K','Tt4=1800 K')
xlim([M0(1) M0(length(M0))])
ylim([0 8*10^-5])




