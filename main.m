clear all;
options = optimset('Display','off');

% TURBOJET DESIGN POINT SIMULATOR
% Real turbojet 
r=287.15; %Specific gas constant
%Flight altitude h=11000 m
h=11000;

%ISA Atmosphere Model
if(h<11000) %Troposphere
 T0=273.15+(15.04 - 0.00649*h);
 p0=101.29*((T0/288.08)^5.256);
elseif(h<25000 && h>=11000) %Lower stratosphere
 T0=273.15-56.46;
 p0=22.65*exp(1.73-0.000157*h); %Upper stratosphere
elseif(h>=25000)
 T0=273.15+(-131.21+0.00299*h);
 p0=2.488*((T0/216.6)^-11.388);
end

%Tabulated values of ISA for this case
% T0=216.5;
% p0=22612;
% a=295.1;

%Design parameters
M0=2;
Tt4=1373;
tau_c=2;

%Loss coeficcients
e_intake=0.075;
e_burner=0.06;
e_nozzle=0.02;

%Politropic efficiencies
n_compressor=0.88;
n_turbine=0.93;

%Nozzle velocity coefficient
phi_nozzle=0.98;

%Turbine cooling bleed
x=0.07;

%TURBOJET CYCLE
%Free stream 0
cp0=r*cp_r(T0,0); 
gamma0=cp0/(cp0-r); % Also we can consider gamma = 1.4 just for this part
a=sqrt(gamma0*r*T0); % Speed of sound at T0

h0=r*h_r(T0,0);
phi0=phi(T0,0);
p_r0=exp(phi0);

v0=M0*a; %Velocity of the air entering the engine
ht0 = h0 + (v0^2)/2; %Definition of total enthalpy

T_t0_guess=390.337; %Initial guess value calculated from the hand-made exercise

fun=@(T) (r*h_r(T,0))-ht0; % Function of the nonlinear equation evaluated at T
Tt0=fsolve(fun,T_t0_guess,options); %This function(fsolve) solves the nonlinear equation of the interpolation that
                                    %must be specified as F(T)=0. Where T
                                    %is the value to be determined and
                                    %T_guess is the initial point to search
                                    %values that the equation accepts
                                    %There may be som issues at convergence
                                    %of the solution, so in order to assure
                                    %that the final value that's been found
                                    %is good, we need to check the
                                    %residuals and see the order of
                                    %magnitude of the errors.
f_sol_ht0=r*h_r(Tt0,0); 
residual_1=f_sol_ht0-ht0; % Residual value from the difference between the computed value and the real value 

phit0=phi(Tt0,0);
p_rt0=exp(phit0);

p_t0=p0*(p_rt0/p_r0);

%Intake 0-2
ht2=ht0; %Total enthalpy is equivalent because it is a duct
Tt2=Tt0;
phit2=phit0;
p_rt2=p_rt0;

p_t2=p_t0*(1-e_intake); %Pressure losses

%Compressor 2-3
Tt3=tau_c*Tt2;

ht3=r*h_r(Tt3,0);
phit3=phi(Tt3,0);
p_rt3=exp(phit3);

p_t3=p_t2*((p_rt3/p_rt2)^n_compressor); %Adiabatic compression

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
ht5_mix=((ht5*(1+alpha)*(1-x))+ht3*x)/(((1+alpha)*(1-x))+x); %Energetic balance of the mixing at the exit of the turbine after bleed air re-incorporation

alpha_eff=alpha*(1-x); %Real fuel air ratio after the bleed re-incorporation

Tt5_guess=1000;
fun=@(T) (r*h_r(T,alpha_eff))-ht5_mix;
Tt5=fsolve(fun,Tt5_guess,options);
f_sol_ht5=r*h_r(Tt5,alpha_eff);
residual_2=f_sol_ht5-ht5_mix;

phit5=phi(Tt5,alpha_eff);
p_rt5=exp(phit5);

p_t5=p_t4*((p_rt5/p_rt4)^(1/n_turbine)); % Adiabatic expansion

%Nozzle 5-9
p9=p0; %Adapted nozzle

ht9=ht5_mix; %Equivalent enthalpy because it is a duct
Tt9=Tt5;
p_rt9=p_rt5;

p_t9=p_t5*(1-e_nozzle);

p_r9i = p_rt9*(p9/p_t9);

T9i_guess=493.62;
fun=@(T) (exp(phi(T,alpha_eff)))-p_r9i;
T9i=fsolve(fun,T9i_guess,options);

f_sol_phi9i=phi(T9i,alpha_eff);
f_sol_p_r9i=exp(f_sol_phi9i);
residual_3=f_sol_p_r9i-p_r9i;

h9i=r*h_r(T9i,alpha_eff);

v9i=sqrt(2*(ht9-h9i));
v9=phi_nozzle*v9i;

%Specific thrust
psi=((1+alpha_eff)*v9) - v0;
%Thermal,propulsive and overall efficiency
delta_K= 0.5*((1+alpha_eff)*(v9^2) - (v0^2)); %Kinetic energy diffrence
q_in=alpha_eff*hf_Tt4; %Specific heat given to the cycle

n_th=delta_K/q_in;

n_pr=(psi*v0)/delta_K;

n_o=(psi*v0)/delta_K;

%Thrust specific fuel consmption
C_TS=alpha_eff/psi;

%Display setup
fprintf('\n')
noms={'Free_stream','Intake','Compressor','Burner','Turbine','Nozzle'};
variables={'ht[kJ/kg]','Tt[K]','Pt[kPa]'};
engine_cycle=table([ht0/1000;Tt0;p_t0],[ht2/1000;Tt2;p_t2],[ht3/1000;Tt3;p_t3],[ht4/1000;Tt4;p_t4],[ht5/1000;Tt5;p_t5],[ht9/1000;Tt9;p_t9],'VariableNames', noms , 'RowNames', variables);
fprintf(['Turbojet engine cycle','\n'])
fprintf(['---------------------------------------','\n'])
disp(engine_cycle)
fprintf(['Turbojet performance parameters','\n'])
fprintf(['---------------------------------------','\n'])
fprintf(['    Specific thrust = ', num2str(psi),' m/s','\n'])
fprintf(['    Thermal efficiency = ', num2str(n_th),'\n'])
fprintf(['    Propulsive efficiency = ', num2str(n_pr),'\n'])
fprintf(['    Overall efficiency = ', num2str(n_o),'\n'])
fprintf(['    Thrust specific fuel consumption = ', num2str(C_TS) ' Kg/s/N','\n'])
fprintf('\n')
fprintf(['fsolver residuals and solutions','\n'])
fprintf(['---------------------------------------','\n'])
fprintf(['    fsolv T_t0 = ', num2str(Tt0),' residual = ',num2str(residual_1),'\n'])
fprintf(['    fsolv T_t5 = ', num2str(Tt5),' residual = ',num2str(residual_2),'\n'])
fprintf(['    fsolv T_9i = ', num2str(T9i),' residual = ',num2str(residual_3)],'\b')
fprintf('\n')

