function [sunSafeties, planetSafeties] = GearStress_PlanetSun(trial_sun, trial_planet,gearbox,motor,duty_cycle,t_div, end_runs)
%% NOTE
%INPUTS:
%pd = diametral pitch (don't confuse this as pitch diameter. Following
%notation from Shigleys)
%ring_n = number of teeth on ring gear
%sun_n = number of teeth on sun gear
%planet_n = number of teeth on planet gear
%
%OUTPUTS:
%sF_s = Factor of Safety from Bending Stress
%sH_s = Factor of Safety from Wear


%% USER-DEFINED CONSTANTS
%Track and Vehicle Variables
end_time = 25; %length of endurance in minutes
Ne = end_runs;%Number of en durance events ran

torque = (motor.torque*duty_cycle)./t_div;%value of average torque retrieved from OptimumLapSim

%Gear Variables


%Y_tab = xlsread('Lewis_Factor.xlsx'); %Lewis Form Factor Table
N_s = motor.speed * end_time * Ne * 3; %load cycle for sun gear
N_p = motor.speed / (trial_planet.teethNum/trial_sun.teethNum) * end_time * Ne; %load cycle for planet gear

%Material/Manufacturing Properties

% v_s %poisson's ratio for sun
% v_p %poisson's ratio for planet
% E_s %Young's modulus of sun
% E_p %Young's modulus of planet


%%INDIVIDUAL VARIABLES NEEDED AFTER THIS POINT
%% STRESS COEFFICIENT CALCULATIONS
V = (trial_sun.pitchDiameter/1000/2)*(motor.speed); %Pitch Line Velocity
Wt = torque*1000/trial_sun.pitchDiameter; %Transmitted Load

%Kv - Dynamic Factor
B = 0.25*(12-gearbox.qualityFactor)^(2/3);
A = 50 + 56*(1-B);
Kv = abs(((A + sqrt(200*V))/A)^B);

%Ks - Size Factor
%{
[~, J_sun]= getGeometryFactor(trial_sun);
[~, J_planetSun]= getGeometryFactor(trial_planet);
%}
[J_planetSun, J_sun] = getGeometryFactor(trial_planet,trial_sun);
%Ks_s = 1.192*(trial_sun.facewidth*sqrt(J_sun)/pd)^0.0535;
%Ks_p = 1.192*(trial_planet.sunStage.facewidth*sqrt(J_planetSun)/pd)^0.0535;

Ks_s =1;
Ks_p =1;

%Load Distribution Factor
KhB_s = getLoadDistributionFactor(true,trial_sun);
KhB_p = getLoadDistributionFactor(true,trial_planet);

%Rim Thickness Factor
Kb=1; %assume thick enough rim

%Overload Factor
Ko = 2;

%Bending Stress Cycle Factor
%{
if N_s < 10^7
    Yn_s = 3.517*N_s^(-0.0817);
else
    Yn_s = 1.3558*N_s.^(-0.0178);
end

if N_p < 10^7
    Yn_p = 3.517.*N_p.^(-0.0817);
else
    Yn_p = 1.3558.*N_p.^(-0.0178);
end
%}
Yn_s = 1;
Yn_p = 1;
%Stress-cycle Factor Equations
%{
if N_s < 10^7
    Zn_s = 1.249.*N_s.^(-0.0138);
else
    Zn_s = 1.4488.*N_s.^(-0.023);
end

if N_p < 10^7
    Zn_p = 1.249.*N_p.^(-0.0138);
else
    Zn_p = 1.4488.*N_p.^(-0.023);
end
%}
Zn_s =1;
Zn_p =1;
%Reliability Factor
Kr = 1; %.99 Reliability

%Temperature Factor 
Kt = 1;

%Surface Condition Factor
Zr = 1; %Subject to change from manufacturer

%Load Sharing Ratio
mN = 1;

%Pitting Resistance Geometry Factor
mG = trial_planet.teethNum/trial_sun.teethNum;
I_sun = (cosd(trial_sun.pressureAngle)*sind(trial_sun.pressureAngle))/(2*mN)*mG/(mG+1);
I_planetSun = (cosd(trial_sun.pressureAngle)*sind(trial_sun.pressureAngle))/(2*mN)/(mG+1);

%Elastic coefficient
%this value can be either found using Eq 14-13 or Table 14-8
v_s = 0.3;
v_p = v_s;
E_s = 2.05 * 10^5; %MPa
E_p = E_s;
%Cp = 2300; %Steel to steel [psi^(1/2)]
Ze = 190;

%AGMA Strength Equations (Assume material same for sun and planet rn)
%Uncomment either grade 1 or grade 2 equations
%Grade 1
St = 0.533*trial_planet.hardness + 88.3; %Allowable bending stress number [MPa]
Sc = 2.22*trial_planet.hardness+ 200; %Allowable contact stress number [MPa]
%Grade 2
%St = 0.703*Hb + 113; %Allowable bending stress number [psi]
%Sc =2.41*trial_sun.hardness + 237; %Allowable contact stress number [psi]

%Hardness Ratio Factor
Zw = 1;

%% STRESS AND ENDURANCE CALCULATIONS
%%Only sun is being analzyed right now as sun is to fail before planet

%%Gear Wear
sigmaH_s = Ze*sqrt(...
                     Wt*Ko*Kv*Ks_s*KhB_s*Zr/I_sun/...
                     (trial_sun.pitchDiameter*trial_sun.facewidth));
sunSafeties.flank = Sc*Zn_s./(Kt*Kr)./sigmaH_s; 

sigmaH_p = Ze*sqrt(...
                     Wt*Ko*Kv*Ks_p*KhB_p*Zr/I_planetSun/...
                     (trial_planet.pitchDiameter*trial_planet.facewidth));
planetSafeties.flank = Sc*Zn_p./(Kt*Kr)./sigmaH_p;

%Bending
sigmaF_s = Wt*Ko*Kv*Ks_s/trial_sun.facewidth*(KhB_s*Kb)/J_sun;
sunSafeties.root = St*Yn_s./(Kt*Kr)./sigmaF_s;

sigmaF_p = Wt*Ko*Kv*Ks_p/trial_planet.facewidth*(KhB_p*Kb)/J_planetSun;
planetSafeties.root = St*Yn_p./(Kt*Kr)./sigmaF_p;

%fprintf('Sun ring has max gear contact stress = %f kpa and a FOS = %f\n',sigmaC_s_metric, sH_s);
%fprintf('Sun ring has max gear bending stress = %f kpa and a FOS = %f\n',sigma_metric, sF_s);
end