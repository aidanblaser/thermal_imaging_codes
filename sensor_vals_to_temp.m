function image_Temp = sensor_vals_to_temp(raw,Emissivity,ObjectDistance,AtmosphericTemperature)
% Converted by: Aidan Blaser (ablaser@ucsd.edu)
% Last Edited: 09/10/21
%
% Sources: https://github.com/detecttechnologies/thermal_base/blob/main/thermal_base/utils.py
%          https://github.com/gtatters/Thermimage/blob/master/R/raw2temp.R

    ReflectedApparentTemperature = AtmosphericTemperature;
    IRWindowTemperature=AtmosphericTemperature;
    IRWindowTransmission=1;
    RelativeHumidity=50;
    PlanckR1=14906.216;
    PlanckB= 1396.5;
    PlanckF=1;
    PlanckO=-7261;
    PlanckR2=0.010956882;
    %"""Convert raw values from the thermographic sensor sensor to temperatures in °C. Tested for Flir cams."""
    % this calculation has been ported to python from https://github.com/gtatters/Thermimage/blob/master/R/raw2temp.R
    % a detailed explanation of what is going on here can be found there

    % constants
    ATA1 = 0.006569;
    ATA2 = 0.01262;
    ATB1 = -0.002276;
    ATB2 = -0.00667;
    ATX = 1.9;

    % transmission through window (calibrated)
    emiss_wind = 1 - IRWindowTransmission;
    refl_wind = 0;
    % transmission through the air
    h2o = (RelativeHumidity / 100) * exp(1.5587 + 0.06939 * (AtmosphericTemperature)...
        - 0.00027816 * (AtmosphericTemperature)^2 ...
        + 0.00000068455 * (AtmosphericTemperature)^3);
    
    tau1 = ATX * exp(-sqrt(ObjectDistance / 2) * (ATA1 + ATB1 * sqrt(h2o))) + (1 - ATX) * ...
        exp(-sqrt(ObjectDistance / 2) * (ATA2 + ATB2 * sqrt(h2o)));

    tau2 = ATX * exp(-sqrt(ObjectDistance / 2) * (ATA1 + ATB1 * sqrt(h2o))) + (1 - ATX) * ...
        exp(-sqrt(ObjectDistance / 2) * (ATA2 + ATB2 * sqrt(h2o)));
    
    % radiance from the environment
    raw_refl1 = PlanckR1 / (PlanckR2 * (exp(PlanckB / (ReflectedApparentTemperature + 273.15)) - PlanckF)) - PlanckO;
    raw_refl1_attn = (1 - Emissivity) / Emissivity * raw_refl1;  % Reflected component

    raw_atm1 = (PlanckR1 / (PlanckR2 * (exp(PlanckB / (AtmosphericTemperature + 273.15)) - PlanckF)) - PlanckO);
    % Emission from atmosphere 1
    raw_atm1_attn = (1 - tau1) / Emissivity / tau1 * raw_atm1;  % attenuation for atmospheric 1 emission

    raw_wind = (PlanckR1 / (PlanckR2 * (exp(PlanckB / (IRWindowTemperature + 273.15)) - PlanckF)) - PlanckO);
    % Emission from window due to its own temp
    raw_wind_attn = (emiss_wind / Emissivity / tau1 / IRWindowTransmission * raw_wind);
    % Componen due to window emissivity

    raw_refl2 = (PlanckR1 / (PlanckR2 * (exp(PlanckB / (ReflectedApparentTemperature + 273.15)) - PlanckF)) - PlanckO);
    % Reflection from window due to external objects
    raw_refl2_attn = (refl_wind / Emissivity / tau1 / IRWindowTransmission * raw_refl2);
    % component due to window reflectivity

    raw_atm2 = (PlanckR1 / (PlanckR2 * (exp(PlanckB / (AtmosphericTemperature + 273.15)) - PlanckF)) - PlanckO);
    % Emission from atmosphere 2
    raw_atm2_attn = ((1 - tau2) / Emissivity / tau1 / IRWindowTransmission / tau2 * raw_atm2);
    % attenuation for atmospheric 2 emission
    
    raw_obj = (raw / Emissivity / tau1 / IRWindowTransmission / tau2 - raw_atm1_attn ...
    - raw_atm2_attn - raw_wind_attn - raw_refl1_attn - raw_refl2_attn);
    
    val_to_log = PlanckR1 ./ (PlanckR2 * (raw_obj + PlanckO)) + PlanckF;
    
    %temperature from radiance
    image_Temp = PlanckB ./ log(val_to_log) - 273.15;
end