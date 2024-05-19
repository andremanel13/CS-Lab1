close all;
clear all;

%% ADC parameters
Vdd  = 1;
Vcm  = Vdd/2;
Vref = 2*Vdd;
num_stages = 3;

gain  = [4.0 4.0 4.0];
nbits = [4 4 6];
nbits_adc = nbits(1) + (nbits(2)-1) + (nbits(3)-1);

%% Component errors
sigma_caps = 0.003;              % capacitor mismatch error - 0.3% of the unit value
sigma_gain = 0.003;              % gain error - 0.3% of the unit value
sigma_comp = 3e-2 / sqrt(2);     % comparator - offset voltage error - 30mV RMS

%% Input signal - sweep transfer function
Vlsb = Vref/2^nbits_adc; 
deltaVin = -Vdd: Vlsb/4 :Vdd;
np = max(size(deltaVin));

%% Monte Carlo Simulation
mc_simulations = 100;

% Accuracy and linearity variables
codes   = zeros(1, mc_simulations);
lin_inl = zeros(1, mc_simulations);
lin_dnl = zeros(1, mc_simulations);

tic
for mc = 1:mc_simulations
    % Simulate pipeline
    [dout] = pipeline_simulation(nbits, gain, np, deltaVin, sigma_caps, sigma_gain, sigma_comp);

    % Measure Vt
    Vt = deltaVin(find(dout(2:end)-dout(1:end-1))+1);
    ncodes = max(size(Vt));

    % INL e DNL
    Vlsb_real = (Vt(end)-Vt(1)) / (ncodes-1);
    inl = (Vt-(0:ncodes-1) * Vlsb_real - Vt(1)) / Vlsb_real;
    dnl = ((Vt(2:end) - Vt(1:end-1)) / Vlsb_real) - 1;

    codes(mc) = ncodes;
    lin_inl(mc) = nbits_adc + log(max(inl)-min(inl)) / log(2);
    lin_dnl(mc) = nbits_adc + log(max(dnl)-min(dnl)) / log(2);
end
toc

%% Plots

if mc_simulations==1
    % Transfer Function
    figure(1)
    plot(deltaVin,dout)
    grid on
    ylabel('Number of Codes')
    xlabel('deltaVin (V)')
    title(['Transfer Function with ' num2str(ncodes) ' codes'])
    
    figure(2)
    subplot(2,1,1)
    plot(inl)
    grid on
    title('INL')
    subplot(2,1,2)
    plot(dnl)
    grid on
    xlabel('Code Number')
    title('DNL')
else
    figure(3)
    subplot(2,1,1)
    hist(lin_inl, 100)
    grid on
    ylabel('Number of Cases')
    title('Linearity INL')
    subplot(2,1,2)
    hist(lin_dnl, 100)
    grid on
    ylabel('Number of Cases')
    title('Linearity DNL')
    
    figure(4)
    hist(2^nbits_adc - codes, 100)
    grid on
    ylabel('Number of Cases')
    title('Number of missing codes')
end


