close all;
clear all;

%% Parâmetros do ADC
Vdd  = 1;
Vcm  = Vdd/2;
Vref = 2*Vdd;
num_stages = 3;

gain  = [4.0 4.0 4.0];
nbits = [4 4 6];
nbits_adc = nbits(1) + (nbits(2)-1) + (nbits(3)-1);

%% Erros dos Componentes
sigma_caps = 0.003;             % erro de incompatibilidade de capacitor - 0.3% do valor unitário
sigma_gain = 0.003;             % erro do ganho - 0.3% do valor unitário
sigma_comp = 3e-2 / sqrt(2);    % comparador - erro de tensão de offset - 30mV RMS

%% Sinal de Entrada - Função de transferência Sweep
Vlsb = Vref/2^nbits_adc; 
deltaVin = -Vdd: Vlsb/4 :Vdd;
np = max(size(deltaVin));

%% Simulações de Monte Carlo
mc_simulations = 100;

% Variáveis de precisão e linearidade
codes   = zeros(1, mc_simulations);
lin_inl = zeros(1, mc_simulations);
lin_dnl = zeros(1, mc_simulations);

tic
for mc = 1:mc_simulations
    % Simular pipeline
    [dout] = pipeline_simulation(nbits, gain, np, deltaVin, sigma_caps, sigma_gain, sigma_comp);

    % Medir Vt
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
    % Função de Transferência
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


