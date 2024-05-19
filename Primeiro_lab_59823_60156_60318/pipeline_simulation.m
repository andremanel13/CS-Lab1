function [dout] = pipeline_simulation (nbits, gain, np, deltaVin, sigma_caps, sigma_gain, sigma_comp)

    %% Parâmetros do ADC
    Vdd  = 1;
    Vref = 2*Vdd;
    Vcm  = Vdd/2;
    num_stages = 3;
    
    gain = [4.0 4.0 4.0];
    nbits = [4 4 6];
    nbits_adc = nbits(1) + (nbits(2)-1) + (nbits(3)-1);

    %% Outras variáveis
    dout = zeros(1,np);
    dout_stages = cell(1,3);
    
    Cup  = cell(1,3);
    Cun  = cell(1,3);
    Gain = cell(1,3);
    Voff = cell(1,3);

    % Geração de componentes considerando erros
    for i = 1:num_stages
        Cup {i} = 1 + randn(1,2^nbits(i))*sigma_caps;
        Cun {i} = 1 + randn(1,2^nbits(i))*sigma_caps;
        Gain{i} = gain(i) + randn*sigma_gain;
        Voff{i} = randn*sigma_comp;
    end

    % Conversão do sinal de entrada - differential to single ended
    Vp = Vcm + deltaVin/2;
    Vn = Vcm - deltaVin/2;

    for i = 1:num_stages
        Vref = Vref/2; % Vref passa para metade de andar em andar
        [dout_stages{i}, Vp, Vn] = stage_simulation(nbits(i), i, np, Vp, Vn, Vref, Vcm, Cup{i}, Cun{i}, Voff{i}, Gain{i});
    end
    
    % Correção digital 
    dout1 = binaryVectorToDecimal(dout_stages{1});
    dout2 = binaryVectorToDecimal(dout_stages{2});
    dout3 = binaryVectorToDecimal(dout_stages{3});

    dout = dout1*(2^8) + dout2*(2^5) -(2^2)*(2^5) + dout3 -(2^4);
end