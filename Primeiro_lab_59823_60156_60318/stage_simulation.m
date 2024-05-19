function [bit, Vresp, Vresn] = stage_simulation (nbits, nstage, np, Vp, Vn, Vref, Vcm, Cup, Cun, Voff, gain)

    bit = zeros(np, nbits);
    Vcp = zeros(np, nbits+1);
    Vcn = zeros(np, nbits+1);

    % Geração de Condensadores Equivalentes
    for i = 0:nbits-1
        Cbinp(i+1) = sum (Cup((2^i):(2^(i+1)-1)));
        Cbinm(i+1) = sum (Cun((2^i):(2^(i+1)-1)));
    end

    Cbinp = fliplr(Cbinp);
    Cbinm = fliplr(Cbinm);
    Cup_total = sum(Cup);
    Cun_total = sum(Cun);

    % Inicialização de tensões
    Vrp = Vcm + Vref/2;
    Vrn = Vcm - Vref/2;

    if nstage == 1
        Vcp(:,1) = 2*Vcm - Vp;
        Vcn(:,1) = 2*Vcm - Vn;
    else
        Vcp(:,1) = Vp;
        Vcn(:,1) = Vn;
    end
    
    for i = 1:nbits
        if nstage == 1
            bit(:,i) = (Vcn(:,i)-Vcp(:,i)) >= Voff;
            Vcp(:,i+1) = Vcp(:,i) - (Cbinp(i)/Cup_total)*Vcm + Vrn*not(bit(:,i))*(Cbinp(i)/Cup_total) + Vrp*bit(:,i)*(Cbinp(i)/Cup_total);
            Vcn(:,i+1) = Vcn(:,i) - (Cbinm(i)/Cun_total)*Vcm + Vrp*not(bit(:,i))*(Cbinm(i)/Cun_total) + Vrn*bit(:,i)*(Cbinm(i)/Cun_total);
        else
            bit(:,i) = (Vcp(:,i)-Vcn(:,i)) >= Voff;
            Vcp(:,i+1) = Vcp(:,i) - (Cbinp(i)/Cup_total)*Vcm + Vrp*not(bit(:,i))*(Cbinp(i)/Cup_total) + Vrn*bit(:,i)*(Cbinp(i)/Cup_total);
            Vcn(:,i+1) = Vcn(:,i) - (Cbinm(i)/Cun_total)*Vcm + Vrn*not(bit(:,i))*(Cbinm(i)/Cun_total) + Vrp*bit(:,i)*(Cbinm(i)/Cun_total);
        end
    end

    if nstage == 1
        Vresp = Vcm - (Vcp(:,nbits+1)-Vcn(:,nbits+1))*gain / 2;
        Vresn = Vcm + (Vcp(:,nbits+1)-Vcn(:,nbits+1))*gain / 2;
    else
        Vresp = Vcm + (Vcp(:,nbits+1)-Vcn(:,nbits+1))*gain / 2;
        Vresn = Vcm - (Vcp(:,nbits+1)-Vcn(:,nbits+1))*gain / 2;
    end
end
