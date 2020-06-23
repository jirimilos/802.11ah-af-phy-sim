function [ secPermutDeint ] = longRunsDeintrlv( data, Nbpsc, Ncols )

%LONGRUNSDEINTRLV - This function provide second permutation due to long runs
%of LSB bits.

    Ncbps = length(data);
    jj = zeros(1, Ncbps);
    secPermutDeint = zeros(1, Ncbps);
    s = max(Nbpsc/2, 1);

    for i = 1:length(data)
        jj(i) = s * fix((i-1) / s) + mod(((i-1) + fix((Ncols*(i-1))/Ncbps)), s);    % to calculate the new bit position
        secPermutDeint(jj(i)+1) = data(i);                                          % move the bit to a new position
    end

end

