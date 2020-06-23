function [ secPermutInt ] = longRunsIntrlv( data, Nbpsc, Ncols )
%LONGRUNSINTRLV - This function provide second permutation due to long runs
%of LSB bits.

    Ncbps = length(data);
    j = zeros(1, Ncbps);
    secPermutInt = zeros(1, Ncbps);
    s = max(Nbpsc/2, 1);

    for i = 1:length(data)
        j(i) = s * fix((i-1) / s) + mod(((i-1) + Ncbps - fix((Ncols*(i-1))/Ncbps)), s);    % to calculate the new bit position
        secPermutInt(j(i)+1) = data(i);                                                    % move the bit to a new position
    end

end

