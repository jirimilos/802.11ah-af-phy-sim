function [ lfsrOutD ] = descrambler( data )
%DESCRAMBLER - Decrambling of the data

    numberOfBits = length(data);

    initSequnece = [1 0 1 1 1 0 1];                    % initial scrambler sequence
    
    lfsrShiftD = zeros(1, 7);
    lfsrOutD = zeros(1, numberOfBits);
    
    for i = 1 : numberOfBits
        if i < 8
            lfsrOutD(i) = xor(data(i), initSequnece(i));
            lfsrShiftD(2:end) = lfsrShiftD(1:end-1);   % shift LFSR
            lfsrShiftD(1) = initSequnece(i);
        else
            lfsrXOR = xor(lfsrShiftD(4), lfsrShiftD(7));
            lfsrOutD(i) = xor(data(i), lfsrXOR);       % xor x7 + x4 + 1
            lfsrShiftD(2:end) = lfsrShiftD(1:end-1);   % shift LFSR
            lfsrShiftD(1) = lfsrXOR;     
        end
    end


end

