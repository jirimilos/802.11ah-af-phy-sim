function [ lfsrOutS ] = scrambler( data )
%SCRAMBLER - Scrambling of the input data

    numberOfBits = length(data);
    
    initSequence = [1 0 1 1 1 0 1];                     % initial scrambler sequence

    lfsrShiftS = zeros(1, 7);
    lfsrOutS = zeros(1, numberOfBits);
    
    for i = 1 : numberOfBits
        if i < 8
            lfsrOutS(i) = xor(data(i), initSequence(i));
            lfsrShiftS(2:end) = lfsrShiftS(1:end-1);    % shift LFSR
            lfsrShiftS(1) = initSequence(i);
        else
            lfsrXOR = xor(lfsrShiftS(4), lfsrShiftS(7));
            lfsrOutS(i) = xor(data(i), lfsrXOR);        % xor x7 + x4 + 1
            lfsrShiftS(2:end) = lfsrShiftS(1:end-1);    % shift LFSR
            lfsrShiftS(1) = lfsrXOR;     
        end
    end
      
end

