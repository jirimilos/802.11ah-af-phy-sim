function [ ] = plotIQdiag( dataMapped, modOrder )
% This function plot constellation diagram.

dataMapped = reshape(dataMapped, 1, []);

scatterplot(dataMapped)

if modOrder == 2
    %title('BPSK, Symbol Mapping')
elseif modOrder == 4
    %title('QPSK, Symbol Mapping')
else
    %title([num2str(modOrder) 'QAM, Symbol Mapping'] )
end
axis([-2 2 -2 2])

end

