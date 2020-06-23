rolloff = 0.25 % betta -> rolloff factor
span = 6; % filter span in symbols
sps = 4; % samples per symbol

% square-toor, raised cosine filter coeffs
b = rcosdesign(rolloff, span, sps); 

%fvtool(b,'impulse')

% vector of bipolar drata
d = 2*randi([0 1],100,1) -1;

% upsample and filter the data for pulse shaping
x = upfirdn(d, b, sps);

%add noise
r = x + randn(size(x))*0.01;

%filter and downsample
y = upfirdn(r,b,1,sps);

subplot(3,1,1)
plot(d)
subplot(3,1,2)
plot(x)
subplot(3,1,3)
plot(y)