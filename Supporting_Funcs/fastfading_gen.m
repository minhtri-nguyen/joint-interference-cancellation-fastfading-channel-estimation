function h = fastfading_gen(Sys)
% This function generates fast fading channel having 1st order Markov
% property: h_(k+1)=a*h(k)+sqrt(1-a^2)*noise
%
% Input:
%   1. Sys  : struct, or System_Setting object
%           Sys.chan_alp        : channel correlation coefficient, 0<a<1
%           Sys.block_length    : length of channel
%           Sys.rx_num          : number of receive antennas
%           Sys.tx_num          : number of transmit antennas
% Output:
%   1.h     : [rx_num x tx_num x block_length] channel matrix
% NOTE:
%   1. all channel coefficients are complex
%   2. h(1), initial value, is Gaussian
% Last modified 30-May-2024

a = Sys.chan_alp;
rx = Sys.rx_num;
tx = Sys.tx_num;
N = Sys.block_length;

if a > 1 || a < 0
    error('Channel correlation coefficient must be in [0,1]')
end

% generate noise
noise = sqrt(1-a^2)*sqrt(1/2)*(randn(rx,tx,N-1)+ 1j*randn(rx,tx,N-1));

% initial channel gain
h(:,:,1) = sqrt(1/2)*(randn(rx,tx)+1j*randn(rx,tx));

for ii = 2:Sys.block_length
    h(:,:,ii) = a*h(:,:,ii-1)+noise(:,:,ii-1);
end
