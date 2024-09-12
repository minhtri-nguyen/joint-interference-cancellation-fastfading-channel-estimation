function x_smap = symbol_detection_smap(Sys,y,h,all_series)
% This function does symbol SMAP detection (see the paper)
% Input:
%   1. Sys      : System_Settings object
%       Sys.rx_num      : number of received antenna
%       Sys.chan_alp    : channel correlation coefficient
%       Sys.data_num    : number of data symbols between two pilots
%       Sys.pilot_num   : number of pilot symbols
%       Sys.mod_ord     : modulation order (4 for QPSK, 4-QAM, etc)
%   2. y        : [rx_num x data_num(pilot_num-1)] received signal
%   3. h        : [rx_num x pilot_num] estimated channels at pilot symbols
%   4. all_series : [mod_ord^data_num x data_num] all possible symbol series
%  Output:
%   1. x_smap   : [data_num(pilot_num-1) x 1] hard-detected symbols
% Last modified: 30-May-2024

if any(size(y) ~= [Sys.rx_num, Sys.data_num*(Sys.pilot_num-1)])
    error('Argument y has wrong size.')
end

if any(size(h) ~= [Sys.rx_num, Sys.pilot_num])
    error('Argument h has wrong size.')
end

if ~isequal(size(all_series),[Sys.mod_ord^Sys.data_num Sys.data_num])
    error('Wrong combinations')
end

x = zeros(Sys.data_num,Sys.pilot_num-1);
a = Sys.chan_alp;

s_n = 1/Sys.snr;
t2 = a/(1-a^2);

S_k = zeros(1,Sys.data_num+1);

% S_k(2) is S_k(1) in the paper, this is for convinient programming

for ii = 1:Sys.data_num
    S_k(ii+1) = 1/(1/s_n+(1+a^2)/(1-a^2)-t2^2*S_k(ii));
end

S_k = S_k(2:end);

gamma_k = ones(Sys.data_num);

for ii = 1:Sys.data_num
    gamma_k(ii,ii:-1:1) = t2.^(0:ii-1).*cumprod([1 S_k(ii-1:-1:1)]);
end

for pp = 1:Sys.pilot_num-1
    
    y_k = y(:,(pp-1)*Sys.data_num+1:pp*Sys.data_num);
    h_h = h(:,pp)*t2;
    h_t = h(:,pp+1)*t2;

    h_h_rep = repmat(h_h,1,Sys.data_num).*repmat(gamma_k(:,1)',Sys.rx_num,1);
    
    h_t_rep = [zeros(Sys.rx_num,Sys.data_num-1) h_t];
    
    F = zeros(Sys.mod_ord^Sys.data_num,1);
    for kk = 1:Sys.mod_ord^Sys.data_num
        y_all = zeros(Sys.rx_num,Sys.data_num);
        for ii = 1:Sys.data_num
            y_all(:,ii) = sum(repmat(all_series(kk,1:ii),Sys.rx_num,1).*y_k(:,1:ii).*gamma_k(ii,1:ii),2)/s_n;
        end
        
        F(kk) = sum(diag((y_all + h_h_rep + h_t_rep)'*(y_all + h_h_rep + h_t_rep)).'./S_k);
    end
    [~,ind_max] = max(F);
    x(:,pp) = all_series(ind_max,:).';
end

x_smap = conj(x(:));
    