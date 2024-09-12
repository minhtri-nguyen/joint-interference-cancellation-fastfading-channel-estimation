function [Estimate, rx_ic] = jicce_algo(Sys,Sig,Interf,pilot_symbols,scenario)
% This function does joint interference cancellation and channel estimation (JICCE)
% This works for independent desired channel coefficients only.
% Input:
%   1. Sys          : System_Settings object, or struct
%       Sys.rx_num      : number of receive antennas
%       Sys.chan_alp    : channel correlation coefficient
%       Sys.data_num    : numer of data symbols between two pilots
%       Sys.pilot_num   : number of pilot symbols
%       Sys.snr         : Signal to Noise Ratio
%   2. Sig          : Signals object, or struct
%       Sig.rx_ip       : [rx_num x block_length] received signal with interference
%       Sig.tx_if       : [rx_num x block_length] received signal without interference
%   3. Interf       : Interference object, or struct
%       Interf.interf_num: interference number
%       Interf.mat      : [rx_num x inter_num x block_length] interfering matrix, if any
%   4. pilot_symbols: [pilot_num x 1] pilot symbols
%   5. scenario     : 0 if there is no interference
%                     1 if interference is removed with average estimate
%
% Output:
%   1. Estimate     : struct
%       Estimate.chan       : rx_num x pilot_num, estimated channels
%       Estimate.eic_avg    : block_length x 1, averaged interference coefficients
%       Estimate.eic        : block_length x pilot_num, local estimate of c
%   2. rx_ic        : [rx_num x block_length] received signal with cancelled interference
%
% NOTES:
%   1. Channel power is fixed to 1. Channel coefficients are independent
%       across different Tx-Rx paths.
%   2. Why don't we include pilot_symbols in Sig object? Because later on
%       it will be tedious and error prone to extract pilot symbols from
%       Sig when running iterative algorithm.
% Last modified: 30-May-2024

%%% Sanitize and shorten variable names
if all([0 1] ~= scenario)
    error('no supported configuration C')
end

N = Sys.pilot_num;
a_p = Sys.chan_alp^(Sys.data_num+1);
Nr = Sys.rx_num;

snr = Sys.snr;
s_n = 1/snr;

%%% Extract necessary variables
if scenario == 1
    y = Sig.rx_ip(:,1:Sys.data_num+1:end); % 1 x pilot_num, rx signal at pilot positions
else
    y = Sig.rx_if(:,1:Sys.data_num+1:end); % 1 x pilot_num, rx signal at pilot positions
end
if ~isequal(size(y),[Sys.rx_num  Sys.pilot_num])
    error('Wrong size')
end

interf_mat = Interf.mat(:,:,1:Sys.data_num+1:end);
if ~isequal(size(interf_mat),[Sys.rx_num  Interf.interf_num  N]) && scenario~=0
    error('Wrong size')
end

%%% pre-allocate matrices
A_nn = zeros(1,N);    % A_n in the paper, with independent channel assumption
d_xSy = zeros(Nr,N);  % the term that goes after A_n in estimate formula of h_n

if scenario == 1
    L =Interf.interf_num;
    eic = zeros(L,N);
    d_xSB = zeros(Nr,L,N);
end

%%% Perform the algorithm
all_ind = 1:N;
for nn = 1:N
    j_n = sign(nn - all_ind);
    temp1_n = snr*(1 - a_p.^(2*(abs(nn - all_ind) - 1)));
    
    x_n = a_p.^(abs(nn - all_ind))./(1 + temp1_n).*pilot_symbols;
    x_n(nn) = pilot_symbols(nn);
    
    beta_n = pilot_symbols.*conj(pilot_symbols(all_ind + j_n))*a_p.*temp1_n./(1 + temp1_n);
    beta_n(nn) = 0;
    
    y_n = y - y(:,all_ind + j_n).*repmat(beta_n,Nr,1);
    S_n = s_n + 1 - a_p.^(2*(abs(nn - all_ind))) - s_n*a_p^2.*temp1_n.^2./(1 + temp1_n);
    S_n(nn) = s_n;
    
    if any(S_n <= 0)
        fprintf('%d\n',nn)
        error('Zero or negative variance.')
    end
    
    %%% for interference matrices
    if scenario == 1
        B_n = interf_mat-interf_mat(:,:,all_ind + j_n).*repmat(permute(beta_n,[3 1 2]),Nr,L,1);
    end
    
    %%% compute matrices
    a_n = 1+sum(abs(x_n).^2./S_n);
    A_nn(nn) = a_n;
    
    d_xSy(:,nn) = sum(y_n.*repmat(conj(x_n)./S_n,Nr,1),2);
    
    if scenario == 1
        D = zeros(L);
        d_BSy = zeros(L,1);
        
        for ii = 1:N
            D = D+B_n(:,:,ii)'*B_n(:,:,ii)/S_n(ii);
            d_xSB(:,:,nn) = d_xSB(:,:,nn)+x_n(ii)'/S_n(ii)*B_n(:,:,ii);
            d_BSy = d_BSy+B_n(:,:,ii)'/S_n(ii)*y_n(:,ii);
        end
        
        D = D-d_xSB(:,:,nn)'*d_xSB(:,:,nn)/a_n;
        
        if all(D == 0)
            eic(:,nn) = zeros(Sys.L,1);
        else
            eic(:,nn) = D\(d_BSy-d_xSB(:,:,nn)'*d_xSy(:,nn)/a_n);
        end
    end
end

%%% compute estimated channel coefficients and eic
if scenario == 0
    % no inteference
    h = zeros(Nr,N);
    for nn = 1:N
        h(:,nn) = d_xSy(:,nn)/A_nn(nn);
    end
    rx_ic = Sig.rx_if;
else
    eic_avg = mean(eic,2);
    Estimate.eic_avg = eic_avg;
    Estimate.eic_all = eic;
    
    % interference was removed by average estimate
    h = zeros(Nr,N);
    for nn = 1:N
        h(:,nn) = (d_xSy(:,nn)-d_xSB(:,:,nn)*eic_avg)/A_nn(nn);
    end
    
    rx_ic = zeros(Sys.rx_num,Sys.block_length);
    for nn = 1:Sys.block_length
        rx_ic(:,nn) = Sig.rx_ip(:,nn) - Interf.mat(:,:,nn)*eic_avg;
    end 
end

Estimate.chan = h;
