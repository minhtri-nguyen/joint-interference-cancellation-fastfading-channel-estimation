function x_imap = symbol_detection_imap(Sys,y,h)
% This function does symbol IMAP detection (see the paper)
% Input:
%   1. Sys      : System_Settings object
%       Sys.rx_num      : number of received antenna
%       Sys.chan_alp    : channel correlation coefficient
%       Sys.data_num    : number of data symbols between two pilots
%       Sys.pilot_num   : number of pilot symbols
%   2. y        : [rx_num x data_num(pilot_num-1)] received signal
%   3. h        : [rx_num x pilot_num] estimated channels at pilot symbols
% Output:
%   1. x_imap  : [data_num(pilot_num-1) x 1] soft detected symbols
% Notes:
% - The IMAP method only works for constant power modulated symbols (i.e., it won't work for 16-QAM).
% - When running iterative method, we can use SMAP instead.
% Last modified: 30-May-2024

%%% Sanitize and shorten variable names
a = Sys.chan_alp;
data_num = Sys.data_num;
pilot_num = Sys.pilot_num;
rx_num = Sys.rx_num;

if any(size(y) ~= [rx_num, data_num*(pilot_num-1)])
    error('Wrong size for y.')
end

if any(size(h) ~= [rx_num, pilot_num])
    error('Wrong size for h.')
end

%%% Perform the algorithm
ind_t = 1:data_num;
a_t = a.^ind_t./(1-a.^(ind_t*2));

ind_head = repmat(a_t,rx_num,pilot_num-1);
ind_tail = repmat(a_t(end:-1:1),rx_num,pilot_num-1);

h_head = kron(h(:,1:end-1),ones(1,data_num));
h_tail = kron(h(:,2:end),ones(1,data_num));

H_comp = h_head.*ind_head + h_tail.*ind_tail;
H_comp = H_comp';

x = diag(H_comp*y);
x_imap = x(:);
x_imap = x_imap./abs(x_imap);