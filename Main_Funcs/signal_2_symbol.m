function Estimate = signal_2_symbol(Sys,Sig,Interf,data_demod,C)
% This function does three things: interference cancellation, channel
% estimation and symbol detection, iteratively.
% Last modified: 30-May-2024


data_ind = 1:Sys.block_length;
data_ind(1:Sys.data_num+1:end)=[];

%%% Perform jicce. We consider all symbols are pilot symbols
Sys_itr = Sys;
Sys_itr.data_num = 0;
Sys_itr.pilot_num = Sys.block_length;

pilot_symbols = Sig.symbols;         % get symbols of the whole block

pilot_symbols(data_ind) = Signals.data_mod(Sys,data_demod);   % replace symbols in data positions by their estimate
[Estimate, rx_ic] = jicce_algo(Sys_itr,Sig,Interf,pilot_symbols,C);

%%% Perform detection
det_method = Sys.detect_method{2};

if isequal(det_method,'imap') || isequal(det_method,'smap')
    % In this case, we use estimated channels before and after every symbol
    % for detection (imap and smap produces same result).
    
    % we need to take care of data_num and pilot_num
    Sys_itr.data_num = 1;
    x_det = zeros(Sys.block_length,1);
    all_series = Signals.data_mod(Sys,[0 Sys.mod_ord]');
    
    % for odd pilot positions
    ind_odd = 1:2:Sys.block_length;
    Sys_itr.pilot_num = length(ind_odd);
    x_det(ind_odd(2:end)-1) = symbol_detection_smap(Sys_itr,rx_ic(:,ind_odd(2:end)-1),Estimate.chan(:,ind_odd),all_series);
    
    % for even positions
    ind_even = 2:2:Sys.block_length;
    Sys_itr.pilot_num = length(ind_even);
    x_det(ind_even(2:end)-1) = symbol_detection_smap(Sys_itr,rx_ic(:,ind_even(2:end)-1),Estimate.chan(:,ind_even),all_series);
    
    
elseif isequal(det_method, 'mrc')
    % In this case, we use estimated channel at every data positions for
    % the detection
    x_det = diag(Estimate.chan'*rx_ic);
    x_det = x_det.'./sum(abs(Estimate.chan).^2,1);  % no need to normalize
else
    error('Invalid detection method in Sys.detect_method{2}.')
end

%%% Post process
x_det(1:Sys.data_num+1:end)=[];
x_det = x_det(:);
x_demod = Signals.symbol_demod(Sys,x_det);

if numel(x_demod) ~= Sys.data_num*(Sys.pilot_num-1)
    error('Wrong size')
end
Estimate.data_demod = x_demod;
Estimate.chan = Estimate.chan(:,1:Sys.data_num+1:end);