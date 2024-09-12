function All_results = main_algo(Sys,Sig,Interf,all_series,C)
% This function does three things: interference cancellation, channel
% estimation and symbol detection.
% Last modified: 30-May-2024


% for the 0-th iteration
pilot_symbols = Sig.symbols(1:Sys.data_num+1:end); % extract symbols at pilot positions
[Estimate, rx_ic] = jicce_algo(Sys,Sig,Interf,pilot_symbols,C);

rx_ic(:,1:Sys.data_num+1:end) = []; % remove signal at pilot positions

if isequal(Sys.detect_method{1},'smap')
    x_det = symbol_detection_smap(Sys,rx_ic,Estimate.chan,all_series);
elseif isequal(Sys.detect_method{1},'imap')
    x_det = symbol_detection_imap(Sys,rx_ic,Estimate.chan);
else
    error('Invalid detection method, see properties of Sys.')
end

Estimate.data_demod = Signals.symbol_demod(Sys,x_det);

All_results{1} = Estimate;

if Sys.itr_max == 0
    return
end

old_data = Estimate.data_demod;

% for the 1-st iteration
Estimate = signal_2_symbol(Sys,Sig,Interf,old_data,C);
All_results{2} = Estimate;
if Sys.itr_max == 1
    return
end

% for the next iterations
new_data = Estimate.data_demod;
ii = 2;

while(any(old_data ~= new_data))
    if ii > Sys.itr_max
        break
    end
    old_data = new_data;
    Estimate = signal_2_symbol(Sys,Sig,Interf,old_data,C);
    All_results{ii+1} = Estimate;
    new_data = Estimate.data_demod;
    ii = ii+1;
end
