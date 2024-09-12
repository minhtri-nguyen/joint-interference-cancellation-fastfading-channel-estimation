classdef Signals
    properties 
        noise_uni
        chan
        rx_noiseless_if
        rx_noiseless_ip
        rx_if
        rx_ip
        data_data  % we only need data at data positions, for error rate computation
        symbols % symbols for the whole block
    end
    
    methods
        function obj = Signals(Sys)
            obj = obj.tx_process(Sys);
        end
        
        function obj = tx_process(obj,Sys)
            data_ = randi([0 Sys.mod_ord-1],[1 Sys.block_length]); % data of the whole block
            obj.symbols = Signals.data_mod(Sys,data_);
            data_(1:Sys.data_num+1:end) = [];   % remove data at pilot positions
            obj.data_data = data_;              % data at data positions
            
            % unit variance noise
            obj.noise_uni = (randn(Sys.rx_num,Sys.block_length)+1j*randn(Sys.rx_num,Sys.block_length))/sqrt(2);
        end
        
        function obj = channel_effect(obj,Sys,Interf)
            % this function add channel effects and interference
            
            % channel
            chan_ = fastfading_gen(Sys);
            obj.chan = squeeze(permute(chan_,[1 3 2]));
    
            % signal with/without interference
            obj.rx_noiseless_if = obj.chan.*repmat(obj.symbols,Sys.rx_num,1);
            obj.rx_noiseless_ip = obj.rx_noiseless_if + Interf.sig;
        end
        
        function obj = add_noise(obj,Sys)
            s_n = 1/Sys.snr;

            obj.rx_if = obj.rx_noiseless_if + sqrt(s_n)*obj.noise_uni;
            obj.rx_ip = obj.rx_noiseless_ip + sqrt(s_n)*obj.noise_uni;
        end
    end
    
    methods (Static)
        function x_mod = data_mod(Sys,data)
            if strcmp(Sys.mod_type, 'QAM')
                x_mod = qammod(data,Sys.mod_ord,'UnitAveragePower',1); % QAM modulated symbols
            elseif strcmp(Sys.mod_type, 'PSK')
                x_mod = pskmod(data,Sys.mod_ord,Sys.phase_offset);     % PSK modulated symbols
            else
                error('Invalid modulation type.')
            end
        end
        
        function x_demod = symbol_demod(Sys,signal)
            if strcmp(Sys.mod_type,'QAM')
                x_demod = qamdemod(signal,Sys.mod_ord,'UnitAveragePower',1); % integers
            elseif strcmp(Sys.mod_type,'PSK')
                x_demod = pskdemod(signal,Sys.mod_ord,Sys.phase_offset);    % intergers
            else
                error('Invalid modulation type.')
            end
        end
    end
end
