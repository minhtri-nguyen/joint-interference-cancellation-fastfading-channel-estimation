classdef System_Settings
    properties
        pilot_num   = 20;       % pilot length
        chan_alp    = 0.99;     % channel correlation coefficient
        tx_num      = 1;        % number of transmit antennas
        rx_num      = 2;        % number of receive antennas
        mod_type    = 'PSK';    % modulation type
        mod_ord     = 16;        % M in M-QAM, takes values of 4,16,64,...
        data_num    = 3;        % number of data symbols between two pilots
        sir         = 1;
        itr_max     = 4;       % maximum number of iterations, for iterative algorithm
        phase_offset = 0;       % for PSK only
        detect_method = {'imap','mrc'}; 
        block_length
        snr
    end
    % Notes on detect_method:
    % - the first one can be 'imap' or 'smap', but 'imap' only works with
    %   constant power modulated symbols
    % - the second one will be used for iterative methods. There are 3
    %   options:
    %   'smap': where estimated channels before and after that symbols are
    %       used
    %   'mrc': where estimated channel at that positions is used
    % Empirically, 'mrc' does not produce good results
    
    methods
        function obj = System_Settings()
            obj = obj.validate();
        end
        
        function obj = validate(obj)
            % This checks errors in system settings, and creates extra fields
            
            if ~strcmp(obj.mod_type, 'QAM') && ~ strcmp(obj.mod_type, 'PSK')
                error('Modulation type must be either QAM or PSK.')
            end
            
            valid_ord = [4 16 64 256 1024 4096];
            if all(valid_ord ~= obj.mod_ord)
                error('Invalid modulation order.')
            end
            
            if obj.chan_alp > 1 || obj.chan_alp <0
                error('Channel correlation coefficient must be in [0,1].')
            end
            obj.block_length = 1 + (obj.data_num+1)*(obj.pilot_num-1);
        end 
    end
end