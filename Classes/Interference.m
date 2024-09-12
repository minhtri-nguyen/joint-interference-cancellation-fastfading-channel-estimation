classdef Interference
    properties (Constant)
        %%% Parameters for desired system
        S_desired = struct('beta',     0.25,...
            'span',     4,...
            'sps',      10,...
            'bandwidth',15e3);
        %%% Parameters for interfering system
        S_interf = struct('beta',     0.25,...
            'span',     4,...
            'sps',      20,...
            'bandwidth',3e4);
        mod_type = 'QAM';
        mod_ord = 4;
    end
    
    properties
        chan
        eic
        interf_num
        sig
        mat
    end
    
    methods
        function obj = Interference(Sys)
            obj = obj.eic_create(Sys);
        end
        
        function obj = eic_create(obj,Sys)
            Ts1 = (1+obj.S_desired.beta)/obj.S_desired.bandwidth;
            df = 0.1/Ts1;
            obj.eic = Interference.xspectrum(obj.S_desired,obj.S_interf,df);
            
            obj.chan = 1/Sys.sir*ones(Sys.rx_num,1); % Sys.rx_num x 1
            obj.interf_num = length(obj.eic);
        end
        
        function obj = interference_create(obj,Sys)
            % interference symbols
            i_int = randi([0 obj.mod_ord-1], [Sys.block_length, obj.interf_num]);
            if strcmp(obj.mod_type, 'QAM')
                t_symbol = qammod(i_int,obj.mod_ord,'UnitAveragePower',1);     % block_length x interf_num
            elseif strcmp(obj.mod_type, 'PSK')
                t_symbol = paskmod(i_int,obj.mod_ord,'UnitAveragePower',1);     % block_length x interf_num
            else
                error('Invalid modulation type.')
            end
            obj.sig = zeros(Sys.rx_num,Sys.block_length);   % interference signal, to add to Rx signal
            obj.mat = zeros(Sys.rx_num,obj.interf_num,Sys.block_length);    % interference matrix B_n in the paper
            for nn = 1:Sys.block_length
                obj.mat(:,:,nn)=obj.chan*t_symbol(nn,:);
                obj.sig(:,nn) = obj.mat(:,:,nn)*obj.eic;
            end
        end
    end
    
    methods (Static)
        function eic = xspectrum(S_desired,S_interf,df)
            % generate normalized eic
            RBW_12 = S_interf.bandwidth/S_desired.bandwidth;
            
            eic = zeros(RBW_12,1);
            
            for tt = 1:RBW_12
                eic(tt) = Interference.channel_coeff_eff(S_desired,S_interf,df,tt-(RBW_12-1)/2);
            end
            eic = eic/norm(eic);
        end
        
        function c = channel_coeff_eff(S_desired,S_interf,df,dt)
            % S_desired and S_interf have fields
            % beta, span, sps, bandwidth
            
            % find effective channel coefficients
            b1 = S_desired.beta;
            b2 = S_interf.beta;
            
            BW1 = S_desired.bandwidth;
            BW2 = S_interf.bandwidth;
            RBW_12 = BW2/BW1;
            
            sps2 = S_interf.sps;
            sps1 = sps2*RBW_12;
            
            span1 = S_desired.span;
            span2 = S_interf.span;
            
            fl1 = rcosdesign(b1,span1,sps1);
            st = (1+b1)/BW1;
            t = (0:length(fl1)-1)*st/sps1;
            
            fl2 = rcosdesign(b2,span2,sps2);
            
            dL = length(fl1)-length(fl2);
            
            fl2 = [zeros(1,dL/2) fl2 zeros(1,dL/2)];
            fl2 = circshift(fl2,[0 dt*sps2]);
            
            fl2 = fl2.*exp(1j*2*pi*df*t);
            
            
            ff = conv(fl1,fl2);
            c = ff((end-1)/2);
        end
    end
end