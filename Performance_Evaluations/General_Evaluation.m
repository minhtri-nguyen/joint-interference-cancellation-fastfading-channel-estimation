%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial clear %%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_path = mfilename('fullpath');
ind_dump = find(file_path == '\',1,'last');

file_path = file_path(1:ind_dump-1);
addpath([parentdir_return(file_path,2) '\Classes']);
addpath([parentdir_return(file_path,2) '\Supporting_Funcs']);
addpath([parentdir_return(file_path,2) '\Main_Funcs']);

%%%%%%%%%%%%%%%%%%%%%%% Simulation settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR_dB = 0:10:40;       % SNR, in dB
SNR = db2pow(SNR_dB);   % SNR, in scalar

SIR_dB = 0;             % SIR, in dB
SIR = db2pow(SIR_dB);   % SIR, in scalar
N_test = 1e3;           % number of tests
N_step = 100;

%%%%%%%%%%%%%%%%%%%%%%% Model setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sys = System_Settings();    % create System_Settings object with default properties
Interf = Interference(Sys); % create Interference object with default properties

%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmse_ip = zeros(length(SNR),Sys.pilot_num);
cmse_if = zeros(length(SNR),Sys.pilot_num);
err_if0 = zeros(length(SNR),1);
err_ip0 = zeros(length(SNR),1);
err_if1 = err_if0;
err_ip1 = err_ip0;

n_data = Sys.data_num*(Sys.pilot_num-1);
all_series = all_series_create(Sys);
L_avg = zeros(numel(SNR_dB),1);

tic
fprintf('Start simulation\n')
for tt = 1:N_test
    
    Sig = Signals(Sys); % generate signals at Tx
    Interf = Interf.interference_create(Sys);       % generate interference
    Sig = Sig.channel_effect(Sys,Interf);           % add channel effects and interference
    chan_pilot = Sig.chan(:,1:Sys.data_num+1:end);  % channel at pilot positions, for CMSE computation
    
    % loop for SNR
    for ss = 1:length(SNR)
        Sys.snr = SNR(ss);
        
        Sig = Sig.add_noise(Sys);   % add noise
        Sys.detect_method{1} = 'imap';
        Allout_ip = main_algo(Sys,Sig,Interf,all_series,1);
        cmse_ip(ss,:) = cmse_ip(ss,:)+sum(abs(Allout_ip{1}.chan-chan_pilot).^2,1)/Sys.rx_num;
        x_demod_ip = Allout_ip{1}.data_demod;
        err_ip0(ss) = err_ip0(ss) + sum(x_demod_ip~=Sig.data_data')/n_data/N_test;
        
        
%         Allout_if = main_algo(Sys,Sig,Interf,all_series,0);
%         cmse_if(ss,:) = cmse_if(ss,:) + sum(abs(Allout_ip{2}.chan-chan_pilot).^2,1)/Sys.rx_num;
        x_demod_if = Allout_ip{2}.data_demod;
        err_if0(ss) = err_if0(ss) + sum(x_demod_if~=Sig.data_data')/n_data/N_test;
        
%         Sys.detect_method{1} = 'smap';
%         Allout_ip = main_algo(Sys,Sig,Interf,all_series,1);
%         cmse_ip(ss,:) = cmse_ip(ss,:)+sum(abs(Allout_ip{1}.chan-chan_pilot).^2,1)/Sys.rx_num;
        x_demod_ip = Allout_ip{end}.data_demod;
        err_ip1(ss) = err_ip1(ss) + sum(x_demod_ip~=Sig.data_data')/n_data/N_test;
        
        L_avg(ss) = L_avg(ss) + length(Allout_ip)/N_test;
    end
    
    if mod(tt,N_step) == 0
        fprintf('Done for %d iterations\n',tt)
    end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%% Post process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmse_ip = sum(cmse_ip,2)/N_test/Sys.pilot_num;
cmse_if = sum(cmse_if,2)/N_test/Sys.pilot_num;

figure
semilogy(SNR_dB,cmse_if)
hold on
semilogy(SNR_dB,cmse_ip)
legend

figure
semilogy(SNR_dB,err_ip0)
hold on
semilogy(SNR_dB,err_if0)
semilogy(SNR_dB,err_ip1)
semilogy(SNR_dB,err_if1)
legend