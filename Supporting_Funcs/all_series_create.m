function all_series = all_series_create(Sys)
% This function creates all possible series based on points in
% constellation
% all_series has size of Sys.mod_ord^Sys.data_num x Sys.data_num,
% where each row is a possible series of symbols

if strcmp(Sys.mod_type,'QAM')
    all_mod = qammod(0:Sys.mod_ord-1,Sys.mod_ord,'UnitAveragePower',1).';
else
    all_mod = pskmod(0:Sys.mod_ord-1,Sys.mod_ord,Sys.phase_offset).';
end

all_series = zeros(Sys.mod_ord^Sys.data_num,Sys.data_num);
for ii = 1:Sys.data_num
    all_series(:,ii) = repmat(kron((1:Sys.mod_ord)',ones(Sys.mod_ord^(Sys.data_num-ii),1)),Sys.mod_ord^(ii-1),1);
end

all_series = all_mod(all_series);
