function All_results = snr_gap_compute(s_upper,s_lower,snr,ser)
% This function finds the SNR gap to achieve the same SER/BER for s_upper (better)
% and s_lower (worse)
% Input:
%  1. s_upper   : ser/ber array [N x 1] or [1 x N]
%  2. s_lower   : ser/ber array [N x 1] or [1 x N]
%  3. snr       : snr array [N x 1] or [1 x N], in dB
%  4. ser       : scalar, the ser/ber that we want to achieve
% Output:
%   All_results : struct
%    .snr_upper : snr for s_upper to achieve the ser
%    .snr_lower : snr for s_lower to achieve the ser
%    .snr_gap   : snr gap

%%% Sanitize
s_upper = s_upper(:);
s_lower = s_lower(:);
snr = snr(:);

if numel(s_upper) ~= numel(s_lower) || numel(s_upper) ~= numel(snr)
    error('s_upper, s_lower, and snr must have the same size.')
end

s_upper = log(s_upper);
s_lower = log(s_lower);
ser = log(ser);

ind_upper_above = find(s_upper>ser,1,'last');    % SNR above
ind_upper_below = find(s_upper<ser,1,'first');   % SNR below

ind_lower_above = find(s_lower>ser,1,'last');    % SNR above
ind_lower_below = find(s_lower<ser,1,'first');   % SNR below

if isempty(ind_upper_above) || isempty(ind_upper_below) ||...
        isempty(ind_lower_above) || isempty(ind_lower_above)
    error('ser is out of range of s_upper or s_lower.')    
end

%%% Perform the interpolation
snr_upper =(snr(ind_upper_below)-snr(ind_upper_above))/...
    (s_upper(ind_upper_below)-s_upper(ind_upper_above))*...
    (ser-s_upper(ind_upper_above))+snr(ind_upper_above);    % "linear" interpolation, in log scale

snr_lower = (snr(ind_lower_below)-snr(ind_lower_above))/...
    (s_lower(ind_lower_below)-s_lower(ind_lower_above))*...
    (ser-s_lower(ind_lower_above))+snr(ind_lower_above);    % "linear" interpolation, in log scale

snr_gap = snr_upper-snr_lower;

All_results.snr_gap = snr_gap;
All_results.snr_upper = snr_upper;
All_results.snr_lower = snr_lower;
