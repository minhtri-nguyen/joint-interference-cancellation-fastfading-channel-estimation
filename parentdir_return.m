function par_dir=parentdir_return(current_dir,C)
% return the parent directory path of current_dir
% C: number of upper level
%  i.e., C=0: no change, C=1: 1 upper level, C=2: 2 upper level
% Input must be a valid path


% stop recursive function
if C == 0
    par_dir = current_dir;
    return
end

% change all forwardslash to backslash
current_dir = strrep(current_dir,'/','\');

% remove backslash, if any
if current_dir(end) == '\'
    current_dir = current_dir(1:end-1);
end

ind = find(current_dir=='\',1,'last');
par_dir = current_dir(1:ind-1);

if isempty(find(par_dir=='\',1,'last'))
    error('Cannot find parent directory of a disk.')
end

par_dir = parentdir_return(par_dir,C-1);