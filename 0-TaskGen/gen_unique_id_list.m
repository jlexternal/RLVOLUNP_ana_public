% gen_unique_id_list

% Type:         Script
% Level:        0
% Group:        TaskGen
% Usage:        Generates the list of unique IDs to be attributed to participants in
%                   the database
% Study:        RLVOLUNP
%
% Note:         To be run only ONE TIME for study
%
% Jun Seok Lee  - <jlexternal@gmail.com> - September 2021

% creates an array of unique ids characterized by the following format:
%   [enumerator_id][random character of length 4][random numeric array of length 6]
%       where enumerator_id has 3 leading zeros 
%

if isfile('unique_id_keys/unique_id_list.mat') 
    error('Unique ID key list already exists!');
end

char_array  = char(65:90);
num_len     = 6;
char_len    = 4;
unique_id_list = strings(num_len,1);

for isubj = 1:100
    unique_id_list(isubj) = sprintf('%03d%s%06d',isubj,randsample(char_array,char_len),randi(999999));
end

if not(isfolder('unique_id_keys'))
    mkdir unique_id_keys
end

save('unique_id_keys/unique_id_list','unique_id_list');
