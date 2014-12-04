% Burrows Wheeler Transformation Analysis
% implemented by Dan Belling

function burrow-wheeler(varargin)

% Part 1 - Burrows-Wheeler Transform
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[desc1, toySeq] = fastaread('sample1.fa');
time_start = cputime;
bwtSeq = bwt(toySeq);
fastawrite('bwtSample1.fa',bwtSeq);
time_end = cputime - time_start;
fprintf('Wrote: %s BWT in %d seconds\n', desc1, time_end);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Part 2 - Inverse Burrows-Wheeler Transform
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[desc2, sapienSeq] = fastaread('bwtHomoSapiens.fa');
time_start = cputime;
seq = reverseBwt(sapienSeq);
fastawrite('rBwtSample2.fa',seq);
time_end = cputime - time_start;
fprintf('Wrote: %s reverse BWT in %d seconds\n', desc2, time_end);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end

% Function - Constructs a Burrows-Wheeler sequence from a given sequence
function bwtSeq = bwt(seq)

% Initialize the appended sequence, declare relevant matrices %
appendSeq = strcat(seq, '$'); % Sequence with special appended char %
len = length(appendSeq); % Length of appended sequence %
bwt_mat = zeros(1,len); % Output of the transform (see II) %
perm = cell(1,len); % Cell array of cyclic permutations (see I)%

% I. Cyclic permutations
% Form cell array of cyclic permutations (See utility function) %
% Alphabetically sort the cell matrix %
[bwCell_sort, suff_array] = suf_con(appendSeq);

% II. BW Transform
% Initialize the transform, and remove the last column %
for j = 1:len
    subSeq = bwCell_sort{j};
    bwt_mat(1,j) = subSeq(len:len);
end
% Cast the integer matrix back to a character array %
bwtSeq = char(bwt_mat);

end

% Function - reverses a Burrows-Wheeler transformation
function seq = reverseBwt(bwtSeq)
% Initialize the cell count vector and tally matrix %

% Declarations - START %
bwtSeq = upper(bwtSeq); % Ensure capitalization of ill-formatted chrY seq %
len = length(bwtSeq); % Length of BWTseq %
outSeq = cell(1,len); % Cell vector for storing the output sequence (see IV) %
index = 1; % Index counter %
% Declarations - END %

% I. Base Count (Utility function) %
% Format is 1: $, 2: A, 3: C, 4: G, 5: T %
% Each slot corresponds to the number of characters less than that char %
count_vec = bse_cnt(bwtSeq);

% II. Occurence Matrix construction (Utility function) %
iCh = bwtSeq(1);
occ_mat = occurence_mat(iCh, bwtSeq);

% III. LF - Mapping (Utility function)
lf_map = last_frnt_map(bwtSeq, occ_mat, count_vec);

% IV. Output sequence construction %
% Reconstruct the original string %
for k = 1:len
    outSeq{k} = bwtSeq(index);
    index = lf_map(1, index);
end
% Flip reverse constructed sequence %
outSeq = flip(outSeq);
seq = cell2mat(outSeq);
end


end


end

% Utility functions %
function occ_mat = occurence_mat(iCh, bwtSeq)
len = length(bwtSeq); % Length of the sequence %
occ_mat = zeros(len,5); % Occurence matrix of base ranks (see II) %
% Initialize occurence matrix %
switch iCh
    case '$'
        occ_mat(1,1) = 1;
    case 'A'
        occ_mat(1,2) = 1;
    case 'C'
        occ_mat(1,3) = 1;
    case 'G'
        occ_mat(1,4) = 1;
    case 'T'
        occ_mat(1,5) = 1;
end
% Populate the occurence matrix %
for i = 2:len
    ch = bwtSeq(i);
    % Add the previous row %
    occ_mat(i,:) = occ_mat(i-1,:);
    switch ch
        case '$'
            occ_mat(i,1) = occ_mat(i,1) + 1;
        case 'A'
            occ_mat(i,2) = occ_mat(i,2) + 1;
        case 'C'
            occ_mat(i,3) = occ_mat(i,3) + 1;
        case 'G'
            occ_mat(i,4) = occ_mat(i,4) + 1;
        case 'T'
            occ_mat(i,5) = occ_mat(i,5) + 1;
    end
end
end

% This function constructs a last-front mapping vector for use in
% reverseing a BWT sequence 
function lf_map = last_frnt_map(bwtSeq, occ_mat, count_vec)
len = length(bwtSeq); % Length of the sequence %
lf_map = zeros(1,len); % Mapping vector for use in final sequence construction (see III) %
% Begin the Last-Front mapping %
for j = 1:len
    iCh = bwtSeq(j);
    switch iCh
        case 'A'
            lf_map(1,j) = occ_mat(j,2) + count_vec(1,2);
        case 'C'
            lf_map(1,j) = occ_mat(j,3) + count_vec(1,3);
        case 'G'
            lf_map(1,j) = occ_mat(j,4) + count_vec(1,4);
        case 'T'
            lf_map(1,j) = occ_mat(j,5) + count_vec(1,5);
        case '$'
            lf_map(1,j) = 1;
    end
end
end

% This function constructs a base count vector which gives the number of
% residues which register as 'alphabetically/lexically smaller' than the
% previous base
function count_vec = bse_cnt(bwtSeq)
bases = basecount(bwtSeq); % Caclculate number of each nucleotide in seq (disregard warning) %
count_vec = zeros(1,5); % Vector of base count references (see I) %
count_vec(1,1) = 0;
count_vec(1,2) = 1;
count_vec(1,3) = count_vec(1,2) + bases.A;
count_vec(1,4) = count_vec(1,3) + bases.C;
count_vec(1,5) = count_vec(1,4) + bases.G;
end

% This function constructs a rank vector which corresponds to residue
% positions within a BWT sequence
function rank_mat = rank_con(bwtSeq, occ_mat)
len = length(bwtSeq); % Length of sequence %
rank_mat = zeros(1,len); % Nucleotide rank vector %
% Fill in the rank vector %
for j = 1:len
    ch = bwtSeq(j);
    switch ch
        case 'A'
            rank_mat(1,j) = occ_mat(j,2);
        case 'C'
            rank_mat(1,j) = occ_mat(j,3);
        case 'G'
            rank_mat(1,j) = occ_mat(j,4);
        case 'T'
            rank_mat(1,j) = occ_mat(j,5);
        case '$'
            rank_mat(1,j) = 1;
    end
end
end

% This function forms and sorts an array of cyclic permutations for use
% within the BWT transform function, also produces a suffix array.
function [bwCell_sort, suff_array] = suf_con(seq)
    len = length(seq); % Length of sequence %
    perm = cell(1,len); % Cell array of cyclic permutations %
    
    % Form cyclic permutations %
    for i = 1:len
        perm{i} = strcat(seq(i:len), seq(1:i-1));
    end

    [bwCell_sort, suff_array] = sort(perm);

end

% This function performs a search given an entry query in a bwtSeq
function [first, last] = query_seq(qry, bwtSeq, count_vec, lf_map)

    % Declarations - START %
    len = length(bwtSeq); % Length of BWT seq %
    qry_len = length(qry); % Length of query sequence %
    first = 1; % Start of query range %
    last = len; % End of query range %
    index = qry_len; % Index counter into the query %
    count_ref = 0; % count_vec reference %
    % Declarations - END %
    
    % Initialize first residue of search %
    ch = qry(index);
        switch ch
            case 'A'
                first = 2;
                last = count_vec(1,3);
            case 'C'
                first = count_vec(1,3) + 1;
                last = count_vec(1,4);
            case 'G'
                first = count_vec(1,4) + 1;
                last = count_vec(1,5);
            case 'T'
                first = count_vec(1,5) + 1;
                last = len;
        end
        
        index = index - 1; % Update index counter %

    while first <= last && index >= 2
         
        first = lf_map(1,first);
        last = lf_map(1,last);
        
        index = index - 1;
    end
end
