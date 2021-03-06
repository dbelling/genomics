function sankoff-tree(varargin)

% Phylogenetic tree construction
% Implemented by: Dan Belling

% Open the sequence files and store the struct contents in an array %
listing = dir('./sequences');
seq_struct = struct('Header', ' ', 'Sequence', ' ' );
fprintf('Opening files in sequences directory..\n');
for i = 3:40
    file_name = listing(i).name;
    current_file = fastaread(strcat('./sequences/',file_name));
    seq_struct(i-2).Header = current_file.Header;
    seq_struct(i-2).Sequence = current_file.Sequence;
end

% ~~~~~~~~~~~~~~~~ PART I ~~~~~~~~~~~~~~~~~~~~~~~~ %
% Perform multi-sequence alignment with multialign %
fprintf('Performing multi-sequence align..\n');
seq_align_struct = multialign(seq_struct);
fprintf('Writing multi sequence alignment results..\n');
fID_one = fopen('multi_align.txt', 'wt');

for j = 1:38
    nbytes = fprintf(fID_one, '%s \t %s\n', seq_align_struct(j).Sequence, seq_align_struct(j).Header);
    if nbytes <= 0
        error('Could not write to sequence file');
    end
end
fclose(fID_one);

% ~~~~~~~~ PART II ~~~~~~~~~ %
% Compute pairwise distances %
fprintf('Computing pairwise distances..\n');
dist_vec = seqpdist(seq_align_struct); 
fID_two = fopen('pair_dist.txt', 'wt');

for k = 1:38
    nbytes = fprintf(fID_two, '%d \t %s\n', dist_vec(k), seq_align_struct(k).Header);
    if nbytes <= 0
        error('Could not write to distance file');
    end
end
fclose(fID_two);

% Tree construction & neighbor joining %
fprintf('Constructing tree..\n');
phylo_tree = seqneighjoin(dist_vec, 'equivar', seq_align_struct);
h = plot(phylo_tree, 'orient', 'top');
ylabel('Evolutionary distance');
set(h.terminalNodeLabels, 'Rotation', 65);

% ~~~~~~~~~~~~~~~~ PART III ~~~~~~~~~~~~~~~~ %
% Hard-coded scoring matrix for toy examples %
score_mat = [0, 3, 4, 9, 8;
            3, 0, 2, 4, 8;
            4, 2, 0, 4, 8;
            9, 4, 4, 0, 8;
            8, 8, 8, 8, 8];

[sankoff_score, pointer_struct] = sankoff(phylo_tree, score_mat, seq_align_struct);
        
% ~~~~~~~~~~~~~~~~ Part IV ~~~~~~~~~~~~~~~~~ %
% Parsimony analysis of ribosomal RNA sequences (handled from part 3) %

end

function [sank_score, pointer_struct] = sankoff(phylo_tree, score_mat, seq_align_struct)
    
    % Leaf + root score declarations
    score_vec = zeros(1,5);
    A_leaf = [0, Inf, Inf, Inf, Inf]; % A leaf score
    U_leaf = [Inf, 0, Inf, Inf, Inf]; % U leaf score
    G_leaf = [Inf, Inf, 0, Inf, Inf]; % G leaf score
    C_leaf = [Inf, Inf, Inf, 0, Inf]; % C leaf score
    Mis_leaf = [Inf, Inf, Inf, Inf, 0]; % - leaf score
      
    % Construct the struct array of pointers and their score vectors %
    % Nodes are already numbered bottom-up from tree construction %
    pointers = get(phylo_tree, 'Pointers'); % Node pairings map %
    len = length(pointers); % Number of pairings %
    nodes = get(phylo_tree, 'NodeNames'); % Node numbering map %
    pointer_struct = struct('Left', ' ', 'Right', ' ', 'LeftIndex', 0, 'RightIndex', 0, 'ScoreVec', [0,0,0,0,0], 'LabelSeq', ' ');
    for m = 1:len
        left_ind = pointers(m,1); % Left index
        right_ind = pointers(m,2); % Right index
        % Resolves discrepencies in orientation of node indices ; correctly
        % assign left & right nodes at each branch 
        if left_ind > right_ind
            pointer_struct(m).Left = nodes(left_ind);
            pointer_struct(m).Right = nodes(right_ind);
            pointer_struct(m).LeftIndex = pointers(m,1);
            pointer_struct(m).RightIndex = pointers(m,2);
        else
            pointer_struct(m).Left = nodes(right_ind);
            pointer_struct(m).Right = nodes(left_ind);
            pointer_struct(m).LeftIndex = pointers(m,2);
            pointer_struct(m).RightIndex = pointers(m,1);
        end
        % Blank score vector holder (for now) %
        pointer_struct(m).ScoreVec = [0, 0, 0, 0, 0]; 
    end
    
    % Label leaves only first, disregard branches for now %
    % Calculate scores via post-order traversal of phylogenic tree %
    for n = 1:len
        left_rna = pointer_struct(n).Left;
        right_rna = pointer_struct(n).Right;
        check_left = strfind(left_rna, 'Branch');
        check_right = strfind(right_rna, 'Branch');
        % Compare the sequences  - if not branches %
        if isempty(check_left{1}) && isempty(check_right{1}) % Not Branches
            % Declarations - START %
            left_scores = zeros(1,5); % Left Scores
            right_scores = zeros(1,5); % Right Scores
            current_scores = zeros(1,5); % Current sequence score vector
            left_seq = seq_align_struct(pointer_struct(n).LeftIndex).Sequence; % Current left sequence 
            right_seq = seq_align_struct(pointer_struct(n).RightIndex).Sequence; % Current right sequence
            seq_len = length(left_seq); % Length of sequence
            char_seq = zeros(1,seq_len); % Output sequence for branch
            % Declarations - END %
            
            % Go through sequence character by character and compare, apply
            % transition scores with helper function %
            for p = 1:seq_len
                left_ch = left_seq(p);
                right_ch = right_seq(p);
                
                % Left character score inputs %
                switch left_ch
                    case 'A'
                        left_scores = A_leaf;
                    case 'U'
                        left_scores = U_leaf;
                    case 'G'
                        left_scores = G_leaf;
                    case 'C'
                        left_scores = C_leaf;
                    case '-'
                        left_scores = Mis_leaf;
                    case 'N'
                        left_scores = Mis_leaf;
                end
                % Right character score inputs %
                switch right_ch
                    case 'A'
                        right_scores = A_leaf;
                    case 'U'
                        right_scores = U_leaf;
                    case 'G'
                        right_scores = G_leaf;
                    case 'C'
                        right_scores = C_leaf;
                    case '-'
                        right_scores = Mis_leaf;
                    case 'N'
                        right_scores = Mis_leaf;
                end
                % Residues different? %
                if ~strcmp(left_ch, right_ch)               
                    
                    this_score = zeros(1,5);
                    for i = 1:5
                        % Scan the left branch
                        left_minimum = Inf;
                        for j = 1:5
                            left_cost = score_mat(i, j) + left_scores(1, j);
                            left_minimum = min(left_minimum, left_cost);
                        end
                        % Scan the right branch
                        right_minimum = Inf;
                        for k = 1:5
                            right_cost = score_mat(i, k) + right_scores(1, k);
                            right_minimum = min(right_minimum, right_cost);
                        end
                        this_score(1, i) = left_minimum + right_minimum;
                    end
                    current_scores = current_scores + this_score;
                    % Which residue is smallest? %
                    [~, index] = min(current_scores);
                    switch index
                        case 1
                            char_seq(1,p) = 'A';
                        case 2
                            char_seq(1,p) = 'U';
                        case 3
                            char_seq(1,p) = 'G';
                        case 4
                            char_seq(1,p) = 'C';
                        case 5
                            char_seq(1,p) = '-'; 
                    end
                % Residues are the same %    
                else
                    char_seq(1,p) = left_ch;
                end
            end
            % Cast vector back to character sequence, add data to struct %
            pointer_struct(n).LabelSeq = char(char_seq);
            pointer_struct(n).ScoreVec = current_scores;
        end           
    end
    % Next handle Branch + Leaf relationships %
    for q = 1:len
        left_rna = pointer_struct(q).Left;
        right_rna = pointer_struct(q).Right;
        check_left = strfind(left_rna, 'Branch');
        check_right = strfind(right_rna, 'Branch');
        % Compare the sequences  - if branches & leaves %
        if xor(~isempty(check_left{1}), ~isempty(check_right{1}))
            % Declarations - START %
            left_scores = zeros(1,5);
            right_scores = zeros(1,5);
            current_scores = zeros(1,5);
            % Left node contains branch
            if ~isempty(check_left{1})
                ind = strrep(left_rna, 'Branch ', '');
                ind = str2num(ind{1}); % Branch index
                left_seq = pointer_struct(ind).LabelSeq; % Left sequence
                right_seq = seq_align_struct(pointer_struct(q).RightIndex).Sequence; % Right sequence
                left_scores = pointer_struct(ind).ScoreVec; % Left scores
                seq_len = length(left_seq); % Length of current sequence
                char_seq = zeros(1,seq_len); % Current branch sequence
                % Declarations - END %
                % Scan sequence char by char
                for r = 1:seq_len
                    left_ch = left_seq(r); % Current left character
                    right_ch = right_seq(r); % Current right character
                    switch right_ch
                        case 'A'
                            right_scores = A_leaf;
                        case 'U'
                            right_scores = U_leaf;
                        case 'G'
                            right_scores = G_leaf;
                        case 'C'
                            right_scores = C_leaf;
                        case '-'
                            right_scores = Mis_leaf;
                        case 'N'
                            right_scores = Mis_leaf;
                    end
                % Residues different? %
                if ~strcmp(left_ch, right_ch)               
                    
                    this_score = zeros(1,5);
                    for i = 1:5
                        % Scan the left branch
                        left_minimum = Inf;
                        for j = 1:5
                            left_cost = score_mat(i, j) + left_scores(1, j);
                            left_minimum = min(left_minimum, left_cost);
                        end
                        % Scan the right branch
                        right_minimum = Inf;
                        for k = 1:5
                            right_cost = score_mat(i, k) + right_scores(1, k);
                            right_minimum = min(right_minimum, right_cost);
                        end
                        this_score(1, i) = left_minimum + right_minimum;
                    end
                    current_scores = current_scores + this_score;
                    % Which residue is smallest? %
                    [~, index] = min(current_scores);
                    switch index
                        case 1
                            char_seq(1,r) = 'A';
                        case 2
                            char_seq(1,r) = 'U';
                        case 3
                            char_seq(1,r) = 'G';
                        case 4
                            char_seq(1,r) = 'C';
                        case 5
                            char_seq(1,r) = '-'; 
                    end
                % Residues are the same %    
                else
                    char_seq(1,r) = left_ch;
                end
                    
                end
            % Right node contains branch
            elseif ~isempty(check_right{1})
                ind = strrep(right_rna, 'Branch ', '');
                ind = str2num(ind{1}); % Current branch index
                left_seq = seq_align_struct(pointer_struct(q).LeftIndex).Sequence; % Left sequence
                right_seq = pointer_struct(ind).LabelSeq; % Right sequence
                right_scores = pointer_struct(ind).ScoreVec; % Right scores
                seq_len = length(left_seq); % Current sequence length
                char_seq = zeros(1,seq_len); % Current branch sequence
                % Declarations - END %
                % Scan sequence char by char
                for t = 1:seq_len
                    left_ch = left_seq(t); % Current left character
                    right_ch = right_seq(t); % Current right character
                    switch left_ch
                        case 'A'
                            left_scores = A_leaf;
                        case 'U'
                            left_scores = U_leaf;
                        case 'G'
                            left_scores = G_leaf;
                        case 'C'
                            left_scores = C_leaf;
                        case '-'
                            left_scores = Mis_leaf;
                        case 'N'
                            left_scores = Mis_leaf;
                    end
                    % Residues different? %
                if ~strcmp(left_ch, right_ch)               
                    
                    this_score = zeros(1,5);
                    for i = 1:5
                        % Scan the left branch
                        left_minimum = Inf;
                        for j = 1:5
                            left_cost = score_mat(i, j) + left_scores(1, j);
                            left_minimum = min(left_minimum, left_cost);
                        end
                        % Scan the right branch
                        right_minimum = Inf;
                        for k = 1:5
                            right_cost = score_mat(i, k) + right_scores(1, k);
                            right_minimum = min(right_minimum, right_cost);
                        end
                        this_score(1, i) = left_minimum + right_minimum;
                    end
                    current_scores = current_scores + this_score;
                    % Which residue is smallest? %
                    [~, index] = min(current_scores);
                    switch index
                        case 1
                            char_seq(1,t) = 'A';
                        case 2
                            char_seq(1,t) = 'U';
                        case 3
                            char_seq(1,t) = 'G';
                        case 4
                            char_seq(1,t) = 'C';
                        case 5
                            char_seq(1,t) = '-'; 
                    end
                % Residues are the same %    
                else
                    char_seq(1,t) = left_ch;
                end
                    
                end
   
            end
            % Cast vector back to character sequence, add data to struct %
            pointer_struct(q).LabelSeq = char(char_seq);
            pointer_struct(q).ScoreVec = current_scores;
        end
    end
    % Finally; add branch/branch relationships to pointers, resolve discrepencies, only after leaves and branches have been handled.
    for w = 1:len
        
        % Handle only empty branch sequences %
        if isempty(pointer_struct(w).LabelSeq)
            left_rna = pointer_struct(w).Left;
            right_rna = pointer_struct(w).Right;
            check_left = strfind(left_rna, 'Branch');
            check_right = strfind(right_rna, 'Branch');
            if xor(~isempty(check_left{1}), ~isempty(check_right{1}))
            % Declarations - START %
            left_scores = zeros(1,5);
            right_scores = zeros(1,5);
            current_scores = zeros(1,5);
            % Left node contains branch
            if ~isempty(check_left{1})
                ind = strrep(left_rna, 'Branch ', '');
                ind = str2num(ind{1}); % Branch index
                left_seq = pointer_struct(ind).LabelSeq; % Left sequence
                right_seq = seq_align_struct(pointer_struct(q).RightIndex).Sequence; % Right sequence
                left_scores = pointer_struct(ind).ScoreVec; % Left scores
                seq_len = length(left_seq); % Length of current sequence
                char_seq = zeros(1,seq_len); % Current branch sequence
                % Declarations - END %
                % Scan sequence char by char
                for r = 1:seq_len
                    left_ch = left_seq(r); % Current left character
                    right_ch = right_seq(r); % Current right character
                    switch right_ch
                        case 'A'
                            right_scores = A_leaf;
                        case 'U'
                            right_scores = U_leaf;
                        case 'G'
                            right_scores = G_leaf;
                        case 'C'
                            right_scores = C_leaf;
                        case '-'
                            right_scores = Mis_leaf;
                        case 'N'
                            right_scores = Mis_leaf;
                    end
                % Residues different? %
                if ~strcmp(left_ch, right_ch)               
                    
                    this_score = zeros(1,5);
                    for i = 1:5
                        % Scan the left branch
                        left_minimum = Inf;
                        for j = 1:5
                            left_cost = score_mat(i, j) + left_scores(1, j);
                            left_minimum = min(left_minimum, left_cost);
                        end
                        % Scan the right branch
                        right_minimum = Inf;
                        for k = 1:5
                            right_cost = score_mat(i, k) + right_scores(1, k);
                            right_minimum = min(right_minimum, right_cost);
                        end
                        this_score(1, i) = left_minimum + right_minimum;
                    end
                    current_scores = current_scores + this_score;
                    % Which residue is smallest? %
                    [~, index] = min(current_scores);
                    switch index
                        case 1
                            char_seq(1,r) = 'A';
                        case 2
                            char_seq(1,r) = 'U';
                        case 3
                            char_seq(1,r) = 'G';
                        case 4
                            char_seq(1,r) = 'C';
                        case 5
                            char_seq(1,r) = '-'; 
                    end
                % Residues are the same %    
                else
                    char_seq(1,r) = left_ch;
                end
                    
                end
            % Right node contains branch
            elseif ~isempty(check_right{1})
                ind = strrep(right_rna, 'Branch ', '');
                ind = str2num(ind{1}); % Current branch index
                left_seq = seq_align_struct(pointer_struct(q).LeftIndex).Sequence; % Left sequence
                right_seq = pointer_struct(ind).LabelSeq; % Right sequence
                right_scores = pointer_struct(ind).ScoreVec; % Right scores
                seq_len = length(left_seq); % Current sequence length
                char_seq = zeros(1,seq_len); % Current branch sequence
                % Declarations - END %
                % Scan sequence char by char
                for t = 1:seq_len
                    left_ch = left_seq(t); % Current left character
                    right_ch = right_seq(t); % Current right character
                    switch left_ch
                        case 'A'
                            left_scores = A_leaf;
                        case 'U'
                            left_scores = U_leaf;
                        case 'G'
                            left_scores = G_leaf;
                        case 'C'
                            left_scores = C_leaf;
                        case '-'
                            left_scores = Mis_leaf;
                        case 'N'
                            left_scores = Mis_leaf;
                    end
                    % Residues different? %
                if ~strcmp(left_ch, right_ch)               
                    
                    this_score = zeros(1,5);
                    for i = 1:5
                        % Scan the left branch
                        left_minimum = Inf;
                        for j = 1:5
                            left_cost = score_mat(i, j) + left_scores(1, j);
                            left_minimum = min(left_minimum, left_cost);
                        end
                        % Scan the right branch
                        right_minimum = Inf;
                        for k = 1:5
                            right_cost = score_mat(i, k) + right_scores(1, k);
                            right_minimum = min(right_minimum, right_cost);
                        end
                        this_score(1, i) = left_minimum + right_minimum;
                    end
                    current_scores = current_scores + this_score;
                    % Which residue is smallest? %
                    [~, index] = min(current_scores);
                    switch index
                        case 1
                            char_seq(1,t) = 'A';
                        case 2
                            char_seq(1,t) = 'U';
                        case 3
                            char_seq(1,t) = 'G';
                        case 4
                            char_seq(1,t) = 'C';
                        case 5
                            char_seq(1,t) = '-'; 
                    end
                % Residues are the same %    
                else
                    char_seq(1,t) = left_ch;
                end
                    
                end
   
            end
            % Cast vector back to character sequence, add data to struct %
            pointer_struct(w).LabelSeq = char(char_seq);
            pointer_struct(w).ScoreVec = current_scores;
            elseif ~isempty(check_left{1}) && ~isempty(check_right{1})
                ind_left = strrep(left_rna, 'Branch ', '');
                ind_left = str2num(ind_left{1}); % Left branch index
                ind_right = strrep(right_rna, 'Branch ', '');
                ind_right = str2num(ind_right{1}); % Right branch index
                right_seq = pointer_struct(ind_right).LabelSeq; % Right sequence
                right_scores = pointer_struct(ind_right).ScoreVec; % Right scores
                left_scores = pointer_struct(ind_left).ScoreVec; % Left scores
                left_seq = pointer_struct(ind_left).LabelSeq % Left Sequence

            end

        end
        
        
        
        
        
    end
    
    
    % Determine minimum score for tree %
    score_vec = pointer_struct(len).ScoreVec;
    sank_score = min(score_vec);
    
end
