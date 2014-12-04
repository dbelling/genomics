function needle-anchor(varargin)

% Needleman-Wunsch algorithm with anchoring %
% Implemented by Dan Belling %

% Variable argument errors %
% Insufficient input data %
if(nargin < 2)
    error('Please provide at least 2 files to perform sequence alignment on')
% Too many inputs %
elseif(nargin > 3)
    error('Too many inputs provided; only 1 match.txt file is possible with this implementation')
end

% Sequence input %
fID1 = fopen(varargin{1});
fID2 = fopen(varargin{2});

% File access errors %
if(fID1 < 0)
    error('Input file1 not found')
elseif(fID2 < 0)
    error('Input file2 not found')
end

% Anchors Away! (provided there is a match.txt)
if(nargin == 3)
    % Commence error checking
    fID3 = fopen(varargin{3});
    if(fID3 < 0)
        error('Match file not found')
    end    
    match = textscan(fID3, '%s');
    numMatch = length(match{1}) / 4;
    
    % Initialize match vectors
    seq1MatchStart = zeros(1, numMatch);
    seq1MatchEnd = zeros(1, numMatch);
    seq2MatchStart = zeros(1 , numMatch);
    seq2MatchEnd = zeros(1, numMatch);
    
    % Fill start/end match positions in match vectors here        
    i = 1;
    while(i < length(match{1}))
        for j = 1:numMatch
            for k = 1:4
                switch k
                    % Sequence 1 match starting location
                    case 1
                        seq1MatchStart(1,j) = str2num(match{1}{i});
                        i = i+1;                        
                    % Sequence 1 match ending location
                    case 2
                        seq1MatchEnd(1,j) = str2num(match{1}{i});
                        i = i+1;
                    % Sequence 2 match starting location
                    case 3
                        seq2MatchStart(1,j) = str2num(match{1}{i});
                        i = i+1;
                    % Sequence 2 match ending location
                    case 4
                        seq2MatchEnd(1,j) = str2num(match{1}{i});
                        i = i+1;
                end
            end
        end
    end
    clearvars match;
    fclose(fID3);
end

% Cell array initialization %
seq1 = textscan(fID1, '%s');
seq2 = textscan(fID2, '%s');
s1 = size(seq1{1});
s2 = size(seq2{1});
size1 = s1(1,1);
size2 = s2(1,1);
fclose(fID1);
fclose(fID2);
clearvars s1 s2 fID1 fID2;

% Sequence Construction/Concatenation %
sequence1 = '';
sequence2 = '';

% Sequence 1 %
for i=2:size1
    sequence1 = strcat(sequence1,seq1{1}{i});
end

% Sequence 2 %
for j=2:size2
    sequence2 = strcat(sequence2,seq2{1}{j});
end

clearvars size1 size2 seq1 seq2;   

length1 = length(sequence1);
length2 = length(sequence2);

% First iteration of needle-wunsch here; permutations to follow
% Perform Needleman-Wunsch here
if(nargin == 3)
    % Anchoring      
    [scoreMat, backtraceMat] = needlemanWunsch(sequence1, sequence2, seq1MatchStart, seq1MatchEnd, seq2MatchStart, seq2MatchEnd);
    % No anchoring
else
    [scoreMat, backtraceMat] = needlemanWunsch(sequence1, sequence2);
end

% Begin backtrace %
Score = backTrace(sequence1, sequence2, scoreMat, backtraceMat);

% sequence done here.
TotScore = zeros(1,10);
iteration = zeros(1,10);
    
for i=1:10
    sequence1 = randoPerm(sequence1);
    sequence2 = randoPerm(sequence2);
    [scoreMat, backtraceMat] = needlemanWunsch(sequence1, sequence2);
    TotScore(1,i) = backTrace(sequence1, sequence2, scoreMat, backtraceMat);
    iteration(1,i) = i;
end

plot(iteration, TotScore);
disp(sequence1);
disp(sequence2);

end

function[scoreMat, backtraceMat] = needlemanWunsch(varargin)

% Score matrix + Gap penalty + Backtrace cell initialization %
sequence1 = varargin{1};
sequence2 = varargin{2};

length1 = length(sequence1);
length2 = length(sequence2);
scoreMat = zeros(length1 + 1, length2 + 1);
backtraceMat = cell(length1 + 1, length2 + 1);
gap = -2;

% Gap assignment for first column and row %
for i=1:length1+1
    scoreMat(i,1) = gap * (i-1);
end

for j=1:length2+1
    scoreMat(1,j) = gap * (j-1);
end

% Perform anchoring here %
if(nargin > 2)
    seq1MatchStart = varargin{3};
    seq1MatchEnd = varargin{4};
    seq2MatchStart = varargin{5};
    seq2MatchEnd = varargin{6};   
    
   for i = 1:length(seq1MatchStart);
       seq1Start = seq1MatchStart(1,i);
       seq1End = seq1MatchEnd(1,i);
       seq2Start = seq2MatchStart(1,i);
       seq2End = seq2MatchEnd(1,i);
       
       for j = seq1Start:seq1End
           for k = seq2Start:seq2End
               % Scoring alignment is local, only consider beginning and
               % end of sequence, 'localized' needle wunsch
               if(sequence1(j - 1) == sequence2(k - 1))
                   score = 1;
               else
                   score = -3;
               end
               F = [scoreMat(j - 1, k - 1) + score, scoreMat(j - 1, k) + gap, scoreMat(j, k - 1) + gap];
               [scoreMax, condition] = max(F);
               scoreMat(j, k) = scoreMax;
               switch condition
                    case 1
                        backtraceMat{j, k} = [j - 1, k - 1];
                    case 2
                        backtraceMat{j, k} = [j - 1, k];
                    case 3
                        backtraceMat{j, k} = [j, k - 1];
               end
                           
           end          
       end                                          
   end      
end

% Vanilla Needlman-Wunsch %
% Score matrix fill
    for i=2:length1 + 1
        for j=2:length2 + 1
            % Found match?
            if(sequence1(i-1) == sequence2(j-1))
                score = 1;
            % Mismatch
            else
                score = -3;
            end
            % Recurrence relation
            F = [scoreMat(i - 1, j - 1) + score, scoreMat(i - 1, j) + gap, scoreMat(i, j - 1) + gap];
            [scoreMax, condition] = max(F);
            scoreMat(i, j) = scoreMat(i, j) + scoreMax;
            switch condition
                case 1
                    backtraceMat{i, j} = [i - 1, j - 1];
                case 2
                    backtraceMat{i, j} = [i - 1, j];
                case 3
                    backtraceMat{i, j} = [i, j - 1];
            end
        end
    end   
end

function [Score] = backTrace(sequence1, sequence2, scoreMat, backtraceMat)

% Backtrace and alignment set
length1 = length(sequence1);
length2 = length(sequence2);

% Start at bottom right of matrix
i = length1+1;
j = length2+1;

% Initialize cell matrix
alignment = cell(1,max([length1, length2])-1);
alignment{1} = [i, j];

Score = scoreMat(length1, length2);
k = 1;
% Perform backtrace from bottom right of matrix until top left, non-edge
while i > 2 && j > 2
   node = backtraceMat{i,j};
   alignment{k} = node;
   k = k + 1;
   i = node(1);
   j = node(2); 
   Score = Score + max([scoreMat(i - 1, j - 1), scoreMat(i - 1, j), scoreMat(i, j - 1)]);
end

% Terminals have no score
iinit = 0;
jinit = 0;
l = 1;

% Character/Amino Acid search in backtrace %
for k=length(alignment):-1:1
    node = alignment{k};
    i = node(1);
    j = node(2);
    % Character match found in sequence 1, set amino acid here
    if i ~= iinit
        residue1 = sequence1(i-1);
    % Sequence1 gap
    else
        residue1 = '_';
    end
    % Character match found in sequence 2, set amino acid here
    if j ~= jinit
        residue2 = sequence2(j-1);
    % Sequence2 gap
    else
        residue2 = '_';
    end
    % Update loop
    iinit = i;
    jinit = j;
    alignedsequence1(l) = residue1;
    alignedsequence2(l) = residue2;
    l = l + 1;
end

% disp(alignedsequence1);
% disp(alignedsequence2);
alength1 = length(alignedsequence1);
Score = 0;

% Score Tally
for i=1:alength1
    % Gap Found?
    if(alignedsequence1(i) == '_' || alignedsequence2(i) == '_')
        Score = Score - 2;        
    % Mismatch found?
    elseif(alignedsequence1(i) ~= alignedsequence2(i))
        Score = Score - 3;
    % Match Found
    else
        Score = Score + 1;
    end        
end


end

function[sequenceOut] = randoPerm(sequenceIn)
% Cast a candidate sequence string to a double ascii array, permute the
% array randomly, and cast back to original sequence

str = double(sequenceIn);
index = randperm(length(str));
permutation = str(index);
sequenceOut = [char(permutation)];

end
