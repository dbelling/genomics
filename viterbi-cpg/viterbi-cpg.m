function viterbi-cpg(varargin)

% Hidden markov models and CpG islands on the human genome %
% Viterbi algorithm implementation by:
% Dan Belling

inpFile = varargin{1};
fID1 = fopen(inpFile);
% Basic error handling
if(fID1 < 0)
    error('Function cannot read FASTA sequence file');
    if(nargin < 1)
        error('Expected input sequence file');
    elseif(nargin > 1)
        error('Can only input 1 sequence file');
    end
end
fclose(fID1);

% Sequence input %
[desc, seq] = fastaread(inpFile);
bases = basecount(seq);
disp('Sequence description:');
disp(desc);
fprintf('It contains: \n %d A nucleotides\n %d C nucleotides\n %d G nucleotides\n %d T nucleotides\n\n', bases.A, bases.C, bases.G, bases.T);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Construct the transition probability matrix
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% A+, C+, G+, T+, A-, C-, G-, T- (Probability ordering for row & column)
trans_mat = zeros(8,8);
emit_mat = zeros(8,4);
p = 0.999;
q = 0.95;

% Simplified transmit matrix (2-state model) %
tranSimp_mat = zeros(2,2);
tranSimp_mat(1,1) = p;
tranSimp_mat(1,2) = 1 - p;
tranSimp_mat(2,1) = 1 - q;
tranSimp_mat(2,2) = q;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initializing the +/+ region
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Row 1
trans_mat(1,1) = 0.180 * p;
trans_mat(1,2) = 0.274 * p;
trans_mat(1,3) = 0.426 * p;
trans_mat(1,4) = 0.120 * p;

% Row 2
trans_mat(2,1) = 0.171 * p;
trans_mat(2,2) = 0.368 * p;
trans_mat(2,3) = 0.274 * p;
trans_mat(2,4) = 0.188 * p;

% Row 3
trans_mat(3,1) = 0.161 * p;
trans_mat(3,2) = 0.339 * p;
trans_mat(3,3) = 0.375 * p;
trans_mat(3,4) = 0.125 * p;

% Row 4
trans_mat(4,1) = 0.079 * p;
trans_mat(4,2) = 0.335 * p;
trans_mat(4,3) = 0.384 * p;
trans_mat(4,4) = 0.182 * p;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initializing the -/- region
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Row 5
trans_mat(5,5) = 0.300 * q;
trans_mat(5,6) = 0.205 * q;
trans_mat(5,7) = 0.285 * q;
trans_mat(5,8) = 0.210 * q;

% Row 6
trans_mat(6,5) = 0.322 * q;
trans_mat(6,6) = 0.298 * q;
trans_mat(6,7) = 0.078 * q;
trans_mat(6,8) = 0.302 * q;

% Row 7
trans_mat(7,5) = 0.248 * q;
trans_mat(7,6) = 0.246 * q;
trans_mat(7,7) = 0.298 * q;
trans_mat(7,8) = 0.208 * q;

% Row 8
trans_mat(8,5) = 0.177 * q;
trans_mat(8,6) = 0.239 * q;
trans_mat(8,7) = 0.292 * q;
trans_mat(8,8) = 0.292 * q;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initializing the +/- region
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

for i = 1:4
    for j = 5:8
        trans_mat(i,j) = (1 - p) / 4;
    end
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initializing the -/+ region
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~

for m = 5:8
    for n = 1:4
        trans_mat(m,n) = (1 - q) / 4;
    end
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Emission probability matrix initialization %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Similar to the transition construction, the emission matrix contains
% (row-wise) each of the nucleotides in/out of the CpG island state while
% emitting the specified nucleotide with probability of unity

% A+ / A- %
emit_mat(1,1) = 1;
emit_mat(5,1) = 1;

% C+ / C- %
emit_mat(2,2) = 1;
emit_mat(6,2) = 1;

% G+ / G- %
emit_mat(3,3) = 1;
emit_mat(7,3) = 1;

% T+ / T- %
emit_mat(4,4) = 1;
emit_mat(8,4) = 1;

% ~~~~~~~~~~~~~~~~~
% Viterbi Algorithm
% ~~~~~~~~~~~~~~~~~

n = length(seq);
log_Viterbi = zeros(2,n);

% Decoding % 
% Transition from CpG-island state (row 1) to non-CpG-island state (row 2)%

% To simplify the calculations, take the log2 of each matrix, to maintain
% cumulative sums, instead of probability products of transitions %
logTrans_mat = log2(trans_mat);

% Assume here that the "start" state occurs off of the chromosome (50:50
% chance) to be on/off a CpG island, while then having an equal probability
% for any of the nucleotides to occur (1/4) culminating with a
% "self-transition" which would discriminate the two initial values
nti = seq(1);

% Initialize first column of log_Viterbi %
if nti == 'A'
    log_Viterbi(1,1) = -3 + logTrans_mat(1,1);
    log_Viterbi(2,1) = -3 + logTrans_mat(5,5);
elseif nti == 'C'
    log_Viterbi(1,1) = -3 + logTrans_mat(2,2);
    log_Viterbi(2,1) = -3 + logTrans_mat(6,6);
elseif nti == 'G'
    log_Viterbi(1,1) = -3 + logTrans_mat(3,3);
    log_Viterbi(2,1) = -3 + logTrans_mat(7,7);
elseif nti == 'T'
    log_Viterbi(1,1) = -3 + logTrans_mat(4,4);
    log_Viterbi(2,1) = -3 + logTrans_mat(8,8);
end

% Buildup the viterbi matrix %
for i = 2:n
    % Check current residue, and previous residue
    ntk = seq(i);
    ntl = seq(i-1);
    pos = log_Viterbi(1,i-1);
    neg = log_Viterbi(2,i-1);
    
    % A Nucleotide transmissions
    if ntl == 'A'
        if ntk == 'A'
            log_Viterbi(1,i) = max(pos + logTrans_mat(1,1), neg + logTrans_mat(5,1));
            log_Viterbi(2,i) = max(pos + logTrans_mat(1,5), neg + logTrans_mat(5,5)); 
        elseif ntk == 'C'
            log_Viterbi(1,i) = max(pos + logTrans_mat(1,2), neg + logTrans_mat(5,2));
            log_Viterbi(2,i) = max(pos + logTrans_mat(1,1), neg + logTrans_mat(5,6));
        elseif ntk == 'G';
            log_Viterbi(1,i) = max(pos + logTrans_mat(1,3), neg + logTrans_mat(5,3));
            log_Viterbi(2,i) = max(pos + logTrans_mat(1,1), neg + logTrans_mat(5,7));
        elseif ntk == 'T'
            log_Viterbi(1,i) = max(pos + logTrans_mat(1,4), neg + logTrans_mat(5,4));
            log_Viterbi(2,i) = max(pos + logTrans_mat(1,1), neg + logTrans_mat(5,8));
        end
        
     
    % C Nucleotide transmissions
    elseif ntl == 'C'
        if ntk == 'A'
            log_Viterbi(1,i) = max(pos + logTrans_mat(2,1), neg + logTrans_mat(6,1));
            log_Viterbi(2,i) = max(pos + logTrans_mat(2,5), neg + logTrans_mat(6,5));
        elseif ntk == 'C'
            log_Viterbi(1,i) = max(pos + logTrans_mat(2,2), neg + logTrans_mat(6,2));
            log_Viterbi(2,i) = max(pos + logTrans_mat(2,6), neg + logTrans_mat(6,6));
        elseif ntk == 'G'
            log_Viterbi(1,i) = max(pos + logTrans_mat(2,3), neg + logTrans_mat(6,3));
            log_Viterbi(2,i) = max(pos + logTrans_mat(2,7), neg + logTrans_mat(6,7));
        elseif ntk == 'T'
            log_Viterbi(1,i) = max(pos + logTrans_mat(2,4), neg + logTrans_mat(6,4));
            log_Viterbi(2,i) = max(pos + logTrans_mat(2,8), neg + logTrans_mat(6,8));
        end
        
    % G Nucleotide transmissions
    elseif ntl == 'T'
        if ntk == 'A'
            log_Viterbi(1,i) = max(pos + logTrans_mat(3,1), neg + logTrans_mat(7,1));
            log_Viterbi(2,i) = max(pos + logTrans_mat(3,5), neg + logTrans_mat(7,5));
        elseif ntk == 'C'
            log_Viterbi(1,i) = max(pos + logTrans_mat(3,2), neg + logTrans_mat(7,2));
            log_Viterbi(2,i) = max(pos + logTrans_mat(3,6), neg + logTrans_mat(7,6));
        elseif ntk == 'G'
            log_Viterbi(1,i) = max(pos + logTrans_mat(3,3), neg + logTrans_mat(7,3));
            log_Viterbi(2,i) = max(pos + logTrans_mat(3,7), neg + logTrans_mat(7,7));
        elseif ntk == 'T'
            log_Viterbi(1,i) = max(pos + logTrans_mat(3,4), neg + logTrans_mat(7,4));
            log_Viterbi(2,i) = max(pos + logTrans_mat(3,8), neg + logTrans_mat(7,8));
        end
        
    % T Nucleotide transmissions
    elseif ntl == 'G'
        if ntk == 'A'
            log_Viterbi(1,i) = max(pos + logTrans_mat(4,1), neg + logTrans_mat(8,1));
            log_Viterbi(2,i) = max(pos + logTrans_mat(4,5), neg + logTrans_mat(8,5));
        elseif ntk == 'C'
            log_Viterbi(1,i) = max(pos + logTrans_mat(4,2), neg + logTrans_mat(8,2));
            log_Viterbi(2,i) = max(pos + logTrans_mat(4,6), neg + logTrans_mat(8,6));
        elseif ntk == 'G'
            log_Viterbi(1,i) = max(pos + logTrans_mat(4,3), neg + logTrans_mat(8,3));
            log_Viterbi(2,i) = max(pos + logTrans_mat(4,7), neg + logTrans_mat(8,7));
        elseif ntk == 'T'
            log_Viterbi(1,i) = max(pos + logTrans_mat(4,4), neg + logTrans_mat(8,4));
            log_Viterbi(2,i) = max(pos + logTrans_mat(4,8), neg + logTrans_mat(8,8));
        end
        
    % N Nucleotides? (Skip)
    else
        log_Viterbi(1,i) = log_Viterbi(1, i - 1);
        log_Viterbi(2,i) = log_Viterbi(2, i - 1);
    end
    
end

% Now we must backtrack the matrix to construct the CpG analysis string %
% CpGSeq = 'CN';
%
vitSeq_mat = zeros(1,n);

for j = 1:n
    mMax = max(log_Viterbi(1,j), log_Viterbi(2,j));
    
    if mMax == log_Viterbi(1,j) % CpG Island higher probability
        vitSeq_mat(1,j) = 1;
    elseif mMax == log_Viterbi(2,j) % Normal region higher probability
        vitSeq_mat(1,j) = -1;
    end
end

% Backtrace computation time was excessive for the long DNA sequence;
% maintained a vector of 1's and -1's during runtime instead to label islands 
% on the DNA strand. Could not apply toy method to long seq.

% Further analysis / motif mapping / gene ID %
% Find islands of length 200+ %

start = 0;
stop = 0;
count = 0;
isle_mat = zeros(3,50); % Found ~50 CpG islands of length 200 +
write = 1;

for k = 1:n
    x = vitSeq_mat(1,k);
    
    if count == 0
        start = k;
    end
        
    if x == 1
        count = count + 1;
    elseif x == -1
        if count ~= 0
            stop = k - 1;
            if count > 200 % Only write 200+ residue islands to mat
                isle_mat(1,write) = start;
                isle_mat(2,write) = stop;
                isle_mat(3,write) = count;
                write = write + 1;
                count = 0;
            end
        end 
        count = 0;
    end          
end

write = write - 1;

% Gene mapping %
% Note: the program expects that the motifs_Chr21.txt file is present in
% the working directory. Genes are hard coded below to prevent repeat hits
chrStart = 43507093;
chrEnd = 46944323;

fID2 = fopen('motifs_Chr21.txt');
motifScan = textscan(fID2, '%s');
fileLen = length(motifScan{1});
motifLen = length(motifScan{1}) / 5;
motifs = cell(motifLen, 1);
motifLocs = zeros(2, motifLen);
motifMap = cell(write, 1);

% Write to accessible cells/matrices %
motifs{1} = motifScan{1}{1};
motifLocs(1,1) = str2num(motifScan{1}{3});
motifLocs(2,1) = str2num(motifScan{1}{4});

for w = 2:motifLen
    motifs{w} = motifScan{1}{5 * (w - 1) + 1};
    motifLocs(1,w) = str2num(motifScan{1}{(5 * w) - 2});
    motifLocs(2,w) = str2num(motifScan{1}{(5 * w) - 1});
end

% Perform mapping of motifs here %
for s = 1:write
    isle_start = isle_mat(1, s);
    isle_stop = isle_mat(2, s);
    for t = 1:motifLen
        mot_start = motifLocs(1,t) - chrStart;
        mot_end = motifLocs(2,t) - chrStart;
        if (mot_start >= isle_start) && (mot_end <= isle_stop)
            motifMap{s} = motifs{t};
            break
        end           
    end
end



% Hard coding this part to remove redundancies from the Chr21 file %
genes = zeros(2,42);
genes(1,1) = 43541266 - chrStart;
genes(2,1) = 43548179 - chrStart;
genes(1,2) = 43621786 - chrStart;
genes(2,2) = 43717352 - chrStart;
genes(1,3) = 43641283 - chrStart;
genes(2,3) = 43641396 - chrStart;
genes(1,4) = 43676459 - chrStart;
genes(2,4) = 43717352 - chrStart;
genes(1,5) = 43720384 - chrStart;
genes(2,5) = 43720455 - chrStart;
genes(1,6) = 43732185 - chrStart;
genes(2,6) = 43735519 - chrStart;
genes(1,7) = 43802835 - chrStart;
genes(2,7) = 43816955 - chrStart;
genes(1,8) = 43824011 - chrStart;
genes(2,8) = 43830697 - chrStart;
genes(1,9) = 43837606 - chrStart;
genes(2,9) = 43837702 - chrStart;
genes(1,10) = 44313391 - chrStart;
genes(2,10) = 44329209 - chrStart;
genes(1,11) = 44337150 - chrStart;
genes(2,11) = 44339420 - chrStart;
genes(1,12) = 44373890 - chrStart;
genes(2,12) = 44376488 - chrStart;
genes(1,13) = 44394643 - chrStart;
genes(2,13) = 44452285 - chrStart;
genes(1,14) = 44473303 - chrStart;
genes(2,14) = 44495964 - chrStart;
genes(1,15) = 44513067 - chrStart;
genes(2,15) = 44527688 - chrStart;
genes(1,16) = 44579355 - chrStart;
genes(2,16) = 44581362 - chrStart;
genes(1,17) = 45079432 - chrStart;
genes(2,17) = 45115960 - chrStart;
genes(1,18) = 45152259 - chrStart;
genes(2,18) = 45156956 - chrStart;
genes(1,19) = 45168246 - chrStart;
genes(2,19) = 45168349 - chrStart;
genes(1,20) = 45209457 - chrStart;
genes(2,20) = 45222360 - chrStart;
genes(1,21) = 45275969 - chrStart;
genes(2,21) = 45276133 - chrStart;
genes(1,22) = 45285116 - chrStart;
genes(2,22) = 45405179 - chrStart;
genes(1,23) = 45339676 - chrStart;
genes(2,23) = 45339782 - chrStart;	
genes(1,24) = 45416203 - chrStart;
genes(2,24) = 45416301 - chrStart;	
genes(1,25) = 45666223 - chrStart;
genes(2,25) = 45682099 - chrStart;	
genes(1,26) = 45710274 - chrStart;	
genes(2,26) = 45718102 - chrStart;	
genes(1,27) = 45725211 - chrStart;	
genes(2,27) = 45747253 - chrStart;	
genes(1,28) = 45749773 - chrStart;	
genes(2,28) = 45754194 - chrStart;	
genes(1,29) = 45751119 - chrStart;	
genes(2,29) = 45751543 - chrStart;	
genes(1,30) = 45876518 - chrStart;	
genes(2,30) = 45877311 - chrStart;	
genes(1,31) = 45894981 - chrStart;	
genes(2,31) = 45895050 - chrStart;	
genes(1,32) = 45917777 - chrStart;
genes(2,32) = 46131495 - chrStart;	
genes(1,33) = 45993624 - chrStart;
genes(2,33) = 45994841 - chrStart;	
genes(1,34) = 46117342 - chrStart;
genes(2,34) = 46117722 - chrStart;
genes(1,35) = 46188956 - chrStart;
genes(2,35) = 46221738 - chrStart;	
genes(1,36) = 46225777 - chrStart;	
genes(2,36) = 46237966 - chrStart;	
genes(1,37) = 46269515 - chrStart;	
genes(2,37) = 46293644 - chrStart;	
genes(1,38) = 46286280 - chrStart;	
genes(2,38) = 46288992 - chrStart;	
genes(1,39) = 46305871 - chrStart;
genes(2,39) = 46340965 - chrStart;	
genes(1,40) = 46355573 - chrStart;
genes(2,40) = 46359828 - chrStart;	
genes(1,41) = 46554641 - chrStart;	
genes(2,41) = 46564597 - chrStart;	
genes(1,42) = 46875424 - chrStart;	
genes(2,42) = 46933633 - chrStart;	

geneStr = {'UMODL1'; 'ABCG1'; 'AP001622.1'; 'ABCG1'; 'AP001623.1'; 'TFF3'; 'TMPRSS3'; 'UBASH3A'; 'AP001624.1'; 'NDUFV3'; 'AP001629.1'; 'AP001630.1'; 'PKNOX1'; 'CBS';
    'U2AF1'; 'AP001631.10'; 'RRP1B'; 'PDXK'; 'AP001052.1'; 'RRP1'; 'MYL6P'; 'AGPAT3'; 'AP001054.1'; 'AB001523.1'; 'DNMT3L'; 'AIRE'; 'PFKL'; 'C21orf2'; 'AP001062.7'; 
    'LRRC3'; 'AP001065.1'; 'C21orf29'; 'KRTAP10-4'; 'KRTAP10-12'; 'UBE2G2'; 'SUMO3'; 'PTTG1IP'; 'AL773603.1'; 'ITGB2'; 'C21orf67'; 'ADARB1'; 'COL18A1'};

len = length(geneStr);


% Print pertinent CpG island info %
fprintf('Total CpG Islands found: %d\n', write);

for l = 1:write
    start = isle_mat(1,l);
    stop = isle_mat(2,l);
    count = isle_mat(3,l);
    found = false;
    for p = 1:len
        diff = genes(1, p) - stop;
        if diff <= 500 % Genes within 500 bp of island
            found = true;
            gene = geneStr{p};
        end
    end
    if found
        if ~isempty(motifMap{l})
            fprintf('CpG Island %d: %d bp (%d - %d): Gene - %s Motif - %s\n', l, count, start, stop, gene, motifMap{l});
        else
            fprintf('CpG Island %d: %d bp (%d - %d): Gene - %s\n', l, count, start, stop, gene);
        end
    else
        fprintf('CpG Island %d: %d bp (%d - %d)\n', l, count, start, stop);
    end       
end

end
