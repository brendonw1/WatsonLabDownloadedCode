function [logQ,alist,B] = BinaryMatrixUniformRnd(SampN,rN,cN,BN,tflagN,oflagN,pflagN)
%function [logQ,alist,B] = BinaryMatrixUniformRnd(N,r,c,Binput,tFlag,oFlag,pflag)
%
% r is mx1 nonnegative integers
% c is 1xn nonnegative integers
% N is a positive integer
%
% Approximate sampling from the uniform distribution over mxn binary matrices
% with row sums r and column sums c.  Generates N independent samples.  
%
% An error is generated if no binary matrix agrees with r and c.
%
% alist stores the locations of the ones in the samples.  
% If d = sum(r) = sum(c), then alist is 2 x d x N.
%
% The 1-entries of the kth matrix are stored as alist(:,:,k).  The
% (row,column) indices are (alist(1,t,k),alist(2,t,k)) for t=1:d.
%
% If the third output B is requested, then the matrices are explicitly
% created via:
%
% B = false(m,n,N); for k = 1:N, for t = 1:d, B(alist(1,t,k),alist(2,t,k),k) = true; end, end
%
% so that B(:,:,k) is the kth random matrix.
%
% logQ(k)=log(probability that algorithm generates B(:,:,k))
%
% OPTIONS: (using [] for an option selects default)
%
% Binput is a mxn binary matrix.  If it is provided, then alist(:,:,1) agrees
% with Binput and logQ(1) computes the corresponding log probability.
% Use Binput = [] to skip this option.
%
% tFlag = {'rows','cols','fast','slow' (default)} determines whether the algorithm
% loops over rows ('rows') or columns ('cols').  'fast' loops over the
% largest dimension (and is fastest).  'slow' does the opposite, but tends
% to have less variable importance weights (but not always).
%
% oFlag = {'same','descend' (default)} determines the selection
% order for the entries of the dimension that is looped over.
%
% pFlag = {'none','rows','cols','both' (default)} determines whether zero margins
% are pruned before using the algorithm.  
%
% If you use this software, please reference:
%
% M Harrison (2009) A dynamic programming approach for approximate uniform 
%  generation of binary matrices with specified margins. arXiv:0906.1004
% (or see if there is an updated published reference)

% Matthew Harrison
% Oct 8, 2014

% Copyright (C) 2014 Matthew T Harrison
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%    See the GNU General Public License for more details.
%
%    See <http://www.gnu.org/licenses/> for a copy of the license.

%------------------------------------------------------%
%----------------  ARGUMENT CHECKING  -----------------%
%------------------------------------------------------%

if nargin < 7 || isempty(pflagN)
    pflagN = 'both';
end
doPRUNEr = false;
doPRUNEc = false;
switch lower(pflagN)
    case 'both'
        doPRUNEr = true;
        doPRUNEc = true;
    case 'rows'
        doPRUNEr = true;
    case 'cols'
        doPRUNEc = true;
    case 'none'
    otherwise
        error('unknown pFlag')
end

if nargin < 6 || isempty(oflagN)
    oflagN = 'descend';
end

if nargin < 5 || isempty(tflagN)
    tflagN = 'slow';
end
doTRANS = false;
switch lower(tflagN)
    case 'fast'
        if numel(rN) > numel(cN), doTRANS = true; end
    case 'slow'
        if numel(cN) > numel(rN), doTRANS = true; end
    case 'rows'
        doTRANS = true;
    case 'cols'
    otherwise
        error('unknown tFlag')
end

if nargin < 4 || isempty(BN)
    doIN = false;
else
    doIN = true;
end

doA = true;
if nargout < 2, doA = false; end

if ~isscalar(SampN) || SampN < 1 || SampN ~= round(SampN), error('N must be a positive integer'), end

Bsize = [numel(rN),numel(cN),SampN];

%------------------------------------------------------%
%--------------- START: PREPROCESSING -----------------%
%------------------------------------------------------%

% transpose the matrix
if doTRANS
    tmpN = rN;
    rN = cN;
    cN = tmpN;
    if doIN, BN = BN.'; end
end

% shaping
rT = rN(:);
cT = cN(:);

% sizing
mT = numel(rT);
nT = numel(cT);

% check for compatible input
if doIN
    if ~isequal([mT,nT],size(BN)) || ~isequal(sum(BN,2),rT) || ~isequal(sum(BN,1).',cT), error('invalid Binput'), end
end

% order the marginals ... ensure all zeros are at the end
[rsort,rndxT] = sort(rT,'descend');
switch lower(oflagN)
    case 'descend'
        [csort,cndx] = sort(cT,'descend');
        cmax = csort(1);
        cmin = csort(nT);
    case 'same'
        [~,cndx] = sort(logical(cT),'descend'); % only move zeros to the end
        csort = cT(cndx);
        cmax = max(csort);
        cmin = min(csort);
    otherwise
        error('unknown oFlag')
end

% generate the inverse index for the row orders to facilitate fast
% sorting during the updating
irndxT = (1:mT).'; irndxT(rndxT) = irndxT;

% basic input checking
if rsort(1) > nT || rsort(mT) < 0 || cmax > mT || cmin < 0 || any(rsort ~= round(rsort)) || any(csort ~= round(csort)), error('marginal entries invalid'), end

% compute the conjugate of c
cconjT = conjugate_local(csort,mT);

% get the running total of number of ones to assign
countT = sum(rsort);

% get the running total of sum of c squared
ccount2T = sum(csort.^2);

% check for compatible marginals
if countT ~= sum(csort) || any(cumsum(rsort) > cumsum(cconjT)), error('marginal sums invalid'), end

% initialize the memory
logQ = zeros(SampN,1);

if doA, AN = SampN; else AN = 1; end
alist = zeros(2,countT,AN);

% pruning
if doPRUNEr
    for k = mT:-1:1
        if rsort(k) == 0, mT=mT-1; else break, end
    end
end
if doPRUNEc
    for k = nT:-1:1
        if csort(k) == 0, nT=nT-1; else break, end
    end
end

% initialize the memory
M = cmax+3; % index 1 corresponds to -1; index 2 corresponds to 0, index 3 corresponds to 1, ..., index M corresponds to c(1)+1
S = zeros(M,mT);
SS = zeros(M,1);

eps0 = eps(0); % used to prevent divide by zero

%------------------------------------------------------%
%--------------- END: PREPROCESSING -------------------%
%------------------------------------------------------%

% loop over the number of samples
for SampLoop = 1:SampN
    
    %--------------- INITIALIZATION -----------------------%
    if doA, ALoop = SampLoop; else ALoop = 1; end
    
    % copy in initialization
    r = rT;
    rndx = rndxT;
    irndx = irndxT;
    
    cconj = cconjT;
    count = countT;
    ccount2 = ccount2T;

    m = mT;
    n = nT;
    
    % initialize
    place = 0; % most recent assigned column in alist
    logq = 0; % running log probability
    
    %------------------------------------------------------%
    %--------------- START: COLUMN-WISE SAMPLING ----------%
    %------------------------------------------------------%
    
    %-------- loop over columns ------------%
    for c1 = 1:nT
        
        %-----------------------------------------------------------------%
        %------------- START: SAMPLE THE NEXT "COLUMN" -------------------%
        %-----------------------------------------------------------------%
        
        % remember the starting point for this columns
        placestart = place + 1;
        
        %--------------------------------
        % sample a col
        %--------------------------------
        
        label = cndx(c1); % current column label
        
        colval = csort(c1); % current column value
        
        if count == 0, break, end
        
        % update the conjugate
        for i = 1:colval
            cconj(i) = cconj(i)-1;
        end
        % update the number of columns remaining
        n = n - 1;
        
        %------------ DP initialization -----------
        
        smin = colval;
        smax = colval;
        cumsums = count;
        % update the count
        count = count - colval;
        % update running total of sum of c squared
        ccount2 = ccount2 - colval^2;
        
        cumconj = count;
        
        SS(colval+3) = 0;
        SS(colval+2) = 1;
        SS(colval+1) = 0;
        
        % get the constants for computing the probabilities
        if (count == 0) || (m*n == count)
            weightA = 0;
        else
            weightA = m*n/(count*(m*n-count));
            weightA = weightA*(1-weightA*(ccount2-count^2/n))/2;
        end
        
        %----------- dynamic programming ----------
        % loop over (remaining and sorted descending) rows in reverse
        for i = m:-1:1
            
            % get the value of this row and use it to compute the
            % probability of a 1 for this row/column pair
            rlabel = rndx(i);
            val = r(rlabel);
            p = val*exp(weightA*(1-2*(val-count/m)));
            p = p/(n+1-val+p);
            q = 1-p;
            
            % update the feasibility constraints
            cumsums = cumsums - val;
            cumconj = cumconj - cconj(i);
            
            sminold = smin;
            smaxold = smax;
            
            % incorporate the feasibility constraints into bounds on the
            % running column sum
            smin = max(0,max(cumsums-cumconj,sminold-1));
            smax = min(smaxold,i-1);
            
            % DP iteration
            SSS = 0;
            
            SS(smin+1) = 0;  % no need to set S(1:smin) = 0, since it is not accessed
            for j = smin+2:smax+2
                a = SS(j)*q;
                b = SS(j+1)*p;
                apb = a + b;
                SSS = SSS + apb;
                SS(j) = apb;
                S(j,i) = b/(apb+eps0);
            end
            SS(smax+3) = 0;  % no need to set S(smax+4:end) = 0, since it is not accessed
            
            % check for impossible
            %if SSS <= 0, error('algorithm error 1'), end
            
            % normalize to prevent overflow/underflow
            for j = smin+2:smax+2
                SS(j) = SS(j) / SSS;
            end
            
        end
        
        %----------- sampling ----------
        j = 2; % running total (offset to match indexing offset)
        jmax = colval + 2;
        if j < jmax % skip assigning anything when colval == 0
            if doIN
                for i = 1:m
                    % get the transition probability of generating a one
                    p = S(j,i);
                    % get the current row
                    rlabel = rndx(i);
                    if BN(rlabel,label)
                        
                        % if we have a generated a one, then decrement the current
                        % row total
                        val = r(rlabel);
                        r(rlabel) = val-1;
                        
                        % record the entry and update the log probability
                        place = place + 1;
                        logq = logq + log(p);

                        alist(1,place,ALoop) = rlabel;
                        alist(2,place,ALoop) = label;
                        j = j + 1;
                        % the next test is not necessary, but seems more efficient
                        % since all the remaining p's must be 0
                        if j == jmax, break, end
                    else
                        logq = logq + log(1-p);
                    end
                end
            else
                for i = 1:m
                    % get the transition probability of generating a one
                    p = S(j,i);
                    if rand <= p
                        
                        % if we have a generated a one, then decrement the current
                        % row total
                        rlabel = rndx(i);
                        val = r(rlabel);
                        r(rlabel) = val-1;
                                                
                        % record the entry and update the log probability
                        place = place + 1;
                        logq = logq + log(p);

                        alist(1,place,ALoop) = rlabel;
                        alist(2,place,ALoop) = label;
                        j = j + 1;
                        % the next test is not necessary, but seems more efficient
                        % since all the remaining p's must be 0
                        if j == jmax, break, end
                    else
                        logq = logq + log(1-p);
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------%
        %------------- END: SAMPLE THE NEXT "COLUMN" ---------------------%
        %-----------------------------------------------------------------%
        
        if count == 0, break, end
        
        %-----------------------------------------------
        % everything is updated except the sorting
        %-----------------------------------------------
        
        %-----------------------------------------------------------------%
        %------------- START: RESORT THE NEW ROW SUMS --------------------%
        %-----------------------------------------------------------------%
        
        % re-sort the assigned rows
        
        % this code block takes each row that was assigned to the list
        % and either leaves it in place or swaps it with the last row
        % that matches its value; this leaves the rows sorted (descending)
        % since each row was decremented by only 1
        
        % looping in reverse ensures that least rows are swapped first
        for j = place:-1:placestart
            % get the row label and its new value (old value -1)
            k = alist(1,j,ALoop);
            val = r(k);
            % find its entry in the sorting index
            irndxk = irndx(k);
            % look to see if the list is still sorted
            irndxk1 = irndxk + 1;
            if irndxk1 <= m && r(rndx(irndxk1)) > val
                % need to re-sort
                % find the first place where k can be inserted
                irndxk1 = irndxk1 + 1;
                while irndxk1 <= m && r(rndx(irndxk1)) > val
                    irndxk1 = irndxk1 + 1;
                end
                irndxk1 = irndxk1 - 1;
                % now swap irndxk and irndxk1
                rndxk1 = rndx(irndxk1);
                rndx(irndxk) = rndxk1;
                rndx(irndxk1) = k;
                irndx(k) = irndxk1;
                irndx(rndxk1) = irndxk;
            end
        end
        
        %-----------------------------------------------------------------%
        %------------- END: RESORT THE NEW ROW SUMS ----------------------%
        %-----------------------------------------------------------------%
        
        % r(rndx(rndx1:rndxm)) is sorted descending and has exactly those
        % unassigned rows
        % rndx(rndx1:rndxm) still gives the labels of those rows
        % rndx(irndx(k)) = k
        %
        % c(c1+1:cn) is sorted descending and has exactly those unassigned columns
        % cndx(c1+1:cn) still gives the labels of those columns
        %
        % m, n, count, ccount2, ccount2c are valid for the remaining rows, cols
        
    end
    
    logQ(SampLoop) = logq;
    doIN = false; % only compute input on the first instance
end

% undo the transpose
if doTRANS && doA, alist = circshift(alist,1,1); end

% construct B
if nargout >= 3
    B = false(Bsize); 
    for k = 1:SampN
        for t = 1:countT
            B(alist(1,t,k),alist(2,t,k),k) = true; 
        end
    end
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%------------------ END OF MAIN FUNCTION ---------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


% helper function (just to keep everything together... not for efficiency,
% since it is only called once)

function cc = conjugate_local(c,n)
% function cc = conjugate(c,n)
%
% let c(:) be nonnegative integers
% cc(k) = sum(c >== k)  for k = 1:n

cc = zeros(n,1);

for j = 1:numel(c)
    k = c(j);
    if k >= n
        cc(n) = cc(n) + 1;
    elseif k >= 1
        cc(k) = cc(k) + 1;
    end
end

s = cc(n);
for j = n-1:-1:1
    s = s + cc(j);
    cc(j) = s;
end

%-----------------------------------


