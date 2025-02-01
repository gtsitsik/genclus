function [U,A,B,eigvals,iters,time,obj_cur] = genclus(X,R,M,thres,...
                       max_iters,constraints,rho,print_type,mtimesx_exists)
% [U,A,B] = GENCLUS(X,R,M) calculates GenClus on the adjacency tensor X
% with a total of R components for the node embeddings and M components for
% the view embeddings. U, A and B are the calculated factor matrices that
% effectively aim to model X based on a PARAFAC with factor matrices U, U
% and A*B.
%
% Further details:
% U - Its i-th row corresponds to the node embedding that captures the
%     structure of the i-th node for all view clusters simultaneously.  
% B - Its i-th row is used to create weighted versions of all rows of U
%     which in turn correspond to node embeddings that capture the
%     structure of the nodes only for the i-th view cluster.
% A - Its i-th row corresponds to the i-th view embedding.
%
% [U,A,B] = GENCLUS(X,R,M,thres,max_iters,...
%                                constraints,rho,print_type,mtimesx_exists)
% specifies additional input parameters as defined below:
% thres           -  Convergence threshold.
% max_iters       -  Maximum iterations allowed.
% constraints     -  A vector with two characters indicating the
%                    constraints on A and B, respectively. The allowed
%                    values of each character are '1','+' and 'U'
%                    indicating non-zero values constrained all to 1,
%                    constrained to be non-negative and being
%                    unconstrained, respectively.  
% rho             -  Please ignore this value and set it equal to an empty
%                    array if necessary.
% print_type      -  Allows increasing levels of message verbosity by
%                    setting it to "none","basic" and "all", respectively. 
% mtimesx_exists  -  A flag indicating the availability of the mtimesx
%                    utility.
%
% Note: Any of the inputs can either be ommited or set equal to an empty
% array [], which will set them to their default value.
%
% [U,A,B,eigvals,iters,time,obj_cur] = GENCLUS(X,R,M,...)
% specifies additional outputs as defined below:
% eigvals  -  The eigenvalues based on which B is constructed. 
% iters    -  The actual number of iterations.
% time     -  The total run duration in seconds.
% obj_cur  -  The final value of the objective function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BSD 3-Clause License
% 
% Copyright (c) 2025, Yorgos Tsitsikas
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if nargin == 0 || isempty(X)
    I = 10;
    K = 4;
    R_true = 2;
    M_true = 2;
    U_true = orth(randn(I,R_true));
    A_true = double(randi(M_true,K,1)==[1:M_true]);
    B_true = double(randi(R_true,M_true,1)==[1:M_true])';
    %     A_true = [1 1 1 0 0 0
    %               0 0 0 1 1 1]';
    %     B_true = [1 1 0 0;
    %               0 0 10 1];
    X = parafac2full(U_true,U_true,A_true*B_true*1+randn(K,R_true)*0);
    
    %         X = ones(10,10);
    noise = rand(size(X))/1e3;
    noise = noise + permute(noise,[2 1 3]);
    X = X + noise;
    %     X = rand(100,100,10);
    X = (X + permute(X,[2 1 3]))/2;
end
if nargin <= 1 || isempty(R)
    R = size(U_true,2);
end
if nargin <= 2 || isempty(M)
    M = size(A_true,2);
end

if nargin <= 3 || isempty(thres)
    thres = 1e-6;
end
if nargin <= 4 || isempty(max_iters)
    max_iters = 50;
end
if nargin <=5 || isempty(constraints)
    constraints='++';
end
if nargin <=6 || isempty(rho)
    rho = 0.00;
end
if nargin<=7 || isempty(print_type)
    print_type ="basic";
end
if nargin<=8 || isempty(mtimesx_exists)
    mtimesx_exists = exist('mtimesx')==3;
end
K = size(X,3);
I = size(X,1);
X_vec = reshape(X,prod(size(X,[1,2])),size(X,3));

% Construct initial U, A and B
switch  constraints(1)
    case '1'
        A = (randi(M,K,1)==[1:M]);
    case '+'
        A = (randi(M,K,1)==[1:M]).*rand(K,1);
    case 'U'
        A = (randi(M,K,1)==[1:M]).*randn(K,1);
    case 'R'
        A = randn(K,M);
        X_3 = reshape(permute(X,[3,1,2]),size(X,3),[])';
    otherwise
        error("Invalid constraint for A")
end
A = double(A);

switch  constraints(2)
    case '1'
        B = ([1:M]'==randi(M,1,R));
    case '+'
        B = ([1:M]'==randi(M,1,R)).*rand(1,R);
    case 'U'
        B = ([1:M]'==randi(M,1,R)).*randn(1,R);
    otherwise
        error("Invalid constraint for B")
end
B = double(B);

U = nan(size(X,1),R);
for i = 1:size(B,1)
    tmp = find(B(i,:));
    U(:,tmp) = orth(randn(size(X,1),numel(tmp)));
end

obj_prev = inf;

obj_cur = obj_eval(X,U,A,B,rho,constraints,mtimesx_exists);
obj_cur_ = obj_cur;
obj_prev_ = obj_prev;

iters = 0;


% norm(X_rec(:)-X_rec_prev(:))/norm(X_rec_prev(:))
eig_cutoff = 0.5;

print_progress("begin")
while (abs(obj_cur-obj_prev)>=thres*obj_prev || (isempty(R) && eig_cuqetoff<0.8) ) && iters<max_iters
    iters = iters+1;
    % Update U and B
    d = [];
    U_all = [];
    B_tmp = [];
    eigvals = [];
    for m = 1:M
        switch constraints(2)
            case '1'
                L = 2*sum(X.*reshape(A(:,m),1,1,[]),3)-norm(A(:,m))^2*eye(size(X,[1,2]));
            case {'U','+'}
                if norm(A(:,m))>eps
                    L = sum(X.*reshape(A(:,m),1,1,[]),3)/norm(A(:,m));
                else
                    L = zeros(size(X,[1,2]));
                end
        end
        %                 D1 = pinv(diag(sqrt(sum(L,2))));
        %                 D2 = pinv(diag(sqrt(sum(L,1))));
        %                 L = D1*L*D2;
        %
        %                 L = (L+L')/2; %makes sure L is still symmetric
        % tmp = norm(L);
        % if tmp>eps
        %        L = L/tmp;
        % end
        [cur_U,cur_D] = eig(L);
        
        cur_d = diag(cur_D)';
        %         if constraints(2)=='+'
        %             cur_d(cur_d<0)=0;
        %         end
        U_all = [U_all cur_U];
        eigvals(:,m)= flip(cur_d);
        d = [d cur_d];
        B_tmp = [B_tmp m*ones(1,numel(cur_d))];
    end
    
    switch constraints(2)
        case {'1','+'}
            if constraints(2)=='+'
                d(d<0)=0;
            end
            [~,inds] = sort(d,'descend');
        case 'U'
            [~,inds] = sort(abs(d),'descend');
        otherwise
            error("Invalid constraint for B")
    end
    
    if ~isempty(R)
        inds = inds(1:R);
    else
        inds = find(d>eig_cutoff);
        if isempty(inds)
            [~,inds] = max(d);
        end
        eig_cutoff = eig_cutoff*1.05;
    end
    U = U_all(:,inds);
    B_tmp = B_tmp(inds);
    B = double([1:M]'==B_tmp);
    if max(constraints(2)==['U','+'])
        A_div=vecnorm(A,2,1)';
        A_div(A_div<eps) = 1; %handles all-zeros columns of A
        B = B.*d(inds)./A_div;
    end
    print_progress("U,B")
    
    % TODO: Normalize A or B after each update to avoid numerical instabilities.
    % Update A
    switch constraints(1)
        case {'1','U','+'}
            C = permute(reshape(parafac2full(U,U,B,mtimesx_exists),prod(size(X,[1,2])),M),[1 3 2]);
            switch constraints(1)
                case '1'
                    dists = vecnorm(X_vec-C,2,1);
                    [~,assgn] = min(dists,[],3);
                case {'U','+'}
                    C_norms = vecnorm(C,2,1);
                    C_norms(C_norms<eps) = 1;
                    in_prods = sum(X_vec.* (C./C_norms),1);
                    if constraints(1) == '+'
                        in_prods(in_prods<0)=0;
                    end
                    [~,assgn] = max(abs(in_prods),[],3);
            end
            
            A = double(assgn'==[1:M]);
            if max(constraints(1)==['U','+'])
                A = A.*permute(in_prods,[2,3,1])./reshape(C_norms(assgn),[],1);
            end
            
            % Updates A as in ComClus. May cause updates of U and B to not be
            % monotonically non-increasing
%         case 'R'
%             S = khatrirao(U,U)*B';
%             A = A.*(abs((2*X_3'*S)./(2*A*S'*S+rho)).^(1/2));
% %             A=rand(size(A));
% %             A =A.*(A==max(A,[],2));
        otherwise
            error("Invalid constraint for A")
    end
    print_progress("A  ")
    
    obj_prev = obj_cur;
    obj_cur = obj_eval(X,U,A,B,rho,constraints,mtimesx_exists);
end
print_progress("end")


    function print_progress(var)
        if print_type~="nothing"
            obj_prev_ = obj_cur_;
            X_rec_ = parafac2full(U,U,A*B,mtimesx_exists);
            obj_cur_ =  obj_eval(X,U,A,B,rho,constraints,mtimesx_exists);
            obj_change_ = (obj_cur_-obj_prev_)/obj_prev_*100;
            
            msg="";
            if max(print_type == ["all","basic"])
                msg = var+" iters: " + iters + " - rec error:   "+num2str(norm(X(:)-X_rec_(:))/norm(X(:))*100) + "   " +  obj_change_;
                if print_type == "all"
                    msg = "GenClus"+constraints+" R:"+R+" M:"+M+" thres:"+thres+" max_iters:"+max_iters +" constraints:"+ constraints +" rho:"+rho+" | "+msg;
                end
            end
            if (obj_change_>-inf) || var=="begin" || var=="end"|| var==""
                disp(msg);
                if var=="end"
                    disp("---------------------------------------------------")
                end
            end
            
            if obj_change_>1e-12% && var =="U,B"
                obj_change_;
            end
        end
    end
time = toc;
end


function L_all = obj_eval(X,U,A,B,rho,constraints,mtimesx_exists)
X_rec = parafac2full(U,U,A*B,mtimesx_exists);
L_all = norm(X(:)-X_rec(:))^2+rho*norm(A(:),1)*(constraints(1)=='R');
end


function X = parafac2full(A,B,C,mtimesx_exists)
% PARAFAC2FULL Compose 3-mode PARAFAC.
%
%    X = PARAFAC2FULL(A,B,C) returns a tensor X which is the composition of
%    a PARAFAC with factor matrices A,B and C.
%
%    X = PARAFAC2FULL(A,B,C,mtimesx_exists) allows avoiding the expensive 
%    exist operation for 'mtimesx'. This is useful when PARAFAC2FULL is 
%    called multiple times.
if nargin==3
    mtimesx_exists = exist('mtimesx','file')==3;
end
if mtimesx_exists 
    X = mtimesx(A.*permute(C,[3 2 1]),B');
else

    X = zeros([size(A,1),size(B,1),size(C,1)]);
    for i =1:size(C,1)
        X(:,:,i) = (A.*C(i,:))*B';
    end
end

% Third way of calculating X_rec using tensor toolbox. Can be very slow.
% X = double(ktensor({A,B,C}));
end
