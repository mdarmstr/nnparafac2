function [Bk,A,Dk,Bs,ssr,timeOut] = nnparafac2(varargin)
%
%Implementation of the Flexible-Coupling (Non-Negative) PARAFAC2, as described by Cohen and Bro:
%
%Cohen, Jeremy E., and Rasmus Bro. "Nonnegative PARAFAC2: A flexible coupling approach." International Conference on Latent Variable Analysis and Signal Separation. Springer, Cham, 2018.
%
%Features estimates for signal-to-noise ratio for each slice to determine an appropriate value for the coupling
%constant
%
%Currently there is a dependency on fcnnls (Fast Combimatorial Non-Negative
%Least Squares) as a non-negative least squares solver:
%https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e5925d25-4a80-11e4-9553-005056977bd0/0b2678a0-fdea-4929-820c-ef5a41b26b46/previews/source/fcnnls.m/index.html
%
%(c) Michael D. Sorochan Armstrong, James J. Harynuk 2021

tic

if size(varargin,2) < 2
    error('Must declare tensor, number of factors');
else
end

Xijkl = varargin{1};
k = varargin{2};

%Turn Xk into a cell array if it isn't already
if ~iscell(Xkl)
    for kl = 1:sz(3)
        Xkc{kl} = Xkl(:,:,kl);
    end
    Xkl = Xkc;
end

%Determining the size of each cell for initialisation purposes
cellsz = cellfun(@size,Xkl,'uni',false);

switch size(varargin,2)
    case 6
        disp('Utilising input initialisations for A,Dk,Bk, and Bs');
        Bsi = varargin{6};
        Bki = varargin{5};
        Dki = varargin{4};
        Ai = varargin{3};
        [Bk,A,Dk,Bs,ssr] = nnparafac2als(Xkl,k,Bki,Ai,Dki,Bsi,1000);
        
    case 5
        Bsi = eye(k,k);
        disp('Utilising input initialisations for A,Dk, and Bk');
        Bki = varargin{5};
        Dki = varargin{4};
        Ai = varargin{3};
        [Bk,A,Dk,Bs,ssr] = nnparafac2als(Xkl,k,Bki,Ai,Dki,Bsi,1000);
    case 4
        Bsi = eye(k,k);
        for kl = 1:sz(3)
            Bki{kl} = rand(cellsz{kl}(1),k);
        end
        disp('Utilising input initialisations for A, and Dk');
        Dki = varargin{4};
        Ai = varargin{3};
        [Bk,A,Dk,Bs,ssr] = nnparafac2als(Xkl,k,Bki,Ai,Dki,Bsi,1000);
        
    case 3
        Bsi = eye(k,k);
        for kl = 1:sz(3)
            Bki{kl} = rand(cellsz{kl}(1),k);
            Dki(:,:,kl) = eye(k,k);
        end
        disp('Utilising input initialisations for A');
        Ai = varargin{3};
        [Bk,A,Dk,Bs,ssr] = nnparafac2als(Xkl,k,Bki,Ai,Dki,Bsi,1000);
    case 2
        disp('Utilising random initialisations - best of 10 initial estimates');
        for re = 1:2 %change this to 10
            Bs = eye(k,k);
            for kl = 1:sz(3)
                Bki{kl} = rand(cellsz{kl}(1),k);
                Dk(:,:,kl) = eye(k,k);
            end
            A = rand(sz(2),k);
            [Bkit,Ait,Dkit,Bsit,ssr] = nnparafac2als(Xkl,k,Bki,A,Dk,Bs,10);
            Bki(:,:,:,re) = Bkit; %This is unresolved currently
            Ai(:,:,re) = Ait;
            Dki(:,:,:,re) = Dkit;
            Bsi(:,:,re) = Bsit;
            ssr_rand(re) = ssr(end);
            
            disp(sprintf('Initialisation %d of %d',re,10)) %#ok
        end
        [~,idx] = min(ssr_rand);
        [Bk,A,Dk,Bs,ssr] = nnparafac2als(Xkl,k,Bki(:,:,:,idx),Ai(:,:,idx),Dki(:,:,:,idx),Bsi(:,:,idx),1000);
end

timeOut = toc;

end

function [Bk,A,Dk,Bs,SSR] = nnparafac2als(Xkl,k,Bki,Ai,Dki,Bsi,maxiter)

sz(3) = size(Xkl,2); sz(2) = size(Xkl{1,1},2);
cellsz = cellfun(@size,Xkl,'uni',false);

Bk = Bki; A = Ai; Dk = Dki; Bs = Bsi;

for kk = 1:sz(3)
    [U, ~, V] = svds(Bk{kk}*Bs,k);
    Pk{kk} = U*V';
    Xh{kk} = Bk{kk}*Dk(:,:,kk)*A';
    mk(kk) = sum(sum(Bk{kk}*Dk(:,:,kk)*A'))./sum(sum(Bk{kk}));
    ssr1(kk) = sum(sum((Xkl{kk} - Xh{kk}).^2)) + mk(kk).*sum(sum((Bk{kk}-Pk{kk}*Bs).^2)); %updated cost function
end

ssr2 = 1;
iter = 1;
eps = 1e-7;
YNorm = vertcat(Xkl{:});
YNorm = sum(YNorm(:).^2);
ssr1 = sum(ssr1);


while abs(ssr1-ssr2)/ssr1 > eps && abs(ssr1 - ssr2) > abs(ssr1)*eps && iter < maxiter
    
    ssr1 = ssr2;
    
    %Pk Estimation
    for kk = 1:sz(3)
        if iter > 1
            [U, ~, V] = svds(Bk{kk}*Bs,k);
            Pk{kk} = U*V';
            %Bs Estimation
            BsT(:,:,kk) = mk(kk).*Pk{kk}'*Bk{kk};
            BkDk{kk} = Bk{kk}*Dk(:,:,kk);
        else
            BsT(:,:,kk) = mk(kk).*Pk{kk}'*Bk{kk};
            BkDk{kk} = Bk{kk}*Dk(:,:,kk);
        end
    end
    
    Bs = 1/(sum(mk))*(sum(BsT,3));
    
    for kk = 1:k
        Bs(:,kk) = Bs(:,kk)./norm(Bs(:,kk));
    end
    
    BkDkIK = vertcat(BkDk{:});
    
    if iter == 1
        Xjki = vertcat(Xkl{:});
    end
    
    for ki = 1:sz(2)
        A1(ki,:) = fcnnls([],[],BkDkIK'*BkDkIK,BkDkIK'*Xjki(:,ki));
    end
    
    if any(sum(A1,2) == 0)
        A1 = 0.9*A + 0.1*A1;
    end
    
    for kk = 1:k
        A1(A1 == 0) = 1e-20;
        A1(:,kk) = A1(:,kk)./norm(A1(:,kk));
    end
    
    %Bk Estimation
    for kk = 1:sz(3)
        for ii = 1:cellsz{kk}(1)
            Bkt(ii,:) = fcnnls([],[],(Dk(:,:,kk)*(A1'*A1)*Dk(:,:,kk) + mk(kk)*eye(k)), (Xkl{kk}(ii,:)*A1*Dk(:,:,kk) + mk(kk)*Pk{kk}(ii,:)*Bs)'); %#ok
        end
        
        Bk{kk} = Bkt;
        
        for rr = 1:k
            tmp = Bk{kk}(:,rr);
            tmp(tmp == 0) = 1e-20;
            Bk{kk}(:,rr) = tmp;
            Bk{kk}(:,rr) = Bk{kk}(:,rr)./norm(Bk{kk}(:,rr));
        end
    end
    
    for kk = 1:sz(3)
        Dktemp = diag((pinv(Bk{kk}'*Bk{kk})*Bk{kk}'*Xkl{kk})/A1');
        Dk(:,:,kk) = diag(Dktemp);
    end
    
    if iter == 1 %If this is the first iteration, define muk
        for kk = 1:sz(3)
            S = svds(Xkl{kk} - mean(Xkl{kk}),2); %perform SVD on mean-centred data as an estimate for SNR.
            SNR = S(1)^2/(S(2)^2);
            mk(kk) = 10^(-SNR/10)*sum((sum(sqrt((Xkl{kk} - Bk{kk}*Dk(:,:,kk)*A1').^2)))/sum(sum(sqrt((Bk{kk} - Pk{kk}*Bs).^2))));
        end
    elseif iter < 10
        for kk = 1:sz(3)
            mk(kk) = min(mk(kk)*1.05,1e12); %Growing mk with each iteration.
        end
    else
        
    end
    
    for kl = 1:sz(3)
        res_mdl(kl) = norm(Xkl{kl} - Bk{kl}*Dk(:,:,kl)*A1','fro')^2;
        res_cpl(kl) = mk(kl).*norm((Bk{kl}-Pk{kl}*Bs),'fro')^2;
    end
    
    ssr2 = sum(res_mdl + res_cpl)/YNorm;
    SSR(iter) = ssr2;
    
    %Disply output, progress of the algorithm
    if iter == 1
        varNames = {'Iteration','Absolute Error','Relative Error','SSR','mk'};
        fprintf(1,'\n%s\t\t%s\t\t%s\t\t%s\t\t%s\n',varNames{:})
        fprintf(1,' \t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n',[iter,abs(ssr2-ssr1),abs(ssr2-ssr1)/abs(ssr2),SSR(iter),mean(mk)]);
    else
        fprintf(1,' \t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n',[iter,abs(ssr2-ssr1),abs(ssr2-ssr1)/abs(ssr2),SSR(iter),mean(mk)]);
    end
    
    A = A1;
    
    iter = iter + 1;
    
end

end
