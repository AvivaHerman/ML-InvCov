function A = getDictionary(dic_type,n,redundancy)
matricesDir = 'Matrices\';
prefix = 'dic';
if strcmp(dic_type,'Rand_ill_cond')
    % Random ill:
    filename = [matricesDir,prefix,num2str(n),'_',num2str(redundancy*n),dic_type,'.mat'];
    if isempty(dir(filename))
        A = randn(n,redundancy*n);
        % manipulation of singular values
        [U,S,V] = svd(A,'econ');
        s = diag(S);
        s = s/s(1);
        p = floor(log(1e-10)/log(s(end)));
        s = s.^p;
        S(:,1:size(A,1)) = diag(s);
        A = U*S*V';
        clear U V;
        % % normalize the dictionary:
        D = 1./(sqrt(sum(A.^2))');
        m = size(A,2);
        A = A*spdiags(D,0,m,m);
        clear D;
        save(filename,'A');
    else
        A = load(filename);
        A = A.A;
    end
    % Random
elseif strcmp(dic_type,'Rand_plain')==1
    filename = [matricesDir,prefix,num2str(n),'_',num2str(redundancy*n),dic_type,'.mat'];
    if isempty(dir(filename))
        A = randn(n,redundancy*n);
        D = 1./(sqrt(sum(A.^2))');
        m = size(A,2);
        A = A*spdiags(D,0,m,m);
        clear D;
        save(filename,'A');
    else
        A = load(filename);
        A = A.A;
    end
elseif strcmp(dic_type,'Rand_blurred')==1
    filename = [matricesDir,prefix,num2str(n),'_',num2str(redundancy*n),dic_type,'.mat'];
    if isempty(dir(filename))
        A = randn(n,redundancy*n);
        A = filter2(ones(7,1)/7,A);
        D = 1./(sqrt(sum(A.^2))');
        m = size(A,2);
        A = A*spdiags(D,0,m,m);
        clear D;
        save(filename,'A');
    else
        A = load(filename);
        A = A.A;
    end
elseif strcmp(dic_type,'K3')==1
    filename = [matricesDir,prefix,num2str(n),'_',num2str(redundancy*n),'K3.mat'];
    if isempty(dir(filename))
        A = randn(n,redundancy*n);
        s = svd(A,'econ');
        K = ceil(0.55*size(A,2));
        v = A(:,K-1);
        v = v/norm(v);
        for k = K:size(A,2)
            nv = randn(size(v));
            nv = nv/norm(nv);
            A(:,k) = v + 0.05*nv;
        end
        [U,S,V] = svd(A,'econ');
        A = U*diag(s)*V';
        clear U V;
        % normalize the dictionary:
        D = 1./(sqrt(sum(A.^2))');
        m = size(A,2);
        A = A*spdiags(D,0,m,m);
        clear D;
       save(filename,'A');
    else
        A = load(filename);
        A = A.A;
    end
elseif strcmp(dic_type,'Rand_sign')==1
    filename = [matricesDir,prefix,num2str(n),'_',num2str(redundancy*n),dic_type,'.mat'];
    if isempty(dir(filename))
        A = sign(randn(n,redundancy*n)-0.5);
        % normalize the dictionary:
        D = 1./(sqrt(sum(A.^2))');
        m = size(A,2);
        A = A*spdiags(D,0,m,m);
        clear D;
        save(filename,'A');
    else
        A = load(filename);
        A = A.A;
    end
elseif strcmp(dic_type,'PartialHadamard')==1
    filename = [matricesDir,prefix,num2str(n),'_',num2str(redundancy*n),dic_type,'.mat'];
    if isempty(dir(filename))
        m = redundancy*n;
        H = hadamard(m);
        p = randperm(m);
        A = H(p(1:n), :);
        % normalize the dictionary:
        D = 1./(sqrt(sum(A.^2))');
        m = size(A,2);
        A = A*spdiags(D,0,m,m);
        clear D;
        save(filename,'A');
    else
        A = load(filename);
        A = A.A;
    end
    
elseif strcmp(dic_type,'PartialDCT')==1
    filename = [matricesDir,prefix,num2str(n),'_',num2str(redundancy*n),dic_type,'.mat'];
    if isempty(dir(filename))
        m = redundancy*n;
        H = dct(eye(m));
        p = randperm(m);
        A = H(p(1:n), :);
        % normalize the dictionary:
        D = 1./(sqrt(sum(A.^2))');
        m = size(A,2);
        A = A*spdiags(D,0,m,m);
        clear D;
        save(filename,'A');
    else
        A = load(filename);
        A = A.A;
    end    
end
return;