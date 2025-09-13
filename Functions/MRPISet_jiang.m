function Fs = MRPISet_jiang(Ak, W, epsilon)
    % Computes an invariant approximation of the 
    % minimal robust positively invariant set for:
    % x^{+} = Ax + w with w \in W
    % according to Algorithm 1 in 'Invariant approximations of
    % the minimal robust positively invariant set' by Rakovic et al. 
    % Requires a matrix A, a Polytope W, and a tolerance 'epsilon'. 
    [nx,~] = size(Ak); 
    s = 0; 
    alpha = 1000;
    Ms = 1000;
    E = eye(nx);
    it = 0;
    while(alpha > epsilon/(epsilon + Ms))
        s = s+1;
        alpha = max(W.support(Ak^s*(W.A)')./W.b);
        mss = zeros(2*nx,1);
        for i = 1:s
            mss = mss+W.support([Ak^i, Ak^i]);
        end
        Ms = max(mss);
        it = it+1;
    end
    fprintf('mPRI approxiamtion: iteration times s = %d, alpha = %f\n.',s,alpha);
    Fs = W;
%     Fs_A = [eye(nx);           
%             -eye(nx)];
    Fs_A = generate_polytope(6,1);
    Fs_b_new = zeros(size(Fs_A,1),1);
    lastStr = '';
    for i =1:s-1
 
        % 当前进度比例
        progress = i/(s-1);
        % 进度条长度
        barLength = 30; 
        % 已完成部分
        numBars = round(progress*barLength);
        % 构造字符串
        bar = [repmat('=',1,numBars), repmat(' ',1,barLength-numBars)];
    
        % 本次要打印的内容
        str = sprintf('[%s] %3.0f%% (step %d/%d)', bar, progress*100, i, s-1);
    
        % 清除上一行
        fprintf(repmat('\b',1,length(lastStr)));
    
        % 打印新进度条
        fprintf('%s', str);
    
        % 记录这次的字符串，方便下次覆盖
        lastStr = str;
        Fs = Fs+Ak^i*W;
        vertices = Fs.V';
        for j=1:size(Fs_A,1)
            [~,idx] = max(Fs_A(j,:)*vertices);
            Fs_b_new_j = Fs_A(j,:)*vertices(:,idx);
            Fs_b_new(j,1) = Fs_b_new_j;
        end
        Fs = Polyhedron('A',Fs_A,'b',Fs_b_new);
    end
    fprintf('\n');
    Fs = (1/(1-alpha))*Fs;
end
