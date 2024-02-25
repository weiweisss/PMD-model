function T = transferMatrix(f, t, given_t, given_f)
    %% 系统参数
    N = 12;
    mean_DGD = 3e-12; % ps
    delta_f = abs(f(2)-f(1));
    delta_t = t(2) - t(1);
    c = 299792458;
    %% 计算T矩阵
    
    Bmatrices = cell(1, N);
    for i = 1 : N
        Bmatrices{i} = BN(f, i, mean_DGD);
    end
    Hmatrices = cell(1, N+1);
    for i = 1 : N+1
        Hmatrices{i} = HN(t);
    end
        
    tindex = ceil(given_t / delta_t);
%   findex = ceil(given_f / delta_f);
    findex = ceil(c/given_f*1e10 - 15300);
    T = Hmatrices{N+1}(:, 2*tindex-1: 2*tindex);
    for i = N : -1 : 1
        T = T * Bmatrices{i}(:, 2*findex-1: 2*findex);
        T = T * Hmatrices{i}(:, 2*tindex-1: 2*tindex);
    end
end
%% 函数
function B = BN(f, n, mean_DGD)
    % f is the independent variable
    % n is the order, which denotes the ith B matrix
    sigma = 1; % controls randmness
    tauN = ( 3*pi / (8*n) ) * mean_DGD * ( 1 + sigma*randn() );
    B = zeros(2, 2*length(f));
    for i = 1 : length(f)
        B(1, 2*i-1) = exp(1j*pi*f(i)*tauN);
        B(2, 2*i) = -exp(1j*pi*f(i)*tauN);
    end
end

function H = HN(t)
    % t is the independent variable
    % initial value
    e_kappa = rand()*2*pi;
    e_alpha = rand()*2*pi;
    e_phi = rand()*2*pi;
    % rotation speeds
    w_total = 1.5e6; % rad/s
    w_kappa = w_total ./ sqrt(3) + randn();
    w_alpha = w_total ./ sqrt(3) + randn();
    w_phi = w_total ./ sqrt(3) + randn();
    % angles
    kappa = w_kappa*t + e_kappa;
    alpha = w_alpha*t + e_alpha;
    phi = w_phi*t + e_phi;
    H = zeros(2, 2*length(t));
    for i = 1 : length(t)
        H(1, 2*i-1) = cos(kappa(i)) * exp(1j*alpha(i));
        H(1, 2*i) = -sin(kappa(i)) * exp(1j*phi(i));
        H(2, 2*i-1) = sin(kappa(i)) * exp(-1j*phi(i));
        H(2, 2*i) = cos(kappa(i)) * exp(-1j*alpha(i));
    end
end