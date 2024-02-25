function [B, H] = JonesMatrix(wavelength, t)
    %% 系统参数
    N = 12;
    c = 299792458;
    f = c./ wavelength;
    w_total = 1.5e6; % rad/s
    mean_DGD = 3e-12; % ps
    %% 计算B, H矩阵
    B = cell(1, N);
    for i = 1 : N
        B{i} = BN(f, i, mean_DGD);
    end
    H = cell(1, N+1);
    for i = 1 : N+1
        H{i} = HN(t, w_total);
    end
end
%% 函数
function B = BN(f, n, mean_DGD)
    % f is the independent variable
    % n is the order, which denotes the ith B matrix
    sigma = 0.5; % controls randmness
    tauN = sqrt( 3*pi / (8*n) ) * mean_DGD * ( 1 + sigma*randn() );
    B = zeros(2, 2*length(f));
    for i = 1 : length(f)
        B(1, 2*i-1) = exp(1j*pi*f(i)*tauN);
        B(2, 2*i) = exp(-1j*pi*f(i)*tauN);
    end
end

function H = HN(t, w_total)
    % t is the independent variable
    % initial value
    e_kappa = rand()*2*pi;
    e_alpha = rand()*2*pi;
    e_phi = rand()*2*pi;
    % rotation speeds
    w_kappa = w_total ./ sqrt(3) ./ 2 + randn()*w_total/10;
    w_alpha = w_total ./ sqrt(3) + randn()*w_total/10;
    w_phi = w_total ./ sqrt(3) + randn()*w_total/10;
    % abs(((4*w_kappa^2 + w_phi^2 + w_alpha^2)^0.5 - w_total)/w_total*100) 误差
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