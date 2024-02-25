%% 时域上
nsfactor = 1e-9;
sampled_time = 1.56250000000000e-12 / nsfactor;
dt = -1000 : sampled_time : 1000; % 时间间隔
td = [86, 26.4, 14, 4.6]; % typical drift time
%% 4.6ns
load('4.6ns');
Input1 = InputPort1.Sampled.Signal(1, [1:4098]);
Input2 = InputPort2.Sampled.Signal(2, [1:4098]);
time = InputPort2.Sampled.Time([1:4098]);
sampled_time = (time(2)-time(1));
given_t = 0;
fx = abs(fft(Input1));
fy = abs(fft(Input2));
for i = 1 : length(Input2)
    given_t = given_t + sampled_time;
    T = transferMatrix(fx, fy, time, given_t, 0, sampled_time, i);
    temp = T * [Input1(i); Input2(i)];
    Input1(i) = temp(1);
    Input2(i) = temp(2);
    i
end
%% 
Ix = Input1;
Iy = Input2;
E_Delta_tau_squared = calculate_EDeltaTauSquared(Ix, Iy, abs(fx(2)-fx(1)));
S4 = computeStokes(Ix, Iy);
seqlength = ceil(length(Ix)/2);
[c,lags] = xcorr(S4(1:seqlength), 'normalized');
plot(lags,c/E_Delta_tau_squared)
%% 画图
figure
hold on
for i = 1 : 4
ACF_PMDoft = ( 1-exp( -abs(dt)/td(i) ) ) ./ (abs(dt)/td(i));
plot(dt, ACF_PMDoft);
title('ACF_{PMD}(\Deltat)')
xlabel('\Deltat')
ylabel('ACF')
end
legend(['typical drift time=',num2str(td(1)),'ns'] ...
    , ['typical drift time=',num2str(td(2)),'ns'] ...
    , ['typical drift time=',num2str(td(3)),'ns'] ...
    , ['typical drift time=',num2str(td(4)),'ns'])
hold off
%% 频域上
psToGhzfactor = 1e-6;
mean_DGD = [1.5, 12.1] * psToGhzfactor;
df = -800: 0.1 :800;
dw = 2*pi * df;
figure
hold on
for i = 1 : 2
ACF_PMDoft = ( 1 - exp( -mean_DGD(i) * dw.^2/3 ) ) ./ ( dw.^2/3 * mean_DGD(i) );
wc(i) = sqrt( 3 /mean_DGD(i) );
fc(i) = wc(i) / (2*pi);
ACFc(i) =  ( 1 - exp( -mean_DGD(i) * wc(i).^2/3 ) ) ./ ( wc(i).^2/3 * mean_DGD(i) );
plot(df, ACF_PMDoft);
plot(fc(i), ACFc(i), 'r*');
line([fc(i), fc(i)],[ACFc(i), 0], 'LineStyle','--')
title('ACF_{PMD}(\Deltaf)')
ylabel('ACF')
xlabel('\Deltaf')
end
legend(['mean DGD=',num2str(mean_DGD(1)/psToGhzfactor),'ps'] ...
    , '', '' , ['mean DGD=',num2str(mean_DGD(2)/psToGhzfactor),'ps'])
hold off
%% 计算Stokes空间向量
function S = computeStokes(X, Y)
    S = zeros(4, length(X));
    % 计算Stokes参数
    S(1, :) = abs(X).^2 + abs(Y).^2;
    S(2, :) = abs(X).^2 - abs(Y).^2;
    S(3, :) = 2*real(X.*conj(Y));
    S(4, :) = 2*imag(X.*conj(Y));
end
%% 计算归一化参数
function E_Delta_tau_squared = calculate_EDeltaTauSquared(X, Y, fs)

    % 计算群时延的差值
    Delta_tau = computeDGD(X, Y, fs);

    % 计算Δτ的平方
    Delta_tau_squared = Delta_tau.^2;

    % 计算E[Δτ^2]
    E_Delta_tau_squared = mean(Delta_tau_squared);
end
