function dgd = computeDGD(signal1, signal2, fs)
    % 傅立叶变换
    FFT1 = fft(signal1);
    FFT2 = fft(signal2);

    % 计算相位响应
    phase1 = angle(FFT1);
    phase2 = angle(FFT2);

    % 计算频率轴
    f = (0:length(signal1)-1)*fs/length(signal1);

    % 计算群延迟
    groupDelay1 = -diff(phase1) ./ diff(f);
    groupDelay2 = -diff(phase2) ./ diff(f);

    % 计算差分群延迟
    dgd = abs(groupDelay1 - groupDelay2);
end
