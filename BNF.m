clear;
clc;
close all;
%%

dt = 1e-9;                                % The time interval     ʱ����    dt
cnt= 1000;                                 % the time points      ʱ���
T = cnt * dt;           

freq =800e9;                             % Frequency range        Ƶ�ʷ�Χ
df=1e10;                              % Frequency resolution      Ƶ�ʷֱ���  df
length_f_rf = freq/df;                    % The number of Frequency points    Ƶ�� 80 ��
f_rf = freq/length_f_rf*(-length_f_rf/2:1:length_f_rf/2);  % frequency interval    Ƶ�ʼ��
omg_0 = 2*pi* 3e8 / (1550e-9);                    % The central optical frequency at 1550nm   ����Ƶ��
omg = 2*pi*f_rf + omg_0;     

len_omg = length(omg);


%%
%   DGD 
N = 12;%����
DGD =10e-12;       %%Channel average DGD                   ƽ��DGD
sigma = 0.5;      %% A parameter that determines the actual value of each segment of DGD
Xn = normrnd(0,1,[1,N]);
tao = sqrt(3*pi/(8*N))*DGD*(1+sigma*Xn);

tao(N+1) = 0;

%%

%%����ʱ���ݻ�
t = 0:1:999;%ʱ��ȡ1000����
len_t = length(t);%1000
kepa = ones(N,len_t);%13*1000
alpha = ones(N,len_t);
phy = ones(N,len_t);%��ʼ��

%%
%13�� ��Ҫ13�����������һ��DGD��taoΪ0
   omg_kepa = 1.3*10^5;%��ʱȡ��ֵ���������ÿ�β�һ��
    omg_alpha = 1.2*10^5;
    omg_phy = 1.1*10^5; 
    
for i = 1:N+1
    e_kepa = rand*2*pi;
    e_alpha= rand*2*pi;
    e_phy= rand*2*pi; 
    
    kepa_end= cnt*omg_kepa*dt+ e_kepa;
    alpha_end= cnt*omg_alpha*dt + e_alpha;
    phy_end= cnt*omg_phy*dt + e_phy; 
    
    kapa = linspace(e_kepa,kepa_end, 1000);
    alpha = linspace(e_alpha,alpha_end, 1000);
    fai = linspace(e_phy,phy_end, 1000) ;
    
    u11=cos(kapa).*exp(-1j*alpha);
    u12=exp(1j*fai).*sin(kapa);
    u21= -exp(-1j*fai).* sin(kapa);
    u22= exp(1j*alpha).*cos(kapa); 
    
    U11(:,i)=u11;                                    %ÿ����������ĸ�Ԫ��
    U12(:,i)=u12;
    U21(:,i)=u21;
    U22(:,i)=u22; 
end

%%
%��ʼ��U���� PMD


%%
for j = 1:1000  
    UU11 = ones(1,len_omg);
UU12 = zeros(1,len_omg);
UU21 = zeros(1,len_omg);
UU22 = ones(1,len_omg);


    RSOP11 = U11(j,:);
    RSOP12 = U12(j,:);
    RSOP21 = U21(j,:);
    RSOP22 = U22(j,:);    
    for i = 1:N
        tao_temp = tao(i);
        dgd_cheng_rsop11 = RSOP11(i)*exp(1j*omg*tao_temp/2);%omg�����飬��֪���ܲ���ֱ���㡪���ܣ�ֱ������1*4����
        dgd_cheng_rsop12 = RSOP12(i)*exp(1j*omg*tao_temp/2);
        dgd_cheng_rsop21 = RSOP21(i)*exp(-1j*omg*tao_temp/2);
        dgd_cheng_rsop22 = RSOP22(i)*exp(-1j*omg*tao_temp/2);
        
        temp11 = dgd_cheng_rsop11.*UU11 + dgd_cheng_rsop12.*UU21;%ֱ�ӳ���Unn
        temp12 = dgd_cheng_rsop11.*UU12 + dgd_cheng_rsop12.*UU22;
        temp21 = dgd_cheng_rsop21.*UU11 + dgd_cheng_rsop22.*UU21;
        temp22 = dgd_cheng_rsop21.*UU12 + dgd_cheng_rsop22.*UU22;
        
        UU11 = temp11;%���Բ�����ʱֵ��_NO����ΪU11�ı��Ӱ�����������
        UU12 = temp12;
        UU21 = temp21;
        UU22 = temp22;
    end 
    U(1,1,:) = UU11;%����ؼ�������ֱ�Ӹ�ֵ�Ϳ��ԣ�����ҪUnn��i��
    U(1,2,:) = UU12;
    U(2,1,:) = UU21;
    U(2,2,:) = UU22;
    
    % ��Q���󣬽�����һ����PMD
    w = 1 : len_omg - 1;
    %��Uw:��Uw���ĸ�Ԫ����ֱ������Ƿֱ���Ui��΢��
    Uw11 = (UU11(w+1)-UU11(w))./(omg(w+1)-omg(w));
    Uw12 = (UU12(w+1)-UU12(w))./(omg(w+1)-omg(w));
    Uw21 = (UU21(w+1)-UU21(w))./(omg(w+1)-omg(w));
    Uw22 = (UU22(w+1)-UU22(w))./(omg(w+1)-omg(w));
    
    Uw(1,1,:) = Uw11;
    Uw(1,2,:) = Uw12;
    Uw(2,1,:) = Uw21;
    Uw(2,2,:) = Uw22;
    
    Q = zeros(2,2,length(Uw11));%300ά
    for ii = 1:length(Uw11)
        Q(:,:,ii) = 1j*Uw(:,:,ii)*U(:,:,ii)';
    end   
 %�ټ���ƫ��ģɫɢ����tao
    tao1 = real(1*(Q(1,1,:) - Q(2,2,:)));
    tao2 = real(1*(Q(1,2,:) + Q(2,1,:)));
    tao3 = real((Q(1,2,:) - Q(2,1,:))*1j);
    
    tao1_squz = squeeze(tao1).';%�Ȱ���άɾ��һά���ٰѾ���ת��Ϊ������
    tao2_squz = squeeze(tao2).';
    tao3_squz = squeeze(tao3).';
    
   delta_tao = sqrt(tao1.^2 + tao2.^2 + tao3.^2);  
   delta_tao_squz = squeeze(delta_tao).'; %�������Ҫ�Ħ�tao DGD
   delta_tao_squz = delta_tao_squz .*1e12;
   mat_delta_tao_squz(j,:) = delta_tao_squz;   %DGDֵ   
end


contourf(mat_delta_tao_squz);
temp = mat_delta_tao_squz(:);

mean_mat_delta_tao_squz = temp;    % DGDֵ����

%%  ��DGD�ֲ� �������˹Τ�ֲ�
figure
[rr_freq, x_centers] = hist(mean_mat_delta_tao_squz,50);       % [The number of statistics for each interval��The center coordinates]
lengthwin=((max(mean_mat_delta_tao_squz(:))-min(mean_mat_delta_tao_squz(:)))/length(rr_freq));
rr_pdf = rr_freq./sum(rr_freq(:))/lengthwin;       %rr_freq is the number of occurrences,sum(rr_pdf)*lengthwin=1
bar(x_centers,rr_pdf);
hold on
x = linspace(0,max(mean_mat_delta_tao_squz),length(rr_freq));
meantao = mean(mean_mat_delta_tao_squz);
y = 32*x.^2.*exp(-(4*x.^2)/(pi*meantao.^2))/(pi.^2*meantao.^3);  %% Theory,References: Polarization Optics in Telecommunications PP-402

y=y./sum(y)/lengthwin;
plot(x,y)