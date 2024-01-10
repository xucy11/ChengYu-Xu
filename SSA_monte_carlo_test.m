function [ interanual_var,decade_var,mutidecade_var,superdecade_var ] = SSA_monte_carlo_test( alpha , L, PN, M, period1, period2, period3)
%Monte Carlo检验
%alpha置信度
%输出各尺度的显著阈值

%判断采样点个数是否是奇数
if mod(L,2)==1
    L=L+1;
end

randn=10000;
noise=randn(randn,L);
Fs=1;
trust_percent=alpha*randn;
period=zeros(1,PN);%存放每个模态的最大特征周期
interanual_var=zeros(1,randn);
decade_var=zeros(1,randn);
mutidecade_var=zeros(1,randn);
superdecade_var=zeros(1,randn);

for ii=1:randn
    for tt=1:PN
        [y,~,~,~]=SSA(noise(ii,:)',M,tt);
        %傅立叶变化判断周期
        Y=squeeze(fft(y));
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        p=1./f;
        [~,index]=max(P1);
        period(tt)=p(index);
    end
    %合成相同周期
    interanual=find(period<period1);
    decade=find(period>=period1&period<period2);
    mutidecade=find(period>=period2&period<period3);
    superdecade=find(period>=period3);
    [~,~,~,interanual_var(ii)]=SSA(noise(ii,:)',M,interanual);
    [~,~,~,decade_var(ii)]=SSA(noise(ii,:)',M,decade);
    [~,~,~,mutidecade_var(ii)]=SSA(noise(ii,:)',M,mutidecade);
    [~,~,~,superdecade_var(ii)]=SSA(noise(ii,:)',M,superdecade);
end
interanual_var=sort(interanual_var);interanual_var=interanual_var(trust_percent);
decade_var=sort(decade_var);decade_var=decade_var(trust_percent);
mutidecade_var=sort(mutidecade_var);mutidecade_var=mutidecade_var(trust_percent);
superdecade_var=sort(superdecade_var);superdecade_var=superdecade_var(trust_percent);
end

