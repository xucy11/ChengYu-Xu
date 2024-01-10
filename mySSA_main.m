%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%By ChengYu Xu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
data=readmatrix('D:\Cygwin64\home\11605\CRU193401-202012_de.csv');
%自己选定数据输入参数，没有的添加到switch分支中
data_num=  6   ;%  1---GPCC 2--GISS 3--cru0.5*0.5 4--NOAA 5--GPCC2*2 6--CRU2*2 7--NOAA2.5*2.5
switch (data_num)
    case 1
        latn=180;
        lonn=360;
        yearn=70;
        Fillvalue= 9.96921e+36;
    case 2
        latn=90;
        lonn=180;
        yearn=70;
        Fillvalue= 9999;
    case 3
        latn=360;
        lonn=720;
        yearn=117;
        Fillvalue=1000000;
    case 4
        latn=36;
        lonn=72;
        yearn=70;
        Fillvalue=999;
    case 5
        latn=90;
        lonn=180;
        yearn=87;
        Fillvalue= -9.96921e+36;
    case 6
        latn=90;
        lonn=180;
        yearn=87;
        Fillvalue=9.96921e+36;
    case 7
        latn=72;
        lonn=144;
        yearn=73;
        Fillvalue=-9.96921e+36;
end

%经纬度范围设置
switch (data_num)
    case 1
        lat=(-89.5:89.5);lon=(-179.5:179.5);%全球1*1
    case 2
        lat=(-89:2:89);lon=(1:2:359);%全球2*2
    case 3
        lat=(-89.75:0.5:89.75);lon=(-179.75:0.5:179.75);%全球0.5*0.5
    case 4
        lat=(-87.5:5:87.5);lon=(2.5:5:357.5);%全球5*5
    case 5
        lat=(-89:2:89);lon=(-179:2:179);%全球2*2
    case 6
        lat=(-89:2:89);lon=(-179:2:179);%全球2*2
    case 7
        lat=(-88.75:2.5:88.75);lon=(1.25:2.5:358.75);%全球2.5*2.5
end

%转换形状，从1维到3维
mydata=permute(reshape(data,[lonn,latn,yearn]),[3,2,1]);

%是否颠倒纬度
if data_num== 1 || data_num== 2
    mydata=mydata(:,end:-1:1,:);
end

%是否在matlab中求距平（一般在前处理中完成）
cal_anno=0;
switch (cal_anno)
    case 0
        %不求距平
        de=mydata;
    case 1
        %求距平
        de=zeros(yearn,latn,lonn);
        for jj=1:latn
            for kk=1:lonn
                if abs(mydata(1,jj,kk))>=abs(Fillvalue)
                    continue %跳过缺测
                end
                de(:,jj,kk)=mydata(:,jj,kk)-mean(mydata(:,jj,kk));
            end
        end
end

%对每个格点进行SSA，取前PN个特征重构序列，超过窗口数M的重构序列视为噪音
%用傅立叶变化计算模态周期，属于一类周期的将其合成。
%方差占比=PN个重构序列中某段年代际的方差/M个重构序列的总方差
L=yearn;
Fs=1;
PN=7;
M=15;

period1=10; %第1段为X<period1的周期
period2=20; %第2段为period1<=X<period2的周期
period3=50; %第3段为period2<=X<period3的周期
            %第4段为X>=period3的周期

test_level=0.9; %蒙特卡洛检验的显著性水平

%%%%%%%%%%%%%%%%%%%%%%主程序%%%%%%%%%%%%%%%%%%%%%%
y=zeros(latn,lonn,PN,yearn);
period=zeros(1,PN);%存放每个模态的最大特征周期
interanual_var=zeros(latn,lonn);
decade_var=zeros(latn,lonn);
mutidecade_var=zeros(latn,lonn);
superdecade_var=zeros(latn,lonn);

for jj=1:latn
    for kk=1:lonn
        if  ( abs( de(1,jj,kk) )>=abs(Fillvalue))
            continue %跳过缺测
        end
        for tt=1:PN
            [y(jj,kk,tt,:),~,~,~]=SSA_function(squeeze(de(:,jj,kk)),M,tt);
            %傅立叶变化判断周期
            Y=squeeze(fft(y(jj,kk,tt,:)));
            P2 = abs(Y/L);
            P1 = P2(1:floor(L/2)+1);
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

        [~,~,~,interanual_var(jj,kk)]=SSA_function(squeeze(de(:,jj,kk)),M,interanual);
        [~,~,~,decade_var(jj,kk)]=SSA_function(squeeze(de(:,jj,kk)),M,decade);
        [~,~,~,mutidecade_var(jj,kk)]=SSA_function(squeeze(de(:,jj,kk)),M,mutidecade);
        [~,~,~,superdecade_var(jj,kk)]=SSA_function(squeeze(de(:,jj,kk)),M,superdecade);
    end
end

%蒙特卡洛检验
[interanual_var_test,decade_var_test,mutidecade_var_test,superdecade_var_test] = SSA_monte_carlo_test(test_level,yearn,PN,M,period1,period2,period3);

%将NaN替换成0，表示该周期为方差占比为0
interanual_var(isnan(interanual_var)) = 0;
decade_var(isnan(decade_var)) = 0;
mutidecade_var(isnan(mutidecade_var)) = 0;
superdecade_var(isnan(superdecade_var)) = 0;

%标记缺测为-1
interanual_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
decade_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
mutidecade_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
superdecade_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
%%%%%%%%%%%%%%%%%%%%%%主程序%%%%%%%%%%%%%%%%%%%%%%

%输出
ssa=zeros(3,latn,lonn);
ssa(1,:,:)=decade_var;ssa(2,:,:)=mutidecade_var;ssa(3,:,:)=superdecade_var;
fid=fopen('D:\Cygwin64\home\11605\xxx.grd','wb');%以二进制数据写入方式打开文件
for i=1:3
    for j=1:latn
        for k=1:lonn
            fwrite(fid,ssa(i,j,k),'float');
        end
    end
end
fclose(fid);

disp (['Test level of < period1: ',num2str(interanual_var_test)])
disp (['Test level of period1 <= XX < period2: ',num2str(decade_var_test)]) 
disp (['Test level of period2 <= XX < period3 ',num2str(mutidecade_var_test)])
disp (['Test level of >= period3: ',num2str(superdecade_var_test)])