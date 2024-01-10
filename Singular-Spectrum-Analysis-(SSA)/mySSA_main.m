%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%By ChengYu Xu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
data=readmatrix('D:\Cygwin64\home\11605\CRU193401-202012_de.csv');
%�Լ�ѡ���������������û�е���ӵ�switch��֧��
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

%��γ�ȷ�Χ����
switch (data_num)
    case 1
        lat=(-89.5:89.5);lon=(-179.5:179.5);%ȫ��1*1
    case 2
        lat=(-89:2:89);lon=(1:2:359);%ȫ��2*2
    case 3
        lat=(-89.75:0.5:89.75);lon=(-179.75:0.5:179.75);%ȫ��0.5*0.5
    case 4
        lat=(-87.5:5:87.5);lon=(2.5:5:357.5);%ȫ��5*5
    case 5
        lat=(-89:2:89);lon=(-179:2:179);%ȫ��2*2
    case 6
        lat=(-89:2:89);lon=(-179:2:179);%ȫ��2*2
    case 7
        lat=(-88.75:2.5:88.75);lon=(1.25:2.5:358.75);%ȫ��2.5*2.5
end

%ת����״����1ά��3ά
mydata=permute(reshape(data,[lonn,latn,yearn]),[3,2,1]);

%�Ƿ�ߵ�γ��
if data_num== 1 || data_num== 2
    mydata=mydata(:,end:-1:1,:);
end

%�Ƿ���matlab�����ƽ��һ����ǰ��������ɣ�
cal_anno=0;
switch (cal_anno)
    case 0
        %�����ƽ
        de=mydata;
    case 1
        %���ƽ
        de=zeros(yearn,latn,lonn);
        for jj=1:latn
            for kk=1:lonn
                if abs(mydata(1,jj,kk))>=abs(Fillvalue)
                    continue %����ȱ��
                end
                de(:,jj,kk)=mydata(:,jj,kk)-mean(mydata(:,jj,kk));
            end
        end
end

%��ÿ��������SSA��ȡǰPN�������ع����У�����������M���ع�������Ϊ����
%�ø���Ҷ�仯����ģ̬���ڣ�����һ�����ڵĽ���ϳɡ�
%����ռ��=PN���ع�������ĳ������ʵķ���/M���ع����е��ܷ���
L=yearn;
Fs=1;
PN=7;
M=15;

period1=10; %��1��ΪX<period1������
period2=20; %��2��Ϊperiod1<=X<period2������
period3=50; %��3��Ϊperiod2<=X<period3������
            %��4��ΪX>=period3������

test_level=0.9; %���ؿ�������������ˮƽ

%%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%
y=zeros(latn,lonn,PN,yearn);
period=zeros(1,PN);%���ÿ��ģ̬�������������
interanual_var=zeros(latn,lonn);
decade_var=zeros(latn,lonn);
mutidecade_var=zeros(latn,lonn);
superdecade_var=zeros(latn,lonn);

for jj=1:latn
    for kk=1:lonn
        if  ( abs( de(1,jj,kk) )>=abs(Fillvalue))
            continue %����ȱ��
        end
        for tt=1:PN
            [y(jj,kk,tt,:),~,~,~]=SSA_function(squeeze(de(:,jj,kk)),M,tt);
            %����Ҷ�仯�ж�����
            Y=squeeze(fft(y(jj,kk,tt,:)));
            P2 = abs(Y/L);
            P1 = P2(1:floor(L/2)+1);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Fs*(0:(L/2))/L;
            p=1./f;
            [~,index]=max(P1);
            period(tt)=p(index);
        end
        %�ϳ���ͬ����
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

%���ؿ������
[interanual_var_test,decade_var_test,mutidecade_var_test,superdecade_var_test] = SSA_monte_carlo_test(test_level,yearn,PN,M,period1,period2,period3);

%��NaN�滻��0����ʾ������Ϊ����ռ��Ϊ0
interanual_var(isnan(interanual_var)) = 0;
decade_var(isnan(decade_var)) = 0;
mutidecade_var(isnan(mutidecade_var)) = 0;
superdecade_var(isnan(superdecade_var)) = 0;

%���ȱ��Ϊ-1
interanual_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
decade_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
mutidecade_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
superdecade_var(abs(de(1,:,:))>=abs(Fillvalue))=-1;
%%%%%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%

%���
ssa=zeros(3,latn,lonn);
ssa(1,:,:)=decade_var;ssa(2,:,:)=mutidecade_var;ssa(3,:,:)=superdecade_var;
fid=fopen('D:\Cygwin64\home\11605\xxx.grd','wb');%�Զ���������д�뷽ʽ���ļ�
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