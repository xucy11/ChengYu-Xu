
function [y,r,vr, percent]=SSA_function(x1,L,I)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SSA generates a trayectory matrix X from the original series x1
% by sliding a window of length L. The trayectory matrix is aproximated 
% using Singular Value Decomposition���ֽ⣩. The last step reconstructs
% the series from the aproximated trayectory matrix. 
%The SSA applications include smoothing, filtering, and trend extraction.
%SSA��Ӧ�ð���ƽ���������Լ����Ʒ���

% x1 Original time series (column vector form) ԭʼʱ������
% L  Window length                             ���ڳ���
% y  Reconstructed time series                 �ع�ʱ������
% r  Residual time series r=x1-y               ����ʱ������
% vr Relative value of the norm of the approximated trajectory matrix with respect to the original trajectory matrix

% The program output is the Singular Spectrum of x1 (must be a column vector),������������
% using a window length L. You must choose the components be used to reconstruct 
%the series in the form [i1,i2:ik,...,iL], based on the Singular Spectrum appearance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %I=('Choose the agrupation of components to reconstruct the series in the form I=[i1,i2:ik,...,iL]  ')
N=length(x1); 
   if L>N/2;L=N-L;end
   
	K=N-L+1; %Ƕ��
   X=zeros(L,K);  
	for i=1:K
	  X(1:L,i)=x1(i:L+i-1); %�õ��켣����
	end
    
% Step 2: SVD

   S=X*X'; 
	[U,autoval]=eig(S);%U������������autoval������ֵ
	[d,i]=sort(-diag(autoval));  
   d=-d;
   U=U(:,i);sev=sum(d); 
	%plot((d./sev)*100),hold on,plot((d./sev)*100,'rx');
    percent=sum(  d(I)./sev  );
	%title('�����׷���');xlabel('Eigenvalue Number');ylabel('Eigenvalue (% Norm of trajectory matrix retained)')
   V=(X')*U; 
   rc=U*V';

% Step 3: Grouping

   Vt=V';
   rca=U(:,I)*Vt(I,:);

% Step 4: Reconstruction

   y=zeros(N,1);  
   Lp=min(L,K);
   Kp=max(L,K);

   for k=0:Lp-2
     for m=1:k+1
      y(k+1)=y(k+1)+(1/(k+1))*rca(m,k-m+2);
     end
   end

   for k=Lp-1:Kp-1
     for m=1:Lp
      y(k+1)=y(k+1)+(1/(Lp))*rca(m,k-m+2);
     end
   end

   for k=Kp:N
      for m=k-Kp+2:N-Kp+1
       y(k+1)=y(k+1)+(1/(N-k))*rca(m,k-m+2);
     end
   end

   %figure;hold on;xlabel('Data plot');ylabel('Original and reconstructed series')
 % plot(y,'k'); grid on;plot(x1,'b--');%�ع���ԭʼͼ��
   r=x1-y;
  %figure;plot(r,'g');xlabel('Data plot');ylabel('Residual series');grid on%�в�����ͼ��
   vr=(sum(d(I))/sev)*100;
end