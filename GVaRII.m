clear all 
clc
tic
%% This code is used to achieve the empirical analysis in paper"Autoregressive models
% of the time series under volatility uncertainty and application to VaR
% model".
%
sp0=load('sp0.txt');% the data of S&P500 from Jan. 4,2010--July, 17,2020
n1=length(sp0);
K=5;
L=10;

P=100;
sp=sp0(n1-L-2499-K-P:n1);% 4300 for Table 3
sp1=diff(sp);
n=length(sp1);
sp2=100*sp1./sp(1:n);
alpha=0.05;

m=n-L+1;
m1=m-K+1;
m2=m1-P;
mean1=zeros(m,1);
var1=zeros(m,1);
uns=zeros(1,m1);
unsm=zeros(1,m1);
unsl=zeros(m1,1);
means=zeros(m1,1);

for i=1:m
 mean1(i)=mean(sp2(i:i+L-1));
 var1(i)=var(sp2(i:i+L-1));
end
for j=1:m1
 uns(j)=min(var1(j:j+K-1))/max(var1(j:j+K-1));   
 unsm(j)=max(var1(j:j+K-1));
 unsl(j)=min(var1(j:j+K-1));
 uns(j)= unsl(j);
 means(j)=mean1(j+K-1);
end
x1=ones(P,1);
Q=0;
meanss=zeros(m2-1,1);
vars=zeros(m2-1,1);
ass=zeros(m2-1,1);
varsk=zeros(m2-1,1);
for p=1:m2-1
A=[x1 uns(p:p+P-1)'];
beta=inv(A'*A)*A'*uns(p+1:p+P)';
A1=[x1 unsm(p:p+P-1)'];
beta1=inv(A1'*A1)*A1'*unsm(p+1:p+P)';
A2=[x1 means(p:p+P-1)];
beta2=inv(A2'*A2)*A2'*means(p+1:p+P);
y=beta(1)+beta(2)*uns(p+P);%measure the undertainty 
y1=beta1(1)+beta1(2)*unsm(p+P);%measure the maximum volatility
y2=beta2(1)+beta2(2)*means(p+P);%measure the return 
ass(p)=sp2(K+L+P+p-1);%the rate of return
vars(p)=y1;% the forecast of maximum variance
meanss(p)=y2;% the forecast of rate of return
varsk(p)=norminv((1/2+sqrt(y)/sqrt(y1)/2)*alpha,y2,sqrt(y1));%the forecast VaR
aa(p)=y;
if ass(p)<varsk(p)
Q=Q+1; % the violation of VaR
end 
line(p)=Q/p;
end
 
[y1 y y2]
%% Testing the G-VaR model Luc and Lcc
 n1=m2-1-Q;
 n2=Q;
 bin=2*(n1*log(n1/(n1+n2)/(1-alpha))+n2*log(n2/(n1+n2)/alpha));
 h=1-chi2cdf(bin, 1);
 
 u01=0;
 u00=0;
 u10=0;
 u11=0;
 for p=1:m2-2
if ass(p)<varsk(p) && ass(p+1)<varsk(p+1)
u11=u11+1;    
end 
if ass(p)<varsk(p) && ass(p+1)>varsk(p+1)
u10=u10+1;    
end   
if ass(p)>varsk(p) && ass(p+1)>varsk(p+1)
u00=u00+1;    
end 
if ass(p)>varsk(p) && ass(p+1)<varsk(p+1)
u01=u01+1;    
end         
 end
 pi=(u01+u11)/(u11+u10+u00+u01);
 pi01=u01/(u00+u01);
 pi11=u11/(u10+u11);
 di=m2-2-(u11+u10+u00+u01);
 if u11==0
 bin1=0;  
 else
 bin1=2*(u00*log(1-pi01)+u01*log(pi01)+u10*log(1-pi11)+u11*log(pi11)-(u00+u10)*log(1-pi)-(u01+u11)*log(pi));
 end
 [u11 u10 u00 u01 m2-1]
 h1=1-chi2cdf(bin1, 1);
 [Q/(m2-1) (m2-1)*alpha   h h1 -mean(varsk)]
 
%% Figure 1 mm1
%  subplot(211);
%  plot(unsm(m1-2499:m1),'b','LineWidth',1.1)
%  hold on
%  plot(unsl(m1-2499:m1),'r','LineWidth',1.1)
%  grid on
% legend('Maximum variance','Minimum variance');
% xlabel('From July 2010 to July 2020')
% 
% subplot(212);
%  plot(unsm(m1-249:m1),'b','LineWidth',1.3)
%  hold on
%  plot(unsl(m1-249:m1),'r','LineWidth',1.3)
%  grid on
% legend('Maximum variance','Minimum variance');
% xlabel('From July 2019 to July 2020')

%% Figure 2 re1
% subplot(311);
% plot(unsm(m1-P-1:m1-2),unsm(m1-P:m1-1),'.')
% y1=beta1(1)+beta1(2).*unsm(m1-P-1:m1-2);%measure the maximum volatility
% hold on
% plot(unsm(m1-P-1:m1-2),y1,'r','LineWidth',1.5)
% 
% legend('Autoregressive relation');
% ylabel('Maximum variance at t')
% xlabel(' Maximum variance at t-1')
% x=unsm(m1-P-1:m1-2);
% y=unsm(m1-P:m1-1);
% Rfenzi1=sum((x-mean(x)).*(y-mean(y)));
% Rfenmu1=sqrt(sum((x-mean(x)).^2).*sum((y-mean(y)).^2));
% R1=Rfenzi1/Rfenmu1; 
% 
% subplot(312);
% plot(uns(m1-P-1:m1-2),uns(m1-P:m1-1),'.')
% y=beta(1)+beta(2).*uns(m1-P-1:m1-2);%measure the undertainty 
% hold on
% plot(uns(m1-P-1:m1-2),y,'r','LineWidth',1.5)
% legend('Autoregressive relation');
% ylabel('Minimum variance at t')
% xlabel('Minimum variance at t-1')
% x=uns(m1-P-1:m1-2);
% y=uns(m1-P:m1-1);
% Rfenzi1=sum((x-mean(x)).*(y-mean(y)));
% Rfenmu1=sqrt(sum((x-mean(x)).^2).*sum((y-mean(y)).^2));
% R2=Rfenzi1/Rfenmu1; 
% 
% 
% subplot(313);
% plot(means(m1-P-1:m1-2),means(m1-P:m1-1),'.')
% hold on
% y2=beta2(1)+beta2(2).*means(m1-P-1:m1-2);%measure the return 
% plot(means(m1-P-1:m1-2),y2,'r','LineWidth',1.5)
% legend('Autoregressive relation');
% ylabel('Return at t')
% xlabel('Return at t-1')
% 
% x=means(m1-P-1:m1-2);
% y=means(m1-P:m1-1);
% Rfenzi1=sum((x-mean(x)).*(y-mean(y)));
% Rfenmu1=sqrt(sum((x-mean(x)).^2).*sum((y-mean(y)).^2));
% R3=Rfenzi1/Rfenmu1; 
% [beta1 beta  beta2]
% [R1 R2 R3]
% [R1 R2 R3].*[R1 R2 R3]
%% Figure 3 var1 

%  subplot(211);
%  plot(ass)
%  hold on
%  plot(varsk,'r','LineWidth',1.3)
%  grid on
% legend('log-return','-G-VaR');
% xlabel('From July 2010 to July 2020')
% 
%  subplot(212);
%   plot(ass(m2-1-249:m2-1))
%   hold on
%  plot(varsk(m2-1-249:m2-1),'r','LineWidth',1.3)
%  grid on
%  legend('log-return','-G-VaR');
% 
% xlabel('From July 2019 to July 2020')
% 

%% Figure 4  ra1

%  subplot(211);
%  plot(line,'b','LineWidth',1.5)
%  grid on
%  legend('The violations');
% 
%  xlabel('From July 2010 to July 2020')
% 
%  subplot(212);
%  plot(line(m2-1-249:m2-1),'b','LineWidth',1.5)
%  grid on
%  legend('The violations')
%  xlabel('From July 2019 to July 2020')
