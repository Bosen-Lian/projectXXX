clear all;
close all;
clc;


%系统
A=[0.906488 0.0816012 -0.0005;
    0.074349 0.90121 -0.000708383;
    0 0 0.132655];

B=[-0.00150808;
    -0.0096;
    0.867345];

D=[0.00951892;
    0.00038373;
    0];



%专家参数
Qe=1*[1 0 0 ;
      0 1 0 ;
      0 0 1 ];
Re=1;
gammae=5;
Pe=dare(A,B,Qe,Re);
Ke=inv(Re+B'*Pe*B)*(B'*Pe*A);
% Le=inv(-gammae^2+D'*Pe*D-D'*Pe*B*inv(Re+B'*Pe*B)*B'*Pe*D)*(D'*Pe*A-D'*Pe*B*inv(Re+B'*Pe*B)*B'*Pe*A);



%学生参数
Q=0*[0.01 0 0 ;
      0 0.01 0 ;
      0 0 0.01];
R=1;
K1=[0.1    0.1    .11];
K=K1;




%迭代初值
x(1:3,1)=[1;-1;-1];
xx(1:3,1)=[1;1;-1];
k=0;
kk=0;
ii=0;
jj=0;
eK=[];
eL=[];
eQ=[];
stopq=0;
alpha=1; %0.75

for i=1:10000
    
    %专家数据
    pe(:,i)=0.00001*(rand(1));
    ue(:,i)=-Ke*x(:,i)+pe(:,i);


    x(:,i+1)=A*x(:,i)+B*ue(:,i);
    
      uue(:,i)=-Ke*xx(:,i);
    xx(:,i+1)=A*xx(:,i)+B*uue(:,i)+pe(:,i);  
    
    
    
    ii=ii+1;
    Hxx(ii,:)=[x(1,i)^2, 2*x(1,i)*x(2,i), 2*x(1,i)*x(3,i), x(2,i)^2, 2*x(2,i)*x(3,i),x(3,i)^2]...
              -[x(1,i+1)^2, 2*x(1,i+1)*x(2,i+1), 2*x(1,i+1)*x(3,i+1), x(2,i+1)^2, 2*x(2,i+1)*x(3,i+1),x(3,i+1)^2];%1-6  
    Hxu(ii,:)=2*kron(x(:,i)'*K1',x(:,i)')+2*kron(ue(:,i)',x(:,i)');%7-9
    Huu(ii,:)=kron(ue(:,i)-K1*x(:,i),ue(:,i)+K1*x(:,i));%10
    r(ii)=x(:,i)'*K1'*R*K1*x(:,i)+x(:,i)'*Q*x(:,i)...
        +alpha*(ue(:,i)+K1*x(:,i))'*R*(ue(:,i)+K1*x(:,i))+(ue(:,i)-K1*x(:,i))'*R*(ue(:,i)+K1*x(:,i));

       Hqxx(ii,:)=[xx(1,i)^2, 2*xx(1,i)*xx(2,i), 2*xx(1,i)*xx(3,i), xx(2,i)^2, 2*xx(2,i)*xx(3,i),xx(3,i)^2]...
              -[xx(1,i+1)^2, 2*xx(1,i+1)*xx(2,i+1), 2*xx(1,i+1)*xx(3,i+1), xx(2,i+1)^2, 2*xx(2,i+1)*xx(3,i+1),xx(3,i+1)^2];%1-6  
%         
        Hxuu(ii,:)=2*kron(xx(:,i)'*K1',xx(:,i)')+2*kron(uue(:,i)',xx(:,i)');%7-9  
        Huuu(ii,:)=kron(uue(:,i)-K1*xx(:,i),uue(:,i)+K1*xx(:,i));%10
        n(ii)=-xx(:,i)'*K1'*R*K1*x(:,i)-(uue(:,i)-K1*xx(:,i))'*R*(K1*xx(:,i)+uue(:,i));
%             n(ii)=xx(:,i)'*Q*xx(:,i)+alpha*(uue(:,i)+K1*xx(:,i))'*R*(uue(:,i)+K1*xx(:,i));

    
    
    
    jj=jj+1;
    if stopq==0
      
       X(jj,:)=[xx(1,i)^2, 2*xx(1,i)*xx(2,i), 2*xx(1,i)*xx(3,i), xx(2,i)^2, 2*xx(2,i)*xx(3,i),xx(3,i)^2];


%           
        Hqxxx(jj,:)=Hqxx(ii,:);
        Hqxu(jj,:)=Hxuu(ii,:);
% 
        Hquu(jj,:)=Huuu(ii,:);

         rq(jj)= n(ii);


    end
    
    rank(X)
    
    if rank(X)==6
        stopq=1;
    end
    
    zbar=[Hxx Hxu Huu];

    m=zbar; 

    a=rank(m)
        
    if a==10
        
       H=m\r';

       K1=inv(H(10))*([H(7) H(8) H(9)]);
       K=[K;K1];

       eK=[eK,norm(K1-Ke)];
  

        Hxx=[];
        Hxu=[];
        Hxw=[];
        Huu=[];
        Huw=[];
        Hww=[];
        zbar=[];
        r=[];
        ii=0;
        m=[];
        
        Q=Q+alpha*(Ke-K1)'*R*(Ke-K1);
           eQ=[eQ,norm(Q)];
%           Hq=Hqxxx*[H(1),H(2),H(3),H(4),H(5),H(6)]'...
%               +Hqxu*[H(7) H(8) H(9)]'+Hquu*H(10)+rq;
%  
% 
%          Hq=rq;
%         HQ=X\Hq';
%         Q=[HQ(1),HQ(2),HQ(3);
%            HQ(2),HQ(4),HQ(5);
%            HQ(3),HQ(5),HQ(6)];
        eQ=[eQ,norm(Q)];
        jj=0;
        stopq=0;
        X=[];
        Hqxxx=[];
        Hqxu=[];
        Hqxw=[];
        Hquu=[];
        Hquw=[];
        Hqww=[];
        Hq=[];
        
    end
    


end
 
figure(1)
plot(eQ,'-*','LineWidth',1)
legend('Q')

figure(2)
plot(eK,'-*','LineWidth',1)
legend('K-Ke')


