clc
clear all
close all
s=tf('s');
%=================>vorudiha<===============================================
kp=0;
ki=1.5 ; 
kd=2 ;

teta=0;

l0=ki*exp(-j*teta);
h0=ki*exp(-j*teta);

l1=0.8+(kp*exp(-j*teta));
h1=1.2+(kp*exp(-j*teta));

l2=kd*exp(-j*teta);
h2=kd*exp(-j*teta);

l3=0.3;
h3=0.7;

l4=0.6;
h4=1;

%H0=input('H0 ra vared konid\n');
H0=[1 3 4];
%==============================>system<====================================
[row col]=size(H0);
n=input('martabe ra vared konid\n');
pow=zeros(1,n+1);

for t=1:n
    
    pow(1,t+1)=input(strcat('alfa ',int2str(t),'\n'));
    
end


pow(1)=0;
alfa0=pow(1);
alfa1=pow(2);
alfa2=pow(3);
alfa3=pow(4);
alfa4=pow(5);

%====================>Rmin,Rmax<===============================
% dar ghesmat haye mR1min,B,D bayad taghirate lazem ra ba tavajoh be
% mojmueye H vared konid...
mF1=[(1/alfa1);(2/alfa2);(3/alfa3);(4/alfa4)];
    F1=max(mF1);
mF2=[(1/(alfa4-alfa3)) (2/(alfa4-alfa2)) (3/(alfa4-alfa1)) (4/(alfa4-0))];
    F2=max(mF2);
    E=min(abs(l0),abs(h0));
 
 A1=max(abs(l1),abs(h1));
 A2=max(abs(l2),abs(h2));
 A3=max(abs(l3),abs(h3));
 A4=max(abs(l4),abs(h4));
mR1min=(A1+A3+A4);
    R1min=min(1,((E/mR1min)^(1/alfa1)));
 
 B=abs(l1)+abs(h1)+abs(l3)+abs(h3)+abs(l4)+abs(h4);
    R2min=[E/(E+max(B))]^F1;
 
 C0=max(abs(l0),abs(h0));
 C1=max(abs(l1),abs(h1));
 C2=max(abs(l2),abs(h2));
 C3=max(abs(l3),abs(h3));
 CC4=min(abs(l4),abs(h4));
 C4=C0/CC4;
 C5=C1/CC4;
 C7=C2/CC4;
 C8=C3/CC4;
 C6=C4+C5+C7+C8;
 
    R1max=max(1,(C6)^(1/(alfa4-alfa3)));
 
 D=[abs(l0) abs(h0) abs(l1) abs(h1) abs(l3) abs(h3)];
    R2max=(1+(max(D)/min(abs(l4),abs(h4))))^F2;
 
    Rmin=max(R1min,R2min)
    Rmax=min(R1max,R2max)






%-----------------------> Hj(r),Hj-prim(c)    <----------------------


f0=floor((alfa0)/2);
f1=floor((alfa1)/2);
f2=floor((alfa2)/2);
f3=floor((alfa3)/2);
f4=floor((alfa4)/2);

frac_f0=((alfa0)/2)-f0;
frac_f1=((alfa1)/2)-f1;
frac_f2=((alfa2)/2)-f2;
frac_f3=((alfa3)/2)-f3;
frac_f4=((alfa4)/2)-f4;

% dakhele Z anhayi k namoayani nadarand ra hazf mikonim (frac haye ozv haye H0 r vared mikonim)
matris=[frac_f0 frac_f1 frac_f2 frac_f3 frac_f4];

for v=1:col
    
   U(1,v)=matris(1,(H0(v)+1)); 
    
end
U1=unique(U);
B0=sort(U1)
[sat sot]=size(B0)

Hj=repmat(0.01,size(H0));
Hj_prim=repmat(0.01,size(H0));

 
  for f=1:col   
     
     for k=1:n+1
      
      if k-1==H0(1,f)&&rem(floor(pow(k)/2),2)==0
          for r=1:sot
            if  B0(r)==[(pow(k)/2)-floor(pow(k)/2)]
      Hj(r)=k-1
      r
            end
          end
      end
      
      if k-1==H0(1,f)&&rem(floor(pow(k)/2),2)~=0
           for r=1:sot
            if  B0(r)==[(pow(k)/2)-floor(pow(k)/2)]
      Hj_prim(r)=k-1
      r
           end
          end
      end
     end
  end        
%===========================>kharitanov<===================================
    
L=[l0 l1 l2 l3 l4];
H=[h0 h1 h2 h3 h4];
% b=1;
%Hj_prim(2)=[];
%Hj_prim(2)=[];
Q=zeros(n,n+1);
Q1=zeros(n,n+1);
for a=1:sot
    
   for h=1:n+1
             
          b=1;
          while b<a&&b<sot
    
    
             if h-1 ==Hj(b)
             %  h==Hj_prim(1)  
    
             Q(a,h)=H(h);
             Q1(a,h)=L(h);
             
             end
                 b=b+1;
                 
          end
              
          b=a;
          while  a-1<b&&b<=sot
    
    
             if h-1==Hj_prim(b)  
    
             Q(a,h)=H(h);
             Q1(a,h)=L(h);
             
             end
             b=b+1;
          end   
             
           
           
             if Q(a,h)~=H(h)
                 
              Q(a,h)=L(h);
              Q1(a,h)=H(h);
             %b=b+1;
             end
 
    end
    
end
Q
Q1
%=============================>sefr ha<=================================

%agar system be forme commenshorit bood in ghesmat ra faal mikonim...

ps=[1 s s^2 s^3 s^4];
pss1=ps*Q(1,:)';
pss2=ps*Q(2,:)';
pss3=ps*Q(3,:)';
pss4=ps*Q(4,:)';

zs1=zero(pss1)
zs2=zero(pss2)
zs3=zero(pss3)
zs4=zero(pss4)


%=====================================================================================
%v=[Rmin:1:Rmax];
%format long
q_positive=[h0 h1 h2 h3 h4];
q_negetive=[l0 l1 l2 l3 l4];
q_d=[q_positive + q_negetive]';


%===============================>Hw<==============================================
% 
% %Hw=zeros(floor(Rmin),floor(Rmax));
%  baze=floor(Rmax-Rmin+1);
% A_prim=zeros(sot,baze) ;
% B_prim=zeros(sot,baze);
% for d=1:sot
%     w=floor(Rmin)+1;
% while w<floor(Rmax)+1
%     for i=w
%        
% 
% D0=[(q_d(1)+(q_d(2)*w.^(alfa1)*exp(j*(alfa1/2)*pi)))+(q_d(3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi)))+(q_d(4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi)))+(q_d(5)*(w.^(alfa4)*exp(j*(alfa4/2)*pi)))];
% 
%  A_prim(d,i)=imag((Q(d,1)*(exp(j*(-B0(d))*pi)))+(Q(d,2)*(w.^(alfa1)*exp(j*(alfa1/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,5)*(w.^(alfa4)*exp(j*(alfa4/2)*pi))*exp(j*(-B0(d))*pi)));
%  B_prim(d,i)=A_prim(d,i)-imag((D0)*exp(j*-B0(d)*pi));
% 
%    
% z1(d,i)=[A_prim(d,i)] ;
% z2(d,i)=[ B_prim(d,i)];
% Z=[z1; z2];
%  
%   Hww(d,i)=max(Z(:,i));   
%    
% 
% % i=i+1;   
% w=w+1;
%     end
% end
% 
% end
% for ba=1:baze-1
%     Hw(1,ba)=max(Hww(:,ba));
% 
% end
% figure
% plot(Hw),grid,title 'tabe test paydarie moghavem'
% xlabel('w')
% ylabel('H(w)')
% 
% Hw=zeros(floor(Rmin),floor(Rmax));


 baze=floor(Rmax-Rmin+1);
A_prim=zeros(sot,baze) ;
B_prim=zeros(sot,baze);
for d=1:sot
    w=floor(Rmin);
while w<floor(Rmax)+1
    for i=1:baze
       

%D0=[(q_d(1)+(q_d(2)*w.^(alfa1)*exp(j*(alfa1/2)*pi)))+(q_d(3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi)))+(q_d(4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi)))];
D0=[(q_d(1)+(q_d(2)*w.^(alfa1)*exp(j*(alfa1/2)*pi)))+(q_d(3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi)))+(q_d(4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi)))+(q_d(5)*(w.^(alfa4)*exp(j*(alfa4/2)*pi)))];
 
 %A_prim(d,i)=imag((Q(d,1)*(exp(j*(-B0(d))*pi)))+(Q(d,2)*(w.^(alfa1)*exp(j*(alfa1/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi))*exp(j*(-B0(d))*pi)));
 A_prim(d,i)=imag((Q(d,1)*(exp(j*(-B0(d))*pi)))+(Q(d,2)*(w.^(alfa1)*exp(j*(alfa1/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi))*exp(j*(-B0(d))*pi))+(Q(d,5)*(w.^(alfa4)*exp(j*(alfa4/2)*pi))*exp(j*(-B0(d))*pi)));
  
 B_prim(d,i)=A_prim(d,i)-imag((D0)*exp(j*-B0(d)*pi));

   
z1(d,i)=[A_prim(d,i)] ;
z2(d,i)=[ B_prim(d,i)];
Z=[z1; z2];


 % Hww(d,i)=max(Z(:,i));   
  

% i=i+1;   
w=w+0.1;
    end
end

end
% for ba=1:baze-1
%     Hw(1,ba)=max(Hww(:,ba));
% 
% end
for lp=1:baze
%  for ni=1:2*sot
     Hw(1,lp)=max(Z(:,lp));
%  end
end
figure
plot(Hw),grid,title 'tabe test paydarie moghavem'
xlabel('w')
ylabel('H(w)')

%===============================>rous<==========================================
%w ra vared konid:
 w=Rmin;
 figure
while w<floor(Rmax)+1


A1_prim=zeros(2*sot,1);
D0=[(q_d(1)+(q_d(2)*w.^(alfa1)*exp(j*(alfa1/2)*pi)))+(q_d(3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi)))+(q_d(4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi)))+(q_d(5)*(w.^(alfa4)*exp(j*(alfa4/2)*pi)))];

for x=1:sot

A1_prim(x,1)=(Q(x,1))+(Q(x,2)*(w.^(alfa1)*exp(j*(alfa1/2)*pi)))+(Q(x,3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi)))+(Q(x,4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi)))+(Q(x,5)*(w.^(alfa4)*exp(j*(alfa4/2)*pi)));
A1_prim(x+sot,1)=(Q1(x,1))+(Q1(x,2)*(w.^(alfa1)*exp(j*(alfa1/2)*pi)))+(Q1(x,3)*(w.^(alfa2)*exp(j*(alfa2/2)*pi)))+(Q1(x,4)*(w.^(alfa3)*exp(j*(alfa3/2)*pi)))+(Q1(x,5)*(w.^(alfa4)*exp(j*(alfa4/2)*pi)));


end
for xy=1:2*sot
if xy>1
plot([real(A1_prim(xy,1)) real(A1_prim(xy-1,1))],[imag(A1_prim(xy,1)) imag(A1_prim(xy-1,1))])
 hold on
end
end
plot([real(A1_prim(2*sot,1)) real(A1_prim(1,1))],[imag(A1_prim(2*sot,1)) imag(A1_prim(1,1))]), grid on
w=w+0.1;
end




%==================================>end<============================================

