
%Initial Conditions

A=1/0.83428;  %1/sqrt(N);
V=2;
m=1;
hbar=1;
knot=(2*m*V/(hbar).^2).^(.5);
d=1;
sigmax=16;
sigmak=1/(2*sigmax);
deltak=(pi)*sigmak/(120*d);
q=1;
deltax=(6)^(1/3)/(10*(3*sigmak+q));    %intput for first block of code
t=0; 
xnot=-50;
deltat=1;
e=2.71828;
k=[(-3*sigmak-q):deltak:(3*sigmak-q)]; 
kt=(knot.^2-k.^2).^(.5);
w=((knot)^2)./(2*i*k.*kt.*coth(kt*d)+2*k.^2-knot^2);
% Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2)-i*k.^2*hbar*t/(2*m));
     %old gauss
     %Gauss=sigmax*e.^(-2*pi*i*k*xnot-2*pi^2*sigmax^2*k.^2-i*k.^2*hbar*t/(2*m)); 
N=0; 

%x=[-6*sigmax:deltax:0*sigmax];
% Gaussx=e.^-((x-xnot).^2/(2*sigmax^2))/(2*pi).^.5;
%plot(x,trapz(k,(abs(Gauss.*e.^(i*2*pi*k*x))).^2);

%Probability of finding particle between B and C
%{
B=-6;
C=0;
t=0;
SUMOPx1=0;
SUMOPx2=0;
SUMOPx3=0;

if all([B<=0,C<=0])   
    fprintf('integrating from I to I'); 
    x=B;
        while x<=C
            ya1=(A.*e.^(i*k*(x))+w.*A.*e.^(-i*k*(x)));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y1= ya1.*Gauss;
            SUMOPx1 = SUMOPx1+(deltax)*(abs(trapz(k,y1.*e.^(i*2*pi*k*x)))).^2;
            x=x+deltax;
        end
    SUMOPx1
end

if all([B>=0,B<d,C>0,C<=d])
    fprintf('integrating from II to II'); 
    x=B;
        while x<=C
            ya2=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*(x));
            yb2=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*(x));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y2=(ya2+yb2).*Gauss;
            SUMOPx2 = SUMOPx2+(deltax)*(abs(trapz(k,y2.*e.^(i*2*pi*k*x)))).^2; 
            x=x+deltax;
        end
    SUMOPx2
end

if all([B>d,C>d])
    fprintf('integrating from III to III'); 
    x=B;
        while x<=C
            ya3=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*d-i*k*d+i*k*(x));
            yb3=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*d-i*k*d+i*k*(x));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y3=(ya3+yb3).*Gauss;
            SUMOPx3 = SUMOPx3+(deltax)*(abs(trapz(k,y3.*e.^(i*2*pi*k*x)))).^2; 
            x=x+deltax;
        end
    SUMOPx3
end

if all([B<0,C>0,C<d])
    fprintf('integrating from I to II'); 
    x=B;
        while x<=0
            ya1=(A.*e.^(i*k*(x))+w.*A.*e.^(-i*k*(x)));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y1= ya1.*Gauss;
            SUMOPx1 = SUMOPx1+(deltax)*(abs(trapz(k,y1.*e.^(i*2*pi*k*x)))).^2;
            x=x+deltax;
        end
    x=0;
        while x<=C
            ya2=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*(x));
            yb2=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*(x));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y2=(ya2+yb2).*Gauss;
            SUMOPx2 = SUMOPx2+(deltax)*(abs(trapz(k,y2.*e.^(i*2*pi*k*x)))).^2; 
            x=x+deltax;
        end
    SUMOPx1+SUMOPx2 
end

if all([B<d,B>=0,C>d])
    fprintf('integrating from II to III'); 
    x=B;
        while x<=d
            ya2=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*(x));
            yb2=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*(x));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y2=(ya2+yb2).*Gauss;
            SUMOPx2 = SUMOPx2+(deltax)*(abs(trapz(k,y2.*e.^(i*2*pi*k*x)))).^2; 
            x=x+deltax;
        end
    x=d;
        while x<=C
            ya3=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*d-i*k*d+i*k*(x));
            yb3=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*d-i*k*d+i*k*(x));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y3=(ya3+yb3).*Gauss;
            SUMOPx3 = SUMOPx3+(deltax)*(abs(trapz(k,y3.*e.^(i*2*pi*k*x)))).^2; 
            x=x+deltax;
        end
    SUMOPx2+SUMOPx3
end
            
if all([B<0,C>d])
    fprintf('integrating from I to III'); 
    x=B;
        while x<=0
            ya1=(A.*e.^(i*k*(x))+w.*A.*e.^(-i*k*(x)));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y1= ya1.*Gauss;
            SUMOPx1 = SUMOPx1+(deltax)*(abs(trapz(k,y1.*e.^(i*2*pi*k*x)))).^2;
            x=x+deltax;
        end
    x=0;
        while x<=d
            ya2=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*(x));
            yb2=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*(x));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y2=(ya2+yb2).*Gauss;
            SUMOPx2 = SUMOPx2+(deltax)*(abs(trapz(k,y2.*e.^(i*2*pi*k*x)))).^2; 
            x=x+deltax;
        end
    x=d;
        while x<=C
            ya3=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*d-i*k*d+i*k*(x));
            yb3=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*d-i*k*d+i*k*(x));
Gauss=e.^(i*sigmax^2*(i*k.^2/2-xnot*k/sigmax^2));
            y3=(ya3+yb3).*Gauss;
            SUMOPx3 = SUMOPx3+(deltax)*(abs(trapz(k,y3.*e.^(i*2*pi*k*x)))).^2; 
            x=x+deltax;
        end
    SUMOPx1+SUMOPx2+SUMOPx3
end
%} 

%{
 %Normalization
 format long
 N=0;
 n=100;
 SUMOPx1=0;
 SUMOPx2=0;
 SUMOPx3=0;
 x=-n*d;
 p=0;
 
 while x<0
   
    ya1=(A*e.^(i*k*x)+w*A.*e.^(-i*k*x)).*e.^(i*k.^2*hbar*t/(2*m));
            
    %Gauss=((1/(2*pi))^.25)*e.^(-k.^2/4);  %for gaussian
    Gauss=e.^-((1/sigmak)^2*((k+q).^2/2-(xnot)*sigmak^2*i*k));
            y1= ya1.*Gauss; 
% SUMOPx1 = SUMOPx1+(deltax)*(abs(trapz(k,(trapz(k,y1)).*e.^(i*k*x)))).^2;
  SUMOPx1 = SUMOPx1+(deltax)*(abs(trapz(k,y1))).^2;
% SUMOPx1 = SUMOPx1+(deltax)*(abs(trapz(k,y1.*((1/(2*pi))^.5).*e.^(i*k*x)))).^2;  %only for gaussian
  x=x+deltax;
  if mod(x,2*n*d/100) <= deltax/(n*d/100) %give fraction of computational time
      round(x/(2*n*d/100)+50)
  end
 end
 
 x=0;
     
 while all([x<d,x>=0])
  
    
             
    ya2=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             yb2=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             Gauss=e.^-((1/sigmak)^2*((k+q).^2/2-(xnot)*sigmak^2*i*k));
             y2=(ya2+yb2).*Gauss;
  SUMOPx2 = SUMOPx2+(deltax)*(abs(trapz(k,y2))).^2;
  
  x=x+deltax;
  if mod(x,2*n*d/100) <= deltax/(n*d/100)
      round(x/(2*n*d/100)+50)
  end
 end
 

  x=d;
     
 while all([x<n*d,x>=d])
  
   ya3=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*d-i*k*d+i*k*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             yb3=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*d-i*k*d+i*k*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             Gauss=e.^-((1/sigmak)^2*((k+q).^2/2-(xnot)*sigmak^2*i*k));
             y3=(ya3+yb3).*Gauss;
 SUMOPx3 = SUMOPx3+(deltax)*(abs(trapz(k,y3))).^2;
 
  x=x+deltax;
  if mod(x,2*n*d/100) <= deltax/(n*d/100)
      round(x/(2*n*d/100)+50)
  end
 end
 
 N=SUMOPx1+SUMOPx2+SUMOPx3;
 (N)^.5
%}
% ktq=(knot.^2-q.^2).^(.5);
% wq=(-(knot)^2)./(-2*i*q.*ktq.*coth(ktq*d)-q.^2+ktq.^2);
% dwq=(2*i*knot^2*((knot^2-2*q^2)*coth(ktq*d)+q*ktq*(d*q*(csc(ktq*d))^2-2*i)))/(ktq*(2*q*ktq*coth(ktq*d)+i*(knot^2-2*q^2))^2);
% subplot(5,1,1)
% plot(k,w);
% subplot(5,1,2)
% plot(k,wq+dwq*(k-q));
% subplot(5,1,3)
% plot(k,i*w);
% subplot(5,1,4)
% plot(k,i*(wq+dwq*(k-q)));
% subplot(5,1,5)
% plot(k,e.^-((1/sigmak)^2*((k-q).^2/2-(xnot)*sigmak^2*i*k)));


%  t=100;
%         x=0;       
% subplot(5,1,2)
% plot(k,-1);
% 
% subplot(5,1,3)
% plot(k,-1+2*i*coth(knot*d)/knot*k)
% subplot(5,1,4)
% plot(k,-1+2*i*coth(knot*d)/knot*k+4*k.^2.*(coth(knot*d)).^2/knot.^2)
% subplot(5,1,5)

%{
%DATA AQUISITION Numerical
q=-1;
B=-100.0;
C= 100.0;
deltax=0.5;
ARR=zeros(20,10,1);                   % (x,t,p)  
z = zeros(20,1);
v=1;
D=100;
while t<=D
x=B;
l=1;
while x<=C
    
if (x<0)  


             ya1=((A.*e.^(i*k*(x))+w*A.*e.^(-i*k*(x)))).*e.^(i*k.^2*hbar*t/(2*m));
            Gauss=e.^-((1/sigmak)^2*((k-q).^2/2+(xnot)*sigmak^2*i*k));
            y1= ya1.*Gauss;
             ARR(l,t+1,1)=((deltax)*(abs(trapz(k,y1))).^2);
            z(l,1) = B+(l-1)*deltax;
end

 if all([x>=0,x<=d])
  
            ya2=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             yb2=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             Gauss=e.^-((1/sigmak)^2*((k-q).^2/2+(xnot)*sigmak^2*i*k));
             y2=(ya2+yb2).*Gauss;
             ARR(l,t+1,1)=((deltax)*(abs(trapz(k,y2))).^2);    
                  z(l,1) = B+(l-1)*deltax;
 end
 
 if (x>d) 
   
             ya3=(A.*(kt-i*k)+A.*w.*(kt+i*k))./(2*kt).*e.^(-kt*d-i*k*d+i*k*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             yb3=(A.*(kt+i*k)+A.*w.*(kt-i*k))./(2*kt).*e.^(kt*d-i*k*d+i*k*(x)).*e.^(i*k.^2*hbar*t/(2*m));
             Gauss=e.^-((1/sigmak)^2*((k-q).^2/2+(xnot)*sigmak^2*i*k));
             y3=(ya3+yb3).*Gauss;
             ARR(l,t+1,1)=((deltax)*(abs(trapz(k,y3))).^2); 
             z(l,1) = B+(l-1)*deltax;
    
 end

 x=x+deltax;
 l=l+1;
end
fprintf('%d\n',t+1); %
 t=t+1;       
end       

 p=0;
while p<D/4
  subplot(5,5,p+1)
        plot(z(:,1),ARR(:,4*p+1)); 
        
        p=p+1;
        
end
%} 

%{
%DATA AQUISITION Analytical

B=-100.0;
C= 100.0;
deltax=0.2;
ARR=zeros(200,100,1);                   % (x,t,p)  
z = zeros(200,1);
D=100;
ktq=(knot.^2-q.^2).^(.5);
wq=(-(knot)^2)./(-2*i*q.*ktq.*coth(ktq*d)-q.^2+ktq.^2);
dwq=(2*i*knot^2*((knot^2-2*q^2)*coth(ktq*d)+q*ktq*(d*q*(csc(ktq*d))^2-2*i)))/(ktq*(2*q*ktq*coth(ktq*d)+i*(knot^2-2*q^2))^2);

while t<=D
x=B;
l=1;
while x<=C
    
if (x<0)  

            Xt=(A*(m)^.5*e.^(i*q*(x-xnot)+(-i*q^2*hbar*t/(2*m))+(-(x-xnot-q*hbar*t/m)^2*m*sigmak^2/(2*(i*hbar*sigmak^2*t+m)))))/(i*hbar*sigmak^2*t+m)^.5;
            Xr=(A*(m)^.5*e.^(-i*q*(x+xnot)+(-i*q^2*hbar*t/(2*m))+(-(x-xnot-q*hbar*t/m)^2*m*sigmak^2/(2*(i*hbar*sigmak^2*t+m))))*(wq-dwq*(i*(x+xnot+q*hbar*t/m))*sigmak^2*m/(i*hbar*t*sigmak^2+m)))/(i*hbar*sigmak^2*t+m)^.5;
            X=Xt+Xr;
            ARR(l,t+1,1)=((deltax)*(abs(X))^2);
            z(l,1) = B+(l-1)*deltax;
            end
% end
% AA=(1+wq-dwq*q);
% BB=i*(-1+wq-dwq*q);
% CC=(q/(knot^2-q^2)^.5-q*knot^2/(knot^2-q^2)^(1.5));
% DD=q^2/(knot^2-q^2)^.5-(q^4+2q^2*(knot^2-q^2))/(knot^2-q^2)^(1.5);
% EE=i*dwq*(q^3+2q*(knot^2-q^2))/(knot^2-q^2)^(1.5)
% if all[x<d,x>0]
% 
% end
%  
% AA=(1+wq-dwq*q);
% BB=i*(-1+wq-dwq*q);
% CC=(q/(knot^2-q^2)^.5-q*knot^2/(knot^2-q^2)^(1.5));
% DD=q^2/(knot^2-q^2)^.5-(q^4+2q^2*(knot^2-q^2))/(knot^2-q^2)^(1.5);
% EE=i*dwq*(q^3+2q*(knot^2-q^2))/(knot^2-q^2)^(1.5) 
% G1=
% G2=
% P1=
% P2=
% F1=
% F2=
% 
%  if (x>d) 
%   
%  end

 x=x+deltax;
 l=l+1;
end
fprintf('%d\n',t+1); %
 t=t+1;  
 x=B;
l=1;
end       

 p=0;
while p<D/4
  subplot(5,5,p+1)
        plot(z(:,1),ARR(:,4*p+1)); 
        
        p=p+1;
        
end
%}        
  

ktq=(knot.^2-q.^2).^(.5);
wq=(-(knot)^2)./(-2*i*q.*ktq.*coth(ktq*d)-2*q.^2+ktq.^2);
dwq=(2*i*knot^2*((knot^2-2*q^2)*coth(ktq*d)+q*ktq*(d*q*(csc(ktq*d))^2-2*i)))/(ktq*(2*q*ktq*coth(ktq*d)+i*(knot^2-2*q^2))^2);

k=[(-3*sigmak+q):deltak:(3*sigmak+q)];  
q=1;
subplot(2,2,1)
plot(k,w);
subplot(2,2,2)
plot(k,wq+dwq*(k-q));        
subplot(2,2,3)
plot(k,i*w);
subplot(2,2,4)
 plot(k,i*(wq+dwq*(k-q)));       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
 %{
 
 while t<100
  
 t
while x < 10

k=[-3*sigma:deltak/10000:3*sigma];

w1= -(knot)^2 * 2*sinh( (knot.^2-k.^2).^(0.5)*d );

w2=(-4*i*k.*(knot.^2-k.^2).^(.5).*cosh((knot.^2-k.^2).^(.5)*d)-2*k.^2.*sinh((knot.^2-k.^2).^(.5)*d)+2*(knot.^2-k.^2).*sinh((knot.^2-k.^2).^(.5)*d)); 

w=w1./w2;

y1=(A*((knot.^2-k.^2).^(0.5)-i*k)+A*w.*((knot.^2-k.^2).^(0.5)+i*k))./(2*(knot.^2-k.^2).^(0.5)).*e.^(-(knot.^2-k.^2).^(0.5)*d-i*k*d+i*k*x);

y2=(A*((knot.^2+k.^2).^(0.5)+i*k)+A*w.*((knot.^2-k.^2).^(0.5)-i*k))./(2*(knot.^2-k.^2).^(0.5)).*e.^((knot.^2-k.^2).^(0.5)*d-i*k*d+i*k*x);

y3=1/(sigma*(2*pi).^(0.5)).*e.^-((k-q).^2/(2*(sigma).^2)+(i*k.^2*hbar*t)/(2*m));

y=(y1+y2).*y3;

SUMofprob = SUMofprob + abs(trapz(k,y)).^2;

x= x+1;
end

x=1;


SUMofprob

t=t+1;


 end


kb=[-3*sigma:deltak/1000:3*sigma];

1b= -(knot)^2 * 2*sinh( (knot.^2-kb.^2).^(0.5)*d );

2b=(-4*i*kb.*(knot.^2-kb.^2).^(.5).*cosh((knot.^2-kb.^2).^(.5)*d)-2*kb.^2.*sinh((knot.^2-kb.^2).^(.5)*d)+2*(knot.^2-kb.^2).*sinh((knot.^2-kb.^2).^(.5)*d)); 

wb=w1b./w2b;

y1b=((A*(knot.^2-kb.^2).^(0.5)-i*kb)+A*wb.*((knot.^2-kb.^2).^(0.5)+i*kb))./(2*(knot.^2-kb.^2).^(0.5)).*e.^(-(knot.^2-kb.^2).^(0.5)*d-i*kb*d+i*kb*x);

y2b=((A*(knot.^2+kb.^2).^(0.5)-i*kb)+A*wb.*((knot.^2-kb.^2).^(0.5)-i*kb))./(2*(knot.^2-kb.^2).^(0.5)).*e.^((knot.^2-kb.^2).^(0.5)*d-i*kb*d+i*kb*x);

y3b=1/(sigma*(2*pi).^(0.5)).*e.^-((kb-q).^2/(2*(sigma).^2)+(i*kb.^2*hbar*t)/(2*m));

yb=y1b+(y2b.*y3b);

trapz(kb,yb)




kc=[-3*sigma:deltak/1000000:3*sigma];

w1c= -(knot)^2 * 2*sinh( (knot.^2-kc.^2).^(0.5)*d );

w2c=(-4*i*kc.*(knot.^2-kc.^2).^(.5).*cosh((knot.^2-kc.^2).^(.5)*d)-2*kc.^2.*sinh((knot.^2-kc.^2).^(.5)*d)+2*(knot.^2-kc.^2).*sinh((knot.^2-kc.^2).^(.5)*d)); 

wc=w1c./w2c;

y1c=((A*(knot.^2-kc.^2).^(0.5)-i*kc)+A*wc.*((knot.^2-kc.^2).^(0.5)+i*kc))./(2*(knot.^2-kc.^2).^(0.5)).*e.^(-(knot.^2-kc.^2).^(0.5)*d-i*kc*d+i*kc*x);

y2c=((A*(knot.^2+kc.^2).^(0.5)-i*kc)+A*wc.*((knot.^2-kc.^2).^(0.5)-i*kc))./(2*(knot.^2-kc.^2).^(0.5)).*e.^((knot.^2-kc.^2).^(0.5)*d-i*kc*d+i*kc*x);

y3c=1/(sigma*(2*pi).^(0.5)).*e.^-((kc-q).^2/(2*(sigma).^2)+(i*kc.^2*hbar*t)/(2*m));

yc=y1c+(y2c.*y3c);

trapz(kc,yc)



kd=[-3*sigma:deltak/10000000:3*sigma];

w1d= -(knot)^2 * 2*sinh( (knot.^2-kd.^2).^(0.5)*d );

w2d=(-4*i*kd.*(knot.^2-kd.^2).^(.5).*cosh((knot.^2-kd.^2).^(.5)*d)-2*kd.^2.*sinh((knot.^2-kd.^2).^(.5)*d)+2*(knot.^2-kd.^2).*sinh((knot.^2-kd.^2).^(.5)*d)); 

wd=w1d./w2d;

y1d=((A*(knot.^2-kd.^2).^(0.5)-i*kd)+A*wd.*((knot.^2-kd.^2).^(0.5)+i*kd))./(2*(knot.^2-kd.^2).^(0.5)).*e.^(-(knot.^2-kd.^2).^(0.5)*d-i*kd*d+i*kd*x);

y2d=((A*(knot.^2+kd.^2).^(0.5)-i*kd)+A*wd.*((knot.^2-kd.^2).^(0.5)-i*kd))./(2*(knot.^2-kd.^2).^(0.5)).*e.^((knot.^2-kd.^2).^(0.5)*d-i*kd*d+i*kd*x);

y3d=1/(sigma*(2*pi).^(0.5)).*e.^-((kd-q).^2/(2*(sigma).^2)+(i*kd.^2*hbar*t)/(2*m));

yd=y1d+(y2d.*y3d);

trapz(kd,yd) 
%} 


