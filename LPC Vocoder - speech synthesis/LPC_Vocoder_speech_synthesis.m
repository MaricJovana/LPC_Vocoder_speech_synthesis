clear all; close all; clc; 

%KREIRANJE POBUDNOG SIGNALA
f0_matrica=load('f0_recenica 33.mat');
f0=f0_matrica.f0;
[x,fs]=audioread("recenica33.wav");
tau=0.03;
preklapanje=0.5;
N=length(f0);
T=ceil(fs*tau);
%zvucnost  zvucni su jedinice, bezvucni 0
zvucnost= ones(1,N);
for br = 1:N
    if(isnan(f0(br)))
      zvucnost(br)= 0;
    end
end
%izdvajamo odbirke
k=1;
for br = 0:T/2:length(x)-T
    x_odb(k,:)=x(br+1:br+T);
    k=k+1;
end
R=T-T*preklapanje;

%pobudni signal
pobudni_signal=zeros(N,R);
%prvi prozor
s=0;
 if(zvucnost(1)==0) %prvi odbirak je bezvucan
      pobudni_signal(1,:)=0.01*randn(1,R);
 else %prvi odbirak je zvucan
     s=s+1;
     pobudni_signal(br,(1:ceil(f0(1)):R))=1;
     lpamti(s)=ceil(mod(R,f0(br)));
 end
%ostali prozori
for br = 2:N
    if(zvucnost(br)==0) %sum
      pobudni_signal(br,:)=0.01*randn(1,R);
    elseif((zvucnost(br)==1) && (zvucnost(br-1)==1))%zvucni pa zvucni
        s=s+1;
        loc=ceil(f0(br)-lpamti(s-1));
        if(loc>0)% t2 je veca od razlike
            pobudni_signal(br,(loc:ceil(f0(br)):R))=1;
            lpamti(s)=ceil(mod(R-loc,f0(br)));
        elseif(loc<0)
            pobudni_signal(br-1,(R+loc))=1;
            pobudni_signal(br,(ceil(f0(br))+loc:ceil(f0(br)):R))=1; 
            lpamti(s)=ceil(mod(R-(ceil(f0(br))+loc),f0(br)));
        else 
            pobudni_signal(br,(1:ceil(f0(br)):R))=1;
            lpamti(s)=ceil(mod(R,f0(br)));
        end
    else %bezvucni pa zvucni
        s=s+1;
        pobudni_signal(br,(1:ceil(f0(br)):R))=1;
        lpamti(s)=ceil(mod(R,f0(br)));
    end
end
pobuda=pobudni_signal';
pobuda=pobuda(:)';
figure,stem((0:length(pobuda)-1)/fs,pobuda); title('Pobudni signal');xlabel('t[s]');

%SINTEZA GOVORA POMOCU LPC ANALIZE
%mnozenje sa hamingovom prozorskom funkcijom
x_odb_niz=x_odb';
x_odb_niz=x_odb_niz(:);
x_h=zeros(N,T);
for br = 0:N-1
    x_w=x_odb_niz(br*T+1:T+br*T).*hamming(T);
    x_h(br+1,:)=x_w;
end
%autolpc
 p=35;
 for br= 1:N
    [A(:,br),G(br)]=autolpc(x_h(br,:),p);
    Gain(br)=G(br)/(sqrt(sum(pobudni_signal(br,:).^2))+0.01); 
    x_filter(br,:)= filter(Gain(br),A(:,br),pobudni_signal(br,:));
 end

 %prikaz
 x_i=x_filter';
 x_i=x_i(:)'; 
 figure, plotyy((0:length(x)-1)/fs,x/abs(max(x)),(0:length(x_i)-1)/fs,x_i/abs(max(x_i))); title('Vremenski oblik izlaznog i originalnog signala');
 legend('originalni','sintetizovan'); xlabel('t[s]');
 sound(x_i,fs);