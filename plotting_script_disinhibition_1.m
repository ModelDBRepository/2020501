clear all
close all
clc
load disinhibition_model_1.mat
%%
f1 = figure(1)

my = 4
mx = 4; 
dt_long = T/kd;
sf = 0.1/dt_long;

dat = REC2(1:sf:kd,1:NE);
dt = 0.1;
time = (1:size(dat,1))*T/(size(dat,1));

tcrit1 = 60; 
tcrit2 = 240;
BD1 = dat((time>tcrit1).*(time<tcrit2)==1,:);
tshort1 = dt*(1:size(BD1,1)) + tcrit1;

tcrit3 = 360;
tcrit4 = 540;
BD2 = dat((time>tcrit3).*(time<tcrit4)==1,:);
tshort2 = dt*(1:size(BD2,1)) + tcrit3;

tcrit5 = 660;  
tcrit6 = 840;
BD3 = dat((time>tcrit5).*(time<tcrit6)==1,:);
tshort3 = dt*(1:size(BD3,1)) + tcrit5;





subplot(my,mx,[9,13])
plot((1:kd)*T/kd,mean(REC2(:,1:NE)'),'LineWidth',2)
ylim([min(mean(REC2(:,1:NE)'))*0.9,max(mean(REC2(:,1:NE)'))*1.1])
%ylim([mean(mean(REC2'))*0.5,max(mean(REC2'))*1.1])
xlabel('Time (s)')
ylabel('Mean of Calcium Traces')
xlim([tcrit1,tcrit2])

subplot(my,mx,[10,14])
plot((1:kd)*T/kd,mean(REC2(:,1:NE)'),'LineWidth',2)
ylim([min(mean(REC2(:,1:NE)'))*0.9,max(mean(REC2(:,1:NE)'))*1.1])
%ylim([mean(mean(REC2'))*0.5,max(mean(REC2'))*1.1])
xlabel('Time (s)')
ylabel('Mean of Calcium Traces')


xlim([tcrit3,tcrit4])
subplot(my,mx,[11,15])
plot((1:kd)*T/kd,mean(REC2(:,1:NE)'),'LineWidth',2)
ylim([min(mean(REC2(:,1:NE)'))*0.9,max(mean(REC2(:,1:NE)'))*1.1])
%ylim([mean(mean(REC2'))*0.5,max(mean(REC2'))*1.1])
xlabel('Time (s)')
ylabel('Mean of Calcium Traces')
xlim([tcrit5,tcrit6])








%%
i1del = find(mean(BD1)==0)';
i2del = find(mean(BD2)==0)';
i3del = find(mean(BD3)==0)';
idel = [i1del;i2del;i3del]; 
BD1(:,idel)=[];
BD2(:,idel)=[];
BD3(:,idel)=[];

BD1 = BD1(:,1:40);
BD2 = BD2(:,1:40);
BD3 = BD3(:,1:40);
N = size(BD1,2);
%%

for i = 1:N 
NBD1(:,i) = BD1(:,i)/sqrt( trapz(tshort1,BD1(:,i).^2));   
end
for i = 1:N 
NBD2(:,i) = BD2(:,i)/sqrt( trapz(tshort2,BD2(:,i).^2));   
end
for i = 1:N 
NBD3(:,i) = BD3(:,i)/sqrt( trapz(tshort3,BD3(:,i).^2));   
end


[u1,s1,v1] = svd(BD1);
[u2,s2,v2] = svd(BD2);
[u3,s3,v3] = svd(BD3);

[m11,i11] = max(abs(u1(:,1)));
[m12,i12] = max(abs(u1(:,2)));
s11 = sign(u1(i11,1));
s12 = sign(u1(i12,2));

[m21,i21] = max(abs(u2(:,1)));
[m22,i22] = max(abs(u2(:,2)));
s21 = sign(u2(i21,1));
s22 = sign(u2(i22,2));

[m31,i31] = max(abs(u3(:,1)));
[m32,i32] = max(abs(u3(:,2)));
s31 = sign(u3(i31,1));
s32 = sign(u3(i32,2));

coeff1 = s12*v1(:,2);
coeff2 = s22*v2(:,2);
coeff3 = s32*v3(:,2);
[c1,ic1] = sort(coeff1);
[c2,ic2] = sort(coeff2); 
[c3,ic3] = sort(coeff3); 





figure(1)
subplot(my,mx,[1,5])
imagesc(tcrit1,(1:(N/10)),BD1(:,ic1)',[0,30])
xlabel('Time (s)')
ylabel('Neuron Index (Calcium Traces)')
set(gca,'Ydir','Normal')
title('Home Cage')
colormap('hot')

subplot(my,mx,[2,6])
imagesc(tcrit2,(1:(N/10)),BD2(:,ic1)',[0,30])
xlabel('Time (s)')
ylabel('Neuron Index (Calcium Traces)')
set(gca,'Ydir','Normal')
title('Up-State')
colormap('hot')

subplot(my,mx,[3,7])
imagesc(tcrit3,(1:(N/10)),BD3(:,ic1)',[0,30])
xlabel('Time (s)')
ylabel('Neuron Index (Calcium Traces)')
set(gca,'Ydir','Normal')
title('Return to HC')
colormap('hot')







subplot(my,mx,[4])
plot(tshort1,s11*u1(:,1),'r','LineWidth',2), hold on 
plot(tshort1,s12*u1(:,2),'b','LineWidth',2)
ylim([-0.1,0.3])
legend('Feature 1','Feature 2')
xlim([tshort1(1),tshort1(end)])
xlabel('Time (s)')
ylabel('Features')
subplot(my,mx,[8])
plot(tshort2,s21*u2(:,1),'r','LineWidth',2), hold on 
plot(tshort2,s22*u2(:,2),'b','LineWidth',2)
ylim([-0.1,0.3])
xlabel('Time (s)')
ylabel('Features')
xlim([tshort2(1),tshort2(end)])
subplot(my,mx,[12])
plot(tshort3,s31*u3(:,1),'r','LineWidth',2), hold on 
plot(tshort3,s32*u3(:,2),'b','LineWidth',2)
ylim([-0.1,0.3])
xlabel('Time (s)')
ylabel('Features')
xlim([tshort3(1),tshort3(end)])
%% 


subplot(my,mx,[16]);

[u1,s1,v1] = svd(BD1);
[u2,s2,v2] = svd(BD2);
[u3,s3,v3] = svd(BD3); 
clear error 

for k = 0:N 
    if k > 0 
s1p=0*s1; s1p(1:k,1:k) = s1(1:k,1:k);
appx = u1*s1p*v1';
    else 
        appx = 0;
    end
error(k+1) = norm(appx-BD1,'fro');
end
plot(0:N,100*error/norm(-BD1,'fro'),'b','LineWidth',2), hold on 
[m1hc,i1hc] = min(abs(error/norm(-BD1,'fro')- 0.1));

clear error 
for k = 0:N
    if k > 0 
s2p=0*s2; s2p(1:k,1:k) = s2(1:k,1:k);
appx = u2*s2p*v2';
    else 
        appx = 0;
    end
error(k+1) = norm(appx-BD2,'fro');
end
plot(0:N,100*error/norm(-BD2,'fro'),'k','LineWidth',2)
[m1us,i1us] = min(abs(error/norm(-BD2,'fro') - 0.1));

clear error
for k = 0:N
    if k >0 
s3p=0*s3; s3p(1:k,1:k) = s3(1:k,1:k);
appx = u3*s3p*v3';
    else 
        appx = 0;
    end
error(k+1) = norm(appx-BD3,'fro');
end
plot(0:N,100*error/norm(-BD3,'fro'),'r','LineWidth',2)
xlabel('Number of Features')
ylabel('% Error in Reconstruction')

 

 
%set(gca,'FontSize',FS)
[m1rhc,i1rhc] = min(abs(error/norm(-BD3,'fro') - 0.1));
title(sprintf('HC - %d, US - %d, HC - %d',i1hc,i1us,i1rhc))
%dim_rel(dan,:) = [i1hc,i1us,i1rhc];
%% 
set(f1,'paperposition',[0,0,15,10])
print(f1,'disinhibition_1.png','-dpng','-r300')
print(f1,'disinhibition_1.svg','-dsvg','-painters')
