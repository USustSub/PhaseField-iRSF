L = 150000 ; 
% eval(['!sed -i ',char(39),'3s/.*/L = ',num2str(L,10),...
%      '; /',char(39),' in/3SolidFib.geo']);
% for ii__ = [11 31 61]
figure
ii = 0 ; 
for l_ = [L/500 L/1000 L/3000 L/6000]
    ii = ii + 1 ; 
    subplot(2,2,ii)
hold on

eval(['load out_results/out' num2str(l_) '_' num2str(0.5) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'r-','LineWidth',1)
eval(['load out_results/out' num2str(l_) '_' num2str(1) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'b-','LineWidth',1)
eval(['load out_results/out' num2str(l_) '_' num2str(2) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'g-','LineWidth',1)
eval(['load out_results/out' num2str(l_) '_' num2str(3) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'c-','LineWidth',1)
eval(['load out_results/out' num2str(l_) '_' num2str(4) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'m-','LineWidth',1)
load out_ref/slip/jj600.mat
plot(ff(:,1),ff(:,2),'k-','LineWidth',2)


pause
end

%%
figure
ii = 0 ; 
for l_ = [L/500 L/1000 L/3000 L/6000]
    ii = ii + 1 ; 
    subplot(2,2,ii)
ff = dlmread(['out_results/out' num2str(l_) '_' num2str(0.5) '/LD.out']);
semilogy(ff(:,1),ff(:,2),'r-','LineWidth',2)
hold on
ff = dlmread(['out_results/out' num2str(l_) '_' num2str(1) '/LD.out']);
semilogy(ff(:,1),ff(:,2),'b-','LineWidth',2)
ff = dlmread(['out_results/out' num2str(l_) '_' num2str(2) '/LD.out']);
semilogy(ff(:,1),ff(:,2),'g-','LineWidth',2)
ff = dlmread(['out_results/out' num2str(l_) '_' num2str(3) '/LD.out']);
semilogy(ff(:,1),ff(:,2),'c-','LineWidth',2)
ff = dlmread(['out_ref/LD.out']);
semilogy(ff(:,1),ff(:,2),'k-','LineWidth',2)

pause
end


%%
cd ../Interface_phasefield_2_2_2
L = 150000 ; 
% eval(['!sed -i ',char(39),'3s/.*/L = ',num2str(L,10),...
%      '; /',char(39),' in/3SolidFib.geo']);
% for ii__ = [11 31 61]
figure
ii = 0 ; 
for l_ = [L/1000 L/3000 L/6000]
    ii = ii + 1 ; 
    subplot(1,3,ii)
hold on

eval(['load out_results/out' num2str(l_) '_' num2str(0.5) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'r-','LineWidth',1)
eval(['load out_results/out' num2str(l_) '_' num2str(1) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'b-','LineWidth',1)
eval(['load out_results/out' num2str(l_) '_' num2str(2) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'g-','LineWidth',1)
eval(['load out_results/out' num2str(l_) '_' num2str(3) '/slip/jj600.mat'])
plot(ff(:,1),ff(:,2),'c-','LineWidth',1)
% eval(['load out_results/out' num2str(l_) '_' num2str(4) '/slip/jj600.mat'])
% plot(ff(:,1),ff(:,2),'m-','LineWidth',1)
load ../Interface_phasefield_2_2/out_ref/slip/jj600.mat
plot(ff(:,1),ff(:,2),'k-','LineWidth',2)


% pause
end
