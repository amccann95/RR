function plot_rp_ecg(ecg_WS,REC,DET,RP,database)

fs = ecg_WS.Samplingfrequency;
lengthECGToPlot = 60; % seconds

%samps = [1:lengthECGToPlot*fs];
%ecgData = ecg_WS.data(samps);
%[R_inds,RR] = detectRR(ecg_WS,samps); 

time = [1/fs:1/fs:lengthECGToPlot];
% figure;
% subplot(2,1,1); plot(time,ecgData(samps)); hold on; plot(time(R_inds),ecgData(R_inds),'*');
% xlimit('time (s)'); limit('ECG (mV)');

% subplot(2,1,1); plot(1/fs:1/fs:length(ecg_WS.data)/fs,ecg_WS.data); 
% hold on; plot(ecg_WS.R_inds,ecg_WS.data(ecg_WS.R_inds));
% subplot(2,1,2); plot(1:length(ecg_WS.RR),ecg_WS.RR); 
% xlabel('IBI interval number'); ylabel('IBI interval duration (ms)');


figure; 
plot(ecg_WS.RR); 
xlabel('RR interval number');
ylabel('RR interval duration (ms)');
xlim([0, length(ecg_WS.RR)]);
ylim([250 1750]);
set(gca, 'fontsize', 12);
set(gca,'linewidth', 1.2);
set(gca,'FontName','carlito');
%set(gca,'box','off');

figure;
imagesc(RP==0);
colormap('gray');
axis('square');
set(gca,'Ydir','normal');
title(strcat(sprintf('REC: %0.2f%%',100*REC),sprintf(' - DET: %0.2f%%',100*DET)));
set(gca, 'fontsize', 12);
set(gca,'FontName','carlito');
set(gca,'linewidth', 1.2);
%set(gca,'box','off');


end