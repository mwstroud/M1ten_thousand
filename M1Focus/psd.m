function [] = psd ()
% PSD plot
file = './output/ecp.h5';
data = h5read(file,'//ecp/data');
fs=10000;
[pxx,f]=periodogram(data,[],length(data),fs);
semilogx(f,10*log10(pxx)); xlim([0 200]);
title("M1 Cortex PSD");
xlabel('Frequency(Hz)');
ylabel('PSD (dB/Hz)');
saveas(gcf, './results/M1_PSD.png')
end