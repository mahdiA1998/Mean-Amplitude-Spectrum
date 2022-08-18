clc
clear
close all
load('S11_Session1.mat');

%% %% Mean Amplitude Spectrums for each Stimulation Frequency_ear-EEG
t_class1 = zeros(1, 50);
t_class2 = zeros(1, 50);
t_class3 = zeros(1, 50);
c = 1;
d = 1;
e = 1;
for i = 1:150
    if ear_cnt.y_dec(i) == 1
       t_class1(c) = ear_cnt.t(i);
       c = c + 1;
    elseif ear_cnt.y_dec(i) == 2
       t_class2(d) = ear_cnt.t(i);
       d = d + 1;
    else
       t_class3(e) = ear_cnt.t(i);
       e = e + 1;
    end
end

class1 = zeros(18, 5000);
class2 = zeros(18, 5000);
class3 = zeros(18, 5000);
for channs = 1:18
    for i=1:50
        class1(channs,:)=class1(channs,:)+ear_cnt.x(channs, t_class1(i)-1000:t_class1(i)+3999);
        class2(channs,:)=class2(channs,:)+ear_cnt.x(channs, t_class2(i)-1000:t_class2(i)+3999);
        class3(channs,:)=class3(channs,:)+ear_cnt.x(channs, t_class3(i)-1000:t_class3(i)+3999);
    end
end
class1_ave = class1/50;
class2_ave = class2/50;
class3_ave = class3/50;
% subplot(2,1,1);
% plot(ear_cnt.x(1, t_class1(1)-1000:t_class1(1)+3999));
% subplot(2,1,2);
% plot(class1_ave(1,:));

channsAve1 = (class1_ave(8,:)+class1_ave(18,:))/2;      % Average for 2 channels(#8 and 18)
channsAve2 = (class2_ave(8,:)+class2_ave(18,:))/2;
channsAve3 = (class3_ave(8,:)+class3_ave(18,:))/2;

for channs = 1:18
     ssvep1 = downsample(channsAve1, 2);
     ssvep2 = downsample(channsAve2, 2);
     ssvep3 = downsample(channsAve3, 2);
end

fs = 250;                     %sampling frequency
[b,a] = butter(5,[6 50]/(fs/2),'bandpass');             % filter (6-50 Hz)
ssvep1 = filtfilt(b,a,ssvep1');
ssvep1 = ssvep1';
ssvep2 = filtfilt(b,a,ssvep2');
ssvep2 = ssvep2';
ssvep3 = filtfilt(b,a,ssvep3');
ssvep3 = ssvep3';

Amp_spec1 = fft(ssvep1);
Amp_spec2 = fft(ssvep2);
Amp_spec3 = fft(ssvep3);

p21 = abs(Amp_spec1/2500);
p11 = p21(1:(2500/2)+1);
p11(2:end-1) = 2*p11(2:end-1);
f = 250*(0:(2500/2))/2500;
subplot(1,3,1);
plot(f, p11);
hold on
title('Amplitude Spectrum for 10Hz flicker');
xlabel('Frequency(Hz)');
ylabel('Amplitude');

p22 = abs(Amp_spec2/2500);
p12 = p22(1:(2500/2)+1);
p12(2:end-1) = 2*p12(2:end-1);
subplot(1,3,2);
plot(f, p12);
hold on
title('Amplitude Spectrum for 8.57Hz flicker');
xlabel('Frequency(Hz)');
ylabel('Amplitude');

p23 = abs(Amp_spec3/2500);
p13 = p23(1:(2500/2)+1);
p13(2:end-1) = 2*p13(2:end-1);
subplot(1,3,3);
plot(f, p13);
hold on
title('Amplitude Spectrum for 7.5Hz flicker');
xlabel('Frequency(Hz)');
ylabel('Amplitude');

%% Mean Amplitude Spectrums for each Stimulation Frequency_scalp-EEG
t_class1 = zeros(1, 50);
t_class2 = zeros(1, 50);
t_class3 = zeros(1, 50);
c = 1;
d = 1;
e = 1;
for i = 1:150
    if cap_cnt.y_dec(i) == 1
       t_class1(c) = cap_cnt.t(i);
       c = c + 1;
    elseif cap_cnt.y_dec(i) == 2
       t_class2(d) = cap_cnt.t(i);
       d = d + 1;
    else
       t_class3(e) = cap_cnt.t(i);
       e = e + 1;
    end
end

class1 = zeros(8, 5000);
class2 = zeros(8, 5000);
class3 = zeros(8, 5000);
for channs = 1:8
    for i=1:50
        class1(channs,:)=class1(channs,:)+cap_cnt.x(channs, t_class1(i)-1000:t_class1(i)+3999);
        class2(channs,:)=class2(channs,:)+cap_cnt.x(channs, t_class2(i)-1000:t_class2(i)+3999);
        class3(channs,:)=class3(channs,:)+cap_cnt.x(channs, t_class3(i)-1000:t_class3(i)+3999);
    end
end
class1_ave = class1/50;
class2_ave = class2/50;
class3_ave = class3/50;

channsAve1 = class1_ave(7,:);
channsAve2 = class2_ave(7,:);
channsAve3 = class3_ave(7,:);

for channs = 1:18
     ssvep1 = downsample(channsAve1, 2);
     ssvep2 = downsample(channsAve2, 2);
     ssvep3 = downsample(channsAve3, 2);
end

fs = 250;                     %sampling frequency
[b,a] = butter(5,[6 50]/(fs/2),'bandpass');
ssvep1 = filtfilt(b,a,ssvep1');
ssvep1 = ssvep1';
ssvep2 = filtfilt(b,a,ssvep2');
ssvep2 = ssvep2';
ssvep3 = filtfilt(b,a,ssvep3');
ssvep3 = ssvep3';

Amp_spec1 = fft(ssvep1);
Amp_spec2 = fft(ssvep2);
Amp_spec3 = fft(ssvep3);

p21 = abs(Amp_spec1/2500);
p11 = p21(1:(2500/2)+1);
p11(2:end-1) = 2*p11(2:end-1);
f = 250*(0:(2500/2))/2500;
subplot(1,3,1);
plot(f, p11);
hold on
title('Amplitude Spectrum for 10Hz flicker');


p22 = abs(Amp_spec2/2500);
p12 = p22(1:(2500/2)+1);
p12(2:end-1) = 2*p12(2:end-1);
subplot(1,3,2);
plot(f, p12);
hold on
title('Amplitude Spectrum for 8.57Hz flicker');

p23 = abs(Amp_spec3/2500);
p13 = p23(1:(2500/2)+1);
p13(2:end-1) = 2*p13(2:end-1);
subplot(1,3,3);
plot(f, p13);
hold on
title('Amplitude Spectrum for 7.5Hz flicker');
