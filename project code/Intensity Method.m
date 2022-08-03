clear all;
close all;
clc ;

load('Sub1.mat')
ppg = Sessiondata(4:end,4);
ppg = cell2mat(ppg);
ppg = ppg/abs(max(ppg));
ppg = smooth(ppg);
figure;
findpeaks(ppg,'MinPeakHeight',0)
[~,locs] = findpeaks(ppg,'MinPeakHeight',0);
x = 0;
g = 0;
 for i = 1:114
     a = locs(i+1);
     b = locs(i);
     d = a-b;
     x = x + d;
     g = g+1;
 end
result = x/g;
Heart_Rate_NeXus_PPG = 128 * 60 / result 

load('Sub1_smrt.mat')
s_ppg = Sessiondata(77:end,2);
s_ppg = cell2mat(s_ppg);
s_ppg = s_ppg/abs(max(s_ppg));
s_ppg = smooth(s_ppg);
figure;
findpeaks(s_ppg)
[~,locs] = findpeaks(s_ppg);
x = 0;
g = 0;
 for i = 1:116
     a = locs(i+1);
     b = locs(i);
     d = a-b;
     x = x + d;
     g = g+1;
 end
 result1 = x/g;
Heart_Rate_SmartPhone_Wrist = 29 * 60 / result1

load('Sub1_smrt.mat')
s_ppg = Sessiondata(77:end,4);
s_ppg = cell2mat(s_ppg
s_ppg = s_ppg/abs(max(s_ppg));
s_ppg = smooth(s_ppg);
figure;
findpeaks(s_ppg)
[pks,locs] = findpeaks(s_ppg);
x = 0;
g = 0;
 for i = 1:92
     a = locs(i+1);
     b = locs(i);
     d = a-b;
     x = x + d;
     g = g+1;
 end
 result1 = x/g;
Heart_Rate_SmartPhones_Fingertip = 15 * 60 / result1