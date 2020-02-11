% MGA852
% LAB 2
% Compteur de pas

% read csv file and create 
% time, x, y, z matrices
clear all
close all
[num, txt, raw] = xlsread('LAB2_Part1_testdata1.csv');
time_matrix     = str2double(txt(2:end, 1));
x_matrix        = str2double(txt(2:end, 2));
y_matrix        = str2double(txt(2:end, 3));
z_matrix        = str2double(txt(2:end, 4));
% plot
% X
figure
plot(time_matrix, x_matrix);
title('X vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on
% Y
figure
plot(time_matrix, y_matrix);
title('Y vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on
% Z
z_matrix_no_g = z_matrix - nanmean(z_matrix);
figure
plot(time_matrix, z_matrix_no_g);
title('Z vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on
% Acceleration Magnitude
acc_mag = sqrt(sum(x_matrix.^2 + y_matrix.^2 + z_matrix_no_g.^2,2));
figure
plot(time_matrix, acc_mag);
title('Total Acceleration vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on

minPeakHeight = nanstd(acc_mag);
[pks,locs] = findpeaks(acc_mag,'MINPEAKHEIGHT',minPeakHeight);
numSteps = numel(pks);
% hold on
% plot(time_matrix(locs), pks, 'r', 'Marker', 'v', 'LineStyle', 'none');
% title('Counting Steps');
% xlabel('Time (s)');
% ylabel('Acceleration Magnitude, No Gravity (m/s^2)');
% hold off


