% MGA852
% LAB 2
% Compteur de pas

% read csv file and create 
% time, x, y, z matrices
clear all
close all
clc
[num, txt, raw] = xlsread('LAB2_Part1_testdata1_mod.csv');
time_matrix     = num(1:end, 1);
x_matrix        = num(1:end, 2);
y_matrix        = num(1:end, 3);
z_matrix        = num(1:end, 4);
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
z_matrix_no_g = z_matrix - nanmean(z_matrix); % substract gravity from acc

indx = 1;
for a=0.0:0.1:4
    minPeakHeight(indx) = a * nanstd(z_matrix_no_g);
    [pks,locs] = findpeaks(z_matrix_no_g,'MINPEAKHEIGHT',minPeakHeight(indx));
    numSteps(indx) = numel(pks);
    indx = indx + 1;
end

figure
plot(minPeakHeight, numSteps)
title('Number of peaks Vs. minimum peak height');
xlabel('Minimum peak height');
ylabel('Number of peaks');
grid on

figure
plot(time_matrix, z_matrix_no_g);
title("Z acceleration without gravity vs Time, Number of steps = " + numSteps(13))
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on
hold on
[pks,locs] = findpeaks(z_matrix_no_g,'MINPEAKHEIGHT',minPeakHeight(13));
plot(time_matrix(locs), pks, 'r', 'Marker', 'v', 'LineStyle', 'none');

% Acceleration Magnitude
acc_mag = sqrt(sum(x_matrix.^2 + y_matrix.^2 + z_matrix_no_g.^2,2));
figure
plot(time_matrix, acc_mag);
title('Total Acceleration vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on


