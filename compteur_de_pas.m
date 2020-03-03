% MGA852
% LAB 2
% Compteur de pas
%
% TABLE OF CONTENTS
% 1 - Read CSV
% 2 - If data contains NaN, replace with average
% 3 - Detect Sampling rate
% 4 - Plot time, x, y and z accelerations
% 5 - Detect distance between peak and zero-crossing
% 6 - Detect number of steps by Total acceleration (mathworks)
% 7 - Detect number of steps by zero-crossing
% 8 - Detect number of steps by filtering data
%
% read csv file and create 
% time, x, y, z matrices
clear all
close all
[num, txt, raw] = xlsread('LAB2_Part1_testdata1.csv');
time_matrix     = str2double(txt(2:end, 1));
x_matrix        = str2double(txt(2:end, 2));
y_matrix        = str2double(txt(2:end, 3));
z_matrix        = str2double(txt(2:end, 4));

% there are NaN values
% interpolate them and set them
% to a number
x_matrix_nan = ~isnan(x_matrix);
x_temp = cumsum(x_matrix_nan'-diff([1,x_matrix_nan'])/2);
x_matrix = interp1(1:nnz(x_matrix_nan'),x_matrix(x_matrix_nan'),x_temp)';

y_matrix_nan = ~isnan(y_matrix);
y_temp = cumsum(y_matrix_nan'-diff([1,y_matrix_nan'])/2);
y_matrix = interp1(1:nnz(y_matrix_nan'),y_matrix(y_matrix_nan'),y_temp)';

z_matrix_nan = ~isnan(z_matrix);
z_temp = cumsum(z_matrix_nan'-diff([1,z_matrix_nan'])/2);
z_matrix = interp1(1:nnz(z_matrix_nan'),z_matrix(z_matrix_nan'),z_temp)';
%
% Detect sampling rate
%
time_between_samples = [0];
for i=2:length(time_matrix)
    time_between_samples = [time_between_samples; (time_matrix(i)-time_matrix(i-1))];
end
sampling_rate = nanmean(time_between_samples);
%
% Plot
% X
figure
plot(time_matrix, x_matrix);
title('X vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on
%
% Y
%
figure
plot(time_matrix, y_matrix);
title('Y vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on
%
% Z
% 
figure
plot(time_matrix, z_matrix);
title('Z acceleration vs Time')
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
grid on
%
% Count Number of Steps 
% Acceleration Magnitude method proposed by Mathworks
%
acc_mag = sqrt(x_matrix.^2 + y_matrix.^2 + z_matrix.^2);
acc_mag = acc_mag - nanmean(acc_mag);
minPeakHeight = nanstd(acc_mag);
[pks,locs] = findpeaks(acc_mag,'MINPEAKHEIGHT', minPeakHeight);
numSteps_mat = numel(pks)
figure
plot(time_matrix, acc_mag, 'blue');
hold on
plot(time_matrix(locs), pks, 'r', 'Marker', 'v', 'LineStyle', 'none');
title(['Total Acceleration without gravity vs Time (Mathworks method), Steps = ',numSteps_mat])
xlabel('Time [s]');
ylabel('Total Acceleration [m/s^2]');
grid on
hold off
%
% Calculate average interval peak-to-peak
% discard points where acceleration is lowest
%
p2p_temp = locs*sampling_rate;
for k=2:length(p2p_temp)
    p2p_diff = p2p_temp(k) - p2p_temp(k-1);
end
p2p_avg = nanmean(p2p_diff);
%
% Count number of Steps
% Zero-Crossing method
%
% get a list of all zero-crossings
% check if a peak is found between 2 zero-crossings
%
zero_crossing = [0];
numSteps_zero_x = 0;
zero_positions = [1];
locs_z_x = [1];
pks_z_x = [0];
for i=(2):length(acc_mag)
    zero_crossing = acc_mag(i) * acc_mag(i-1) < 0;
    if zero_crossing 
       zero_positions = [zero_positions; (i)];
    end
end
for j=2:length(zero_positions)
    for k=1:length(locs)
        if (locs(k) >= zero_positions(j-1)) && (locs(k) <= zero_positions(j))
            pks_z_x = [pks_z_x; pks(k)];
            locs_z_x = [locs_z_x; locs(k)];
            numSteps_zero_x = numSteps_zero_x + 1;
            break;
        end
    end
end
numSteps_zero_x
figure
plot(time_matrix, acc_mag, 'blue');
hold on
plot(time_matrix(locs_z_x), pks_z_x, 'r', 'Marker', 'v', 'LineStyle', 'none');
title(['Total Acceleration without gravity vs Time (Zero-Crossing method), Steps = ',numSteps_zero_x])
xlabel('Time [s]');
ylabel('Total Acceleration [m/s^2]');
grid on
hold off
%
% Filter method
% Use the Pythagorean theorem to calculate the magnitude 
% of the acceleration vector of each sample from the accelerometer. 
% Low-pass filter the magnitude signal to remove high frequency noise 
% and then look for peaks and valleys in the filtered signal. 
% You may need to add additional requirements to remove false positives. 
% This is by far the simplest way to detect steps, it is also the 
% way that most if not all ordinary pedometers of the sort 
% that you can buy from a sports store work.
%
freq_spectrum = fft(acc_mag);
L = 8263; % length of signal
Fs = 1000; % sampling frequency
f = Fs*(0:(L/2))/L; % define frequency domain
P2 = abs(freq_spectrum/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
%
figure 
plot(f, P1);
title('Power Vs Frequency')
xlabel('Frequency [Hz]');
ylabel('Power []');
grid on
%
% Construct lowpass filter cut-frequency 8.593Hz
% where a = T/tau, T = the time between samples, and tau is the filter
% time constant
%
a = sampling_rate/(1/8.593);
acc_mag_filt = filter(a, [1 a-1], acc_mag);
% 
minPeakHeight_filt = nanstd(acc_mag_filt);
[pks2,locs2] = findpeaks(acc_mag_filt,'MINPEAKHEIGHT', minPeakHeight_filt);
numSteps_filt = numel(pks2)
%
figure
plot(time_matrix, acc_mag, 'blue');
hold on
plot(time_matrix(locs2), pks2, 'r', 'Marker', 'v', 'LineStyle', 'none');
hold on
plot(time_matrix, acc_mag_filt, 'red');
title('Filtered Total Acceleration vs Time, Steps = ')
xlabel('Time [s]');
ylabel('Filtered Total Acceleration [m/s^2]');
grid on
hold off