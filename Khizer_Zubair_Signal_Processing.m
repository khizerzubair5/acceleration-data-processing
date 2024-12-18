clc
clear all
close all

%% Foot Acceleration Ensemble Averaging
% Load files into tables
footTrial1 = readtable('Khizer_Trial_1(foot).txt');
footTrial2 = readtable('Khizer_Trial_2(foot).txt');

% Foot Trial 1
timevector1 = second(footTrial1.time);
dt1 = (timevector1(end)-timevector1(1)) / length(timevector1);
t1 = 0:dt1:(timevector1(end)-timevector1(1))-dt1;   % We won't use this
trial1_acc_matrix = [footTrial1.AccX_g_ footTrial1.AccY_g_ footTrial1.AccZ_g_];

% Foot Trial 2
timevector2 = second(footTrial2.time);
dt2 = (timevector2(end)-timevector2(1)) / length(timevector2);
t2 = 0:dt2:(timevector2(end)-timevector2(1))-dt2;
trial2_acc_matrix = [footTrial2.AccX_g_ footTrial2.AccY_g_ footTrial2.AccZ_g_];

% Shortening Trial 1
trial1_acc_matrix = trial1_acc_matrix(1:end-2,:);

% Making a combined matrix for both trials in xyz directions
combined_acc_matrix = [trial1_acc_matrix trial2_acc_matrix];

%Ensemble Averaging
ensemble_averageX = (combined_acc_matrix(:,1) + combined_acc_matrix(:,4)) / 2;  % DIVIDE BY HOW MANY TRIALS YOU HAVE
ensemble_averageY = (combined_acc_matrix(:,2) + combined_acc_matrix(:,5)) / 2;
ensemble_averageZ = (combined_acc_matrix(:,3) + combined_acc_matrix(:,6)) / 2;
ensemble_averaged_data_matrix = [ensemble_averageX ensemble_averageY ensemble_averageZ];


figure;
hold on;
plot(t2,ensemble_averaged_data_matrix(:,1),LineWidth=2) % Ensemble Averaged Data in x-direction
plot(t2,ensemble_averaged_data_matrix(:,2),LineWidth=2) % y-direction
plot(t2,ensemble_averaged_data_matrix(:,3),LineWidth=2) % z-direction
hold off

xlabel('Time(s)')
ylabel('Acceleration(g)')
title('Ensemble Averaged Foot Acceleration Data')
legend('AccX','AccY','AccZ')

%% Thigh Acceleration Ensemble Averaging
thighTrial1 = readtable('Khizer_Trial_1(thigh).txt');
thighTrial2 = readtable('Khizer_Trial_2(thigh).txt');

% Thigh Trial 1
datetime_vector1_thigh = second(thighTrial1.time);
dt1_thigh = (datetime_vector1_thigh(end)-datetime_vector1_thigh(1)) / length(datetime_vector1_thigh);
t1_thigh = 0:dt1_thigh:(datetime_vector1_thigh(end)-datetime_vector1_thigh(1))-dt1_thigh;   % We won't use this
thigh_trial1_acc_matrix = [thighTrial1.AccX_g_ thighTrial1.AccY_g_ thighTrial1.AccZ_g_];

% Thigh Trial 2
datetime_vector2_thigh = second(thighTrial2.time);
dt2_thigh = (datetime_vector2_thigh(end)-datetime_vector2_thigh(1)) / length(datetime_vector2_thigh);
t2_thigh = 0:dt2_thigh:(datetime_vector2_thigh(end)-datetime_vector2_thigh(1))-dt2_thigh;
thigh_trial2_acc_matrix = [thighTrial2.AccX_g_ thighTrial2.AccY_g_ thighTrial2.AccZ_g_];

% Shortening Trial 1
thigh_trial1_acc_matrix = thigh_trial1_acc_matrix(1:end-1,:);     % Taking off 1 row from trial 1

% Making a combined matrix for both trials in xyz directions
thigh_combined_acc_matrix = [thigh_trial1_acc_matrix thigh_trial2_acc_matrix];

%Ensemble Averaging
thigh_ensemble_averageX = (thigh_combined_acc_matrix(:,1) + thigh_combined_acc_matrix(:,4)) / 2;  % DIVIDE BY HOW MANY TRIALS YOU HAVE
thigh_ensemble_averageY = (thigh_combined_acc_matrix(:,2) + thigh_combined_acc_matrix(:,5)) / 2;
thigh_ensemble_averageZ = (thigh_combined_acc_matrix(:,3) + thigh_combined_acc_matrix(:,6)) / 2;
thigh_ensemble_averaged_data_matrix = [thigh_ensemble_averageX thigh_ensemble_averageY thigh_ensemble_averageZ];



figure;
hold on;
plot(t2_thigh,thigh_ensemble_averaged_data_matrix(:,1),LineWidth=2) % Ensemble Averaged Data in x-direction
plot(t2_thigh,thigh_ensemble_averaged_data_matrix(:,2),LineWidth=2) % y-direction
plot(t2_thigh,thigh_ensemble_averaged_data_matrix(:,3),LineWidth=2) % z-direction
hold off

xlabel('Time(s)')
ylabel('Acceleration(g)')
title('Ensemble Averaged Thigh Acceleration Data')
legend('AccX','AccY','AccZ')

%% Foot Acceleration Data Filtering
fs = length(trial1_acc_matrix) / (timevector2(end)-timevector2(1)); %sampling freq
fc = 3;   % cutoff freq
normal_fc = fc / (fs/2);
[b,a] = butter(4,normal_fc,'low');    % Filter design

%Filter Application
filtered_foot_data_X = filtfilt(b,a,ensemble_averageX);
filtered_foot_data_Y = filtfilt(b,a,ensemble_averageY);
filtered_foot_data_Z = filtfilt(b,a,ensemble_averageZ);

figure;
%subplot(3,1,1);
plot(t2,filtered_foot_data_X,LineWidth=2)
hold on;
plot(t2,filtered_foot_data_Y,LineWidth=2)
plot(t2,filtered_foot_data_Z,LineWidth=2)
legend('AccX','AccY','AccZ');
xlabel('Time (s)');
ylabel('Acceleration (g)');
title('Khizer Filtered Foot Acceleration versus Time');

%% Thigh Acceleration Data Filtering
fs = length(thigh_trial1_acc_matrix) / (datetime_vector2_thigh(end)-datetime_vector2_thigh(1)); %sampling freq
fc = 3;   % cutoff freq
normal_fc = fc / (fs/2);
[b,a] = butter(4,normal_fc,'low');    % Filter design

%Filter Application
filtered_thigh_data_X = filtfilt(b,a,thigh_ensemble_averageX);
filtered_thigh_data_Y = filtfilt(b,a,thigh_ensemble_averageY);
filtered_thigh_data_Z = filtfilt(b,a,thigh_ensemble_averageZ);

figure;
plot(t2_thigh,filtered_thigh_data_X,LineWidth=2)
hold on;
plot(t2_thigh,filtered_thigh_data_Y,LineWidth=2)
plot(t2_thigh,filtered_thigh_data_Z,LineWidth=2)
legend('AccX','AccY','AccZ');
xlabel('Time (s)');
ylabel('Acceleration (g)');
title('Khizer Filtered Thigh Acceleration versus Time');