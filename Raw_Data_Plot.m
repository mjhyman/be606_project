%import txt files
clear; clc; close all;
% 'Exercise A,B' = exercise 2 on/off. Velocity = 13.4 Km/hr
% 'Exercise C,D' = exercise 3 on/off. velocity = 14.4 Km/hr
% 'Exercise E,F' = exercise 4 on/off. velocity = 15.7 Km/hr
wkdir = 'S1_Dataset';

% Exercise 2
a = importdata(fullfile(wkdir, 'Exercise A.txt'));
b = importdata(fullfile(wkdir, 'Exercise B.txt'));
% Exercise 3
c = importdata(fullfile(wkdir, 'Exercise C.txt'));
d = importdata(fullfile(wkdir, 'Exercise D.txt'));
% Exercise 4
e = importdata(fullfile(wkdir, 'Exercise E.txt'));
f = importdata(fullfile(wkdir, 'Exercise F.txt'));

%resize matrix to remove zeros
at = a(1:962,1);
ahr = a(1:962,2);

bt = b(1:705,1);
bhr = b(1:705,2);

ct = c(1:922,1);
chr = c(1:922,2);

dt = d(1:734,1);
dhr = d(1:734,2);

et = e(1:933,1);
ehr = e(1:933,2);

ft = f(1:786,1);
fhr = f(1:786,2);

%add time from exercise time vector to recovery vector
btnew = bt + at(end,1);
dtnew = dt+ct(end,1);
ftnew = ft + et(end,1);

%append time and heart rate vectors
ex2t = [at; btnew];
ex2hr = [ahr; bhr];
ex3t = [ct; dtnew];
ex3hr = [chr; dhr];
ex4t = [et; ftnew];
ex4hr = [ehr; fhr];

%plot figures
figure
plot(ex2t, ex2hr);
xlabel('Time (s)');
ylabel('hr (beats/min)');

figure 
plot(ex3t,ex3hr);
xlabel('Time (s)');
ylabel('hr (beats/min)');

figure
plot(ex4t,ex4hr);
xlabel('Time (s)');
ylabel('hr (beats/min)');

