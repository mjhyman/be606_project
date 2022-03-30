%import txt files
a = importdata('Exercise A.txt');
b = importdata('Exercise B.txt');
c = importdata('Exercise C.txt');
d = importdata('Exercise D.txt');
e = importdata('Exercise E.txt');
f = importdata('Exercise F.txt');
g = importdata('Exercise G.txt');
h = importdata('Exercise H.txt');

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

gt = g(1:860,1);
ghr = g(1:860,2);

ht = h(1:748,1);
hhr = h(1:748,2);

%add time from exercise time vector to recovery vector
btnew = bt + at(end,1);
dtnew = dt+ct(end,1);
ftnew = ft + et(end,1);
htnew = ht + gt(end,1);

%append time and heart rate vectors
ex2t = [at; btnew];
ex2hr = [ahr; bhr];
ex3t = [ct; dtnew];
ex3hr = [chr; dhr];
ex4t = [et; ftnew];
ex4hr = [ehr; fhr];
ex5t = [gt; htnew];
ex5hr = [ghr; hhr];

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

figure
plot(ex5t,ex5hr);
xlabel('Time (s)');
ylabel('hr (beats/min)');



