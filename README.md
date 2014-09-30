Phase-I
=======

%Hodgkin Huxley Model 

%Question1
%Set up the constants
g_K = 36;  %maximum conductance of potassium
g_Na = 120; %maximum conductance of sodium
g_L = 0.3; %maximum conductance
Ek = -12; % potassium voltage across membrane
Ena = 115; % sodium voltage across membrane
El = 10.6; %Voltage across membrance
Vrest = -70; %resting potential
Cm = 1.0; % membrance capacitance
I =0; %initialize current. steay-state --> 0

%initialize membrane potential 
Vm = 0; %function of the membrane currents
stepsize = 0.01; 

X = [0:0.01:100]; %initialize x-axis for part(a)
V_Vec = [Vm]; %initialize y-axis for part (a)

%Gating variables
am = 0.1 * ((25-Vm)/(exp((25-Vm)/10)-1));
bm = 4 * exp((-Vm/18));
an = 0.01 * ((10-Vm)/ (exp((10-Vm)/10)-1));
bn = 0.125 * exp(-Vm/80);
ah = 0.07 * exp(-Vm/20);
bh = 1/(exp((30-Vm)/10)+1);

%initial conditions for m, n, and h.(when t= 0)
m = am / (am + bm); %activation probability
n = an / (an + bn); %activation probability
h = ah / (ah + bh); %inactivation probability

%initialize sodium and potassium conductance
gNa = (m.^3) * h * g_Na; %equation for sodium conductance
gK = (n.^4) * g_K; %equation for potassium conductance
gNa_Vec = [gNa];%initialize conductance of potassium for part(b)
gK_Vec = [gK];%initialize conductance of sodium for part(b)

%initialize the time
t = 0;
%LOOP
while t <= 100

%Gating variables    
am = 0.1 * ((25-Vm)/(exp((25-Vm)/10)-1));
bm = 4 * exp((-Vm/18));
an = 0.01 * ((10-Vm)/ (exp((10-Vm)/10)-1));
bn = 0.125 * exp(-Vm/80);
ah = 0.07 * exp(-Vm/20);
bh = 1/(exp((30-Vm)/10)+1);

%Currents
Ina = (m^3) * g_Na * h * (Vm-Ena);
Ik = (n^4) * g_K *(Vm-Ek);
Il = g_L * (Vm - El);
Iion = I-Ik-Ina - Il;

%New voltage value 
Vm = Vm + stepsize*(Iion / Cm);
%building the vector for y-axis (adding the new voltage value)
V_Vec = [V_Vec Vm];

%next step probability
m = m + stepsize.*((am *(1-m)) - (bm*m));
n = n + stepsize.*((an* (1-n)) - (bn*n));
h = h + stepsize.*((ah *(1-h))- (bh*h));

%calculating conductance with next values
gNa =  (m.^3) * h * g_Na;
gK = (n.^4) * g_K;
gNa_Vec = [gNa_Vec gNa];
gK_Vec = [gK_Vec gK];

t = t + 0.01;

end

%subtracting voltage to set to the resting potential value
V_Vec = V_Vec -70;
%plotting the voltage vs. time 
plot (X,V_Vec);
title('membrane potential');
axis([0,100,-100,40])
xlabel('time(ms)');
ylabel('Voltage (mV)');

%plotting conductances vs. time
figure
plot(X,gNa_Vec, 'r');
hold on
plot(X,gK_Vec, 'b');
legend ('gNa', 'gK');
title('gK and gNa');
axis([0,100,0,40]);
xlabel('Time(ms)');
ylabel('Conductance(S/cm2)');

%% Question2
clear
%Set up the constants
g_K = 36;  %maximum conductance of potassium
g_Na = 120; %maximum conductance of sodium
g_L = 0.3; %maximum conductance
Ek = -12; % potassium voltage across membrane
Ena = 115; % sodium voltage across membrane
El = 10.6; %Voltage across membrance
Vrest = -70; %resting potential
Cm = 1.0; % membrance capacitance
I = 0;

%initialize membrane potential 
Vm = 0; %function of the membrane currents
stepsize = 0.01; 

X = [0:0.01:100]; %initialize x-axis for part(a)
V_Vec = [Vm]; %initialize y-axis for part (a)

%Gating variables
am = 0.1 * ((25-Vm)/(exp((25-Vm)/10)-1));
bm = 4 * exp((-Vm/18));
an = 0.01 * ((10-Vm)/ (exp((10-Vm)/10)-1));
bn = 0.125 * exp(-Vm/80);
ah = 0.07 * exp(-Vm/20);
bh = 1/(exp((30-Vm)/10)+1);

%initial conditions for m, n, and h.(when t= 0)
m = am / (am + bm); %activation probability
n = an / (an + bn); %activation probability
h = ah / (ah + bh); %inactivation probability

%initialize sodium and potassium conductance
gNa = (m.^3) * h * g_Na; 
gK = (n.^4) * g_K;
gNa_Vec = [gNa];%initialize conductance of potassium for part(b)
gK_Vec = [gK];%initialize conductance of sodium for part(b)

t = 0

while t <= 100
if t < 0.5
    I = 5;
else 
    I = 0;
end

%Gating variables    
am = 0.1 * ((25-Vm)/(exp((25-Vm)/10)-1));
bm = 4 * exp((-Vm/18));
an = 0.01 * ((10-Vm)/ (exp((10-Vm)/10)-1));
bn = 0.125 * exp(-Vm/80);
ah = 0.07 * exp(-Vm/20);
bh = 1/(exp((30-Vm)/10)+1);

%Currents
Ina = (m^3) * g_Na * h * (Vm-Ena);
Ik = (n^4) * g_K *(Vm-Ek);
Il = g_L * (Vm - El);
Iion = I-Ik-Ina - Il;

%New voltage value 
Vm = Vm + stepsize*(Iion / Cm);
%building the vector for y-axis (adding the new voltage value)
V_Vec = [V_Vec Vm];

%next step probability
m = m + stepsize.*((am *(1-m)) - (bm*m));
n = n + stepsize.*((an* (1-n)) - (bn*n));
h = h + stepsize.*((ah *(1-h))- (bh*h));

%calculating conductance with next values
gNa =  (m.^3) * h * g_Na;
gK = (n.^4) * g_K;
gNa_Vec = [gNa_Vec gNa];
gK_Vec = [gK_Vec gK];
 
t = t + 0.01;

end

%subtracting voltage to set to the resting potential value
V_Vec = V_Vec -70;
%plotting the voltage vs. time 
plot (X,V_Vec);
title('membrane potential');
axis([0,100,-100,40])
xlabel('time(ms)');
ylabel('Voltage (mV)');

%plotting conductances vs. time
figure
plot(X,gNa_Vec, 'r');
hold on
plot(X,gK_Vec, 'b');
legend ('gNa', 'gK');
title('gK and gNa');
axis([0,100,0,40]);
xlabel('Time(ms)');
ylabel('Conductance(S/cm2)');

%% Question3
clear

%Question1
%Set up the constants
g_K = 36;  %maximum conductance of potassium
g_Na = 120; %maximum conductance of sodium
g_L = 0.3; %maximum conductance
Ek = -12; % potassium voltage across membrane
Ena = 115; % sodium voltage across membrane
El = 10.6; %Voltage across membrance
Vrest = -70; %resting potential
Cm = 1.0; % membrance capacitance
I =5; %initialize current. steay-state --> 0

%initialize membrane potential 
Vm = 0; %function of the membrane currents
stepsize = 0.01; 

X = [0:0.01:100]; %initialize x-axis for part(a)
V_Vec = [Vm]; %initialize y-axis for part (a)

%Gating variables
am = 0.1 * ((25-Vm)/(exp((25-Vm)/10)-1));
bm = 4 * exp((-Vm/18));
an = 0.01 * ((10-Vm)/ (exp((10-Vm)/10)-1));
bn = 0.125 * exp(-Vm/80);
ah = 0.07 * exp(-Vm/20);
bh = 1/(exp((30-Vm)/10)+1);

%initial conditions for m, n, and h.(when t= 0)
m = am / (am + bm); %activation probability
n = an / (an + bn); %activation probability
h = ah / (ah + bh); %inactivation probability

%initialize sodium and potassium conductance
gNa = (m.^3) * h * g_Na; 
gK = (n.^4) * g_K;
gNa_Vec = [gNa];%initialize conductance of potassium for part(b)
gK_Vec = [gK];%initialize conductance of sodium for part(b)

t = 0;

%LOOP
while t <= 100

%Gating variables    
am = 0.1 * ((25-Vm)/(exp((25-Vm)/10)-1));
bm = 4 * exp((-Vm/18));
an = 0.01 * ((10-Vm)/ (exp((10-Vm)/10)-1));
bn = 0.125 * exp(-Vm/80);
ah = 0.07 * exp(-Vm/20);
bh = 1/(exp((30-Vm)/10)+1);

%Currents
Ina = (m^3) * g_Na * h * (Vm-Ena);
Ik = (n^4) * g_K *(Vm-Ek);
Il = g_L * (Vm - El);
Iion = I-Ik-Ina - Il;

%New voltage value 
Vm = Vm + stepsize*(Iion / Cm);
%building the vector for y-axis (adding the new voltage value)
V_Vec = [V_Vec Vm];

%next step probability
m = m + stepsize.*((am *(1-m)) - (bm*m));
n = n + stepsize.*((an* (1-n)) - (bn*n));
h = h + stepsize.*((ah *(1-h))- (bh*h));

%calculating conductance with next values
gNa =  (m.^3) * h * g_Na;
gK = (n.^4) * g_K;
gNa_Vec = [gNa_Vec gNa];
gK_Vec = [gK_Vec gK];

t = t + 0.01;

end

%subtracting voltage to set to the resting potential value
V_Vec = V_Vec -70;
%plotting the voltage vs. time 
plot (X,V_Vec);
title('membrane potential');
axis([0,100,-100,40])
xlabel('time(ms)');
ylabel('Voltage (mV)');

%plotting conductances vs. time
figure
plot(X,gNa_Vec, 'r');
hold on
plot(X,gK_Vec, 'b');
legend ('gNa', 'gK');
title('gK and gNa');
axis([0,100,0,40]);
xlabel('Time(ms)');
ylabel('Conductance(S/cm2)');





