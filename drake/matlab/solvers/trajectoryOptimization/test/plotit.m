
dvar = load('variational-rms.mat');
dfo = load('fo-rms.mat');

figure(1)
subplot(1,2,1);
plot(dvar.N,dvar.rms_pos,'b.-');
hold on;
plot(dfo.N,dfo.rms_pos,'r.-');
hold off;

subplot(1,2,2);
plot(dvar.N,dvar.rms_vel,'b.-');
hold on;
plot(dfo.N,dfo.rms_vel,'r.-');
hold off;
