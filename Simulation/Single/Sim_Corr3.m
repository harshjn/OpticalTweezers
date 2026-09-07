clc
clear all
close all
tMat = readtable(['trajectory_N3_alpha0.035_FoverFc0.1_nums100000.csv']);
plot(tMat.particle_1_position)
X1 = tMat.particle_1_position(2100:end);
X2 = tMat.particle_2_position(2100:end);
X3 = tMat.particle_3_position(2100:end);
%X4 = tMat.particle_4_position;

%%
plot(X1); hold on; plot(X2); plot(X3)

Vx1 = diff(X1);
Vx2 = diff(X2);

%%
Vx1_ = Vx1(1:end) - mean(Vx1);
Vx2_ = Vx2(1:end) - mean(Vx2);

% --- Parameters ---
maxShift = 200;
shifts = -maxShift:maxShift;

C = zeros(size(shifts));

parfor k = 1:length(shifts)
    N = shifts(k);
    
    if N >= 0
        v1 = Vx1_(1:end-N);
        v2 = Vx2_(1+N:end);
    else
        v1 = Vx1_(1-N:end);
        v2 = Vx2_(1:end+N);
    end
    
    % Compute truncated variances
    C1 = mean(v1.^2);
    C2 = mean(v2.^2);
    
    % Normalized correlation
    C(k) = mean(v1 .* v2) / sqrt(C1 * C2);
    k
end

% --- Plot ---
figure;
plot(shifts, C, 'k-', 'LineWidth', 1.5);
hold on;
yline(0,'k:');  % zero reference
xlabel('Shift (frames)');
ylabel('C_{xx}(\tau)');
set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','in');
box on;
title('Large Trap, alpha = input, Fc/F = input')