%% Method 1
clear;
close all;
clc;

a_L = 0.015; % 1.5cm = 0.015 m
tort = 0.8; %tortuosity
K = 17; %m/day
n = 0.38;
i = 0.0033;
vx = K/n*i; %m/day

disp(vx);

D0_TCE = 10.1 * 10^-6 *3600*24*10^-4; % cm^2/s --> m^2/day
D0_cDCE = 11.4 * 10^-6 *3600*24*10^-4; % cm^2/s --> m^2/day
D0_VC = 13.3 * 10^-6 *3600*24*10^-4; % cm^2/s --> m^2/day

disp(D0_TCE);
disp(D0_cDCE);
disp(D0_VC);

Dx_TCE = vx*a_L + n*tort*D0_TCE ; 
Dx_cDCE = vx*a_L + n*tort*D0_cDCE;
Dx_VC = vx*a_L + n*tort*D0_VC;

disp(Dx_TCE);
disp(Dx_cDCE);
disp(Dx_VC);

disp(2*Dx_TCE/vx);
disp(2*Dx_cDCE/vx);
disp(2*Dx_VC/vx);

dx = 0.03; % dx < 2Dx/vx
dt = 0.03; % dt < (dx)^2  /   2Dx

L1 = 35; % m
L2 = 30; % m 

x1 = 0:dx:L1;
x2 = 0:dx:L2;

Nx1 = length(x1);
Nx2 = length(x2);

Nxi1 = Nx1 +1;
Nxi2 = Nx2 +1;

t = 0:dt:700;
Nt = length(t);

KT_A_V = 7.8 * 10^-4; % day^-1

C1 = zeros(Nxi1,Nt);

C2_T = zeros(Nxi2,Nt);
C2_c = zeros(Nxi2,Nt);
C2_V = zeros(Nxi2,Nt);

T = zeros(1,Nt);
Cw = zeros(1,Nt);

C1(1:Nx1,1) = 10; % C0 = 10000 miug/L = 10 mg/L
C1(1,1) = 0;

C2_T(:,1) = 0;
C2_c(:,1) = 0;
C2_V(:,1) = 0;

C1(Nxi1,1) = C1(Nx1-1,1); % Imaginary node
C2_T(Nxi2,1) = C2_T(Nx2-1,1);
C2_c(Nxi2,1) = C2_c(Nx2-1,1);
C2_V(Nxi2,1) = C2_V(Nx2-1,1);

%Cw(1)=C2_T(Nx2,1)/10;
%T(1)=0;
miu_max_T = 3.42; %day^-1
miu_max_c = 0.16; 
miu_max_V = 0.52;
X_cell = 2*10^9; % cells/L
Ks_T = 10 * 131.38; % mMol/L --> Mol/L --> g/L --> mg/L
Ks_c = 3.3 * 96.94;
Ks_V = 2.6 * 62.5;
Y = 5.2 * 10^8 / 35.45; %cells/mmol --> cells/mg

%Volume of cells
A = 10.67 *35; %m^2
Vcell = n * A * dx * 1000;

V_zone1 = n * A * L1 * 1000;
M_zone1 = zeros(1,Nt); 
M_zone1(1) = 10^9;
source = true;

%% Zone 1 mass transfer and Zone 2 biodegradation
for j =1:(Nt-1)

    C1(Nxi1,j) = C1(Nx1-1,j);

    for i = 2:Nx1
        adv1 = -vx*(C1(i+1,j)-C1(i-1,j))/(2*dx);
        disp1 = Dx_TCE* (C1(i+1,j)-2*C1(i,j)+C1(i-1,j))/((dx)^2);

        if source
            R1 = KT_A_V * (1300 - C1(i,j));
        else
            R1 = 0;
        end

        C1(i,j+1) = C1(i,j)+ dt* (adv1 + disp1+R1);
    end

    
    C1(1,j+1) = 0;
    C1(Nxi1,j+1)=C1(Nx1-1,j+1);


    if source
        dM_dt = V_zone1 * R1;
        M_zone1(j+1) = M_zone1(j)-dM_dt * dt;
        
        if M_zone1(j+1) < 0
            M_zone1(j+1)=0;
        end

        if M_zone1(j+1) <= 0.001*M_zone1(1)
            source = false;
            fprintf('Source stops releasing at t = %.2f day (step %d)\n', t(j+1), j+1);
        end
    else
        M_zone1(j+1) = M_zone1(j);
    end


    C2_T(1,j) = C1(Nx1,j);
    C2_T(1,j+1) = C1(Nx1,j+1);

    C2_c(1,j)   = 0;
    C2_V(1,j)   = 0;
    C2_c(1,j+1) = 0;
    C2_V(1,j+1) = 0;
    
    C2_T(Nxi2,j) = C2_T(Nx2-1,j); % Imaginary node for C2
    C2_c(Nxi2,j) = C2_c(Nx2-1,j);
    C2_V(Nxi2,j) = C2_V(Nx2-1,j);

    for i = 2:Nx2

        adv2_T = -vx*(C2_T(i+1,j)-C2_T(i-1,j))/(2*dx);
        adv2_c = -vx*(C2_c(i+1,j)-C2_c(i-1,j))/(2*dx);
        adv2_V = -vx*(C2_V(i+1,j)-C2_V(i-1,j))/(2*dx);

        disp2_T = Dx_TCE * (C2_T(i+1,j)-2*C2_T(i,j)+C2_T(i-1,j))/((dx)^2);
        disp2_c = Dx_cDCE * (C2_c(i+1,j)-2*C2_c(i,j)+C2_c(i-1,j))/((dx)^2);
        disp2_V = Dx_VC * (C2_V(i+1,j)-2*C2_V(i,j)+C2_V(i-1,j))/((dx)^2);
        
        bio_T = -miu_max_T*C2_T(i,j)*X_cell/((Ks_T+C2_T(i,j))*Y);
        bio_c = -miu_max_c*C2_c(i,j)*X_cell/((Ks_c+C2_c(i,j))*Y);
        bio_V = -miu_max_V*C2_V(i,j)*X_cell/((Ks_V+C2_V(i,j))*Y);
        
        C2_T(i,j+1) = C2_T(i,j) + dt * (adv2_T + disp2_T + bio_T);
        C2_c(i,j+1) = C2_c(i,j) + dt * (adv2_c + disp2_c - bio_T + bio_c);
        C2_V(i,j+1) = C2_V(i,j) + dt * (adv2_V + disp2_V - bio_c + bio_V);

    end
    
    C2_T(Nxi2,j+1) = C2_T(Nx2-1,j+1);
    C2_c(Nxi2,j+1) = C2_c(Nx2-1,j+1);
    C2_V(Nxi2,j+1) = C2_V(Nx2-1,j+1);

end

%% Concentration vs distance
x_full = [x1, L1 + x2(2:end)]; 

t_plot = [10, 50, 200, 600];

for k = 1:length(t_plot)
    [~, jidx] = min(abs(t - t_plot(k)));

    C_T_full = [C1(1:Nx1, jidx); C2_T(2:Nx2, jidx)];
    C_c_full = [zeros(Nx1,1);     C2_c(2:Nx2, jidx)];  % cDCE only in zone 2
    C_V_full = [zeros(Nx1,1);     C2_V(2:Nx2, jidx)];  % VC only in zone 2

    figure;
    plot(x_full, C_T_full, 'r-',  ...
         x_full, C_c_full, 'b--', ...
         x_full, C_V_full, 'g-.', 'LineWidth', 3);
    xlabel('x (m)');
    ylabel('C (mg/L)');
    legend({'TCE','cDCE','VC'}, 'Location','best');
    grid on;
end


%% Concentration vs time: one point in Zone 1 and one point in Zone 2

% Choose a point in zone 1
x1_target = L1 / 2;                      
[~, i1] = min(abs(x1 - x1_target));      
C_T_z1_time = C1(i1,:);                  

% Choose a point in zone 2
x2_target_1 = 2;                      
[~, i2] = min(abs(x2 - x2_target_1));      
C_T_z2_time_1 = C2_T(i2,:);                
C_c_time_1    = C2_c(i2,:);                
C_V_time_1    = C2_V(i2,:);                

x2_target = L2 / 2;                      
[~, i3] = min(abs(x2 - x2_target));      
C_T_z2_time = C2_T(i3,:);                
C_c_time    = C2_c(i3,:);                
C_V_time    = C2_V(i3,:);       

x2_target_2 = 28;                      
[~, i4] = min(abs(x2 - x2_target_2));      
C_T_z2_time_2 = C2_T(i4,:);                
C_c_time_2    = C2_c(i4,:);                
C_V_time_2    = C2_V(i4,:);       

% TCE in zone 1
figure;
plot(t, C_T_z1_time, 'c-',  'LineWidth', 3, ...
    'DisplayName', ['TCE, Zone 1, x = ' num2str(x1(i1),'%.1f') ' m']);
xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
grid on;

% TCE / cDCE / VC in zone 2
% In zone 2 at x=2m
figure;
hold on;
plot(t, C_T_z2_time_1, 'r-',  'LineWidth', 3, ...
    'DisplayName', ['TCE, Zone 2, x = ' num2str(x2(i2),'%.1f') ' m']);
plot(t, C_c_time_1,    'b-', 'LineWidth', 3, ...
    'DisplayName', ['cDCE, Zone 2, x = ' num2str(x2(i2),'%.1f') ' m']);
plot(t, C_V_time_1,    'g-', 'LineWidth', 3, ...
    'DisplayName', ['VC, Zone 2, x = ' num2str(x2(i2),'%.1f') ' m']);
hold off;


xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
grid on;

% In zone 2 at x=15m
figure;
hold on;
plot(t, C_T_z2_time, 'r-',  'LineWidth', 3, ...
    'DisplayName', ['TCE, Zone 2, x = ' num2str(x2(i3),'%.1f') ' m']);
plot(t, C_c_time,    'b-', 'LineWidth', 3, ...
    'DisplayName', ['cDCE, Zone 2, x = ' num2str(x2(i3),'%.1f') ' m']);
plot(t, C_V_time,    'g-', 'LineWidth', 3, ...
    'DisplayName', ['VC, Zone 2, x = ' num2str(x2(i3),'%.1f') ' m']);
hold off;

xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
grid on;

% In zone 2 at x=28m
figure;
hold on;
plot(t, C_T_z2_time_2, 'r-',  'LineWidth', 3, ...
    'DisplayName', ['TCE, Zone 2, x = ' num2str(x2(i4),'%.1f') ' m']);
plot(t, C_c_time_2,    'b-', 'LineWidth', 3, ...
    'DisplayName', ['cDCE, Zone 2, x = ' num2str(x2(i4),'%.1f') ' m']);
plot(t, C_V_time_2,    'g-', 'LineWidth', 3, ...
    'DisplayName', ['VC, Zone 2, x = ' num2str(x2(i4),'%.1f') ' m']);
hold off;

xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
grid on;
%% Mass vs time in Zone 1
figure;
plot(t, M_zone1/1000, 'r-', 'LineWidth', 3);   % mg -> g
xlabel('Time (day)');
ylabel('Total mass in Zone 1 (g)');
grid on;

%% Method 2
clear;
close all;
clc;

a_L = 0.015; % 1.5cm = 0.015 m
tort = 0.8; %tortuosity
K = 17; %m/day
n = 0.38;
i = 0.0033;
vx = K/n*i; %m/day

disp(vx);

D0_TCE = 10.1 * 10^-6 *3600*24*10^-4; % cm^2/s --> m^2/day
D0_cDCE = 11.4 * 10^-6 *3600*24*10^-4; % cm^2/s --> m^2/day
D0_VC = 13.3 * 10^-6 *3600*24*10^-4; % cm^2/s --> m^2/day

disp(D0_TCE);
disp(D0_cDCE);
disp(D0_VC);

Dx_TCE = vx*a_L + n*tort*D0_TCE ; 
Dx_cDCE = vx*a_L + n*tort*D0_cDCE;
Dx_VC = vx*a_L + n*tort*D0_VC;

disp(Dx_TCE);
disp(Dx_cDCE);
disp(Dx_VC);

disp(2*Dx_TCE/vx);
disp(2*Dx_cDCE/vx);
disp(2*Dx_VC/vx);

dx = 0.03; % dx < 2Dx/vx
dt = 0.05; % dt < (dx)^2  /   2Dx

L1 = 35; % m
L2 = 20; % m 
L3 = 10;

x1 = 0:dx:L1;
x2 = 0:dx:L2;
x3 = 0:dx:L3;

Nx1 = length(x1);
Nx2 = length(x2);
Nx3 = length(x3);

Nxi1 = Nx1 +1;
Nxi2 = Nx2 +1;
Nxi3 = Nx3 +1;

t = 0:dt:5000;
Nt = length(t);

KT_A_V = 7.8 * 10^-4; % day^-1

C1 = zeros(Nxi1,Nt);

C2_T = zeros(Nxi2,Nt);

C3_T = zeros(Nxi3,Nt);
C3_c = zeros(Nxi3,Nt);
C3_V = zeros(Nxi3,Nt);

T = zeros(1,Nt);
Cw = zeros(1,Nt);

C1(1:Nx1,1) = 10; % C0 = 10000 miug/L = 10 mg/L
C1(1,1) = 0;

C2_T(:,1) = 0;

C3_T(:,1) = 0;
C3_c(:,1) = 0;
C3_V(:,1) = 0;

C1(Nxi1,1) = C1(Nx1-1,1); % Imaginary node
C2_T(Nxi2,1) = C2_T(Nx2-1,1);

C3_T(Nxi3,1) = C3_T(Nx3-1,1);
C3_c(Nxi3,1) = C3_c(Nx3-1,1);
C3_V(Nxi3,1) = C3_V(Nx3-1,1);

%Cw(1)=C2_T(Nx2,1)/10;
%T(1)=0;
miu_max_T = 3.42; %day^-1
miu_max_c = 0.16; 
miu_max_V = 0.52;
X_cell = 100.0 * 10^8; % cells/L
Ks_T = 10 * 131.38; % mMol/L --> Mol/L --> g/L --> mg/L
Ks_c = 3.3 * 96.94;
Ks_V = 2.6 * 62.5;
Y = 5.2 * 10^8 / 35.45; %cells/mmol --> cells/mg

%Volume of cells
A = 10.67 *35; %m^2
Vcell = n * A * dx * 1000;
M_source0 = 10^9;
M_source = M_source0 ; %mg
source = true;
M_source_record = zeros(1, Nt);
M_source_record(1) = M_source0;

%% Zone 1 mass transfer and Zone 2 biodegradation
for j =1:(Nt-1)

    C1(Nxi1,j) = C1(Nx1-1,j);
    d_mass_release = 0;

    for i = 2:Nx1
        adv1 = -vx*(C1(i+1,j)-C1(i-1,j))/(2*dx);
        disp1 = Dx_TCE* (C1(i+1,j)-2*C1(i,j)+C1(i-1,j))/((dx)^2);

        if source
        
            R1 = KT_A_V * (1300 - C1(i,j));
        else
            R1 = 0;
        end

        C1(i,j+1) = C1(i,j)+ dt* (adv1 + disp1+R1);

        if source && (R1>0)
            d_mass_release = d_mass_release + R1*Vcell*dt;
        end


    end
    
    C1(1,j+1) =0;
    C1(Nxi1,j+1)=C1(Nx1-1,j+1);

    if source
        M_source = M_source - d_mass_release;
        if M_source <= 0.001 * M_source0
            source = false;
            fprintf('Source stops releasing at t = %.2f day (step %d)\n', t(j+1), j+1);
        end
    end

    if source
        M_source_record(j+1) = M_source;
    else
        M_source_record(j+1) = M_source_record(j);
    end

    C2_T(1,j) = C1(Nx1,j);
    C2_T(1,j+1) = C1(Nx1,j+1);

    C2_T(Nxi2,j) = C2_T(Nx2-1,j);

    for i = 2:Nx2
        adv2_T  = -vx * (C2_T(i+1,j) - C2_T(i-1,j)) / (2*dx);
        disp2_T =  Dx_TCE * (C2_T(i+1,j) - 2*C2_T(i,j) + C2_T(i-1,j)) / dx^2;

        C2_T(i,j+1) = C2_T(i,j) + dt * (adv2_T + disp2_T);
    end

    C2_T(Nxi2,j+1) = C2_T(Nx2-1,j+1);

    %zone 3
    C3_T(1,j)   = C2_T(Nx2,j);
    C3_T(1,j+1) = C2_T(Nx2,j+1);

    C3_c(1,j)   = 0;
    C3_V(1,j)   = 0;
    C3_c(1,j+1) = 0;
    C3_V(1,j+1) = 0;
    
    C3_T(Nxi3,j) = C3_T(Nx3-1,j); % Imaginary node for C2
    C3_c(Nxi3,j) = C3_c(Nx3-1,j);
    C3_V(Nxi3,j) = C3_V(Nx3-1,j);

    for i = 2:Nx3

        adv3_T = -vx*(C3_T(i+1,j)-C3_T(i-1,j))/(2*dx);
        adv3_c = -vx*(C3_c(i+1,j)-C3_c(i-1,j))/(2*dx);
        adv3_V = -vx*(C3_V(i+1,j)-C3_V(i-1,j))/(2*dx);

        disp3_T = Dx_TCE * (C3_T(i+1,j)-2*C3_T(i,j)+C3_T(i-1,j))/((dx)^2);
        disp3_c = Dx_cDCE * (C3_c(i+1,j)-2*C3_c(i,j)+C3_c(i-1,j))/((dx)^2);
        disp3_V = Dx_VC * (C3_V(i+1,j)-2*C3_V(i,j)+C3_V(i-1,j))/((dx)^2);
        
        bio_T = -miu_max_T*C3_T(i,j)*X_cell/((Ks_T+C3_T(i,j))*Y);
        bio_c = -miu_max_c*C3_c(i,j)*X_cell/((Ks_c+C3_c(i,j))*Y);
        bio_V = -miu_max_V*C3_V(i,j)*X_cell/((Ks_V+C3_V(i,j))*Y);
        
        C3_T(i,j+1) = C3_T(i,j) + dt * (adv3_T + disp3_T + bio_T);
        C3_c(i,j+1) = C3_c(i,j) + dt * (adv3_c + disp3_c - bio_T + bio_c);
        C3_V(i,j+1) = C3_V(i,j) + dt * (adv3_V + disp3_V - bio_c + bio_V);

    end
    
    C3_T(Nxi3,j+1) = C3_T(Nx3-1,j+1);
    C3_c(Nxi3,j+1) = C3_c(Nx3-1,j+1);
    C3_V(Nxi3,j+1) = C3_V(Nx3-1,j+1);

end


%% 
x_full = [ ...
    x1, ...
    L1 + x2(2:Nx2), ...
    L1 + L2 + x3(2:Nx3) ...
];

% Times to plot
t_plot = [100, 200, 300, 400,500,620];

for k = 1:length(t_plot)
    [~, jidx] = min(abs(t - t_plot(k)));

    % TCE over three zones
    C_T_full = [ ...
        C1(1:Nx1, jidx); ...
        C2_T(2:Nx2, jidx); ...
        C3_T(2:Nx3, jidx) ...
    ];

    % cDCE: only in Zone 3
    C_c_full = [ ...
        zeros(Nx1 + Nx2 - 1, 1); ...
        C3_c(2:Nx3, jidx) ...
    ];

    % VC: only in Zone 3
    C_V_full = [ ...
        zeros(Nx1 + Nx2 - 1, 1); ...
        C3_V(2:Nx3, jidx) ...
    ];

    figure;
    plot(x_full, C_T_full, 'r-',  ...
         x_full, C_c_full, 'b--', ...
         x_full, C_V_full, 'g-.', 'LineWidth', 3);
    xlabel('x (m)');
    ylabel('C (mg/L)');
    legend({'TCE','cDCE','VC'}, 'Location','best');
    title(['t = ' num2str(t(jidx)) ' day']);
    grid on;
end

%% Concentration vs time at selected points

% One point in Zone 1 (middle)
x1_target = L1 / 2;
[~, i1] = min(abs(x1 - x1_target));
C_T_z1_time = C1(i1, :);

figure;
plot(t, C_T_z1_time, 'c-', 'LineWidth', 3, ...
    'DisplayName', ['TCE, Zone 1, x = ' num2str(x1(i1),'%.1f') ' m']);
xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
title('TCE in Zone 1 (middle point)');
xlim([0 700]);
grid on;

% One point in Zone 2 (middle)
x2_target = L2 / 2;
[~, i2] = min(abs(x2 - x2_target));
C_T_z2_time = C2_T(i2, :);

figure;
plot(t, C_T_z2_time, 'm-', 'LineWidth', 3, ...
    'DisplayName', ['TCE, Zone 2, x = ' num2str(x2(i2),'%.1f') ' m']);
xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
title('TCE in Zone 2 (middle point)');
xlim([0 700]);
grid on;

% Zone 3: near inlet, e.g. x = 2 m
x3_target_1 = 0.3;
[~, i3] = min(abs(x3 - x3_target_1));
C_T_z3_time_1 = C3_T(i3, :);
C_c_time_1    = C3_c(i3, :);
C_V_time_1    = C3_V(i3, :);

figure;
hold on;
plot(t, C_T_z3_time_1, 'r-', 'LineWidth', 3, ...
    'DisplayName', ['TCE']);
plot(t, C_c_time_1,    'b-', 'LineWidth', 3, ...
    'DisplayName', ['cDCE']);
plot(t, C_V_time_1,    'g-', 'LineWidth', 3, ...
    'DisplayName', ['VC']);
hold off;
xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
title('Species in Zone 3 near inlet');
xlim([0 700]);

grid on;

% Zone 3: middle point, x = L3/2
x3_target_2 = L3 / 2;
[~, i4] = min(abs(x3 - x3_target_2));
C_T_z3_time_2 = C3_T(i4, :);
C_c_time_2    = C3_c(i4, :);
C_V_time_2    = C3_V(i4, :);

figure;
hold on;
plot(t, C_T_z3_time_2, 'r-', 'LineWidth', 3, ...
    'DisplayName', ['TCE']);
plot(t, C_c_time_2,    'b-', 'LineWidth', 3, ...
    'DisplayName', ['cDCE']);
plot(t, C_V_time_2,    'g-', 'LineWidth', 3, ...
    'DisplayName', ['VC']);
hold off;
xlabel('Time (day)');
ylabel('C (mg/L)');
legend('Location','best');
title('Species in Zone 3 (middle point)');
xlim([0 700]);

grid on;

%% Source mass vs time in zone 1

figure;
plot(t, M_source_record / 1000, 'r-', 'LineWidth', 3); % mg -> g
xlabel('Time (day)');
ylabel('Source mass (g)');
title('Source mass vs time');
xlim([0 500]);

grid on;