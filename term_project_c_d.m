%loc = [0, 0.3, 0.7, 1.1, 1.2, 2.2, 2.3, 2.9, 3.1, 3.3, 3.4, 4, 4.2, 4.3, 4.9, 6, 6.5, 6.7, 6.9, 7.1, 7.6];
%m = [0, 20, 25, 25, 20, 25, 70, 1, 10, 15, 1, 10, 120, 1, 50, 5, 1, 5, 5, 10, 10];
clear; clc; close all;

loc = [0, 0.5, 6.5, 11.5, 13, 15]; 
mass = [0, 20, 120, 1, 50, 5]; 

dim = [1.9, 6.0, 0.8, 1.0, 0.2]; 

tolerance_dim = 0.002; 
tolerance_mass = 0.1; 

total_mass_nom = sum(mass);
CG_nominal = sum(loc .* mass) / total_mass_nom;

L_total = sum(dim);
L_nose = dim(1);
CP_nose = 0.67 * L_nose;
CN_nose = 2.0;
CP_body = L_total / 2;
CN_body = 0.5;
CP_fins = L_total - 0.15;
CN_fins = 3.5;
CP_nominal = (CN_nose * CP_nose + CN_body * CP_body + CN_fins * CP_fins) / ...
             (CN_nose + CN_body + CN_fins);

fprintf('=== NOMINAL ROCKET PARAMETERS ===\n');
fprintf('Total length: %.2f m\n', L_total);
fprintf('Total mass: %.2f kg\n', total_mass_nom);
fprintf('CG location: %.4f m from nose\n', CG_nominal);
fprintf('CP location: %.4f m from nose\n', CP_nominal);
fprintf('Static stability margin: %.4f m\n\n', CP_nominal - CG_nominal);

fprintf('=== PART (c): Monte Carlo Trajectory Simulation Model ===\n');
fprintf('Assumption: CG location remains CONSTANT (no change due to fuel consumption)\n');
fprintf('Using ASE 367K Flight Dynamics trajectory equations\n\n');

g0 = 9.807; 
R_earth = 6371e3;
mu = 3.986e14; 
Omega_earth = 7.29e-5;
latitude = 30.6280; 
target_altitude = 35000;

Isp = 250; 
thrust = 50000;
m_fuel = 120; 
burn_time = 60; 
burn_rate = m_fuel / burn_time; 

Cd = 0.4; 
diameter = 1.0; 
A_ref = pi * (diameter/2)^2; 

rho_sealevel = 1.225; 
beta_atm = 1/9042; 

n_simulations = 5000;
fprintf('Running %d Monte Carlo simulations...\n', n_simulations);

max_altitudes = zeros(n_simulations, 1);
max_velocities = zeros(n_simulations, 1);
max_accelerations = zeros(n_simulations, 1);
max_dynamic_pressure = zeros(n_simulations, 1);
CG_samples = zeros(n_simulations, 1);
CP_samples = zeros(n_simulations, 1);
stability_margins = zeros(n_simulations, 1);

fprintf('Progress: ');
progress_step = floor(n_simulations / 20);

dt = 0.1; 
t_max = 400; 
time_vec = 0:dt:t_max;
n_time = length(time_vec);

for sim = 1:n_simulations
    dim_varied = dim + (tolerance_dim/3) * randn(size(dim));
    mass_varied = mass + (tolerance_mass/3) * randn(size(mass));
    mass_varied(mass_varied < 0) = 0; 
    
    total_mass_varied = sum(mass_varied);
    CG_samples(sim) = sum(loc .* mass_varied) / total_mass_varied;
    
    L_total_var = sum(dim_varied);
    L_nose_var = dim_varied(1);
    CP_nose_var = 0.67 * L_nose_var;
    CP_body_var = L_total_var / 2;
    CP_fins_var = L_total_var - 0.15;
    CP_samples(sim) = (CN_nose * CP_nose_var + CN_body * CP_body_var + CN_fins * CP_fins_var) / ...
                      (CN_nose + CN_body + CN_fins);
    
    stability_margins(sim) = CP_samples(sim) - CG_samples(sim);
    
    m_total = total_mass_varied;
    
    x = zeros(n_time, 1);
    z = zeros(n_time, 1);
    vx = zeros(n_time, 1);
    vz = zeros(n_time, 1);
    altitude = zeros(n_time, 1);
    velocity = zeros(n_time, 1);
    acceleration = zeros(n_time, 1);
    gamma_vec = zeros(n_time, 1);
    dynamic_pressure = zeros(n_time, 1);
    
    x(1) = 0;
    z(1) = 0;
    vx(1) = Omega_earth * R_earth * cosd(latitude); 
    vz(1) = 0;
    m = m_total;
    
    gamma_rad = pi/2;
    
    for i = 1:n_time-1
        t = time_vec(i);
        h = z(i); 
        
        if h < 0
            rho = rho_sealevel;
        else
            rho = rho_sealevel * exp(-beta_atm * h);
        end
        
        r = sqrt((R_earth + h)^2);
        
        g_grav = mu / r^2;
        
        g_centripetal_x = Omega_earth^2 * (R_earth + h) * cosd(latitude);
        
        V = sqrt(vx(i)^2 + vz(i)^2);
        velocity(i) = V;
        
        if V > 0.1
            gamma_rad = atan2(vz(i), vx(i));
        end
        gamma_vec(i) = gamma_rad * 180/pi;
        
        dynamic_pressure(i) = 0.5 * rho * V^2;
        
        if t <= burn_time
            F_thrust_x = thrust * cos(gamma_rad);
            F_thrust_z = thrust * sin(gamma_rad);
            m = m_total - burn_rate * t;
        else
            F_thrust_x = 0;
            F_thrust_z = 0;
            m = m_total - m_fuel;
        end
        
        if V > 0.1
            F_drag_x = -0.5 * rho * V * vx(i) * Cd * A_ref;
            F_drag_z = -0.5 * rho * V * vz(i) * Cd * A_ref;
        else
            F_drag_x = 0;
            F_drag_z = 0;
        end
        
        F_grav_x = -m * g_centripetal_x;
        F_grav_z = -m * g_grav;
        
        Fx = F_thrust_x + F_drag_x + F_grav_x;
        Fz = F_thrust_z + F_drag_z + F_grav_z;
        
        ax = Fx / m;
        az = Fz / m;
        acceleration(i) = sqrt(ax^2 + az^2);
        
        vx(i+1) = vx(i) + ax * dt;
        vz(i+1) = vz(i) + az * dt;
        
        x(i+1) = x(i) + vx(i) * dt;
        z(i+1) = z(i) + vz(i) * dt;
        
        if z(i+1) < 0 && i > 10
            altitude(i+1:end) = 0;
            velocity(i+1:end) = 0;
            acceleration(i+1:end) = 0;
            gamma_vec(i+1:end) = gamma_vec(i);
            dynamic_pressure(i+1:end) = 0;
            break;
        end
        
        altitude(i+1) = z(i+1);
        
    end
    
    max_altitudes(sim) = max(altitude);
    max_velocities(sim) = max(velocity);
    max_accelerations(sim) = max(acceleration);
    max_dynamic_pressure(sim) = max(dynamic_pressure);
    
    if sim == 1
        time_ex = time_vec;
        alt_ex = altitude;
        vel_ex = velocity;
        acc_ex = acceleration;
        gamma_ex = gamma_vec;
        q_ex = dynamic_pressure;
    end
    
    if mod(sim, progress_step) == 0
        fprintf('█');
    end
end

fprintf('=== PART (b): Wind Output at given Input Altitude ===\n');
function [ug, vg, wg] = wind_turb(alt, sigma_gust, V_fs, dt, t_max)
    %alt units in feet
    time_vec = 0:dt:t_max;
    n_time = length(time_vec);
    ug = zeros(n_time, 1);
    vg = zeros(n_time, 1);
    wg = zeros(n_time, 1);
    alt = 0.3048*alt; %feet to meter
    
    for i = 1:n_time-1
        h = alt(i);
        if h <= 1000
            Lw = h; 
            Lu = h / (0.177 + 0.000823*h)^1.2; 
            Lv = Lu;
            sigmaw = sigma_gust;
            sigmau = sigma_gust / (0.177 + 0.000823*h)^0.4;
            sigmav = sigmau;
        elseif h >= 2000
            Lw = 1750; 
            Lu = Lw; 
            Lv = Lw;
            sigmaw = sigma_gust;
            sigmau = sigmaw;
            sigmav = sigmau;
        elseif h > 1000 && h < 2000
            h1 = 1000; 
            Lw1     = h1;
            Lu1     = h1 / (0.177 + 0.000823*h1)^1.2;
            Lv1     = Lu1;
            sigmaw1 = sigma_gust;
            sigmau1 = sigma_gust / (0.177 + 0.000823*h1)^0.4;
            sigmav1 = sigmau1;

            h2 = 2000;
            Lw2     = 1750;
            Lu2     = 1750;
            Lv2     = 1750;        
            sigmaw2 = sigma_gust;
            sigmau2 = sigma_gust;
            sigmav2 = sigma_gust;

            slope = (h - h1) / (h2 - h1);   % ranges from 0 to 1
        
            Lw = Lw1 + slope*(Lw2 - Lw1);
            Lu = Lu1 + slope*(Lu2 - Lu1);
            Lv = Lv1 + slope*(Lv2 - Lv1);
        
            sigmaw = sigmaw1 + slope*(sigmaw2 - sigmaw1);
            sigmau = sigmau1 + slope*(sigmau2 - sigmau1);
            sigmav = sigmav1 + slope*(sigmav2 - sigmav1);
        else
            warning('Height is not in range for wind model');
        end
        V_dt = V_fs*dt;

        ug(i+1) = (1 - V_dt/Lu)*ug(i) + sigmau * sqrt(2*V_dt/Lu) * randn; 
        vg(i+1) = (1 - V_dt/Lv)*vg(i) + sigmav * sqrt(2*V_dt/Lv) * randn; 
        wg(i+1) = (1 - V_dt/Lw)*wg(i) + sigmaw * sqrt(2*V_dt/Lw) * randn;
    end
end
sigma_gust = 0.44704*46*rand; 
%mph to m/s multiplied by 46mph multiplied by rand 
%gives random wind speed between "light air" to "gale"
%source - www.weather.gov
[ug, vg, wg] = wind_turb(altitude, sigma_gust, velocity, dt, t_max);

vx = vx - ug;
vy = -vg;
vz = vz - wg;
V = sqrt(vx.^2 + vy.^2 + vz.^2);
velocity = V;

figure;
wind_gust = sqrt(ug.^2 + vg.^2 + wg.^2);
plot(altitude, wind_gust, 'b-', 'LineWidth', 2);
ylabel('Altitude (km)');
xlabel('Wind Gust (km/s)');
title('Part B: Wind Gusts varying with Altitude');
grid on;

fprintf(' Done!\n\n');

fprintf('TRAJECTORY RESULTS:\n');
fprintf('  Maximum Altitude:\n');
fprintf('    Mean: %.2f km\n', mean(max_altitudes)/1000);
fprintf('    Std Dev: %.2f km\n', std(max_altitudes)/1000);
fprintf('    Range: [%.2f, %.2f] km\n', min(max_altitudes)/1000, max(max_altitudes)/1000);
fprintf('    Target: %.2f km\n', target_altitude/1000);

success_count = sum(max_altitudes >= target_altitude);
success_rate = success_count / n_simulations * 100;
fprintf('    Success Rate: %.2f%% (%d/%d reach target)\n\n', ...
    success_rate, success_count, n_simulations);

fprintf('  Maximum Velocity:\n');
fprintf('    Mean: %.2f m/s\n', mean(max_velocities));
fprintf('    Range: [%.2f, %.2f] m/s\n\n', min(max_velocities), max(max_velocities));

fprintf('  Maximum Acceleration:\n');
fprintf('    Mean: %.2f g\n', mean(max_accelerations)/g0);
fprintf('    Max: %.2f g\n\n', max(max_accelerations)/g0);

fprintf('  Maximum Dynamic Pressure:\n');
fprintf('    Mean: %.2f kPa\n', mean(max_dynamic_pressure)/1000);
fprintf('    Max: %.2f kPa\n\n', max(max_dynamic_pressure)/1000);

figure('Name', 'Part (c) - Monte Carlo Trajectory Results', 'Position', [50 50 1400 900]);

subplot(3,3,1);
histogram(max_altitudes/1000, 50, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
hold on;
xline(target_altitude/1000, 'r--', 'LineWidth', 2.5, 'Label', 'Target');
xlabel('Maximum Altitude (km)'); 
ylabel('Frequency');
title('Distribution of Maximum Altitudes');
grid on;

subplot(3,3,2);
sorted_alt = sort(max_altitudes/1000);
cdf_y = (1:n_simulations)' / n_simulations * 100;
plot(sorted_alt, cdf_y, 'b-', 'LineWidth', 2);
hold on;
xline(target_altitude/1000, 'r--', 'LineWidth', 2);
xlabel('Maximum Altitude (km)'); 
ylabel('Cumulative Probability (%)');
title('CDF of Maximum Altitude');
grid on;

subplot(3,3,3);
histogram(max_velocities, 50, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');
xlabel('Maximum Velocity (m/s)'); 
ylabel('Frequency');
title('Distribution of Maximum Velocity');
grid on;

subplot(3,3,4);
histogram(max_accelerations/g0, 50, 'FaceColor', [0.6 0.2 0.8], 'EdgeColor', 'none');
xlabel('Maximum Acceleration (g)'); 
ylabel('Frequency');
title('Distribution of Maximum Acceleration');
grid on;

subplot(3,3,5);
histogram(max_dynamic_pressure/1000, 50, 'FaceColor', [0.2 0.8 0.6], 'EdgeColor', 'none');
xlabel('Max Dynamic Pressure (kPa)'); 
ylabel('Frequency');
title('Distribution of Max q (Dynamic Pressure)');
grid on;

subplot(3,3,6);
histogram(stability_margins, 50, 'FaceColor', [0.2 0.8 0.4], 'EdgeColor', 'none');
hold on;
xline(0, 'r--', 'LineWidth', 2, 'Label', 'Unstable');
xlabel('Stability Margin (CP - CG) [m]'); 
ylabel('Frequency');
title('Static Stability Margin Distribution');
grid on;

subplot(3,3,7);
plot(time_ex, alt_ex/1000, 'b-', 'LineWidth', 2);
ylabel('Altitude (km)');
xlabel('Time (s)');
title('Example Trajectory - Altitude Profile');
grid on;

subplot(3,3,8);
yyaxis left
plot(time_ex, vel_ex, 'r-', 'LineWidth', 2);
ylabel('Velocity (m/s)');
yyaxis right
plot(time_ex, gamma_ex, 'g-', 'LineWidth', 2);
ylabel('Flight Path Angle (deg)');
xlabel('Time (s)');
title('Example Trajectory - Velocity & Flight Path Angle');
legend('Velocity', 'γ', 'Location', 'best');
grid on;

subplot(3,3,9);
plot(q_ex/1000, alt_ex/1000, 'b-', 'LineWidth', 2);
xlabel('Dynamic Pressure q (kPa)');
ylabel('Altitude (km)');
title('Example Trajectory - Dynamic Pressure Profile');
grid on;

fprintf('=== PART (d): Monte Carlo Stability Analysis During Ascent ===\n');
fprintf('Considering fuel consumption and CG movement during flight\n\n');

n_stability_sims = 1000;
n_fuel_steps = 100; 

fprintf('Running %d stability simulations with fuel consumption...\n', n_stability_sims);

unstable_count = 0;
min_margins_all_sims = zeros(n_stability_sims, 1);
time_to_instability = zeros(n_stability_sims, 1);

for sim = 1:n_stability_sims
    dim_varied = dim + (tolerance_dim/3) * randn(size(dim));
    mass_varied = mass + (tolerance_mass/3) * randn(size(mass));
    mass_varied(mass_varied < 0) = 0;
    
    L_total_var = sum(dim_varied);
    L_nose_var = dim_varied(1);
    CP_nose_var = 0.67 * L_nose_var;
    CP_body_var = L_total_var / 2;
    CP_fins_var = L_total_var - 0.15;
    CP_sim = (CN_nose * CP_nose_var + CN_body * CP_body_var + CN_fins * CP_fins_var) / ...
             (CN_nose + CN_body + CN_fins);
    
    fuel_levels = linspace(m_fuel, 0, n_fuel_steps);
    stability_margin_sim = zeros(n_fuel_steps, 1);
    
    for step = 1:n_fuel_steps
        mass_current = mass_varied;
        mass_current(3) = fuel_levels(step); 
        
        total_mass_current = sum(mass_current);
        if total_mass_current > 0
            CG_sim = sum(loc .* mass_current) / total_mass_current;
        else
            CG_sim = 0;
        end
        
        stability_margin_sim(step) = CP_sim - CG_sim;
    end
    
    min_margin = min(stability_margin_sim);
    min_margins_all_sims(sim) = min_margin;
    
    if min_margin < 0
        unstable_count = unstable_count + 1;
        unstable_idx = find(stability_margin_sim < 0, 1);
        time_to_instability(sim) = unstable_idx / n_fuel_steps * burn_time;
    else
        time_to_instability(sim) = NaN; 
    end
end

instability_rate = unstable_count / n_stability_sims * 100;

fprintf('STABILITY ANALYSIS RESULTS:\n');
fprintf('  Simulations with instability: %d/%d (%.2f%%)\n', ...
    unstable_count, n_stability_sims, instability_rate);
fprintf('  Minimum stability margin:\n');
fprintf('    Mean: %.4f m\n', mean(min_margins_all_sims));
fprintf('    Worst case: %.4f m\n', min(min_margins_all_sims));

if unstable_count > 0
    fprintf('\n  ⚠ WARNING: Rocket becomes UNSTABLE in %.2f%% of cases!\n', instability_rate);
    fprintf('  Action required: Add ballast weight at nose\n\n');
    
    desired_margin = 0.5;
    target_CG = CP_nominal - desired_margin;
    masses_empty = mass;
    masses_empty(3) = 0;
    CG_empty = sum(loc .* masses_empty) / sum(masses_empty);
    
    if CG_empty > target_CG
        total_moment = sum(loc .* masses_empty);
        total_mass_empty = sum(masses_empty);
        required_weight = (total_moment - target_CG * total_mass_empty) / target_CG;
        required_weight = max(0, required_weight);
    else
        required_weight = 0;
    end
    
    fprintf('  RECOMMENDATION:\n');
    fprintf('    Add %.2f kg of ballast at the nose (location 0 m)\n', required_weight);
    fprintf('    This will ensure stability margin > 0.5 m in worst case\n');
else
    fprintf('\n  ✓ Rocket remains STABLE in all simulations!\n');
    fprintf('  No additional weight needed.\n');
end

figure('Name', 'Part (d) - Stability During Ascent', 'Position', [100 100 1400 600]);

subplot(1,2,1);
n_plot = 50; 
for sim = 1:n_plot
    dim_varied = dim + (tolerance_dim/3) * randn(size(dim));
    mass_varied = mass + (tolerance_mass/3) * randn(size(mass));
    mass_varied(mass_varied < 0) = 0;
    
    L_total_var = sum(dim_varied);
    L_nose_var = dim_varied(1);
    CP_nose_var = 0.67 * L_nose_var;
    CP_body_var = L_total_var / 2;
    CP_fins_var = L_total_var - 0.15;
    CP_sim = (CN_nose * CP_nose_var + CN_body * CP_body_var + CN_fins * CP_fins_var) / ...
             (CN_nose + CN_body + CN_fins);
    
    fuel_levels = linspace(m_fuel, 0, n_fuel_steps);
    stability_profile = zeros(n_fuel_steps, 1);
    
    for step = 1:n_fuel_steps
        mass_current = mass_varied;
        mass_current(3) = fuel_levels(step);
        total_mass_current = sum(mass_current);
        CG_sim = sum(loc .* mass_current) / total_mass_current;
        stability_profile(step) = CP_sim - CG_sim;
    end
    
    if min(stability_profile) < 0
        plot(fuel_levels, stability_profile, 'r-', 'LineWidth', 0.5, 'Color', [1 0 0 0.3]);
    else
        plot(fuel_levels, stability_profile, 'g-', 'LineWidth', 0.5, 'Color', [0 0.7 0 0.3]);
    end
    hold on;
end
plot([m_fuel 0], [0 0], 'k--', 'LineWidth', 2);
xlabel('Fuel Remaining (kg)');
ylabel('Stability Margin (CP - CG) [m]');
title(sprintf('Stability Margin During Fuel Consumption\n(%d sample trajectories)', n_plot));
grid on;

subplot(1,2,2);
histogram(min_margins_all_sims, 50, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'none');
hold on;
xline(0, 'r--', 'LineWidth', 2.5, 'Label', 'Unstable');
xline(mean(min_margins_all_sims), 'b--', 'LineWidth', 2, 'Label', 'Mean');
xlabel('Minimum Stability Margin (m)');
ylabel('Frequency');
title('Distribution of Minimum Stability Margins');
grid on;
legend('Location', 'best');

fprintf('\n');
fprintf('================================================================\n');
fprintf('                    FINAL SUMMARY REPORT                        \n');
fprintf('================================================================\n\n');

fprintf('PART (c) - Monte Carlo Trajectory Analysis:\n');
fprintf('  • Model: Round, rotating Earth with atmosphere\n');
fprintf('  • Simulations run: %d\n', n_simulations);
fprintf('  • Mean max altitude: %.2f km (target: 35 km)\n', mean(max_altitudes)/1000);
fprintf('  • Success rate: %.2f%%\n', success_rate);
fprintf('  • Mean max velocity: %.2f m/s\n', mean(max_velocities));
fprintf('  • Mean max acceleration: %.2f g\n', mean(max_accelerations)/g0);
fprintf('  • Assumption: Constant CG location\n\n');

fprintf('PART (d) - Stability Analysis with Fuel Consumption:\n');
fprintf('  • Simulations run: %d\n', n_stability_sims);
fprintf('  • Instability rate: %.2f%%\n', instability_rate);
fprintf('  • Worst stability margin: %.4f m\n', min(min_margins_all_sims));

if unstable_count > 0
    fprintf('  • STATUS: ⚠ UNSTABLE - Action Required\n');
    fprintf('  • SOLUTION: Add %.2f kg ballast at nose\n', required_weight);
else
    fprintf('  • STATUS: ✓ STABLE - No action needed\n');
end
fprintf('\n================================================================\n');

