% 麦克纳姆轮8字轨迹跟踪（模型预测控制）
clear; clc; close all;

%% ====================== 1. 核心参数优化 =======================
% 麦克纳姆轮物理参数
r = 0.05;           % 车轮半径 (m)
lx = 0.2;           % 车身x方向半长
ly = 0.2;           % 车身y方向半宽
umax = 20;          % 最大转速
dt = 0.1;           % 采样时间 (s)
tau = 0.1;          % 响应时间常数

% 8字轨迹参数
omega_ref = 0.2; % 8 字轨迹的 "角频率"（控制轨迹周期）        
n_cycles = 2;           
T_cycle = 2*pi/omega_ref;  %1 个 8 字轨迹的 "周期"（转 1 圈的时间）
T_total = n_cycles * T_cycle;  
steps = round(T_total/dt);     
time = 0:dt:T_total;           

if length(time) ~= steps + 1
    time = linspace(0, T_total, steps + 1);
end

% 麦克纳姆轮运动学矩阵
compensation = 1.05;  
B_kinematic = compensation * r/4 * [
    1,  1,  1,  1;          % vx
   -1,  1,  1, -1;          % vy
    1/(lx+ly), -1/(lx+ly), 1/(lx+ly), -1/(lx+ly)  % omega
];

%% ====================== 2. 状态空间模型 =======================
A = [
    1 0 0 dt 0 0 0.5*dt^2 0         0;
    0 1 0 0  dt 0 0         0.5*dt^2 0;
    0 0 1 0  0  dt 0         0         0.5*dt^2;
    0 0 0 1  0  0  dt        0         0;
    0 0 0 0  1  0  0         dt        0;
    0 0 0 0  0  1  0         0         dt;
    0 0 0 -dt/tau 0  0  1   -dt/tau    0;
    0 0 0 0  -dt/tau 0 0     1       -dt/tau;
    0 0 0 0  0  -dt/tau 0    0         1
];

B = zeros(9,4); 
B(7:9,:) = (1/tau) * B_kinematic;

C = eye(9); 
D = zeros(9,4);
sys = ss(A,B,C,D,dt);

%% ====================== 3. MPC控制器优化 =======================
prediction_horizon = 25;  
control_horizon = 12;     
mpc_obj = mpc(sys, dt, prediction_horizon, control_horizon);

mpc_obj.Model.Disturbance = [];
mpc_obj.Model.Noise = [];

mpc_obj.Weights.OutputVariables = [150 150 30 30 30 30 20 20 20];
mpc_obj.Weights.ManipulatedVariables = [0.2 0.2 0.2 0.2];
mpc_obj.Weights.ManipulatedVariablesRate = [1 1 1 1];

for i = 1:4
    mpc_obj.MV(i).Min = -umax;
    mpc_obj.MV(i).Max = umax;
    mpc_obj.MV(i).RateMin = -50;
    mpc_obj.MV(i).RateMax = 50;
end

setEstimator(mpc_obj, 'custom');%不使用 MPC 默认的状态估计器（卡尔曼滤波器），而是后续代码中定义状态的更新逻辑
mpc_state = mpcstate(mpc_obj);%创建一个与 MPC 控制器 mpc_obj 对象匹配的 "状态对象"，用于存储和跟踪系统的当前状态

%% ====================== 4. 仿真初始化 ======================
initial_state = zeros(9,1);
mpc_state.Plant = initial_state;

% 预计算参考轨迹缓存（所有时间步的参考轨迹已提前算好）
ref_cache = zeros(9, steps+1);
for k = 1:steps+1
    t = (k-1)*dt;
    ref_cache(:,k) = ref_generator(t);
end

ref0 = ref_cache(:,1);
ref_hist = zeros(9, steps+1);
ref_hist(:,1) = ref0;

state_hist = zeros(9, steps+1);
w_hist = zeros(4, steps);
pos_errors = zeros(1, steps+1);
state_hist(:,1) = initial_state;
pos_errors(1) = norm(initial_state(1:2) - ref0(1:2));

%% ====================== 5. 动画初始化（固定Y轴范围至±2）======================
figure('Position', [100, 100, 1000, 700]); hold on; grid on; axis equal;
xlabel('X 位置 (m)', 'FontSize', 11); ylabel('Y 位置 (m)', 'FontSize', 11);
title('麦克纳姆轮8字轨迹跟踪（MPC）', 'FontSize', 12);

% ========== 关键修改1：参考轨迹直接用预计算的完整ref_cache绘制，而非仅初始点 ==========
% h_ref：完整参考轨迹（一开始就全画出来，红色虚线）
h_ref = plot(ref_cache(1,:), ref_cache(2,:), 'r--', 'LineWidth', 2.5, 'DisplayName', '参考轨迹');
% h_real：实际轨迹（随仿真逐步更新，蓝色实线）
h_real = plot(state_hist(1,1), state_hist(2,1), 'b-', 'LineWidth', 2, 'DisplayName', '实际轨迹');
% h_pred：预测轨迹（随仿真逐步更新，青色点线）
h_pred = plot(nan, nan, 'c.-', 'LineWidth', 1.5, 'MarkerSize', 10, 'DisplayName', '预测轨迹');

% 固定Y轴范围为±2，X轴范围根据轨迹自动调整并保持比例
y_min = -2;    % 固定Y轴最小值
y_max = 2;     % 固定Y轴最大值

% 计算X轴范围以保持坐标轴比例
x_vals = ref_cache(1,:);
x_range = max(x_vals) - min(x_vals);
y_range = y_max - y_min;

% 确保X轴范围足够显示所有X数据并保持与Y轴的比例
x_margin = 0.1;  % 10%的余量
x_center = (max(x_vals) + min(x_vals)) / 2;
required_x_half_range = max(x_range/2 * (1 + x_margin), y_range/2);

x_min = x_center - required_x_half_range;
x_max = x_center + required_x_half_range;

axis([x_min, x_max, y_min, y_max]);

% 调整文本位置
h_info = text(x_min + 0.2, y_max - 0.3, sprintf('步数: 0/%d  误差: %.4fm', steps, pos_errors(1)), ...
             'BackgroundColor', 'w', 'FontSize', 10);

legend('Location', 'northwest', 'FontSize', 10);

%% ====================== 6. 主循环 =======================
animation_speed = 0.05;

for k = 1:steps
    t = time(k);
    current_state = mpc_state.Plant;

    ref = ref_cache(:,k);
    ref_hist(:,k) = ref;

    noise = 0.005;
    measured_acc = ref(7:9) + noise*randn(3,1);
    current_state(7:9) = 0.9 * current_state(7:9) + 0.1 * measured_acc;%模拟加速度计数据融合

    Yref = zeros(prediction_horizon, 9);
    for j = 1:prediction_horizon
        idx = k + j;
        if idx <= steps+1
            Yref(j,:) = ref_cache(:,idx)';
        else
            t_fut = t + j*dt;
            Yref(j,:) = ref_generator(t_fut)';
        end
    end
    [u_opt, info] = mpcmove(mpc_obj, mpc_state, current_state, Yref);

    w_hist(:,k) = u_opt;
    state_hist(:,k+1) = current_state;
    pos_errors(1,k+1) = norm(current_state(1:2) - ref(1:2));

    next_state = A * current_state + B * u_opt;
    mpc_state.Plant = next_state;

    % ========== 关键修改2：删除参考轨迹的更新代码（因为一开始已画完整） ==========
    % 仅更新实际轨迹和预测轨迹
    set(h_real, 'XData', state_hist(1,1:k+1), 'YData', state_hist(2,1:k+1));
    
    pred_traj = zeros(2, prediction_horizon);
    for j = 1:prediction_horizon
        pred_traj(:,j) = info.Xopt(j,1:2)';
    end
    set(h_pred, 'XData', pred_traj(1,:), 'YData', pred_traj(2,:));
    
    set(h_info, 'String', sprintf('步数: %d/%d  时间: %.1fs  误差: %.4fm', ...
                                 k, steps, t, pos_errors(k+1)));
    
    drawnow;
    pause(animation_speed);
end

%% ====================== 7. 仿真结束处理 ======================
ref_hist(:,end) = ref_cache(:,end);
state_hist(:,end) = mpc_state.Plant;

if length(pos_errors) ~= length(time)
    pos_errors = pos_errors(1:length(time));
end

%% ====================== 8. 结果分析 =======================
fprintf('\n8字轨迹%d圈跟踪结果：\n', n_cycles);
fprintf('最终位置误差: %.4f m\n', pos_errors(end));
fprintf('平均位置误差: %.4f m\n', mean(pos_errors));
fprintf('均方根误差: %.4f m\n', rms(pos_errors));

% 车轮转速图
figure('Position', [200, 200, 1000, 400]);
plot(time(1:steps), w_hist', 'LineWidth', 1.5);
title('车轮转速变化', 'FontSize', 12);
xlabel('时间 (s)', 'FontSize', 11); ylabel('转速 (rad/s)', 'FontSize', 11);
ylim([-umax-2, umax+2]);
yline(umax, 'r--', '最大转速', 'FontSize', 10);
yline(-umax, 'r--', '最小转速', 'FontSize', 10);
legend('车轮1','车轮2','车轮3','车轮4', 'Location', 'best', 'FontSize', 10);
grid on;

% 位置误差图
figure('Position', [200, 300, 1000, 400]);
plot(time, pos_errors, 'LineWidth', 1.5, 'Color', [0.8 0 0]);
title('位置跟踪误差', 'FontSize', 12);
xlabel('时间 (s)', 'FontSize', 11); ylabel('误差 (m)', 'FontSize', 11);
ylim([0, max(pos_errors)*1.2]);
grid on;

%% ====================== 9. 参考点函数 =======================
function ref = ref_generator(t)
    a = 3; omega = 0.2; dt = 0.1;

    denom = 1 + cos(omega*t)^2;
    x = a*sin(omega*t)/denom;
    y = a*sin(omega*t)*cos(omega*t)/denom;

    dx = a*omega*(cos(omega*t)*denom + 2*sin(omega*t)^2*cos(omega*t))/denom^2;
    dy = a*omega*(cos(2*omega*t)*denom + sin(omega*t)^2*cos(omega*t)^2)/denom^2;
    theta = atan2(dy, dx);

    speed = 0.3 + 0.1*sin(0.3*t);
    vx = speed * cos(theta);
    vy = speed * sin(theta);

    ddx_term1 = sin(omega*t)*denom*(1+3*cos(omega*t)^2);
    ddx_term2 = 4*cos(omega*t)^2*sin(omega*t)*denom;
    ddx = -a*omega^2*(ddx_term1 + ddx_term2)/denom^3;

    ddy_term1 = sin(2*omega*t)*denom*(3+cos(2*omega*t));
    ddy_term2 = 2*cos(omega*t)*sin(omega*t)^3;
    ddy = -a*omega^2*(ddy_term1 + ddy_term2)/denom^3;

    curvature = (dx*ddy - dy*ddx)/max((dx^2 + dy^2), 1e-6)^(3/2);
    omega_ref = speed * curvature;

    persistent prev_vx prev_vy prev_omega last_t
    if isempty(prev_vx)
        prev_vx = vx; prev_vy = vy; prev_omega = omega_ref; last_t = t;
    end
    ax = (vx - prev_vx) / (t - last_t + 1e-5);
    ay = (vy - prev_vy) / (t - last_t + 1e-5);
    alpha = (omega_ref - prev_omega) / (t - last_t + 1e-5);

    prev_vx = vx; prev_vy = vy; prev_omega = omega_ref; last_t = t;

    ref = [x; y; theta; vx; vy; omega_ref; ax; ay; alpha];
end