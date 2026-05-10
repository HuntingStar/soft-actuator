%% nonuniform_pneunet_contact_model.m
% 不等高相邻气动腔室接触建模
% 单位体系：
%   长度：mm
%   力：N
%   压力/模量：MPa = N/mm^2
%   弯矩：N*mm
%
% 理论框架：
%   1. 腔室侧壁视为初始平膜；
%   2. 鼓涨轮廓近似为圆弧；
%   3. 膜材料采用不可压缩 neo-Hookean；
%   4. 两个不等高腔室通过自由鼓涨曲线叠加判断初始接触；
%   5. 接触区域近似为椭圆；
%   6. 接触力 Fp = p Ac；
%   7. 接触力方向由高度差诱导的接触面倾角确定；
%   8. 驱动器弯曲用等效弯矩 + 常曲率梁近似展示。

clear; clc; close all;

%% ===================== 1. 参数设置 =====================
par = struct();

% -------- 腔室几何参数 --------
% 原文中膜的高度方向尺寸为 2a
% 这里两个腔室底部齐平，但是高度不同：
% H1 = 2a1, H2 = 2a2
par.a1 = 6.0;          % 高腔室半高，mm
par.a2 = 5.0;          % 矮腔室半高，mm

% 腔室长度方向半尺寸，原文为 b，对应膜尺寸 2b
par.b1 = 8.0;          % 高腔室半长度，mm
par.b2 = 8.0;          % 矮腔室半长度，mm

% 相邻腔室半间隙
% 两个未变形膜之间总间隙为 2hx
par.hx = 0.75;          % mm

% 膜厚
par.tm1 = 1.5;         % mm
par.tm2 = 1.5;         % mm

% 材料剪切模量
% DS20 可先取 mu = 0.136 MPa；DS30 可先取 mu = 0.198 MPa
par.mu1 = 0.198;       % MPa
par.mu2 = 0.198;       % MPa

% -------- actuator cover / 底部梁参数 --------
par.Hcover = 5.0;      % cover 厚度，mm
par.Bcover = 16;     % cover 宽度，mm
par.Lcover = 20.0;     % 驱动器有效长度，mm

% cover 等效弹性模量
% 若近似不可压缩橡胶，E ≈ 3mu
par.Ecover = 3 * 0.136;    % MPa

% 一个气室单元的实体长度参数
% 原文中 hL 表示 air chamber length along X axis
par.hL = 4.2;          % mm

% -------- 绘图压力 --------
p_plot_kPa = 10;       % 单个示意图使用的压力，kPa

% -------- 扫描压力范围 --------
p_sweep_kPa = linspace(1, 70, 160);

%% ===================== 2. 单压力计算 =====================
res = evalOnePressure(par, p_plot_kPa);

disp('================ 单压力计算结果 ================');
fprintf('pressure p = %.2f kPa\n', p_plot_kPa);
fprintf('theta_mf1 = %.4f rad, theta_mf2 = %.4f rad\n', res.theta1, res.theta2);
fprintf('Rf1 = %.4f mm, Rf2 = %.4f mm\n', res.R1, res.R2);
fprintf('Y* = %.4f mm\n', res.Ystar);
fprintf('rho1* = %.4f mm, rho2* = %.4f mm\n', res.rho1, res.rho2);
fprintf('theta_h_init = %.4f deg\n', rad2deg(res.theta_h_init));

if res.hasContact
    fprintf('Contact: YES\n');
    fprintf('Y- = %.4f mm, Y+ = %.4f mm\n', res.Yminus, res.Yplus);
    fprintf('c1 = %.4f mm, c2 = %.4f mm, cy = %.4f mm, cz = %.4f mm\n', ...
        res.c1, res.c2, res.cy, res.cz);
    fprintf('Ac = %.4f mm^2\n', res.Ac);
    fprintf('Fp = %.4f N\n', res.Fp);
    fprintf('Fx = %.4f N, FY = %.4f N\n', res.Fx, res.FY);
    fprintf('theta_h = %.4f deg\n', rad2deg(res.theta_h));
    fprintf('Mm = %.4f N*mm\n', res.Mm);
    fprintf('Mpi = %.4f N*mm\n', res.Mpi);
    fprintf('kappa = %.6f 1/mm\n', res.kappa);
    fprintf('tip angle = %.4f deg\n', rad2deg(res.theta_tip));
else
    fprintf('Contact: NO\n');
end

%% ===================== 3. 压力扫描计算 =====================
sweep = evalPressureSweep(par, p_sweep_kPa);

%% ===================== 4. 绘图 =====================
plotContactGeometry(par, res);
plotVariableSweep(sweep);
plotEquivalentBending(par, res);


%% ========================================================================
%% ========================== 局部函数区 ==================================
%% ========================================================================

function res = evalOnePressure(par, p_kPa)
    % 单个压力点下的全部计算

    p = p_kPa * 1e-3;   % kPa -> MPa = N/mm^2

    res = struct();
    res.p_kPa = p_kPa;
    res.p = p;

    % -------- 1. 求两个腔室的自由鼓涨角 theta_mf --------
    theta1 = solveThetaMF(p, par.a1, par.mu1, par.tm1);
    theta2 = solveThetaMF(p, par.a2, par.mu2, par.tm2);

    res.theta1 = theta1;
    res.theta2 = theta2;

    if isnan(theta1) || isnan(theta2)
        res.valid = false;
        res.hasContact = false;
        fillInvalidFields();
        return;
    end

    res.valid = true;

    % -------- 2. 求自由鼓涨圆弧半径 --------
    R1 = par.a1 / sin(theta1);
    R2 = par.a2 / sin(theta2);

    res.R1 = R1;
    res.R2 = R2;

    % -------- 3. 求不等高腔室的初始接触点 Y* --------
    % Y* = (R2*a1 + R1*a2)/(R1 + R2)
    Ystar = (R2 * par.a1 + R1 * par.a2) / (R1 + R2);

    % 由于矮腔只存在于 [0, 2a2]，接触高度必须在共同高度范围内
    YminDomain = 0;
    YmaxDomain = min(2 * par.a1, 2 * par.a2);

    if Ystar < YminDomain
        Ystar = YminDomain;
    elseif Ystar > YmaxDomain
        Ystar = YmaxDomain;
    end

    res.Ystar = Ystar;
    res.rho1 = Ystar - par.a1;
    res.rho2 = Ystar - par.a2;

    % 初始接触倾角，由高度差诱导
    arg = (par.a1 - par.a2) / (R1 + R2);
    arg = max(min(arg, 1), -1);
    res.theta_h_init = asin(arg);

    % -------- 4. 判断是否接触 --------
    deltaFun = @(Y) uMem(Y, par.a1, R1, theta1) + ...
                    uMem(Y, par.a2, R2, theta2) - 2 * par.hx;

    Phi = deltaFun(Ystar);
    res.Phi = Phi;

    if Phi <= 0
        res.hasContact = false;

        res.Yminus = NaN;
        res.Yplus  = NaN;
        res.Yc     = Ystar;
        res.c1 = 0;
        res.c2 = 0;
        res.cy = 0;
        res.cz = 0;
        res.Ac = 0;
        res.Fp = 0;
        res.Fx = 0;
        res.FY = 0;
        res.theta_h = res.theta_h_init;
        res.Mm = 0;
        res.Mpi = 0;
        res.kappa = 0;
        res.theta_tip = 0;
        [res.beamX, res.beamY] = beamCurve(par.Lcover, 0);

        return;
    end

    res.hasContact = true;

    % -------- 5. 求接触边界 Y-, Y+ --------
    Ygrid = linspace(YminDomain, YmaxDomain, 2000);
    Dgrid = deltaFun(Ygrid);

    [~, imax] = max(Dgrid);

    % 左侧根
    idxL = find(Dgrid(1:imax) <= 0, 1, 'last');
    if isempty(idxL)
        Yminus = YminDomain;
    else
        Yminus = fzero(deltaFun, [Ygrid(idxL), Ygrid(idxL + 1)]);
    end

    % 右侧根
    idxR_rel = find(Dgrid(imax:end) <= 0, 1, 'first');
    if isempty(idxR_rel)
        Yplus = YmaxDomain;
    else
        idxR = imax + idxR_rel - 1;
        if idxR <= imax
            Yplus = YmaxDomain;
        else
            Yplus = fzero(deltaFun, [Ygrid(idxR - 1), Ygrid(idxR)]);
        end
    end

    res.Yminus = Yminus;
    res.Yplus  = Yplus;

    % 接触宽度
    cy = Yplus - Yminus;
    c1 = Ystar - Yminus;
    c2 = Yplus - Ystar;

    c1 = max(c1, 0);
    c2 = max(c2, 0);

    res.c1 = c1;
    res.c2 = c2;
    res.cy = cy;

    % -------- 6. 计算 z 方向接触宽度 cz --------
    % 原文形式：
    % cz = cy * b * (a - rho*) / a^2
    % 不等高推广时，两腔分别给出一个可接触长度，实际取较小值
    cz1 = cy * par.b1 * max(par.a1 - abs(res.rho1), 0) / par.a1^2;
    cz2 = cy * par.b2 * max(par.a2 - abs(res.rho2), 0) / par.a2^2;
    cz = min(cz1, cz2);

    res.cz = cz;

    % -------- 7. 接触面积 --------
    Ac = pi / 4 * cy * cz;
    res.Ac = Ac;

    % -------- 8. 接触区域质心高度 Yc --------
    YcGrid = linspace(Yminus, Yplus, 800);
    Dpos = max(deltaFun(YcGrid), 0);

    if trapz(YcGrid, Dpos) > 1e-12
        Yc = trapz(YcGrid, YcGrid .* Dpos) / trapz(YcGrid, Dpos);
    else
        Yc = Ystar;
    end

    res.Yc = Yc;

    % -------- 9. 有限接触后的等效接触面倾角 --------
    du1 = duMem(Yc, par.a1, R1);
    du2 = duMem(Yc, par.a2, R2);

    sc = 0.5 * (du1 - du2);
    theta_h = atan(sc);

    res.theta_h = theta_h;

    % -------- 10. 接触力 --------
    Fp = p * Ac;      % MPa * mm^2 = N
    Fx = Fp * cos(theta_h);
    FY = -Fp * sin(theta_h);

    res.Fp = Fp;
    res.Fx = Fx;
    res.FY = FY;

    % -------- 11. 接触力产生的弯矩 --------
    % 完整形式：M = |r x F|
    % 近似保留横向偏心 hx 和竖向偏心 Hcover/2 + Yc
    rx = par.hx;
    rY = par.Hcover / 2 + Yc;

    Mm = Fp * (rY * cos(theta_h) + rx * sin(theta_h));
    Mm = max(Mm, 0);

    res.Mm = Mm;

    % -------- 12. 等效驱动弯矩 --------
    Mpi = Mm * sqrt(2 * par.hx / (2 * par.hx + par.hL));
    res.Mpi = Mpi;

    % -------- 13. 等效弯曲曲率 --------
    I = par.Bcover * par.Hcover^3 / 12;
    kappa = Mpi / (par.Ecover * I);

    res.kappa = kappa;
    res.theta_tip = kappa * par.Lcover;

    [res.beamX, res.beamY] = beamCurve(par.Lcover, kappa);

    function fillInvalidFields()
        names = {'R1','R2','Ystar','rho1','rho2','theta_h_init','Phi', ...
                 'Yminus','Yplus','Yc','c1','c2','cy','cz','Ac', ...
                 'Fp','Fx','FY','theta_h','Mm','Mpi','kappa','theta_tip'};
        for k = 1:numel(names)
            res.(names{k}) = NaN;
        end
        res.beamX = NaN;
        res.beamY = NaN;
    end
end


function sweep = evalPressureSweep(par, pList_kPa)
    % 压力扫描

    n = numel(pList_kPa);

    sweep.p_kPa = pList_kPa(:);

    fields = {'theta1','theta2','R1','R2','Ystar','rho1','rho2', ...
              'theta_h_init','Yminus','Yplus','Yc','c1','c2','cy','cz', ...
              'Ac','Fp','Fx','FY','theta_h','Mm','Mpi','kappa','theta_tip','Phi'};

    for f = 1:numel(fields)
        sweep.(fields{f}) = NaN(n, 1);
    end

    sweep.hasContact = false(n, 1);
    sweep.valid = false(n, 1);

    for i = 1:n
        res = evalOnePressure(par, pList_kPa(i));

        sweep.valid(i) = res.valid;
        sweep.hasContact(i) = res.hasContact;

        for f = 1:numel(fields)
            if isfield(res, fields{f})
                sweep.(fields{f})(i) = res.(fields{f});
            end
        end
    end
end


function theta = solveThetaMF(p, a, mu, tm)
    % 求解自由鼓涨角 theta_mf
    %
    % 方程：
    % sin(theta) - sin(theta)^7/theta^6 = p a / (2 mu tm)

    if p <= 0
        theta = 0;
        return;
    end

    rhs = p * a / (2 * mu * tm);

    fun = @(th) sin(th) - sin(th).^7 ./ th.^6 - rhs;

    thGrid = linspace(1e-5, pi - 1e-5, 4000);
    val = fun(thGrid);

    idx = find(val(1:end-1) .* val(2:end) <= 0, 1, 'first');

    if isempty(idx)
        theta = NaN;
        warning('No solution for theta_mf. Try lower pressure or higher stiffness/thickness.');
        return;
    end

    theta = fzero(fun, [thGrid(idx), thGrid(idx + 1)]);
end


function u = uMem(Y, a, R, theta)
    % 自由鼓涨膜轮廓
    %
    % u(Y) = sqrt(R^2 - (Y-a)^2) - R cos(theta)

    inside = R^2 - (Y - a).^2;
    inside = max(inside, 0);

    u = sqrt(inside) - R * cos(theta);
end


function du = duMem(Y, a, R)
    % 自由鼓涨轮廓导数 du/dY

    den = sqrt(max(R^2 - (Y - a).^2, 1e-12));
    du = -(Y - a) ./ den;
end


function [x, y] = beamCurve(L, kappa)
    % 常曲率梁中心线
    %
    % 若 kappa = 0:
    %   x = s, y = 0
    %
    % 若 kappa ~= 0:
    %   x = sin(kappa*s)/kappa
    %   y = (1 - cos(kappa*s))/kappa

    s = linspace(0, L, 400);

    if abs(kappa) < 1e-12
        x = s;
        y = zeros(size(s));
    else
        x = sin(kappa * s) / kappa;
        y = (1 - cos(kappa * s)) / kappa;
    end
end


function plotContactGeometry(par, res)
    % 绘制不等高腔室接触轮廓图

    figure('Color', 'w', 'Name', 'Nonuniform chamber contact geometry');
    hold on; box on; axis equal;

    cyan = [0.0, 0.75, 0.80];
    darkCyan = [0.0, 0.45, 0.50];

    % 左/高腔室
    Y1 = linspace(0, 2 * par.a1, 800);
    u1 = uMem(Y1, par.a1, res.R1, res.theta1);
    xL_free = -par.hx + u1;

    % 右/矮腔室
    Y2 = linspace(0, 2 * par.a2, 800);
    u2 = uMem(Y2, par.a2, res.R2, res.theta2);
    xR_free = par.hx - u2;

    % 未变形侧壁
    plot([-par.hx, -par.hx], [0, 2 * par.a1], 'k:', 'LineWidth', 1.2);
    plot([ par.hx,  par.hx], [0, 2 * par.a2], 'k:', 'LineWidth', 1.2);

    % 自由鼓涨轮廓，虚线
    plot(xL_free, Y1, '--', 'Color', cyan, 'LineWidth', 1.6);
    plot(xR_free, Y2, '--', 'Color', cyan, 'LineWidth', 1.6);

    % 实际接触后的轮廓
    xL_actual = xL_free;
    xR_actual = xR_free;

    if res.hasContact
        idxL = Y1 >= res.Yminus & Y1 <= res.Yplus;
        idxR = Y2 >= res.Yminus & Y2 <= res.Yplus;

        xContactL = contactX(par, res, Y1(idxL));
        xContactR = contactX(par, res, Y2(idxR));

        xL_actual(idxL) = xContactL;
        xR_actual(idxR) = xContactR;
    end

    plot(xL_actual, Y1, '-', 'Color', cyan, 'LineWidth', 3.0);
    plot(xR_actual, Y2, '-', 'Color', cyan, 'LineWidth', 3.0);

    % 接触线
    if res.hasContact
        YcLine = linspace(res.Yminus, res.Yplus, 200);
        xCLine = contactX(par, res, YcLine);
        plot(xCLine, YcLine, 'k-', 'LineWidth', 4.0);

        % 初始接触点 Cp
        xStar = contactX(par, res, res.Ystar);
        plot(xStar, res.Ystar, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        text(xStar + 0.25, res.Ystar + 0.15, '$C_p$', ...
            'Interpreter', 'latex', 'FontSize', 12);

        % 接触质心
        xCent = contactX(par, res, res.Yc);
        plot(xCent, res.Yc, 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
        text(xCent + 0.25, res.Yc, '$Y_c$', ...
            'Interpreter', 'latex', 'FontSize', 12);

        % c1, c2 标注
        xMark = xStar + 0.55;
        plot([xMark, xMark], [res.Yminus, res.Ystar], 'r-', 'LineWidth', 1.5);
        plot([xMark, xMark], [res.Ystar, res.Yplus], 'r-', 'LineWidth', 1.5);
        text(xMark + 0.15, 0.5 * (res.Yminus + res.Ystar), '$c_1$', ...
            'Interpreter', 'latex', 'Color', 'r', 'FontSize', 12);
        text(xMark + 0.15, 0.5 * (res.Ystar + res.Yplus), '$c_2$', ...
            'Interpreter', 'latex', 'Color', 'r', 'FontSize', 12);

        % 接触力方向
        scaleF = 1.8;
        quiver(xCent, res.Yc, ...
            scaleF * cos(res.theta_h), ...
           -scaleF * sin(res.theta_h), ...
            0, 'r', 'LineWidth', 2.2, 'MaxHeadSize', 0.8);

        text(xCent + scaleF * cos(res.theta_h) + 0.15, ...
             res.Yc - scaleF * sin(res.theta_h), ...
             '$\mathbf{F}_p$', ...
             'Interpreter', 'latex', 'Color', 'r', 'FontSize', 13);

        % 接触面倾角
        text(xCent + 0.4, res.Yc - 0.8, ...
            sprintf('$\\theta_h=%.2f^\\circ$', rad2deg(res.theta_h)), ...
            'Interpreter', 'latex', 'FontSize', 12);
    end

    % 腔室中心线
    yline(par.a1, ':', '$a_1$', 'Interpreter', 'latex', 'LineWidth', 1.0);
    yline(par.a2, ':', '$a_2$', 'Interpreter', 'latex', 'LineWidth', 1.0);

    % 底部齐平标注
    plot([-par.hx - 1.2, par.hx + 1.2], [0, 0], 'k-', 'LineWidth', 1.0);
    text(-par.hx - 1.0, -0.45, 'bottom aligned', 'FontSize', 10);

    % 坐标轴箭头
    quiver(0, 2 * par.a1 + 0.6, 0, 0.9, 0, 'k', 'LineWidth', 1.2);
    quiver(0, 2 * par.a1 + 0.6, 0.9, 0, 0, 'k', 'LineWidth', 1.2);
    text(0.08, 2 * par.a1 + 1.5, '$Y_0^\prime$', 'Interpreter', 'latex', 'FontSize', 12);
    text(0.95, 2 * par.a1 + 0.55, '$X_0$', 'Interpreter', 'latex', 'FontSize', 12);

    xlabel('$X$ / mm', 'Interpreter', 'latex');
    ylabel('$Y$ / mm', 'Interpreter', 'latex');

    title(sprintf('Nonuniform chamber contact, p = %.1f kPa', res.p_kPa), ...
        'Interpreter', 'latex');

    xlim([-par.hx - 1.5, par.hx + 3.0]);
    ylim([-0.8, 2 * par.a1 + 1.8]);

    legend({'undeformed wall', 'undeformed wall', ...
            'free profile', 'free profile', ...
            'actual profile', 'actual profile'}, ...
            'Location', 'bestoutside');

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end


function xC = contactX(par, res, Y)
    % 接触线位置，取两个自由轮廓的中面

    u1 = uMem(Y, par.a1, res.R1, res.theta1);
    u2 = uMem(Y, par.a2, res.R2, res.theta2);

    xLeft  = -par.hx + u1;
    xRight =  par.hx - u2;

    xC = 0.5 * (xLeft + xRight);
end


function plotVariableSweep(sweep)
    % 绘制随压力变化的变量

    figure('Color', 'w', 'Name', 'Pressure sweep variables');
    tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    p = sweep.p_kPa;

    % 1. 自由鼓涨角
    nexttile; hold on; box on;
    plot(p, rad2deg(sweep.theta1), 'LineWidth', 1.8);
    plot(p, rad2deg(sweep.theta2), 'LineWidth', 1.8);
    xlabel('p / kPa');
    ylabel('\theta_{mf} / deg');
    title('Free inflation angle');
    legend('\theta_{mf,1}', '\theta_{mf,2}', 'Location', 'best');

    % 2. 初始接触点
    nexttile; hold on; box on;
    plot(p, sweep.Ystar, 'LineWidth', 1.8);
    xlabel('p / kPa');
    ylabel('Y^* / mm');
    title('Initial contact height');

    % 3. 接触倾角
    nexttile; hold on; box on;
    plot(p, rad2deg(sweep.theta_h_init), '--', 'LineWidth', 1.5);
    plot(p, rad2deg(sweep.theta_h), 'LineWidth', 1.8);
    xlabel('p / kPa');
    ylabel('\theta_h / deg');
    title('Contact normal angle');
    legend('initial', 'finite contact', 'Location', 'best');

    % 4. 接触面积
    nexttile; hold on; box on;
    plot(p, sweep.Ac, 'LineWidth', 1.8);
    xlabel('p / kPa');
    ylabel('A_c / mm^2');
    title('Contact area');

    % 5. 接触力
    nexttile; hold on; box on;
    plot(p, sweep.Fp, 'k-', 'LineWidth', 2.0);
    plot(p, sweep.Fx, '--', 'LineWidth', 1.6);
    plot(p, sweep.FY, '--', 'LineWidth', 1.6);
    xlabel('p / kPa');
    ylabel('Force / N');
    title('Contact force components');
    legend('F_p', 'F_x', 'F_Y', 'Location', 'best');

    % 6. 接触宽度
    nexttile; hold on; box on;
    plot(p, sweep.c1, 'LineWidth', 1.6);
    plot(p, sweep.c2, 'LineWidth', 1.6);
    plot(p, sweep.cy, 'k-', 'LineWidth', 1.8);
    plot(p, sweep.cz, '--', 'LineWidth', 1.6);
    xlabel('p / kPa');
    ylabel('width / mm');
    title('Contact widths');
    legend('c_1', 'c_2', 'c_y', 'c_z', 'Location', 'best');

    % 7. 弯矩
    nexttile; hold on; box on;
    plot(p, sweep.Mm, 'LineWidth', 1.8);
    plot(p, sweep.Mpi, '--', 'LineWidth', 1.8);
    xlabel('p / kPa');
    ylabel('Moment / N mm');
    title('Moment');
    legend('M_m', 'M_{pi}', 'Location', 'best');

    % 8. 曲率
    nexttile; hold on; box on;
    plot(p, sweep.kappa, 'LineWidth', 1.8);
    xlabel('p / kPa');
    ylabel('\kappa / mm^{-1}');
    title('Equivalent curvature');

    % 9. 末端角
    nexttile; hold on; box on;
    plot(p, rad2deg(sweep.theta_tip), 'LineWidth', 1.8);
    xlabel('p / kPa');
    ylabel('\theta_{tip} / deg');
    title('Equivalent tip angle');

    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 10);
end


function plotEquivalentBending(par, res)
    % 绘制等效弯曲形状

    figure('Color', 'w', 'Name', 'Equivalent bending of actuator');
    hold on; box on; axis equal;

    % 未变形中心线
    plot([0, par.Lcover], [0, 0], 'k--', 'LineWidth', 1.2);

    % 变形中心线
    plot(res.beamX, res.beamY, 'b-', 'LineWidth', 3.0);

    % 在末端画切线方向
    xTip = res.beamX(end);
    yTip = res.beamY(end);
    quiver(xTip, yTip, ...
        8 * cos(res.theta_tip), ...
        8 * sin(res.theta_tip), ...
        0, 'r', 'LineWidth', 2.0, 'MaxHeadSize', 0.7);

    % 等效弯矩标注
    text(0.05 * par.Lcover, max(res.beamY) + 3, ...
        sprintf('$M_{pi}=%.3f\\ \\mathrm{N\\,mm}$', res.Mpi), ...
        'Interpreter', 'latex', 'FontSize', 12);

    text(0.05 * par.Lcover, max(res.beamY) + 0.8, ...
        sprintf('$\\kappa=%.5f\\ \\mathrm{mm^{-1}}$', res.kappa), ...
        'Interpreter', 'latex', 'FontSize', 12);

    text(0.05 * par.Lcover, max(res.beamY) - 1.4, ...
        sprintf('$\\theta_{tip}=%.2f^\\circ$', rad2deg(res.theta_tip)), ...
        'Interpreter', 'latex', 'FontSize', 12);

    xlabel('$X$ / mm', 'Interpreter', 'latex');
    ylabel('$Y$ / mm', 'Interpreter', 'latex');
    title('Equivalent bending caused by chamber contact', 'Interpreter', 'latex');

    legend({'undeformed centerline', 'deformed centerline', 'tip tangent'}, ...
        'Location', 'best');

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
end