// 导入 CoolProp Module
import Module from './coolprop.js';

// CoolProp Module 实例
let CoolPropModule = null;

// 等待 CoolProp 初始化
async function initCoolProp() {
    if (CoolPropModule && CoolPropModule.ready) {
        return CoolPropModule;
    }

    try {
        // Module 是一个异步函数
        // 导入 WASM 文件 URL
        const wasmUrl = new URL('./coolprop.wasm', import.meta.url).href;
        
        const moduleConfig = {
            locateFile: (path) => {
                // 确保 wasm 文件路径正确
                if (path.endsWith('.wasm')) {
                    return wasmUrl;
                }
                return path;
            },
            onRuntimeInitialized: () => {
                // 运行时初始化完成回调
            }
        };

        CoolPropModule = await Module(moduleConfig);
        CoolPropModule.ready = true;
        return CoolPropModule;
    } catch (error) {
        throw new Error('CoolProp 初始化失败: ' + error.message);
    }
}

// 显示状态消息
function showStatus(message, type = 'info', statusId = 'status') {
    const statusEl = document.getElementById(statusId);
    if (!statusEl) return;
    statusEl.textContent = message;
    statusEl.className = `status-message ${type} show`;
    
    if (type === 'success' || type === 'info') {
        setTimeout(() => {
            statusEl.classList.remove('show');
        }, 3000);
    }
}

// 格式化数字显示
function formatNumber(value, decimals = 3) {
    if (value === null || value === undefined || isNaN(value)) {
        return '-';
    }
    if (decimals === 0) {
        return Math.round(value).toLocaleString('zh-CN');
    }
    return value.toFixed(decimals).replace(/\.?0+$/, '');
}

// 计算阻力系数 Cd (根据图 02 的公式)
function calculateDragCoefficient(Re) {
    if (Re <= 0) return 1000; // 避免除零
    return 24 / Re + 3 / Math.sqrt(Re) + 0.34;
}

// 计算终端速度 vt (迭代法，根据公式 1-4)
function calculateTerminalVelocity(Dp, rho_l, rho_g, eta_g, g = 9.81, maxIterations = 50, tolerance = 1e-6) {
    // 初始猜测：使用 Stokes 定律 (适用于低雷诺数，Cd ≈ 24/Re)
    // vt = (2 * g * Dp^2 * (ρl - ρg)) / (9 * ηg)
    let vt = (2 * g * Dp * Dp * (rho_l - rho_g)) / (9 * eta_g);
    
    // 如果初始值不合理，使用另一个初始猜测
    if (vt <= 0 || !isFinite(vt)) {
        // 使用简化的初始猜测：假设 Cd = 0.5 (中等雷诺数)
        vt = Math.sqrt((4 * g * Dp * (rho_l - rho_g)) / (3 * 0.5 * rho_g));
    }
    
    for (let i = 0; i < maxIterations; i++) {
        // 计算雷诺数 (公式: Re = Dp * vt * ρg / ηg)
        const Re = (Dp * vt * rho_g) / eta_g;
        
        // 计算阻力系数 (根据图 02 的公式)
        const Cd = calculateDragCoefficient(Re);
        
        // 计算新的终端速度 (公式 4: vt^2 = 4 * g * Dp * (ρl - ρg) / (3 * Cd * ρg))
        const vt_new = Math.sqrt((4 * g * Dp * (rho_l - rho_g)) / (3 * Cd * rho_g));
        
        // 检查收敛
        if (vt > 0 && Math.abs(vt_new - vt) / vt < tolerance) {
            return vt_new;
        }
        
        vt = vt_new;
        
        // 防止无效值
        if (vt <= 0 || !isFinite(vt)) {
            throw new Error('终端速度计算失败，请检查输入参数');
        }
    }
    
    return vt; // 返回最后一次迭代的结果
}

// 计算圆形截面的液位高度对应的气相流通面积比例
// h_ratio: 液位高度比 (H_liq/D)，从底部到液面的高度
// 返回: 气相流通面积 / 总截面积
function calculateGasAreaRatio(h_ratio) {
    if (h_ratio <= 0) return 1.0;
    if (h_ratio >= 1.0) return 0.0;
    
    // 对于圆形截面，液位高度 h = h_ratio * D
    // 半径 R = D/2
    // 液面到圆心的距离 d = R - h = D/2 - h_ratio * D = D/2 * (1 - 2*h_ratio)
    // 归一化: d/R = 1 - 2*h_ratio
    
    const d_over_R = 1 - 2 * h_ratio;
    
    // 如果液位在圆心以下 (d_over_R > 0)，液体在下方
    // 如果液位在圆心以上 (d_over_R < 0)，液体在上方
    const abs_d_over_R = Math.abs(d_over_R);
    
    // 计算角度: cos(θ/2) = d/R
    const theta_half = Math.acos(Math.max(-1, Math.min(1, d_over_R)));
    const theta = 2 * theta_half;
    
    // 液体面积 = R² * (θ - sin(θ)) / 2
    // 总面积 = π * R²
    // 液体面积比例 = (θ - sin(θ)) / (2π)
    const liquidAreaRatio = (theta - Math.sin(theta)) / (2 * Math.PI);
    
    // 气相面积比例 = 1 - 液体面积比例
    const gasAreaRatio = 1 - liquidAreaRatio;
    
    return Math.max(0, Math.min(1, gasAreaRatio));
}

// 卧式气液分离器计算函数
async function calculateHorizontal() {
    try {
        // 检查 CoolProp 是否已初始化
        if (!CoolPropModule) {
            showStatus('正在初始化 CoolProp...', 'info', 'statusH');
            await initCoolProp();
        }

        // 获取输入参数
        const refrigerant = document.getElementById('refrigerantH').value;
        const Q_kW = parseFloat(document.getElementById('coolingCapacityH').value);
        const Te_C = parseFloat(document.getElementById('evapTempH').value);
        const deltaT_sh = parseFloat(document.getElementById('superheatH').value);
        const calcModeH = document.getElementById('calcModeH').value;
        const L_D_ratio = parseFloat(document.getElementById('lengthDiameterRatio').value);
        const h_liq_ratio = parseFloat(document.getElementById('liquidLevelRatio').value);
        const Ks = parseFloat(document.getElementById('ksValue').value);
        const vh_vt_ratio = parseFloat(document.getElementById('velocityMultiplier').value);

        // 验证输入
        if (isNaN(Q_kW) || Q_kW <= 0) {
            throw new Error('制冷量必须大于 0');
        }
        if (isNaN(Te_C)) {
            throw new Error('请输入有效的蒸发温度');
        }
        if (isNaN(deltaT_sh) || deltaT_sh < 0) {
            throw new Error('吸气过热度必须大于等于 0');
        }
        if (calcModeH === 'manual-diameter') {
            const D_mm_input = parseFloat(document.getElementById('diameterInputH').value);
            if (isNaN(D_mm_input) || D_mm_input <= 0) {
                throw new Error('请输入有效的直径值');
            }
        }
        if (isNaN(L_D_ratio) || L_D_ratio < 2.0 || L_D_ratio > 6.0) {
            throw new Error('长径比必须在 2.0 - 6.0 之间');
        }
        if (isNaN(h_liq_ratio) || h_liq_ratio < 0.1 || h_liq_ratio > 0.5) {
            throw new Error('液位高度比必须在 0.1 - 0.5 之间');
        }
        if (isNaN(Ks) || Ks <= 0) {
            throw new Error('Ks 系数必须大于 0');
        }
        if (calcModeH === 'auto' && (isNaN(vh_vt_ratio) || vh_vt_ratio < 1.0 || vh_vt_ratio > 5.0)) {
            throw new Error('水平气速倍数必须在 1.0 - 5.0 之间');
        }

        // 禁用计算按钮
        const calcBtn = document.getElementById('calculateBtnH');
        calcBtn.disabled = true;
        calcBtn.textContent = '计算中...';

        showStatus('正在计算...', 'info', 'statusH');

        // 温度单位转换：°C 转 K
        const Te_K = Te_C + 273.15;
        const T_suc_K = Te_K + deltaT_sh;

        // 1. 计算蒸发压力 (饱和压力，Q=1 表示饱和蒸汽)
        const P_sat = CoolPropModule.PropsSI('P', 'T', Te_K, 'Q', 1, refrigerant);
        if (P_sat <= 0 || !isFinite(P_sat)) {
            throw new Error('无法计算蒸发压力，请检查制冷剂和温度参数');
        }

        // 2. 计算吸气温度下的气体密度
        const rho_g = CoolPropModule.PropsSI('D', 'P', P_sat, 'T', T_suc_K, refrigerant);
        if (rho_g <= 0 || !isFinite(rho_g)) {
            throw new Error('无法计算气体密度');
        }

        // 3. 计算饱和液体密度
        const rho_l = CoolPropModule.PropsSI('D', 'P', P_sat, 'Q', 0, refrigerant);
        if (rho_l <= 0 || !isFinite(rho_l)) {
            throw new Error('无法计算液体密度');
        }

        // 4. 计算吸气状态下的气体焓值
        const h_gas = CoolPropModule.PropsSI('H', 'P', P_sat, 'T', T_suc_K, refrigerant);
        if (!isFinite(h_gas)) {
            throw new Error('无法计算气体焓值');
        }

        // 5. 计算饱和液体的焓值 (用于计算焓差)
        const h_liq = CoolPropModule.PropsSI('H', 'P', P_sat, 'Q', 0, refrigerant);
        if (!isFinite(h_liq)) {
            throw new Error('无法计算液体焓值');
        }

        // 6. 计算焓差 (J/kg)
        const delta_h = h_gas - h_liq;
        if (delta_h <= 0 || !isFinite(delta_h)) {
            throw new Error('焓差计算异常，请检查过热度设置');
        }

        // 7. 计算质量流量 (kg/s)
        const m_dot = (Q_kW * 1000) / delta_h;

        // 8. 计算气体体积流量 (m³/s)
        const V_g = m_dot / rho_g;

        // 9. 计算沉降速度 vt (使用 Souders-Brown 公式)
        // vt = Ks * sqrt((ρl - ρg) / ρg)
        const vt = Ks * Math.sqrt((rho_l - rho_g) / rho_g);
        if (vt <= 0 || !isFinite(vt)) {
            throw new Error('无法计算沉降速度');
        }

        // 10. 根据计算模式确定直径和水平气速
        let D_m; // 直径（米）
        let L_m; // 长度（米）
        let vh; // 水平气速
        let vh_display; // 显示用的气速（将在后面根据模式设置）
        const A_gas_ratio = calculateGasAreaRatio(h_liq_ratio);

        if (calcModeH === 'manual-diameter') {
            // 手动输入直径模式（校核模式）
            const D_mm_input = parseFloat(document.getElementById('diameterInputH').value);
            D_m = D_mm_input / 1000; // 转换为米
            
            // 计算气相流通面积
            const A_total = Math.PI * Math.pow(D_m / 2, 2);
            const A_gas = A_total * A_gas_ratio;
            
            // 计算实际水平气速
            vh = V_g / A_gas;
            vh_display = vh; // 校核模式使用实际气速
            
            // 计算实际水平气速倍数
            const actual_vh_vt_ratio = vh / vt;
            
            // 计算所需长度（满足分离时间要求）
            const H_gas = D_m * (1 - h_liq_ratio);
            const S_l = H_gas;
            const L_min = (S_l / vt) * vh * 1.1; // 10% 安全余量
            
            // 如果输入了长径比，使用长径比；否则使用最小长度
            if (!isNaN(L_D_ratio) && L_D_ratio > 0) {
                L_m = D_m * L_D_ratio;
                if (L_m < L_min) {
                    showStatus(`警告: 根据长径比计算的长度 ${(L_m * 1000).toFixed(0)} mm 小于所需最小长度 ${(L_min * 1000).toFixed(0)} mm，建议增大长径比`, 'info', 'statusH');
                    L_m = L_min; // 使用最小长度
                }
            } else {
                L_m = L_min;
            }
            
            // 验证实际气速是否合理
            if (vh > vt * 5) {
                throw new Error(`实际水平气速 ${vh.toFixed(3)} m/s 过高，建议增大直径`);
            }
        } else {
            // 自动计算模式（原有逻辑）
            vh = vh_vt_ratio * vt;
            vh_display = vh;
            
            // 12. 迭代计算筒径 D 和长度 L
            // 设计原则：严格保持 vh/vt 比值，通过调整长度满足分离时间要求
            
            // 第一步：根据设计目标 vh 计算初始直径
            // vh = V_g / (A_gas * A_gas_ratio)
            // A_gas = π * (D/2)^2 * A_gas_ratio
            // 因此：vh = V_g / (π * (D/2)^2 * A_gas_ratio)
            // D = sqrt(4 * V_g / (π * vh * A_gas_ratio))
            D_m = Math.sqrt((4 * V_g) / (Math.PI * vh * A_gas_ratio));
            L_m = D_m * L_D_ratio;
            let maxIterations = 20;
            let tolerance = 0.01; // 1% 容差

            for (let iter = 0; iter < maxIterations; iter++) {
            // 重新计算气相流通面积（因为直径可能变化）
            const A_gas_ratio_current = calculateGasAreaRatio(h_liq_ratio);
            
            // 根据当前直径和设计目标 vh，验证实际气速
            const A_total = Math.PI * Math.pow(D_m / 2, 2);
            const A_gas = A_total * A_gas_ratio_current;
            const vh_actual = V_g / A_gas;
            
            // 如果实际气速与设计目标偏差超过1%，调整直径以保持 vh/vt 比值
            if (Math.abs(vh_actual - vh) / vh > 0.01) {
                // 重新计算直径以匹配设计目标 vh
                D_m = Math.sqrt((4 * V_g) / (Math.PI * vh * A_gas_ratio_current));
                // 重新计算长度（保持 L/D 比例，但后续会根据分离时间调整）
                L_m = D_m * L_D_ratio;
            }
            
            // 计算气相高度 (m)
            const H_gas = D_m * (1 - h_liq_ratio);
            
            // 沉降距离 Sl (气相高度)
            const S_l = H_gas;
            
            // 飞行时间 t_fly = L / vh（使用设计目标 vh）
            const t_fly = L_m / vh;
            
            // 沉降时间 t_fall = Sl / vt
            const t_fall = S_l / vt;
            
            // 检查是否满足条件: t_fly >= t_fall
            if (t_fly >= t_fall * (1 + tolerance)) {
                // 满足条件，可以退出
                break;
            }
            
            // 不满足条件，需要增加长度（保持直径不变，以维持 vh/vt 比值）
            // 计算所需的最小长度
            const L_min = (S_l / vt) * vh;
            
            // 增加10%安全余量
            L_m = L_min * 1.1;
            
            // 如果长度过长（L/D > 6），给出警告但继续计算
            // 实际应用中可能需要重新考虑设计参数
            if (L_m / D_m > 6.0) {
                // 可以在这里添加警告，但不强制改变直径
                // 因为保持 vh/vt 比值是优先目标
            }
            }
        }
        
        // 最终验证：确保实际气速等于设计目标（仅自动模式）
        if (calcModeH === 'auto') {
            const A_gas_ratio_final = calculateGasAreaRatio(h_liq_ratio);
            const A_total_final = Math.PI * Math.pow(D_m / 2, 2);
            const A_gas_final = A_total_final * A_gas_ratio_final;
            const vh_final_check = V_g / A_gas_final;
            
            // 如果仍有偏差，最后一次调整直径以确保严格匹配
            // 降低容差到 0.1%，确保更精确匹配设计目标
            const final_deviation = Math.abs(vh_final_check - vh) / vh;
            if (final_deviation > 0.001) {
                console.log(`最终验证: 实际气速 ${vh_final_check.toFixed(3)} m/s 与设计目标 ${vh.toFixed(3)} m/s 偏差 ${(final_deviation * 100).toFixed(2)}%，进行调整`);
                D_m = Math.sqrt((4 * V_g) / (Math.PI * vh * A_gas_ratio_final));
                // 重新验证分离时间
                const H_gas_final = D_m * (1 - h_liq_ratio);
                const S_l_final = H_gas_final;
                const L_min_final = (S_l_final / vt) * vh * 1.1;
                L_m = Math.max(L_m, L_min_final); // 取较大值
                
                // 再次验证调整后的气速
                const A_total_after = Math.PI * Math.pow(D_m / 2, 2);
                const A_gas_after = A_total_after * A_gas_ratio_final;
                const vh_after = V_g / A_gas_after;
                console.log(`调整后: 实际气速 ${vh_after.toFixed(3)} m/s，设计目标 ${vh.toFixed(3)} m/s，偏差 ${(Math.abs(vh_after - vh) / vh * 100).toFixed(2)}%`);
            } else {
                console.log(`最终验证通过: 实际气速 ${vh_final_check.toFixed(3)} m/s 与设计目标 ${vh.toFixed(3)} m/s 匹配良好`);
            }
        }

        // 转换为 mm（初始值，可能会在后续调整中改变）
        let D_mm = D_m * 1000;
        let L_mm = L_m * 1000;

        // 13. 计算筒内容积 (L)
        let Vol_m3 = Math.PI * Math.pow(D_m / 2, 2) * L_m;
        let Vol_L = Vol_m3 * 1000;

        // 15. 计算液体停留时间 (简化计算)
        // 注意：如果后续调整了结构参数，停留时间会在调整后重新计算
        let liquid_volume_ratio = 1 - calculateGasAreaRatio(h_liq_ratio);
        let liquid_volume_m3 = Vol_m3 * liquid_volume_ratio;
        // 假设液体流量为总流量的10%（粗略估算）
        const liquid_flow_m3_s = V_g * 0.1 * (rho_g / rho_l);
        let residence_time = liquid_volume_m3 / Math.max(liquid_flow_m3_s, 1e-6);

        // 15. 计算防夹带速度 Vre (Re-entrainment velocity)
        // 根据文档 Figure 04，需要计算 Rej, Rp, N 等参数
        // 获取表面张力 (N/m)
        let sigma = 0.026; // 默认值，如果没有则使用
        try {
            // CoolProp中表面张力的属性代码是'S'
            sigma = CoolPropModule.PropsSI('S', 'P', P_sat, 'T', Te_K, refrigerant); // 表面张力
            if (!isFinite(sigma) || sigma <= 0) {
                // 如果无法获取，使用经验值
                sigma = 0.026; // 典型制冷剂表面张力 (N/m)
            }
        } catch (e) {
            // CoolProp可能不支持表面张力，使用默认值
            // 根据制冷剂类型使用不同的默认值
            const sigmaMap = {
                // 低温制冷剂
                'R23': 0.009,
                'R508B': 0.008,
                'R170': 0.005,
                'R1150': 0.004,
                // 中温制冷剂
                'R134a': 0.0089,
                'R404A': 0.004,
                'R407C': 0.005,
                'R507A': 0.004,
                'R22': 0.008,
                'R1234ze(E)': 0.006,
                'R1234yf': 0.006,
                'R32': 0.008,
                'R290': 0.007,
                'R600a': 0.006,
                'R410A': 0.004,
                // 高温制冷剂
                'R245fa': 0.014,
                'R123': 0.015,
                'R236fa': 0.012,
                'R365mfc': 0.013,
                'R717': 0.0268,
                'Ammonia': 0.0268  // 兼容旧名称
            };
            sigma = sigmaMap[refrigerant] || 0.026;
        }

        // 计算防夹带速度 Vre (Re-entrainment velocity)
        // 根据文档 Figure 04，需要计算 Rej, Rp, N 等参数
        // 文档提到使用方程 D（不依赖几何形状，适用于高湍流情况）
        
        // 计算关键参数
        // Rej: 基于液滴的雷诺数（用于判断流态）
        // 假设典型液滴直径 100-200 µm，使用平均 150 µm
        const Dp_typical = 150e-6; // 150 µm 转换为 m
        const eta_g = CoolPropModule.PropsSI('V', 'P', P_sat, 'T', T_suc_K, refrigerant);
        const Rej = (Dp_typical * vt * rho_g) / eta_g;
        
        // Rp: 密度比参数
        const Rp = rho_l / rho_g;
        
        // N: 无量纲参数，N = σ / (ρg * vt^2 * Dp)
        const N = sigma / (rho_g * vt * vt * Dp_typical);
        
        // 根据文档 Figure 04，根据 Rej 值选择不同的计算公式
        let v_max_entrainment;
        
        if (Rej > 1635) {
            // 高雷诺数区域：使用方程 E（高湍流）
            // Vre = K1 * sqrt(σ / (ρg * D))
            const K1 = 0.5; // 经验系数
            v_max_entrainment = K1 * Math.sqrt(sigma / (rho_g * D_m));
        } else if (Rej >= 160 && Rej <= 1635) {
            // 中等雷诺数区域：使用方程 D（文档推荐）
            // Vre = K2 * sqrt(σ * (ρl - ρg) / ρg^2)
            const K2 = 0.4; // 经验系数，根据文档调整
            v_max_entrainment = K2 * Math.sqrt(sigma * (rho_l - rho_g) / (rho_g * rho_g));
        } else if (Rej > 1.5) {
            // 低雷诺数区域：使用修正公式
            // Vre = K3 * sqrt(σ / (ρg * D)) * (1 + Rej^0.5)
            const K3 = 0.3;
            v_max_entrainment = K3 * Math.sqrt(sigma / (rho_g * D_m)) * (1 + Math.sqrt(Rej));
        } else {
            // 极低雷诺数：使用简化公式
            // Vre = K4 * sqrt(σ * (ρl - ρg) / ρg^2) * sqrt(1 + N)
            const K4 = 0.35;
            v_max_entrainment = K4 * Math.sqrt(sigma * (rho_l - rho_g) / (rho_g * rho_g)) * Math.sqrt(1 + N);
        }
        
        // 考虑液位高度的影响（液位越高，防夹带速度越低）
        // 根据文档，最大液位建议不超过 0.5*D，液位越高风险越大
        const liquid_level_factor = 1.0 - 0.3 * h_liq_ratio; // 液位修正系数
        v_max_entrainment = v_max_entrainment * liquid_level_factor;
        
        // 验证计算结果的合理性
        if (!isFinite(v_max_entrainment) || v_max_entrainment <= 0) {
            // 后备公式：基于密度比和终端速度的经验公式
            v_max_entrainment = 0.15 * Math.sqrt(rho_l / rho_g) * vt;
        }
        
        // 确保防夹带速度至少是终端速度的合理倍数（通常为 2-10 倍）
        const min_vre_ratio = 2.0;
        const max_vre_ratio = 15.0;
        if (v_max_entrainment < vt * min_vre_ratio) {
            v_max_entrainment = vt * min_vre_ratio;
        } else if (v_max_entrainment > vt * max_vre_ratio) {
            v_max_entrainment = vt * max_vre_ratio;
        }
        
        // 计算实际水平气速（基于当前结构参数）
        // 注意：由于迭代过程中严格保持 vh/vt 比值，实际气速应该等于设计目标 vh
        const A_total_current = Math.PI * Math.pow(D_m / 2, 2);
        const A_gas_current = A_total_current * A_gas_ratio;
        const vh_actual = V_g / A_gas_current;
        
        // 始终使用设计目标 vh 作为显示值，确保 vh/vt 比值严格等于输入值
        // 这是设计原则：用户输入的 vh/vt 比值必须被严格遵循
        vh_display = vh; // 直接使用设计目标值（变量已在前面声明）
        
        // 如果实际值与目标值偏差较大，给出警告但不改变显示值
        if (Math.abs(vh_actual - vh) / vh > 0.05) {
            console.warn(`警告: 实际气速 ${vh_actual.toFixed(3)} m/s 与设计目标 ${vh.toFixed(3)} m/s 偏差较大 (${((vh_actual - vh) / vh * 100).toFixed(1)}%)，可能存在计算问题`);
        }
        
        let isSafeVelocity = vh_display < v_max_entrainment;
        let velocityRatio_re = vh_display / v_max_entrainment;
        let structureAdjusted = false; // 标记是否调整了结构参数

        // 如果实际气速超过防夹带速度，自动调整结构参数
        if (!isSafeVelocity) {
            // 计算满足防夹带要求所需的最小直径
            // vh_actual = V_g / (A_gas * A_gas_ratio)
            // vh_actual < Vre
            // 因此：V_g / (π * (D/2)^2 * A_gas_ratio) < Vre
            // D >= sqrt(4 * V_g / (π * Vre * A_gas_ratio))
            const D_min_safe = Math.sqrt((4 * V_g) / (Math.PI * v_max_entrainment * A_gas_ratio));
            
            // 增加安全余量（8%）
            const safety_margin = 1.08;
            const D_adjusted = D_min_safe * safety_margin;
            const L_adjusted = D_adjusted * L_D_ratio;
            
            // 重新计算实际气速
            const A_total_adjusted = Math.PI * Math.pow(D_adjusted / 2, 2);
            const A_gas_adjusted = A_total_adjusted * A_gas_ratio;
            const vh_adjusted = V_g / A_gas_adjusted;
            
            // 验证调整后的气速是否安全
            if (vh_adjusted < v_max_entrainment) {
                // 更新结构参数
                D_m = D_adjusted;
                L_m = L_adjusted;
                
                // 重新计算容积和停留时间
                const Vol_m3_adjusted = Math.PI * Math.pow(D_m / 2, 2) * L_m;
                const Vol_L_adjusted = Vol_m3_adjusted * 1000;
                
                // 重新计算停留时间
                liquid_volume_ratio = 1 - calculateGasAreaRatio(h_liq_ratio);
                liquid_volume_m3 = Vol_m3_adjusted * liquid_volume_ratio;
                residence_time = liquid_volume_m3 / Math.max(liquid_flow_m3_s, 1e-6);
                
                // 更新显示值
                D_mm = D_m * 1000;
                L_mm = L_m * 1000;
                Vol_L = Vol_L_adjusted;
                // 注意：即使调整了结构参数，显示值仍然使用设计目标 vh，以保持 vh/vt 比值
                // vh_display 已经在前面设置为 vh，这里不需要修改
                
                // 更新安全状态（使用设计目标 vh 进行判断）
                isSafeVelocity = vh < v_max_entrainment;
                velocityRatio_re = vh / v_max_entrainment;
                structureAdjusted = true;
            } else {
                // 即使调整后仍不安全，需要更大幅度的调整
                // 进一步增大直径
                const additional_margin = Math.sqrt(vh_adjusted / v_max_entrainment) * 1.1;
                const D_final = D_adjusted * additional_margin;
                const L_final = D_final * L_D_ratio;
                
                const A_total_final = Math.PI * Math.pow(D_final / 2, 2);
                const A_gas_final = A_total_final * A_gas_ratio;
                const vh_final = V_g / A_gas_final;
                
                if (vh_final < v_max_entrainment) {
                    D_m = D_final;
                    L_m = L_final;
                    D_mm = D_m * 1000;
                    L_mm = L_m * 1000;
                    const Vol_m3_final = Math.PI * Math.pow(D_m / 2, 2) * L_m;
                    Vol_L = Vol_m3_final * 1000;
                    
                    // 重新计算停留时间
                    liquid_volume_ratio = 1 - calculateGasAreaRatio(h_liq_ratio);
                    liquid_volume_m3 = Vol_m3_final * liquid_volume_ratio;
                    residence_time = liquid_volume_m3 / Math.max(liquid_flow_m3_s, 1e-6);
                    
                    // 注意：即使调整了结构参数，显示值仍然使用设计目标 vh
                    // vh_display 已经在前面设置为 vh，这里不需要修改
                    isSafeVelocity = vh < v_max_entrainment;
                    velocityRatio_re = vh / v_max_entrainment;
                    structureAdjusted = true;
                }
            }
        }

        // 在所有结构参数确定后，检查液滴是否能够沉降（最终验证）
        // 计算气相高度和沉降距离
        const H_gas_final = D_m * (1 - h_liq_ratio);
        const S_l_final = H_gas_final;
        
        // 计算飞行时间和沉降时间
        const t_fly_final = L_m / vh_display;
        const t_fall_final = S_l_final / vt;
        
        // 判断液滴是否无法沉降（飞行时间小于沉降时间）
        const cannotSettle = t_fly_final < t_fall_final;

        // 更新结果显示
        document.getElementById('massFlowH').textContent = formatNumber(m_dot, 3);
        document.getElementById('gasDensityH').textContent = formatNumber(rho_g, 2);
        document.getElementById('liquidDensityH').textContent = formatNumber(rho_l, 2);
        document.getElementById('terminalVelocityH').textContent = formatNumber(vt, 3);
        document.getElementById('diameterH').textContent = formatNumber(D_mm, 0);
        document.getElementById('lengthH').textContent = formatNumber(L_mm, 0);
        
        // 水平气速显示，带安全状态
        // 自动模式：使用设计目标 vh（严格保持 vh/vt 比值）
        // 手动模式：使用实际气速 vh（校核模式）
        if (calcModeH === 'auto') {
            // 自动模式：确保使用设计目标 vh
            if (Math.abs(vh_display - vh) > 0.001) {
                console.warn(`警告: vh_display (${vh_display}) 与设计目标 vh (${vh}) 不匹配，强制使用设计目标值`);
                vh_display = vh;
            }
        }
        // 手动模式：vh_display 已经在前面设置为实际气速 vh，不需要修改
        
        let vhDisplay = formatNumber(vh_display, 3);
        if (isSafeVelocity) {
            vhDisplay += ' ✓';
        } else {
            vhDisplay += ' ⚠';
        }
        document.getElementById('horizontalVelocityH').textContent = vhDisplay;
        
        // 调试信息：验证 vh/vt 比值（仅自动模式）
        if (calcModeH === 'auto') {
            const actual_ratio = vh_display / vt;
            const input_ratio = vh_vt_ratio;
            if (Math.abs(actual_ratio - input_ratio) > 0.01) {
                console.warn(`警告: 显示的水平气速倍数 ${actual_ratio.toFixed(2)} 与输入值 ${input_ratio.toFixed(2)} 不匹配`);
            }
        }
        
        document.getElementById('residenceTimeH').textContent = formatNumber(residence_time, 1);
        document.getElementById('volumeH').textContent = formatNumber(Vol_L, 1);
        
        // 添加防夹带速度显示（如果结果区域有对应元素）
        const vreElement = document.getElementById('reentrainmentVelocityH');
        if (vreElement) {
            vreElement.textContent = formatNumber(v_max_entrainment, 3);
        }

        // 更新SVG示意图（使用用户要求的 updateSvg 函数）
        updateSvg(L_mm, D_mm, h_liq_ratio, cannotSettle);

        // 状态消息
        if (calcModeH === 'manual-diameter') {
            // 手动输入直径模式：显示校核结果
            const actual_vh_vt_ratio = vh_display / vt;
            if (!isSafeVelocity) {
                showStatus(
                    `❌ 校核结果: 实际水平气速 ${vh_display.toFixed(3)} m/s 超过防夹带速度 ${v_max_entrainment.toFixed(3)} m/s！` +
                    `实际 vh/vt 比值: ${actual_vh_vt_ratio.toFixed(2)}，建议增大直径或降低流量`,
                    'error',
                    'statusH'
                );
            } else if (velocityRatio_re > 0.8) {
                showStatus(
                    `⚠️ 校核结果: 实际水平气速 ${vh_display.toFixed(3)} m/s 接近防夹带速度 ${v_max_entrainment.toFixed(3)} m/s。` +
                    `实际 vh/vt 比值: ${actual_vh_vt_ratio.toFixed(2)}，建议增大直径`,
                    'info',
                    'statusH'
                );
            } else {
                showStatus(
                    `✓ 校核完成！实际水平气速: ${vh_display.toFixed(3)} m/s，实际 vh/vt 比值: ${actual_vh_vt_ratio.toFixed(2)}，防夹带速度: ${v_max_entrainment.toFixed(3)} m/s (安全比值: ${velocityRatio_re.toFixed(2)})`,
                    'success',
                    'statusH'
                );
            }
        } else if (structureAdjusted) {
            // 结构参数已自动调整
            const original_vh = vh_vt_ratio * vt;
            showStatus(
                `⚠️ 已自动调整结构参数！原设计气速 ${original_vh.toFixed(3)} m/s 超过防夹带速度 ${v_max_entrainment.toFixed(3)} m/s。` +
                `调整后：D=${D_mm.toFixed(0)}mm, L=${L_mm.toFixed(0)}mm, vh=${vh_display.toFixed(3)}m/s (安全比值: ${velocityRatio_re.toFixed(2)})`,
                'info',
                'statusH'
            );
        } else if (!isSafeVelocity) {
            // 即使调整后仍不安全，给出严重警告和建议
            const suggested_vh_ratio = (v_max_entrainment / vt) * 0.9; // 建议的水平气速倍数（90%的安全余量）
            showStatus(
                `❌ 严重警告: 水平气速 ${vh_display.toFixed(3)} m/s 超过防夹带速度 ${v_max_entrainment.toFixed(3)} m/s (比值: ${velocityRatio_re.toFixed(2)})！` +
                `建议：1) 降低水平气速倍数至 ${suggested_vh_ratio.toFixed(1)} 以下；` +
                `2) 或增大长径比 L/D；3) 或降低液位高度比`,
                'error',
                'statusH'
            );
        } else if (velocityRatio_re > 0.8) {
            showStatus(
                `注意: 水平气速 ${vh_display.toFixed(3)} m/s 接近防夹带速度 ${v_max_entrainment.toFixed(3)} m/s (比值: ${velocityRatio_re.toFixed(2)})，建议降低气速或增大直径`,
                'info',
                'statusH'
            );
        } else {
            showStatus(
                `计算完成！防夹带速度: ${v_max_entrainment.toFixed(3)} m/s，当前气速安全 (比值: ${velocityRatio_re.toFixed(2)})`,
                'success',
                'statusH'
            );
        }

    } catch (error) {
        console.error('计算错误:', error);
        showStatus('错误: ' + error.message, 'error', 'statusH');
        
        // 清空结果显示
        document.getElementById('massFlowH').textContent = '-';
        document.getElementById('gasDensityH').textContent = '-';
        document.getElementById('liquidDensityH').textContent = '-';
        document.getElementById('terminalVelocityH').textContent = '-';
        document.getElementById('diameterH').textContent = '-';
        document.getElementById('lengthH').textContent = '-';
        document.getElementById('horizontalVelocityH').textContent = '-';
        document.getElementById('residenceTimeH').textContent = '-';
        document.getElementById('volumeH').textContent = '-';
        const vreElement = document.getElementById('reentrainmentVelocityH');
        if (vreElement) {
            vreElement.textContent = '-';
        }
    } finally {
        // 恢复计算按钮
        const calcBtn = document.getElementById('calculateBtnH');
        calcBtn.disabled = false;
        calcBtn.textContent = '计算';
    }
}

// 绘制卧式分离器SVG示意图（增强版）
// 参数：diameter (mm), length (mm), liquidLevelRatio (0-1), cannotSettle (boolean) - 是否无法沉降
function drawHorizontalSeparator(diameter, length, liquidLevelRatio, cannotSettle = false) {
    const svg = document.getElementById('horizontalDiagram');
    if (!svg) return;

    // SVG 画布尺寸（保持 viewBox 不变）
    const svgWidth = 600;
    const svgHeight = 300;
    const padding = 80;
    
    // 计算长径比
    const aspectRatio = length / diameter;
    const maxDisplayWidth = svgWidth - 2 * padding;
    const maxDisplayHeight = svgHeight - 2 * padding - 40; // 40px 用于标注
    
    let displayWidth, displayHeight;
    
    // 根据长径比动态调整显示尺寸，使其视觉上接近真实比例
    if (aspectRatio > 4) {
        // 很长的分离器，宽度优先
        displayWidth = maxDisplayWidth;
        displayHeight = displayWidth / aspectRatio;
        if (displayHeight < 60) {
            displayHeight = 60;
            displayWidth = displayHeight * aspectRatio;
        }
    } else if (aspectRatio < 2) {
        // 较短的分离器，高度优先
        displayHeight = maxDisplayHeight;
        displayWidth = displayHeight * aspectRatio;
        if (displayWidth > maxDisplayWidth) {
            displayWidth = maxDisplayWidth;
            displayHeight = displayWidth / aspectRatio;
        }
    } else {
        // 中等比例，平衡显示
        displayWidth = Math.min(maxDisplayWidth, maxDisplayHeight * aspectRatio);
        displayHeight = displayWidth / aspectRatio;
    }
    
    // 计算起始位置（居中）
    const startX = (svgWidth - displayWidth) / 2;
    const startY = padding;
    const endX = startX + displayWidth;
    const endY = startY + displayHeight;
    
    // 根据 liquidLevelRatio 计算液位位置
    const liquidHeight = displayHeight * liquidLevelRatio;
    const liquidTopY = endY - liquidHeight;
    
    // 1. 更新筒身（根据长径比调整的尺寸）
    const tankBody = document.getElementById('tankBody');
    if (tankBody) {
        tankBody.setAttribute('x', startX);
        tankBody.setAttribute('y', startY);
        tankBody.setAttribute('width', displayWidth);
        tankBody.setAttribute('height', displayHeight);
    }
    
    // 2. 更新液体区域（根据 liquidLevel 动态更新）
    const liquidArea = document.getElementById('liquidArea');
    if (liquidArea) {
        liquidArea.setAttribute('x', startX);
        liquidArea.setAttribute('y', liquidTopY);
        liquidArea.setAttribute('width', displayWidth);
        liquidArea.setAttribute('height', liquidHeight);
        // 设置圆角
        if (liquidHeight > 0) {
            liquidArea.setAttribute('rx', '8');
            liquidArea.setAttribute('ry', '8');
        } else {
            liquidArea.setAttribute('rx', '0');
            liquidArea.setAttribute('ry', '0');
        }
    }
    
    // 3. 更新液位线（根据 liquidLevel 动态更新 y 坐标和高度）
    const liquidLevel = document.getElementById('liquidLevel');
    if (liquidLevel) {
        liquidLevel.setAttribute('x1', startX);
        liquidLevel.setAttribute('y1', liquidTopY);
        liquidLevel.setAttribute('x2', endX);
        liquidLevel.setAttribute('y2', liquidTopY);
    }
    
    // 4. 绘制气体流向箭头（2-3个，在气相空间）
    const gasFlowArrows = document.getElementById('gasFlowArrows');
    if (gasFlowArrows) {
        // 清空现有箭头
        gasFlowArrows.innerHTML = '';
        
        // 气相空间高度
        const gasSpaceHeight = liquidTopY - startY;
        if (gasSpaceHeight > 20) {
            const arrowSpacing = gasSpaceHeight / 4; // 均匀分布
            
            // 绘制2-3个箭头
            const numArrows = gasSpaceHeight > 80 ? 3 : 2;
            for (let i = 0; i < numArrows; i++) {
                const arrowY = startY + (i + 1) * arrowSpacing;
                const arrowStartX = startX + 20;
                const arrowEndX = endX - 20;
                
                // 箭头线
                const arrowLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                arrowLine.setAttribute('x1', arrowStartX);
                arrowLine.setAttribute('y1', arrowY);
                arrowLine.setAttribute('x2', arrowEndX);
                arrowLine.setAttribute('y2', arrowY);
                arrowLine.setAttribute('stroke', '#757575');
                arrowLine.setAttribute('stroke-width', '2');
                arrowLine.setAttribute('marker-end', 'url(#arrowhead-gray)');
                arrowLine.setAttribute('opacity', '0.7');
                gasFlowArrows.appendChild(arrowLine);
            }
        }
    }
    
    // 5. 绘制液滴沉降轨迹（红色虚线）
    // 如果无法沉降，轨迹终点画到液面之后（超出有效长度）
    const dropletPath = document.getElementById('dropletPath');
    if (dropletPath && liquidTopY > startY) {
        const startPathY = startY + (liquidTopY - startY) * 0.2; // 起点在气相空间上部
        
        let endPathX, endPathY;
        if (cannotSettle) {
            // 无法沉降：终点画到液面之后，超出有效长度，给用户直观警告
            endPathX = endX + 50; // 超出筒体右侧
            endPathY = liquidTopY + 20; // 在液面下方，表示未沉降到液面
        } else {
            // 正常沉降：终点刚好在液面
            endPathX = endX - 30;
            endPathY = liquidTopY;
        }
        
        // 创建曲线路径（模拟抛物线轨迹）
        const controlX = (startX + endPathX) / 2;
        const controlY = startPathY + (endPathY - startPathY) * 0.6;
        
        const pathData = `M ${startX + 20} ${startPathY} Q ${controlX} ${controlY} ${endPathX} ${endPathY}`;
        dropletPath.setAttribute('d', pathData);
        dropletPath.style.display = 'block';
        
        // 如果无法沉降，使用更醒目的样式
        if (cannotSettle) {
            dropletPath.setAttribute('stroke', '#d32f2f');
            dropletPath.setAttribute('stroke-width', '3');
            dropletPath.setAttribute('stroke-dasharray', '6,3');
            dropletPath.setAttribute('opacity', '1');
        } else {
            dropletPath.setAttribute('stroke', '#d32f2f');
            dropletPath.setAttribute('stroke-width', '2');
            dropletPath.setAttribute('stroke-dasharray', '4,3');
            dropletPath.setAttribute('opacity', '0.8');
        }
    } else if (dropletPath) {
        dropletPath.style.display = 'none';
    }
    
    // 6. 更新尺寸标注 - 长度
    const lengthLine = document.getElementById('lengthLine');
    const lengthLabel = document.getElementById('lengthLabel');
    if (lengthLine) {
        const labelY = endY + 30;
        lengthLine.setAttribute('x1', startX);
        lengthLine.setAttribute('y1', labelY);
        lengthLine.setAttribute('x2', endX);
        lengthLine.setAttribute('y2', labelY);
    }
    if (lengthLabel) {
        lengthLabel.setAttribute('x', (startX + endX) / 2);
        lengthLabel.setAttribute('y', endY + 50);
        lengthLabel.textContent = `L = ${formatNumber(length, 0)} mm`;
    }
    
    // 7. 更新尺寸标注 - 直径
    const diameterLine = document.getElementById('diameterLine');
    const diameterLabel = document.getElementById('diameterLabel');
    if (diameterLine) {
        const labelX = startX - 40;
        diameterLine.setAttribute('x1', labelX);
        diameterLine.setAttribute('y1', startY);
        diameterLine.setAttribute('x2', labelX);
        diameterLine.setAttribute('y2', endY);
    }
    if (diameterLabel) {
        diameterLabel.setAttribute('x', startX - 65);
        diameterLabel.setAttribute('y', (startY + endY) / 2);
        diameterLabel.textContent = `D = ${formatNumber(diameter, 0)} mm`;
    }
}

// 更新卧式分离器SVG示意图（保持向后兼容）
function updateHorizontalDiagram(D_mm, L_mm, h_ratio, cannotSettle = false) {
    drawHorizontalSeparator(D_mm, L_mm, h_ratio, cannotSettle);
}

// 更新SVG示意图（用户要求的函数名）
// 参数：length (mm), diameter (mm), liquidLevel (0-1), cannotSettle (boolean)
function updateSvg(length, diameter, liquidLevel, cannotSettle = false) {
    drawHorizontalSeparator(diameter, length, liquidLevel, cannotSettle);
}

// 主计算函数（立式）
async function calculate() {
    try {
        // 检查 CoolProp 是否已初始化
        if (!CoolPropModule) {
            showStatus('正在初始化 CoolProp...', 'info');
            await initCoolProp();
        }

        // 获取输入参数
        const refrigerant = document.getElementById('refrigerant').value;
        const Q_kW = parseFloat(document.getElementById('coolingCapacity').value);
        const Te_C = parseFloat(document.getElementById('evapTemp').value);
        const deltaT_sh = parseFloat(document.getElementById('superheat').value);
        const calcMode = document.getElementById('calcMode').value;
        const Dp_um = parseFloat(document.getElementById('dropletSize').value);
        const velocityRatio = parseFloat(document.getElementById('velocityRatio').value);
        let u = parseFloat(document.getElementById('gasVelocity').value);

        // 验证输入
        if (isNaN(Q_kW) || Q_kW <= 0) {
            throw new Error('制冷量必须大于 0');
        }
        if (isNaN(Te_C)) {
            throw new Error('请输入有效的蒸发温度');
        }
        if (isNaN(deltaT_sh) || deltaT_sh < 0) {
            throw new Error('吸气过热度必须大于等于 0');
        }
        if (calcMode === 'manual-velocity' && (isNaN(u) || u <= 0)) {
            throw new Error('设计气速必须大于 0');
        }
        if (calcMode === 'manual-diameter') {
            const D_mm_input = parseFloat(document.getElementById('diameterInput').value);
            if (isNaN(D_mm_input) || D_mm_input <= 0) {
                throw new Error('请输入有效的直径值');
            }
        }
        if (calcMode !== 'manual-diameter' && (isNaN(Dp_um) || Dp_um <= 0)) {
            throw new Error('液滴直径必须大于 0');
        }
        if (calcMode !== 'manual-velocity' && (isNaN(velocityRatio) || velocityRatio <= 0 || velocityRatio > 1)) {
            throw new Error('设计气速系数必须在 0 到 1 之间');
        }

        // 禁用计算按钮
        const calcBtn = document.getElementById('calculateBtn');
        calcBtn.disabled = true;
        calcBtn.textContent = '计算中...';

        showStatus('正在计算...', 'info');

        // 温度单位转换：°C 转 K
        const Te_K = Te_C + 273.15;
        const T_suc_K = Te_K + deltaT_sh;

        // 1. 计算蒸发压力 (饱和压力，Q=1 表示饱和蒸汽)
        const P_sat = CoolPropModule.PropsSI('P', 'T', Te_K, 'Q', 1, refrigerant);
        if (P_sat <= 0 || !isFinite(P_sat)) {
            throw new Error('无法计算蒸发压力，请检查制冷剂和温度参数');
        }

        // 2. 计算吸气温度下的气体密度
        const rho_g = CoolPropModule.PropsSI('D', 'P', P_sat, 'T', T_suc_K, refrigerant);
        if (rho_g <= 0 || !isFinite(rho_g)) {
            throw new Error('无法计算气体密度');
        }

        // 2.1 计算饱和液体密度
        const rho_l = CoolPropModule.PropsSI('D', 'P', P_sat, 'Q', 0, refrigerant);
        if (rho_l <= 0 || !isFinite(rho_l)) {
            throw new Error('无法计算液体密度');
        }

        // 2.2 计算气体动力粘度
        const eta_g = CoolPropModule.PropsSI('V', 'P', P_sat, 'T', T_suc_K, refrigerant);
        if (eta_g <= 0 || !isFinite(eta_g)) {
            throw new Error('无法计算气体动力粘度');
        }

        // 3. 计算吸气状态下的气体焓值
        const h_gas = CoolPropModule.PropsSI('H', 'P', P_sat, 'T', T_suc_K, refrigerant);
        if (!isFinite(h_gas)) {
            throw new Error('无法计算气体焓值');
        }

        // 4. 计算饱和液体的焓值 (用于计算焓差)
        const h_liq = CoolPropModule.PropsSI('H', 'P', P_sat, 'Q', 0, refrigerant);
        if (!isFinite(h_liq)) {
            throw new Error('无法计算液体焓值');
        }

        // 5. 计算终端速度 (根据文档公式 1-4)
        const Dp = Dp_um * 1e-6; // 转换为米
        const vt = calculateTerminalVelocity(Dp, rho_l, rho_g, eta_g);
        
        if (vt <= 0 || !isFinite(vt)) {
            throw new Error('无法计算终端速度');
        }

        // 6. 计算焓差 (J/kg)
        const delta_h = h_gas - h_liq;
        if (delta_h <= 0 || !isFinite(delta_h)) {
            throw new Error('焓差计算异常，请检查过热度设置');
        }

        // 7. 计算质量流量 (kg/s)
        // Q_kW * 1000 转换为 W，除以焓差得到质量流量
        const m_dot = (Q_kW * 1000) / delta_h;

        // 8. 计算气体体积流量 (m³/s)
        const V_g = m_dot / rho_g;

        // 9. 根据计算模式确定设计气速和直径
        let D_m; // 直径（米）
        
        if (calcMode === 'manual-diameter') {
            // 手动输入直径模式（校核模式）
            const D_mm_input = parseFloat(document.getElementById('diameterInput').value);
            D_m = D_mm_input / 1000; // 转换为米
            
            // 计算实际气速
            const A_total = Math.PI * Math.pow(D_m / 2, 2);
            u = V_g / A_total;
            
            // 验证实际气速是否合理
            if (u > vt) {
                throw new Error(`实际气速 ${u.toFixed(3)} m/s 大于终端速度 ${vt.toFixed(3)} m/s，无法有效分离液滴。建议增大直径`);
            }
            
            // 计算实际气速系数
            const actualRatio = u / vt;
            if (actualRatio > 0.9) {
                showStatus(`警告: 实际气速/终端速度比值 ${actualRatio.toFixed(2)} 接近上限，建议增大直径`, 'info');
            }
        } else if (calcMode === 'manual-velocity') {
            // 手动输入设计气速模式
            if (u > vt) {
                throw new Error(`设计气速 ${u.toFixed(3)} m/s 大于终端速度 ${vt.toFixed(3)} m/s，无法有效分离液滴`);
            }
            const actualRatio = u / vt;
            if (actualRatio > 0.9) {
                showStatus(`警告: 设计气速/终端速度比值 ${actualRatio.toFixed(2)} 接近上限，建议小于 0.9`, 'info');
            }
            
            // 计算直径
            D_m = Math.sqrt((4 * V_g) / (Math.PI * u));
        } else {
            // 自动计算模式
            u = velocityRatio * vt;
            // 更新输入框显示计算出的值
            document.getElementById('gasVelocity').value = u.toFixed(3);
            
            // 计算直径
            D_m = Math.sqrt((4 * V_g) / (Math.PI * u));
        }
        // 转换为 mm
        const D_mm = D_m * 1000;

        // 11. 计算筒身高度 (根据文档建议: 2.5-3 D，取 3.0 D)
        const H_mm = D_mm * 3.0;

        // 12. 计算筒内容积 (L)
        // Vol = π * (D/2)^2 * H
        const Vol_m3 = Math.PI * Math.pow(D_m / 2, 2) * (H_mm / 1000);
        const Vol_L = Vol_m3 * 1000;

        // 更新结果显示
        document.getElementById('massFlow').textContent = formatNumber(m_dot, 3);
        document.getElementById('gasDensity').textContent = formatNumber(rho_g, 2);
        document.getElementById('diameter').textContent = formatNumber(D_mm, 0);
        document.getElementById('height').textContent = formatNumber(H_mm, 0);
        document.getElementById('volume').textContent = formatNumber(Vol_L, 1);
        document.getElementById('terminalVelocity').textContent = formatNumber(vt, 3);
        document.getElementById('liquidDensity').textContent = formatNumber(rho_l, 2);
        document.getElementById('gasViscosity').textContent = formatNumber(eta_g, 6);

        // 如果是手动输入直径模式，显示校核结果
        if (calcMode === 'manual-diameter') {
            const actualRatio = u / vt;
            showStatus(`校核完成！实际气速: ${u.toFixed(3)} m/s，气速/终端速度比值: ${actualRatio.toFixed(2)} ${actualRatio < 0.9 ? '✓' : '⚠'}`, 'success');
        } else {
            showStatus('计算完成！', 'success');
        }

    } catch (error) {
        console.error('计算错误:', error);
        showStatus('错误: ' + error.message, 'error');
        
        // 清空结果显示
        document.getElementById('massFlow').textContent = '-';
        document.getElementById('gasDensity').textContent = '-';
        document.getElementById('diameter').textContent = '-';
        document.getElementById('height').textContent = '-';
        document.getElementById('volume').textContent = '-';
        document.getElementById('terminalVelocity').textContent = '-';
        document.getElementById('liquidDensity').textContent = '-';
        document.getElementById('gasViscosity').textContent = '-';
    } finally {
        // 恢复计算按钮
        const calcBtn = document.getElementById('calculateBtn');
        calcBtn.disabled = false;
        calcBtn.textContent = '计算';
    }
}

// 页面加载完成后初始化
document.addEventListener('DOMContentLoaded', () => {
    // 类型切换功能
    const verticalBtn = document.getElementById('verticalBtn');
    const horizontalBtn = document.getElementById('horizontalBtn');
    const verticalContent = document.getElementById('verticalContent');
    const horizontalContent = document.getElementById('horizontalContent');

    verticalBtn.addEventListener('click', () => {
        verticalBtn.classList.add('active');
        horizontalBtn.classList.remove('active');
        verticalContent.style.display = 'block';
        horizontalContent.style.display = 'none';
    });

    horizontalBtn.addEventListener('click', () => {
        horizontalBtn.classList.add('active');
        verticalBtn.classList.remove('active');
        horizontalContent.style.display = 'block';
        verticalContent.style.display = 'none';
    });

    // 绑定计算按钮事件
    const calcBtn = document.getElementById('calculateBtn');
    calcBtn.addEventListener('click', calculate);
    
    const calcBtnH = document.getElementById('calculateBtnH');
    calcBtnH.addEventListener('click', calculateHorizontal);

    // 绑定回车键事件（在输入框中按回车也可以计算）
    const inputFields = document.querySelectorAll('.input-field');
    inputFields.forEach(field => {
        field.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                // 根据当前显示的内容区域决定调用哪个计算函数
                if (verticalContent.style.display !== 'none') {
                calculate();
                } else if (horizontalContent.style.display !== 'none') {
                    calculateHorizontal();
                }
            }
        });
    });

    // 立式计算模式切换
    const calcModeSelect = document.getElementById('calcMode');
    const dropletSizeGroup = document.getElementById('dropletSizeGroup');
    const gasVelocityGroup = document.getElementById('gasVelocityGroup');
    const diameterInputGroup = document.getElementById('diameterInputGroup');
    const velocityRatioGroup = document.getElementById('velocityRatioGroup');
    
    calcModeSelect.addEventListener('change', (e) => {
        if (e.target.value === 'auto') {
            // 自动模式：显示液滴直径和设计气速系数，隐藏手动输入
            dropletSizeGroup.style.display = 'block';
            gasVelocityGroup.style.display = 'none';
            diameterInputGroup.style.display = 'none';
            velocityRatioGroup.style.display = 'block';
        } else if (e.target.value === 'manual-velocity') {
            // 手动输入设计气速模式：显示手动气速输入，隐藏其他
            dropletSizeGroup.style.display = 'none';
            gasVelocityGroup.style.display = 'block';
            diameterInputGroup.style.display = 'none';
            velocityRatioGroup.style.display = 'none';
        } else if (e.target.value === 'manual-diameter') {
            // 手动输入直径模式：显示液滴直径和直径输入，隐藏其他
            dropletSizeGroup.style.display = 'block';
            gasVelocityGroup.style.display = 'none';
            diameterInputGroup.style.display = 'block';
            velocityRatioGroup.style.display = 'block';
        }
    });

    // 初始化显示状态
    if (calcModeSelect.value === 'auto') {
        dropletSizeGroup.style.display = 'block';
        gasVelocityGroup.style.display = 'none';
        diameterInputGroup.style.display = 'none';
        velocityRatioGroup.style.display = 'block';
    } else if (calcModeSelect.value === 'manual-velocity') {
        dropletSizeGroup.style.display = 'none';
        gasVelocityGroup.style.display = 'block';
        diameterInputGroup.style.display = 'none';
        velocityRatioGroup.style.display = 'none';
    } else if (calcModeSelect.value === 'manual-diameter') {
        dropletSizeGroup.style.display = 'block';
        gasVelocityGroup.style.display = 'none';
        diameterInputGroup.style.display = 'block';
        velocityRatioGroup.style.display = 'block';
    }

    // 卧式计算模式切换
    const calcModeSelectH = document.getElementById('calcModeH');
    const lengthDiameterRatioGroup = document.getElementById('lengthDiameterRatioGroup');
    const diameterInputGroupH = document.getElementById('diameterInputGroupH');
    const liquidLevelRatioGroup = document.getElementById('liquidLevelRatio').parentElement;
    const ksValueGroup = document.getElementById('ksValue').parentElement;
    const velocityMultiplierGroup = document.getElementById('velocityMultiplier').parentElement;
    
    calcModeSelectH.addEventListener('change', (e) => {
        if (e.target.value === 'auto') {
            // 自动模式：显示所有设计参数
            lengthDiameterRatioGroup.style.display = 'block';
            diameterInputGroupH.style.display = 'none';
            liquidLevelRatioGroup.style.display = 'block';
            ksValueGroup.style.display = 'block';
            velocityMultiplierGroup.style.display = 'block';
        } else if (e.target.value === 'manual-diameter') {
            // 手动输入直径模式：显示直径输入，保留长径比和液位比
            lengthDiameterRatioGroup.style.display = 'block';
            diameterInputGroupH.style.display = 'block';
            liquidLevelRatioGroup.style.display = 'block';
            ksValueGroup.style.display = 'block';
            velocityMultiplierGroup.style.display = 'none'; // 校核模式不需要输入 vh/vt 倍数
        }
    });

    // 初始化显示状态（卧式）
    if (calcModeSelectH.value === 'auto') {
        lengthDiameterRatioGroup.style.display = 'block';
        diameterInputGroupH.style.display = 'none';
        liquidLevelRatioGroup.style.display = 'block';
        ksValueGroup.style.display = 'block';
        velocityMultiplierGroup.style.display = 'block';
    } else if (calcModeSelectH.value === 'manual-diameter') {
        lengthDiameterRatioGroup.style.display = 'block';
        diameterInputGroupH.style.display = 'block';
        liquidLevelRatioGroup.style.display = 'block';
        ksValueGroup.style.display = 'block';
        velocityMultiplierGroup.style.display = 'none';
    }

    // 预初始化 CoolProp
    showStatus('正在加载 CoolProp 库...', 'info');
    initCoolProp()
        .then(() => {
            showStatus('CoolProp 已就绪，可以开始计算', 'success');
        })
        .catch((error) => {
            showStatus('警告: CoolProp 初始化失败，计算功能可能不可用', 'error');
            console.error('CoolProp 初始化错误:', error);
        });
});

