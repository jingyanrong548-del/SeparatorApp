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
        // 如果过热度为 0，使用 0.001 进行计算
        const deltaT_sh_calc = deltaT_sh === 0 ? 0.001 : deltaT_sh;
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
        const T_suc_K = Te_K + deltaT_sh_calc;

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
        document.getElementById('volumeFlowH').textContent = formatNumber(V_g, 6);
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

        // 评估设计并显示结论
        // 在校核模式下，获取用户输入的目标防夹带速度比值
        let targetVelocityRatio = null;
        if (calcModeH === 'manual-diameter') {
            const targetInput = document.getElementById('targetVelocityRatioH');
            if (targetInput) {
                const targetValue = parseFloat(targetInput.value);
                if (!isNaN(targetValue) && targetValue > 0 && targetValue < 1.0) {
                    targetVelocityRatio = targetValue;
                }
            }
        }
        const evaluationH = evaluateHorizontalDesign(vh_display, vt, v_max_entrainment, D_mm, L_mm, L_D_ratio, h_liq_ratio, Vol_L, m_dot, rho_g, cannotSettle, calcModeH, targetVelocityRatio);
        const conclusionHTMLH = generateConclusionHTML(evaluationH);
        
        const conclusionSectionH = document.getElementById('conclusionSectionH');
        const conclusionContentH = document.getElementById('conclusionContentH');
        if (conclusionSectionH && conclusionContentH) {
            conclusionContentH.innerHTML = conclusionHTMLH;
            conclusionSectionH.style.display = 'block';
        }

        // 根据评估结果显示状态消息（不合理的设计不会自动隐藏）
        if (calcModeH === 'manual-diameter') {
            const actual_vh_vt_ratio = vh_display / vt;
            if (evaluationH.overallStatus === 'unqualified') {
                showStatus(`❌ 校核结果：设计不合格！请查看下方结论与建议进行调整。`, 'error', 'statusH');
            } else if (evaluationH.overallStatus === 'warning') {
                showStatus(`⚠️ 校核结果：设计需注意！请查看下方结论与建议。`, 'error', 'statusH');
            } else {
                showStatus(`✓ 校核完成！设计合格。`, 'success', 'statusH');
            }
        } else {
            if (evaluationH.overallStatus === 'unqualified') {
                showStatus('❌ 计算完成，但设计不合格！请查看下方结论与建议进行调整。', 'error', 'statusH');
            } else if (evaluationH.overallStatus === 'warning') {
                showStatus('⚠️ 计算完成，但设计需注意！请查看下方结论与建议。', 'error', 'statusH');
            } else {
                showStatus('✓ 计算完成！设计合格。', 'success', 'statusH');
            }
        }

    } catch (error) {
        console.error('计算错误:', error);
        showStatus('错误: ' + error.message, 'error', 'statusH');
        
        // 隐藏结论部分
        const conclusionSectionH = document.getElementById('conclusionSectionH');
        if (conclusionSectionH) {
            conclusionSectionH.style.display = 'none';
        }
        
        // 清空结果显示
        document.getElementById('massFlowH').textContent = '-';
        document.getElementById('volumeFlowH').textContent = '-';
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

// 评估立式分离器设计
// 参数：u (设计气速), vt (终端速度), D_mm (直径mm), H_mm (高度mm), Vol_L (容积L), 
//      m_dot (质量流量kg/s), rho_g (气体密度kg/m³), calcMode (计算模式), velocityRatio (气速系数，仅自动模式)
function evaluateVerticalDesign(u, vt, D_mm, H_mm, Vol_L, m_dot, rho_g, calcMode, velocityRatio = null) {
    const evaluation = {
        overallStatus: 'qualified', // 'qualified', 'warning', 'unqualified'
        items: [],
        suggestions: []
    };

    // 1. 气速比值评估
    const velocityRatio_actual = u / vt;
    let velocityStatus = 'qualified';
    let velocityStatusText = '';
    
    if (velocityRatio_actual > 1.0) {
        velocityStatus = 'unqualified';
        velocityStatusText = '危险：气速超过终端速度，无法有效分离液滴';
    } else if (velocityRatio_actual > 0.95) {
        velocityStatus = 'unqualified';
        velocityStatusText = '危险：气速接近终端速度，分离效果差';
    } else if (velocityRatio_actual > 0.90) {
        velocityStatus = 'warning';
        velocityStatusText = '警告：气速接近上限，建议降低';
    } else if (velocityRatio_actual >= 0.85) {
        velocityStatus = 'qualified';
        velocityStatusText = '良好：气速在推荐范围内';
    } else if (velocityRatio_actual >= 0.75) {
        velocityStatus = 'qualified';
        velocityStatusText = '优秀：气速在最佳范围内';
    } else {
        velocityStatus = 'warning';
        velocityStatusText = '注意：气速较低，设计较保守，可考虑优化';
    }
    
    evaluation.items.push({
        name: '气速比值 (u/vt)',
        value: velocityRatio_actual.toFixed(3),
        status: velocityStatus,
        description: velocityStatusText
    });
    
    if (velocityStatus === 'unqualified') {
        evaluation.overallStatus = 'unqualified';
    } else if (velocityStatus === 'warning' && evaluation.overallStatus === 'qualified') {
        evaluation.overallStatus = 'warning';
    }

    // 2. 直径标准化检查
    const standardDiameters = [100, 150, 200, 250, 300, 400, 500, 600, 800, 1000];
    let closestStandard = null;
    let minDiff = Infinity;
    
    for (const stdD of standardDiameters) {
        const diff = Math.abs(D_mm - stdD);
        if (diff < minDiff) {
            minDiff = diff;
            closestStandard = stdD;
        }
    }
    
    const diameterDeviation = (minDiff / closestStandard) * 100;
    let diameterStatus = 'qualified';
    let diameterStatusText = '';
    
    if (diameterDeviation <= 5) {
        diameterStatus = 'qualified';
        diameterStatusText = `标准直径：${closestStandard} mm（偏差 ${diameterDeviation.toFixed(1)}%）`;
    } else if (diameterDeviation <= 10) {
        diameterStatus = 'warning';
        diameterStatusText = `接近标准直径：${closestStandard} mm（偏差 ${diameterDeviation.toFixed(1)}%）`;
    } else {
        diameterStatus = 'warning';
        diameterStatusText = `非标准直径，建议使用：${closestStandard} mm`;
    }
    
    evaluation.items.push({
        name: '直径标准化',
        value: `${D_mm.toFixed(0)} mm`,
        status: diameterStatus,
        description: diameterStatusText,
        suggestedValue: diameterDeviation > 5 ? `${closestStandard} mm` : null
    });
    
    if (diameterStatus === 'warning' && evaluation.overallStatus === 'qualified') {
        evaluation.overallStatus = 'warning';
    }

    // 3. 高度范围验证
    const H_D_ratio = H_mm / D_mm;
    let heightStatus = 'qualified';
    let heightStatusText = '';
    
    if (H_D_ratio >= 2.5 && H_D_ratio <= 3.0) {
        heightStatus = 'qualified';
        heightStatusText = '推荐：高度在最佳范围内';
    } else if (H_D_ratio >= 2.0 && H_D_ratio < 2.5) {
        heightStatus = 'warning';
        heightStatusText = '可接受：高度略低，建议 H/D = 2.5-3.0';
    } else if (H_D_ratio > 3.0 && H_D_ratio <= 3.5) {
        heightStatus = 'warning';
        heightStatusText = '可接受：高度略高，建议 H/D = 2.5-3.0';
    } else if (H_D_ratio < 2.0) {
        heightStatus = 'warning';
        heightStatusText = '不推荐：高度过低，建议 H/D ≥ 2.5';
    } else {
        heightStatus = 'warning';
        heightStatusText = '不推荐：高度过高，建议 H/D ≤ 3.5';
    }
    
    const suggestedHeight = D_mm * 2.75; // 推荐使用 2.75 倍直径
    
    evaluation.items.push({
        name: '高度范围 (H/D)',
        value: H_D_ratio.toFixed(2),
        status: heightStatus,
        description: heightStatusText,
        suggestedValue: (H_D_ratio < 2.5 || H_D_ratio > 3.0) ? `${suggestedHeight.toFixed(0)} mm` : null
    });
    
    if (heightStatus === 'warning' && evaluation.overallStatus === 'qualified') {
        evaluation.overallStatus = 'warning';
    }

    // 4. 容积合理性评估
    // 基于质量流量估算最小停留时间（假设需要至少 5 秒停留时间）
    const minResidenceTime = 5; // 秒
    const minVolume_L = m_dot * minResidenceTime * 0.1; // 粗略估算，假设液体占10%
    
    let volumeStatus = 'qualified';
    let volumeStatusText = '';
    
    if (Vol_L >= minVolume_L * 2) {
        volumeStatus = 'qualified';
        volumeStatusText = '充足：容积满足停留时间要求';
    } else if (Vol_L >= minVolume_L) {
        volumeStatus = 'qualified';
        volumeStatusText = '合理：容积基本满足要求';
    } else {
        volumeStatus = 'warning';
        volumeStatusText = `注意：容积可能偏小，建议至少 ${(minVolume_L * 1.5).toFixed(1)} L`;
    }
    
    evaluation.items.push({
        name: '容积合理性',
        value: `${Vol_L.toFixed(1)} L`,
        status: volumeStatus,
        description: volumeStatusText
    });

    // 5. 生成具体建议
    if (velocityRatio_actual > 0.90) {
        if (calcMode === 'auto') {
            evaluation.suggestions.push(`降低设计气速系数至 ${(0.85 * vt / vt).toFixed(2)} 以下（当前：${velocityRatio?.toFixed(2) || velocityRatio_actual.toFixed(2)}）`);
        } else if (calcMode === 'manual-velocity') {
            evaluation.suggestions.push(`降低设计气速至 ${(0.85 * vt).toFixed(3)} m/s 以下（当前：${u.toFixed(3)} m/s）`);
        } else {
            const V_g = m_dot / rho_g;
            evaluation.suggestions.push(`增大直径至 ${(Math.sqrt((4 * V_g) / (Math.PI * 0.85 * vt)) * 1000).toFixed(0)} mm 以上（当前：${D_mm.toFixed(0)} mm）`);
        }
    }
    
    if (diameterDeviation > 5) {
        evaluation.suggestions.push(`使用标准直径 ${closestStandard} mm（当前：${D_mm.toFixed(0)} mm）`);
    }
    
    if (H_D_ratio < 2.5 || H_D_ratio > 3.0) {
        evaluation.suggestions.push(`调整高度至 ${suggestedHeight.toFixed(0)} mm（当前：${H_mm.toFixed(0)} mm，H/D = ${H_D_ratio.toFixed(2)}）`);
    }
    
    if (Vol_L < minVolume_L) {
        const suggestedD = Math.sqrt((minVolume_L * 1.5 * 1000) / (Math.PI * (H_mm / 1000) * 0.25));
        evaluation.suggestions.push(`增大直径或高度以增加容积（建议容积：${(minVolume_L * 1.5).toFixed(1)} L）`);
    }

    return evaluation;
}

// 生成结论HTML内容
function generateConclusionHTML(evaluation) {
    let html = '';
    
    // 总体状态
    const statusClass = evaluation.overallStatus === 'qualified' ? 'conclusion-qualified' : 
                       evaluation.overallStatus === 'warning' ? 'conclusion-warning' : 
                       'conclusion-unqualified';
    const statusText = evaluation.overallStatus === 'qualified' ? '✓ 设计合格' : 
                      evaluation.overallStatus === 'warning' ? '⚠ 设计需注意' : 
                      '❌ 设计不合格';
    
    html += `<div class="conclusion-overall ${statusClass}">${statusText}</div>`;
    
    // 评估项目
    html += '<div class="conclusion-items">';
    evaluation.items.forEach(item => {
        const itemStatusClass = item.status === 'qualified' ? 'conclusion-item-qualified' : 
                               item.status === 'warning' ? 'conclusion-item-warning' : 
                               'conclusion-item-unqualified';
        html += `<div class="conclusion-item ${itemStatusClass}">`;
        html += `<div class="conclusion-item-header">`;
        html += `<span class="conclusion-item-name">${item.name}</span>`;
        html += `<span class="conclusion-item-value">${item.value}</span>`;
        html += `</div>`;
        html += `<div class="conclusion-item-desc">${item.description}`;
        if (item.suggestedValue) {
            html += ` <strong>建议值：${item.suggestedValue}</strong>`;
        }
        html += `</div>`;
        html += `</div>`;
    });
    html += '</div>';
    
    // 具体建议
    if (evaluation.suggestions.length > 0) {
        html += '<div class="conclusion-suggestions">';
        html += '<h4>调整建议：</h4>';
        html += '<ul>';
        evaluation.suggestions.forEach(suggestion => {
            html += `<li>${suggestion}</li>`;
        });
        html += '</ul>';
        html += '</div>';
    }
    
    return html;
}

// 评估卧式分离器设计
// 参数：vh (水平气速), vt (沉降速度), v_max_entrainment (防夹带速度), D_mm (直径mm), L_mm (长度mm), 
//      L_D_ratio (长径比), h_liq_ratio (液位高度比), Vol_L (容积L), m_dot (质量流量kg/s), 
//      rho_g (气体密度), cannotSettle (是否无法沉降), calcMode (计算模式), targetVelocityRatio (目标防夹带速度比值，可选)
function evaluateHorizontalDesign(vh, vt, v_max_entrainment, D_mm, L_mm, L_D_ratio, h_liq_ratio, Vol_L, m_dot, rho_g, cannotSettle, calcMode, targetVelocityRatio = null) {
    const evaluation = {
        overallStatus: 'qualified',
        items: [],
        suggestions: []
    };

    // 1. 防夹带速度评估（最重要）
    const velocityRatio_re = vh / v_max_entrainment;
    let entrainmentStatus = 'qualified';
    let entrainmentStatusText = '';
    
    // 如果提供了目标值，使用目标值进行评估
    const targetRatio = targetVelocityRatio !== null ? targetVelocityRatio : 0.8;
    const deviationFromTarget = velocityRatio_re - targetRatio;
    
    if (velocityRatio_re >= 1.0) {
        entrainmentStatus = 'unqualified';
        entrainmentStatusText = '危险：水平气速超过防夹带速度，会发生液滴夹带';
    } else if (velocityRatio_re >= 0.9) {
        entrainmentStatus = 'unqualified';
        entrainmentStatusText = '危险：水平气速接近防夹带速度，存在夹带风险';
    } else if (velocityRatio_re >= 0.8) {
        entrainmentStatus = 'warning';
        if (targetVelocityRatio !== null && velocityRatio_re > targetRatio) {
            entrainmentStatusText = `警告：水平气速超过目标值 ${targetRatio.toFixed(2)}，当前 ${velocityRatio_re.toFixed(3)}，建议降低`;
        } else {
            entrainmentStatusText = '警告：水平气速接近防夹带速度，建议降低';
        }
    } else if (velocityRatio_re >= 0.6) {
        entrainmentStatus = 'qualified';
        if (targetVelocityRatio !== null) {
            if (velocityRatio_re <= targetRatio) {
                entrainmentStatusText = `良好：水平气速满足目标值 ${targetRatio.toFixed(2)}（当前：${velocityRatio_re.toFixed(3)}）`;
            } else {
                entrainmentStatusText = `注意：水平气速略高于目标值 ${targetRatio.toFixed(2)}（当前：${velocityRatio_re.toFixed(3)}）`;
                entrainmentStatus = 'warning';
            }
        } else {
            entrainmentStatusText = '良好：水平气速在安全范围内';
        }
    } else {
        entrainmentStatus = 'qualified';
        if (targetVelocityRatio !== null && velocityRatio_re <= targetRatio) {
            entrainmentStatusText = `优秀：水平气速满足目标值 ${targetRatio.toFixed(2)}，有充足安全余量（当前：${velocityRatio_re.toFixed(3)}）`;
        } else {
            entrainmentStatusText = '优秀：水平气速安全，有充足余量';
        }
    }
    
    evaluation.items.push({
        name: '防夹带速度比值 (vh/Vre)',
        value: velocityRatio_re.toFixed(3),
        status: entrainmentStatus,
        description: entrainmentStatusText
    });
    
    if (entrainmentStatus === 'unqualified') {
        evaluation.overallStatus = 'unqualified';
    } else if (entrainmentStatus === 'warning' && evaluation.overallStatus === 'qualified') {
        evaluation.overallStatus = 'warning';
    }

    // 2. 水平气速倍数评估 (vh/vt)
    const vh_vt_ratio = vh / vt;
    let vhVtStatus = 'qualified';
    let vhVtStatusText = '';
    
    if (vh_vt_ratio > 5.0) {
        vhVtStatus = 'warning';
        vhVtStatusText = '注意：水平气速倍数较高，建议控制在 1.0-5.0 之间';
    } else if (vh_vt_ratio >= 1.0 && vh_vt_ratio <= 5.0) {
        vhVtStatus = 'qualified';
        vhVtStatusText = '合理：水平气速倍数在推荐范围内';
    } else {
        vhVtStatus = 'warning';
        vhVtStatusText = '注意：水平气速倍数较低，设计较保守';
    }
    
    evaluation.items.push({
        name: '水平气速倍数 (vh/vt)',
        value: vh_vt_ratio.toFixed(2),
        status: vhVtStatus,
        description: vhVtStatusText
    });

    // 3. 液滴沉降能力评估
    let settleStatus = 'qualified';
    let settleStatusText = '';
    
    if (cannotSettle) {
        settleStatus = 'unqualified';
        settleStatusText = '不合格：液滴无法在有效长度内沉降到液面';
    } else {
        settleStatus = 'qualified';
        settleStatusText = '合格：液滴可以在有效长度内沉降到液面';
    }
    
    evaluation.items.push({
        name: '液滴沉降能力',
        value: cannotSettle ? '无法沉降' : '可以沉降',
        status: settleStatus,
        description: settleStatusText
    });
    
    if (settleStatus === 'unqualified') {
        evaluation.overallStatus = 'unqualified';
    }

    // 4. 直径标准化检查
    const standardDiameters = [100, 150, 200, 250, 300, 400, 500, 600, 800, 1000];
    let closestStandard = null;
    let minDiff = Infinity;
    
    for (const stdD of standardDiameters) {
        const diff = Math.abs(D_mm - stdD);
        if (diff < minDiff) {
            minDiff = diff;
            closestStandard = stdD;
        }
    }
    
    const diameterDeviation = (minDiff / closestStandard) * 100;
    let diameterStatus = 'qualified';
    let diameterStatusText = '';
    
    if (diameterDeviation <= 5) {
        diameterStatus = 'qualified';
        diameterStatusText = `标准直径：${closestStandard} mm（偏差 ${diameterDeviation.toFixed(1)}%）`;
    } else if (diameterDeviation <= 10) {
        diameterStatus = 'warning';
        diameterStatusText = `接近标准直径：${closestStandard} mm（偏差 ${diameterDeviation.toFixed(1)}%）`;
    } else {
        diameterStatus = 'warning';
        diameterStatusText = `非标准直径，建议使用：${closestStandard} mm`;
    }
    
    evaluation.items.push({
        name: '直径标准化',
        value: `${D_mm.toFixed(0)} mm`,
        status: diameterStatus,
        description: diameterStatusText,
        suggestedValue: diameterDeviation > 5 ? `${closestStandard} mm` : null
    });

    // 5. 长径比评估
    let ldRatioStatus = 'qualified';
    let ldRatioStatusText = '';
    
    if (L_D_ratio >= 2.0 && L_D_ratio <= 6.0) {
        if (L_D_ratio >= 2.5 && L_D_ratio <= 4.0) {
            ldRatioStatus = 'qualified';
            ldRatioStatusText = '推荐：长径比在最佳范围内';
        } else {
            ldRatioStatus = 'qualified';
            ldRatioStatusText = '可接受：长径比在允许范围内';
        }
    } else if (L_D_ratio < 2.0) {
        ldRatioStatus = 'warning';
        ldRatioStatusText = '不推荐：长径比过小，建议 L/D ≥ 2.0';
    } else {
        ldRatioStatus = 'warning';
        ldRatioStatusText = '不推荐：长径比过大，建议 L/D ≤ 6.0';
    }
    
    evaluation.items.push({
        name: '长径比 (L/D)',
        value: L_D_ratio.toFixed(2),
        status: ldRatioStatus,
        description: ldRatioStatusText
    });

    // 6. 液位高度比评估
    let liquidLevelStatus = 'qualified';
    let liquidLevelStatusText = '';
    
    if (h_liq_ratio >= 0.1 && h_liq_ratio <= 0.5) {
        if (h_liq_ratio >= 0.2 && h_liq_ratio <= 0.4) {
            liquidLevelStatus = 'qualified';
            liquidLevelStatusText = '推荐：液位高度比在最佳范围内';
        } else {
            liquidLevelStatus = 'qualified';
            liquidLevelStatusText = '可接受：液位高度比在允许范围内';
        }
    } else {
        liquidLevelStatus = 'warning';
        liquidLevelStatusText = '不推荐：液位高度比超出推荐范围 (0.1-0.5)';
    }
    
    evaluation.items.push({
        name: '液位高度比 (H_liq/D)',
        value: h_liq_ratio.toFixed(2),
        status: liquidLevelStatus,
        description: liquidLevelStatusText
    });

    // 生成具体建议
    if (velocityRatio_re >= 0.8 || (targetVelocityRatio !== null && velocityRatio_re > targetVelocityRatio)) {
        const targetRatio = targetVelocityRatio !== null ? targetVelocityRatio : 0.75;
        const suggested_vh = v_max_entrainment * targetRatio; // 使用目标比值或75%的安全余量
        
        if (calcMode === 'auto') {
            // 自动模式：告诉用户调整水平气速倍数输入框
            const current_vh_vt_ratio = vh / vt;
            const suggested_vh_vt_ratio = suggested_vh / vt;
            evaluation.suggestions.push(`【调整方法】在左侧输入参数区域，找到"水平气速倍数 (vh/vt)"输入框，将其值调整为 ${suggested_vh_vt_ratio.toFixed(2)} 以下（当前值：${current_vh_vt_ratio.toFixed(2)}），然后重新计算`);
            evaluation.suggestions.push(`调整后目标水平气速：${suggested_vh.toFixed(3)} m/s（当前：${vh.toFixed(3)} m/s）`);
            evaluation.suggestions.push(`【替代方案】或增大长径比 L/D 以增加分离长度，或降低液位高度比 H_liq/D 以增加气相流通面积`);
        } else {
            // 手动直径模式（校核模式）：建议增大直径
            if (targetVelocityRatio !== null) {
                evaluation.suggestions.push(`当前防夹带速度比值 ${velocityRatio_re.toFixed(3)} 超过目标值 ${targetVelocityRatio.toFixed(2)}，建议增大直径以降低实际气速`);
                // 计算气相流通面积比例（使用与计算函数相同的逻辑）
                // 简化估算：使用当前液位高度比计算气相面积比例
                const d_over_R = 1 - 2 * h_liq_ratio;
                const theta_half = Math.acos(Math.max(-1, Math.min(1, d_over_R)));
                const theta = 2 * theta_half;
                const liquidAreaRatio = (theta - Math.sin(theta)) / (2 * Math.PI);
                const gasAreaRatio = Math.max(0, Math.min(1, 1 - liquidAreaRatio));
                
                const V_g = m_dot / rho_g;
                const suggested_vh_target = v_max_entrainment * targetVelocityRatio;
                const suggested_D = Math.sqrt((4 * V_g) / (Math.PI * suggested_vh_target * gasAreaRatio)) * 1000;
                evaluation.suggestions.push(`建议直径：${suggested_D.toFixed(0)} mm 以上（当前：${D_mm.toFixed(0)} mm），以达到目标比值 ${targetVelocityRatio.toFixed(2)}`);
            } else {
                evaluation.suggestions.push(`降低水平气速至 ${suggested_vh.toFixed(3)} m/s 以下（当前：${vh.toFixed(3)} m/s）`);
                evaluation.suggestions.push(`增大直径以降低实际气速`);
            }
        }
    }
    
    if (cannotSettle) {
        if (calcMode === 'auto') {
            const current_vh_vt_ratio = vh / vt;
            const suggested_vh_vt_ratio_for_settle = (vt * 0.8) / vt; // 降低气速以增加沉降时间
            evaluation.suggestions.push(`降低"水平气速倍数 (vh/vt)"至 ${suggested_vh_vt_ratio_for_settle.toFixed(2)} 以下（当前：${current_vh_vt_ratio.toFixed(2)}），或增大长径比 L/D，确保液滴能在有效长度内沉降`);
        } else {
            evaluation.suggestions.push(`增大长度或降低水平气速，确保液滴能在有效长度内沉降`);
        }
    }
    
    if (diameterDeviation > 5) {
        evaluation.suggestions.push(`使用标准直径 ${closestStandard} mm（当前：${D_mm.toFixed(0)} mm）`);
    }
    
    if (L_D_ratio < 2.0 || L_D_ratio > 6.0) {
        const suggestedL = D_mm * 3.0;
        evaluation.suggestions.push(`调整长度至 ${suggestedL.toFixed(0)} mm（当前：${L_mm.toFixed(0)} mm，L/D = ${L_D_ratio.toFixed(2)}）`);
    }

    return evaluation;
}

// 评估水蒸汽立式分离器设计
// 参数：v_max (最大允许气速), D_mm (直径mm), H_mm (高度mm), H_D_ratio (高度/直径比), 
//      A_min (最小流通面积m²), efficiency (分离效率%), K (Souders-Brown系数)
function evaluateSteamDesign(v_max, D_mm, H_mm, H_D_ratio, A_min, efficiency, K) {
    const evaluation = {
        overallStatus: 'qualified',
        items: [],
        suggestions: []
    };

    // 1. 分离效率评估
    let efficiencyStatus = 'qualified';
    let efficiencyStatusText = '';
    
    if (efficiency >= 99.0) {
        efficiencyStatus = 'qualified';
        efficiencyStatusText = '优秀：分离效率很高';
    } else if (efficiency >= 95.0) {
        efficiencyStatus = 'qualified';
        efficiencyStatusText = '良好：分离效率满足要求';
    } else if (efficiency >= 90.0) {
        efficiencyStatus = 'warning';
        efficiencyStatusText = '可接受：分离效率较低，建议提高至 95% 以上';
    } else {
        efficiencyStatus = 'warning';
        efficiencyStatusText = '不推荐：分离效率过低';
    }
    
    evaluation.items.push({
        name: '分离效率',
        value: `${efficiency.toFixed(2)}%`,
        status: efficiencyStatus,
        description: efficiencyStatusText
    });

    // 2. Souders-Brown 系数评估
    let kStatus = 'qualified';
    let kStatusText = '';
    
    if (K >= 0.10 && K <= 0.20) {
        kStatus = 'qualified';
        kStatusText = '推荐：K值在标准范围内（适用于丝网除沫器）';
    } else if (K >= 0.05 && K <= 0.30) {
        kStatus = 'qualified';
        kStatusText = '可接受：K值在允许范围内';
    } else {
        kStatus = 'warning';
        kStatusText = '注意：K值超出常规范围，请确认适用性';
    }
    
    evaluation.items.push({
        name: 'Souders-Brown 系数 K',
        value: K.toFixed(3),
        status: kStatus,
        description: kStatusText
    });

    // 3. 直径标准化检查
    const standardDiameters = [100, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 1200, 1500];
    let closestStandard = null;
    let minDiff = Infinity;
    
    for (const stdD of standardDiameters) {
        const diff = Math.abs(D_mm - stdD);
        if (diff < minDiff) {
            minDiff = diff;
            closestStandard = stdD;
        }
    }
    
    const diameterDeviation = (minDiff / closestStandard) * 100;
    let diameterStatus = 'qualified';
    let diameterStatusText = '';
    
    if (diameterDeviation <= 5) {
        diameterStatus = 'qualified';
        diameterStatusText = `标准直径：${closestStandard} mm（偏差 ${diameterDeviation.toFixed(1)}%）`;
    } else if (diameterDeviation <= 10) {
        diameterStatus = 'warning';
        diameterStatusText = `接近标准直径：${closestStandard} mm（偏差 ${diameterDeviation.toFixed(1)}%）`;
    } else {
        diameterStatus = 'warning';
        diameterStatusText = `非标准直径，建议使用：${closestStandard} mm`;
    }
    
    evaluation.items.push({
        name: '直径标准化',
        value: `${D_mm.toFixed(0)} mm`,
        status: diameterStatus,
        description: diameterStatusText,
        suggestedValue: diameterDeviation > 5 ? `${closestStandard} mm` : null
    });

    // 4. 高度范围验证
    let heightStatus = 'qualified';
    let heightStatusText = '';
    
    if (H_D_ratio >= 2.5 && H_D_ratio <= 3.0) {
        heightStatus = 'qualified';
        heightStatusText = '推荐：高度在最佳范围内';
    } else if (H_D_ratio >= 2.0 && H_D_ratio < 2.5) {
        heightStatus = 'warning';
        heightStatusText = '可接受：高度略低，建议 H/D = 2.5-3.0';
    } else if (H_D_ratio > 3.0 && H_D_ratio <= 3.5) {
        heightStatus = 'warning';
        heightStatusText = '可接受：高度略高，建议 H/D = 2.5-3.0';
    } else if (H_D_ratio < 2.0) {
        heightStatus = 'warning';
        heightStatusText = '不推荐：高度过低，建议 H/D ≥ 2.5';
    } else {
        heightStatus = 'warning';
        heightStatusText = '不推荐：高度过高，建议 H/D ≤ 3.5';
    }
    
    const suggestedHeight = D_mm * 2.75;
    
    evaluation.items.push({
        name: '高度范围 (H/D)',
        value: H_D_ratio.toFixed(2),
        status: heightStatus,
        description: heightStatusText,
        suggestedValue: (H_D_ratio < 2.5 || H_D_ratio > 3.0) ? `${suggestedHeight.toFixed(0)} mm` : null
    });

    // 5. 最大允许气速合理性
    let velocityStatus = 'qualified';
    let velocityStatusText = '';
    
    if (v_max > 0.5 && v_max < 5.0) {
        velocityStatus = 'qualified';
        velocityStatusText = '合理：最大允许气速在典型范围内';
    } else if (v_max <= 0.5) {
        velocityStatus = 'warning';
        velocityStatusText = '注意：最大允许气速较低，可能需要较大直径';
    } else {
        velocityStatus = 'warning';
        velocityStatusText = '注意：最大允许气速较高，请确认计算正确性';
    }
    
    evaluation.items.push({
        name: '最大允许气速',
        value: `${v_max.toFixed(3)} m/s`,
        status: velocityStatus,
        description: velocityStatusText
    });

    // 生成具体建议
    if (efficiency < 95.0) {
        evaluation.suggestions.push(`提高分离效率至 95% 以上（当前：${efficiency.toFixed(2)}%）`);
    }
    
    if (diameterDeviation > 5) {
        evaluation.suggestions.push(`使用标准直径 ${closestStandard} mm（当前：${D_mm.toFixed(0)} mm）`);
    }
    
    if (H_D_ratio < 2.5 || H_D_ratio > 3.0) {
        evaluation.suggestions.push(`调整高度至 ${suggestedHeight.toFixed(0)} mm（当前：${H_mm.toFixed(0)} mm，H/D = ${H_D_ratio.toFixed(2)}）`);
    }
    
    if (K < 0.10 || K > 0.20) {
        evaluation.suggestions.push(`建议使用标准 K 值范围 0.10-0.20（当前：${K.toFixed(3)}），适用于丝网除沫器`);
    }

    return evaluation;
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
        // 如果过热度为 0，使用 0.001 进行计算
        const deltaT_sh_calc = deltaT_sh === 0 ? 0.001 : deltaT_sh;
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
        const T_suc_K = Te_K + deltaT_sh_calc;

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

        // 13. 计算液体停留时间
        // 立式分离器：假设液体占底部10-20%的高度（取15%作为典型值）
        const liquid_height_ratio = 0.15; // 液体高度占筒身高度的比例
        const liquid_volume_m3 = Vol_m3 * liquid_height_ratio;
        // 估算液体流量：假设液体流量为总质量流量的10%（粗略估算）
        const liquid_mass_flow_kg_s = m_dot * 0.1;
        const liquid_flow_m3_s = liquid_mass_flow_kg_s / rho_l;
        const residence_time = liquid_volume_m3 / Math.max(liquid_flow_m3_s, 1e-6);

        // 更新结果显示
        document.getElementById('massFlow').textContent = formatNumber(m_dot, 3);
        document.getElementById('volumeFlow').textContent = formatNumber(V_g, 6);
        document.getElementById('gasDensity').textContent = formatNumber(rho_g, 2);
        document.getElementById('diameter').textContent = formatNumber(D_mm, 0);
        document.getElementById('height').textContent = formatNumber(H_mm, 0);
        document.getElementById('volume').textContent = formatNumber(Vol_L, 1);
        document.getElementById('actualVelocity').textContent = formatNumber(u, 3);
        document.getElementById('terminalVelocity').textContent = formatNumber(vt, 3);
        document.getElementById('liquidDensity').textContent = formatNumber(rho_l, 2);
        document.getElementById('gasViscosity').textContent = formatNumber(eta_g, 6);
        document.getElementById('residenceTime').textContent = formatNumber(residence_time, 1);

        // 评估设计并显示结论
        const currentVelocityRatio = calcMode === 'auto' ? velocityRatio : null;
        const evaluation = evaluateVerticalDesign(u, vt, D_mm, H_mm, Vol_L, m_dot, rho_g, calcMode, currentVelocityRatio);
        const conclusionHTML = generateConclusionHTML(evaluation);
        
        const conclusionSection = document.getElementById('conclusionSection');
        const conclusionContent = document.getElementById('conclusionContent');
        if (conclusionSection && conclusionContent) {
            conclusionContent.innerHTML = conclusionHTML;
            conclusionSection.style.display = 'block';
        }

        // 根据评估结果显示状态消息（不合理的设计不会自动隐藏）
        if (calcMode === 'manual-diameter') {
            const actualRatio = u / vt;
            if (evaluation.overallStatus === 'unqualified') {
                showStatus(`❌ 校核结果：设计不合格！实际气速: ${u.toFixed(3)} m/s，气速/终端速度比值: ${actualRatio.toFixed(2)}，请查看下方结论与建议进行调整。`, 'error');
            } else if (evaluation.overallStatus === 'warning') {
                showStatus(`⚠️ 校核结果：设计需注意！实际气速: ${u.toFixed(3)} m/s，气速/终端速度比值: ${actualRatio.toFixed(2)}，请查看下方结论与建议。`, 'error');
            } else {
                showStatus(`✓ 校核完成！实际气速: ${u.toFixed(3)} m/s，气速/终端速度比值: ${actualRatio.toFixed(2)}，设计合格。`, 'success');
            }
        } else {
            if (evaluation.overallStatus === 'unqualified') {
                showStatus('❌ 计算完成，但设计不合格！请查看下方结论与建议进行调整。', 'error');
            } else if (evaluation.overallStatus === 'warning') {
                showStatus('⚠️ 计算完成，但设计需注意！请查看下方结论与建议。', 'error');
            } else {
                showStatus('✓ 计算完成！设计合格。', 'success');
            }
        }

    } catch (error) {
        console.error('计算错误:', error);
        showStatus('错误: ' + error.message, 'error');
        
        // 隐藏结论部分
        const conclusionSection = document.getElementById('conclusionSection');
        if (conclusionSection) {
            conclusionSection.style.display = 'none';
        }
        
        // 清空结果显示
        document.getElementById('massFlow').textContent = '-';
        document.getElementById('volumeFlow').textContent = '-';
        document.getElementById('gasDensity').textContent = '-';
        document.getElementById('diameter').textContent = '-';
        document.getElementById('height').textContent = '-';
        document.getElementById('volume').textContent = '-';
        document.getElementById('actualVelocity').textContent = '-';
        document.getElementById('terminalVelocity').textContent = '-';
        document.getElementById('liquidDensity').textContent = '-';
        document.getElementById('gasViscosity').textContent = '-';
        document.getElementById('residenceTime').textContent = '-';
    } finally {
        // 恢复计算按钮
        const calcBtn = document.getElementById('calculateBtn');
        calcBtn.disabled = false;
        calcBtn.textContent = '计算';
    }
}

// 获取K值对应的不锈钢丝网规格信息
function getMeshSpecByK(K) {
    const meshSpecs = {
        '0.05-0.08': {
            range: '0.05 - 0.08',
            mesh: '100目/英寸',
            wireDiameter: '0.10 mm',
            thickness: '50-100 mm',
            material: '304/316不锈钢',
            application: '低负荷，精细分离',
            porosity: '约85%'
        },
        '0.08-0.12': {
            range: '0.08 - 0.12',
            mesh: '80目/英寸',
            wireDiameter: '0.12 mm',
            thickness: '100-150 mm',
            material: '304/316不锈钢',
            application: '中等负荷，标准分离',
            porosity: '约88%'
        },
        '0.12-0.15': {
            range: '0.12 - 0.15',
            mesh: '60目/英寸',
            wireDiameter: '0.15 mm',
            thickness: '100-150 mm',
            material: '304/316不锈钢',
            application: '标准应用，推荐值',
            porosity: '约90%'
        },
        '0.15-0.18': {
            range: '0.15 - 0.18',
            mesh: '50目/英寸',
            wireDiameter: '0.18 mm',
            thickness: '100-200 mm',
            material: '304/316不锈钢',
            application: '高负荷，粗分离',
            porosity: '约92%'
        },
        '0.18-0.25': {
            range: '0.18 - 0.25',
            mesh: '40目/英寸',
            wireDiameter: '0.20 mm',
            thickness: '150-200 mm',
            material: '304/316不锈钢',
            application: '极高负荷，粗分离',
            porosity: '约94%'
        }
    };
    
    if (K <= 0.08) {
        return meshSpecs['0.05-0.08'];
    } else if (K <= 0.12) {
        return meshSpecs['0.08-0.12'];
    } else if (K <= 0.15) {
        return meshSpecs['0.12-0.15'];
    } else if (K <= 0.18) {
        return meshSpecs['0.15-0.18'];
    } else {
        return meshSpecs['0.18-0.25'];
    }
}

// 更新K值对应的丝网规格提示
function updateMeshSpecInfo(K) {
    const specInfo = document.getElementById('meshSpecInfo');
    const specDetails = document.getElementById('meshSpecDetails');
    
    if (!specInfo || !specDetails) return;
    
    if (K >= 0.05 && K <= 0.3) {
        const spec = getMeshSpecByK(K);
        specDetails.innerHTML = `
            <div><strong>K值范围：</strong>${spec.range}</div>
            <div><strong>丝网规格：</strong>${spec.mesh}</div>
            <div><strong>丝径：</strong>${spec.wireDiameter}</div>
            <div><strong>推荐厚度：</strong>${spec.thickness}</div>
            <div><strong>材质：</strong>${spec.material}</div>
            <div><strong>孔隙率：</strong>${spec.porosity}</div>
            <div style="margin-top: 8px; color: #1976d2;"><strong>应用：</strong>${spec.application}</div>
        `;
        specInfo.style.display = 'block';
    } else {
        specInfo.style.display = 'none';
    }
}

// 水蒸汽立式计算函数
async function calculateSteam() {
    try {
        // 检查 CoolProp 是否已初始化
        if (!CoolPropModule) {
            showStatus('正在初始化 CoolProp...', 'info', 'statusSteam');
            await initCoolProp();
        }

        // 获取输入参数
        const P_in_MPa = parseFloat(document.getElementById('inletPressure').value);
        const x_in = parseFloat(document.getElementById('inletDryness').value);
        const m_total_kg_h = parseFloat(document.getElementById('massFlowRate').value);
        const efficiency = parseFloat(document.getElementById('separationEfficiency').value);
        const K = parseFloat(document.getElementById('kFactor').value);

        // 验证输入
        if (isNaN(P_in_MPa) || P_in_MPa <= 0 || P_in_MPa > 22.1) {
            throw new Error('入口压力必须在 0.1 - 22.1 MPa 之间（水的临界压力为 22.1 MPa）');
        }
        if (isNaN(x_in) || x_in < 0 || x_in > 1) {
            throw new Error('入口干度必须在 0.0 - 1.0 之间');
        }
        if (isNaN(m_total_kg_h) || m_total_kg_h <= 0) {
            throw new Error('总质量流量必须大于 0');
        }
        if (isNaN(efficiency) || efficiency < 90 || efficiency > 99.99) {
            throw new Error('分离效率必须在 90% - 99.99% 之间');
        }
        if (isNaN(K) || K <= 0 || K > 0.5) {
            throw new Error('Souders-Brown 系数必须在 0.05 - 0.5 之间');
        }

        // 禁用计算按钮
        const calcBtn = document.getElementById('calculateBtnSteam');
        calcBtn.disabled = true;
        calcBtn.textContent = '计算中...';

        showStatus('正在计算...', 'info', 'statusSteam');

        // 单位转换：MPa 转 Pa
        const P_in_Pa = P_in_MPa * 1e6;

        // 使用 CoolProp 计算水/蒸汽属性
        // CoolProp 中水的名称是 "Water" 或 "H2O"
        const fluid = 'Water';

        // 1. 计算饱和温度（基于压力）
        const T_sat_K = CoolPropModule.PropsSI('T', 'P', P_in_Pa, 'Q', 1, fluid);
        if (T_sat_K <= 0 || !isFinite(T_sat_K)) {
            throw new Error('无法计算饱和温度，请检查压力参数');
        }
        const T_sat_C = T_sat_K - 273.15;

        // 2. 计算饱和液体密度
        const rho_l = CoolPropModule.PropsSI('D', 'P', P_in_Pa, 'Q', 0, fluid);
        if (rho_l <= 0 || !isFinite(rho_l)) {
            throw new Error('无法计算液体密度');
        }

        // 3. 计算饱和蒸汽密度
        const rho_v = CoolPropModule.PropsSI('D', 'P', P_in_Pa, 'Q', 1, fluid);
        if (rho_v <= 0 || !isFinite(rho_v)) {
            throw new Error('无法计算蒸汽密度');
        }

        // 4. 质量平衡计算
        // m_vapor = m_total * x_in (入口蒸汽质量流量)
        // m_liquid = m_total * (1 - x_in) (入口液体质量流量)
        const m_vapor_in_kg_h = m_total_kg_h * x_in;
        const m_liquid_in_kg_h = m_total_kg_h * (1 - x_in);

        // 4.1 计算进口体积流量 (m³/s)
        // 总进口体积流量 = 蒸汽体积流量 + 液体体积流量
        const m_vapor_in_kg_s = m_vapor_in_kg_h / 3600;
        const m_liquid_in_kg_s = m_liquid_in_kg_h / 3600;
        const V_vapor_in = m_vapor_in_kg_s / rho_v; // 蒸汽体积流量
        const V_liquid_in = m_liquid_in_kg_s / rho_l; // 液体体积流量
        const V_total_in = V_vapor_in + V_liquid_in; // 总进口体积流量

        // 5. 分离效率计算
        // 假设分离效率应用于液体分离（即有多少液体被分离出来）
        // 分离后的液体流量 = 入口液体流量 + 部分蒸汽冷凝（根据效率）
        // 简化模型：假设分离效率主要影响液体分离
        const m_liquid_separated_kg_h = m_liquid_in_kg_h + m_vapor_in_kg_h * (1 - efficiency / 100);
        const m_steam_out_kg_h = m_total_kg_h - m_liquid_separated_kg_h;

        // 6. Souders-Brown 方程计算
        // v_max = K * sqrt((rho_l - rho_v) / rho_v)
        const v_max = K * Math.sqrt((rho_l - rho_v) / rho_v);
        if (v_max <= 0 || !isFinite(v_max)) {
            throw new Error('无法计算最大允许气速');
        }

        // 7. 计算最小流通面积
        // 基于出口蒸汽流量计算
        // A_min = (m_steam / 3600) / (v_max * rho_v)
        // 注意：质量流量单位从 kg/h 转换为 kg/s
        const m_steam_kg_s = m_steam_out_kg_h / 3600;
        const A_min = m_steam_kg_s / (v_max * rho_v);
        if (A_min <= 0 || !isFinite(A_min)) {
            throw new Error('无法计算最小流通面积');
        }

        // 8. 计算推荐容器直径
        // D = sqrt(4 * A_min / PI)
        const D_m = Math.sqrt((4 * A_min) / Math.PI);
        const D_mm = D_m * 1000;

        // 8.1 计算推荐气分高度（通常为2.5-3倍直径）
        const H_D_ratio = 2.75; // 推荐长径比
        const H_mm = D_mm * H_D_ratio;

        // 8.2 计算实际气速（基于出口蒸汽流量和容器直径）
        // 实际气速 = 出口蒸汽体积流量 / 容器截面积
        const V_steam_out = m_steam_kg_s / rho_v; // 出口蒸汽体积流量 (m³/s)
        const A_vessel = Math.PI * Math.pow(D_m / 2, 2); // 容器截面积 (m²)
        const v_actual = V_steam_out / A_vessel; // 实际气速 (m/s)

        // 8.3 计算液体停留时间
        // 计算容器总容积
        const Vol_m3 = Math.PI * Math.pow(D_m / 2, 2) * (H_mm / 1000);
        // 假设液体占底部10-20%的高度（取15%作为典型值）
        const liquid_height_ratio = 0.15;
        const liquid_volume_m3 = Vol_m3 * liquid_height_ratio;
        // 计算液体流量（分离后的液体质量流量转换为体积流量）
        const m_liquid_separated_kg_s = m_liquid_separated_kg_h / 3600;
        const liquid_flow_m3_s = m_liquid_separated_kg_s / rho_l;
        const residence_time = liquid_volume_m3 / Math.max(liquid_flow_m3_s, 1e-6);

        // 9. 更新结果显示
        document.getElementById('saturationTemp').textContent = formatNumber(T_sat_C, 2);
        document.getElementById('volumeFlowSteam').textContent = formatNumber(V_total_in, 6);
        document.getElementById('liquidDensitySteam').textContent = formatNumber(rho_l, 2);
        document.getElementById('vaporDensitySteam').textContent = formatNumber(rho_v, 4);
        document.getElementById('separatedLiquidFlow').textContent = formatNumber(m_liquid_separated_kg_h, 2);
        document.getElementById('drySteamFlow').textContent = formatNumber(m_steam_out_kg_h, 2);
        document.getElementById('actualVelocitySteam').textContent = formatNumber(v_actual, 3);
        document.getElementById('maxVelocity').textContent = formatNumber(v_max, 3);
        document.getElementById('vesselDiameter').textContent = formatNumber(D_mm, 0);
        const heightElement = document.getElementById('vesselHeight');
        if (heightElement) {
            heightElement.textContent = formatNumber(H_mm, 0);
        }
        document.getElementById('minArea').textContent = formatNumber(A_min, 4);
        document.getElementById('residenceTimeSteam').textContent = formatNumber(residence_time, 1);

        // 10. 更新 SVG 示意图
        updateSteamDiagram(D_mm, H_mm, x_in, m_liquid_separated_kg_h / m_total_kg_h);

        // 11. 评估设计并显示结论
        const evaluationSteam = evaluateSteamDesign(v_max, D_mm, H_mm, H_D_ratio, A_min, efficiency, K);
        const conclusionHTMLSteam = generateConclusionHTML(evaluationSteam);
        
        const conclusionSectionSteam = document.getElementById('conclusionSectionSteam');
        const conclusionContentSteam = document.getElementById('conclusionContentSteam');
        if (conclusionSectionSteam && conclusionContentSteam) {
            conclusionContentSteam.innerHTML = conclusionHTMLSteam;
            conclusionSectionSteam.style.display = 'block';
        }

        // 12. 根据评估结果显示状态消息（不合理的设计不会自动隐藏）
        if (evaluationSteam.overallStatus === 'unqualified') {
            showStatus('❌ 计算完成，但设计不合格！请查看下方结论与建议进行调整。', 'error', 'statusSteam');
        } else if (evaluationSteam.overallStatus === 'warning') {
            showStatus('⚠️ 计算完成，但设计需注意！请查看下方结论与建议。', 'error', 'statusSteam');
        } else {
            showStatus(`✓ 计算完成！推荐容器直径: ${formatNumber(D_mm, 0)} mm，最大允许气速: ${formatNumber(v_max, 3)} m/s，设计合格。`, 'success', 'statusSteam');
        }

    } catch (error) {
        console.error('计算错误:', error);
        showStatus('错误: ' + error.message, 'error', 'statusSteam');
        
        // 隐藏结论部分
        const conclusionSectionSteam = document.getElementById('conclusionSectionSteam');
        if (conclusionSectionSteam) {
            conclusionSectionSteam.style.display = 'none';
        }
        
        // 清空结果显示
        document.getElementById('saturationTemp').textContent = '-';
        document.getElementById('volumeFlowSteam').textContent = '-';
        document.getElementById('liquidDensitySteam').textContent = '-';
        document.getElementById('vaporDensitySteam').textContent = '-';
        document.getElementById('separatedLiquidFlow').textContent = '-';
        document.getElementById('drySteamFlow').textContent = '-';
        document.getElementById('actualVelocitySteam').textContent = '-';
        document.getElementById('maxVelocity').textContent = '-';
        document.getElementById('vesselDiameter').textContent = '-';
        const heightElement = document.getElementById('vesselHeight');
        if (heightElement) {
            heightElement.textContent = '-';
        }
        document.getElementById('minArea').textContent = '-';
        const residenceTimeSteamElement = document.getElementById('residenceTimeSteam');
        if (residenceTimeSteamElement) {
            residenceTimeSteamElement.textContent = '-';
        }
    } finally {
        // 恢复计算按钮
        const calcBtn = document.getElementById('calculateBtnSteam');
        calcBtn.disabled = false;
        calcBtn.textContent = '计算';
    }
}

// 更新水蒸汽立式 SVG 示意图
function updateSteamDiagram(diameter_mm, height_mm, inletDryness, liquidRatio) {
    const svg = document.getElementById('steamDiagram');
    if (!svg) return;

    const svgWidth = 400;
    const svgHeight = 500;
    const tankX = 150;
    const tankY = 50;
    const tankWidth = 100;
    const tankHeight = 400;

    // 根据液体比例计算液位高度
    // liquidRatio: 分离后的液体占总流量的比例
    const liquidHeight = tankHeight * Math.min(1.0, Math.max(0.0, liquidRatio));
    const liquidTopY = tankY + tankHeight - liquidHeight;

    // 更新容器主体（如果需要调整尺寸）
    const tankBody = document.getElementById('steamTankBody');
    if (tankBody) {
        // 保持容器尺寸不变，只更新位置
    }

    // 更新液体区域
    const liquidArea = document.getElementById('steamLiquidArea');
    if (liquidArea) {
        liquidArea.setAttribute('x', tankX);
        liquidArea.setAttribute('y', liquidTopY);
        liquidArea.setAttribute('width', tankWidth);
        liquidArea.setAttribute('height', liquidHeight);
        if (liquidHeight > 0) {
            liquidArea.setAttribute('rx', '5');
        } else {
            liquidArea.setAttribute('rx', '0');
        }
    }

    // 更新液位线
    const liquidLevel = document.getElementById('steamLiquidLevel');
    if (liquidLevel) {
        liquidLevel.setAttribute('x1', tankX);
        liquidLevel.setAttribute('y1', liquidTopY);
        liquidLevel.setAttribute('x2', tankX + tankWidth);
        liquidLevel.setAttribute('y2', liquidTopY);
        if (liquidHeight <= 0) {
            liquidLevel.style.display = 'none';
        } else {
            liquidLevel.style.display = 'block';
        }
    }

    // 更新直径标注（水平标注）
    const diameterLabel = document.getElementById('steamDiameterLabel');
    if (diameterLabel) {
        diameterLabel.textContent = `内径 D = ${formatNumber(diameter_mm, 0)} mm`;
    }
    
    // 更新气分高度标注
    const heightLabel = document.getElementById('steamHeightLabel');
    if (heightLabel) {
        heightLabel.textContent = `推荐高度 H = ${formatNumber(height_mm, 0)} mm`;
    }
}

// 页面加载完成后初始化
document.addEventListener('DOMContentLoaded', () => {
    // 类型切换功能
    const verticalBtn = document.getElementById('verticalBtn');
    const horizontalBtn = document.getElementById('horizontalBtn');
    const steamBtn = document.getElementById('steamBtn');
    const verticalContent = document.getElementById('verticalContent');
    const horizontalContent = document.getElementById('horizontalContent');
    const steamContent = document.getElementById('steamContent');

    verticalBtn.addEventListener('click', () => {
        verticalBtn.classList.add('active');
        horizontalBtn.classList.remove('active');
        steamBtn.classList.remove('active');
        if (verticalContent) {
            verticalContent.style.display = 'block';
            verticalContent.style.visibility = 'visible';
        }
        if (horizontalContent) {
            horizontalContent.style.display = 'none';
        }
        if (steamContent) {
            steamContent.style.display = 'none';
        }
    });

    horizontalBtn.addEventListener('click', () => {
        horizontalBtn.classList.add('active');
        verticalBtn.classList.remove('active');
        steamBtn.classList.remove('active');
        if (horizontalContent) {
            horizontalContent.style.display = 'block';
            horizontalContent.style.visibility = 'visible';
        }
        if (verticalContent) {
            verticalContent.style.display = 'none';
        }
        if (steamContent) {
            steamContent.style.display = 'none';
        }
    });

    steamBtn.addEventListener('click', () => {
        steamBtn.classList.add('active');
        verticalBtn.classList.remove('active');
        horizontalBtn.classList.remove('active');
        if (steamContent) {
            steamContent.style.display = 'block';
            steamContent.style.visibility = 'visible'; /* 确保可见 */
        }
        if (verticalContent) {
            verticalContent.style.display = 'none';
        }
        if (horizontalContent) {
            horizontalContent.style.display = 'none';
        }
    });

    // 绑定计算按钮事件
    const calcBtn = document.getElementById('calculateBtn');
    calcBtn.addEventListener('click', calculate);
    
    const calcBtnH = document.getElementById('calculateBtnH');
    calcBtnH.addEventListener('click', calculateHorizontal);
    
    const calcBtnSteam = document.getElementById('calculateBtnSteam');
    if (calcBtnSteam) {
        calcBtnSteam.addEventListener('click', calculateSteam);
    }

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
                } else if (steamContent && steamContent.style.display !== 'none') {
                    calculateSteam();
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
    
    const targetVelocityRatioGroupH = document.getElementById('targetVelocityRatioGroupH');
    
    calcModeSelectH.addEventListener('change', (e) => {
        if (e.target.value === 'auto') {
            // 自动模式：显示所有设计参数
            lengthDiameterRatioGroup.style.display = 'block';
            diameterInputGroupH.style.display = 'none';
            liquidLevelRatioGroup.style.display = 'block';
            ksValueGroup.style.display = 'block';
            velocityMultiplierGroup.style.display = 'block';
            if (targetVelocityRatioGroupH) {
                targetVelocityRatioGroupH.style.display = 'none';
            }
        } else if (e.target.value === 'manual-diameter') {
            // 手动输入直径模式：显示直径输入，保留长径比和液位比
            lengthDiameterRatioGroup.style.display = 'block';
            diameterInputGroupH.style.display = 'block';
            liquidLevelRatioGroup.style.display = 'block';
            ksValueGroup.style.display = 'block';
            velocityMultiplierGroup.style.display = 'none'; // 校核模式不需要输入 vh/vt 倍数
            if (targetVelocityRatioGroupH) {
                targetVelocityRatioGroupH.style.display = 'block'; // 校核模式显示目标防夹带速度比值
            }
        }
    });

    // 初始化显示状态（卧式）
    if (calcModeSelectH.value === 'auto') {
        lengthDiameterRatioGroup.style.display = 'block';
        diameterInputGroupH.style.display = 'none';
        liquidLevelRatioGroup.style.display = 'block';
        ksValueGroup.style.display = 'block';
        velocityMultiplierGroup.style.display = 'block';
        if (targetVelocityRatioGroupH) {
            targetVelocityRatioGroupH.style.display = 'none';
        }
    } else if (calcModeSelectH.value === 'manual-diameter') {
        lengthDiameterRatioGroup.style.display = 'block';
        diameterInputGroupH.style.display = 'block';
        liquidLevelRatioGroup.style.display = 'block';
        ksValueGroup.style.display = 'block';
        velocityMultiplierGroup.style.display = 'none';
        if (targetVelocityRatioGroupH) {
            targetVelocityRatioGroupH.style.display = 'block';
        }
    }

    // K值输入框变化时更新丝网规格提示
    const kFactorInput = document.getElementById('kFactor');
    if (kFactorInput) {
        // 初始化时显示默认K值的规格
        updateMeshSpecInfo(parseFloat(kFactorInput.value) || 0.15);
        
        // 监听输入变化
        kFactorInput.addEventListener('input', (e) => {
            const K = parseFloat(e.target.value);
            if (!isNaN(K)) {
                updateMeshSpecInfo(K);
            }
        });
        
        // 监听值变化（包括通过步进按钮改变）
        kFactorInput.addEventListener('change', (e) => {
            const K = parseFloat(e.target.value);
            if (!isNaN(K)) {
                updateMeshSpecInfo(K);
            }
        });
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

