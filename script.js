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
        const moduleConfig = {
            locateFile: (path) => {
                // 确保 wasm 文件路径正确
                if (path.endsWith('.wasm')) {
                    return './coolprop.wasm';
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
function showStatus(message, type = 'info') {
    const statusEl = document.getElementById('status');
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

// 主计算函数
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
        if (calcMode === 'manual' && (isNaN(u) || u <= 0)) {
            throw new Error('设计气速必须大于 0');
        }
        if (isNaN(Dp_um) || Dp_um <= 0) {
            throw new Error('液滴直径必须大于 0');
        }
        if (isNaN(velocityRatio) || velocityRatio <= 0 || velocityRatio > 1) {
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

        // 6. 根据计算模式确定设计气速
        if (calcMode === 'auto') {
            // 自动模式：根据终端速度和系数计算设计气速 (公式 8)
            u = velocityRatio * vt;
            // 更新输入框显示计算出的值
            document.getElementById('gasVelocity').value = u.toFixed(3);
        } else {
            // 手动模式：验证输入的气速是否合理
            if (u > vt) {
                throw new Error(`设计气速 ${u.toFixed(3)} m/s 大于终端速度 ${vt.toFixed(3)} m/s，无法有效分离液滴`);
            }
            const actualRatio = u / vt;
            if (actualRatio > 0.9) {
                showStatus(`警告: 设计气速/终端速度比值 ${actualRatio.toFixed(2)} 接近上限，建议小于 0.9`, 'info');
            }
        }

        // 7. 计算焓差 (J/kg)
        const delta_h = h_gas - h_liq;
        if (delta_h <= 0 || !isFinite(delta_h)) {
            throw new Error('焓差计算异常，请检查过热度设置');
        }

        // 8. 计算质量流量 (kg/s)
        // Q_kW * 1000 转换为 W，除以焓差得到质量流量
        const m_dot = (Q_kW * 1000) / delta_h;

        // 9. 计算气体体积流量 (m³/s)
        const V_g = m_dot / rho_g;

        // 10. 计算筒身直径 (m) - 根据文档公式 9
        // D = sqrt(4 * V_g / (π * u))
        const D_m = Math.sqrt((4 * V_g) / (Math.PI * u));

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

        showStatus('计算完成！', 'success');

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
    // 绑定计算按钮事件
    const calcBtn = document.getElementById('calculateBtn');
    calcBtn.addEventListener('click', calculate);

    // 绑定回车键事件（在输入框中按回车也可以计算）
    const inputFields = document.querySelectorAll('.input-field');
    inputFields.forEach(field => {
        field.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                calculate();
            }
        });
    });

    // 计算模式切换
    const calcModeSelect = document.getElementById('calcMode');
    const dropletSizeGroup = document.getElementById('dropletSizeGroup');
    const gasVelocityGroup = document.getElementById('gasVelocityGroup');
    
    calcModeSelect.addEventListener('change', (e) => {
        if (e.target.value === 'auto') {
            // 自动模式：显示液滴直径，隐藏手动气速输入
            dropletSizeGroup.style.display = 'block';
            gasVelocityGroup.style.display = 'none';
        } else {
            // 手动模式：显示手动气速输入，隐藏液滴直径
            dropletSizeGroup.style.display = 'none';
            gasVelocityGroup.style.display = 'block';
        }
    });

    // 初始化显示状态
    if (calcModeSelect.value === 'auto') {
        dropletSizeGroup.style.display = 'block';
        gasVelocityGroup.style.display = 'none';
    } else {
        dropletSizeGroup.style.display = 'none';
        gasVelocityGroup.style.display = 'block';
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

