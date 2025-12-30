# 立式气液分离器选型计算器

基于 Web 的立式气液分离器（Vertical Suction Accumulator）选型计算器，用于 HVAC 工程中的气液分离器设计计算。

## 功能特性

- ✅ 支持多种制冷剂：R134a, R245fa, R1234ze(E), R410A, Ammonia, R290
- ✅ 基于 CoolProp 的精确物性计算
- ✅ 自动计算液滴终端速度
- ✅ 两种计算模式：自动计算 / 手动输入
- ✅ 符合工程标准的计算公式
- ✅ 现代化的工业风格界面

## 技术栈

- **Vite** - 现代化构建工具
- HTML5
- CSS3 (Flexbox/Grid 布局)
- Vanilla JavaScript (ES6 模块)
- CoolProp.js (本地 WASM)

## 计算方法

### 核心公式

1. **终端速度计算**（迭代法）：
   - `vt² = 4 × g × Dp × (ρl - ρg) / (3 × Cd × ρg)`
   - `Cd = 24/Re + 3/√Re + 0.34`
   - `Re = Dp × vt × ρg / ηg`

2. **设计气速**：
   - `v = (0.75 to 0.9) × vt`

3. **筒身直径**：
   - `D = √(4 × V / (π × v))`

4. **筒身高度**：
   - `H = 3.0 × D` (推荐范围: 2.5-3 D)

## 使用方法

### 开发环境

```bash
# 安装依赖
npm install

# 启动开发服务器（支持热更新）
npm run dev

# 构建生产版本
npm run build

# 预览生产构建
npm run preview
```

开发服务器会自动在 http://localhost:3000 打开。

### 在线访问

访问 [Vercel 部署地址](https://separator-app.vercel.app)（部署后更新）

## 项目结构

```
SeparatorApp/
├── index.html          # 主页面（Vite 入口）
├── vite.config.js      # Vite 配置文件
├── package.json        # 项目配置和依赖
├── vercel.json         # Vercel 部署配置
├── src/
│   ├── script.js      # 计算逻辑
│   └── style.css      # 样式表
├── public/
│   ├── coolprop.js    # CoolProp 库
│   └── coolprop.wasm  # WASM 文件
└── README.md          # 项目说明
```

## 输入参数

- **制冷剂类型**：选择制冷剂
- **制冷量 Q**：kW（默认 100）
- **蒸发温度 Te**：°C（默认 5）
- **吸气过热度 ΔT_sh**：K（默认 5）
- **计算模式**：
  - 自动：基于液滴终端速度自动计算设计气速
  - 手动：手动输入设计气速
- **液滴直径 Dp**：µm（默认 152，范围 100-200）
- **设计气速系数**：0.75-0.9（默认 0.75）

## 输出结果

- 质量流量 (kg/s)
- 吸气气体密度 (kg/m³)
- 筒身内径 D (mm)
- 筒身高度 H (mm)
- 筒内容积 (L)
- 终端速度 vt (m/s)
- 液体密度 (kg/m³)
- 气体动力粘度 (Pa·s)

## 参考文档

计算逻辑基于 Alfa Laval 分离器设计指南和 HVAC 工程标准。

## License

MIT

