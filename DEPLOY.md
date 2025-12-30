# 部署指南

## GitHub 部署步骤

### 1. 创建 GitHub 仓库（如果还没有）

1. 访问 https://github.com/new
2. 仓库名称：`SeparatorApp`
3. 选择 **Public** 或 **Private**
4. **不要** 初始化 README、.gitignore 或 license（我们已经有了）
5. 点击 "Create repository"

### 2. 推送代码到 GitHub

如果仓库已存在，直接运行：

```bash
cd /Users/jingyanrong/Desktop/MyGitHubProjects/SeparatorApp
git push -u origin main
```

如果仓库不存在或需要重新设置远程地址：

```bash
# 检查当前远程地址
git remote -v

# 如果地址不对，设置正确的 SSH 地址
git remote set-url origin git@github.com:jingyanrong548-del/SeparatorApp.git

# 推送代码
git push -u origin main
```

### 3. 验证推送

访问 https://github.com/jingyanrong548-del/SeparatorApp 确认代码已上传。

---

## Vercel 部署步骤

### 方法 1：通过 Vercel Dashboard（推荐）

1. **登录 Vercel**
   - 访问 https://vercel.com
   - 使用 GitHub 账号登录

2. **导入项目**
   - 点击 "Add New..." → "Project"
   - 选择 `jingyanrong548-del/SeparatorApp` 仓库
   - 点击 "Import"

3. **配置项目**
   - **Framework Preset**: Other
   - **Root Directory**: `./` (默认)
   - **Build Command**: 留空（静态网站，无需构建）
   - **Output Directory**: 留空
   - 点击 "Deploy"

4. **等待部署完成**
   - 通常需要 1-2 分钟
   - 部署完成后会显示访问链接

### 方法 2：通过 Vercel CLI

```bash
# 安装 Vercel CLI（如果还没有）
npm i -g vercel

# 在项目目录中登录
cd /Users/jingyanrong/Desktop/MyGitHubProjects/SeparatorApp
vercel login

# 部署
vercel

# 生产环境部署
vercel --prod
```

### 方法 3：通过 GitHub 集成（自动部署）

1. 在 Vercel Dashboard 中导入项目后
2. 每次推送到 `main` 分支会自动触发部署
3. 无需手动操作

---

## 注意事项

### WASM 文件配置

Vercel 已经配置了正确的 HTTP 头来支持 WASM 文件：
- `Cross-Origin-Embedder-Policy: require-corp`
- `Cross-Origin-Opener-Policy: same-origin`

### 文件大小限制

- `coolprop.js` 和 `coolprop.wasm` 文件较大
- Vercel 免费版支持最大 100MB 的文件
- 如果超过限制，考虑使用 CDN 或优化文件

### 自定义域名（可选）

1. 在 Vercel Dashboard 中进入项目设置
2. 点击 "Domains"
3. 添加你的自定义域名

---

## 故障排除

### SSH 连接问题

如果推送时遇到权限问题：

```bash
# 测试 SSH 连接
ssh -T git@github.com

# 如果显示 "Permission denied"，检查 SSH 密钥
ls -la ~/.ssh

# 如果没有 id_rsa 或 id_ed25519，生成新的 SSH 密钥
ssh-keygen -t ed25519 -C "your_email@example.com"

# 添加 SSH 密钥到 ssh-agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

# 复制公钥到剪贴板
pbcopy < ~/.ssh/id_ed25519.pub

# 然后到 GitHub Settings → SSH and GPG keys → New SSH key 添加
```

### Vercel 部署失败

1. 检查 `vercel.json` 配置是否正确
2. 查看 Vercel Dashboard 中的构建日志
3. 确保所有必需文件都已提交到 GitHub

### 浏览器控制台错误

如果部署后出现 `import.meta` 错误：
- 确保使用 ES6 模块（已配置）
- 检查 Vercel 是否正确设置了 HTTP 头

---

## 更新部署

每次代码更新后：

```bash
# 提交更改
git add .
git commit -m "更新描述"

# 推送到 GitHub
git push origin main

# Vercel 会自动部署（如果已配置 GitHub 集成）
# 或手动运行：vercel --prod
```

