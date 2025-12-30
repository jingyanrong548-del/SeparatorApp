import { defineConfig } from 'vite';

export default defineConfig({
  // 公共基础路径
  base: './',
  
  // 开发服务器配置
  server: {
    port: 3000,
    open: true,
    // 支持 WASM 文件
    headers: {
      'Cross-Origin-Embedder-Policy': 'require-corp',
      'Cross-Origin-Opener-Policy': 'same-origin'
    }
  },
  
  // 构建配置
  build: {
    outDir: 'dist',
    assetsDir: 'assets',
    // 确保 WASM 文件被正确处理
    rollupOptions: {
      output: {
        // 保持 WASM 文件在根目录
        assetFileNames: (assetInfo) => {
          if (assetInfo.name && assetInfo.name.endsWith('.wasm')) {
            return '[name][extname]';
          }
          return 'assets/[name]-[hash][extname]';
        }
      }
    }
  },
  
  // 优化配置
  optimizeDeps: {
    exclude: ['coolprop.js']
  },
  
  // 插件配置（如果需要）
  plugins: []
});

