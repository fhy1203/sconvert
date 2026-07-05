# sconvert - 单细胞数据格式转换工具

实现 RDS (Seurat) 和 H5AD (AnnData) 格式的相互转换。

## 功能特点

- **H5AD → RDS**: 将 Scanpy/AnnData 格式转换为 Seurat 格式
- **RDS → H5AD**: 将 Seurat 格式转换为 Scanpy/AnnData 格式
- **中间文件策略**: 使用 CSV 作为中间格式，确保转换的可靠性和可追溯性
- **完整保留数据**: 保留表达矩阵和细胞元数据

## 环境要求

- Python 3.x
- R 4.x
- 必要的 Python 包: `scanpy`, `anndata`, `pandas`, `numpy`
- 必要的 R 包: `Seurat`, `SeuratObject`

## 安装

使用 conda 创建环境：

```bash
conda create -n sconvert -c conda-forge python r-seurat scanpy anndata pandas numpy
```

## 使用方法

### 1. H5AD 转换为 RDS

```bash
# 使用默认 X 矩阵
python sconvert.py h5ad2rds input.h5ad output.rds

# 指定要转换的 layer（建议指定原始表达矩阵）
python sconvert.py h5ad2rds input.h5ad output.rds --layer counts
```

**转换流程：**
1. Python 检查 H5AD 文件中的所有 layers，并显示每个 layer 的形状
2. 用户可以选择使用默认 X 矩阵或指定特定的 layer（建议指定原始表达矩阵）
3. Python 读取选定的表达矩阵，导出为 CSV
4. R 读取 CSV 文件，创建 Seurat 对象并保存为 RDS

**重要提示：**
- 建议只转换原始表达矩阵（Raw counts），因为 R 和 Python 的归一化、log1p 等处理逻辑可能不同
- 使用 `--layer` 参数明确指定哪个 layer 是原始表达矩阵
- 如果不指定，默认使用 `X` 主矩阵

### 2. RDS 转换为 H5AD

```bash
python sconvert.py rds2h5ad input.rds output.h5ad
```

**转换流程：**
1. R 读取 RDS 文件，导出表达矩阵和元数据为 CSV
2. Python 读取 CSV 文件，创建 AnnData 对象并保存为 H5AD

### 3. 指定临时目录

```bash
python sconvert.py h5ad2rds input.h5ad output.rds --temp-dir /path/to/temp
```

默认情况下，临时文件会保存在系统临时目录中，转换完成后自动清理。

## 示例

```bash
# 激活 conda 环境
conda activate sconvert

# 进入工作目录
cd /home/fhy/work/sconvert

# 示例 1: 检查 H5AD 文件中的 layers
python sconvert.py h5ad2rds singcell/pbmc3k_with_layers.h5ad test.rds
# 这会显示所有可用的 layers 及其数据类型

# 示例 2: H5AD 转 RDS (使用默认 X 矩阵)
python sconvert.py h5ad2rds singcell/pbmc3k.h5ad pbmc3k.rds

# 示例 3: H5AD 转 RDS (指定使用 counts layer)
python sconvert.py h5ad2rds singcell/pbmc3k_with_layers.h5ad pbmc3k_counts.rds --layer counts

# 示例 4: H5AD 转 RDS (指定使用 log1p layer)
python sconvert.py h5ad2rds singcell/pbmc3k_with_layers.h5ad pbmc3k_log1p.rds --layer log1p

# 示例 5: RDS 转 H5AD
python sconvert.py rds2h5ad pbmc3k.rds pbmc3k_converted.h5ad
```

## 转换说明

### 数据保留

- **表达矩阵**: 完整保留，包括稀疏矩阵的转换
- **细胞元数据**: 保留在 `obs` 中（H5AD）或 `meta.data` 中（RDS）
- **基因信息**: 保留基因名称

### 注意事项

1. **表达矩阵选择**:
   - 建议转换原始表达矩阵（如 counts），因为 R 和 Python 的归一化、log1p 等处理逻辑可能不同
   - 使用 `--layer` 参数可以明确指定要转换的数据层
   - 如果不指定，默认使用 `X` 主矩阵

2. **内存使用**: 转换过程中会将稀疏矩阵转换为密集矩阵，大文件可能需要较多内存

3. **基因名称**: Seurat 不支持基因名中包含下划线 `_`，会自动替换为短横线 `-`

4. **Seurat 版本**: 支持 Seurat 5.x 版本，使用 `layer` 参数而非旧的 `slot` 参数

5. **元数据**: 如果原始数据的元数据为空，转换后会创建空的元数据框

## 技术细节

### 中间文件格式

转换使用 CSV 作为中间格式：
- `matrix.csv`: 表达矩阵（细胞×基因 或 基因×细胞）
- `metadata.csv`: 细胞元数据
- `variables.csv`: 基因信息（仅 H5AD→RDS 方向）

### 矩阵方向

- **H5AD (AnnData)**: 表达矩阵为 (细胞 × 基因)
- **RDS (Seurat)**: 表达矩阵为 (基因 × 细胞)

转换过程中会自动处理矩阵转置。

## 测试

使用测试数据验证转换：

```bash
# 测试 H5AD -> RDS
python sconvert.py h5ad2rds singcell/pbmc3k.h5ad test_output.rds

# 测试 RDS -> H5AD
python sconvert.py rds2h5ad test_output.rds test_output.h5ad

# 验证转换结果
python -c "
import scanpy as sc
import numpy as np

orig = sc.read_h5ad('singcell/pbmc3k.h5ad')
conv = sc.read_h5ad('test_output.h5ad')

X_orig = orig.X.toarray() if hasattr(orig.X, 'toarray') else orig.X
X_conv = conv.X

print(f'形状匹配: {X_orig.shape == X_conv.shape}')
print(f'数值匹配: {np.allclose(X_orig, X_conv)}')
"
```

## 故障排除

### 问题：找不到 Rscript

确保 R 在 conda 环境中：
```bash
conda install -c conda-forge r-base
```

### 问题：内存不足

对于大数据集，考虑：
- 增加系统内存
- 使用更强大的服务器
- 分批处理数据

### 问题：Seurat API 变更

脚本已适配 Seurat 5.x，如果使用旧版本 Seurat，可能需要调整代码。

## 作者

凡python

## 许可证

MIT License
