#!/bin/bash
# sconvert 快速使用示例

# 设置路径
SCONVERT_DIR="/home/fhy/work/sconvert"
CONDA_ENV="sconvert"
PYTHON="/home/fhy/biosoft/miniforge/envs/${CONDA_ENV}/bin/python"

cd "$SCONVERT_DIR"

echo "=========================================="
echo "sconvert 使用示例"
echo "=========================================="
echo

# 示例 1: H5AD -> RDS
echo "示例 1: H5AD 转换为 RDS"
echo "命令: python sconvert.py h5ad2rds singcell/pbmc3k.h5ad pbmc3k.rds"
$PYTHON sconvert.py h5ad2rds singcell/pbmc3k.h5ad pbmc3k.rds
echo

# 示例 2: RDS -> H5AD
echo "示例 2: RDS 转换为 H5AD"
echo "命令: python sconvert.py rds2h5ad pbmc3k.rds pbmc3k_back.h5ad"
$PYTHON sconvert.py rds2h5ad pbmc3k.rds pbmc3k_back.h5ad
echo

# 示例 3: 验证转换
echo "示例 3: 验证转换结果"
$PYTHON -c "
import scanpy as sc
import numpy as np

orig = sc.read_h5ad('singcell/pbmc3k.h5ad')
conv = sc.read_h5ad('pbmc3k_back.h5ad')

X_orig = orig.X.toarray() if hasattr(orig.X, 'toarray') else orig.X
X_conv = conv.X

print(f'原始数据: {orig.n_obs} 细胞 x {orig.n_vars} 基因')
print(f'转换后数据: {conv.n_obs} 细胞 x {conv.n_vars} 基因')
print(f'形状匹配: {X_orig.shape == X_conv.shape}')
print(f'数值匹配: {np.allclose(X_orig, X_conv)}')
print(f'最大差异: {np.max(np.abs(X_orig - X_conv))}')
"

echo
echo "=========================================="
echo "转换完成！"
echo "=========================================="
echo
echo "=========================================="
echo "空间转录组转换示例"
echo "=========================================="
echo

# 示例 4: 空间转录组 H5AD -> RDS
echo "示例 4: 空间转录组 H5AD 转换为 RDS"
echo "命令: python sconvert.py h5ad2rds spatial/visuim_test.h5ad spatial_output.rds --spatial"
$PYTHON sconvert.py h5ad2rds spatial/visuim_test.h5ad spatial_output.rds --spatial
echo

# 示例 5: 空间转录组 RDS -> H5AD
echo "示例 5: 空间转录组 RDS 转换为 H5AD"
echo "命令: python sconvert.py rds2h5ad spatial_output.rds spatial_back.h5ad --spatial"
$PYTHON sconvert.py rds2h5ad spatial_output.rds spatial_back.h5ad --spatial
echo

# 示例 6: 验证空间转录组转换
echo "示例 6: 验证空间转录组转换结果"
$PYTHON -c "
import scanpy as sc
import numpy as np

orig = sc.read_h5ad('spatial/visuim_test.h5ad')
conv = sc.read_h5ad('spatial_back.h5ad')

print(f'原始数据: {orig.n_obs} 细胞 x {orig.n_vars} 基因')
print(f'转换后数据: {conv.n_obs} 细胞 x {conv.n_vars} 基因')

# 检查空间坐标
if 'spatial' in orig.obsm and 'spatial' in conv.obsm:
    orig_coords = orig.obsm['spatial']
    conv_coords = conv.obsm['spatial']
    print(f'原始坐标: {orig_coords.shape}')
    print(f'转换后坐标: {conv_coords.shape}')
    print(f'坐标匹配: {np.allclose(orig_coords, conv_coords)}')
else:
    print('警告: 空间坐标缺失')

# 检查缩放因子
if 'spatial' in orig.uns and 'spatial' in conv.uns:
    print('✓ 缩放因子保留完整')
else:
    print('警告: 缩放因子缺失')
"

echo
echo "=========================================="
echo "所有示例完成！"
echo "=========================================="
