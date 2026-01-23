#!/bin/bash

# 下载人类 AlphaFold 数据库
# Human (Homo sapiens) - UniProt Proteome: UP000005640

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}======================================"
echo "下载人类 AlphaFold 数据库"
echo "======================================${NC}"

# 创建目录
BASE_DIR="/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads"
HUMAN_DIR="$BASE_DIR/Human"
mkdir -p "$HUMAN_DIR"

# 人类数据URL
HUMAN_URL="https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v6.tar"
TAR_FILE="$HUMAN_DIR/UP000005640_9606_HUMAN_v6.tar"

echo -e "\n${YELLOW}步骤 1: 下载人类蛋白质组数据${NC}"
echo "URL: $HUMAN_URL"
echo "目标: $TAR_FILE"

if [ -f "$TAR_FILE" ]; then
    echo -e "${YELLOW}文件已存在，跳过下载${NC}"
else
    echo "开始下载..."
    wget -c -O "$TAR_FILE" "$HUMAN_URL"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ 下载完成${NC}"
    else
        echo -e "${RED}✗ 下载失败${NC}"
        exit 1
    fi
fi

echo -e "\n${YELLOW}步骤 2: 解压文件${NC}"
echo "解压到: $HUMAN_DIR"

tar -xf "$TAR_FILE" -C "$HUMAN_DIR"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ 解压完成${NC}"
else
    echo -e "${RED}✗ 解压失败${NC}"
    exit 1
fi

echo -e "\n${YELLOW}步骤 3: 统计文件信息${NC}"
PDB_COUNT=$(find "$HUMAN_DIR" -name "*.pdb.gz" -o -name "*.pdb" | wc -l)
CIF_COUNT=$(find "$HUMAN_DIR" -name "*.cif.gz" -o -name "*.cif" | wc -l)
echo "PDB 文件数量: $PDB_COUNT"
echo "CIF 文件数量: $CIF_COUNT"

echo -e "\n${GREEN}======================================"
echo "下载完成！"
echo "数据保存在: $HUMAN_DIR"
echo "======================================${NC}"

# 询问是否删除tar文件
read -p "是否删除原始tar文件以节省空间? (yes/no): " response
if [ "$response" = "yes" ] || [ "$response" = "y" ]; then
    rm "$TAR_FILE"
    echo -e "${GREEN}✓ 已删除tar文件${NC}"
fi
