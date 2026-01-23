#!/bin/bash

# AlphaFold Database 下载和解压脚本
# 使用方法: ./download_and_extract.sh

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 输入文件
LIST_FILE="af2db_list.txt"

# 检查列表文件是否存在
if [ ! -f "$LIST_FILE" ]; then
    echo -e "${RED}错误: 找不到 $LIST_FILE 文件${NC}"
    exit 1
fi

# 创建下载目录
DOWNLOAD_DIR="af2db_downloads"
mkdir -p "$DOWNLOAD_DIR"

echo -e "${GREEN}开始下载和解压 AlphaFold 数据库...${NC}"
echo "======================================"

# 读取文件并处理每一行
while IFS=', ' read -r name url || [ -n "$name" ]; do
    # 跳过空行
    [ -z "$name" ] && continue
    
    echo -e "\n${YELLOW}处理: $name${NC}"
    echo "URL: $url"
    
    # 提取文件名
    filename=$(basename "$url")
    filepath="$DOWNLOAD_DIR/$filename"
    
    # 下载文件
    if [ -f "$filepath" ]; then
        echo -e "${YELLOW}文件已存在，跳过下载: $filename${NC}"
    else
        echo "正在下载 $filename ..."
        wget -c -P "$DOWNLOAD_DIR" "$url"
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ 下载完成${NC}"
        else
            echo -e "${RED}✗ 下载失败${NC}"
            continue
        fi
    fi
    
    # 解压文件
    extract_dir="$DOWNLOAD_DIR/${name}"
    mkdir -p "$extract_dir"
    
    echo "正在解压到 $extract_dir ..."
    tar -xf "$filepath" -C "$extract_dir"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ 解压完成${NC}"
    else
        echo -e "${RED}✗ 解压失败${NC}"
    fi
    
done < "$LIST_FILE"

echo -e "\n${GREEN}======================================"
echo "所有任务完成！"
echo "文件保存在: $DOWNLOAD_DIR"
echo -e "======================================${NC}"
