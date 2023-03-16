#!/bin/bash
#这是一个简单的运行脚本，具体的可以直接看代码

shopt -s  extglob

clear(){
    rm -rf !(CMakeLists.txt|.vscode|CodeCraft_zip.sh|main|main.cpp|run.sh)
    echo "Clear Done!"
}

if [ "$1" = "clear" ]
then 
    clear
elif [ "$1" = "compile" ]
then
    cmake CMakeLists.txt
    make
else
    clear
    cmake CMakeLists.txt
    make
    clear
    echo "All Work Done!"
fi