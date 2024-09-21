#!bin/sh
clear
g++ -std=c++11 main.cpp -o main.out -I/$MOSEK_HOME/$MOSEK_VERSION/tools/platform/linux64x86/h/ -L/$MOSEK_HOME/$MOSEK_VERSION/tools/platform/linux64x86/bin/ -Wl,-rpath-link,/$MOSEK_HOME/$MOSEK_VERSION/tools/platform/linux64x86/bin/ -Wl,-rpath=/$MOSEK_HOME/$MOSEK_VERSION/tools/platform/linux64x86/bin/   -lmosek64 -lfusion64
./main.out