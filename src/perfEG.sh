#!/bin/bash

NAME=".././perfEG"
CPU=3
METRICA="FLOPS_DP"

# Aumenta a frequência do processador
echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

# Exercuta o programa e imprime a os MFLOPS para cada região
out_text=$(likwid-perfctr  -C ${CPU} -g ${METRICA} -m ${NAME} < $1)
echo "$out_text" | awk '/^-+$/ {count++} count <= 3 {print} count == 3 {exit}'
echo "$out_text"| grep -i "FLOPS\| \ DP \[MFLOP/s\]"

# Diminui a frequência do processador
echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor
