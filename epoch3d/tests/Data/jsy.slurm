#!/bin/bash
#SBATCH --job-name=EPOCH3D
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --nodelist=cluster001
#SBATCH --ntasks=256
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# 打印作业信息
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"

# 进入提交任务的目录
cd $SLURM_SUBMIT_DIR

# 加载运行时所需的模块
module load gnu12/12.2.0 openmpi/4.1.4 phdf5/1.10.8

# 启动 EPOCH3D 程序
# 输入 "." 作为输出目录，程序会自动找到当前目录下的 input.deck
# 并将模拟结果输出到当前目录
echo "." | mpirun -np 256 ../../bin/epoch3d

# 作业完成信息
echo "End time: $(date)"
