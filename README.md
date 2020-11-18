GEOS_CHEM_ADJOINT
Adjoint: 伴随文件夹
OCO2_CO2_Litev9_ND_Observation_Error_Covariance： OCO-2天地模式观测数据

#输入数据
初始浓度数据 通量数据 可以从GEOS-Chem数据目录下载

#编译环境
编译器: gnu
linux系统: ubuntu 16.03

#运行过程
1.在rundirs/merra2_2x25_CO2/目录下运行make命令，得到可执行程序geos
2.准备观测数据
2.先运行rundirs/merra2_2x25_CO2/BEC/bec.sh，这个脚本用来计算背景误差协方差
3.在gcadj.run.bac里面设置同化时间
4.运行gcadj.run.bac进行同化
