%%
% 1. 模块化
% 2. 结构体
% 3. 功能封装为子函数
clc;
clear;
%% data:  data, model
global data model;

%% read data
filename = 'testdata_33bus.xlsx';
read_IEEE33();%读取网络参数
h_filename = 'testdata_33bus.xlsx';
read_heat51(h_filename);

%% initializing
Initialize();

%% modeling
model = [];
model.st = [];
model.objective = [];

%% define variables
DistFlow();
Heat_Net();

%% define obj
define_objective();
Pareto_NNC();%NNC求前沿
CERC_topsis();%对Pareto上的点进行分析