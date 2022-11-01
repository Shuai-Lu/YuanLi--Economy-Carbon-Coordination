function read_IEEE33()
global data model

%% 读取网络参数SB=10MVA,VB=12.66kV
data.mpc.bus=xlsread('testdata_33bus','mpc.bus','C4:O36');
data.mpc.branch=xlsread('testdata_33bus','mpc.branch','C4:P35');
data.mpc.device=xlsread('testdata_33bus','mpc.device ','C4:R13');
data.mpc.load=xlsread('testdata_33bus','profiles','F4:AC36');
data.mpc.cost=xlsread('testdata_33bus','mpc.cost','C18:Z19');
base=xlsread('testdata_33bus','profiles','E4:E36');
Base=base*ones(1,24);
P_wt_b=xlsread('testdata_33bus','profiles','F37:AC40');
base_wt=xlsread('testdata_33bus','profiles','E37:E40');
data.P_wt=P_wt_b.*(base_wt*ones(1,24))/10000;



%% 设置参数
data.E_net.num_node=size(data.mpc.bus);
data.E_net.num_node=data.E_net.num_node(1,1);%节点数33
data.E_net.num_branch=size(data.mpc.branch);
data.E_net.num_branch=data.E_net.num_branch(1,1);%网络中的支路数32
data.mpc.pload2=sum(data.mpc.load)/2;%各时刻系统中的总负荷数
%% 规模变量
%时间范围
data.Horizon = 24;
%机组数
data.num.Nunits_WT = 4;%风力发电机组
data.num.Nunits_BT = 4;%蓄电池
data.num.Nunits_CHP = 1;%CHP机组
data.num.Nunits_EB = 1;%电锅炉
%% 
data.mpc.Pload= data.mpc.load.*Base/10000;
data.Eload=sum(data.mpc.Pload);
data.mpc.Qload=xlsread('testdata_33bus','profiles','F47:AC79');
