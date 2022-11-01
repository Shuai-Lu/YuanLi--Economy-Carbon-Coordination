function read_heat51(h_filename)
global data model

h_filename = 'testdata_33bus.xlsx';

data.heat.node_info=xlsread(h_filename,'heatingnet.node','C4:K54');%热网节点信息
data.heat.pipe_info=xlsread(h_filename,'heatingnet.pipe','D4:N53');%热网管道信息
data.heat.build_info=xlsread(h_filename,'buildings','C4:G29');%负荷信息
data.heat.T_out=xlsread(h_filename,'profiles','F41:AC41');%负荷信息

data.heat.num_node=size(data.heat.node_info);
data.heat.num_node=data.heat.num_node(1,1);%节点数51
data.heat.num_pipe=size(data.heat.pipe_info);
data.heat.num_pipe=data.heat.num_pipe(1,1);%供回水管道数50（*2）
data.heat.num_load=size(data.heat.build_info);
data.heat.num_load=data.heat.num_load(1,1);%热用户数目26
