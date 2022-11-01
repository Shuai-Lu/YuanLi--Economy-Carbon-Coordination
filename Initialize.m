function Initialize()
clc;
clear;
global data model

%% 设参
data.r = data.mpc.branch(:,4)*10/(12.66^2);
data.x = data.mpc.branch(:,5)*10/(12.66^2);

data.upstream = zeros(data.E_net.num_node,data.E_net.num_branch);
data.dnstream = zeros(data.E_net.num_node,data.E_net.num_branch);
for i = 1:data.E_net.num_branch
    data.upstream(i+1,i) = 1;
end
for i=[1:17,19:21,23:24,26:32]
    data.dnstream(i,i) = 1;
end
data.dnstream(2,18) = 1;
data.dnstream(3,22) = 1;
data.dnstream(6,25) = 1;

data.Vmax = [ones(1,data.Horizon)
        1.1^2*ones(data.E_net.num_node-1,data.Horizon)];
data.Vmin = [ones(1,data.Horizon)
        0.9^2*ones(data.E_net.num_node-1,data.Horizon)];
data.Pgmax=[ones(1,data.Horizon);zeros(4,data.Horizon);ones(1,data.Horizon);...
    zeros(3,data.Horizon);ones(1,data.Horizon);zeros(7,data.Horizon);...
    ones(1,data.Horizon);zeros(5,data.Horizon);ones(1,data.Horizon);...
    zeros(7,data.Horizon);ones(1,data.Horizon);zeros(1,data.Horizon)];
data.Qgmax=[ones(1,data.Horizon);zeros(4,data.Horizon);ones(1,data.Horizon);...
    zeros(3,data.Horizon);ones(1,data.Horizon);zeros(7,data.Horizon);...
    ones(1,data.Horizon);zeros(5,data.Horizon);ones(1,data.Horizon);...
    zeros(7,data.Horizon);ones(1,data.Horizon);zeros(1,data.Horizon)];
