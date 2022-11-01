%,Sb=1MVA,Vb=12.66KV
function DistFlow()

global data model

%% 设参
r = data.mpc.branch(:,4)*10/(12.66^2);
x = data.mpc.branch(:,5)*10/(12.66^2);

upstream = zeros(data.E_net.num_node,data.E_net.num_branch);
dnstream = zeros(data.E_net.num_node,data.E_net.num_branch);
for i = 1:data.E_net.num_branch
    upstream(i+1,i) = 1;
end
for i=[1:17,19:21,23:24,26:32]
    dnstream(i,i) = 1;
end
dnstream(2,18) = 1;
dnstream(3,22) = 1;
dnstream(6,25) = 1;

Vmax = [ones(1,data.Horizon)
        1.1^2*ones(data.E_net.num_node-1,data.Horizon)];
Vmin = [ones(1,data.Horizon)
        0.9^2*ones(data.E_net.num_node-1,data.Horizon)];

Pgmax=[ones(1,data.Horizon);zeros(4,data.Horizon);ones(1,data.Horizon);...
    zeros(3,data.Horizon);ones(1,data.Horizon);zeros(7,data.Horizon);...
    ones(1,data.Horizon);zeros(5,data.Horizon);ones(1,data.Horizon);...
    zeros(7,data.Horizon);ones(1,data.Horizon);zeros(1,data.Horizon)];
Qgmax=[ones(1,data.Horizon);zeros(4,data.Horizon);ones(1,data.Horizon);...
    zeros(3,data.Horizon);ones(1,data.Horizon);zeros(7,data.Horizon);...
    ones(1,data.Horizon);zeros(5,data.Horizon);ones(1,data.Horizon);...
    zeros(7,data.Horizon);ones(1,data.Horizon);zeros(1,data.Horizon)];

%% 2.设变量
model.V_2 = sdpvar(data.E_net.num_node,data.Horizon);%电压的平方
model.I_2 = sdpvar(data.E_net.num_branch,data.Horizon);%电流的平方
model.P_line = sdpvar(data.E_net.num_branch,data.Horizon);%线路有功
model.Q_line = sdpvar(data.E_net.num_branch,data.Horizon);%线路无功
model.P_gen = sdpvar(data.E_net.num_node,data.Horizon);%发电机有功
model.Q_gen = sdpvar(data.E_net.num_node,data.Horizon);%发电机无功

model.p_wt=sdpvar(data.num.Nunits_WT,data.Horizon,'full');%风机出力

model.judge.P_BT_charge = sdpvar(data.num.Nunits_BT,data.Horizon,'full');%ESS充电功率
model.judge.P_BT_discharge = sdpvar(data.num.Nunits_BT,data.Horizon,'full');%ESS放电功率
model.judge.U_BT_ch = binvar(data.num.Nunits_BT,data.Horizon,'full');%ESS充电状态
model.judge.U_BT_dis = binvar(data.num.Nunits_BT,data.Horizon,'full');%ESS放电状态
model.judge.SOC = sdpvar(data.num.Nunits_BT,data.Horizon+1,'full');
%ESS的电量，这个25的原因要搞懂才能理解储能首末功率相等的意思


model.judge.P_CHP=sdpvar(data.num.Nunits_CHP ,data.Horizon,'full');%chp电功率出力
model.judge.H_CHP=sdpvar(data.num.Nunits_CHP ,data.Horizon,'full');%chp热功率出力
model.judge.U_CHP=binvar(data.num.Nunits_CHP ,data.Horizon,'full');%chp开停机标记data.num.Nunits_CHP ,data.Horizon,'full');%chp电功率出力

model.judge.Pbuy=sdpvar(1,24,'full');%从电网购电电量
Psell=sdpvar(1,24,'full');%向电网售电电量
model.judge.Pnet=sdpvar(1,24,'full');%交换功率
model.judge.Temp_net=binvar(1,24,'full'); % 购|售电标志

%% 3.设约束
%% 与电网进行交易，买卖电
for i=1:data.Horizon
    model.st = [model.st, 0<=model.judge.Pnet(i)<=1000/10000,...
        0<=model.judge.Pbuy(i)<=0.1,...
        0<=Psell(i)<=0]; %主网功率交换约束
    model.st = [model.st, implies(model.judge.Temp_net(i),...
        [model.judge.Pnet(i)>=0,model.judge.Pbuy(i)==model.judge.Pnet(i),Psell(i)==0])]; %购电情况约束
    model.st = [model.st, implies(1-model.judge.Temp_net(i),...
        [model.judge.Pnet(i)<=0,Psell(i)==model.judge.Pnet(i),model.judge.Pbuy(i)==0])]; %售电情况约束   
end
model.P_buy=[model.judge.Pnet;zeros(9,data.Horizon); zeros(23,data.Horizon);];
%% 风力发电机组，约束：风电实际出力小于风力预测
model.st = [model.st, 0 <= model.p_wt,model.p_wt<= data.P_wt*1];
P_wt_all=[zeros(5,data.Horizon);model.p_wt(1,:);zeros(11,data.Horizon);...
    model.p_wt(2,:);zeros(5,data.Horizon);model.p_wt(3,:);...
    zeros(7,data.Horizon);model.p_wt(4,:);zeros(1,data.Horizon)];
%% 储能装置（ESS）约束
%充放电状态约束
model.st = [model.st, model.judge.U_BT_dis + model.judge.U_BT_ch <= 1];%表示充电，放电，不充不放三种状态
%功率约束
for k=1:data.num.Nunits_BT
model.st = [model.st, 0 <= model.judge.P_BT_discharge(k,:) <= model.judge.U_BT_dis(k,:)*0.05];
model.st = [model.st, 0 <= model.judge.P_BT_charge(k,:) <= model.judge.U_BT_ch(k,:)*0.05];
end
%容量约束
for t = 1:data.Horizon
model.st = [model.st, model.judge.SOC(:,t+1) == ...
    model.judge.SOC(:,t) + 0.92*model.judge.P_BT_charge(:,t) - 1.08*model.judge.P_BT_discharge(:,t)]; 
end
model.st = [model.st, model.judge.SOC(:,1) ==model.judge.SOC(:,25)];
model.st = [model.st, model.judge.SOC(1,1) ==0.02];model.st = [model.st, model.judge.SOC(2,1) ==0.02];
model.st = [model.st, model.judge.SOC(3,1) ==0];model.st = [model.st, model.judge.SOC(4,1) ==0];
model.st = [model.st, 0<=model.judge.SOC(1,:) <= 0.1];
model.st = [model.st, 0<=model.judge.SOC(2,:) <= 0.1];
model.st = [model.st, 0 <=model.judge.SOC(3,:) <=0.1];
model.st = [model.st, 0 <=model.judge.SOC(4,:)<=0.1];

%投入节点选择
model.P_dch = [zeros(5,data.Horizon);model.judge.P_BT_discharge(1,:);zeros(11,data.Horizon);...
    model.judge.P_BT_discharge(2,:);zeros(5,data.Horizon);model.judge.P_BT_discharge(3,:);...
    zeros(7,data.Horizon);model.judge.P_BT_discharge(4,:);zeros(1,data.Horizon)];
model.P_ch = [zeros(5,data.Horizon);model.judge.P_BT_charge(1,:);zeros(11,data.Horizon);...
    model.judge.P_BT_charge(2,:);zeros(5,data.Horizon);model.judge.P_BT_charge(3,:);...
    zeros(7,data.Horizon);model.judge.P_BT_charge(4,:);zeros(1,data.Horizon)];
%disp(model.P_dch);
%% CHP机组约束
%热电出力上下限约束
 for t = 1:data.Horizon
     model.st = [model.st,((0 <= model.judge.P_CHP(:,t)<= 1):'hahahahahahaha')];
     model.st = [model.st,0 <= model.judge.H_CHP(:,t) ];
 end
 for t = 2:data.Horizon
     model.st = [model.st,-0.25 <= model.judge.P_CHP(:,t)-model.judge.P_CHP(:,t-1)<=0.2500];%爬坡约束
 end
model.P_chp_all=[zeros(9,data.Horizon); (model.judge.P_CHP(1,:));zeros(23,data.Horizon);];

%% 潮流约束
%松弛后的节点功率约束
Pin = -upstream*model.P_line + upstream*(model.I_2.*(r*ones(1,data.Horizon))) + dnstream*model.P_line;%节点注入有功
Qin = -upstream*model.Q_line + upstream*(model.I_2.*(x*ones(1,data.Horizon))) + dnstream*model.Q_line;%节点注入无功
%电平衡
model.st = [model.st, Pin+data.mpc.Pload- model.P_chp_all-P_wt_all-model.P_buy - model.P_dch + model.P_ch == 0];
model.P_gen=model.P_buy+model.P_chp_all+P_wt_all + model.P_dch - model.P_ch ;
model.st = [model.st, Qin+data.mpc.Qload-model.Q_gen == 0];

%松弛后的电压平衡约束
model.st = [model.st, model.V_2(data.mpc.branch(:,3),:) == ...
    model.V_2(data.mpc.branch(:,2),:) - 2*(r*ones(1,data.Horizon)).*model.P_line - ...
    2*(x*ones(1,data.Horizon)).*model.Q_line + ((r.^2+x.^2)*ones(1,data.Horizon)).*model.I_2];
%支路潮流二阶锥约束
model.st = [model.st, model.V_2(data.mpc.branch(:,2),:).*model.I_2 >= model.P_line.^2 + model.Q_line.^2];

%% 普通约束
%节点电压约束
model.st = [model.st, Vmin <= model.V_2,model.V_2 <= Vmax];
%发电机功率约束
model.st = [model.st, 0 <= abs(model.P_gen) <= Pgmax,-Qgmax <= model.Q_gen,model.Q_gen <= Qgmax];
%支路电流约束
 model.st = [model.st, 0 <= model.I_2,model.I_2 <= 1.1];
model.P_in=Pin;

