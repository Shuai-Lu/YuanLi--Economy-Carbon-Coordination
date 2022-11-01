%% 热网建模
%利用'testdata_6bus.xlsx'中六节点的热网模型
function Heat_Net()
global data model

%% 参数设置
C_water=4200;%定义水的比热容
density_w=1;%水的密度,t/m3
n_H_EX=0.98;%换热设备平均效率
delta_t=1;%单位h

%管道相关参数
pipe_fnode=data.heat.pipe_info(:,1)+1;%管道首节点
pipe_tnode=data.heat.pipe_info(:,2)+1;%管道末端节点
L_pipe=data.heat.pipe_info(:,3);%管道长度
D_pipe=data.heat.pipe_info(:,4);%管道内径
rough_pipe=data.heat.pipe_info(:,5);%管道粗糙度
R_pipe=data.heat.pipe_info(:,6)/1000;%conductivity传导率，热损失系数,管道单位长度热阻kW/m*C
q_pipe=data.heat.pipe_info(:,9);%管道流量,单位转换为t/h
t_bs_in=data.heat.pipe_info(:,10);%供水管道的首端进水初始温度
t_br_in=data.heat.pipe_info(:,11);%回水管道的首端进水初始温度

%节点相关参数
node_type=data.heat.node_info(:,2);%节点类型，0-发电机，1-连接节点，2-负荷节点
T_s_border=data.heat.node_info(:,3:4);%节点供水温度上下限，T_s_border（：,1）表示下限，T_s_border（：,2）表示上限
T_r_border=data.heat.node_info(:,5:6);%节点回水温度上下限，T_r_border（：,1）表示下限，T_r_border（：,2）表示上限
Pl_min=50000;%
T_out=-25*data.heat.T_out';%24*1,24h室外温度变化

%% 热网决策变量
data.num.Nunits_HS=1;
model.judge.W_HS = sdpvar(data.num.Nunits_HS,data.Horizon+1 ,'full');%HS在t时段内,总的蓄/放热量
model.judge.H_HS_in = sdpvar(data.num.Nunits_HS,data.Horizon,'full');
model.judge.U_HS_in = binvar(data.num.Nunits_HS,data.Horizon,'full');%蓄热槽充热
model.judge.H_HS_out = sdpvar(data.num.Nunits_HS,data.Horizon,'full');
model.judge.U_HS_out = binvar(data.num.Nunits_HS,data.Horizon,'full');%蓄热槽放热

model.heat_judge.T_s_in=sdpvar(data.heat.num_pipe,data.Horizon,'full');%5*24，供水管道首端温度
model.heat_judge.T_r_in=sdpvar(data.heat.num_pipe,data.Horizon,'full');%5*24，回水管道首端温度
model.heat_judge.T_s_out=sdpvar(data.heat.num_pipe,data.Horizon,'full');%5*24，供水管道末端温度
model.heat_judge.T_r_out=sdpvar(data.heat.num_pipe,data.Horizon,'full');%5*24，回水管道末端温度
model.heat_judge.T_s_delay_out=sdpvar(data.heat.num_pipe,data.Horizon,'full');%5*24，供水管道末端温度
model.heat_judge.T_r_delay_out=sdpvar(data.heat.num_pipe,data.Horizon,'full');%5*24，回水管道末端温度
model.heat_judge.T_s_node=sdpvar(data.heat.num_node,data.Horizon,'full');%6*24，供水管道节点温度
model.heat_judge.T_r_node=sdpvar(data.heat.num_node,data.Horizon,'full');%6*24，回水管道节点温度



model.heat_judge.pd=sdpvar(data.Horizon,data.heat.num_load,'full');%24*26热负荷大小
model.heat_judge.t_in=sdpvar(data.heat.num_load,data.Horizon+1,'full');%26*24建筑物室温
%% 加入燃气锅炉
model.judge.H_GB=sdpvar(1 ,data.Horizon,'full');%燃气锅炉出力
model.judge.U_GB=binvar(1 ,data.Horizon,'full');%燃气锅炉状态
for t = 1:data.Horizon
     model.st = [model.st,((0 <= model.judge.H_GB(:,t) <= 10000):'hahahahahahaha')];
     model.st = [model.st,0 <= model.judge.U_GB(:,t)<=1];
end
for t = 2:data.Horizon
     model.st = [model.st,-2500 <= model.judge.H_GB(:,t)-model.judge.H_GB(:,t-1)<=2500];%爬坡约束
end
%% 储热设备约束
model.judge.W_HS(1,1)=1500;
Self_h_loss=0.03; n_in=0.95;%HS的能量自损系数%HS储热效率
n_out=0.95;%HS放热效率

%状态标志位满足互斥约束g
model.st = [model.st, model.judge.U_HS_in+model.judge.U_HS_out <= 1];
%蓄放热功率约束
for k=1:data.num.Nunits_HS
    model.st = [model.st, 0 <= model.judge.H_HS_out(k,:) <= model.judge.U_HS_out(k,:)*1500];
    model.st = [model.st, 0 <= model.judge.H_HS_in(k,:) <=model.judge.U_HS_in(k,:)*2000];
end
%容量约束
for t = 2:data.Horizon+1
    model.st = [model.st, model.judge.W_HS(:,t) == ...
        model.judge.W_HS(:,t-1)*(1-Self_h_loss) + (n_in*model.judge.H_HS_in(:,t-1) - model.judge.H_HS_out(:,t-1)/n_out)];   
end
model.st = [model.st,model.judge.W_HS(:,1) == model.judge.W_HS(:,25)];
model.st = [model.st,300 <= model.judge.W_HS <=4000];

%% 节点法
%相关参数计算
r_p=ceil((density_w*(pi/4*D_pipe.^2).*L_pipe)./(q_pipe*delta_t))-1;
R_p=(1+r_p).*q_pipe*delta_t;
k_p=(R_p-density_w*(pi/4*D_pipe.^2).*L_pipe)./(q_pipe*delta_t);%
n_p=1-exp(-R_pipe*(delta_t*3600).*(r_p+3/2-k_p)./(density_w*C_water*(pi/4*D_pipe.^2)));%单位是(m-1)
% disp(r_p);
for t=1:data.Horizon
    for j=1:data.heat.num_pipe
        %热网管道的传输延时
        time_delay1=t-r_p(j,1)-1;  
        time_delay2=t-r_p(j,1);
        %判断model.heat_judge.T_s_in(j,t-1-r_p(j,1))中的第二个索引是不是负数
        if time_delay1<=0
            T_s_in1=t_bs_in(j,1);
            T_r_in1=t_br_in(j,1);
%             disp(T_s_in1);
%             disp(T_r_in1);
        end
        if time_delay1>0
            T_s_in1=model.heat_judge.T_s_in(j,t-1-r_p(j,1));
            T_r_in1=model.heat_judge.T_r_in(j,t-1-r_p(j,1));
        end
        if time_delay2<=0
            T_s_in2=t_bs_in(j,1);
            T_r_in2=t_br_in(j,1);
        end
        if time_delay2>0
            T_s_in2=model.heat_judge.T_s_in(j,t-r_p(j,1));
            T_r_in2=model.heat_judge.T_r_in(j,t-r_p(j,1));
        end 
        model.heat_judge.T_s_delay_out(j,t)=(1-k_p(j,1))*T_s_in1+k_p(j,1)*T_s_in2;
        model.heat_judge.T_r_delay_out(j,t)=(1-k_p(j,1))*T_r_in1+k_p(j,1)*T_r_in2;   
    end
end

%延时加上热损
model.st = [model.st ,((model.heat_judge.T_s_out ==...
    ((1-n_p)*ones(1,data.Horizon)).*(model.heat_judge.T_s_delay_out-ones(data.heat.num_pipe,1)*T_out')+ones(data.heat.num_pipe,1)*T_out'):'DELAY1')];

% model.st = [model.st , ((value(model.heat_judge.T_r_out)==(1-n_p).*(model.heat_judge.T_r_delay_out(:,t)-T_out(t,1))+T_out(t,1)): 'T_r_out')];
model.st = [model.st , ((model.heat_judge.T_r_out==...
    ((1-n_p)*ones(1,data.Horizon)).*(model.heat_judge.T_r_delay_out-ones(data.heat.num_pipe,1)*T_out')+ones(data.heat.num_pipe,1)*T_out'):'delay2')];

%节点处静态热平衡方程
%找出热源节点，交汇节点和负荷节点所在位置集合
source_node=find(node_type==0);%源
load_node=find(node_type==2);%荷
in_node=find(node_type==1);%网络中节点
num_source=size(source_node);
num_source=num_source(1,1);
num_load=size(load_node);
num_load=num_load(1,1);
num_inode=size(in_node);
num_inode=num_inode(1,1);
%对于每一个网络中的节点
for k=1:num_inode
    p=in_node(k,1);%记录该节点的编号
    %找到每个节点流向它的管道
    up_pipe=find(pipe_tnode==p);%管道编号
    num_uppipe=size(up_pipe);
    num_uppipe=num_uppipe(1,1);
    %找到每个节点流出它的管道
    down_pipe=find(pipe_fnode==p);
    num_downpipe=size(down_pipe);
    num_downpipe=num_downpipe(1,1); 
    sum_qout_r=0;
    %生成一个存放该节点流向他的管道的质量流量矩阵,供水     
    m_all=zeros(data.heat.num_pipe,data.Horizon);%5*24
    %回水
    up_pipe_r_m=zeros(num_downpipe,1);  
    m_r_all=zeros(data.heat.num_pipe,data.Horizon);%5*24
    %节点处的热平衡方程    
    for i=1:num_uppipe
        up_pipe_m=q_pipe(up_pipe(i,1),1);
        up_pipe_M=up_pipe_m*ones(1,data.Horizon);%行数为流向该节点的管道数，列数为24
        m_all(up_pipe(i,1),:)=up_pipe_M;
%         %供水
        %回水的时候原来流向它的管道，现在相当于是从这个节点流出，所以要计算总的流出的流量
        sum_qout_r = sum_qout_r +q_pipe(up_pipe(i,1),1);%2*1      
    end
    for i=1:num_downpipe
        up_pipe_r_m=q_pipe(down_pipe(i,1),1);
        up_pipe_r_M=up_pipe_r_m*ones(1,data.Horizon);%行数为流向该节点的管道数，列数为24
        m_r_all(down_pipe(i,1),:)=up_pipe_r_M;
    end
    h_s_all=m_all.*model.heat_judge.T_s_out;
    h_s_all=sum(h_s_all);
    h_r_all=m_r_all.*model.heat_judge.T_r_out;
    h_r_all=sum(h_r_all);
    model.st = [model.st ,((model.heat_judge.T_s_node(p,:)*sum_qout_r == h_s_all): 'T_s_node')];
    model.st = [model.st ,((model.heat_judge.T_r_node(p,:)*sum_qout_r == h_r_all): 'T_r_node')];
    %从该节点流出的管道的首端温度与节点温度相同
    %供水管道
    for i=1:num_downpipe
        model.st = [model.st , ((model.heat_judge.T_s_node(p,:) == model.heat_judge.T_s_in(down_pipe(i,1),:)): 'T_s_node=pipe_out')];
    end
    %回水管道
    for i=1:num_uppipe
        model.st = [model.st ,(( model.heat_judge.T_r_node(p,:) == model.heat_judge.T_r_in(up_pipe(i,1),:)): 'T_r_node=pipe_out')];
    end
end
%对于热源节点来说
for t = 1:data.Horizon
     model.st = [model.st,((0 <= model.judge.H_CHP(:,t) ):'chp_st')];
     %热平衡
     model.st = model.st+[(( sum(sum(model.heat_judge.pd(t,:)))+280 >= sum(model.judge.H_CHP(:,t))+model.judge.H_GB(:,t)+...
         sum(model.judge.H_HS_out(:,t))- sum(model.judge.H_HS_in(:,t)) ...
         >= sum(sum(model.heat_judge.pd(t,:))) ):'chp_load')];
end

%热点比2.3：1
model.st = [model.st, model.judge.H_CHP == 2.3*model.judge.P_CHP*10000];
%对于热源节点，节点的供水温度等于下层供水管道的首端温度，节点回水温度等于下面连接回水管道的末端温度
 
for s=1:num_source
    p=find(pipe_fnode==source_node(s,1));
    model.heat_judge.T_s_node(source_node(s,1),:)=model.heat_judge.T_s_in(p(s),:);
    model.heat_judge.T_r_node(source_node(s,1),:)=model.heat_judge.T_r_out(p,:);
    model.st = [model.st ,(((model.judge.H_CHP(s,:) +model.judge.H_GB(s,:)+ model.judge.H_HS_out(s,:)-model.judge.H_HS_in(s,:) ) == ...
        C_water*q_pipe(source_node(s,1),1)/3600 ...
        *(model.heat_judge.T_s_node(source_node(s,1),:)-model.heat_judge.T_r_node(source_node(s,1),:))): 'H_source=c*m*delta_T')];
end

%对于负荷节点来说
 %对于负荷节点，节点的供水温度等于上层供水管道的末端温度，节点回水温度等于上层连接回水管道的首端温度
 
for l=1:num_load
    p=find(pipe_tnode==load_node(l,1));
    model.heat_judge.T_s_node(load_node(l,1),:)=model.heat_judge.T_s_out(p,:);
    model.heat_judge.T_r_node(load_node(l,1),:)=model.heat_judge.T_r_in(p,:);
    model.st = [model.st , ((model.heat_judge.pd(:,l)' == ...
        C_water*q_pipe(p,1)/3600 ...
        *(model.heat_judge.T_s_node(load_node(l,1),:)-model.heat_judge.T_r_node(load_node(l,1),:))): 'H_load=c*m*delta_T')];
end

%供回水上下限约束
model.st = [model.st , (( 65<= model.heat_judge.T_s_in <=95 ): 'min<T_s<max')];
model.st = [model.st , (( 65<= model.heat_judge.T_s_out <=95 ): 'min<T_s<max')];%70-95
model.st = [model.st , (( 35<= model.heat_judge.T_r_in <=60 ): 'min<T_r<max')];%35-60
model.st = [model.st , (( 35<= model.heat_judge.T_r_out <=60 ): 'min<T_r<max')];

%% 热负荷运行约束
N=data.heat.build_info(:,5);%3*1,
h_load_place=data.heat.build_info(:,2);%3*1,连接节点编号
C_air=data.heat.build_info(:,3).*N;%26*1,空气的比热容
R_s=data.heat.build_info(:,4)./N;%26*1,
a=exp((-delta_t*ones(26,1))./(R_s(:,1).*C_air(:,1)));
disp(a);

for t=2:data.Horizon
      model.st = [model.st ,model.heat_judge.t_in(:,t) == exp((-1*delta_t*ones(26,1))./(R_s(:,1).*C_air(:,1))).*model.heat_judge.t_in(:,t-1) + ...
          (1-exp((-1*delta_t*ones(26,1))./(R_s(:,1).*C_air(:,1)))).*(R_s(:,1).*(model.heat_judge.pd(t-1,:))'+T_out(t-1,1))];
end
 model.st = [model.st ,model.heat_judge.t_in(:,25) == exp((-1*delta_t*ones(26,1))./(R_s(:,1).*C_air(:,1))).*model.heat_judge.t_in(:,24) + ...
          (1-exp((-1*delta_t*ones(26,1))./(R_s(:,1).*C_air(:,1)))).*(R_s(:,1).*(model.heat_judge.pd(24,:))'+T_out(24,1))];
    
%用户热舒适度约束
model.st = [model.st , (( 17 <= model.heat_judge.t_in <=27 ): 'min<T_r<max')];
% model.st = [model.st , (( 21 <= model.heat_judge.t_in <=21 ): 'min<T_r<max')];
