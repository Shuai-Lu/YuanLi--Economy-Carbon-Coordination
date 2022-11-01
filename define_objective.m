function define_objective()
global data model
price_buy=data.mpc.cost(1,:);
price_sell=0.49;%0.49

model.cost.Cgrid_buy=sum((10000*model.judge.Pbuy).*price_buy);
model.cost.Cgrid=model.cost.Cgrid_buy;
%% calculate CO2
%购电部分产生的CO2
model.CO2_buy_grid=0.5810*sum(sum(10000*model.judge.Pbuy));%
%data.E_G表示六台发电机组的机组节点碳势，碳排放等于节点碳势乘以出力乘以时间
P_g=[10000*model.p_wt(1,:); 10000*model.judge.P_CHP; 10000*model.p_wt(2,:); 10000*model.p_wt(3,:); 10000*model.p_wt(4,:)];
%5*24,10是chp，其他的都是风力发电机和蓄电池
p_H=[model.judge.H_CHP;model.judge.H_GB];%1*24,6/51节点热网管道都只有一个热源点
data.E_G=([0; 443; 0; 0; 0]/1000)';
data.E_H=([443.0;516]/1000)';
model.CO2_all=0;
CO2_G=sum(sum(data.E_G*P_g));
CO2_H=sum(sum(data.E_H*p_H));
model.CO2_all=CO2_G + CO2_H +model.CO2_buy_grid;
model.CO2_e=CO2_G +model.CO2_buy_grid;
%%%%%%%%%%%
model.CO2_chp=CO2_G;
model.CO2_grid=model.CO2_buy_grid;
model.CO2_GB=sum(sum(model.judge.H_GB*0.516));
model.CO2_CHP_all=sum(sum(model.judge.H_CHP*0.443+10000*model.judge.P_CHP*0.443));

%% calculat cost
model.cost.CHP.coal = 0;%CHP机组煤耗成本
c_gas=0.575;%0.571
model.cost.CHP.coal=c_gas*(sum(model.judge.P_CHP*10000)/0.3+sum(model.judge.H_GB)/0.99);
model.cost.CHP.total = model.cost.CHP.coal ;%CHP机组总成本
%% IES运行成本
model.cost.IES_run = 0;%运行成本
K_ESS=0.01685;%运行维护费用
K_CHP=0.0480;%CHP机组运行维护
K_WT=0.02;%风机运行维护费用
K_TST=0.013;
K_GB=0.0457;
model.cost.IES_run = ...
    (K_WT*(sum(sum(model.p_wt*10000))) +...
    K_CHP*(sum(sum(model.judge.P_CHP*10000))) + K_CHP*(sum(sum(model.judge.H_CHP)))+...
    K_ESS*(sum(sum(abs(model.judge.P_BT_charge*10000)))+sum(sum(abs(model.judge.P_BT_discharge*10000))))+...
    K_TST*(sum(sum(abs(model.judge.H_HS_in)))+sum(sum(abs(model.judge.H_HS_out)))) +...
    +K_GB*(sum(sum(model.judge.H_GB))));

model.cost_obj=model.cost.IES_run+model.cost.CHP.total+model.cost.Cgrid;
model.CO2_obj=model.CO2_all;
model.obj1=model.cost_obj;
model.obj2=model.CO2_obj;

%% 探究CHP局部碳减排对整体的影响
% %CHP
% m1=20;
% F2_new=zeros(20,1);
% F1_new=zeros(20,1);
% delta_F2=zeros(20,1);delta_F1=zeros(20,1);
% delta_CHP=zeros(20,1);
% TCHP=[66160.1645700000;68976.7318400000;71645.5313900000;74103.2038300000;76521.0783900000;78679.2305400000;80742.9631200000;82804.0476000000;84911.6212500000;86869.0288500000;88655.3394400000;90413.9354200000;92105.2402400000;93780.6945400000;95433.8139400000;97083.8130600000;98890.4039000000;100980.876600000;102750.906500000;104150.187000000];
% F1=[182397.865916657;182524.454718706;182675.316393106;182843.134267172;183027.782553574;183228.697592839;183438.388926405;183656.462514397;183880.874738848;184110.915454502;184345.403398763;184584.078607678;184827.138496165;185074.596797007;185326.364934565;185582.768660740;185844.357913700;186129.989065402;186480.836239024;186930.763229854];
% F2=[141803.691569005;141317.147058091;140864.257716798;140434.915506760;140028.866895398;139645.413181390;139274.106101138;138914.439585869;138563.569474724;138220.510988072;137883.603197424;137552.507168951;137227.498830555;136908.614797601;136595.665226767;136289.151004832;135989.846523783;135723.884964778;135548.397844598;135513.839067998];
% %循环求CHP的CLCER
% for j=1:m1
%     model.st_New_CHP1= [model.st, model.CO2_CHP_all ==TCHP(j,1)-200];%CHP源节点10的碳排放约束
%     model.st_New_CHP1= [model.st_New_CHP1, model.obj2== F2(j,1)];
%     model.ops=sdpsettings('verbose', 0, 'solver', 'cplex');
%     model.sol=optimize(model.st_New_CHP1,model.obj1,model.ops);
%     F1_new(j,1)=value(model.obj1);
%     delta_F1(j,1)=F1_new(j,1)-F1(j,1);
%     delta_CHP(j,1)=200;
%     if model.sol.problem == 0
%         disp('succcessful solved');disp(j);
%     else
%         disp('error');disp(j);
%         yalmiperror(model.sol.problem)
%     end
% end
% CLCER_CHP=delta_F1./delta_CHP;
% %循环求CHP的SLCER
% for j=1:m1
%     model.st_New_CHP2= [model.st, model.CO2_CHP_all ==TCHP(j,1)-200];%CHP源节点10的碳排放约束
%     model.st_New_CHP2= [model.st_New_CHP2, model.obj1== F1(j,1)];
%     model.ops=sdpsettings('verbose', 0, 'solver', 'cplex');
%     model.sol=optimize(model.st_New_CHP2,model.obj2,model.ops);
%     F2_new(j,1)=value(model.obj2);
%     delta_F2(j,1)=F2_new(j,1)-F2(j,1);
%     delta_CHP(j,1)=200;
%     if model.sol.problem == 0
%         disp('succcessful solved');disp(j);
%     else
%         disp('error');disp(j);
%         yalmiperror(model.sol.problem)
%     end
% end
% SCER_CHP=delta_F2./delta_CHP;
%% 探究GB局部碳减排对整体的影响
% %GB
% m1=20;
% F2_new=zeros(20,1);F1_new=zeros(20,1);
% delta_F2=zeros(20,1);delta_F1=zeros(20,1);
% delta_GB=zeros(20,1);
% TGB=[61699.5287800000;59410.2323600000;57243.6428100000;55248.3810600000;53294.7231400000;51542.0055800000;49862.4904900000;48185.0340100000;46474.2228200000;44885.6794100000;43435.5165300000;42008.6072800000;40636.0495000000;39276.5557900000;37934.4916300000;36594.9716600000;35127.6683500000;33428.4465200000;32050.6095000000;31094.1468600000];
% F1=[182397.865916657;182524.454718706;182675.316393106;182843.134267172;183027.782553574;183228.697592839;183438.388926405;183656.462514397;183880.874738848;184110.915454502;184345.403398763;184584.078607678;184827.138496165;185074.596797007;185326.364934565;185582.768660740;185844.357913700;186129.989065402;186480.836239024;186930.763229854];
% F2=[141803.691569005;141317.147058091;140864.257716798;140434.915506760;140028.866895398;139645.413181390;139274.106101138;138914.439585869;138563.569474724;138220.510988072;137883.603197424;137552.507168951;137227.498830555;136908.614797601;136595.665226767;136289.151004832;135989.846523783;135723.884964778;135548.397844598;135513.839067998];
% %循环求GB的CLCER
% for j=1:m1
%     model.st_New_GB1= [model.st, model.CO2_GB ==TGB(j,1)-100];%CHP源节点10的碳排放约束
%     model.st_New_GB1= [model.st_New_GB1, model.obj2== F2(j,1)];
%     model.ops=sdpsettings('verbose', 0, 'solver', 'cplex');
%     model.sol=optimize(model.st_New_GB1,model.obj1,model.ops);
%     F1_new(j,1)=value(model.obj1);
%     delta_F1(j,1)=F1_new(j,1)-F1(j,1);
%     delta_GB(j,1)=100;
%     if model.sol.problem == 0
%         disp('succcessful solved');disp(j);
%     else
%         disp('error');disp(j);
%         yalmiperror(model.sol.problem)
%     end
% end
% CLCER_GB=delta_F1./delta_GB;
% %循环求GB的SLCER
% for j=1:m1
%     model.st_New_GB2= [model.st, model.CO2_GB ==TGB(j,1)-100];%CHP源节点10的碳排放约束
%     model.st_New_GB2= [model.st_New_GB2, model.obj1== F1(j,1)];
%     model.ops=sdpsettings('verbose', 0, 'solver', 'cplex');
%     model.sol=optimize(model.st_New_GB2,model.obj2,model.ops);
%     F2_new(j,1)=value(model.obj2);
%     delta_F2(j,1)=F2_new(j,1)-F2(j,1);
%     delta_GB(j,1)=100;
%     if model.sol.problem == 0
%         disp('succcessful solved');disp(j);
%     else
%         disp('error');disp(j);
%         yalmiperror(model.sol.problem)
%     end
% end
% SCERP_GB=delta_F2./delta_GB;

%% 准备求Pareto
model.ops=sdpsettings('verbose', 0, 'solver', 'cplex');
model.sol=optimize(model.st,model.obj1,model.ops);
if model.sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(model.sol.problem)
end
data.G1_x1=value(model.obj1);
data.G2_x1=value(model.obj2);
model.ops=sdpsettings('verbose', 0, 'solver', 'cplex');
model.sol=optimize(model.st,model.obj2,model.ops);
data.G1_x2=value(model.obj1);
data.G2_x2=value(model.obj2);
if model.sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(model.sol.problem)
end

