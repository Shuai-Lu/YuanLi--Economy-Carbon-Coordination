function Pareto_NNC()

global data model

data.L1=data.G1_x2-data.G1_x1;
data.L2=data.G2_x1-data.G2_x2;
G1_jun=(model.cost_obj-data.G1_x1)/data.L1;
G2_jun=(model.CO2_obj-data.G2_x2)/data.L2;
m1=20;
G_jun=[G1_jun ,G2_jun];
model.objective=G2_jun;

%% 
G1=[data.G1_x1,data.G2_x1];
G2=[data.G1_x2,data.G2_x2];
% 乌托邦点
G_u=[data.G1_x1,data.G2_x2];


% 乌托邦线的方向
u1=[0,1];
u2=[1,0];
N=[1,-1];

delta=1/(m1-1);
o=1;
for i = 0:delta:1
    a1(1,o)=i;
    o=o+1;
end
m=ones(1,m1);
a2=m-a1;
for j=1:m1
    x_pj(j,:)=a1(1,j)*u2+a2(1,j)*u1;
end

model.ops=sdpsettings( 'solver', 'cplex','verbose', 0);

%% 绘制Pareto前沿
% for j=19:19
for j=1:m1
    model.st_NNC = model.st+((G_jun-x_pj(j,:))*N'<=0);
    model.sol=optimize(model.st_NNC,model.objective, model.ops);
    if model.sol.problem == 0
        disp('succcessful solved');disp(j);
    else
        disp('error');
        yalmiperror(model.sol.problem)
    end
    G=[value(G1_jun)*data.L1+data.G1_x1, value(G2_jun)*data.L2 + data.G2_x2];
    model.Pareto(j,:)=G;
end

