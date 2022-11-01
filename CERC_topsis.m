%有n个评价对象，2个评价指标（第一列为成本，升序排列；第二列是CO2排放量）
function CERC_topsis()
global model
D=model.Pareto;
[n,m]=size(D);

%%
%成本和CO2排放都是极小型，所以需要进行正向化处理
D2(:,1)=1./D(:,1);
D2(:,2)=1./D(:,2);
%对正向化后的矩阵进行标准化
Z=D2./repmat(sum(D2.*D2).^0.5,n,1);

%计算权重
wj=(sum(Z)/sum(sum(Z)));

%最优解向量
[max_P]=max(Z,[],1);
[min_P]=min(Z,[],1);
for i=1:n
    D_P(i,1)=(wj(1,1)*((max_P(1,1)-Z(i,1))^2)+wj(1,2)*((max_P(1,2)-Z(i,2))^2))^0.5;
    D_N(i,1)=(wj(1,1)*((min_P(1,1)-Z(i,1))^2)+wj(1,2)*((min_P(1,2)-Z(i,2))^2))^0.5;
end
S = D_N ./ (D_P+D_N);    % 未归一化的得分
[sorted_S,index] = sort(S ,'descend');

CO2_acost=zeros(20,1);
m=(D(2,1)-D(1,1))/(D(2,2)-D(1,2));
for i=19:-1:2
    CO2_acost(i,1)=(D(i+1,1)-D(i-1,1))/(D(i+1,2)-D(i-1,2));
end

