% V12 + ICCM有限差分暂态水网 + 双向hP/hM-SOC普通水管独立版。
% 普通水管旧mita/Taylor SCA已完全移除；水泵功率FPP-SCA与气网SOC惩罚保留。
% function [cons,W_cons,G_cons,Cost,yita,p_W,y_it]=robust_IES(yita)
%% 参数录入
tic;
clear;clc
yalmip clear;
output_dir=getenv('CEWGS_OUTPUT_DIR');
if isempty(output_dir), output_dir=pwd; end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
sv=str2double(getenv('CEWGS_INITIAL_GAS_PRESSURE_SCALE'));
if ~(isfinite(sv)&&sv>0),sv=0.3;end%低压初始状态：普通气网节点取额定压力的0.3倍
gas_soc_penalty_weight=str2double(getenv('CEWGS_GAS_SOC_PENALTY_WEIGHT'));
if ~(isfinite(gas_soc_penalty_weight)&&gas_soc_penalty_weight>=0),gas_soc_penalty_weight=0.1;end
enforce_gas_node_min_pressure_flag=strtrim(getenv('CEWGS_ENFORCE_GAS_NODE_MIN_PRESSURE'));
enforce_gas_node_min_pressure=strcmpi(enforce_gas_node_min_pressure_flag,'1')|| ...
    strcmpi(enforce_gas_node_min_pressure_flag,'true');
gas_customer_inactive_pressure_upper_flag=strtrim(getenv('CEWGS_GAS_CUSTOMER_INACTIVE_PRESSURE_UPPER'));
gas_customer_inactive_pressure_upper=strcmpi(gas_customer_inactive_pressure_upper_flag,'1')|| ...
    strcmpi(gas_customer_inactive_pressure_upper_flag,'true');
water_friction_formulation='soc_auxiliary';%独立程序固定采用hP/hM-SOC，不保留旧SCA分支
water_soc_penalty_weight=str2double(getenv('CEWGS_WATER_SOC_PENALTY_WEIGHT'));
if ~(isfinite(water_soc_penalty_weight)&&water_soc_penalty_weight>=0),water_soc_penalty_weight=0.1;end
pipe_reverse_sca_flag=strtrim(getenv('CEWGS_PIPE_REVERSE_SCA'));
enable_pipe_reverse_sca=strcmpi(pipe_reverse_sca_flag,'1')||strcmpi(pipe_reverse_sca_flag,'true');
pipe_sca_tol=str2double(getenv('CEWGS_PIPE_SCA_TOL'));
if ~(isfinite(pipe_sca_tol)&&pipe_sca_tol>=0),pipe_sca_tol=0.1;end
pipe_sca_penalty_start=str2double(getenv('CEWGS_PIPE_SCA_PENALTY_START'));
if ~(isfinite(pipe_sca_penalty_start)&&pipe_sca_penalty_start>0),pipe_sca_penalty_start=1;end
pipe_sca_reference_damping=str2double(getenv('CEWGS_PIPE_SCA_REFERENCE_DAMPING'));
if ~isfinite(pipe_sca_reference_damping),pipe_sca_reference_damping=1;end
assert(pipe_sca_reference_damping>0&&pipe_sca_reference_damping<=1, ...
    'CEWGS_PIPE_SCA_REFERENCE_DAMPING must be in (0,1].');
sw=str2double(getenv('CEWGS_INITIAL_WATER_PRESSURE_SCALE'));
if ~(isfinite(sw)&&sw>0),sw=1;end%与V12基线一致时为1
case_file=strtrim(getenv('CEWGS_CASE_FILE'));
if isempty(case_file)
    case_file=fullfile(fileparts(mfilename('fullpath')), ...
        '33_20_15systems_increased_water_storage.mat');
end
assert(exist(case_file,'file')==2,'CEWGS_CASE_FILE does not identify an existing MAT file.');
supdata=importdata(case_file);%系统数据录入；允许隔离测试覆盖案例文件而不覆盖正式基线
[main_line]=supdata.main_line;%电系统支路数据
% main_line(32:36,:)=[];%去掉联络线
[load]=supdata.load;%电系统节点数据
load(:,2)=load(:,2)*0.5;
load(:,3)=load(:,3)*0.5;
% GT同时作为黑启动根节点和实际发电设备，模型不再读取独立DG数据。
[WT]=supdata.WT;%风电参数
wind_capacity_boost_MW=str2double(getenv('CEWGS_WT_CAPACITY_BOOST_MW'));
if ~isfinite(wind_capacity_boost_MW),wind_capacity_boost_MW=0;end
wind_capacity_boost_unit=str2double(getenv('CEWGS_WT_CAPACITY_BOOST_UNIT'));
if ~isfinite(wind_capacity_boost_unit),wind_capacity_boost_unit=1;end
assert(wind_capacity_boost_MW>=0,'CEWGS_WT_CAPACITY_BOOST_MW must be nonnegative.');
assert(wind_capacity_boost_unit>=1&&wind_capacity_boost_unit<=size(WT,1)&&wind_capacity_boost_unit==floor(wind_capacity_boost_unit), ...
    'CEWGS_WT_CAPACITY_BOOST_UNIT must identify an existing wind unit.');
WT(wind_capacity_boost_unit,3)=WT(wind_capacity_boost_unit,3)+wind_capacity_boost_MW;%单变量风电容量灵敏度测试
[GT]=supdata.GT;%燃气轮机参数
% GT数据表直接使用实际额定容量(MW)和爬坡率(MW/min)，
% 不再保留旧程序的经验放大倍数。

[GDS_compessor]=supdata.GDS_compessor;%燃气压缩机参数
[GDS_node]=supdata.GDS_node;%天然气系统节点
[GDS_pipe]=supdata.GDS_pipe;%天然气系统管道
[GDS_station]=supdata.GDS_station;%天然气站
gas_node1_source_capacity_kcfh=str2double(getenv('CEWGS_GAS_NODE1_SOURCE_CAPACITY_KCFH'));
if ~(isfinite(gas_node1_source_capacity_kcfh)&&gas_node1_source_capacity_kcfh>0)
    gas_node1_source_capacity_kcfh=40;
end
gas_node1_source_station=find(GDS_station(:,2)==1);
assert(isscalar(gas_node1_source_station), ...
    'The gas system must contain exactly one source station at node 1.');
GDS_station(gas_node1_source_station,3)=gas_node1_source_capacity_kcfh;
% 压力驱动服务区间只属于用气/燃机节点。定压气源保留自身物理
% 压力边界；恒流/限流气源与其余连通节点统一采用10 bar最低压力，
% 避免流量边界与源端遗留15 bar下界共同制造虚假的Weymouth冲突。
gas_customer_min_pressure_bar=str2double(getenv('CEWGS_GAS_CUSTOMER_MIN_PRESSURE_BAR'));
if ~(isfinite(gas_customer_min_pressure_bar)&&gas_customer_min_pressure_bar>0)
    gas_customer_min_pressure_bar=10;
end
gas_flow_source_node=str2double(getenv('CEWGS_GAS_FLOW_SOURCE_NODE'));
if ~isfinite(gas_flow_source_node),gas_flow_source_node=0;end
% 定压源保留自身压力边界；恒流/限流源的压力由网络方程决定，不能再
% 使用数据表中15 bar的负荷压力门槛，否则给定10 kcf/h会与Weymouth
% 等式隐含的流量发生冲突。它与普通连通节点统一采用10 bar物理下限。
[GDS_node,gas_source_data_nodes,gas_pressure_source_data_nodes,gas_uniform_pressure_nodes]= ...
    configure_gas_pressure_boundaries(GDS_node,GDS_station,gas_flow_source_node, ...
    gas_customer_min_pressure_bar);
[GDS_storage]=supdata.GDS_storage;%天然气储气罐
[WDS_pump]=supdata.WDS_pump;%水泵参数
[WDS_load]=supdata.WDS_load;%配水系统负荷
[WDS_node]=supdata.WDS_node;%配水系统节点
[WDS_pipe]=supdata.WDS_pipe;%配水系统管道
[WDS_tank]=supdata.WDS_tank;%配水系统水箱
water_pipe_flow_capacity_text=strtrim(getenv('CEWGS_WATER_PIPE_FLOW_CAPACITY_M3S'));
if isempty(water_pipe_flow_capacity_text)
    water_pipe_flow_capacity_m3s=WDS_pipe(:,7);
else
    water_pipe_flow_capacity_m3s=sscanf( ...
        strrep(water_pipe_flow_capacity_text,',',' '),'%f');
    if isscalar(water_pipe_flow_capacity_m3s)
        water_pipe_flow_capacity_m3s=repmat( ...
            water_pipe_flow_capacity_m3s,size(WDS_pipe,1),1);
    end
    assert(numel(water_pipe_flow_capacity_m3s)==size(WDS_pipe,1), ...
        ['CEWGS_WATER_PIPE_FLOW_CAPACITY_M3S must provide either one shared ', ...
        'limit or one limit per water pipe.']);
    assert(all(isfinite(water_pipe_flow_capacity_m3s))&& ...
        all(water_pipe_flow_capacity_m3s>0), ...
        'All water-pipe flow limits must be finite and positive.');
    WDS_pipe(:,7)=min(WDS_pipe(:,7),water_pipe_flow_capacity_m3s);
end
n2=size(GT,1);%微型燃气轮机数量
n_root=n2;%GT同时作为黑启动源，孤岛根节点数等于GT数量
detT=1/4;%最小时间间隔(h)
RU=GT(:,5)*60*detT;%MW/min -> MW/时段
RD=GT(:,6)*60*detT;%MW/min -> MW/时段；数据表中下爬坡为负值
n3=size(WT,1);%风电机组数量
p1=size(WDS_pump,1);%水泵数
cp1=size(GDS_compessor,1);%天然气压缩机数
Kt=6;%开关控制上限
% T=ceil((size(load,1)-n_root)/Kt);%时间断面数（以GT作为黑启动源）
% 真实时间总时段与窗口
TT = 16;   % 4 h / 15 min = 16 个实际执行时段
test_steps=str2double(getenv('CEWGS_TEST_STEPS'));
if isfinite(test_steps)
    assert(test_steps>=1&&test_steps<=32&&test_steps==floor(test_steps));
    TT=test_steps;
end
Th = 2;    % 预测时域2；配合控制时域1的0-1变量移动阻塞保证可计算性
requested_horizon=str2double(getenv('CEWGS_ROLLING_HORIZON'));
if isfinite(requested_horizon)&&requested_horizon>=1&&requested_horizon<=TT&&requested_horizon==floor(requested_horizon)
    Th=requested_horizon;
end
executed_steps=str2double(getenv('CEWGS_EXECUTE_STEPS'));
if ~isfinite(executed_steps),executed_steps=TT;end
assert(executed_steps>=1&&executed_steps<=TT&&executed_steps==floor(executed_steps), ...
    'CEWGS_EXECUTE_STEPS must be an integer from 1 to TT.');
lookahead_start_step=str2double(getenv('CEWGS_LOOKAHEAD_START_STEP'));
if ~isfinite(lookahead_start_step),lookahead_start_step=1;end
assert(lookahead_start_step>=1&&lookahead_start_step==floor(lookahead_start_step), ...
    'CEWGS_LOOKAHEAD_START_STEP must be a positive integer.');
move_blocking_flag=strtrim(getenv('CEWGS_BINARY_MOVE_BLOCKING'));
if isempty(move_blocking_flag)
    binary_move_blocking=true;
else
    binary_move_blocking=strcmpi(move_blocking_flag,'1');
end
topology_move_blocking_flag=strtrim(getenv('CEWGS_TOPOLOGY_MOVE_BLOCKING'));
if isempty(topology_move_blocking_flag)
    topology_move_blocking=binary_move_blocking;%向后兼容：旧总开关仍包含电网拓扑
else
    topology_move_blocking=strcmpi(topology_move_blocking_flag,'1')|| ...
        strcmpi(topology_move_blocking_flag,'true');
end
flow_direction_move_blocking_flag=strtrim(getenv('CEWGS_FLOW_DIRECTION_MOVE_BLOCKING'));
if isempty(flow_direction_move_blocking_flag)
    flow_direction_move_blocking=binary_move_blocking;%向后兼容：旧总开关仍包含水/气流向
else
    flow_direction_move_blocking=strcmpi(flow_direction_move_blocking_flag,'1')|| ...
        strcmpi(flow_direction_move_blocking_flag,'true');
end
water_region_flag=strtrim(getenv('CEWGS_WATER_REGION_MOVE_BLOCKING'));
if isempty(water_region_flag)
    water_region_move_blocking=binary_move_blocking;%默认保持既有行为，供隔离测试显式释放
else
    water_region_move_blocking=strcmpi(water_region_flag,'1')||strcmpi(water_region_flag,'true');
end
gt_inactive_pressure_upper_flag=strtrim(getenv('CEWGS_GT_INACTIVE_PRESSURE_UPPER'));
if isempty(gt_inactive_pressure_upper_flag)
    gt_inactive_pressure_upper=false;%物理基线：高压表示GT可用，不反向强制GT开机
else
    gt_inactive_pressure_upper=strcmpi(gt_inactive_pressure_upper_flag,'1')|| ...
        strcmpi(gt_inactive_pressure_upper_flag,'true');
end
sca_warm_flag=strtrim(getenv('CEWGS_SCA_ROLLING_WARM_START'));
if isempty(sca_warm_flag)
    sca_rolling_warm_start=false;%默认保持既有结果；显式启用后沿用上一执行点作为下窗SCA初值
else
    sca_rolling_warm_start=strcmpi(sca_warm_flag,'1')||strcmpi(sca_warm_flag,'true');
end
sca_warm_start_mode=lower(strtrim(getenv('CEWGS_SCA_WARM_START_MODE')));
if isempty(sca_warm_start_mode)
    if sca_rolling_warm_start,sca_warm_start_mode='all';else,sca_warm_start_mode='none';end
end
assert(any(strcmp(sca_warm_start_mode,{'none','pump'})), ...
    'CEWGS_SCA_WARM_START_MODE must be none or pump in the water-SOC program.');
sca_rolling_warm_start=~strcmp(sca_warm_start_mode,'none');
sca_warm_start_weight=str2double(getenv('CEWGS_SCA_WARM_START_WEIGHT'));
if ~isfinite(sca_warm_start_weight),sca_warm_start_weight=1;end
assert(sca_warm_start_weight>=0&&sca_warm_start_weight<=1, ...
    'CEWGS_SCA_WARM_START_WEIGHT must be in [0,1].');
enforce_pump_flag=strtrim(getenv('CEWGS_ENFORCE_PUMP_PHYSICS'));
if isempty(enforce_pump_flag)
    enforce_pump_physics=true;%正式对比统一验收原始泵曲线及功率等式
else
    enforce_pump_physics=strcmpi(enforce_pump_flag,'1')||strcmpi(enforce_pump_flag,'true');
end
pump_power_sca_mode=lower(strtrim(getenv('CEWGS_PUMP_POWER_SCA_MODE')));
if isempty(pump_power_sca_mode),pump_power_sca_mode='full';end
assert(any(strcmp(pump_power_sca_mode,{'partial','full'})), ...
    'CEWGS_PUMP_POWER_SCA_MODE must be partial or full.');
pump_speed_model=lower(strtrim(getenv('CEWGS_PUMP_SPEED_MODEL')));
if isempty(pump_speed_model),pump_speed_model='projection';end
assert(any(strcmp(pump_speed_model,{'sca','projection'})), ...
    'CEWGS_PUMP_SPEED_MODEL must be sca or projection.');
capture_solver_trace_flag=strtrim(getenv('CEWGS_CAPTURE_SOLVER_TRACE'));
capture_solver_trace=strcmpi(capture_solver_trace_flag,'1')||strcmpi(capture_solver_trace_flag,'true');
sca_mip_warm_start_flag=strtrim(getenv('CEWGS_SCA_MIP_WARM_START'));
if isempty(sca_mip_warm_start_flag)
    sca_mip_warm_start=true;%仅提供上一SCA解作为MIP start，不固定变量或改变可行域
else
    sca_mip_warm_start=strcmpi(sca_mip_warm_start_flag,'1')||strcmpi(sca_mip_warm_start_flag,'true');
end
sca_integer_stabilization=lower(strtrim(getenv('CEWGS_SCA_INTEGER_STABILIZATION')));
if isempty(sca_integer_stabilization),sca_integer_stabilization='none';end
assert(any(strcmp(sca_integer_stabilization, ...
    {'none','after_first','network_after_first','electric_plan_after_first'})), ...
    ['CEWGS_SCA_INTEGER_STABILIZATION must be none, after_first, ', ...
    'network_after_first, or electric_plan_after_first.']);
load(9,2:3)=load(8,2:3);
load(8,2:3)=0;
electric_load_reduction_MW=str2double(getenv('CEWGS_ELECTRIC_LOAD_REDUCTION_MW'));
if ~isfinite(electric_load_reduction_MW),electric_load_reduction_MW=0;end
electric_load_before_reduction_MW=sum(load(:,2));
assert(electric_load_reduction_MW>=0&&electric_load_reduction_MW<electric_load_before_reduction_MW, ...
    'CEWGS_ELECTRIC_LOAD_REDUCTION_MW must be nonnegative and smaller than the total electric load.');
electric_load_scale=(electric_load_before_reduction_MW-electric_load_reduction_MW)/electric_load_before_reduction_MW;
load(:,2:3)=electric_load_scale*load(:,2:3);%各节点P/Q同比例缩放，保持原功率因数与负荷空间分布
electric_load_after_reduction_MW=sum(load(:,2));
p_rate=[load(:,2),zeros(size(load(:,2),1),TT-1)]*ones(TT,TT);%用电用户有功需求
q_rate=[load(:,3),zeros(size(load(:,3),1),TT-1)]*ones(TT,TT);%用电用户无功需求
sto=size(GDS_storage,1);%储气罐数量

%% 电系统关键决策变量
l0=size(load,1);%节点数目
k0=0;%k0 为系统中的用户数目
A1=zeros(size(load,1),length(load(load(:,2)~=0)));
for i=1:size(load,1)
    if load(i,2)~=0
        k0=k0+1;
        A1(i,k0)=1;
    end
end
%% 气系统关键决策变量
ng1=size(GDS_node,1);%气网节点数
l1=size(GDS_pipe,1);%气网管道数
s1=size(GDS_station,1);%气源数
%负荷节点关联矩阵
A2=eye(ng1);
A2(GDS_node(:,6)==0,:)=[];
kg=size(A2,1);%气系统用户数独立于电/水用户数，避免数据列切换时维数隐式耦合
kw=size(WDS_node,1)-1;%水系统用户对应除水源外的节点
assert(kw==k0,'Water and electric customer counts must match the shared priority vector.');
line=main_line;
%% 计算所需初值
Mf=length(load)*5;%虚拟潮流大M
v_imin=ones(size(load,1),TT)*144.65;%%最小电压
v_imax=ones(size(load,1),TT)*176.71;
S_ijmax=repmat(main_line(:,7),1,TT);%线路容量
Mv=176.71+4*max(main_line(:,7))*max(main_line(:,4)+main_line(:,5));%线路容量大M值
% 压缩机功率因数
fi=repmat(tan(acos(GDS_compessor(:,6))),1,TT);%功率因数换算为Q/P
% 燃气轮机发电效率
sigma=repmat(GDS_compessor(:,5),1,TT);
fa=repmat(tan(acos(WDS_pump(:,5))),1,TT);%功率因数换算为Q/P
Pp_max=repmat(WDS_pump(:,7),1,TT);%水泵功率上限
Pcomp_max=repmat(GDS_compessor(:,10),1,TT);%压缩机功率上限
%%......................天然气系统初值...............................%%
% 天然气管道损失系数kcf/bar
Cij=repmat(GDS_pipe(:,4),1,TT);
% 天然气压缩机最大压缩率：直接读取数据表，不再附加经验修正。
gama=repmat(GDS_compessor(:,7),1,TT);
% 天然气管道容量上限
G_ijmax=repmat(GDS_pipe(:,6),1,TT);
% 天然气节点压力限制：最低/最高压力均直接读取数据表；满供压力
% 在下文统一定义为最低压力的3倍。
kesai_imin=repmat(GDS_node(:,3),1,TT);
kesai_imax=repmat(GDS_node(:,4),1,TT);
assert(all(GDS_node(:,4)>=3*GDS_node(:,3)), ...
    'Each gas-node maximum pressure must cover three times its minimum pressure.');
G_smax=GDS_station(:,3);%气源流量上限
%管线包存气量系数
K_mn=repmat(GDS_pipe(:,7),1,TT);
fai=0.279;%天然气电转化系数279kW/kcf
effi=GT(:,7);%燃气轮机效率
% Sl_tmin=repmat(Sl_t0,1,TT);%气罐最低储量
Sl_tmin=zeros(sto,TT);%气罐最低储量
Sl_tmax=repmat(GDS_storage(:,3),1,TT);%气罐最大储量
g_soutmax=repmat(GDS_storage(:,6),1,TT);%流出流量上限
g_sinmax=repmat(GDS_storage(:,5),1,TT);%流入流量上限
%% 储气量
% Sl_t0=zeros(sto,1);%初始储气罐存量
%%......................水系统初值参数....................................%%
pipe=WDS_pipe;
W_imax=repmat(WDS_load(:,2),1,TT);%用户水负荷需求
R_w=repmat(WDS_pipe(:,8),1,TT);% 管道阻尼参数
alpha=repmat(WDS_pump(:,2),1,TT);% 水泵参数1
beta=repmat(WDS_pump(:,3),1,TT);% 水泵参数2
W_ijmax=repmat(pipe(:,7),1,TT);%管道流量上限m^3/s
water_source_node=1;
water_source_max=sum(pipe(pipe(:,2)==water_source_node|pipe(:,3)==water_source_node,7));%水源固定上限取源节点相连管道容量
H_imin=repmat(WDS_node(:,3),1,TT);%节点水头下限
H_imax=repmat(WDS_node(:,4),1,TT);%节点水头上限
eff=repmat(WDS_pump(:,4),1,TT);%水泵效率
g=9.8;%重力加速度
Sq=5;%水箱底面积m^2
ramp_ta=repmat(WDS_tank(:,4),1,TT)/3600;%水源流量爬坡m^3/s^2
Tm=3600*detT;%时间段长度s
K_ij=repmat(WDS_pipe(:,9),1,TT);%水管存压缩系数
M_ij=repmat(WDS_pipe(:,10),1,TT);%水锤压力系数
%%......................配电系统矩阵参数计算.................................%%
A=zeros(size(load,1),size(main_line,1));
B=zeros(size(load,1),size(main_line,1));
for i=1:size(main_line,1)
    B(main_line(i,2),i)=1;
    A(main_line(i,3),i)=1;
end
R=[main_line(:,4),zeros(size(main_line(:,4),1),TT-1)]*ones(TT,TT);%电阻矩阵
X=[main_line(:,5),zeros(size(main_line(:,5),1),TT-1)]*ones(TT,TT);%电抗矩阵
%燃气轮机作为黑启动源关联矩阵
C=zeros(size(load,1),1);
k=1;
for i=1:n_root
    C(GT(i,2),k)=1;
    k=k+1;
end
%非根节点位置关联矩阵
D=eye(size(load,1));
D(:,GT(:,2))=[];
%微型燃气轮机关联矩阵
E=zeros(size(load,1),1);
k=1;
for i=1:n2
    E(GT(i,2),k)=1;
    k=k+1;
end
%风电机组关联矩阵
F=zeros(size(load,1),1);
k=1;
for i=1:n3
    F(WT(i,2),k)=1;
    k=k+1;
end
%天然气压缩机关联矩阵
EtG_c=zeros(size(load,1),cp1);
k=1;
for i=1:cp1
    EtG_c(GDS_compessor(i,9),k)=1;
    k=k+1;
end
%水泵关联矩阵
EtW_p=zeros(size(load,1),p1);
k=1;
for i=1:p1
    EtW_p(WDS_pump(i,6),k)=1;
    k=k+1;
end
%%.....................配气系统矩阵参数计算................................%%
Ag1=zeros(ng1,l1);%首节点管道关联矩阵
for i=1:l1
    Ag1(GDS_pipe(i,2),i)=1;
end
Ag2=zeros(ng1,l1);%末节点管道关联矩阵
for i=1:l1
    Ag2(GDS_pipe(i,3),i)=1;
end
G_imax=repmat(GDS_node(:,6),1,TT);%节点天然气需求
% 压缩机管道关联矩阵
Bg1=zeros(cp1,l1);
for i=1:cp1
    Bg1(i,GDS_compessor(i,1))=1;
end
% 普通管道关联矩阵
Bg2=eye(l1);
Bg2(GDS_compessor(:,1),:)=[];
%微型燃气轮机关联矩阵
E=zeros(size(load,1),1);
k=1;
for i=1:n2
    E(GT(i,2),k)=1;
    k=k+1;
end

%微型燃气轮机节点关联矩阵
Cg=zeros(n2,ng1);
for i=1:n2
    Cg(i,GT(i,3))=1;
end
%存储节点关联矩阵
St1=zeros(sto,ng1);
for i=1:sto
    St1(i,GDS_storage(i,2))=1;
end
M1=Bg1*Ag1'*kesai_imax*10;%天然气系统压强大M
%% 初始气管存
% X_mn0=Bg2*GDS_pipe(:,9);
% X_mn0=zeros(l1-cp1,1);
kesai_ireq=3*kesai_imin;
Dg=zeros(s1,ng1);% 气源节点关联矩阵
for i=1:s1
    Dg(i,GDS_station(i,2))=1;
end
gas_source_nodes=GDS_station(:,2);
gas_source_pressure_bar=str2double(getenv('CEWGS_GAS_SOURCE_PRESSURE_BAR'));
if ~(isfinite(gas_source_pressure_bar)&&gas_source_pressure_bar>0)
    gas_source_pressure_bar=20;
end
gas_flow_source_rate_kcfh=str2double(getenv('CEWGS_GAS_FLOW_SOURCE_RATE_KCFH'));
if ~isfinite(gas_flow_source_rate_kcfh),gas_flow_source_rate_kcfh=0;end
gas_flow_source_mode=lower(strtrim(getenv('CEWGS_GAS_FLOW_SOURCE_MODE')));
if isempty(gas_flow_source_mode),gas_flow_source_mode='fixed';end
assert(any(strcmp(gas_flow_source_mode,{'fixed','upper_bound'})), ...
    'CEWGS_GAS_FLOW_SOURCE_MODE must be fixed or upper_bound.');
if gas_flow_source_node>0
    gas_flow_source_station=find(gas_source_nodes==gas_flow_source_node);
    assert(isscalar(gas_flow_source_station), ...
        'The constant-flow gas-source node must identify exactly one station.');
else
    gas_flow_source_station=[];
end
gas_pressure_source_stations=setdiff((1:s1).',gas_flow_source_station);
gas_pressure_source_nodes=gas_source_nodes(gas_pressure_source_stations);
gas_source_pressure=gas_source_pressure_bar*ones(numel(gas_pressure_source_stations),1);%定压气源边界
assert(all(gas_source_pressure<=GDS_node(gas_pressure_source_nodes,4)), ...
    'Fixed gas-source pressure exceeds the source-node maximum pressure.');
if ~isempty(gas_flow_source_station)
    assert(gas_flow_source_rate_kcfh>=0&& ...
        gas_flow_source_rate_kcfh<=G_smax(gas_flow_source_station), ...
        'Constant gas-source flow must lie within the station capacity.');
end
%%.....................配水系统矩阵参数计算................................%%
Dw=zeros(size(load,1),size(WDS_pump,1));%水泵电系统耦合关联矩阵
for i=1:size(WDS_pump,1)
    Dw(WDS_pump(i,6),i)=1;
end
Tw=zeros(size(WDS_tank,1),size(WDS_node,1));%配水系统水箱节点关联矩阵
for i=1:size(WDS_tank,1)
    Tw(i,WDS_tank(i,2))=1;
end
Aw1=zeros(size(WDS_node,1),size(WDS_pipe,1));%配水系统首节点管道关联矩阵
Aw2=zeros(size(WDS_node,1),size(WDS_pipe,1));%配水系统末节点管道关联矩阵
pipe=WDS_pipe;
for i=1:size(pipe,1)
    Aw1(pipe(i,2),i)=1;
    Aw2(pipe(i,3),i)=1;
end
Aw=Aw2-Aw1;
% 水泵配水系统耦合关联矩阵
Bw=zeros(size(WDS_pump,1),size(WDS_pipe,1));
for i=1:size(WDS_pump,1)
    Bw(i,WDS_pump(i,1))=1;
end
Cw=eye(size(WDS_pipe,1),size(WDS_pipe,1));
Cw(WDS_pump(:,1),:)=[];% 管道关联矩阵
M=2500;%水头大M
source_head_value=str2double(getenv('CEWGS_WATER_SOURCE_HEAD'));
if ~isfinite(source_head_value)
    source_head_value=20.6;
end
assert(source_head_value>0,'CEWGS_WATER_SOURCE_HEAD must be positive.');
H_source_initial=source_head_value;%各节点地理高程假设为0；默认20.6 m，可用环境变量做隔离灵敏度测试
H_source_max=H_source_initial;
source_head_flag=strtrim(getenv('CEWGS_FIX_WATER_SOURCE_HEAD'));
if isempty(source_head_flag)
    fix_water_source_head=true;%经16时段测试确认恒压20.6 m可行，作为共同算例边界
else
    fix_water_source_head=strcmpi(source_head_flag,'true')||strcmp(source_head_flag,'1');
end
service_model=lower(strtrim(getenv('CEWGS_SERVICE_MODEL')));
if isempty(service_model),service_model='pressure_driven';end
assert(any(strcmp(service_model,{'pressure_driven','fixed_pressure'})), ...
    'CEWGS_SERVICE_MODEL must be pressure_driven or fixed_pressure.');
storage_model=lower(strtrim(getenv('CEWGS_STORAGE_MODEL')));
if isempty(storage_model),storage_model='dynamic_storage';end
assert(any(strcmp(storage_model,{'dynamic_storage','no_storage'})), ...
    'CEWGS_STORAGE_MODEL must be dynamic_storage or no_storage.');
water_tank_model=lower(strtrim(getenv('CEWGS_WATER_TANK_MODEL')));
if isempty(water_tank_model),water_tank_model=storage_model;end
assert(any(strcmp(water_tank_model,{'dynamic_storage','no_storage'})), ...
    'CEWGS_WATER_TANK_MODEL must be dynamic_storage or no_storage.');
gas_tank_model=lower(strtrim(getenv('CEWGS_GAS_TANK_MODEL')));
if isempty(gas_tank_model),gas_tank_model=storage_model;end
assert(any(strcmp(gas_tank_model,{'dynamic_storage','no_storage'})), ...
    'CEWGS_GAS_TANK_MODEL must be dynamic_storage or no_storage.');
water_network_model=lower(strtrim(getenv('CEWGS_WATER_NETWORK_MODEL')));
if isempty(water_network_model)
    if strcmp(storage_model,'no_storage'),water_network_model='steady_state';else,water_network_model='transient';end
end
assert(any(strcmp(water_network_model,{'transient','steady_state'})), ...
    'CEWGS_WATER_NETWORK_MODEL must be transient or steady_state.');
gas_network_model=lower(strtrim(getenv('CEWGS_GAS_NETWORK_MODEL')));
if isempty(gas_network_model)
    if strcmp(storage_model,'no_storage'),gas_network_model='steady_state';else,gas_network_model='transient';end
end
assert(any(strcmp(gas_network_model,{'transient','steady_state'})), ...
    'CEWGS_GAS_NETWORK_MODEL must be transient or steady_state.');
fixed_pressure_formulation=lower(strtrim(getenv('CEWGS_FIXED_PRESSURE_FORMULATION')));
if isempty(fixed_pressure_formulation),fixed_pressure_formulation='legacy';end
assert(any(strcmp(fixed_pressure_formulation,{'legacy','compact_exact'})), ...
    'CEWGS_FIXED_PRESSURE_FORMULATION must be legacy or compact_exact.');
pressure_driven_formulation=lower(strtrim(getenv('CEWGS_PRESSURE_DRIVEN_FORMULATION')));
if isempty(pressure_driven_formulation),pressure_driven_formulation='legacy';end
assert(any(strcmp(pressure_driven_formulation,{'legacy','compact_exact'})), ...
    'CEWGS_PRESSURE_DRIVEN_FORMULATION must be legacy or compact_exact.');
topology_formulation=lower(strtrim(getenv('CEWGS_TOPOLOGY_FORMULATION')));
if isempty(topology_formulation),topology_formulation='legacy';end
assert(any(strcmp(topology_formulation,{'legacy','strengthened'})), ...
    'CEWGS_TOPOLOGY_FORMULATION must be legacy or strengthened.');
topology_solution_mode=lower(strtrim(getenv('CEWGS_TOPOLOGY_SOLUTION_MODE')));
if isempty(topology_solution_mode),topology_solution_mode='direct';end
assert(any(strcmp(topology_solution_mode,{'direct','relax_project_polish','mcb_prototype'})), ...
    'CEWGS_TOPOLOGY_SOLUTION_MODE must be direct, relax_project_polish, or mcb_prototype.');
rpp_relax_time_limit=str2double(getenv('CEWGS_RPP_RELAX_TIME_LIMIT'));
if ~(isfinite(rpp_relax_time_limit)&&rpp_relax_time_limit>0),rpp_relax_time_limit=60;end
rpp_projection_time_limit=str2double(getenv('CEWGS_RPP_PROJECTION_TIME_LIMIT'));
if ~(isfinite(rpp_projection_time_limit)&&rpp_projection_time_limit>0),rpp_projection_time_limit=10;end
rpp_polish_time_limit=str2double(getenv('CEWGS_RPP_POLISH_TIME_LIMIT'));
if ~(isfinite(rpp_polish_time_limit)&&rpp_polish_time_limit>0),rpp_polish_time_limit=60;end
mcb_max_iterations=str2double(getenv('CEWGS_MCB_MAX_ITERATIONS'));
if ~(isfinite(mcb_max_iterations)&&mcb_max_iterations>=1&&mcb_max_iterations==floor(mcb_max_iterations)),mcb_max_iterations=3;end
mcb_relative_gap=str2double(getenv('CEWGS_MCB_RELATIVE_GAP'));
if ~(isfinite(mcb_relative_gap)&&mcb_relative_gap>=0&&mcb_relative_gap<1),mcb_relative_gap=0.05;end
mcb_relax_time_limit=str2double(getenv('CEWGS_MCB_RELAX_TIME_LIMIT'));
if ~(isfinite(mcb_relax_time_limit)&&mcb_relax_time_limit>0),mcb_relax_time_limit=30;end
mcb_integer_time_limit=str2double(getenv('CEWGS_MCB_INTEGER_TIME_LIMIT'));
if ~(isfinite(mcb_integer_time_limit)&&mcb_integer_time_limit>0),mcb_integer_time_limit=30;end
if strcmp(topology_formulation,'strengthened')
    Mf=length(load)-n_root;%多根森林中任一支路最多承载全部非根节点的单位商品流
    topology_cycle_matrix=restoration_cycle_matrix(main_line(:,2:3));
else
    topology_cycle_matrix=zeros(0,size(main_line,1));
end
V_max=WDS_tank(:,3);%水箱中的最大水量
omiga=[2,2,100,10,10,2,100,2,100,2,2,10,10,2];%负荷优先级权重
if kg==numel(omiga)
    omiga_g=omiga;%正式14用户案例保持原优先级不变
else
    omiga_g=(sum(omiga)/kg)*ones(1,kg);%列选择诊断：气用户等权且保持系统总权重不变
end

% 可选的单因素诊断：沿用已执行的基线0-1拓扑，只比较连续调度对风机增容的响应。
% 默认关闭，不参与正式恢复优化。
diagnostic_plan_file=strtrim(getenv('CEWGS_DIAGNOSTIC_PLAN_FILE'));
diagnostic_plan_active=~isempty(diagnostic_plan_file);
diagnostic_fix_scope=lower(strtrim(getenv('CEWGS_DIAGNOSTIC_FIX_SCOPE')));
if isempty(diagnostic_fix_scope),diagnostic_fix_scope='topology_and_pumps';end
assert(any(strcmp(diagnostic_fix_scope,{'topology','topology_and_pumps'})), ...
    'CEWGS_DIAGNOSTIC_FIX_SCOPE must be topology or topology_and_pumps.');
if diagnostic_plan_active
    assert(isfile(diagnostic_plan_file),'CEWGS_DIAGNOSTIC_PLAN_FILE does not exist.');
    if strcmp(diagnostic_fix_scope,'topology')
        diagnostic_plan=builtin('load',diagnostic_plan_file,'y_it_r','a_ij_r');
        assert(size(diagnostic_plan.y_it_r,2)>=TT&&size(diagnostic_plan.a_ij_r,2)>=TT, ...
            'Diagnostic topology plan is shorter than the requested trajectory.');
    else
        diagnostic_plan=builtin('load',diagnostic_plan_file,'y_it_r','a_ij_r','Xm_r');
        assert(size(diagnostic_plan.y_it_r,2)>=TT&&size(diagnostic_plan.a_ij_r,2)>=TT&&size(diagnostic_plan.Xm_r,2)>=TT, ...
            'Diagnostic baseline plan is shorter than the requested trajectory.');
    end
end

force_full_water_step=str2double(getenv('CEWGS_FORCE_FULL_WATER_STEP'));
if ~isfinite(force_full_water_step),force_full_water_step=0;end
assert(force_full_water_step==0||(force_full_water_step>=1&&force_full_water_step<=TT&&force_full_water_step==round(force_full_water_step)), ...
    'CEWGS_FORCE_FULL_WATER_STEP must be 0 or an integer within the simulated trajectory.');
force_water_node_text=strtrim(getenv('CEWGS_FORCE_FULL_WATER_NODES'));
if isempty(force_water_node_text)
    force_water_node_ids=(2:size(WDS_node,1))';
else
    force_water_tokens=regexp(force_water_node_text,'[,;\s]+','split');
    force_water_node_ids=str2double(force_water_tokens(:));
    assert(all(isfinite(force_water_node_ids))&&all(force_water_node_ids==round(force_water_node_ids))&& ...
        all(force_water_node_ids>=2)&&all(force_water_node_ids<=size(WDS_node,1)), ...
        'CEWGS_FORCE_FULL_WATER_NODES must contain water-node IDs from 2 to the final node.');
end
switch_action_cost=str2double(getenv('CEWGS_SWITCH_ACTION_COST'));
if ~isfinite(switch_action_cost),switch_action_cost=0;end
assert(switch_action_cost>=0,'CEWGS_SWITCH_ACTION_COST must be nonnegative.');
tie_switch_multiplier=str2double(getenv('CEWGS_TIE_SWITCH_MULTIPLIER'));
if ~isfinite(tie_switch_multiplier),tie_switch_multiplier=2;end
assert(tie_switch_multiplier>=1,'CEWGS_TIE_SWITCH_MULTIPLIER must be at least 1.');
switch_action_weight=ones(size(main_line,1),1);
tie_line_rows=max(1,size(main_line,1)-4):size(main_line,1);%本算例最后5条为联络线（原程序亦按32:36识别）
switch_action_weight(tie_line_rows)=tie_switch_multiplier;

%% 滚动优化所需初值
g_i0_0=zeros(kg,1);
w_i0_0=zeros(kw,1);
p_G0=zeros(n2,1);%微型燃气轮机初始功率
H_i0=3*H_imin(:,1)*sw;%初始节点压力
H_i0(1)=H_source_initial;%水源节点初始压力水头
kesai_i0=kesai_ireq(:,1)*sv;%普通气网节点采用低压系数构造初始压力
kesai_i0(gas_pressure_source_nodes)=gas_source_pressure;%仅定压气源覆盖初始压力；恒流气源保留网络初压
Sl_t0=kesai_i0(GDS_storage(:,2))*3.2/1.64/100;%初始储气罐存量与节点初压一致
X_mn0=Bg2*(K_mn(:,1).*((Ag1'+Ag2')*kesai_i0)/2);%初始气管存采用修正后的端点平均压力
X_ij0=Cw*(K_ij(:,1).*(Aw1'*H_i0+Aw2'*H_i0)./2);%初始水管存，与瞬态状态方程一致
H_0=3*H_imin(WDS_tank(:,2),1)*sw;%水箱中的初始水压
W_ij0=zeros(size(WDS_pipe,1),1);%初始管道首端流量
W_ji0=zeros(size(WDS_pipe,1),1);%初始管道末端反向流量
y_it0=load(:,5);%初始节点接入状态
a_ij0=zeros(size(main_line,1),1);%初始线路连接状态
sca_flow_seed=[];%上一执行时段的普通管平均流量/泵流量，仅用于SCA线性化初值
sca_speed_seed=[];%上一执行时段的泵速，仅用于SCA线性化初值
gas_flow_seed=zeros(l1-cp1,1);%气管反向SCA的上一执行平均流量参考点
gas_low_pressure_seed=min(Bg2*Ag1'*kesai_i0,Bg2*Ag2'*kesai_i0);%气管反向SCA低压端参考点
%% 变量输出
    w_i0_r=zeros(kw,TT);H_tank_r=zeros(size(WDS_tank,1),TT);P_pump_r=zeros(p1,TT);omima_r=zeros(p1,TT);
    y_it_r=zeros(size(load,1),TT);r_i_r=zeros(k0,TT);g_i0_r=zeros(kg,TT);W_WR_r=zeros(1,TT);V_t_r=zeros(size(WDS_tank,1),TT);
    W_ij_r=zeros(size(WDS_pipe,1),TT);H_i_r=zeros(size(WDS_node,1),TT);L_ij_r=zeros(size(WDS_pipe,1)-p1,TT);
    W_ji_r=zeros(size(WDS_pipe,1),TT);X_ij_r=zeros(size(WDS_pipe,1)-p1,TT);W_ijave_r=zeros(size(WDS_pipe,1),TT);
    P_comp_r=zeros(cp1,TT);
    Sl_t_r=zeros(sto,TT);
    X_mn_r=zeros(l1-cp1,TT);
    p_G_r=zeros(n2,TT);
    q_G_r=zeros(n2,TT);
    g_GS_r=zeros(s1,TT);
    p_W_r=zeros(n3,TT);
    a_ij_r=zeros(size(main_line,1),TT);
    F_ij_r=zeros(size(main_line,1),TT);
    v_i_r=zeros(size(load,1),TT);
    g_GT_r=zeros(n2,TT);
    kesai_i_r=zeros(ng1,TT);
    kesai_ijH_r=zeros(l1-cp1,TT);
    kesai_ijL_r=zeros(l1-cp1,TT);
    Xc_r=zeros(cp1,TT);
    Xm_r=zeros(size(WDS_pump,1),TT);
    I_ij_r=zeros(l1-cp1,TT);
    G_ij_r=zeros(l1,TT);
    G_ji_r=zeros(l1,TT);
    G_ijave_r=zeros(l1,TT);
    Gp_ij_r=zeros(size(WDS_pump,1),TT);
    g_sin_r=zeros(sto,TT); 
    g_sout_r=zeros(sto,TT); 
    q_W_r=zeros(n3,TT);
    wind_lower_r=zeros(n3,TT);wind_mean_r=zeros(n3,TT);wind_sigma_r=zeros(n3,TT);
    water_hP_r=zeros(size(WDS_pipe,1)-p1,TT);water_hM_r=zeros(size(WDS_pipe,1)-p1,TT);
    water_qP_r=zeros(size(WDS_pipe,1)-p1,TT);water_qM_r=zeros(size(WDS_pipe,1)-p1,TT);
    water_soc_gap_r=zeros(1,TT);gas_soc_gap_r=zeros(1,TT);
    water_pipe_sca_slack_r=zeros(1,TT);gas_pipe_sca_slack_r=zeros(1,TT);
    sca_iter_r=zeros(1,TT);sca_slack_r=zeros(1,TT);pump_slack_r=zeros(1,TT);solve_time_r=zeros(1,TT);sca_plan_locked_r=false(1,TT);
    sca_ref_flow_rms_r=zeros(1,TT);sca_ref_pump_rms_r=zeros(1,TT);
    solver_problem_r=zeros(1,TT);friction_res_r=zeros(1,TT);linepack_res_r=zeros(1,TT);gas_dynamic_res_r=zeros(1,TT);pump_curve_res_r=zeros(1,TT);pump_power_res_r=zeros(1,TT);switch_action_count_r=zeros(1,TT);switch_action_weighted_r=zeros(1,TT);
    sca_solve_time_trace_r=cell(1,TT);solver_output_trace_r=cell(1,TT);rpp_trace_r=cell(1,TT);mcb_trace_r=cell(1,TT);
    sca_objective_trace_r=cell(1,TT);sca_ref_flow_rms_trace_r=cell(1,TT);sca_ref_pump_rms_trace_r=cell(1,TT);
    sca_water_soc_gap_trace_r=cell(1,TT);sca_gas_soc_gap_trace_r=cell(1,TT);sca_pipe_residual_trace_r=cell(1,TT);sca_pump_slack_trace_r=cell(1,TT);
    sca_water_pipe_slack_trace_r=cell(1,TT);sca_gas_pipe_slack_trace_r=cell(1,TT);
    sca_pump_curve_residual_trace_r=cell(1,TT);sca_pump_power_residual_trace_r=cell(1,TT);
%% ============== 真实时间滚动循环 ==============
for t = 1:executed_steps
    step_tic=tic;
    if t<lookahead_start_step
        T=1;%先用已验证的单时段模型建立可行滚动状态，再启用前瞻窗口
    else
        T = min(Th, TT - t + 1);  % 本轮窗口长度
    end
    %% 策略输入
    if binary_move_blocking&&T>1
        r_i=repmat(binvar(k0,1,'full'),1,T);
    else
        r_i=binvar(k0,T,'full');
    end
    r_i0=A1*r_i;%用电负荷恢复比例
    %% 水系统关键决策变量
    % w_i0=importdata('策略.XLSX').w_i;
    %% 配电系统中间变量
    P_ij=sdpvar(size(main_line,1),T,'full');%线路传输有功
    Q_ij=sdpvar(size(main_line,1),T,'full');%线路传输无功
    v_i=sdpvar(size(load,1),T,'full');%节点电压平方
    p_i=sdpvar(size(load,1),T,'full');%节点有功
    q_i=sdpvar(size(load,1),T,'full');%节点无功
    p_G=sdpvar(n2,T,'full');%微型燃气轮机有功(电网侧)
    q_G=sdpvar(n2,T,'full');%微型燃气轮机无功(电网侧)
    p_W=sdpvar(n3,T,'full');%风电机组有功标幺
    q_W=sdpvar(n3,T,'full');%风电机组无功标幺
    P_comp=sdpvar(cp1,T,'full');%压缩机功率(电网侧)
    q_comp=sdpvar(cp1,T,'full');%压缩机无功(电网侧)
    P_pump=sdpvar(p1,T,'full');%水泵功率(电网侧)
    q_pump=sdpvar(p1,T,'full');%压缩机无功(电网侧)
    %%................................虚拟流控制变量...........................................%%
    R_i=sdpvar(n_root,T,'full');
    F_ij=sdpvar(size(main_line,1),T,'full');
    %% 配气系统中间变量
    G_ij=sdpvar(l1,T,'full');%ij方向管道流量
    G_ji=sdpvar(l1,T,'full');%ji方向管道流量
    G_ijave=sdpvar(l1,T,'full');%管道平均气流量
    g_GS=sdpvar(s1,T,'full');%气源输气量
    kesai_i=sdpvar(ng1,T,'full');%节点气压
    if binary_move_blocking&&T>1
        Xc=repmat(binvar(cp1,1,'full'),1,T);%控制时域1：窗口内保持启停状态
    else
        Xc=binvar(cp1,T,'full');%压缩机启停状态
    end
    X_mn=sdpvar(l1-cp1,T,'full');%气管存
    g_GT=sdpvar(n2,T,'full');%燃气轮机天然气需求
    g_sin=sdpvar(sto,T,'full');%流入储气罐流量
    g_sout=sdpvar(sto,T,'full');%流出储气罐流量
    Sl_t=sdpvar(sto,T,'full');%储气罐存量
    %% 配水系统中间变量
    W_ij=sdpvar(size(WDS_pipe,1),T,'full');%正向管道流量
    W_ji=sdpvar(size(WDS_pipe,1),T,'full');%反向管道流量
    X_ij=sdpvar(size(WDS_pipe,1)-size(WDS_pump,1),T,'full');%管道内存水
    W_WR=[sdpvar(1,T,'full');zeros(size(WDS_node,1)-1,T)];%水源
    Gp_ij=sdpvar(size(WDS_pump,1),T,'full');%水泵扬程
    H_i=sdpvar(size(WDS_node,1),T,'full');%节点水压
    if binary_move_blocking&&T>1
        Xm=repmat(binvar(size(WDS_pump,1),1,'full'),1,T);%控制时域1：窗口内保持启停状态
    else
        Xm=binvar(size(WDS_pump,1),T,'full');%水泵启停
    end
    V_t=sdpvar(size(WDS_tank,1),T,'full');%水箱水量
    H_tank=sdpvar(size(WDS_tank,1),T,'full');%水箱压头辅助变量
    L_ij=sdpvar(size(WDS_pipe,1)-size(WDS_pump,1),T,'full');%管道水头损失
    W_ijave=sdpvar(size(WDS_pipe,1),T+1,'full');%端点净流量的空间平均，用于ICCM后向差分


    %% 优化模型
    if strcmp(pump_speed_model,'sca')
        omima=sdpvar(size(WDS_pump,1),T,'full');%原模型：水泵转速作为SCA变量
    else
        omima=[];%投影模型：转速由流量和扬程在求解后唯一恢复，不再作为优化变量
    end
    % 两个对比算例使用相同的连续恢复比例。fixed-pressure仅改变
    % “供给是否要求额定压力”的关系，不能把客户负荷额外改成全取或全弃。
    g_i0=[g_i0_0,sdpvar(kg,T)];%配气恢复比例
    w_i0=[w_i0_0,sdpvar(kw,T)];%配水恢复比例
    fixed_pressure_g_available=[];
    fixed_pressure_w_available=[];
    pressure_driven_g_active=[];
    pressure_driven_w_active=[];
    pressure_driven_gt_active=[];
    if strcmp(service_model,'fixed_pressure')&&strcmp(fixed_pressure_formulation,'compact_exact')
        if binary_move_blocking&&T>1
            fixed_pressure_g_available=repmat(binvar(kg,1,'full'),1,T);
        else
            fixed_pressure_g_available=binvar(kg,T,'full');
        end
        if water_region_move_blocking&&T>1
            fixed_pressure_w_available=repmat(binvar(kw,1,'full'),1,T);
        else
            fixed_pressure_w_available=binvar(kw,T,'full');
        end
    end
    if strcmp(service_model,'pressure_driven')&&strcmp(pressure_driven_formulation,'compact_exact')
        if binary_move_blocking&&T>1
            pressure_driven_g_active=repmat(binvar(kg,1,'full'),1,T);
            pressure_driven_gt_active=repmat(binvar(n2,1,'full'),1,T);
        else
            pressure_driven_g_active=binvar(kg,T,'full');
            pressure_driven_gt_active=binvar(n2,T,'full');
        end
        if water_region_move_blocking&&T>1
            pressure_driven_w_active=repmat(binvar(kw,1,'full'),1,T);
        else
            pressure_driven_w_active=binvar(kw,T,'full');
        end
    end
    g_i=A2'*g_i0;
    if topology_move_blocking&&T>1
        y_it=[y_it0,repmat(binvar(size(load,1),1,'full'),1,T)];%控制时域1、预测时域2
        a_ij=[a_ij0,repmat(binvar(size(main_line,1),1,'full'),1,T)];
    else
        y_it=[y_it0,binvar(size(load,1),T,'full')];%节点切除标志
        a_ij=[a_ij0,binvar(size(main_line,1),T,'full')];%支路切除标志
    end
    V_0=H_0*Sq;%水箱中的初始水量
    cons=[];%主问题约束
    if diagnostic_plan_active
        if topology_move_blocking&&T>1
            target_y=repmat(round(diagnostic_plan.y_it_r(:,t)),1,T);
            target_a=repmat(round(diagnostic_plan.a_ij_r(:,t)),1,T);
        else
            target_cols=t:(t+T-1);
            target_y=round(diagnostic_plan.y_it_r(:,target_cols));
            target_a=round(diagnostic_plan.a_ij_r(:,target_cols));
        end
        if strcmp(diagnostic_fix_scope,'topology_and_pumps')
            if binary_move_blocking&&T>1
                target_Xm=repmat(round(diagnostic_plan.Xm_r(:,t)),1,T);
            else
                target_cols=t:(t+T-1);
                target_Xm=round(diagnostic_plan.Xm_r(:,target_cols));
            end
        end
        cons=[cons,y_it(:,2:T+1)==target_y,a_ij(:,2:T+1)==target_a];
        if strcmp(diagnostic_fix_scope,'topology_and_pumps')
            cons=[cons,Xm==target_Xm];
        end
    end
    % 受控反事实诊断：仅在指定的全局时段强制所有用水负荷满供。
    % 默认关闭；用于判断未满供是物理不可行，还是滚动/SCA轨迹与目标权衡所致。
    forced_water_col=force_full_water_step-t+1;
    if force_full_water_step>0&&forced_water_col>=1&&forced_water_col<=T
        cons=[cons,w_i0(force_water_node_ids-1,forced_water_col+1)==1];
    end
    vec=@(z)reshape(z,[],1);%强制逐元素不等式，避免方阵被YALMIP解释为半正定约束
    %% 约束条件
    %............电系统约束条件................%%
    cons=[cons,vec(-repmat(WT(:,3),1,T).*p_W.*repmat(tan(acos(WT(:,4))),1,T)+q_W)<=0];%风电机组无功上界
    cons=[cons,vec(-q_W-repmat(WT(:,3),1,T).*p_W.*repmat(tan(acos(WT(:,4))),1,T))<=0];%风电机组无功下界
    %微型燃气轮机功率约束
    if binary_move_blocking&&T>1
        u_it=repmat(binvar(3,1,'full'),1,T);%窗口内不重复启停
    else
        u_it=binvar(3,T,'full');%燃气轮机启停
    end
    cons=[cons,vec(p_G-u_it.*repmat(GT(:,4),1,T))<=0];%燃气轮机接入有功约束
    cons=[cons,vec(-p_G)<=0];%燃气轮机接入有功约束
    cons=[cons,vec(q_G-u_it.*repmat(GT(:,4),1,T)*3)<=0];%燃气轮机接入无功约束
    cons=[cons,vec(-q_G-u_it.*repmat(GT(:,4),1,T)*3)<=0];%燃气轮机接入无功约束
    cons=[cons,vec(p_G-[p_G0,p_G(:,1:T-1)])<=vec(repmat(RU,1,T))];%燃气轮机上爬坡约束
    cons=[cons,vec(-p_G+[p_G0,p_G(:,1:T-1)])<=vec(-repmat(RD,1,T))];%燃气轮机下爬坡约束
    cons=[cons,A1'*(y_it(:,2:T+1))<=r_i];%负荷接入逻辑
    %传输线线性交流约束
    cons=[cons,v_i(line(:,2),:)-v_i(line(:,3),:)<=Mv*(1-a_ij(:,2:T+1))+2*(P_ij.*R(:,1:T)+Q_ij.*X(:,1:T))];%正确方向的线路电压降
    cons=[cons,-Mv*(1-a_ij(:,2:T+1))+2*(P_ij.*R(:,1:T)+Q_ij.*X(:,1:T))<=v_i(line(:,2),:)-v_i(line(:,3),:)];%正确方向的线路电压降
    cons=[cons,-a_ij(:,2:T+1)*Mf<=F_ij];%虚拟流容量上边界
    cons=[cons,F_ij<=a_ij(:,2:T+1)*Mf];%虚拟流容量下边界
    cons=[cons,D'*(y_it(:,2:T+1).*v_imin(:,1:T))<=D'*reshape(v_i,l0,[])];%节点电压约束
    cons=[cons,D'*reshape(v_i,l0,[])<=D'*(y_it(:,2:T+1).*v_imax(:,1:T))];
    %线路容量约束：strengthened为二进制开关集合的精确凸包SOC，整数可行域不变。
    if strcmp(topology_formulation,'strengthened')
        for l=1:size(main_line,1)
            for tau=1:T
                cons=[cons,cone([P_ij(l,tau);Q_ij(l,tau)],S_ijmax(l,tau)*a_ij(l,tau+1))];
            end
        end
    else
        cons=[cons,P_ij.^2+Q_ij.^2<=S_ijmax(:,1:T).^2.*a_ij(:,2:(T+1))];%原始较弱连续松弛
    end
    % %线路容量约束(线性化处理后)
    % cons=[cons,P_ij<=S_ijmax.*a_ij(:,2:(T+1))/sqrt(2)];%非线性约束
    % cons=[cons,P_ij>=-S_ijmax.*a_ij(:,2:(T+1))/sqrt(2)];%非线性约束
    % cons=[cons,Q_ij<=S_ijmax.*a_ij(:,2:(T+1))/sqrt(2)];%非线性约束Q_ij
    % cons=[cons,Q_ij>=-S_ijmax.*a_ij(:,2:(T+1))/sqrt(2)];%非线性约束Q_ij
    %网络黑启动约束
    cons=[cons,sum(a_ij(:,2:T+1)-a_ij(:,1:T))<=ones(1,T)*Kt];%单位时间恢复线路数约束
    cons=[cons,-y_it(line(:,2),2:(T+1))<=-a_ij(line(:,1),2:(T+1))];%线路恢复两端节点一定恢复
    cons=[cons,-y_it(line(:,3),2:(T+1))<=-a_ij(line(:,1),2:(T+1))];
    cons=[cons,-a_ij(line(:,1),2:(T+1))<=-a_ij(line(:,1),1:T)];%恢复稳定性约束
    cons=[cons,vec(-R_i)<=0];%黑启动源虚拟流约束
    %功率平衡约束
    cons=[cons,A*P_ij+p_i==B*P_ij];%有功平衡
    cons=[cons,p_i==F*(repmat(WT(:,3),1,T).*p_W)+E*p_G-r_i0.*p_rate(:,t:t+T-1)-EtG_c*P_comp/1000-EtW_p*P_pump/1000];%节点有功；GT按P_bus接入
    cons=[cons,A*Q_ij+q_i==B*Q_ij];%无功平衡
    cons=[cons,q_i==F*q_W+E*q_G-r_i0.*q_rate(:,t:t+T-1)-EtG_c*q_comp/1000-EtW_p*q_pump/1000];%节点无功；GT按P_bus接入
    cons=[cons,C'*v_i==ones(n_root,T)*160.2756];%GT黑启动根节点电压
    cons=[cons,q_comp==fi(:,1:T).*P_comp];%压缩机无功
    cons=[cons,q_pump==fa(:,1:T).*P_pump];%水泵无功
    cons=[cons,vec(P_pump)>=0,vec(P_pump)<=vec(EtW_p'*y_it(:,2:end).*Pp_max(:,1:T))];%水泵负荷拾取约束
    cons=[cons,vec(EtW_p'*y_it(:,2:end))>=vec(Xm)];
    Pcomp_cap=EtG_c'*y_it(:,2:end).*Pcomp_max(:,1:T);
    cons=[cons,vec(P_comp)>=0,vec(P_comp)<=vec(Pcomp_cap)];%压缩机负荷拾取约束
    cons=[cons,vec(EtG_c'*y_it(:,2:end))>=vec(Xc)];
    %虚拟流约束
    cons=[cons,sum(a_ij(:,2:(T+1)))==size(load,1)*ones(1,T)-sum(1-y_it(:,2:(T+1)))-n_root*ones(1,T)];%以GT为根的孤岛生成树约束
    cons=[cons,D'*(B*F_ij)+D'*y_it(:,2:(T+1))==D'*(A*F_ij)];%黑启动节点连通性约束
    cons=[cons,C'*(B*F_ij-A*F_ij)==R_i];%黑启动源连通性约束
    if strcmp(topology_formulation,'strengthened')
        for cycle_index=1:size(topology_cycle_matrix,1)
            cycle_rows=find(topology_cycle_matrix(cycle_index,:));
            cons=[cons,sum(a_ij(cycle_rows,2:T+1),1)<=numel(cycle_rows)-1];%基本环割，仅强化LP松弛
        end
    end
    %%.........天然气系统约束条件....................%%
    if strcmp(gas_tank_model,'dynamic_storage')
        cons=[cons,-Ag1*G_ij-Ag2*G_ji+Dg'*g_GS-g_i(:,2:end).*G_imax(:,t:t+T-1)-Cg'*g_GT+St1'*(g_sout-g_sin)==0];%含储气罐的气体流量平衡
    else
        cons=[cons,-Ag1*G_ij-Ag2*G_ji+Dg'*g_GS-g_i(:,2:end).*G_imax(:,t:t+T-1)-Cg'*g_GT==0];%无储能稳态基线
    end
    cons=[cons,fai*g_GT.*repmat(effi,1,T)==p_G];%燃气轮机耗气量关系GT(:,4)./effi/fai
    cons=[cons,(0<=g_i0(:,2:end))&(g_i0(:,2:end)<=1)];%仅约束当前决策列，避免上一时段容差值造成假不可行
    % cons=[cons,g_i0(:,2:end)-g_i0(:,1:end-1)>=0];%恢复后不可短时中断
    cons=[cons,vec(g_GS)>=0];
    cons=[cons,vec(g_GS)<=vec(repmat(G_smax,1,T))];%气源仅受固定容量上限约束，不设源端爬坡
    if enforce_gas_node_min_pressure
        % 当前算例不含气网管道停运/节点隔离变量，全部气网节点保持连通，
        % 因而中转节点也必须满足数据表中的物理最低压力。
        cons=[cons,kesai_imin(:,1:T)<=kesai_i];
    else
        cons=[cons,0<=kesai_i];%旧结果复现：允许无负荷中转节点降至0 bar
    end
    if ~isempty(gas_pressure_source_stations)
        cons=[cons,Dg(gas_pressure_source_stations,:)*kesai_i== ...
            repmat(gas_source_pressure,1,T)];%至少保留一个定压参考气源
    end
    if ~isempty(gas_flow_source_station)
        if strcmp(gas_flow_source_mode,'fixed')
            cons=[cons,g_GS(gas_flow_source_station,:)==gas_flow_source_rate_kcfh];%正式算例：节点8恒流气源
        else
            cons=[cons,g_GS(gas_flow_source_station,:)<=gas_flow_source_rate_kcfh];%仅供稳态容量诊断：10 kcf/h为可调上限
        end
    end
    % cons=[cons,kesai_imin(GT(:,3),:).*y_it(GT(:,2),2:(T+1))<=kesai_i(GT(:,3),:)];%GT气压下限
    cons=[cons,kesai_i<=kesai_imax(:,1:T)];%节点气压上限
    cons=[cons,Bg2*G_ijave==0.5*Bg2*(G_ij-G_ji)];%管道平均流量
    % 压缩机流入和流出节点气压关系
    cons=[cons,-Bg1*G_ji==Bg1*G_ij];%压缩机流入与流出流量关系
    cons=[cons,vec(-M1(:,1:T).*(ones(cp1,T)-Xc)+(Bg1*Ag1'*kesai_i))<=vec(Bg1*Ag2'*kesai_i)];
    cons=[cons,vec(M1(:,1:T).*(ones(cp1,T)-Xc)+gama(:,1:T).*(Bg1*Ag1'*kesai_i))>=vec(Bg1*Ag2'*kesai_i)];
    cons=[cons,vec(Bg1*G_ij)<=vec(Xc.*(Bg1*G_ijmax(:,1:T)))];%压缩机传输流量约束
    cons=[cons,vec(Bg1*G_ij)>=0];
    % 气网管道模型与显式储气罐分开控制。稳态气网仅取消管线包，
    % 而 gas_tank_model 独立决定是否保留调节气罐。
    if strcmp(gas_network_model,'transient')
        cons=[cons,Bg2*(G_ij+G_ji)*detT==X_mn-[X_mn0,X_mn(:,1:end-1)]];%每条管道的气体流量传输平衡
        cons=[cons,X_mn==Bg2*(K_mn(:,1:T).*((Ag1'+Ag2')*kesai_i*0.5))];%管线包与平均节点压力成正比
        cons=[cons,X_mn>=0];
    else
        cons=[cons,Bg2*(G_ij+G_ji)==0];%无管存基线：每条普通管道两端流量瞬时守恒
        cons=[cons,X_mn==0];%保留统一输出维度，但不允许管道储气参与调度
    end
    % %,,,,,,,,,,,,,,,,,weymouth 方程(双向凸松弛)....................%
    if flow_direction_move_blocking&&T>1
        I_ij=repmat(binvar(l1-cp1,1,'full'),1,T);%30 min预测窗内气流方向保持不变
    else
        I_ij=binvar(l1-cp1,T,'full');%管道方向判别0-1变量
    end
    kesai_ijH=sdpvar(l1-cp1,T,'full');%i-j中高压节点气压
    kesai_ijL=sdpvar(l1-cp1,T,'full');%i-j中低压节点气压
    % %凸约束（线性化后）
    % theta = linspace(0,2*pi,4);
    % for k = 1:length(theta)
    %     u = [cos(theta(k)); sin(theta(k))];
    %     cons=[cons,u(1)*Bg2*(G_ijave./Cij) + u(2)*kesai_ijL <= kesai_ijH];
    % end
    cons=[cons,sqrt(Bg2*(G_ijave./Cij(:,1:T)).^2+kesai_ijL.^2)<=kesai_ijH];%凸约束（线性化前）
    cons=[cons,-I_ij.*(Bg2*G_ijmax(:,1:T))<=Bg2*G_ijave];
    cons=[cons,(1-I_ij).*(Bg2*G_ijmax(:,1:T))>=Bg2*G_ijave];
    %I_ij=0时首节点气压大于末节点
    cons=[cons,-I_ij.*(Bg2*Ag1'*kesai_imax(:,1:T))<=kesai_ijH-Bg2*(Ag1'*kesai_i)];
    cons=[cons,I_ij.*(Bg2*Ag1'*kesai_imax(:,1:T))>=kesai_ijH-Bg2*(Ag1'*kesai_i)];
    cons=[cons,-I_ij.*(Bg2*Ag2'*kesai_imax(:,1:T))<=kesai_ijL-Bg2*(Ag2'*kesai_i)];
    cons=[cons,I_ij.*(Bg2*Ag2'*kesai_imax(:,1:T))>=kesai_ijL-Bg2*(Ag2'*kesai_i)];
    %I_ij=1时首节点气压小于末节点
    cons=[cons,(I_ij-1).*(Bg2*Ag2'*kesai_imax(:,1:T))<=kesai_ijH-Bg2*(Ag2'*kesai_i)];
    cons=[cons,(1-I_ij).*(Bg2*Ag2'*kesai_imax(:,1:T))>=kesai_ijH-Bg2*(Ag2'*kesai_i)];
    cons=[cons,(I_ij-1).*(Bg2*Ag1'*kesai_imax(:,1:T))<=kesai_ijL-Bg2*(Ag1'*kesai_i)];
    cons=[cons,(1-I_ij).*(Bg2*Ag1'*kesai_imax(:,1:T))>=kesai_ijL-Bg2*(Ag1'*kesai_i)];
    %....................................................%
    % %压缩机功率相关约束
    cons=[cons,P_comp==sigma(:,1:T).*(Bg1*G_ij*28.3168)];
    %% 储气罐相关约束
    if strcmp(gas_tank_model,'dynamic_storage')
        cons=[cons,Sl_t-[Sl_t0,Sl_t(:,1:end-1)]==(g_sin-g_sout).*detT];%储气罐存量约束
        % cons=[cons,Sl_t<=Sl_tmax];%存量上限约束
        % cons=[cons,Sl_t>=Sl_tmin];%最低存量约束
        cons=[cons,St1*kesai_i==1.64*Sl_t./3.2*100];
        cons=[cons,vec(g_sout)>=0,vec(g_sout)<=vec(g_soutmax(:,1:T))];%流率上限不重复乘时间步长
        cons=[cons,vec(g_sin)>=0,vec(g_sin)<=vec(g_sinmax(:,1:T))];%流率上限不重复乘时间步长
    else
        cons=[cons,Sl_t==0,g_sin==0,g_sout==0];%无储能基线：移除调节气罐
    end
    %%.........水系统约束条件....................%%
    if strcmp(water_tank_model,'dynamic_storage')
        cons=[cons,(-Aw1*W_ij-Aw2*W_ji)*Tm+(W_WR-[zeros(1,T);w_i0(:,2:end)].*(W_imax(:,1:T)/3600))*Tm-Tw'*(V_t-[V_0,V_t(:,1:T-1)])==0];%含水箱储水的流量平衡
    else
        cons=[cons,-Aw1*W_ij-Aw2*W_ji+W_WR-[zeros(1,T);w_i0(:,2:end)].*(W_imax(:,1:T)/3600)==0];%无储能稳态水量平衡
    end
    cons=[cons,0<=W_WR(1,:),W_WR(1,:)<=water_source_max];%水源仅受固定容量上限约束，不设源端爬坡
    %配水系统负荷约束
    cons=[cons,(0<=w_i0(:,2:end))&(w_i0(:,2:end)<=1)];%仅约束当前决策列，避免上一时段容差值造成假不可行
    % cons=[cons,w_i0(:,2:end)-w_i0(:,1:end-1)>=0];%恢复后不可短时中断
    % 水管模型与显式水箱分开控制。稳态水网取消管内储水和惯性项，
    % water_tank_model 独立决定是否保留调节水箱。
    if strcmp(water_network_model,'transient')
        cons=[cons,X_ij==Cw*(K_ij(:,1:T).*(Aw1'*H_i+Aw2'*H_i)./2)];%节点水头与管道存水状态
        % ICCM：时间导数采用后向差分，空间两端水头取中心值；
        % 端点流量取当前时刻，不再对相邻两个时刻作梯形平均。
        cons=[cons,Cw*(W_ij+W_ji)*Tm==X_ij-[X_ij0,X_ij(:,1:end-1)]];%跨时段管道存水更新
        cons=[cons,X_ij>=0];
    else
        cons=[cons,Cw*(W_ij+W_ji)==0,X_ij==0];%无储能稳态基线：普通水管无储水
    end
    %水箱相关约束
    H_i0(1)=H_source_max;
    if fix_water_source_head
        cons=[cons,H_i(1,:)==H_source_initial];%恒压水源隔离测试
    else
        cons=[cons,H_i(1,:)<=H_source_max];%源端后续压力上限；不是恒压等式
    end
    if strcmp(water_tank_model,'dynamic_storage')
        cons=[cons,H_tank==V_t./Sq];
        cons=[cons,H_tank==Tw*H_i];
        cons=[cons,vec(V_t)<=vec(repmat(V_max,1,T))];
        cons=[cons,vec(V_t)>=0];
        cons=[cons,vec(H_tank)>=0];
        cons=[cons,vec(V_t-[V_0,V_t(:,1:T-1)])<=vec(ramp_ta(:,1:T)*Tm)];%水箱流速限制
        cons=[cons,vec(V_t-[V_0,V_t(:,1:T-1)])>=vec(-ramp_ta(:,1:T)*Tm)];
    else
        cons=[cons,V_t==0,H_tank==0];%无储能稳态基线：移除调节水箱
    end
    %水泵管道水压关系
    cons=[cons,Bw*Aw'*H_i==Gp_ij];
    %% 水头损失松弛约束
    if strcmp(water_network_model,'transient')
        cons=[cons,L_ij==-Cw*Aw'*H_i];%ICCM空间差分：采用当前时刻节点水头
        cons=[cons,W_ijave==([W_ij0,W_ij]-[W_ji0,W_ji])*0.5];%当前/前一时刻的管道空间平均流量
    else
        cons=[cons,L_ij==-Cw*Aw'*H_i];%稳态水头差
        cons=[cons,W_ijave(:,1)==W_ij0,W_ijave(:,2:end)==W_ij];%统一输出维度
    end
    %......................双向流量凸松弛..............................%
    if flow_direction_move_blocking&&T>1
        deta=repmat(binvar(size(WDS_pipe,1)-size(WDS_pump,1),1),1,T);%预测窗内水流方向保持不变
    else
        deta=binvar(size(WDS_pipe,1)-size(WDS_pump,1),T);%管道流向判断
    end
    if strcmp(water_network_model,'transient')
        W_ijbar=W_ijave(:,2:end);%ICCM：摩阻采用当前时刻的管道空间平均流量
        D_ij=L_ij+Cw*(M_ij(:,1:T).*(W_ijave(:,1:end-1)-W_ijave(:,2:end))/Tm);%水头差与惯性项
    else
        W_ijbar=W_ij;%稳态管流量
        D_ij=L_ij;%稳态模型无管存惯性项
    end
    cons=[cons,Cw*W_ijbar>=-Cw*W_ijmax(:,1:T).*deta];%deta=0为正向
    cons=[cons,Cw*W_ijbar<= Cw*W_ijmax(:,1:T).*(1-deta)];%deta=1为反向
    ordinary_water_kappa=Cw*R_w(:,1:T);
    ordinary_water_qmax=Cw*W_ijmax(:,1:T);
    % 双向水管SOC辅助压降模型：q=qP-qM，h=hP-hM。
    % hP/hM的线性上包络与旋转SOC下界共同给出凸包；目标函数
    % 最小化hP+hM以推动SOC松弛达到原始Darcy-Weisbach等式。
    ordinary_water_count=size(Cw,1);
    water_qP=sdpvar(ordinary_water_count,T,'full');
    water_qM=sdpvar(ordinary_water_count,T,'full');
    water_hP=sdpvar(ordinary_water_count,T,'full');
    water_hM=sdpvar(ordinary_water_count,T,'full');
    cons=[cons,Cw*W_ijbar==water_qP-water_qM];
    cons=[cons,0<=water_qP,water_qP<=ordinary_water_qmax.*(1-deta)];
    cons=[cons,0<=water_qM,water_qM<=ordinary_water_qmax.*deta];
    cons=[cons,0<=water_hP,water_hP<=ordinary_water_kappa.*ordinary_water_qmax.*water_qP];
    cons=[cons,0<=water_hM,water_hM<=ordinary_water_kappa.*ordinary_water_qmax.*water_qM];
    cons=[cons,D_ij==water_hP-water_hM];
    for water_pipe_index=1:ordinary_water_count
        for water_time_index=1:T
            kappa_value=ordinary_water_kappa(water_pipe_index,water_time_index);
            cons=[cons,cone([sqrt(2*kappa_value)*water_qP(water_pipe_index,water_time_index); ...
                water_hP(water_pipe_index,water_time_index)-0.5], ...
                water_hP(water_pipe_index,water_time_index)+0.5)];
            cons=[cons,cone([sqrt(2*kappa_value)*water_qM(water_pipe_index,water_time_index); ...
                water_hM(water_pipe_index,water_time_index)-0.5], ...
                water_hM(water_pipe_index,water_time_index)+0.5)];
        end
    end
    %流量上下限约束
    cons=[cons,-W_ijmax(:,1:T)<=W_ij];
    cons=[cons,W_ij<=W_ijmax(:,1:T)];
    cons=[cons,-W_ijmax(:,1:T)<=W_ji];
    cons=[cons,W_ji<=W_ijmax(:,1:T)];
    %压力上下限约束
    % cons=[cons,H_imin<=H_i];
    cons=[cons,0<=H_i];
    cons=[cons,H_i<=H_imax(:,1:T)];
    %........................负荷压力驱动约束...................%
    %水负荷压力驱动
    H_ireq=3*H_imin;
    if strcmp(service_model,'fixed_pressure')&&strcmp(fixed_pressure_formulation,'compact_exact')
        % 优化固定压力需求：恢复比例仍可连续调度，但只要供给为正，
        % 对应可用状态就必须达到额定压力。该紧凑投影删除了三段压力区的冗余二进制。
        cons=[cons,w_i0(:,2:end)<=fixed_pressure_w_available];
        cons=[cons,H_i(2:end,:)>=H_ireq(2:end,1:T).*fixed_pressure_w_available];
        % 完整描述额定压力分区：低于阈值时available=0，高于阈值时available=1。
        % 该反向界不改变(H,w)投影可行域，只消除高压零供给时available可取0/1的对称分支。
        cons=[cons,H_i(2:end,:)<=H_ireq(2:end,1:T)+ ...
            (H_imax(2:end,1:T)-H_ireq(2:end,1:T)).*fixed_pressure_w_available];
    elseif strcmp(service_model,'pressure_driven')&&strcmp(pressure_driven_formulation,'compact_exact')
        % 三段压力驱动曲线的精确单二进制表述。active=0对应零供给低压区；
        % active=1同时覆盖线性压力驱动区和额定供给平台，平台由w<=1直接截断。
        water_pressure_span=H_ireq(2:end,1:T)-H_imin(2:end,1:T);
        cons=[cons,w_i0(:,2:end)<=pressure_driven_w_active];
        cons=[cons,H_i(2:end,:)>=H_imin(2:end,1:T).*pressure_driven_w_active];
        cons=[cons,H_i(2:end,:)<=H_imin(2:end,1:T)+ ...
            (H_imax(2:end,1:T)-H_imin(2:end,1:T)).*pressure_driven_w_active];
        cons=[cons,w_i0(:,2:end)<=(H_i(2:end,:)-H_imin(2:end,1:T))./water_pressure_span+ ...
            (H_imin(2:end,1:T)./water_pressure_span).*(1-pressure_driven_w_active)];
    else
        if water_region_move_blocking&&T>1
            vw1=repmat(binvar(kw,1,'full'),1,T);vw2=repmat(binvar(kw,1,'full'),1,T);vw3=repmat(binvar(kw,1,'full'),1,T);
        else
            vw1=binvar(kw,T,'full');vw2=binvar(kw,T,'full');vw3=binvar(kw,T,'full');%状态分类
        end
        cons=[cons,vw1+vw2+vw3==1];
        cons=[cons,vw2.*H_imin(2:end,1:T)+vw3.*H_ireq(2:end,1:T)<=H_i(2:end,:)];%压力区分类
        cons=[cons,H_i(2:end,:)<=vw1.*H_imin(2:end,1:T)+vw2.*H_ireq(2:end,1:T)+vw3.*H_imax(2:end,1:T)];%压力区分类
        if strcmp(service_model,'fixed_pressure')
            cons=[cons,0<=w_i0(:,2:end),w_i0(:,2:end)<=vw3];%原始扩展表述
        else
            cons=[cons,(w_i0(:,2:end)>=0)&(w_i0(:,2:end)<=1-vw1)];%零负荷区
            cons=[cons,(w_i0(:,2:end)>=0)&(w_i0(:,2:end)<=(H_i(2:end,:)-H_imin(2:end,1:T))./(H_ireq(2:end,1:T)-H_imin(2:end,1:T))+1-vw2)];%压力驱动区
        end
    end
    %气负荷压力驱动
    kesai_ireq=3*kesai_imin;
    kesai_isaf=1.2*kesai_imin;
    g_imin=0.1*ones(kg,T);%最小安全流量
    if strcmp(service_model,'fixed_pressure')&&strcmp(fixed_pressure_formulation,'compact_exact')
        % 与水网相同：恢复比例连续，二进制变量只表示额定压力可用状态。
        cons=[cons,g_i0(:,2:end)<=fixed_pressure_g_available];
        cons=[cons,A2*kesai_i>=A2*kesai_ireq(:,1:T).*fixed_pressure_g_available];
        if gas_customer_inactive_pressure_upper
            cons=[cons,A2*kesai_i<=A2*kesai_ireq(:,1:T)+ ...
                (A2*kesai_imax(:,1:T)-A2*kesai_ireq(:,1:T)).*fixed_pressure_g_available];
        end
    elseif strcmp(service_model,'pressure_driven')&&strcmp(pressure_driven_formulation,'compact_exact')
        gas_pressure_span=A2*kesai_ireq(:,1:T)-A2*kesai_imin(:,1:T);
        cons=[cons,g_i0(:,2:end)<=pressure_driven_g_active];
        cons=[cons,g_i0(:,2:end)>=g_imin.*pressure_driven_g_active];
        cons=[cons,A2*kesai_i>=A2*kesai_isaf(:,1:T).*pressure_driven_g_active];
        if gas_customer_inactive_pressure_upper
            % 仅用于旧结果复现。物理恢复模型允许节点已经升压而用户阀门仍关闭，
            % 因而默认不再用“停供”状态反向限制节点压力。
            cons=[cons,A2*kesai_i<=A2*kesai_isaf(:,1:T)+ ...
                (A2*kesai_imax(:,1:T)-A2*kesai_isaf(:,1:T)).*pressure_driven_g_active];
        end
        cons=[cons,g_i0(:,2:end)<=(A2*kesai_i-A2*kesai_imin(:,1:T))./gas_pressure_span+ ...
            ((A2*kesai_imin(:,1:T))./gas_pressure_span).*(1-pressure_driven_g_active)];
    else
        if binary_move_blocking&&T>1
            vg1=repmat(binvar(kg,1,'full'),1,T);vg2=repmat(binvar(kg,1,'full'),1,T);vg3=repmat(binvar(kg,1,'full'),1,T);
        else
            vg1=binvar(kg,T,'full');vg2=binvar(kg,T,'full');vg3=binvar(kg,T,'full');%状态分类
        end
        cons=[cons,vg1+vg2+vg3==1];
        cons=[cons,vg2.*(A2*kesai_isaf(:,1:T))+vg3.*(A2*kesai_ireq(:,1:T))<=A2*kesai_i];%压力区分类
        cons=[cons,A2*kesai_i<=vg1.*(A2*kesai_isaf(:,1:T))+vg2.*(A2*kesai_ireq(:,1:T))+vg3.*(A2*kesai_imax(:,1:T))];%压力区分类
        if strcmp(service_model,'fixed_pressure')
            cons=[cons,0<=g_i0(:,2:end),g_i0(:,2:end)<=vg3];%原始扩展表述
        else
            cons=[cons,(g_i0(:,2:end)>=0)&(g_i0(:,2:end)<=1-vg1)];%零负荷区
            cons=[cons,(g_i0(:,2:end)>=g_imin.*(1-vg1))&(g_i0(:,2:end)<=(A2*kesai_i-A2*kesai_imin(:,1:T))./(A2*kesai_ireq(:,1:T)-A2*kesai_imin(:,1:T))+1-vg2)];%安全压力驱动区
        end
    end
    %燃气轮机压力驱动
    kssk=repmat(GT(:,4)./effi/fai,1,T);%逆转化率
    % 与普通气负荷共用同一压力标尺：非气源节点统一采用10--30 bar；
    % 气源节点没有恢复负荷，其压力由物理边界及定压/定流源条件决定。
    kesai_ireqGT=3*kesai_imin(GT(:,3),1:T);
    kesai_isafGT=1.2*kesai_imin(GT(:,3),1:T);
    g_GTmin=0.1*ones(3,T);%最小安全流量
    if strcmp(service_model,'fixed_pressure')&&strcmp(fixed_pressure_formulation,'compact_exact')
        % GT已有启停变量u_it；直接复用为压力可用状态，避免为同一开机事件
        % 再引入一套独立二进制并产生对称/冲突组合。
        fixed_pressure_gt_available=u_it;
        cons=[cons,vec(g_GT./kssk)>=vec(g_GTmin.*fixed_pressure_gt_available), ...
            vec(g_GT./kssk)<=vec(fixed_pressure_gt_available)];
        cons=[cons,vec(kesai_i(GT(:,3),:))>=vec(kesai_ireqGT.*fixed_pressure_gt_available)];
        if gt_inactive_pressure_upper
            % 与Proposed使用相同的“压力满足则GT激活”开关。fixed-pressure
            % 的激活门槛为额定压力Psi_req，而不是压力驱动模型的安全压力Psi_safe。
            cons=[cons,vec(kesai_i(GT(:,3),:))<=vec(kesai_ireqGT+ ...
                (kesai_imax(GT(:,3),1:T)-kesai_ireqGT).*fixed_pressure_gt_available)];
        end
    elseif strcmp(service_model,'pressure_driven')&&strcmp(pressure_driven_formulation,'compact_exact')
        gt_pressure_span=kesai_ireqGT-kesai_imin(GT(:,3),1:T);
        gt_service=g_GT./kssk;
        cons=[cons,vec(gt_service)<=vec(pressure_driven_gt_active)];
        cons=[cons,vec(gt_service)>=vec(g_GTmin.*pressure_driven_gt_active)];
        cons=[cons,vec(kesai_i(GT(:,3),:))>=vec(kesai_isafGT.*pressure_driven_gt_active)];
        if gt_inactive_pressure_upper
            % 旧式反向界：停机时强制GT节点不超过安全压力。该界会把
            % “管网已经升压”反向解释为“GT必须开机”，仅用于旧结果复现。
            cons=[cons,vec(kesai_i(GT(:,3),:))<=vec(kesai_isafGT+ ...
                (kesai_imax(GT(:,3),1:T)-kesai_isafGT).*pressure_driven_gt_active)];
        end
        cons=[cons,vec(gt_service)<=vec((kesai_i(GT(:,3),:)-kesai_imin(GT(:,3),1:T))./gt_pressure_span+ ...
            (kesai_imin(GT(:,3),1:T)./gt_pressure_span).*(1-pressure_driven_gt_active))];
    else
        if binary_move_blocking&&T>1
            vg1GT=repmat(binvar(n2,1,'full'),1,T);vg2GT=repmat(binvar(n2,1,'full'),1,T);vg3GT=repmat(binvar(n2,1,'full'),1,T);
        else
            vg1GT=binvar(n2,T,'full');vg2GT=binvar(n2,T,'full');vg3GT=binvar(n2,T,'full');%状态分类
        end
        cons=[cons,vg1GT+vg2GT+vg3GT==1];
        cons=[cons,vec(vg2GT.*kesai_isafGT+vg3GT.*kesai_ireqGT)<=vec(kesai_i(GT(:,3),:))];%压力区分类
        cons=[cons,vec(kesai_i(GT(:,3),:))<=vec(vg1GT.*kesai_isafGT+vg2GT.*kesai_ireqGT+vg3GT.*kesai_imax(GT(:,3),1:T))];%压力区分类
        if strcmp(service_model,'fixed_pressure')
            cons=[cons,vec(g_GT./kssk)>=0,vec(g_GT./kssk)<=vec(vg3GT)];
        else
            cons=[cons,vec(g_GT./kssk)>=0,vec(g_GT./kssk)<=vec(1-vg1GT)];%零负荷区
            cons=[cons,vec(g_GT./kssk)>=vec(g_GTmin.*(1-vg1GT)),vec(g_GT./kssk)<=vec((kesai_i(GT(:,3),:)-kesai_imin(GT(:,3),1:T))./(kesai_ireqGT-kesai_imin(GT(:,3),1:T))+1-vg2GT)];%安全压力驱动区
        end
    end
    %........................凸松弛.............................%
    %水泵约束
    cons=[cons,Bw*W_ji==-Bw*W_ij];%水泵不储水，两个端点流量守恒
    % cons=[cons,-M*(1-Xm)+(alpha.*(Bw*W_ij).*omima+beta.*omima.^2)<=Gp_ij];%凸约束非凸约束
    % cons=[cons,M*(1-Xm)-alpha.^2./(4*beta).*(Bw*W_ij).*(Bw*W_ij_temp)+beta.*(omima+0.5*alpha./beta.*(Bw*W_ij)).*...
    %     (omima_temp+0.5*alpha./beta.*(Bw*W_ij_temp))+sasa>=Gp_ij];%非凸约束
    cons=[cons,vec(Bw*W_ij)>=0];
    cons=[cons,vec(Bw*W_ij)<=vec((Bw*W_ijmax(:,1:T)).*Xm)];
    if strcmp(pump_speed_model,'projection')
        omega_lower=zeros(size(WDS_pump,1),T);
        omega_upper=ones(size(WDS_pump,1),T);
        [pump_head_lower,pump_head_upper]=pump_head_projection_bounds( ...
            alpha(:,1:T),beta(:,1:T),Bw*W_ij,omega_lower,omega_upper);
        cons=[cons,vec(Gp_ij)>=vec(pump_head_lower-M*(1-Xm))];
        cons=[cons,vec(Gp_ij)<=vec(pump_head_upper+M*(1-Xm))];
        cons=[cons,vec(P_pump)<=vec(Pp_max(:,1:T).*Xm)];%停泵功率严格为0
    else
        cons=[cons,vec(omima)>=0,vec(omima)<=vec(Xm)];%停泵时转速必须为0，避免无物理意义的自由转速振荡
    end
    % cons=[cons,P_pump==g*(1./eff).*Gp_ij.*W_ij_temp];%水泵功率
    %% 风电机会约束的确定性等价形式
    set0=importdata('盒式不确定集.txt');
    % 延长恢复诊断超过原预测长度时保持末端预测值，不额外增加风电。
    wind_index=min(t:t+T-1,size(set0.data,2));
    base_set=set0.data(1:3,wind_index);%日前风电预测结果
    max_set=set0.data(4,wind_index);%不确定集上界
    min_set=set0.data(5,wind_index);%不确定集下界
    epsilon_wind=0.05;%单侧机会约束允许的违约概率
    box_coverage=0.95;%将盒式上下界解释为中心预测区间的覆盖率
    [wind_lower,wind_mean,wind_sigma,wind_z]=wind_chance_limit_from_box(...
        base_set,min_set,max_set,epsilon_wind,box_coverage);
    ops = sdpsettings('verbose',0,'solver','gurobi');
    ops.gurobi.Seed=0;
    requested_threads=str2double(getenv('CEWGS_THREADS'));
    if ~isfinite(requested_threads), requested_threads=1; end
    ops.gurobi.Threads=requested_threads;
    ops.gurobi.OutputFlag=0;
    ops.savesolveroutput=capture_solver_trace;
    ops.savesolverinput=capture_solver_trace;
    requested_mip_gap=str2double(getenv('CEWGS_MIP_GAP'));
    if isfinite(requested_mip_gap)
        assert(requested_mip_gap>=0&&requested_mip_gap<1,'CEWGS_MIP_GAP must be in [0,1).');
        ops.gurobi.MIPGap=requested_mip_gap;
    end
    requested_time_limit=str2double(getenv('CEWGS_SOLVE_TIME_LIMIT'));
    if isfinite(requested_time_limit)
        assert(requested_time_limit>0,'CEWGS_SOLVE_TIME_LIMIT must be positive.');
        ops.gurobi.TimeLimit=requested_time_limit;
    end
    rolling_window_time_limit=str2double(getenv('CEWGS_ROLLING_WINDOW_TIME_LIMIT'));
    if isfinite(rolling_window_time_limit)
        assert(rolling_window_time_limit>0, ...
            'CEWGS_ROLLING_WINDOW_TIME_LIMIT must be positive.');
    end
    % ops.gurobi.MIPGap = 0.0005;
    % ops.gurobi.TimeLimit = 40;
    cons=[cons,vec(p_W)>=0,vec(p_W)<=vec(wind_lower.*(F'*y_it(:,2:T+1)))];
    %%................内层问题初始化.................%%
    %% 附加变量
    if strcmp(pump_speed_model,'sca')
        sasa1=sdpvar(size(WDS_pump,1),T,'full');%原转速SCA的水泵非凸松弛变量
        sasa2=sdpvar(size(WDS_pump,1),T,'full');
    else
        sasa1=[];sasa2=[];%转速投影后不再需要泵曲线SCA松弛变量
    end
    if strcmp(pump_power_sca_mode,'full')
        power_sasa=sdpvar(size(WDS_pump,1),T,'full');%参考Zamzam等FPP-SCA：功率等式两侧共用一个松弛，单位kW
    else
        power_sasa=[];
    end
    if enable_pipe_reverse_sca
        % 隔离测试：给普通水管摩阻等式及普通气管Weymouth等式补上
        % 非凸反向侧的FPP-SCA内逼近。默认关闭，不改变四方案正式基线。
        water_pipe_sasaP=sdpvar(ordinary_water_count,T,'full');
        water_pipe_sasaM=sdpvar(ordinary_water_count,T,'full');
        gas_pipe_sasa=sdpvar(l1-cp1,T,'full');
    else
        water_pipe_sasaP=[];water_pipe_sasaM=[];gas_pipe_sasa=[];
    end
    sca_penalty_start=str2double(getenv('CEWGS_SCA_PENALTY_START'));
    if ~(isfinite(sca_penalty_start)&&sca_penalty_start>0)
        sca_penalty_start=1;%单因素测试确认可减少达到物理精度所需的迭代
    end
    gigo3=sca_penalty_start;gigo4=sca_penalty_start;%仅保留水泵FPP-SCA惩罚因子
    water_pipe_sca_penalty=pipe_sca_penalty_start;
    gas_pipe_sca_penalty=pipe_sca_penalty_start;
    pump_penalty_multiplier=2;%水泵FPP-SCA惩罚快速更新
    pipe_penalty_multiplier=1.5;%水管/气管反向SOC-SCA惩罚平缓更新
    gigo3max=1000;gigo4max=1000;%惩罚因子最大值
    pipe_sca_penalty_max=1000;
    sca_tol_ref=str2double(getenv('CEWGS_SCA_FLOW_REF_TOL'));
    if ~(isfinite(sca_tol_ref)&&sca_tol_ref>0),sca_tol_ref=1e-2;end%水泵流量参考点阈值，单位m^3/s
    sca_iter_max=str2double(getenv('CEWGS_SCA_ITER_MAX'));
    if ~(isfinite(sca_iter_max)&&sca_iter_max>=1&&sca_iter_max==floor(sca_iter_max))
        sca_iter_max=40;
    end
    sca_pump_ref_tol=str2double(getenv('CEWGS_SCA_PUMP_REF_TOL'));
    if ~(isfinite(sca_pump_ref_tol)&&sca_pump_ref_tol>0)
        sca_pump_ref_tol=1;%泵扬程参考点仅用于迭代稳定性判断，默认容差1 m
    end
    sca_pump_tol=str2double(getenv('CEWGS_SCA_PUMP_TOL'));
    if ~(isfinite(sca_pump_tol)&&sca_pump_tol>=0),sca_pump_tol=0.1;end%原始水泵扬程等式验收容差，单位m
    pump_power_tol=str2double(getenv('CEWGS_PUMP_POWER_TOL_KW'));
    if ~(isfinite(pump_power_tol)&&pump_power_tol>=0),pump_power_tol=2;end
    W_ij_temp0=zeros(size(W_ijmax(:,1:T)));%仅水泵行用作泵功率FPP-SCA参考点
    Gp_ij_temp0=max(0,repmat(Bw*Aw'*H_i0,1,T));%功率完整SCA的上一执行水头参考点
    if strcmp(pump_speed_model,'sca')
        omima_temp0=ones(size(WDS_pump,1),T);
    else
        omima_temp0=[];
        W_ij_temp0(WDS_pump(:,1),:)=repmat(max(0,W_ij0(WDS_pump(:,1))),1,T);%转速消元后用上一执行泵流量作参考
    end
    if sca_rolling_warm_start&&~isempty(sca_flow_seed)
        pump_seed=repmat(sca_flow_seed(WDS_pump(:,1)),1,T);
        W_ij_temp0(WDS_pump(:,1),:)=sca_warm_start_weight*pump_seed+ ...
            (1-sca_warm_start_weight)*W_ij_temp0(WDS_pump(:,1),:);%仅泵流量采用阻尼热启动
        if strcmp(pump_speed_model,'sca')
            omima_temp0=sca_warm_start_weight*repmat(sca_speed_seed,1,T)+(1-sca_warm_start_weight)*omima_temp0;%泵速阻尼热启动
        end
    end
    W_ij_temp=W_ij_temp0;
    Gp_ij_temp=Gp_ij_temp0;
    omima_temp=omima_temp0;
    water_q_initial=Cw*(W_ij0-W_ji0)*0.5;
    water_qP_ref=repmat(max(water_q_initial,0),1,T);
    water_qM_ref=repmat(max(-water_q_initial,0),1,T);
    gas_q_ref=repmat(gas_flow_seed,1,T);
    gas_pL_ref=repmat(gas_low_pressure_seed,1,T);
    % Qs=[0,alpha(1)/2;
    %     alpha(1)/2,beta(1)];
    % [Vs,Ds]=eig(Qs);%计算A的特征值对角阵D和特征向量V，使AV=VD成立
    % lambda_P=Ds;
    % lambda_N=-Ds;
    % lambda_P(find(Ds<0))=0;
    % lambda_N(find(Ds>0))=0;
    % Ps=Vs*lambda_P*Vs';
    % Ns=Vs*lambda_N*Vs';
    sca_iter=0;
    plan_lock_cons=[];
    plan_locked=false;
    rpp_topology_cons=[];
    rpp_topology_plan=[];
    while(1)
        sca_iter=sca_iter+1;
        cons_in=[];
        %...........................................%
        if strcmp(pump_speed_model,'sca')
            cons_in=[cons_in,vec(M*(1-Xm)-alpha(:,1:T).^2./(4*beta(:,1:T)).*(Bw*W_ij).^2+beta(:,1:T).*(omima+0.5*alpha(:,1:T)./beta(:,1:T).*(Bw*W_ij)).*...
                (omima_temp+0.5*alpha(:,1:T)./beta(:,1:T).*(Bw*W_ij_temp))+sasa1)>=vec(Gp_ij)];%原泵曲线SCA非凸侧
            cons_in=[cons_in,vec(-M*(1-Xm)-alpha(:,1:T).^2./(4*beta(:,1:T)).*(Bw*W_ij).*(Bw*W_ij_temp)+beta(:,1:T).*(omima+0.5*alpha(:,1:T)./beta(:,1:T).*(Bw*W_ij)).^2-sasa2)<=vec(Gp_ij)];
        end
        %...............................................%
        %     cons_in=[cons_in,M*(1-Xm)-alpha.^2./(4*beta).*(Bw*W_ij).*(Bw*W_ij_temp)+beta.*(omima+0.5*alpha./beta.*(Bw*W_ij)).*...
        %         (omima_temp+0.5*alpha./beta.*(Bw*W_ij_temp))+sasa1>=Gp_ij];%非凸约束
        %     cons_in=[cons_in,-M*(1-Xm)+(alpha.*(Bw*W_ij).*omima_temp+beta.*omima.^2)-sasa2<=Gp_ij];%凸约束非凸约束
        if strcmp(pump_power_sca_mode,'full')
            pump_q_scale=max(Bw*W_ijmax(:,1:T),1e-3);
            pump_H_scale=max(H_imax(:,1:T),[],'all')*ones(size(WDS_pump,1),T);
            [power_ge_residual,power_le_residual]=pump_power_sca_dc_residuals( ...
                g,eff(:,1:T),P_pump,Gp_ij,Bw*W_ij,Gp_ij_temp,Bw*W_ij_temp, ...
                pump_q_scale,pump_H_scale);
            pump_power_scale=g.*(1./eff(:,1:T)).*pump_q_scale.*pump_H_scale;%归一化功率基值，kW
            power_projection_big_m=10;%归一化后的泵功率DC残差大M
            cons_in=[cons_in,vec(power_ge_residual)<=vec(power_sasa./pump_power_scale+power_projection_big_m*(1-Xm))];%P>=gHq/eta：保留凸二次项，线性化右侧凸项
            cons_in=[cons_in,vec(power_le_residual)<=vec(power_sasa./pump_power_scale+power_projection_big_m*(1-Xm))];%P<=gHq/eta：保留凸二次项，线性化右侧凸项
            cons_in=[cons_in,vec(power_sasa)>=0];
        else
            cons_in=[cons_in,P_pump==g*(1./eff(:,1:T)).*Gp_ij.*(Bw*W_ij_temp)];%V12原单边功率近似
        end
        if strcmp(pump_speed_model,'sca')
            cons_in=[cons_in,vec(sasa1)>=0];
            cons_in=[cons_in,vec(sasa2)>=0];
        end
        if enable_pipe_reverse_sca
            % h <= k q^2 的一阶下切平面加非负松弛；与已有
            % h >= k q^2 的旋转SOC共同逼近水管摩阻等式。
            water_reverse_rhsP=2*ordinary_water_kappa.*water_qP_ref.*water_qP- ...
                ordinary_water_kappa.*water_qP_ref.^2;
            water_reverse_rhsM=2*ordinary_water_kappa.*water_qM_ref.*water_qM- ...
                ordinary_water_kappa.*water_qM_ref.^2;
            cons_in=[cons_in,vec(water_hP)<=vec(water_reverse_rhsP+water_pipe_sasaP)];
            cons_in=[cons_in,vec(water_hM)<=vec(water_reverse_rhsM+water_pipe_sasaM)];
            cons_in=[cons_in,vec(water_pipe_sasaP)>=0,vec(water_pipe_sasaM)>=0];

            % pH^2 <= pL^2+(q/C)^2 的凸右侧在上一参考点处作
            % 一阶下估计；与已有SOC相反方向共同逼近Weymouth等式。
            gas_scaled_flow=Bg2*(G_ijave./Cij(:,1:T));
            gas_reverse_rhs=2*gas_pL_ref.*kesai_ijL-gas_pL_ref.^2+ ...
                2*gas_q_ref.*gas_scaled_flow-gas_q_ref.^2;
            cons_in=[cons_in,vec(kesai_ijH.^2)<=vec(gas_reverse_rhs+gas_pipe_sasa)];
            cons_in=[cons_in,vec(gas_pipe_sasa)>=0];
        end
        %%.................目标函数.............................%%
        switch_action_expr=sum(repmat(switch_action_weight,1,T).*(a_ij(:,2:T+1)-a_ij(:,1:T)),'all');%联络线合闸含额外操作权重
        Fs=sum(omiga*(1-r_i),'all')+sum(omiga_g*(1-g_i0(:,2:end)),'all')+sum(omiga*(1-w_i0(:,2:end)),'all')+...
            switch_action_cost*switch_action_expr;
        water_soc_pressure_drop=sum(water_hP+water_hM,'all');
        Fs=Fs+water_soc_penalty_weight*water_soc_pressure_drop;
        % Zhang et al. (IEEE TSG 2024), Eq. (66): minimizing the pressure
        % drop of ordinary gas pipes promotes equality in the relaxed SOC.
        gas_soc_pressure_drop=sum(kesai_ijH-kesai_ijL,'all');
        Fs=Fs+gas_soc_penalty_weight*gas_soc_pressure_drop;
        if strcmp(pump_speed_model,'sca')
            Fs=Fs+gigo3*sum(sasa1(:))+gigo4*sum(sasa2(:));
        end
        if strcmp(pump_power_sca_mode,'full')
            Fs=Fs+max(gigo3,gigo4)*sum(power_sasa(:));
        end
        if enable_pipe_reverse_sca
            Fs=Fs+water_pipe_sca_penalty*sum(water_pipe_sasaP(:)+water_pipe_sasaM(:))+ ...
                gas_pipe_sca_penalty*sum(gas_pipe_sasa(:));
        end
        sca_solve_tic=tic;
        if isfinite(rolling_window_time_limit)
            rolling_window_remaining=rolling_window_time_limit-toc(step_tic);
            if rolling_window_remaining<=0
                error('ChanceRestoration:RollingWindowTimeLimit', ...
                    'Rolling step %d exceeded its cumulative %.1f s time limit before SCA iteration %d.', ...
                    t,rolling_window_time_limit,sca_iter);
            end
            if isfinite(requested_time_limit)
                ops.gurobi.TimeLimit=max(0.1,min(requested_time_limit,rolling_window_remaining));
            else
                ops.gurobi.TimeLimit=max(0.1,rolling_window_remaining);
            end
        end
        ops.usex0=sca_mip_warm_start&&sca_iter>1;
        % 本机YALMIP在sdpsettings构造后不会由usex0字段回填warmstart；
        % 两个字段同时设置，才能在yalmip2gurobi中生成model.start。
        ops.warmstart=ops.usex0;
        if capture_solver_trace
            ops.gurobi.OutputFlag=1;
            ops.gurobi.LogFile=fullfile(output_dir,sprintf('roll_%02d_sca_%02d_gurobi.log',t,sca_iter));
        end
        if strcmp(topology_solution_mode,'relax_project_polish')&&isempty(rpp_topology_plan)
            rpp_context=struct('line',main_line,'roots',GT(:,2), ...
                'y_previous',y_it0,'a_previous',a_ij0,'switch_budget',Kt, ...
                'commodity_bound',Mf);
            rpp_context.node_benefit=zeros(size(load,1),1);
            rpp_context.node_benefit(find(sum(A1,2)>0))=omiga(:);
            rpp_context.protected=struct('variable',{Xm,Xc}, ...
                'buses',{WDS_pump(:,6),GDS_compessor(:,9)});
            rpp_options=struct('relax_time_limit',rpp_relax_time_limit, ...
                'projection_time_limit',rpp_projection_time_limit, ...
                'polish_time_limit',rpp_polish_time_limit,'mip_gap',requested_mip_gap, ...
                'threads',requested_threads,'seed',0,'output_flag',0);
            [sol,rpp_trace,rpp_topology_plan]=solve_topology_relax_project_polish( ...
                [cons,plan_lock_cons,cons_in],Fs,ops,y_it(:,2:T+1), ...
                a_ij(:,2:T+1),rpp_context,rpp_options);
            if sol.problem==0
                rpp_topology_cons=[y_it(:,2:T+1)==rpp_topology_plan.y, ...
                    a_ij(:,2:T+1)==rpp_topology_plan.a];
                rpp_trace_r{t}=rpp_trace;
            end
        elseif strcmp(topology_solution_mode,'mcb_prototype')&&isempty(rpp_topology_plan)
            mcb_context=struct('line',main_line,'roots',GT(:,2), ...
                'y_previous',y_it0,'a_previous',a_ij0,'switch_budget',Kt, ...
                'commodity_bound',Mf);
            graph_context=struct('line',main_line,'roots',GT(:,2), ...
                'y_previous',y_it0,'a_previous',a_ij0, ...
                'customer_nodes',find(sum(A1,2)>0),'priority',omiga(:), ...
                'pump_buses',WDS_pump(:,6), ...
                'water_pump_paths',restoration_water_pump_paths( ...
                    WDS_pipe,WDS_pump,(2:size(WDS_node,1))'), ...
                'compressor_buses',GDS_compessor(:,9), ...
                'gas_demand',GDS_node(GDS_node(:,6)~=0,6), ...
                'gas_base_capacity',1.5, ...
                'compressor_capacity_gain',[3;zeros(max(cp1-1,0),1)], ...
                'switch_budget',Kt);
            mcb_context.graph=graph_context;
            mcb_options=struct('max_iterations',mcb_max_iterations, ...
                'relative_gap',mcb_relative_gap,'global_lower_bound',0, ...
                'master_time_limit',10,'relax_time_limit',mcb_relax_time_limit, ...
                'integer_time_limit',mcb_integer_time_limit, ...
                'mip_gap',requested_mip_gap,'threads',requested_threads,'seed',0);
            [sol,mcb_trace,rpp_topology_plan]=solve_topology_mcb_prototype( ...
                [cons,plan_lock_cons,cons_in],Fs,ops,y_it(:,2:T+1), ...
                a_ij(:,2:T+1),mcb_context,mcb_options);
            if sol.problem==0
                rpp_topology_cons=[y_it(:,2:T+1)==rpp_topology_plan.y, ...
                    a_ij(:,2:T+1)==rpp_topology_plan.a];
                mcb_trace_r{t}=mcb_trace;
            end
        else
            sol=solvesdp([cons,plan_lock_cons,rpp_topology_cons,cons_in],Fs,ops);%恢复方案稳定后仅抛光连续SCA物理
        end
        sca_solve_time_trace_r{t}(sca_iter)=toc(sca_solve_tic);
        if capture_solver_trace
            solver_output_trace_r{t}{sca_iter}=sol.solveroutput;
        end
        if sol.problem~=0
            if capture_solver_trace
                failure_file=fullfile(output_dir,sprintf('failed_roll_%02d_sca_%02d.mat',t,sca_iter));
                failure_solver_output=sol.solveroutput;
                if isfield(sol,'solverinput'),failure_solver_input=sol.solverinput;else,failure_solver_input=[];end
                save(failure_file,'failure_solver_output','failure_solver_input','t','sca_iter', ...
                    'service_model','storage_model','water_tank_model','gas_tank_model', ...
                    'water_network_model','gas_network_model', ...
                    'fixed_pressure_formulation','topology_formulation');
            end
            error('ChanceRestoration:SolverFailure', ...
                'Rolling step %d, SCA iteration %d failed with YALMIP code %d: %s', ...
                t,sca_iter,sol.problem,sol.info);
        end
        if isfinite(rolling_window_time_limit)&&toc(step_tic)>rolling_window_time_limit
            error('ChanceRestoration:RollingWindowTimeLimit', ...
                'Rolling step %d exceeded its cumulative %.1f s time limit during SCA iteration %d.', ...
                t,rolling_window_time_limit,sca_iter);
        end
        display([sol.info,' 目标函数值：',num2str(value(Fs))])
        W_ij_temp_before=W_ij_temp;
        Gp_ij_temp_before=Gp_ij_temp;
        if strcmp(pump_speed_model,'sca'),omima_temp_before=omima_temp;end
        W_ij_temp(WDS_pump(:,1),:)=value(W_ij(WDS_pump(:,1),:));
        Gp_ij_temp=value(Gp_ij);
        if strcmp(pump_speed_model,'sca'),omima_temp=value(omima);end
        water_soc_gap_values=[value(water_hP-ordinary_water_kappa.*water_qP.^2); ...
            value(water_hM-ordinary_water_kappa.*water_qM.^2)];
        water_soc_gap=max(max(water_soc_gap_values,0),[],'all');
        gas_scaled_flow_values=value(Bg2*(G_ijave./Cij(:,1:T)));
        gas_soc_gap_values=value(kesai_ijH).^2-value(kesai_ijL).^2-gas_scaled_flow_values.^2;
        gas_soc_gap=max(max(gas_soc_gap_values,0),[],'all');
        if enable_pipe_reverse_sca
            water_pipe_slack=max([value(water_pipe_sasaP(:));value(water_pipe_sasaM(:))]);
            gas_pipe_slack=max(value(gas_pipe_sasa(:)));
            water_qP_ref=(1-pipe_sca_reference_damping)*water_qP_ref+ ...
                pipe_sca_reference_damping*value(water_qP);
            water_qM_ref=(1-pipe_sca_reference_damping)*water_qM_ref+ ...
                pipe_sca_reference_damping*value(water_qM);
            gas_q_ref=(1-pipe_sca_reference_damping)*gas_q_ref+ ...
                pipe_sca_reference_damping*gas_scaled_flow_values;
            gas_pL_ref=(1-pipe_sca_reference_damping)*gas_pL_ref+ ...
                pipe_sca_reference_damping*value(kesai_ijL);
        else
            water_pipe_slack=0;gas_pipe_slack=0;
        end
        pump_slack_values=[];
        if strcmp(pump_speed_model,'sca')
            pump_slack_values=[pump_slack_values;value(sasa1(:));value(sasa2(:))];
        end
        if strcmp(pump_power_sca_mode,'full')
            pump_slack_values=[pump_slack_values;value(power_sasa(:))];
        end
        if isempty(pump_slack_values),pump_slack=0;else,pump_slack=max(pump_slack_values);end
        sca_slack=pump_slack;
        pump_flow_reference_change=W_ij_temp(WDS_pump(:,1),:)-W_ij_temp_before(WDS_pump(:,1),:);
        ref_flow_rms=sqrt(sum(pump_flow_reference_change(:).^2)/numel(pump_flow_reference_change));
        if strcmp(pump_speed_model,'sca')
            ref_speed_rms=sqrt(sum((omima_temp(:)-omima_temp_before(:)).^2)/numel(omima_temp));
        else
            ref_speed_rms=sqrt(sum((Gp_ij_temp(:)-Gp_ij_temp_before(:)).^2)/numel(Gp_ij_temp));%投影模型用扬程参考点替代转速参考点
        end
        reference_converged=(ref_flow_rms<sca_tol_ref)&&(ref_speed_rms<sca_pump_ref_tol);
        qbar_sca=value(Cw*W_ijbar);
        D_sca=value(D_ij);
        R_sca=Cw*R_w(:,1:T);
        pipe_equation_residual=max(abs(D_sca-R_sca.*qbar_sca.*abs(qbar_sca)),[],'all');
        % 水管SOC精度由惩罚权重控制并作为诊断量输出，不参与水泵SCA停止。
        % 它不依赖参考点更新，重复求解不会自动收紧该凸松弛。
        pump_flow_sca=value(Bw*W_ij);
        pump_head_sca=value(Gp_ij);
        if strcmp(pump_speed_model,'sca')
            pump_speed_sca=value(omima);
            pump_on=value(Xm)>0.5;
            curve_error=abs(pump_head_sca-alpha(:,1:T).*pump_flow_sca.*pump_speed_sca-beta(:,1:T).*pump_speed_sca.^2).*pump_on;
            speed_bound_error=max(max(-pump_speed_sca,0),max(pump_speed_sca-1,0)).*pump_on;
            pump_curve_equation_residual=max([curve_error(:);speed_bound_error(:)]);
        else
            % The monotone pump curve has already been projected exactly onto
            % the (flow, head) plane. Speed is neither optimized nor recovered.
            pump_speed_sca=nan(size(pump_head_sca));
            pump_curve_equation_residual=0;
        end
        pump_power_sca=value(P_pump);
        pump_power_equation_residual=max(abs(pump_power_sca-g*(1./eff(:,1:T)).*pump_head_sca.*pump_flow_sca),[],'all');
        if strcmp(pump_power_sca_mode,'full')&&strcmp(pump_speed_model,'projection')
            pump_slack_tolerance=pump_power_tol;%此时泵松弛量仅包含kW功率松弛
        else
            pump_slack_tolerance=sca_pump_tol;
        end
        if strcmp(pump_speed_model,'projection')
            % With speed projected out, s_p (= power_sasa, kW) is the FPP-SCA
            % feasibility gap.  Use this gap alone for convergence; the
            % reconstructed pump-power residual remains a reported diagnostic.
            pump_physics_converged=(pump_slack<=pump_slack_tolerance);
        else
            pump_physics_converged=(pump_slack<=pump_slack_tolerance)&& ...
                (pump_curve_equation_residual<=sca_pump_tol)&&(pump_power_equation_residual<=pump_power_tol);
        end
        water_pipe_physics_converged=(~enable_pipe_reverse_sca)||(water_soc_gap<=pipe_sca_tol);
        gas_pipe_physics_converged=(~enable_pipe_reverse_sca)||(gas_soc_gap<=pipe_sca_tol);
        sca_physics_accepted=((~enforce_pump_physics)||pump_physics_converged)&& ...
            water_pipe_physics_converged&&gas_pipe_physics_converged;
        if enable_pipe_reverse_sca
            fprintf(['SCA %02d: water gap/slack %.6g/%.6g m, gas gap/slack ' ...
                '%.6g/%.6g bar^2, pump gap %.6g kW.\n'],sca_iter, ...
                water_soc_gap,water_pipe_slack,gas_soc_gap,gas_pipe_slack,pump_slack);
        end
        sca_objective_trace_r{t}(sca_iter)=value(Fs);
        sca_ref_flow_rms_trace_r{t}(sca_iter)=ref_flow_rms;
        sca_ref_pump_rms_trace_r{t}(sca_iter)=ref_speed_rms;
        sca_water_soc_gap_trace_r{t}(sca_iter)=water_soc_gap;
        sca_gas_soc_gap_trace_r{t}(sca_iter)=gas_soc_gap;
        sca_pipe_residual_trace_r{t}(sca_iter)=pipe_equation_residual;
        sca_pump_slack_trace_r{t}(sca_iter)=pump_slack;
        sca_water_pipe_slack_trace_r{t}(sca_iter)=water_pipe_slack;
        sca_gas_pipe_slack_trace_r{t}(sca_iter)=gas_pipe_slack;
        sca_pump_curve_residual_trace_r{t}(sca_iter)=pump_curve_equation_residual;
        sca_pump_power_residual_trace_r{t}(sca_iter)=pump_power_equation_residual;
        if any(strcmp(sca_integer_stabilization, ...
                {'after_first','network_after_first','electric_plan_after_first'}))&& ...
                sca_iter==1&&~sca_physics_accepted
            % 首个MISOCP按给定MIPGap确定恢复计划；后续SCA仅修正连续水力物理。
            % 这不会把连续水/气恢复比例四舍五入，也不会改变首轮整数优化问题。
            plan_lock_cons=[plan_lock_cons, ...
                r_i==round(value(r_i)), ...
                y_it(:,2:end)==round(value(y_it(:,2:end))), ...
                a_ij(:,2:end)==round(value(a_ij(:,2:end)))];
            if strcmp(sca_integer_stabilization,'after_first')
                plan_lock_cons=[plan_lock_cons, ...
                    u_it==round(value(u_it)),Xc==round(value(Xc)),Xm==round(value(Xm)), ...
                    I_ij==round(value(I_ij)),deta==round(value(deta))];
            elseif strcmp(sca_integer_stabilization,'network_after_first')
                plan_lock_cons=[plan_lock_cons, ...
                    u_it==round(value(u_it)),I_ij==round(value(I_ij)),deta==round(value(deta))];
            end
            lock_service_regions=~strcmp(sca_integer_stabilization,'electric_plan_after_first');
            if lock_service_regions&&strcmp(service_model,'fixed_pressure')&&strcmp(fixed_pressure_formulation,'compact_exact')
                plan_lock_cons=[plan_lock_cons, ...
                    fixed_pressure_g_available==round(value(fixed_pressure_g_available)), ...
                    fixed_pressure_w_available==round(value(fixed_pressure_w_available)), ...
                    fixed_pressure_gt_available==round(value(fixed_pressure_gt_available))];
            elseif lock_service_regions&&strcmp(service_model,'pressure_driven')&&strcmp(pressure_driven_formulation,'compact_exact')
                plan_lock_cons=[plan_lock_cons, ...
                    pressure_driven_g_active==round(value(pressure_driven_g_active)), ...
                    pressure_driven_w_active==round(value(pressure_driven_w_active)), ...
                    pressure_driven_gt_active==round(value(pressure_driven_gt_active))];
            elseif lock_service_regions
                plan_lock_cons=[plan_lock_cons, ...
                    vg1==round(value(vg1)),vg2==round(value(vg2)),vg3==round(value(vg3)), ...
                    vw1==round(value(vw1)),vw2==round(value(vw2)),vw3==round(value(vw3)), ...
                    vg1GT==round(value(vg1GT)),vg2GT==round(value(vg2GT)),vg3GT==round(value(vg3GT))];
            end
            plan_locked=true;
        end
        if sca_physics_accepted
            break
        elseif sca_iter>=sca_iter_max
            if enable_pipe_reverse_sca
                ordinary_water_pipe_rows=setdiff((1:size(WDS_pipe,1))',WDS_pump(:,1));
                ordinary_gas_pipe_rows=setdiff((1:size(GDS_pipe,1))',GDS_compessor(:,1));
                sca_diagnostic=struct;
                sca_diagnostic.rolling_step=t;
                sca_diagnostic.sca_iteration=sca_iter;
                sca_diagnostic.water_gap=water_soc_gap_values;
                sca_diagnostic.gas_gap=gas_soc_gap_values;
                sca_diagnostic.water_slack=[value(water_pipe_sasaP);value(water_pipe_sasaM)];
                sca_diagnostic.gas_slack=value(gas_pipe_sasa);
                sca_diagnostic.pump_slack=pump_slack_values;
                sca_diagnostic.water_pipe_rows=ordinary_water_pipe_rows;
                sca_diagnostic.water_pipe_data=WDS_pipe(ordinary_water_pipe_rows,:);
                sca_diagnostic.water_head=value(H_i);
                sca_diagnostic.water_head_min=H_imin(:,1:T);
                sca_diagnostic.water_head_max=H_imax(:,1:T);
                sca_diagnostic.water_qP=value(water_qP);
                sca_diagnostic.water_qM=value(water_qM);
                sca_diagnostic.water_hP=value(water_hP);
                sca_diagnostic.water_hM=value(water_hM);
                sca_diagnostic.water_D=value(D_ij);
                sca_diagnostic.water_L=value(L_ij);
                sca_diagnostic.water_flow_mean=value(Cw*W_ijbar);
                sca_diagnostic.water_flow_limit=ordinary_water_qmax;
                sca_diagnostic.water_linepack_change=value(X_ij-[X_ij0,X_ij(:,1:end-1)]);
                sca_diagnostic.water_source=value(W_WR);
                sca_diagnostic.water_service=value(w_i0(:,2:end));
                sca_diagnostic.gas_pipe_rows=ordinary_gas_pipe_rows;
                sca_diagnostic.gas_pipe_data=GDS_pipe(ordinary_gas_pipe_rows,:);
                sca_diagnostic.gas_pressure=value(kesai_i);
                sca_diagnostic.gas_pressure_min=kesai_imin(:,1:T);
                sca_diagnostic.gas_pressure_max=kesai_imax(:,1:T);
                sca_diagnostic.gas_pressure_high=value(kesai_ijH);
                sca_diagnostic.gas_pressure_low=value(kesai_ijL);
                sca_diagnostic.gas_scaled_flow=gas_scaled_flow_values;
                sca_diagnostic.gas_flow_mean=value(Bg2*G_ijave);
                sca_diagnostic.gas_flow_limit=Bg2*G_ijmax(:,1:T);
                sca_diagnostic.gas_linepack_change=value(X_mn-[X_mn0,X_mn(:,1:end-1)]);
                sca_diagnostic.gas_source=value(g_GS);
                sca_diagnostic.gas_service=value(g_i0(:,2:end));
                sca_diagnostic.water_pipe_sca_penalty=water_pipe_sca_penalty;
                sca_diagnostic.gas_pipe_sca_penalty=gas_pipe_sca_penalty;
                diagnostic_file=fullfile(output_dir,sprintf( ...
                    'sca_nonconvergence_roll_%02d_iter_%02d.mat',t,sca_iter));
                save(diagnostic_file,'sca_diagnostic');
                fprintf('Saved SCA nonconvergence diagnostic: %s\n',diagnostic_file);
            end
            error('ChanceRestoration:SCANonconvergence', ...
                ['Rolling step %d reached the FPP-SCA iteration cap (reference RMS %.3g/%.3g, ' ...
                'water/gas SOC gap %.3g/%.3g, water friction %.3g, pump slack/curve/power %.3g/%.3g/%.3g).'], ...
                t,ref_flow_rms,ref_speed_rms,water_soc_gap,gas_soc_gap,pipe_equation_residual, ...
                pump_slack,pump_curve_equation_residual,pump_power_equation_residual);
        else
            gigo3=min(gigo3max,pump_penalty_multiplier*gigo3);
            gigo4=min(gigo4max,pump_penalty_multiplier*gigo4);
            water_pipe_sca_penalty=min(pipe_sca_penalty_max,pipe_penalty_multiplier*water_pipe_sca_penalty);
            gas_pipe_sca_penalty=min(pipe_sca_penalty_max,pipe_penalty_multiplier*gas_pipe_sca_penalty);
        end
    end
    if sol.problem~=0
        error('ChanceRestoration:SolverFailure', ...
            'Rolling step %d failed with YALMIP code %d: %s',t,sol.problem,sol.info);
    end
    qbar_eval=value(Cw*W_ijbar);
    D_eval=value(D_ij);
    R_eval=Cw*R_w(:,1:T);
    friction_res_r(t)=max(abs(D_eval-R_eval.*qbar_eval.*abs(qbar_eval)),[],'all');
    hP_eval=value(water_hP);hM_eval=value(water_hM);
    qP_eval=value(water_qP);qM_eval=value(water_qM);
    water_soc_gap_r(t)=max([hP_eval-ordinary_water_kappa.*qP_eval.^2; ...
        hM_eval-ordinary_water_kappa.*qM_eval.^2],[],'all');
    gas_soc_gap_r(t)=gas_soc_gap;
    water_pipe_sca_slack_r(t)=water_pipe_slack;
    gas_pipe_sca_slack_r(t)=gas_pipe_slack;
    if strcmp(water_network_model,'transient')
        linepack_eval=value(Cw*(W_ij+W_ji)*Tm-(X_ij-[X_ij0,X_ij(:,1:end-1)]));
    else
        linepack_eval=value(Cw*(W_ij+W_ji));
    end
    linepack_res_r(t)=max(abs(linepack_eval),[],'all');
    if strcmp(gas_network_model,'transient')
        gas_dynamic_eval=value(Bg2*(G_ij+G_ji)*detT-(X_mn-[X_mn0,X_mn(:,1:end-1)]));
    else
        gas_dynamic_eval=value(Bg2*(G_ij+G_ji));
    end
    gas_dynamic_res_r(t)=max(abs(gas_dynamic_eval),[],'all');
    pump_curve_res_r(t)=pump_curve_equation_residual;
    pump_power_res_r(t)=pump_power_equation_residual;
    sca_iter_r(t)=sca_iter;sca_slack_r(t)=sca_slack;pump_slack_r(t)=pump_slack;solver_problem_r(t)=sol.problem;sca_plan_locked_r(t)=plan_locked;
    sca_ref_flow_rms_r(t)=ref_flow_rms;sca_ref_pump_rms_r(t)=ref_speed_rms;
    solve_time_r(t)=toc(step_tic);
    if sca_rolling_warm_start
        sca_flow_seed=zeros(size(WDS_pipe,1),1);
        sca_flow_seed(WDS_pump(:,1))=value(W_ij(WDS_pump(:,1),1));
        sca_speed_seed=pump_speed_sca(:,1);
    end
    if enable_pipe_reverse_sca
        gas_flow_seed=gas_scaled_flow_values(:,1);
        gas_low_pressure_seed=value(kesai_ijL(:,1));
    end
    % ====== 把窗口第1时段的结果“执行并回写为下一轮初始” ======
    % 状态类
    % --- 水/气/电 动态项：都要用第1列 ---
    H_i0   = value(H_i(:,1));          % 节点水压(用于下轮 L_ij 动力项里的 [H_i0, H_i(:,1:end-1)])
    H_0    = value(H_tank(:,1)); % 若你确实需要 H_0，可用 H_tank 第1列；如果你的 H_0 是选择水箱节点的水头，请按你之前定义一致地取第1列
    V_0=H_0*Sq;
    % kesai_i0 = value(kesai_i(:,1));    % 气网节点压力（用于气体动态方程的上一时刻）
    Sl_t0    = value(Sl_t(:,1));       % 储气罐库存（与 Sl_t 递推约束对齐）
    X_mn0    = value(X_mn(:,1));       % 气管存（与 X_mn 递推约束对齐）
    X_ij0    = value(X_ij(:,1));       % 水管存（与 X_ij 递推约束对齐）
    W_ij0    = value(W_ij(:,1));       % 管道首端流量
    W_ji0    = value(W_ji(:,1));       % 管道末端反向流量
    p_G0     = value(p_G(:,1));        % 燃气轮机电功率（用于 RU/RD）

    % --- 带“初始列”的连通/投切状态：下一时刻状态=第2列 ---
    y_it0    = round(value(y_it(:,2))); % 二进制状态按求解容差投影
    executed_switch_action=max(round(value(a_ij(:,2)))-a_ij0,0);
    switch_action_count_r(t)=sum(executed_switch_action);
    switch_action_weighted_r(t)=switch_action_weight'*executed_switch_action;
    a_ij0    = round(value(a_ij(:,2))); % 二进制状态按求解容差投影

    % --- 负荷恢复比例：也是“状态”，下一时刻应取第2列（保持单调不降） ---
    w_i0_0   = min(max(value(w_i0(:,2)),0),1);%连续状态投影到物理边界
    g_i0_0   = min(max(value(g_i0(:,2)),0),1);%连续状态投影到物理边界






    fprintf('Rolling step %d/%d done.\n', t, executed_steps);
    w_i0_r(:,t)=value(w_i0(:,2));H_tank_r(:,t)=value(H_tank(:,1));P_pump_r(:,t)=value(P_pump(:,1));omima_r(:,t)=pump_speed_sca(:,1);
    y_it_r(:,t)=value(y_it(:,2));r_i_r(:,t)=value(r_i(:,1));g_i0_r(:,t)=value(g_i0(:,2));W_WR_r(:,t)=value(W_WR(1,1));V_t_r(:,t)=value(V_t(:,1));
    W_ij_r(:,t)=value(W_ij(:,1));W_ji_r(:,t)=value(W_ji(:,1));X_ij_r(:,t)=value(X_ij(:,1));
    W_ijave_r(:,t)=value(W_ijbar(:,1));H_i_r(:,t)=value(H_i(:,1));L_ij_r(:,t)=value(L_ij(:,1));P_comp_r(:,t)=value(P_comp(:,1));
    water_hP_r(:,t)=value(water_hP(:,1));water_hM_r(:,t)=value(water_hM(:,1));
    water_qP_r(:,t)=value(water_qP(:,1));water_qM_r(:,t)=value(water_qM(:,1));
    Sl_t_r(:,t)=value(Sl_t(:,1));X_mn_r(:,t)=value(X_mn(:,1));p_G_r(:,t)=value(p_G(:,1));q_G_r(:,t)=value(q_G(:,1));g_GS_r(:,t)=value(g_GS(:,1));
    p_W_r(:,t)=value(p_W(:,1));a_ij_r(:,t)=value(a_ij(:,2));
    F_ij_r(:,t)=value(F_ij(:,1));v_i_r(:,t)=value(v_i(:,1));g_GT_r(:,t)=value(g_GT(:,1));kesai_i_r(:,t)=value(kesai_i(:,1));
    kesai_ijH_r(:,t)=value(kesai_ijH(:,1));kesai_ijL_r(:,t)=value(kesai_ijL(:,1));Xc_r(:,t)=value(Xc(:,1));Xm_r(:,t)=value(Xm(:,1));I_ij_r(:,t)=value(I_ij(:,1));
    G_ij_r(:,t)=value(G_ij(:,1));G_ji_r(:,t)=value(G_ji(:,1));G_ijave_r(:,t)=value(G_ijave(:,1));Gp_ij_r(:,t)=value(Gp_ij(:,1));
    g_sin_r(:,t)=value(g_sin(:,1)); g_sout_r(:,t)=value(g_sout(:,1)); q_W_r(:,t)=value(q_W(:,1));
    wind_lower_r(:,t)=wind_lower(:,1);wind_mean_r(:,t)=wind_mean(:,1);wind_sigma_r(:,t)=wind_sigma(:,1);
end
%% 求解

toc;
filename = fullfile(output_dir,['results_' service_model '_' water_network_model '_' ...
    gas_network_model '_' storage_model '_chance_' ...
    char(datetime('now','Format','yyyyMMdd_HHmmss')) '.mat']);

% 将变量统一保存到文件中
save(filename, ...
    'w_i0_r','H_tank_r', 'P_pump_r', 'omima_r', 'y_it_r', 'r_i_r', ...
    'g_i0_r', 'W_WR_r', 'V_t_r', 'W_ij_r', 'W_ji_r', 'X_ij_r', ...
    'W_ijave_r', 'H_i_r',  'L_ij_r', ...
    'P_comp_r', 'Sl_t_r', 'X_mn_r', 'p_G_r', 'q_G_r', 'g_GS_r', ...
    'p_W_r', 'a_ij_r', 'F_ij_r', 'v_i_r',  'g_GT_r',  ...
    'kesai_i_r', 'kesai_ijH_r', 'kesai_ijL_r', 'Xc_r', 'Xm_r', ...
    'I_ij_r', 'G_ij_r', 'G_ji_r', ...
    'G_ijave_r', 'Gp_ij_r', 'g_sin_r', 'g_sout_r', 'q_W_r', ...
    'wind_lower_r', 'wind_mean_r', 'wind_sigma_r', 'epsilon_wind', ...
    'box_coverage', 'wind_z', 'TT', 'Th', 'detT', 'case_file', 'H_source_initial', 'fix_water_source_head', 'sv', 'sw', 'gas_customer_min_pressure_bar', 'gas_source_data_nodes', 'gas_uniform_pressure_nodes', 'gas_soc_penalty_weight', 'enforce_gas_node_min_pressure', 'gas_customer_inactive_pressure_upper', 'water_friction_formulation', 'water_soc_penalty_weight', 'water_pipe_flow_capacity_m3s', 'enable_pipe_reverse_sca', 'pipe_sca_tol', 'pipe_sca_penalty_start', 'pipe_sca_reference_damping', 'kg', 'kw', 'omiga_g', 'water_hP_r', 'water_hM_r', 'water_qP_r', 'water_qM_r', 'water_soc_gap_r', 'gas_soc_gap_r', 'water_pipe_sca_slack_r', 'gas_pipe_sca_slack_r', 'kesai_i0', 'gas_source_nodes', 'gas_pressure_source_nodes', 'gas_source_pressure', 'gas_node1_source_capacity_kcfh', 'gas_flow_source_node', 'gas_flow_source_rate_kcfh', 'gas_flow_source_mode', 'G_smax', 'water_source_node', 'water_source_max', 'service_model', 'storage_model', 'water_tank_model', 'gas_tank_model', 'water_network_model', 'gas_network_model', 'fixed_pressure_formulation', 'pressure_driven_formulation', 'topology_formulation', 'topology_solution_mode', 'executed_steps', 'wind_capacity_boost_MW', 'wind_capacity_boost_unit', 'electric_load_reduction_MW', 'electric_load_before_reduction_MW', 'electric_load_after_reduction_MW', 'electric_load_scale', 'diagnostic_plan_file', 'diagnostic_plan_active', 'diagnostic_fix_scope', 'force_full_water_step', 'force_water_node_ids', 'switch_action_cost', 'tie_switch_multiplier', 'switch_action_count_r', 'switch_action_weighted_r', 'lookahead_start_step', 'binary_move_blocking', 'topology_move_blocking', 'flow_direction_move_blocking', 'water_region_move_blocking', 'gt_inactive_pressure_upper', 'sca_rolling_warm_start', 'sca_warm_start_mode', 'sca_warm_start_weight', 'pump_power_sca_mode', 'pump_speed_model', 'enforce_pump_physics', 'sca_iter_r', 'sca_slack_r', 'pump_slack_r', 'pump_curve_res_r', 'pump_power_res_r', 'sca_plan_locked_r', 'sca_tol_ref', 'sca_pump_ref_tol', 'sca_ref_flow_rms_r', 'sca_ref_pump_rms_r', 'sca_pump_tol', 'pump_power_tol', 'sca_penalty_start', 'solve_time_r', ...
    'solver_problem_r', 'friction_res_r', 'linepack_res_r', 'gas_dynamic_res_r', ...
    'sca_objective_trace_r', 'sca_ref_flow_rms_trace_r', 'sca_ref_pump_rms_trace_r', ...
    'sca_water_soc_gap_trace_r', 'sca_gas_soc_gap_trace_r', 'sca_pipe_residual_trace_r', 'sca_pump_slack_trace_r', ...
    'sca_water_pipe_slack_trace_r', 'sca_gas_pipe_slack_trace_r', ...
    'sca_pump_curve_residual_trace_r', 'sca_pump_power_residual_trace_r', ...
    'capture_solver_trace', 'sca_mip_warm_start', 'sca_integer_stabilization', 'sca_solve_time_trace_r', 'solver_output_trace_r', 'rpp_trace_r', 'mcb_trace_r');

disp(['✅ 所有变量已保存至文件：' filename]);
