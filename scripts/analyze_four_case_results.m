function outputs=analyze_four_case_results(result_root,case_file,output_dir,expected_horizon)
%ANALYZE_FOUR_CASE_LOW_PRESSURE_H1 Matched low-pressure comparison.
% The plotted recovery metric is priority weighted, matching the restoration
% objective. Supply trajectories and demand-weighted recovery are exported as
% source data so that the visual comparison remains auditable.

if nargin<3||isempty(output_dir),output_dir=pwd;end
if nargin<4||isempty(expected_horizon),expected_horizon=1;end
if ~exist(output_dir,'dir'),mkdir(output_dir);end

case_keys={'proposed','fixed_pressure','water_steady','gas_steady'};
case_labels={'Proposed','Fixed-pressure','Steady-state water','Steady-state gas'};
expected_service={'pressure_driven','fixed_pressure','pressure_driven','pressure_driven'};
expected_water={'transient','transient','steady_state','transient'};
expected_gas={'transient','transient','transient','steady_state'};
expected_water_tank={'dynamic_storage','dynamic_storage','no_storage','dynamic_storage'};
expected_gas_tank={'dynamic_storage','dynamic_storage','dynamic_storage','no_storage'};
expected_binary_linked=false(1,4);
expected_topology_linked=true(1,4);
expected_flow_direction_linked=true(1,4);
expected_water_region_linked=false(1,4);
colors=[0.90 0.31 0.10;0.00 0.38 0.70;0.32 0.58 0.25;0.48 0.31 0.68];
styles={'-','--',':','-.'};
priority=[2,2,100,10,10,2,100,2,100,2,2,10,10,2]';
ncase=numel(case_keys);
if isscalar(expected_horizon)
    expected_horizon=repmat(expected_horizon,1,ncase);
else
    expected_horizon=reshape(expected_horizon,1,[]);
    assert(numel(expected_horizon)==ncase, ...
        'Expected horizon must be scalar or contain one value per case.');
end
d=cell(1,ncase);result_files=cell(1,ncase);

for c=1:ncase
    listing=dir(fullfile(result_root,case_keys{c},'results_*.mat'));
    assert(~isempty(listing),'Missing result for %s.',case_keys{c});
    [~,idx]=max([listing.datenum]);
    result_files{c}=fullfile(listing(idx).folder,listing(idx).name);
    d{c}=load(result_files{c});
    assert(d{c}.TT==16&&d{c}.executed_steps==16&&d{c}.Th==expected_horizon(c), ...
        '%s is not a completed H=%d/16-period result.',case_labels{c},expected_horizon(c));
    assert(all(d{c}.solver_problem_r==0),'%s contains a failed window.',case_labels{c});
    assert(strcmp(d{c}.service_model,expected_service{c}));
    assert(strcmp(d{c}.water_network_model,expected_water{c}));
    assert(strcmp(d{c}.gas_network_model,expected_gas{c}));
    assert(strcmp(d{c}.water_tank_model,expected_water_tank{c}));
    assert(strcmp(d{c}.gas_tank_model,expected_gas_tank{c}));
    assert(abs(d{c}.sw-0.4)<1e-12&&abs(d{c}.sv-0.4)<1e-12);
    assert(abs(d{c}.water_soc_penalty_weight-0.001)<1e-12);
    assert(abs(d{c}.gas_soc_penalty_weight-0.1)<1e-12);
    assert(strcmp(d{c}.gas_flow_source_mode,'upper_bound'));
    assert(~d{c}.gt_inactive_pressure_upper);
    if expected_horizon(c)>1
        assert(logical(d{c}.binary_move_blocking)==expected_binary_linked(c)&& ...
            logical(d{c}.topology_move_blocking)==expected_topology_linked(c)&& ...
            logical(d{c}.flow_direction_move_blocking)==expected_flow_direction_linked(c), ...
            '%s does not use the matched H>1 hybrid binary contract.',case_labels{c});
        assert(logical(d{c}.water_region_move_blocking)==expected_water_region_linked(c), ...
            '%s has an unexpected water-pressure-region linking setting.',case_labels{c});
    end
end

casedata=load(case_file);
electric_case=casedata.load;
electric_case(:,2:3)=0.5*electric_case(:,2:3);
electric_case(9,2:3)=electric_case(8,2:3);
electric_case(8,2:3)=0;
electric_demand=electric_case(electric_case(:,2)~=0,2);
water_demand=casedata.WDS_load(2:end,2);
gas_demand=casedata.GDS_node(casedata.GDS_node(:,6)~=0,6);
gas_demand_nodes=find(casedata.GDS_node(:,6)~=0);
water_req=3*casedata.WDS_node(2:end,3);
gas_req=30*ones(numel(gas_demand_nodes),1);

n=16;time_h=(0:n)*0.25;
recovery_priority=zeros(ncase,3,n+1);
recovery_demand=zeros(ncase,3,n+1);
supply=zeros(ncase,3,n+1);
pressure_ready=zeros(ncase,2,n+1);
underpressure_service=zeros(ncase,2,n+1);

for c=1:ncase
    recovery_priority(c,1,2:end)=reshape(100*(priority'*d{c}.r_i_r)/sum(priority),1,1,[]);
    recovery_priority(c,2,2:end)=reshape(100*(priority'*d{c}.w_i0_r)/sum(priority),1,1,[]);
    recovery_priority(c,3,2:end)=reshape(100*(priority'*d{c}.g_i0_r)/sum(priority),1,1,[]);
    recovery_demand(c,1,2:end)=reshape(100*(electric_demand'*d{c}.r_i_r)/sum(electric_demand),1,1,[]);
    recovery_demand(c,2,2:end)=reshape(100*(water_demand'*d{c}.w_i0_r)/sum(water_demand),1,1,[]);
    recovery_demand(c,3,2:end)=reshape(100*(gas_demand'*d{c}.g_i0_r)/sum(gas_demand),1,1,[]);
    supply(c,1,2:end)=reshape(sum(d{c}.p_G_r,1)+sum(casedata.WT(:,3).*d{c}.p_W_r,1),1,1,[]);
    supply(c,2,2:end)=reshape(d{c}.W_WR_r,1,1,[]);
    supply(c,3,2:end)=reshape(sum(d{c}.g_GS_r,1),1,1,[]);
    pressure_ready(c,1,2:end)=reshape(100*(priority'*(d{c}.H_i_r(2:end,:)>=water_req))/sum(priority),1,1,[]);
    pressure_ready(c,2,2:end)=reshape(100*(priority'*(d{c}.kesai_i_r(gas_demand_nodes,:)>=gas_req))/sum(priority),1,1,[]);
    underpressure_service(c,1,2:end)=reshape(100*(priority'*(d{c}.w_i0_r.* ...
        (d{c}.H_i_r(2:end,:)<water_req)))/sum(priority),1,1,[]);
    underpressure_service(c,2,2:end)=reshape(100*(priority'*(d{c}.g_i0_r.* ...
        (d{c}.kesai_i_r(gas_demand_nodes,:)<gas_req)))/sum(priority),1,1,[]);
end

source=table;
for c=1:ncase
    block=table(repmat(string(case_labels{c}),n+1,1),time_h.', ...
        squeeze(supply(c,1,:)),squeeze(recovery_priority(c,1,:)),squeeze(recovery_demand(c,1,:)), ...
        squeeze(supply(c,2,:)),squeeze(recovery_priority(c,2,:)),squeeze(recovery_demand(c,2,:)), ...
        squeeze(supply(c,3,:)),squeeze(recovery_priority(c,3,:)),squeeze(recovery_demand(c,3,:)), ...
        squeeze(pressure_ready(c,1,:)),squeeze(pressure_ready(c,2,:)), ...
        squeeze(underpressure_service(c,1,:)),squeeze(underpressure_service(c,2,:)), ...
        'VariableNames',{'case_name','time_h','electric_supply_MW', ...
        'electric_priority_recovered_pct','electric_demand_recovered_pct', ...
        'water_supply_m3s','water_priority_recovered_pct','water_demand_recovered_pct', ...
        'gas_supply_kcfh','gas_priority_recovered_pct','gas_demand_recovered_pct', ...
        'water_full_pressure_ready_pct','gas_full_pressure_ready_pct', ...
        'water_underpressure_service_pct','gas_underpressure_service_pct'});
    source=[source;block]; %#ok<AGROW>
end
writetable(source,fullfile(output_dir,'four_case_low_pressure_source.csv'));

summary=table;
for c=1:ncase
    normal_g=setdiff((1:size(casedata.GDS_pipe,1))',casedata.GDS_compessor(:,1));
    cij=casedata.GDS_pipe(normal_g,4);
    gas_gap=d{c}.kesai_ijH_r.^2-d{c}.kesai_ijL_r.^2- ...
        (d{c}.G_ijave_r(normal_g,:)./cij).^2;
    t95_e=time_to_threshold(time_h,squeeze(recovery_priority(c,1,:)),95);
    t95_w=time_to_threshold(time_h,squeeze(recovery_priority(c,2,:)),95);
    t95_g=time_to_threshold(time_h,squeeze(recovery_priority(c,3,:)),95);
    row=table(string(case_labels{c}), ...
        recovery_priority(c,1,end),recovery_priority(c,2,end),recovery_priority(c,3,end), ...
        mean(recovery_priority(c,1,2:end)),mean(recovery_priority(c,2,2:end)),mean(recovery_priority(c,3,2:end)), ...
        mean(underpressure_service(c,1,2:end)),max(underpressure_service(c,1,2:end)), ...
        mean(underpressure_service(c,2,2:end)),max(underpressure_service(c,2,2:end)), ...
        sum(abs(diff(squeeze(recovery_priority(c,2,:))))), ...
        sum(abs(diff(squeeze(recovery_priority(c,3,:))))), ...
        t95_e,t95_w,t95_g,sum(d{c}.solve_time_r),max(d{c}.solve_time_r),sum(d{c}.sca_iter_r), ...
        mean(d{c}.friction_res_r),max(d{c}.friction_res_r), ...
        max(d{c}.linepack_res_r),max(d{c}.gas_dynamic_res_r), ...
        mean(abs(gas_gap),'all'),max(abs(gas_gap),[],'all'),max(d{c}.pump_power_res_r), ...
        'VariableNames',{'case_name','final_electric_priority_pct','final_water_priority_pct', ...
        'final_gas_priority_pct','mean_electric_priority_pct','mean_water_priority_pct', ...
        'mean_gas_priority_pct','mean_water_underpressure_service_pct', ...
        'max_water_underpressure_service_pct','mean_gas_underpressure_service_pct', ...
        'max_gas_underpressure_service_pct','water_total_variation_pp', ...
        'gas_total_variation_pp','time_to_95_electric_h','time_to_95_water_h', ...
        'time_to_95_gas_h','total_solve_time_s','max_window_time_s','total_sca_iterations', ...
        'mean_window_max_water_friction_residual_m','max_water_friction_residual_m', ...
        'max_water_storage_residual','max_gas_dynamic_residual', ...
        'mean_gas_soc_gap_bar2','max_gas_soc_gap_bar2','max_pump_power_residual_kW'});
    summary=[summary;row]; %#ok<AGROW>
end
writetable(summary,fullfile(output_dir,'four_case_low_pressure_summary.csv'));

fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18.3 14.8]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
ylabels={'Electric load recovered (%)','Water load recovered (%)','Gas load recovered (%)'};
panels={'(a)','(b)','(c)'};handles=gobjects(1,ncase);
for p=1:3
    ax=nexttile(tl);hold(ax,'on');box(ax,'on');grid(ax,'on');
    for c=1:ncase
        handles(c)=plot(ax,time_h,squeeze(recovery_priority(c,p,:)),styles{c}, ...
            'Color',colors(c,:),'LineWidth',1.7,'Marker','o','MarkerIndices',1:2:n+1, ...
            'MarkerSize',3.2,'MarkerFaceColor','w');
    end
    ylabel(ax,ylabels{p});ylim(ax,[0 105]);yticks(ax,0:25:100);
    xlim(ax,[0 4]);xticks(ax,0:0.5:4);
    text(ax,0.015,0.94,panels{p},'Units','normalized','FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    if p==1
        lg=legend(ax,handles,case_labels,'Box','off','Orientation','horizontal', ...
            'Location','northoutside','FontSize',8.4);lg.NumColumns=4;
    end
end
xlabel(ax,'Time (h)');style_axes(fig);
save_figure(fig,output_dir,'four_case_low_pressure_recovery');close(fig);

figp=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18.3 11.5]);
tl=tiledlayout(figp,2,1,'TileSpacing','compact','Padding','compact');
pressure_labels={'Water below threshold (%)', ...
    'Gas below threshold (%)'};
for p=1:2
    ax=nexttile(tl);hold(ax,'on');box(ax,'on');grid(ax,'on');
    for c=1:ncase
        handles(c)=plot(ax,time_h,squeeze(underpressure_service(c,p,:)),styles{c}, ...
            'Color',colors(c,:),'LineWidth',1.7,'Marker','o','MarkerIndices',1:2:n+1, ...
            'MarkerSize',3.2,'MarkerFaceColor','w');
    end
    pressure_upper=max(50,10*ceil(max(underpressure_service(:,p,:),[],'all')/10));
    ylabel(ax,pressure_labels{p});ylim(ax,[0 pressure_upper]);yticks(ax,0:10:pressure_upper);
    xlim(ax,[0 4]);xticks(ax,0:0.5:4);
    text(ax,0.015,0.94,panels{p},'Units','normalized','FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    if p==1
        lg=legend(ax,handles,case_labels,'Box','off','Orientation','horizontal', ...
            'Location','northoutside','FontSize',8.4);lg.NumColumns=4;
    end
end
xlabel(ax,'Time (h)');style_axes(figp);
save_figure(figp,output_dir,'four_case_low_pressure_underpressure_service');close(figp);

fig2=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18.3 15.4]);
tl=tiledlayout(fig2,3,1,'TileSpacing','compact','Padding','compact');
supply_labels={'Total generation (MW)','Water source supply (m^3/s)','Gas source supply (kcf/h)'};
line_handles=gobjects(1,ncase);bar_handles=gobjects(1,ncase);
for p=1:3
    ax=nexttile(tl);hold(ax,'on');box(ax,'on');grid(ax,'on');
    yyaxis(ax,'left');
    bars=bar(ax,time_h,squeeze(supply(:,p,:)).',0.78,'grouped','EdgeColor','none','FaceAlpha',0.35);
    for c=1:ncase,bars(c).FaceColor=colors(c,:);bar_handles(c)=bars(c);end
    ylabel(ax,supply_labels{p});
    yyaxis(ax,'right');
    for c=1:ncase
        line_handles(c)=plot(ax,time_h,squeeze(recovery_priority(c,p,:)),styles{c}, ...
            'Color',colors(c,:),'LineWidth',1.55,'Marker','o','MarkerIndices',1:2:n+1, ...
            'MarkerSize',3,'MarkerFaceColor','w');
    end
    ylabel(ax,ylabels{p});ylim(ax,[0 105]);yticks(ax,0:25:100);
    xlim(ax,[0 4]);xticks(ax,0:0.5:4);
    text(ax,0.015,0.94,panels{p},'Units','normalized','FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','top');
    if p==1
        lg=legend(ax,[line_handles,bar_handles], ...
            {'Proposed recovery','Fixed-pressure recovery','Steady-water recovery','Steady-gas recovery', ...
            'Proposed supply','Fixed-pressure supply','Steady-water supply','Steady-gas supply'}, ...
            'Box','off','Orientation','horizontal','Location','northoutside','FontSize',7.3);lg.NumColumns=4;
    end
end
xlabel(ax,'Time (h)');style_axes(fig2);
save_figure(fig2,output_dir,'four_case_low_pressure_recovery_supply');close(fig2);

save(fullfile(output_dir,'four_case_low_pressure_data.mat'),'result_files','case_keys', ...
    'case_labels','time_h','recovery_priority','recovery_demand','supply','pressure_ready', ...
    'underpressure_service','summary');
outputs=struct('summary',summary,'source',source,'result_files',{result_files},'output_dir',output_dir);
disp(summary);
fprintf('[FOUR-CASE-LOW-PRESSURE] output=%s H=%s excluded=0 periods=%d cases=%d\n', ...
    output_dir,mat2str(expected_horizon),n,ncase);
end

function value=time_to_threshold(time_h,series,threshold)
idx=find(series>=threshold-1e-8,1,'first');
if isempty(idx),value=NaN;else,value=time_h(idx);end
end

function style_axes(fig)
set(findall(fig,'Type','axes'),'FontName','Times New Roman','FontSize',9, ...
    'LineWidth',0.8,'TickDir','in','GridAlpha',0.16);
set(findall(fig,'Type','text'),'FontName','Times New Roman');
end

function save_figure(fig,output_dir,stem)
savefig(fig,fullfile(output_dir,[stem '.fig']));
exportgraphics(fig,fullfile(output_dir,[stem '.png']),'Resolution',600);
exportgraphics(fig,fullfile(output_dir,[stem '.pdf']),'ContentType','vector');
print(fig,fullfile(output_dir,[stem '.svg']),'-dsvg');
end

