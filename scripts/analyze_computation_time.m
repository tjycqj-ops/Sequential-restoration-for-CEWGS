function outputs=analyze_computation_time(result_root,output_dir)
%ANALYZE_H2_COMPUTATION_TIME Export rolling-window and SCA timing records.

if nargin<2||isempty(output_dir),output_dir=pwd;end
if ~exist(output_dir,'dir'),mkdir(output_dir);end
case_keys={'proposed','fixed_pressure','water_steady','gas_steady'};
case_labels={'Proposed','Fixed-pressure','Steady-state water','Steady-state gas'};
colors=[0.90 0.31 0.10;0.00 0.38 0.70;0.32 0.58 0.25;0.48 0.31 0.68];
styles={'-','--',':','-.'};

window_table=table;
sca_table=table;
data=cell(size(case_keys));
for c=1:numel(case_keys)
    listing=dir(fullfile(result_root,case_keys{c},'results_*.mat'));
    assert(~isempty(listing),'Missing result for %s.',case_keys{c});
    [~,idx]=max([listing.datenum]);
    data{c}=load(fullfile(listing(idx).folder,listing(idx).name));
    assert(data{c}.Th==2&&data{c}.executed_steps==16);
    for t=1:data{c}.executed_steps
        sca_times=data{c}.sca_solve_time_trace_r{t};
        sca_times=sca_times(isfinite(sca_times));
        solver_sum=sum(sca_times);
        row=table(string(case_labels{c}),t,data{c}.solve_time_r(t), ...
            data{c}.sca_iter_r(t),solver_sum,data{c}.solve_time_r(t)-solver_sum, ...
            'VariableNames',{'case_name','rolling_window','window_time_s', ...
            'sca_iterations','sum_solver_time_s','modeling_overhead_s'});
        window_table=[window_table;row]; %#ok<AGROW>
        for k=1:numel(sca_times)
            detail=table(string(case_labels{c}),t,k,sca_times(k), ...
                'VariableNames',{'case_name','rolling_window','sca_iteration','solver_time_s'});
            sca_table=[sca_table;detail]; %#ok<AGROW>
        end
    end
end
writetable(window_table,fullfile(output_dir,'h2_window_timing.csv'));
writetable(sca_table,fullfile(output_dir,'h2_sca_iteration_timing.csv'));

summary=table;
for c=1:numel(case_keys)
    mask=window_table.case_name==string(case_labels{c});
    rows=window_table(mask,:);
    row=table(string(case_labels{c}),sum(rows.window_time_s), ...
        mean(rows.window_time_s),max(rows.window_time_s), ...
        sum(rows.sca_iterations),mean(rows.sca_iterations), ...
        'VariableNames',{'case_name','total_time_s','mean_window_time_s', ...
        'max_window_time_s','total_sca_iterations','mean_sca_iterations'});
    summary=[summary;row]; %#ok<AGROW>
end
writetable(summary,fullfile(output_dir,'h2_computation_summary.csv'));

fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 18.3 11.2]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);hold(ax,'on');box(ax,'on');grid(ax,'on');
handles=gobjects(1,numel(case_keys));
for c=1:numel(case_keys)
    mask=window_table.case_name==string(case_labels{c});
    rows=window_table(mask,:);
    handles(c)=plot(ax,rows.rolling_window,rows.window_time_s,styles{c}, ...
        'Color',colors(c,:),'LineWidth',1.7,'Marker','o','MarkerSize',3.2, ...
        'MarkerFaceColor','w');
end
yline(ax,120,'--','120-s limit','Color',[0.35 0.35 0.35], ...
    'LineWidth',1.0,'LabelHorizontalAlignment','right');
ylabel(ax,'Window time (s)');xlim(ax,[1 16]);xticks(ax,1:16);
ylim(ax,[0 125]);
text(ax,0.015,0.94,'(a)','Units','normalized','FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top');
lg=legend(ax,handles,case_labels,'Box','off','Orientation','horizontal', ...
    'Location','northoutside','FontSize',8.2);lg.NumColumns=4;

ax=nexttile(tl);hold(ax,'on');box(ax,'on');grid(ax,'on');
proposed=window_table(window_table.case_name=="Proposed",:);
bar(ax,proposed.rolling_window,proposed.sca_iterations,0.68, ...
    'FaceColor',colors(1,:),'EdgeColor','none','FaceAlpha',0.78);
ylabel(ax,'Proposed SCA iterations');xlabel(ax,'Rolling window');
xlim(ax,[1 16]);xticks(ax,1:16);ylim(ax,[0 max(proposed.sca_iterations)+1]);
text(ax,0.015,0.94,'(b)','Units','normalized','FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','top');

set(findall(fig,'Type','axes'),'FontName','Times New Roman','FontSize',9, ...
    'LineWidth',0.8,'TickDir','in','GridAlpha',0.16);
set(findall(fig,'Type','text'),'FontName','Times New Roman');
savefig(fig,fullfile(output_dir,'h2_computation_time.fig'));
exportgraphics(fig,fullfile(output_dir,'h2_computation_time.png'),'Resolution',600);
exportgraphics(fig,fullfile(output_dir,'h2_computation_time.pdf'),'ContentType','vector');
print(fig,fullfile(output_dir,'h2_computation_time.svg'),'-dsvg');
close(fig);

outputs=struct('summary',summary,'window_table',window_table, ...
    'sca_table',sca_table,'output_dir',output_dir);
disp(summary);
end

