function outputs = plot_four_case_source_csv(source_csv,output_dir)
%PLOT_FOUR_CASE_SOURCE_CSV Plot saved four-case source data without optimization.

validateattributes(source_csv,{'char','string'},{'scalartext'});
validateattributes(output_dir,{'char','string'},{'scalartext'});
source_csv = char(source_csv);
output_dir = char(output_dir);
assert(exist(source_csv,'file')==2,'Missing plotting source data: %s',source_csv);
if ~exist(output_dir,'dir'), mkdir(output_dir); end

source = readtable(source_csv,'TextType','string');
required = {'case_name','time_h','electric_supply_MW', ...
    'electric_priority_recovered_pct','water_supply_m3s', ...
    'water_priority_recovered_pct','gas_supply_kcfh', ...
    'gas_priority_recovered_pct'};
assert(all(ismember(required,source.Properties.VariableNames)), ...
    'The plotting source CSV does not contain the required columns.');

case_names = ["Proposed","Fixed-pressure","Steady-state water","Steady-state gas"];
case_labels = {'Proposed','Fixed-pressure','Steady-water','Steady-gas'};
colors = [0.90 0.31 0.10;0.00 0.38 0.70;0.32 0.58 0.25;0.48 0.31 0.68];
styles = {'-','--',':','-.'};
supply_fields = {'electric_supply_MW','water_supply_m3s','gas_supply_kcfh'};
recovery_fields = {'electric_priority_recovered_pct', ...
    'water_priority_recovered_pct','gas_priority_recovered_pct'};
supply_labels = {'Total generation (MW)','Water source supply (m^3/s)', ...
    'Gas source supply (kcf/h)'};
recovery_labels = {'Electric load recovered (%)','Water load recovered (%)', ...
    'Gas load recovered (%)'};
panel_labels = {'(a)','(b)','(c)'};

series = cell(1,numel(case_names));
for case_index = 1:numel(case_names)
    rows = source.case_name == case_names(case_index);
    assert(any(rows),'Missing case in plotting source: %s',case_names(case_index));
    series{case_index} = sortrows(source(rows,:),'time_h');
end
time_h = series{1}.time_h;
assert(all(cellfun(@(block) isequal(block.time_h,time_h),series)), ...
    'The four cases do not share the same time grid.');

fig = figure('Visible','off','Color','w','Units','centimeters', ...
    'Position',[2 2 18.3 15.4]);
layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
line_handles = gobjects(1,numel(case_names));
bar_handles = gobjects(1,numel(case_names));
for panel_index = 1:3
    ax = nexttile(layout); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    supply = zeros(numel(time_h),numel(case_names));
    for case_index = 1:numel(case_names)
        supply(:,case_index) = series{case_index}.(supply_fields{panel_index});
    end
    yyaxis(ax,'left');
    bars = bar(ax,time_h,supply,0.78,'grouped','EdgeColor','none','FaceAlpha',0.35);
    for case_index = 1:numel(case_names)
        bars(case_index).FaceColor = colors(case_index,:);
        bar_handles(case_index) = bars(case_index);
    end
    ylabel(ax,supply_labels{panel_index});

    yyaxis(ax,'right');
    for case_index = 1:numel(case_names)
        line_handles(case_index) = plot(ax,time_h, ...
            series{case_index}.(recovery_fields{panel_index}),styles{case_index}, ...
            'Color',colors(case_index,:),'LineWidth',1.55,'Marker','o', ...
            'MarkerIndices',1:2:numel(time_h),'MarkerSize',3, ...
            'MarkerFaceColor','w');
    end
    ylabel(ax,recovery_labels{panel_index}); ylim(ax,[0 105]); yticks(ax,0:25:100);
    xlim(ax,[min(time_h),max(time_h)]); xticks(ax,min(time_h):0.5:max(time_h));
    text(ax,0.015,0.94,panel_labels{panel_index},'Units','normalized', ...
        'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top');
    if panel_index == 1
        legend(ax,[line_handles,bar_handles], ...
            [strcat(case_labels,' recovery'),strcat(case_labels,' supply')], ...
            'Box','off','Orientation','horizontal','Location','northoutside', ...
            'FontSize',7.3,'NumColumns',4);
    end
end
xlabel(ax,'Time (h)');
set(findall(fig,'Type','axes'),'FontName','Times New Roman','FontSize',9, ...
    'LineWidth',0.8,'TickDir','in','GridAlpha',0.16);
set(findall(fig,'Type','text'),'FontName','Times New Roman');

stem = 'four_case_low_pressure_recovery_supply';
fig_file = fullfile(output_dir,[stem '.fig']);
png_file = fullfile(output_dir,[stem '.png']);
tiff_file = fullfile(output_dir,[stem '.tiff']);
pdf_file = fullfile(output_dir,[stem '.pdf']);
svg_file = fullfile(output_dir,[stem '.svg']);
savefig(fig,fig_file);
exportgraphics(fig,png_file,'Resolution',600);
exportgraphics(fig,tiff_file,'Resolution',600);
exportgraphics(fig,pdf_file,'ContentType','vector');
print(fig,svg_file,'-dsvg');
close(fig);

outputs = struct('source_csv',source_csv,'output_dir',output_dir, ...
    'fig',fig_file,'png',png_file,'tiff',tiff_file,'pdf',pdf_file,'svg',svg_file);
fprintf('[PLOT-ONLY] Source: %s\n[PLOT-ONLY] Output: %s\n',source_csv,output_dir);
end
