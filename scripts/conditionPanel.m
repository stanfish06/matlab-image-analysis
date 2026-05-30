% some thoughts:
% condition panel is a grid with markers
% if some treatments are used, hightlight their markers and connect them with a Line
% treatment can also have dose and time
% for dose, we can use a color gradient
% for time, we can fill the disk with different angular fractions
% format for condition string: condition(dose,time)
% no unit for dose and time for now
% separate conditions by space?
function h = conditionPanel(conditions, opts)
% first extract individual conditions
conditions_list = cellfun(@(x)split(x),conditions,'UniformOutput',false);
conditions_dict = dictionary('dummy',struct);
% for each condition, extract its dose and time save
for i = 1:length(conditions_list)
    clear arr_temp
    for j = 1:length(conditions_list{i})
        condition_info = regexp(conditions_list{i}{j},'(\w+)\(([\d\-\.]+)([a-zA-Z]*),([\d\.\-]+)([a-zA-Z]*)\)','tokens');
        if isempty(condition_info)
            % if the condition is not in the right format, skip it
            continue;
        end
        name = condition_info{1}{1};
        dose = str2double(condition_info{1}{2});
        dose_unit = condition_info{1}{3};
        time = str2double(condition_info{1}{4});
        time_unit = condition_info{1}{5};
        if isnan(dose)
            dose = '';
        end
        if isnan(time)
            time = '';
        end
        arr_temp(j) = struct(...
            'name',name,...
            'dose',dose,...
            'time',time,...
            'dose_unit',dose_unit,...
            'time_unit',time_unit...
            );
        if ~conditions_dict.isKey(name)
            conditions_dict(name) = struct('dose',dose,'time',time);
        else
            % if the condition already exists, update the dose and time
            existing = conditions_dict(name);
            existing.dose = [existing.dose,dose];
            existing.time = [existing.time,time];
            conditions_dict(name) = existing;
        end
    end
    conditions_list{i} = arr_temp;
end
conditions_dict('dummy') = [];
colors = orderedcolors('gem12');
line_color = [59 59 59]/255;
negative_color = [233 233 233]/255;
n_treats = length(conditions_dict.keys);
n_conditions = length(conditions_list);
h = hggroup;
treats_order = sort(conditions_dict.keys);
disk_radius = 0.35;
num_points = 100;
min_alpha = 0.15;
theta = linspace(0, 2*pi, num_points);
dx = disk_radius * cos(theta);
dy = disk_radius * sin(theta);
y = 1:n_treats;
for i = 1:n_conditions
    mask = ismember(treats_order,{conditions_list{i}.name});
    % plot(repelem(i,length(n_treats)),1:n_treats,'o','Parent',h,'LineWidth',3,'MarkerFaceColor',negative_color, 'MarkerEdgeColor',negative_color,'MarkerSize',20);
    for j = 1:n_treats
        patch((i-1)*opts.aspect+1+dx*1.25+opts.hshift,n_treats-j+1+dy*1.25,negative_color,'Parent',h,'EdgeColor',negative_color);
    end
    plot(repelem((i-1)*opts.aspect+1+opts.hshift,sum(mask)),n_treats-y(mask)+1,'-','Parent',h,'LineWidth',opts.lineWidth,'Color',line_color);
    for j = y(mask)
        alpha_patch = 1 - (1 - min_alpha) * (1 - 1 / sum(conditions_list{i}(cellfun(@(x)strcmp(x,treats_order{j}),{conditions_list{i}.name})).dose <= unique(conditions_dict(treats_order{j}).dose)));
        if alpha_patch == Inf
            alpha_patch = 1;
        end
        patch((i-1)*opts.aspect+1+1.25*dx+opts.hshift,n_treats-j+1+1.25*dy,line_color,'Parent',h,'EdgeColor',line_color,'FaceAlpha',alpha_patch,'EdgeAlpha',alpha_patch);
        patch((i-1)*opts.aspect+1+dx+opts.hshift,n_treats-j+1+dy,colors(j,:),'Parent',h,'EdgeColor',colors(j,:),'FaceAlpha',alpha_patch,'EdgeAlpha',alpha_patch);

        text((i-1)*opts.aspect+1-0.85*disk_radius+opts.hshift,n_treats-j+1,[num2str(conditions_list{i}(cellfun(@(x)strcmp(x,treats_order{j}),{conditions_list{i}.name})).dose) conditions_list{i}(cellfun(@(x)strcmp(x,treats_order{j}),{conditions_list{i}.name})).dose_unit],'FontSize',opts.fontSizeCondition,'FontWeight','bold','Color','k');
        patch([(i-1)*opts.aspect+1+0.5*disk_radius+opts.hshift (i-1)*opts.aspect+1+0.5*disk_radius+opts.hshift (i-1)*opts.aspect+1+1.85*disk_radius+opts.hshift (i-1)*opts.aspect+1+1.85*disk_radius+opts.hshift],[n_treats-j+1-0.5*disk_radius n_treats-j+1+0.5*disk_radius n_treats-j+1+0.5*disk_radius n_treats-j+1-0.5*disk_radius],line_color,'Parent',h,'EdgeColor',line_color);
        text((i-1)*opts.aspect+1+0.5*disk_radius+opts.hshift,n_treats-j+1,[num2str(conditions_list{i}(cellfun(@(x)strcmp(x,treats_order{j}),{conditions_list{i}.name})).time) conditions_list{i}(cellfun(@(x)strcmp(x,treats_order{j}),{conditions_list{i}.name})).time_unit],'FontSize',opts.fontSizeTime,'FontWeight','bold','Color','white');
    end
end
axis equal;
xlim([1-opts.leftmargin 1+opts.aspect*(n_conditions-1)+opts.rightmargin]);
ylim([0.25 n_treats+0.75]);
xticklabels([]);
yticklabels(treats_order(end:-1:1));
ax = ancestor(h, 'axes');
ax.YAxis.FontSize = opts.fontSizeAxis;
ax.YAxis.FontWeight = 'bold';
ax.XTick = [];
ax.YTick = 1:n_treats;
box off;
set(get(ax, 'XAxis'), 'Visible', 'off');
ax.LineWidth = 2;  % Set axis line width
ax.XAxis.FontWeight = 'bold';  % x-axis tick labels
ax.YAxis.FontWeight = 'bold';  % y-axis tick labels
ax.XAxis.FontSize = opts.fontSizeAxis;  % x-axis tick labels
ax.YAxis.FontSize = opts.fontSizeAxis;  % y-axis tick labels