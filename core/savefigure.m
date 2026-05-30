function savefigure(fname,f)
%SAVEFIGURE Save a figure to file, defaulting to a dated figures/ folder.
%   SAVEFIGURE(FNAME) saves the current figure (gcf). SAVEFIGURE(FNAME, F)
%   saves figure handle F. If FNAME has no extension, both a .png and a .pdf
%   are written; if it has no path, output goes under ./figures/<yymmdd>/.
%   PDF output uses export_fig.

if ~exist('f','var')
    f = gcf;
    gcfflag = true;
else
    gcfflag = false;
end

[filepath,~,ext] = fileparts(fname);
if isempty(ext) || strcmp(ext,"")
    %if no extension is specified, save a png and pdf
    defaultflag = true;
    pdfname = strcat(fname,'.pdf');
    fname = strcat(fname,'.png');
else
    defaultflag = false;
end
if isempty(filepath) || length(filepath) == 1 %#ok<ISCL>
    %if no path is specified, save to a default folder under today's date
    baseDir = pwd; %current working folder
    % c = clock;
    t = datetime("today");
    %format = figures_yyyymmdd
    % filepath = fullfile(baseDir,sprintf('figures_%d%.2d%.2d',c(1),c(2),c(3)));

    %format = figures\yymmdd
    filepath = fullfile(baseDir,'figures',sprintf('%.2d%.2d%.2d',mod(year(t),1000),month(t),day(t)));
    
    %if no such folder exists, make it
    if ~exist(filepath,'dir'), mkdir(filepath); end
    
    fname = fullfile(filepath,fname);
end
%if the specified folder does not exist, make it
if ~exist(filepath,'dir'), mkdir(filepath); end

if strcmp(ext,'.pdf')
    export_fig(fname,'-nocrop')
else
    saveas(f,fname);
end

if gcfflag && defaultflag
    export_fig(pdfname,'-nocrop')
end

end