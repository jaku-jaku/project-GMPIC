%% Helper Functions %%
classdef helper
    methods(Static)
        %% Logs:
        function logtitle(message)
            fprintf("\n[ ===== %s ===== ]:\n", message);
        end
        function loginfo(message)
            fprintf("[INFO] %s\n", message);
        end
        function logerr(message)
            fprintf("[-ERR] %s\n", message);
        end
        function logwrn(message)
            fprintf("[WARN] %s\n", message);
        end
        function logutil(message)
            fprintf("[HELPER]--> %s\n", message);
        end
        %% SYS:
        function createFolder(path, clear_folder)
            global gFigureIndex;
            gFigureIndex = 1;
            if ~exist(path)
                mkdir(path)
                helper.logutil("Folder created!");
            else
                if ~isempty(path) & clear_folder
                    rmdir(path, 's');
                    mkdir(path);
                    helper.logutil("Folder is emptied!");
                else
                    helper.logutil("Folder already existed, skipped!");
                end
            end 
        end
        function val=getTimeStr()
            val = datestr(datetime('now'), 'mm-dd_HH-MM-SS');
        end
        function directory = declareSection(section, subsection, ...
                if_record, clear_folder, close_window)
            if close_window
                close all;
                global gFigureIndex;
                gFigureIndex = 1;
            end
            directory = sprintf("%s/%s", section, subsection);
            fprintf("=========== [%s / %s] (init) ===========\n",  section, subsection);
            helper.createFolder(sprintf("output/%s", directory), clear_folder);
            if if_record
                helper.createDiary(directory, helper.getTimeStr());
            end
            fprintf("=========== [%s / %s] (start) ===========\n",  section, subsection);
        end
        function endSection(close_all)
            fprintf("=========== [ End of Section ] ===========\n")
            if close_all
                close all;
                global gFigureIndex;
                gFigureIndex = 1;
            end
            helper.stopDiary();
        end
        %% Utils:
        function rgb=cmap(i)
            color=[[0 0.4470 0.7410]
            [0.8500 0.3250 0.0980]
            [0.9290 0.6940 0.1250]
            [0.4940 0.1840 0.5560]
            [0.4660 0.6740 0.1880]
            [0.3010 0.7450 0.9330]
            [0.6350 0.0780 0.1840]];
            
            ind=mod(i, length(color))+1;
            rgb=color(ind,:);
        end
        function rgb=cmap_jet(ind,N)
            colormap jet
            color = colormap;
            rgb=color(ind*N,:);
        end
        function rgb=cmap_hsv(ind,N)
            colormap hsv
            color = colormap;
            rgb=color(ind*N,:);
        end
        %% Plot:
        function f=newFigure(index)
            global gFigureIndex;
            f = nan;
            if index < 0
                gFigureIndex = gFigureIndex+1;
                f = figure(gFigureIndex);
            elseif index == 1
                gFigureIndex = 1;
                f = figure(index);
            else
                f = figure(index);
            end
        end
        function plotLineTime(x,y,c,style)
            l = plot(x,y,style,'Linewidth',1,'color',c);
            % scatter(x,y,60,linspace(1,250,length(x)),"filled");
        end
        function annotate(fig,xy,offset,text,c)
            ha = annotation('textarrow','String',text, 'Color', c);
            ha.Parent=fig.CurrentAxes;
            ha.X = [xy(1)+offset(1) xy(1)];
            ha.Y = [xy(2)+offset(2) xy(2)];
        end
        %% Save:
        function saveFigure(DIMENSION, FOLDER, FILE_NAME)
            set(gcf,'units','points','position',[0, 0, DIMENSION(1), DIMENSION(2)]);
            EXP_PATH = sprintf('output/%s/%s.png', FOLDER, FILE_NAME);
            exportgraphics(gcf,EXP_PATH,'BackgroundColor','white');
            helper.logutil("PLOT SAVED @: "+ EXP_PATH);
        end
        function saveData(DATA, FOLDER, FILE_NAME)
            EXP_PATH = sprintf('output/%s/%s.txt', FOLDER, FILE_NAME);
            writematrix(DATA, EXP_PATH)
            helper.logutil("DATA SAVED @: "+ EXP_PATH);
        end
        %% Casadi:
        function includeCasadi()
            addpath('~/JX-Platform/casadi-3');
            helper.logutil("Included Casadi-3!");
        end
        %% Diary:
        function createDiary(FOLDER, FILE_NAME)
            OUTPUT_FILE = sprintf('output/%s/diary_%s.txt', FOLDER, FILE_NAME);
            helper.logutil("Outputting Console to: "+ OUTPUT_FILE);
            diary(OUTPUT_FILE)
        end
        function stopDiary()
            diary  off;
            helper.logutil("Diary OFF -x")
        end
        %% Latex:
        function latexcode = m2latex(matrix)
            if ~isa(matrix,'sym')
                matrix = sym(matrix);
            end
            latexcode = latex(vpa(simplify(matrix)));
            if numel(latexcode)>=6 && strcmp(latexcode(6),'(') && strcmp(latexcode(end),')')
                latexcode(6) = '[';
                latexcode(end) = ']';
            end
            % clipboard('copy',latexcode);
        end
    end
end