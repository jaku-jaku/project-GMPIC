%% Helper Functions %%
classdef helper
    methods(Static)
        %% Logs:
        function logtitle(message)
            fprintf("\n[ ===== %s ===== ]:\n", message);
        end
        function loginfo(message)
            global glevel;
            if glevel < 1
                return;
            end
            fprintf("[INFO] %s\n", message);
        end
        function logerr(message)
            global glevel;
            if glevel < 2
                return;
            end
            fprintf("[-ERR] %s\n", message);
        end
        function logdebug(message)
            global glevel;
            if glevel < 3
                return
            end
            fprintf("[DEBUG] %s\n", message);
        end
        function logwrn(message)
            fprintf("[WARN] %s\n", message);
        end
        function logutil(message)
            fprintf("[HELPER]--> %s\n", message);
        end
        function setLogLevel(level)
            global glevel;
            switch level
                case "all"
                    glevel = 4;
                case "debug"
                    glevel = 3;
                case "error"
                    glevel = 2;
                case "info"
                    glevel = 1;
                otherwise
                    glevel = -1;
            end
            helper.logwrn(sprintf("> Setting Log Level @ [%s:%d]",level,glevel));
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
            global gdirectory;
            directory = sprintf("%s/%s", section, subsection);
            gdirectory = directory;
            fprintf("=========== [%s / %s] (init) ===========\n",  section, subsection);
            helper.createFolder(sprintf("output/%s", directory), clear_folder);
            if if_record
                helper.createDiary(directory, helper.getTimeStr());
            end
            fprintf("=========== [%s / %s] (start) ===========\n",  section, subsection);
        end
        function endSection(close_all)
            global gdirectory;
            fprintf("=========== [ End of Section : %s ] ===========\n", gdirectory)
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
        %% Video Rec:
        function createRecorder(FOLDER, FILE_NAME, QUALITY_, FPS_)
            global vObj;
            global vObj_exist;
            vObj_exist = true;
            EXP_PATH = sprintf('output/%s/%s.avi', FOLDER, FILE_NAME);
            vObj = VideoWriter(EXP_PATH);
            vObj.Quality = QUALITY_;
            vObj.FrameRate = FPS_;
            open(vObj);
        end
        function recordRecorder()
            global vObj;
            global vObj_exist;
            if vObj_exist
                writeVideo(vObj, getframe(gcf));
            end
        end
        function terminateRecorder()
            global vObj;
            global vObj_exist;
            if vObj_exist
                close(vObj);
                vObj_exist = false;
            end
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
            usr = getenv("USER");
            if usr == "jaku"
                addpath('~/JX-Platform/casadi-3');
                helper.logutil("Included Casadi-3!");
            else
                import casadi.*; % by default
            end
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
            % array or matrix to string for best latex
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
        function str = a2str(name, a)
            % array or matrix to string for best console logs
            if ~isa(a,'sym')
                a = sym(a);
            end
            a=vpa(simplify(a), 8);
            [n,m] = size(a);
            if m>1
                str = sprintf('\n%s ', join(string(a), ', '));
            else
                str = sprintf('%s ', join(string(a), ', '));
            end
            str = sprintf("%s=[ %s ]\n",name,str);
        end
        %% Progress bar:
        function textprogressbar(c)
            % REF: https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar
            % This function creates a text progress bar. It should be called with a 
            % STRING argument to initialize and terminate. Otherwise the number correspoding 
            % to progress in % should be supplied.
            % INPUTS:   C   Either: Text string to initialize or terminate 
            %                       Percentage number to show progress 
            % OUTPUTS:  N/A
            % Example:  Please refer to demo_textprogressbar.m
            % Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
            % Version: 1.0
            % Changes tracker:  29.06.2010  - First version
            % Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/
            % Initialization
            persistent strCR;           %   Carriage return pesistent variable
            % Vizualization parameters
            strPercentageLength = 10;   %   Length of percentage string (must be >5)
            strDotsMaximum      = 10;   %   The total number of dots in a progress bar
            % Main 
            if isempty(strCR) && ~ischar(c),
                % Progress bar must be initialized with a string
                error('The text progress must be initialized with a string');
            elseif isempty(strCR) && ischar(c),
                % Progress bar - initialization
                fprintf('%s',c);
                strCR = -1;
            elseif ~isempty(strCR) && ischar(c),
                % Progress bar  - termination
                strCR = [];  
                fprintf([c '\n']);
            elseif isnumeric(c)
                % Progress bar - normal progress
                c = floor(c);
                percentageOut = [num2str(c) '%%'];
                percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
                nDots = floor(c/100*strDotsMaximum);
                dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
                strOut = [percentageOut dotOut];
                
                % Print it on the screen
                if strCR == -1
                    % Don't do carriage return during first run
                    fprintf(strOut);
                else
                    % Do it during all the other runs
                    fprintf([strCR strOut]);
                end
                
                % Update carriage return
                strCR = repmat('\b',1,length(strOut)-1);
                
            else
                % Any other unexpected input
                error('Unsupported argument type');
            end
        end
    end
end