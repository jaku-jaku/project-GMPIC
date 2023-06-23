%% Helper Functions %%
classdef helper
    methods(Static)
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
        function createFolder(path, clear_folder)
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
        function sisoPlot(L_TF, DIMENSION, FOLDER, TAG)
            % Plot
            figure()

            subplot(2, 2, [1,3])
            margin(L_TF)
            grid on
            
            subplot(2, 2, 2)
            rlocus(L_TF)
            
            subplot(2, 2, 4)
            G_TF = minreal(L_TF/(L_TF + 1));
            step(G_TF)
            grid on

            helper.saveFigure(DIMENSION, FOLDER, sprintf("siso_plot_%s", TAG))
        end
        function SO = run_simulink(model_name, duration)
            %% Run simulation model
            MODEL_FILE_PATH = "Models/" + model_name + ".slx";
            open(MODEL_FILE_PATH); % opens the model file
            set_param(model_name, 'StopTime', int2str(duration));%set stop time to 4 s
            SO = sim(MODEL_FILE_PATH); % runs the simulation
            helper.logutil("Running Model ["+ model_name +"] for "+ duration +" s");
        end
        function includeCasadi()
            addpath('~/JX-Platform/casadi-3');
            helper.logutil("Included Casadi-3!");
        end
        function createDiary(FOLDER, FILE_NAME)
            OUTPUT_FILE = sprintf('output/%s/diary_%s.txt', FOLDER, FILE_NAME);
            helper.logutil("Outputting Console to: "+ OUTPUT_FILE);
            diary(OUTPUT_FILE)
        end
    end
end