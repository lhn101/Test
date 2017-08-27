function Config_file_creator
module=1;
disp([ num2str(module) ' - Enter the analysis reference name [max 25 characters without spaces]'])
check=1;
while check
    name=input('        file name >> ','s');
    if ~ischar(name)
        disp('     wrong input -> not a character string')
        continue
    elseif length(name)>25
        disp('     wrong input -> string > 25 characters')
        continue
    elseif any(isspace(name))
        disp('     wrong input -> string contains spaces')
        continue
    else
        check=0;
    end
end

disp(' ')
module=module+1;
check=1;
while check
    nb_obs=input([ num2str(module) ' - How many time-series do you want to analyze simultanously? >> ']);
    if ischar(nb_obs)||rem(nb_obs,1)~=0
        disp('     wrong input -> not an integer number')
        continue
    else
        check=0;
    end
end

disp(' ')
module=module+1;
disp([ num2str(module) ' - Enter the time series reference names [max 10 characters]'])
for i=1:nb_obs
    check=1;
    while check                 
        data_labels{1,i}=input(['    time serie #' num2str(i) ' >> '],'s');
        if ~ischar(data_labels{1,i})
            disp('     wrong input -> not a character string')
            continue
        elseif length(data_labels{1,i})>10
            disp('     wrong input -> string > 10 characters ')
            continue
        else
            check=0;
        end
    end
end

if nb_obs>1
    disp(' ')
    module=module+1;
    disp([ num2str(module) ' - Identifies dependence between time series; use [] to indicate no dependence'])
    for i=1:nb_obs
        check=1;
        while check         
            comp_ic{1,i}=input(['    time serie #' num2str(i) ' depends on time series # >> ']);
            if ischar(comp_ic{1,i})
                disp('     wrong input -> should be integers')
                continue
            elseif length(comp_ic{1,i})>nb_obs-1
                disp('     wrong input -> invalid input')
                continue
            else
                check=0;
            end
        end
    end
else
    comp_ic={[]};
end

disp(' ')
module=module+1;
check=1;
while check
    nb_models=input([ num2str(module) ' - How many model classes do you want for each time-series? >> ']);
    if ischar(nb_models)||rem(nb_models,1)~=0
        disp('     wrong input -> not an integer number')
        continue
    else
        check=0;
    end
end
disp(' ')
disp('    --------------------------------------------------------')
disp('     BDLM Component reference numbers')
disp('    --------------------------------------------------------')
disp('     11: Local level')
disp('     12: Local trend')
disp('     13: Local acceleration')
disp('     21: Local level compatible with local trend')
disp('     22: Local level compatible with local acceleration')
disp('     23: Local trend compatible with local acceleration')
disp('     31: Perodic')
disp('     41: Autoregressive')
disp('     51: Dynamic regression with hidden component')
disp('     61: LSTM/Deep learning correction component')
disp('    --------------------------------------------------------')
disp(' ')

module=module+1;
disp([ num2str(module) ' - Identify components for each model class and observation; e.g. [11 31 41]'])
for j=1:nb_models
    if nb_models>1   
                disp(['    Model class #' num2str(j)])
    end
    for i=1:nb_obs
        check=1;
        while check
            comp{j}{i}=input(['     time serie #' num2str(i) ' >> ']);
            if ischar(comp{j}{i})
                disp('     wrong input -> should be integers')
                continue
            elseif j>1&& length(comp{j}{i})~=length(comp{j-1}{i})
                disp('     wrong input -> all model classes must have the same number of components')
                continue
            else
                check=0;
            end
        end
    end
end
disp(' ')
if nb_models>1
    module=module+1;
    disp([ num2str(module) ' - Identify shared parameters between the components of the model class #1; e.g. [0 1 1]'])
    for j=2:nb_models
        disp(['    Model class #' num2str(j)])
        for i=1:nb_obs
            check=1;
            while check
                const{j}{i}=input(['     time serie #' num2str(i) ' >> ']);
                if ischar(const{j}{i})
                    disp('     wrong input -> should be integers')
                    continue
                elseif length(const{j}{i})~=length(comp{j}{i})
                    disp('     wrong input -> the number of constraint parameters must be the same as the number of components')
                    continue
                else
                    check=0;
                end
            end
        end
    end
end
disp(' ')
module=module+1;
disp([ num2str(module) ' - Identify the number of days for the start and end of the training period'])
check=1;
while check          
    training(1)=input('            start >> ');
    training(2)=input('              end >> ');
    
    if ischar(training(1))
        disp('     wrong input -> should be integers')
        continue
    else
        check=0;
    end
end
cd config_files
fileID=fopen(['CFG_' name '.m'],'w');

fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fileID,'%%%% Project file_data reference\n');
fprintf(fileID,'option.name=''%s'';\n',name);
fprintf(fileID,'option.training_period=[%d,%d];\n\n',training);
fprintf(fileID,'%%%% Data\n');
fprintf(fileID,'data.labels={');
for i=1:nb_obs
    fprintf(fileID,'''%s''', data_labels{i});
    if i<nb_obs&&nb_obs>1
        fprintf(fileID,',');
    end
end
fprintf(fileID,'};\n');
fprintf(fileID,'data.values= %%Enter your dataset here; each column is a dataset; missing data -> NaN \n');
fprintf(fileID,'data.timestamps= %%Enter the time stamps (Unix format); column vector; \n\n');

fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fileID,'%% BDLM Component reference numbers\n');
fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fileID,'%% 11: Local level\n');
fprintf(fileID,'%% 12: Local trend\n');
fprintf(fileID,'%% 13: Local acceleration\n');
fprintf(fileID,'%% 21: Local level compatible with local trend\n');
fprintf(fileID,'%% 22: Local level compatible with local acceleration\n');
fprintf(fileID,'%% 23: Local trend compatible with local acceleration\n');
fprintf(fileID,'%% 31: Perodic\n');
fprintf(fileID,'%% 41: Autoregressive\n');
fprintf(fileID,'%% 51: Dynamic regression with hidden component\n');
fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
fprintf(fileID,'%%%% Model components | Choose which component to employ to build the model\n');
for j=1:nb_models
    fprintf(fileID,'model.components.block{%d}={',j);
    for i=1:nb_obs
        fprintf(fileID,'[');
        fprintf(fileID,'%d ',comp{j}{i});
        fprintf(fileID,'] ');
    end
    fprintf(fileID,'};\n');
end
fprintf(fileID,'\n');

fprintf(fileID,'%%%% Model component constrains | Take the same parameter as model class ##1\n');
for j=2:nb_models
    fprintf(fileID,'model.components.const{%d}={',j);
    for i=1:nb_obs
        fprintf(fileID,'[');
        fprintf(fileID,'%d ',const{j}{i});
        fprintf(fileID,'] ');
    end
    fprintf(fileID,'};\n\n');
end

fprintf(fileID,'%%%% Model inter-components dependence | {[components form dataset_i depends on components form dataset_j]_i,[...]}\n');
fprintf(fileID,'model.components.ic={');
for j=1:nb_obs
    fprintf(fileID,'[');
    fprintf(fileID,'%d ',comp_ic{1,j});
    fprintf(fileID,'] ');
end
fprintf(fileID,'};\n\n');

cd ..
disp(' ')
disp(['    The new configuration file was saved in .../BDLM_V2/config_files/' 'CFG_' name '.m'])

