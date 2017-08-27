function file_data=load_saved_files

files=dir(strcat(cd,'/saved_files/*'));
files_2=zeros(size(files,1),1);
for i_1=1:size(files,1)
    files_2(i_1) = strcmp(files(i_1).name,'file_data.mat'); %identify if file_data.mat exists
end
if any(files_2)
    load(strcat(cd,'/saved_files/','file_data.mat'));      %Choose file to load
    if size(file_data,1)>1
        for loop=2:size(file_data,1)
            add_char=repmat(' ',[1,max(66-length([file_data{loop,1},'-',file_data{loop,2}]),0)]);
            if loop<=10
                disp([num2str(loop-1),'  -> ',file_data{loop,1},add_char, ' /' file_data{loop,2}]);
            else
                disp([num2str(loop-1),' -> ',file_data{loop,1},add_char,' /' file_data{loop,2}]);
            end
            
        end
    end
else
    disp('file_data.mat does not exists -> a new empty one was created')
    file_data=cell(1,2);
end