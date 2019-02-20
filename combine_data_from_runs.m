function varargout = combine_data_from_runs(base_dir, varargin)

if nargin > 1
    output_dir = varargin{1};
else
    output_dir = base_dir;
end

m = what(base_dir); % gets all matlab-related files (looking for .mat)

n_target_files = 0;
for i_file = 1:length(m.mat)
    file_name = lower(m.mat{i_file});
    if isequal(file_name(1), 's')
        %this is a subject data file
        n_target_files = n_target_files + 1;
    end
end

if n_target_files > 0
    file_list = cell(1, n_target_files);
    k_mat = 1;
    for i_file = 1:n_target_files
        file_name = lower(m.mat{i_file});
        if isequal(file_name(1), 's')
            %this is a subject data file
            file_list{k_mat} = file_name;
            k_mat = k_mat + 1;
        end
    end

    sub = file_list{1}(1:4);
    
    Data_main = struct('MT', [], 'RT', [], 'Succ', [], 'pPT', [], ...
        'Type', [], 'ViewTime', [], 'Kinematics', [], 'Target', []);
    
    for i_file = 1:n_target_files
        try
            load([base_dir, file_list{i_file}])
            Data1 = Data;
            clear Data;

            Data_main.MT = cat(1, Data_main.MT, Data1.MT);
            Data_main.RT = cat(1, Data_main.RT, Data1.RT);
            Data_main.Succ = cat(1, Data_main.Succ, Data1.Succ);
            Data_main.pPT = cat(1, Data_main.pPT, Data1.pPT);
            Data_main.Type = cat(1, Data_main.Type, Data1.Type);
            Data_main.ViewTime = cat(1, Data_main.ViewTime, Data1.ViewTime);
            Data_main.Kinematics = cat(1, Data_main.Kinematics, Data1.Kinematics);
            Data_main.Target = cat(1, Data_main.Target(:), Data1.Target(:));
            
            varargout{1} = 1;
            varargout{2} = [];
        catch err
            varargout{1} = -1;
            varargout{2} = err;
        end
        Data = Data_main;
        save([output_dir, filesep, sub], 'Data');
    end
else
    varargout{1} = 0;
    varargout{2} = [];
end