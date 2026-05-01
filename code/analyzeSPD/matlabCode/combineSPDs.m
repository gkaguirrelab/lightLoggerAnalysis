function combineSPDs(spds, output_path, options)
    arguments 
        spds
        output_path
        options.overwrite_existing = false 
        options.verbose = true; 
    end     

    % If given a path to the SPDs instead of the struct itself, 
    % read it in 
    if(isstring(spds) || ischar(spds))
        spds = load(spds).spds; 
    end 

    %{
    # Input dictionary should be of the fomr
    # {subject: {activity: 
    #               {color_mode: 
    #                           {projection_type: path}
    #                }
    #           }
    # }

    %}

    % Next, we will iterate over the subjects
    % and activities
    subjects = fieldnames(spds); 
    for ss = 1:numel(subjects)
        % Retrieve the subjet ID
        subject_id = subjects{ss}; 
        subject_struct = spds.(subject_id); 

        % Retrieve the activities to terate over 
        activity_names = fieldnames(subject_struct);
        for aa = 1:numel(activity_names)
            % Retrieve the activity name 
            activity_name = activity_names{aa}; 
            activity_struct = subject_struct.(activity_name); 
            
            % At this point, we should plot the SPD 
            % for each colormode on the same plot, including 
            % all of the available projection types







        end 




    end 


    

end


% Local function to plot the SPDs 
function localPlotSPD(activity_struct)


end