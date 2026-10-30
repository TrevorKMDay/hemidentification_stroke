% DATA FROM NATALYA 
% 
% file = '/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/data/StrokeRHFCpaper_gruffalos_20260506.mat' ;
% 
% data = load(file) ;
% 
% % fc1 = data.gruffalos.fc{1, 1} ;
% 
% for i = 1:size(data.gruffalos, 1)
% 
%     id = data.gruffalos.subname{i} ;
% 
%     unpacked = data.gruffalos.fc{i} ;
%     writematrix(unpacked, ...
%         strcat("/Users/tkmd/MyDrive/Projects/", ...
%                 "hemisphere_fingerprinting_stroke/data/", ...
%                 "individual_matrices/sub-", id, ".csv"))
% 
% end

% DATA FROM ANDREW

% Original file
% file = '/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/data/brainnetome274_scfc.mat' ;

% Fixed file from Kyle,
file = '/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/data/brainnetome274_scfc_20260922.mat' ;
load(file)

fcs = m3(m3.Unit == "fc", :) ;

lesions = m3(m3.Unit == "lesion", ["SubjectCode", "ScanDate", "Data"]) ;

writetable(lesions, ...
    strcat("/Users/tkmd/MyDrive/Projects/", ...
           "hemisphere_fingerprinting_stroke/data/lesions")) ;

fc_date = fcs(:, ["SubjectCode", "ScanDate"]) ;
writetable(fc_date, ...
            strcat("/Users/tkmd/MyDrive/Projects/", ...
                   "hemisphere_fingerprinting_stroke/data/", ...
                   "scan_dates.csv")) ;

output_dir = strcat("/Users/tkmd/MyDrive/Projects/", ...
    "hemisphere_fingerprinting_stroke/data/", ...
    "individual_matrices_andrew/") ;

mkdir(output_dir) ;

for i = 1:size(fcs, 1)

    id = fcs.SubjectCode{i} ;
    patient = fcs.PatientStatus(i) ;
    date = erase(fcs.ScanDate(i), "-") ;

    unpacked = fcs.Data{i} ;
    writematrix(unpacked, ...
        strcat(output_dir, "sub-", id, ...
                "_date-", date, "_patient-", ... 
                num2str(patient), "-r.csv"))

end