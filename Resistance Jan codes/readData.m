function patients = readData(fileName)

    data = readtable(fileName);
    
    NoPatients = size(data,2)/2;
    
    patients = cell(1,NoPatients);
    
    for i = 1:NoPatients
     
       pat.time = data{:,2*(i-1)+1};
       pat.data = data{:,2*(i-1)+2};
       
       indx = pat.data < 0;
       pat.time(indx) = [];
       pat.data(indx) = [];
       
       patients{i} = pat;
    end

end

