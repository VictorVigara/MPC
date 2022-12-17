function cell_data = cell_data2Plot(T,X,sim_name,cell_data)
    idx_mass = size(cell_data,1) + 1;
    cell_data{idx_mass,1} = [T];
    cell_data{idx_mass,2} = [X];
    cell_data{idx_mass,3} = sim_name;
end