clc;clear all; close all
monkey={'Quincy','Michel'};
area={'PITd','LIP'};
load full_stats

for mm=1:length(monkey)
    for aa=1:length(area)
        [RT_data]= get_RT_from_cell_list(monkey{mm},area{aa});
        full_stats(aa,mm).RT_data=RT_data;        
    end
end

save full_stats full_stats