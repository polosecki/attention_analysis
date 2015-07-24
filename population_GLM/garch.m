monkeys={'Quincy','Michel'};
areas={'PITd','LIP'};
high_res_output=0;
noise_model='poisson';
for monkey=1:2
    for area=1:2
        save_many_GLMs(monkeys{monkey},areas{area},high_res_output,noise_model)
    end
end