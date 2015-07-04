function bad_files=bad_files(monkey,area)

if strcmp(monkey,'Quincy') && strcmp(area,'PITd')
    bad_files=[4 8 22 23 24 32 33 42:46 52];
elseif strcmp(monkey,'Quincy') && strcmp(area,'LIP')
    bad_files=[2 12:14 17 26];
elseif strcmp(monkey,'Michel') && strcmp(area,'LIP')
    bad_files=[1 2 4 5 21 31 32 33];
elseif strcmp(monkey,'Michel') && strcmp(area,'PITd')
    bad_files=[6 10:12 16 19 32 38 40 48];
else
    error('You suck, dude')
end

