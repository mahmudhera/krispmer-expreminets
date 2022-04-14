num_targets = 105

for i in range(1, num_targets+1):
    i_str = str(i)
    f = open( "target" + i_str + "_guidescan_out", 'r' )
    lines = f.readlines()[8:]
    for line in lines:
        mer23 = line.split()[9]
        if 'NGG' in mer23:
            print(mer23[:-3])
        else:
            print(mer23[3:])
    print('done for ' + i_str)
    break
