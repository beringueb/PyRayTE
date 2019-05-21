



for sensi in [0,1,2]:
    for fsky in [0.1,0.2,0.4]:
        with open('../ini/SO.ini',"r") as f:
            text = f.read()
        new_text = text.replace('fsky = 0.4','fsky = {:3.1f}'.format(fsky))
        new_new_text = new_text.replace('noise_file = ../data/noise_spectra/SO_2_0.4.dat','noise_file = ../data/noise_spectra/SO_{:d}_{:3.1f}.dat'.format(sensi,fsky))

        with open('../ini/SO/SO_{:d}_{:3.1f}.ini'.format(sensi,fsky),"w") as f:
            f.write(new_text)
        
        
        