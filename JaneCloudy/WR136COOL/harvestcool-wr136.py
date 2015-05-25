import numpy as N

#Tstar_range = [5.0e4, 1.0e5, 1.5e5]
Trange = N.arange(2.0, 9.1, 0.1)
phirange = [9.0, 10.0, 11.0]
hdenrange = [0.0, 1.0, 2.0]

ionized = "photo"
abunid = "ngc6888"

#for Tstar in Tstar_range:
for loghden in hdenrange:
    for logphi in phirange:
        outfile = file("coolfunc-%s-wr136-phi%.2f-%s-n%.2f.dat" % (
                ionized, logphi, abunid, loghden), "w")
        outfile.write('# Temperature\tLambda (erg cm3/s)\tNp\tNe\tL (erg/cm3/s)\tH (erg/cm3/s)\tTop Cooler\tfrac\tTop Heater\tfrac\n')
        for logT in Trange:
            filename = "%s-wr136-phi%.2f-%s-n%.2f-T%.2f" % (
                ionized, logphi, abunid, loghden, logT)
            print "Trying " + filename
            try: 
                f = file("out/%s.cool" % filename, 'r')
                fields = f.readlines()[-1].split('\t')
                T = 10**logT
                Ctot = float(fields[3])
                coolants = '\t'.join(fields[4:6])

                f = file("out/%s.heat" % filename, 'r')
                fields = f.readlines()[-1].split('\t')
                Htot = float(fields[2])
                heatants = '\t'.join(fields[4:6])

                f = file("out/%s.ovr" % filename, 'r')
                fields = f.readlines()[-1].split('\t')
                pden = (10.0**float(fields[3]))*(10.0**float(fields[7]))
                hden = 10.0**float(fields[3])
                eden = 10.0**float(fields[4])
                cool = Ctot/(hden*eden)
            except:
                continue

            outfile.write("%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%s\t%s\n" % (
                    T, cool, pden, eden, Ctot, Htot, coolants, heatants))
