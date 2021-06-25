#Calculates the distance to the nearest massive galaxy

def distance_to_nearest_host(data):
    distances = []
    hostrvirs = []
    for i in range(len(data)):
        s = data['sim'].tolist()[i]
       
        if s=='h148' or s=='h229' or s=='h242' or s=='h329': # if sat simulation, find distance to halo 1
            h1dist = data['h1dist'].tolist()[i]*0.6776942783267969
            distances.append(h1dist)
       
            h1rvir = data['Rvir'][(data.sim==s) & (data.haloid==1)].tolist()[0]*0.6776942783267969
            hostrvirs.append(h1rvir)
           
        else: # if field simulation, find distance to nearest massive DM halo (currently > 0.5e12 Msol)
            if s=='cptmarvel':
                path = '/home/akinshol/Data/Sims/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

            if s=='elektra':
                path = '/home/akinshol/Data/Sims/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096.dir/elektra.cosmo25cmb.4096g5HbwK1BH.004096'

            if s=='rogue':
                path = '/home/akinshol/Data/Sims/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096.dir/rogue.cosmo25cmb.4096g5HbwK1BH.004096'

            if s=='storm':
                path = '/home/akinshol/Data/Sims/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
           
            coords = []
            with open(path+'.coords','rb') as f:
                while True:
                    try:
                        coords.append(pickle.load(f,encoding='latin1'))
                    except EOFError:
                        break
            coords = pd.DataFrame(coords)
           
            threshold = 5*10**(11) # this threshold can be adjusted,
            # i tried to pick something similar to the virial masses of the host in the JL simulations
           
            coords = coords[coords.mass > threshold]
           
            halocoords = np.array([data['Xc'].tolist()[i],data['Yc'].tolist()[i],data['Zc'].tolist()[i]])

            x = np.array(coords['Xc'])
            y = np.array(coords['Yc'])
            z = np.array(coords['Zc'])
            Rvir = np.array(coords['Rv'])

           
            c = np.array([x,y,z])
            c = np.transpose(c)
            dist = np.sqrt(np.sum((halocoords-c)**2, axis=1))*0.6776942783267969
            distances.append(np.min(dist))
            hostrvirs.append(Rvir[np.argmin(dist)]*0.6776942783267969)
           
    return np.array(distances),np.array(hostrvirs)
