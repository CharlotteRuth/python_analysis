import pandas as pd
import numpy as np
import sys
import os 
import tqdm

import pynbody

pynbody.config['halo-class-priority'] =  [pynbody.halo.ahf.AHFCatalogue,
                                          pynbody.halo.GrpCatalogue,
                                          pynbody.halo.AmigaGrpCatalogue,
                                          pynbody.halo.legacy.RockstarIntermediateCatalogue,
                                          pynbody.halo.rockstar.RockstarCatalogue,
                                          pynbody.halo.subfind.SubfindCatalogue, pynbody.halo.hop.HOPCatalogue]


def get_shortname_base(sim):
    if sim=='h148mint':
        shortname='h148.cosmo50PLK.6144g'
        base='h148.cosmo50PLK.6144g3HbwK1BH'
    if sim=='h329mint':
        shortname='h329.cosmo50PLK.6144g'
        base='h329.cosmo50PLK.6144g5HbwK1BH'
    if sim=='h148'
        shortname='h148.cosmo50PLK.3072g'
        base='h148.cosmo50PLK.3072g3HbwK1BH'
    if sim=='h229':
        shortname='h229.cosmo50PLK.3072g'
        base='h229.cosmo50PLK.3072gst5HbwK1BH'
    if sim=='h242':
        shortname='h242.cosmo50PLK.3072g'
        base='h242.cosmo50PLK.3072gst5HbwK1BH'
    if sim=='h329':
        shortname='h329.cosmo50PLK.3072g'
        base='h329.cosmo50PLK.3072gst5HbwK1BH'
    return shortname, base



#new_keys = ['h148-12','h148-27','h148-34','h148-38','h148-55','h148-65','h148-251','h148-249','h148-282','h229-14','h229-18','h229-20','h229-22','h229-49','h242-21','h242-38','h242-69','h329-29','h329-117']
#new_keys = ['h329-7','h329-12','h329-18','h329-41','h329-64']# Mint
#new_keys = ['h148-2','h148-3','h148-4','h148-7','h148-11','h148-16','h148-17','h148-19','h148-26','h148-30','h148-32','h148-34','h148-39','h148-49','h148-64','h148-72','h148-133','h148-170','h148-181','h148-195','h148-262','h148-348'] #Mint
new_keys =  ['h148-2','h148-3','h148-4','h148-7','h148-11','h148-16','h148-17','h148-19','h148-26','h148-30','h148-32','h148-34','h148-39','h148-49','h148-64','h148-72','h148-133','h148-170','h148-181','h148-195','h148-262','h148-348', 'h329-7','h329-12','h329-18','h329-41','h329-64']

#new_keys = np.append(new_keys, ['h148-1','h229-1','h242-1','h329-1'])
#new_keys = np.append(new_keys, ['h329-1'])
#new_keys = np.append(new_keys, ['h148-1'])
new_keys = np.append(new_keys, ['h148-1','h329-1'])



new_keys = ['h148-12','h148-27','h148-34','h148-38','h148-55','h148-65','h148-251','h148-249','h148-282',
    'h229-14','h229-18','h229-20','h229-22','h229-49','h242-21','h242-38','h242-69','h329-29','h329-117']

new_keys = np.append(new_keys, ['h148-1','h229-1','h242-1','h329-1'])

haloids_h148 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h148')], 'h148-', '').astype(int)
haloids_h148 = dict(zip(haloids_h148, [[]]*len(haloids_h148)))

haloids_h329 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h329')], 'h329-', '').astype(int)
haloids_h329 = dict(zip(haloids_h329, [[]]*len(haloids_h329)))

haloids_h229 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h229')], 'h229-', '').astype(int)
haloids_h229 = dict(zip(haloids_h229, [[]]*len(haloids_h229)))

haloids_h242 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h242')], 'h242-', '').astype(int)
haloids_h242 = dict(zip(haloids_h242, [[]]*len(haloids_h242)))

haloids_h329 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h329')], 'h329-', '').astype(int)
haloids_h329 = dict(zip(haloids_h329, [[]]*len(haloids_h329)))

for key in tqdm.tqdm(new_keys):
    sim = key[:4]
    z0haloid = int(key[5:])

    shortname, base = get_shortname_base(sim)


    os.environ['TANGOS_DB_CONNECTION'] = '/data2/REPOSITORY/tangos.db'
    os.environ['TANGOS_SIMULATION_FOLDER'] = f'/home/christenc/REPOSITORY/e12Gals/{base}'

    os.environ['TANGOS_DB_CONNECTION'] = '/home/akinshol/Tangos-2022/databases/JL.db'
    os.environ['TANGOS_SIMULATION_FOLDER'] = f'/home/christenc/Data/Sims/{shortname}/{base}/'
    import tangos

    print(f'Computing merger trees for {sim}-{z0haloid}')

    # Using the tangos merger tree
    #tangos_halo = tangos.get_halo(f'snapshots_200crit_{sim}mint/{base}.004096/halo_{z0haloid}')
    #time, haloid = tangos_halo.calculate_for_progenitors("t()","halo_number()")

    # Using the halo_trace merger tree
    #print(haloid)
    #print(f'/home/christenc/REPOSITORY/e12Gals/{base}/{base}.004096.M200.trace_back.hdf5')
    
    #trace = pd.read_hdf(f'/home/christenc/REPOSITORY/e12Gals/{base}/{base}.004096.M200.trace_back.hdf5')
    #haloid = [z0haloid] + trace.loc[z0haloid].tolist()
    tangos_halo = tangos.get_halo(f'snapshots_200crit_{sim}/{base}.004096/halo_{z0haloid}')
    time, haloid = tangos_halo.calculate_for_progenitors("t()","halo_number()")

import glob
#for sim in ['h148mint', 'h329mint']: # mint
for sim in ['h148','h229','h242','h329']: 
    shortname, base = get_shortname_base(sim)
    filepaths = np.array(glob.glob(f'/home/christenc/REPOSITORY/e12Gals/{base}/snapshots_200crit_{sim}mint/{base}.00????')) # mint
    filepaths = np.array(glob.glob(f'/home/christenc/Data/Sims/{shortname}/{base}/snapshots_200crit_{sim}/{base}.00????'))
    s = np.argsort([-int(f[-4:]) for f in filepaths])
    filepaths = filepaths[s]
    exec(f'filepaths_{sim} = filepaths')




output = dict(
    filepaths = dict(
        h148 = filepaths_h148, 
        h229 = filepaths_h229, 
        h242 = filepaths_h242, 
        h329 = filepaths_h329
    ),
    haloids = dict(
        h148 = haloids_h148, 
        h229 = haloids_h229, 
        h242 = haloids_h242,
        h329 = haloids_h329
    )
)

print(output)

import pickle
# with open('/home/christenc/Code/Datafiles/filepaths_haloids_2022.pickle','wb') as f:
with open('filepaths_haloids_2022.pickle','wb') as f:
    pickle.dump(output, f)
    



#assert 1==2
#tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + base + f'.004096/halo_{haloid}')
#time, a, id_main, mvir_main, rvir_main = tangos_halo.calculate_for_progenitors("t()", "a()", "halo_number()",
#"Mvir", "Rvir")




