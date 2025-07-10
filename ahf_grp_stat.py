import numpy as np
import pynbody
from pathlib import Path
import sys

def main(filename):
    f = pynbody.load(filename)
    h = f.halos(halo_numbers='v1')
    if not Path(filename+'.amiga.grp').is_file():
        f['amiga.grp'] = h.get_group_array()
        f['amiga.grp'].write()
        
    my_properties = h.get_properties_all_halos()
    halo_ids = np.arange(1,len(my_properties['ID']) + 1)

    # Find all halos with less than contam_threshold mass in low-res DM particles
    contam_threshold = 0.1
    contam_cut = my_properties['fMhires'] > 1 - contam_threshold

    header = ('Grp', 'N_tot', 'N_gas', 'N_star', 'N_dark', 'Mvir(M_sol)', 'Rvir(kpc)', 'GasMass(M_sol)', 'StarMass(M_sol)', 'DarkMass(M_sol)', 'V_max', 'R@V_max', 'VelDisp', 'Xc', 'Yc', 'Zc', 'VXc', 'VYc', 'VZc', 'Contam', 'Satellite?') #, 'N_BH', 'False?')
    grp_arr = (halo_ids[contam_cut]).astype(int)
    n_tot_arr = (my_properties['npart'][contam_cut]).astype(int)
    n_gas_arr = (my_properties['n_gas'][contam_cut]).astype(int)
    n_star_arr = (my_properties['n_star'][contam_cut]).astype(int)
    n_dark_arr = n_tot_arr-n_star_arr-n_gas_arr
    mvir_arr = my_properties['Mvir'][contam_cut]/f.properties['h']
    rvir_arr = my_properties['Rvir'][contam_cut]/f.properties['h']
    gmass_arr = my_properties['M_gas'][contam_cut]/f.properties['h']
    smass_arr = my_properties['M_star'][contam_cut]/f.properties['h']
    dmass_arr = rvir_arr-gmass_arr-smass_arr
    vmax_arr = my_properties['Vmax'][contam_cut]
    rvmax_arr = my_properties['Rmax'][contam_cut]
    vdisp_arr = my_properties['sigV'][contam_cut]
    xc_arr = my_properties['Xc'][contam_cut]/f.properties['h']
    yc_arr = my_properties['Yc'][contam_cut]/f.properties['h']
    zc_arr = my_properties['Zc'][contam_cut]/f.properties['h']
    vxc_arr = my_properties['VXc'][contam_cut]
    vyc_arr = my_properties['VYc'][contam_cut]
    vzc_arr = my_properties['VZc'][contam_cut]
    contam_arr = my_properties['fMhires'][contam_cut]<1
    sat_arr = my_properties['hostHalo'][contam_cut]

    # Stack all arrays column-wise
    data = np.column_stack([
        grp_arr,
        n_tot_arr,
        n_gas_arr,
        n_star_arr,
        n_dark_arr,
        mvir_arr,
        rvir_arr,
        gmass_arr,
        smass_arr,
        dmass_arr,
        vmax_arr,
        rvmax_arr,
        vdisp_arr,
        xc_arr,
        yc_arr,
        zc_arr,
        vxc_arr,
        vyc_arr,
        vzc_arr,
        contam_arr,
        sat_arr
    ])
    
    # Convert header tuple to a string (space-separated)
    header_str = ' '.join(header)
    
    # Save to file
    np.savetxt(
        filename + '.amiga.stat',   # change to your desired filename
        data,
        delimiter='\t',
        header=header_str,
        fmt='%s'  # use '%s' to handle mixed data types (e.g. integers, floats, booleans)
    )


    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print ('Usage: python ahf_grp_stat.py <sim>')
        sys.exit()
    else:
        filename = str(sys.argv[1])

    main(filename)

    
