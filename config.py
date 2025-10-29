import numpy as np
from geometry import *
from util.h5py_util import *
import yaml,sys

WATERSPEED = 1481.44134805
WATERDENST = 999.6150851557516
REFFREQUENCY = 1000000.0

class structtype():
    pass

def parse_config(filename = 'config.yaml'):

    g = structtype()
    
    with open(filename, 'r') as f:
        
        cfg = yaml.safe_load(f)
        
        g.enhance = 1.1
        g.n = 3
        g.golden_ratio = 1
        g.grid_ratio = 0.25
        g.m = 0
        g.hankel_kind = 2
        
        g.rj = [WATERDENST]     # 1 kg/m^3
        g.cj = [WATERSPEED]                # 1490 m/s
        g.frequency_reference = REFFREQUENCY # 1 MHz frequency
        
        g.fac0 = np.array([1,- 1])
        g.fac = np.array([1,- 1])
        g.amp = np.array([1,1])
        g.mud = np.array([[1,- 1],[1,- 1]])
        g.tra = np.array([1,1])        
        
        g.convergemod = cfg['CONVE_MOD']
        g.model_scale = cfg['MODEL_SCL']
        g.att = cfg['ATTEN_RAT']
        g.tra = cfg['TRANS_RAT']
              
        g.ifu = int(cfg['FREQU_INI'])
        g.ffu = int(cfg['FREQU_FIN'])
        
        g.nff = int(cfg['FREQF_DEL'])
        g.iff = cfg['FREQF_INI']
        g.fff = cfg['FREQF_FIN']

        g.fr  = g.model_scale * g.frequency_reference
        g.sfr = g.model_scale * g.frequency_reference * np.linspace(g.iff,g.fff,g.nff)

        # REFERENCE IMPEDANCE AND WAVELENGTH
        g.zed = np.prod((g.cj,g.rj))
        g.wav = np.divide(g.cj,g.fr)

        g.piston_radius = g.model_scale * cfg['PISTO_RAD']
        g.piston__pitch = g.model_scale * cfg['PISTO_VCT']
        g.piston__catch = g.model_scale * cfg['PISTO_VCR']
        g.piston_distan = g.model_scale * cfg['PISTO_DST']
        g.interf_centre = g.model_scale * cfg['INTER_CEN']
        g.piston_centre = g.model_scale * cfg['PISTO_CEN']
        g.tank_rad = g.model_scale * cfg['TANK_RAD']
        g.scale = 4

        GX = ((cfg['GRIDX_INI']),(cfg['GRIDX_FIN']),int(cfg['GRIDX_DEL']))
        GY = ((cfg['GRIDY_INI']),(cfg['GRIDY_FIN']),int(cfg['GRIDY_DEL']))
        CX = ((cfg['CENTX_INI']),(cfg['CENTX_FIN']),int(cfg['CENTX_DEL']))
        CY = ((cfg['CENTY_INI']),(cfg['CENTY_FIN']),int(cfg['CENTY_DEL']))
        SP = ((cfg['RATSP_INI']),(cfg['RATSP_FIN']),int(cfg['RATSP_DEL']))
        DS = ((cfg['RATDS_INI']),(cfg['RATDS_FIN']),int(cfg['RATDS_DEL']))
        PS = ((cfg['RATPS_INI']),(cfg['RATPS_FIN']),int(cfg['RATPS_DEL'])) 
 
        g.grid_x = g.model_scale * np.linspace(GX[0],GX[1],GX[2])
        g.grid_y = g.model_scale * np.linspace(GY[0],GY[1],GY[2])
        g.cent_x = g.model_scale * np.linspace(CX[0],CX[1],CX[2])
        g.cent_y = g.model_scale * np.linspace(CY[0],CY[1],CY[2])

        g.skr = np.linspace(SP[0],SP[1],SP[2])
        g.sdr = np.linspace(DS[0],DS[1],DS[2])
        g.spr = np.linspace(PS[0],PS[1],PS[2])
        
        g.radia_siz = int(cfg['RADIA_SIZ'])
        g.eleme_wav = int(cfg['ELEME_WAV'])

        g.piston_angle = int(cfg['PISTO_ANG'])
        g.interf_angle = int(cfg['INTER_ANG'])

    return g

def create_configfile(fn,filename, output_path): 
    
    g = fn(filename)
       
    nRD = g.radia_siz
    
    Neltoverlambda = g.eleme_wav
       
    RD = nRD*g.wav
    
    max_k0 = 2 * np.pi* nRD * g.sfr[-1] / g.fr
    
    max_k = np.real(max_k0 * g.skr[-1] / 1000)
    
    max_k0 = np.real(max_k0)
       
    # SET DIMENSIONS FOR RECTANGULAR GRID
    lambda0 = 2 * np.pi/ max_k0
   
    # COORDINATE CENTRE INTERFACE AND PISTON
    height_interf_centre = g.interf_centre
    height_piston_centre = g.piston_centre
    
    width_interf_centre = g.piston_distan / 2
    width_piston_centre = g.piston_distan
    
    interf_centre = np.array((width_interf_centre,height_interf_centre))
    piston_centre = np.array((width_piston_centre,height_piston_centre))
    
    interf_angle = np.array(g.interf_angle)
    piston_angle = np.array(g.piston_angle)
    ##
    # CREATE SURROUNDING SURFACE, TRANSMITER TRANSDUCER SURFACE, RECEPTOR TRANSDUCER (SURFACE PROBE)
    ##

    if g.tank_rad:    
        interf_centre[0] = interf_centre[0] * 2
        S = fn_discretize_geometry_plane(4 * g.tank_rad, interf_centre , - np.pi / 2 ,Neltoverlambda / 100, 0, list(reversed(fn_rotation())))
        T = fn_discretize_geometry_circular(g.piston__pitch,[piston_centre[1], piston_centre[0]], np.pi ,Neltoverlambda / (100), 0,  fn_rotation(), g.tank_rad)
        R = fn_discretize_geometry_circular(g.piston__catch, [piston_centre[1], piston_centre[0]], 0,Neltoverlambda / (100), 0,  fn_rotation(), g.tank_rad)
    else:
        S = fn_discretize_geometry_plane(2 * piston_centre[0], interf_centre, - np.pi / 2 ,Neltoverlambda / 100, 0, list(reversed(fn_rotation())))
        T = fn_discretize_geometry_plane(g.piston__pitch,[0,piston_centre[1]], 0 ,Neltoverlambda / (100), piston_angle )
        R = fn_discretize_geometry_plane(g.piston__catch,piston_centre, np.pi,Neltoverlambda / (100), 0 )
               
    S.c = S.x + 1j * S.z

    ##
    # RADIAL SHIFT BASED ON RECTANGULAR PACKING
    ##

    fn_surface_rectangular_stacking(S)
    ##
    # SCALE SURROUNDING SURFACE WITH SCALING FACTOR
    ##    
    
    S.a = S.a * RD
    S.x = S.x * RD
    S.z = S.z * RD
    S.y = S.y * RD
    S.c = S.c * RD
    S.ci = S.ci * RD
    S.co = S.co * RD

    ##    
    # SLICE TRANSMITER SURFACE WITH SURROUNDING SURFACE
    ##
    
    _,ndx0 = fn_enclosure_rectan([T.x,T.z],[ 1, 1],interf_centre,0)
    _,ndx1 = fn_enclosure_rectan([T.x,T.z],[ 1,-1],interf_centre,0)

    T.ndx  = [ndx0,ndx1]

    ##
    # SCALE TRANSMITTER SURFACE WITH SCALING FACTOR
    ##    
        
    T.a = T.a * RD * (T.ndx[0] + T.ndx[1]*np.sqrt(g.tra))
    T.x = T.x * RD
    T.z = T.z * RD
    T.y = T.y * RD
    T.c = T.x + 1j * T.z  
    
    ##
    # SLICE RECEPTOR SURFACE WITH SURROUNDING SURFACE
    ##
    
    _,ndx0 = fn_enclosure_rectan([R.x, R.z],[ 1, 1],interf_centre,0)
    _,ndx1 = fn_enclosure_rectan([R.x, R.z],[ 1,-1],interf_centre,0)

    R.ndx  = [ndx0,ndx1]
    
    ##
    # SCALE RECEPTOR SURFACE WITH SCALING FACTOR
    ##

    R.a = R.a * RD * (R.ndx[0] + R.ndx[1]*np.sqrt(g.tra))
    R.x = R.x * RD
    R.z = R.z * RD
    R.y = R.y * RD
    R.c = R.x + 1j * R.z

    ##
    # CREATE DOMAIN MESH (FIELD PROBE)
    ##
    
    D = fn_discretize_geometry_domain([g.grid_x,g.grid_y],[g.cent_x,g.cent_y],g.grid_ratio)

    ##
    # REMOVE DOMAIN CLOSE TO SURROUNDING SURFACE
    ##

    _,ndx0 = fn_enclosure_rectan([D.x, D.z],[ 1, 1],interf_centre,.3)
    _,ndx1 = fn_enclosure_rectan([D.x, D.z],[ 1,-1],interf_centre,.3)
    
    D.ndx = np.logical_and(ndx0,ndx1)
    fn_copy_filter(D,D)
    D.ndx0 = np.logical_and(ndx0,ndx1)

    ##
    # SLICE DOMAIN WITH SURROUNDING SURFACE
    ##
    
    _,ndx0 = fn_enclosure_rectan([D.x,D.z],[ 1, 1],interf_centre,0)
    _,ndx1 = fn_enclosure_rectan([D.x,D.z],[ 1, -1],interf_centre,0)
    
    D.ndx = [ndx0,ndx1]
    
    ##
    # SCALE DOMAIN SURFACE WITH SCALING FACTOR
    ##

    D.a = D.a * RD
    D.x = D.x * RD
    D.z = D.z * RD
    D.y = D.z * 0
    D.c = D.x + 1j * D.z
    
    keys = ['T','S','D','R','Neltoverlambda','nRD','g']
    values = [T,S,D,R,Neltoverlambda,nRD,g]
    dict = {key: value for key, value in zip(keys, values)}

    configfile = 'P' + g.convergemod + '_' + str(g.nff) + '_' + str(int(g.iff*2*g.model_scale*100)) + '_' + str(int(g.fff*2*g.model_scale*100))
    save_dict_to_hdf5(dict, output_path, configfile)

    return T,S,D,R,Neltoverlambda,nRD,g

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} config_file.yaml'.format(sys.argv[0]))

    input_file = sys.argv[1]
    output_path = sys.argv[2]
    
    create_configfile(parse_config, input_file, output_path)