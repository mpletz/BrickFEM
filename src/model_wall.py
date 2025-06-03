# functions for the optimizing a Lego wall design

import numpy as np
import matplotlib.pyplot as plt
import subprocess,os,shutil

#from itertools import *
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from PIL import Image

DIR0 = os.path.abspath('')

def make_dir(dir_name, if_change=0, if_clear=0):
    """Create folder and possibly change the work directory to it

    Args:
    dir_name (str): String that defines the directory name
    if_change (bool): If the directory should be set as the work directory
    if_clear (bool): If the directory already exists, should its content (including all sub-folders) be deleted 

    Returns:
    dir_abs (str): The absolute path of the created directory
    """
    dir_abs = os.path.abspath('')

    # if path does not exist: create
    if os.path.exists(dir_name) == 0:
        os.mkdir(dir_name)
    else:
        # if it exists: clear if if_clear==1
        if if_clear:
            shutil.rmtree(dir_name)
            os.mkdir(dir_name)
    dir1 = dir_abs + "/" + dir_name

    # change into the dir_name directory
    if if_change:
        os.chdir(dir1)
    return dir_abs

# Lego colors from bricklink (only solid colors): https://www.bricklink.com/catalogColors.asp
cm_lego = {
'White': 'ffffff', 'Very Light Gray': 'e5e5e5', 'Very Light Bluish Gray': 'e0e5e5',
'Light Bluish Gray': 'a5acbf', 'Light Gray': '919191', 'Dark Gray': '615050',
'Dark Bluish Gray': '4e5355', 'Black': '1e1e1e', 'Dark Red': '631116',
'Red': 'b10b0f', 'Coral': 'ff7669', 'Dark Salmon': 'ff592a', 'Salmon': 'ff7256',
'Light Salmon': 'ffc0af', 'Sand Red': 'c28276', 'Dark Brown': '48302a', 'Brown': '633721',
'Light Brown': '925b39', 'Medium Brown': '9b613d', 'Reddish Brown': '7b3a27',
'Fabuland Brown': 'af5e47', 'Dark Tan': 'b28d61', 'Medium Tan': 'd6bd8c', 'Tan': 'eed39e',
'Light Nougat': 'ffc5a9', 'Nougat': 'ffa576', 'Medium Nougat': 'e49556',
'Dark Nougat': 'cd6e3f', 'Fabuland Orange': 'f4852c', 'Earth Orange': 'e97c28',
'Dark Orange': 'b04a17', 'Rust': 'af401c', 'Neon Orange': 'ff5041', 'Orange': 'ff7324',
'Medium Orange': 'ff9a38', 'Bright Light Orange': 'ffbe2c', 'Light Orange': 'ffb23f',
'Very Light Orange': 'ffd69e', 'Dark Yellow': 'de8c34', 'Yellow': 'ffda32',
'Light Yellow': 'ffe39a', 'Bright Light Yellow': 'ffec89', 'Neon Yellow': 'fff938',
'Neon Green': 'd5ef5b', 'Light Lime': 'e9eab8', 'Yellowish Green': 'e3efa2',
'Medium Lime': 'dcd931', 'Fabuland Lime': 'a1ca41', 'Lime': 'bbd930', 'Dark Olive Green': '6b693b',
'Olive Green': 'a29e50', 'Dark Green': '254b3c', 'Green': '00863c', 'Bright Green': '00c23b',
'Medium Green': '7ed986', 'Light Green': 'cfebcc', 'Sand Green': '95b69a',
'Dark Turquoise': '009794', 'Light Turquoise': '00bdb4', 'Aqua': 'afe1d7', 'Light Aqua': 'c5ece7',
'Dark Blue': '1e314c', 'Blue': '004f99', 'Dark Azure': '0096d8', 'Little Robots Blue': '43b7de',
'Maersk Blue': '69b9d1', 'Medium Azure': '50c7da', 'Sky Blue': '76cedc', 'Medium Blue': '71a4d0',
'Bright Light Blue': 'b1cbe9', 'Light Blue': 'bfd4dc', 'Sand Blue': '7b8fa0',
'Dark Blue-Violet': '1731a2', 'Violet': '2a4296', 'Blue-Violet': '4065e7', 'Lilac': '6d5bc3',
'Medium Violet': '8788dc', 'Light Lilac': 'c6c6ea', 'Light Violet': 'c1c3dd',
'Dark Purple': '572675', 'Purple': '72257f', 'Light Purple': 'ab3087', 'Medium Lavender': 'c280d0',
'Lavender': 'ceb6dd', 'Clikits Lavender': 'dfa2d2', 'Sand Purple': 'b07399', 'Magenta': 'b52469',
'Dark Pink': 'f554a7', 'Medium Dark Pink': 'fd7ca6', 'Bright Pink': 'fab5d3', 'Pink': 'f7c7d0',
'Light Pink': 'f3cdcb'}

col_dict = {1:'Light Bluish Gray',2:'Bright Green',3:'Yellow',4:'Medium Blue',5:'Light Lime',
            6:'Medium Orange',7:'Magenta',8:'Medium Lavender'}

col_set = {i:'#'+cm_lego[j] for i,j in col_dict.items()}

# brick functions
# ---------------------------------------------------------------------------

def make_full_from_sym(wall_des):
    # get full array of 0,1 entries from the right side of it
    return np.append(np.array(wall_des)[:,:0:-1],np.array(wall_des),axis=1)

def make_bricks(wall):
    def get_bricks(row):
        # function for getting brick lengths from 0,1 list (one row)
        brick_list = []
        i_temp = 1
        for is_con in row:
            if is_con == 1:
                i_temp += 1
            else:
                brick_list += [i_temp]
                i_temp = 1
        brick_list += [i_temp]
        return brick_list

    wall_bricks = []
    for row in wall:
        wall_bricks += [get_bricks(row)]
    return wall_bricks

def plot_wall(wall_bricks,out_name='',if_show=1):
    b0, h0 = 8.,9.6

    fig, ax = plt.subplots(figsize=(4,4))

    for i_row,row in enumerate(wall_bricks[::-1]):
        x0 = 0
        for brick_i in row:
            ax.add_patch(plt.Rectangle((x0,i_row*h0),brick_i*b0,h0,
                        facecolor=col_set[brick_i],edgecolor='k'))
            x0 += brick_i*b0

    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()
    if out_name != '':
        plt.savefig(out_name+'.png',dpi=600)
    if if_show:
        plt.show()
    plt.close(fig)

# model functions
# ---------------------------------------------------------------------------

def create_ass_dict_plates(run_dir,model_name,wall_bricks,nx=6,ny=5):
    # make dict and wall Python file
    #
    t_step = 0.0006

    ass_str = """# automatically generated brickFEM wall model

import sys
sys.path.append("..")

from brickfem import make_model

assembly = {'name':'%s',
'bricks':{0:{'type':'base-plate', 'nx':%i, 'nz':1},
1:{'type':'plate', 'nx':1, 'nz':1},
2:{'type':'plate', 'nx':2, 'nz':1},
3:{'type':'plate', 'nx':3, 'nz':1},
4:{'type':'plate', 'nx':4, 'nz':1},
5:{'type':'plate', 'nx':5, 'nz':1},
6:{'type':'plate', 'nx':6, 'nz':1},
7:{'type':'plate', 'nx':7, 'nz':1},
8:{'type':'plate', 'nx':8, 'nz':1}},
    """%(model_name,nx)

    ass_end_str = """
'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
'loads_rp':{},
'mesh_size':0.75, 'mu':0.2
}

e_par_wall = {'t_step': %f, 'mass_scale_t': 1.2e-7, 'is_acc': 0,
                'load_str': '',
                'loads_rigid': {1: {'shape': 'sphere', 'radius': 5., 'loc': (%f*8, 3.2*%f, 4+5),
                                    'u': (0, 0, -12/0.0006*%f)}}}

make_model(assembly, e_par_wall)
    """%(t_step, (nx-1)/2., ny*0.5, t_step)

    b0, h0 = 8.,3.2 #9.67

    part_str = """'parts': {
    1:{'brick_id':0, 'loc':(0,0,0), 'c':'Brown'},
    """

    i_part = 2

    for i_row,row in enumerate(wall_bricks[::-1]):
        x0 = 0
        y0 = i_row*h0
        for brick_i in row:
            part_str += "%i:{'brick_id':%i,'loc':(%.2f,%.2f,0),'c':'%s'},\n  "%(i_part,
                                    brick_i,x0,y0,col_dict[brick_i])
            i_part += 1
            x0 += brick_i*b0
    part_str = part_str[:-4]+'\n},'

    with open(run_dir+'/temp_wall.py','w') as f:
        f.write(ass_str+part_str+ass_end_str)
    # egal, wird ja eh geloescht
    # with open(model_name+'-expl-mesh075mm/temp_wall.py','w') as f:
    #    f.write(ass_str+part_str+ass_end_str)    
    return ass_str+part_str+ass_end_str


def old_create_ass_dict(run_dir,model_name,wall_bricks,nx=6,ny=5,mesh_size=0.75):
    # make dict and wall Python file
    #
    if 'msize' in model_name:
        t_step = 0.0002 #
    else:
        t_step = 0.0006 #

    ass_str = """# automatically generated brickFEM wall model

import sys
sys.path.append("..")

from brickfem import make_model

assembly = {'name':'%s',
'bricks':{0:{'type':'base-plate', 'nx':%i, 'nz':1},
1:{'type':'regular', 'nx':1, 'nz':1},
2:{'type':'regular', 'nx':2, 'nz':1},
3:{'type':'regular', 'nx':3, 'nz':1},
4:{'type':'regular', 'nx':4, 'nz':1},
5:{'type':'regular', 'nx':5, 'nz':1},
6:{'type':'regular', 'nx':6, 'nz':1},
7:{'type':'regular', 'nx':7, 'nz':1},
8:{'type':'regular', 'nx':8, 'nz':1}},
    """%(model_name,nx)

    ass_end_str = """
'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
'loads_rp':{},
'mesh_size':%f, 'mu':0.2
}

e_par_wall = {'t_step': %f, 'mass_scale_t': 1.2e-7, 'is_acc': 0,
                'load_str': '',
                'loads_rigid': {1: {'shape': 'sphere', 'radius': 5., 'loc': (%f*8, 9.6*%f, 4+5),
                                    'u': (0, 0, -12/0.0006*%f)}}}

make_model(assembly, e_par_wall)
    """%(mesh_size,t_step, (nx-1)/2., ny*0.5, t_step)

    b0, h0 = 8.,9.6

    part_str = """'parts': {
    1:{'brick_id':0, 'loc':(0,0,0), 'c':'Brown'},
    """

    i_part = 2

    for i_row,row in enumerate(wall_bricks[::-1]):
        x0 = 0
        y0 = i_row*h0
        for brick_i in row:
            part_str += "%i:{'brick_id':%i,'loc':(%.2f,%.2f,0),'c':'%s'},\n  "%(i_part,
                                    brick_i,x0,y0,col_dict[brick_i])
            i_part += 1
            x0 += brick_i*b0
    part_str = part_str[:-4]+'\n},'

    with open(run_dir+'/temp_wall.py','w') as f:
        f.write(ass_str+part_str+ass_end_str)
    # egal, wird ja eh geloescht
    # with open(model_name+'-expl-mesh075mm/temp_wall.py','w') as f:
    #    f.write(ass_str+part_str+ass_end_str)    
    return ass_str+part_str+ass_end_str

def create_ass_dict(model_name,wall_bricks,mesh_size,t_step,v,m=1.5,is_new=1,nx=6,ny=5):
    # make dict and wall Python file (m: mass in g)

    ass_str = """# automatically generated brickFEM wall model

import sys
sys.path.append("..")

from brickfem import make_model

assembly = {'name':'%s',
'bricks':{0:{'type':'base-plate', 'nx':%i, 'nz':1},
1:{'type':'regular', 'nx':1, 'nz':1},
2:{'type':'regular', 'nx':2, 'nz':1},
3:{'type':'regular', 'nx':3, 'nz':1},
4:{'type':'regular', 'nx':4, 'nz':1},
5:{'type':'regular', 'nx':5, 'nz':1},
6:{'type':'regular', 'nx':6, 'nz':1},
7:{'type':'regular', 'nx':7, 'nz':1},
8:{'type':'regular', 'nx':8, 'nz':1}},
    """%(model_name,nx)

    if m != 0:
        ass_end_str = """
'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
'loads_rp':{},
'mesh_size':%f, 'mu':0.2
}

e_par_wall = {'t_step': %f, 'mass_scale_t': 1.2e-7, 'is_acc': 0,
                'load_str': '%s',
                'loads_rigid': {1: {'shape': 'sphere', 'radius': 5, 'loc': (%f*8, 9.6*%f, 4+5), 'm':%f,
                                    'v0': (0, 0, -%f)}}}

make_model(assembly, e_par_wall, is_new=%i)
        """%(mesh_size, t_step, 'mass'+str(int(m*10)), (nx-1)/2., ny*0.5,m*1e-6, v*1000, is_new)
    else:
            ass_end_str = """
'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
'loads_rp':{},
'mesh_size':%f, 'mu':0.2
}

e_par_wall = {'t_step': %f, 'mass_scale_t': 1.2e-7, 'is_acc': 0,
                'load_str': '',
                'loads_rigid': {1: {'shape': 'sphere', 'radius': 5, 'loc': (%f*8, 9.6*%f, 4+5),
                                    'u': (0, 0, %f)}}}

make_model(assembly, e_par_wall, is_new=%i)
        """%(mesh_size, t_step, (nx-1)/2., ny*0.5, -t_step*v*1000, is_new)
    
    b0, h0 = 8., 9.6

    part_str = """'parts': {
    1:{'brick_id':0, 'loc':(0,0,0), 'c':'Brown'},
    """

    i_part = 2

    for i_row,row in enumerate(wall_bricks[::-1]):
        x0 = 0
        y0 = i_row*h0
        for brick_i in row:
            part_str += "%i:{'brick_id':%i,'loc':(%.2f,%.2f,0),'c':'%s'},\n  "%(i_part,
                                    brick_i,x0,y0,col_dict[brick_i])
            i_part += 1
            x0 += brick_i*b0
    part_str = part_str[:-4]+'\n},'

    with open('temp_wall.py','w') as f:
        f.write(ass_str+part_str+ass_end_str)

    return ass_str+part_str+ass_end_str


def get_fu_curve0(res_dir):
    #
    loc_files = os.listdir(res_dir)
    # name
    res_name = res_dir.split('/')[-1][:-15]
    #
    ho_file = [i for i in loc_files if i.endswith('-ass.dat')][0]
    rp_file = [i for i in loc_files if i.endswith('-node00.dat')][0]
    #
    # RF3, RF2, RF1, U1, U3, U2, V1, V2, V3, time
    rp_res = np.loadtxt(res_dir+'/'+rp_file,delimiter=',',skiprows=1)

    # time, ALLSE, ALLKE, ALLWK, ALLFD
    ho_res = np.loadtxt(res_dir+'/'+ho_file,delimiter=',',skiprows=1)

    # name, ho_res, t_rp, z_rp, vz_rp, rfz_rp
    return res_name, ho_res, rp_res[:,-1], rp_res[:,4], rp_res[:,-2], rp_res[:,0]

def get_energies(ho_res,t,z,rfz):
    # get the total energies of one Lego wall model
    
    # compute dissipated energy of sphere (J)
    #e_diss_sphere = np.trapz(rfz,z)/1000

    fz_max = max(rfz)
    se_max = max(ho_res[:,1])

    # obtain maximum frictional dissipation (J)
    e_fric = ho_res[-1,-1]/1000

    # obtain kinetic energy at the end (J)
    e_kin = ho_res[-1,-3]/1000

    # obtain work at the end (J)
    e_wk = ho_res[-1,-2]/1000

    return e_wk, e_kin, e_fric

def make_energy_plot(res_dir,wall_bricks,m=0):

    name, ho_data, t_rp, z_rp, vz_rp, rfz_rp = get_fu_curve0(res_dir)

    e_wk, e_kin, e_fric = get_energies(ho_data,t_rp,z_rp,rfz_rp)

    # time, ALLSE, ALLKE, ALLWK, ALLFD
    #fig, (ax0,ax) = plt.subplots(figsize=(5,3.5),nrows=2,sharex=True)
    fig, (ax0,ax) = plt.subplots(figsize=(5,4.5),nrows=2,sharex=True)

    ax.plot(ho_data[:,0]*1000,ho_data[:,1]/1000,label='ALLSE')
    ax.plot(ho_data[:,0]*1000,ho_data[:,2]/1000,label='ALLKE')
    ax.plot(ho_data[:,0]*1000,ho_data[:,3]/1000,label='ALLWK')
    ax.plot(ho_data[:,0]*1000,ho_data[:,4]/1000,label='ALLFD')

    if m == 0.:
        ax0.plot(t_rp*1000,-rfz_rp/1000,label='RFZ (kN)')
        ax0.set_ylim(ymin=0)
    else:
        ax0.plot(t_rp*1000,-vz_rp/1000,label=r'$v_\mathrm{z}$ (m/s)')
        ax0.set_ylim(-21,21)

    ax0.set_title(name+r' ($E_\mathrm{fric}$ = '+str(round(e_fric,4))+' Nm)')

    ax0.set_xlim(0,max(t_rp*1000))
    ax.set_ylim(ymin=-0.01)

    # Define the position for the image axes
    #axins = fig.add_axes([0.55,0.4,0.47, 0.45])
    axins = fig.add_axes([0.7,0.62,0.47*0.5, 0.45*0.5])

    b0, h0 = 8.,9.6

    for i_row,row in enumerate(wall_bricks[::-1]):
        x0 = 0
        for brick_i in row:
            axins.add_patch(plt.Rectangle((x0,i_row*h0),brick_i*b0,h0,
                            facecolor=col_set[brick_i],edgecolor='k',alpha=0.7))
            x0 += brick_i*b0

    axins.axis('equal')
    axins.axis('off')
    
    ax.set_xlabel('time $t$ (ms)')
    if m == 0.:
        ax0.set_ylabel('reaction force (kN)')
    else:
        ax0.set_ylabel('sphere velocity (m/s)')
    
    ax.set_ylabel('energies (J)')

    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(res_dir+'/'+name+'_res.png',dpi=600)
    plt.close(fig)
    #plt.show()
    return e_wk, e_kin, e_fric

def make_wall_model(model_name,con_arr=(),if_run=1,m=0,is_new=1,t_step=0.0006,v=20,mesh_size=0.75):
    """Run the BrickFEM model for a Lego wall. 

    Args:
        model_name (str): Name for the wall model
        con_arr (tuple, optional): Definition of the bricks in the Lego wall with width of nx(studs) and height of ny(studs). The size of con_arr  ((nx-1) x ny list): 1 if bricks are together at that position and 0 for separate bricks. For example, the bottom row can be (0,0,1,1,1), which corresponds to these bricks 1x1, 1x1, 1x4. If con_arr is (), a random design is created. Defaults to ().
        if_run (int, optional): If the model should be run by starting Abaqus from the command line. Defaults to 1.
        m (float): mass of the sphere. If 0, it moves with a constant speed v


    Returns:
        float: Total dissipated frictional energy in the model
    """
    model_dir = model_name+'-expl-mesh'+str(int(mesh_size*100)).zfill(3)+'mm'
    make_dir(model_dir, if_clear=1)
    print(model_dir)
    
    # create brick array, write image
    ny, nx = np.shape(con_arr)
    nx += 1
    print(nx,ny)

    # create wall bricks from bool array
    wall_bricks = make_bricks(con_arr)

    # create the Abaqus/Python file
    create_ass_dict(model_name,wall_bricks,mesh_size,t_step,v,m,is_new=is_new,nx=nx,ny=ny)

    # run the model as well
    if if_run:
        subprocess.call('abaqus cae nogui=temp_wall.py', shell=True)
        np.savetxt(model_dir+'/wall-arr.dat',con_arr)
        plot_wall(wall_bricks,model_dir+'/wall_design',0)
    
    # TODO: test -- load the results, particularly the energies 
    e_wk, e_kin, e_fric = make_energy_plot(model_dir,wall_bricks,m)
    print('...')
    return e_wk, e_kin, e_fric


# run this file using cmd: python model_wall.py
if __name__ == "__main__":
   
    wall_alternate1x4 = np.array([[0, 1, 1, 1], [1, 1, 0, 1], [0, 1, 1, 1], [1, 1, 0, 1],
                                  [0, 1, 1, 1], [1, 1, 0, 1],[0, 1, 1, 1]])
    wall_alternate1x4 = make_full_from_sym(wall_alternate1x4)
    
    wall_opt = np.array([[1, 0, 1, 0], [1, 0, 0, 1], [1, 1, 0, 1], [0, 0, 1, 0],
                         [1, 1, 0, 1], [0, 0, 1, 0],[0, 0, 0, 0]])
    wall_opt = make_full_from_sym(wall_opt)

    wall_singles = np.zeros((7,7))
    wall_full = np.ones((7,7))
    
    make_wall_model('wall-singles', wall_singles,if_run=1,is_new=1)
    make_wall_model('wall-full', wall_full,if_run=1,is_new=1)
    make_wall_model('wall-altern1x4', wall_alternate1x4,if_run=1,is_new=1)
    
    #make_wall_model('wall-optimized', wall_opt,if_run=1,is_new=1)

    # simple setup for mesh study (3x4 wall)
    #make_wall_model('wall-m0', con_arr=np.array(((1,1),(1,0),(1,1),(0,1))),  
    #                if_run=1,mesh_size=0.75,m=0,t_step=0.0004)
    