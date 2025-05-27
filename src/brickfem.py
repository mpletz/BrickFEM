"""Creating a simple Lego(r) model in Abaqus
# -----------------------------------------------------------------
# Martin Pletz, 2022-11

Use this lego model in a separate Abaqus-Python file that loads the function make_model using `from brickfem import make_model`. Then call the abaqus model similar as in the test cases provided in `model_case_1.py`, etc. If you want to run the model, write `abaqus cae nogui=model_case_1.py` in the command line. Note that the calculations usually crash if you run them interactively in the opened Abaqus CAE or in the command line using `abaqus cae script=lego_model.py`.
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
import numpy as np
import os, random, json, copy, shutil, subprocess

TOL = 1e-3
DIR0 = os.path.abspath('')

# default general Lego brick geometry parameters

lego_geom = {'b, b_gap, b_wall': (8., 0.1, 1.6),
                'h, h_stud, h_top': (9.6, 1.7, 1.5),
                'inside small': {'r': 1.6, 't_rib':1.1, 'h_rib':8.1},
                'inside big': {'r': 3.3, 't': 0.9, 't_rib':0.9, 'h_rib': 6.8},
                'delta_r': 0.05,
                'E, nu, dens': (2200., 0.35, 1e-9)}

# General functions -----------------------------------------------------------------

def remove_files(dir0,type_list=('.log', '.msg')):
    """Remove all files in the directory dir that end with the strings
    given in type_list and additional predefined strings. Some only with
    'ms.mdl' because they are needed for the first, implicit model to
    load the results"""
    type0_list = ('.com', '.sim', '.SMABulk', '.pac', '.abq',
                  '.dmp', '.exception', '.simdir', 'ms.mdl', 'ms.prt',
                  'ms.res', 'ms.res', 'ms.sel', 'ms.stt')
    
    # all files in the folder
    file_list = [i for i in os.listdir(dir0) if '.' in i]

    # number-Types
    number_list = tuple('.'+str(i) for i in range(1,21))

    # select files
    for file_name in file_list:
        for end_str in (type0_list + type_list + number_list):
            if file_name.endswith(end_str):
                try:
                    os.remove(file_name)
                except:
                    None
    return


def run_model(model, job_name, n_proc=4):
    """Run Abaqus model with a job called job_name using n_proc processors. Use double precision.
    """
    # if job_name is empty, just use the model name
    if job_name == '':
        job_name = model.name
    
    # Save the cae, create and run job. Watch out: job.waitForCompletion() usually crashes if the script is not run using `abaqus cae nogui=...` in the command line
    mdb.saveAs(pathName=model.name + '.cae')
    job = mdb.Job(model=model.name, name=job_name, type=ANALYSIS,
                  multiprocessingMode=THREADS, numCpus=n_proc,
                  explicitPrecision=DOUBLE_PLUS_PACK, 
                  nodalOutputPrecision=FULL, numDomains=n_proc)
    job.submit(consistencyChecking=OFF)
    # waitForCompletion() crashes, if run interactively
    #  --> try opening window to manually check when job is finished
    try:
        getWarningReply('Wait until job has finished (see .sta file), then press Yes.',buttons=(YES,NO))
    except:
        job.waitForCompletion()
    return

# Draw the Lego(r) parts and widen their cavities -----------------------------------

def make_middle_pos(x_arr, z_arr):
    """Create mid positions between the points in a regular x-z grid and mid positions every other stud for the bottom walls in the cavities.

    Args:
    x_arr (np.ndarray): Array of x positions of studs
    z_arr (np.ndarray): Array of z positions of studs
    
    Returns:
    x_mid_arr (np.array): Mid positions in x direction, len(x_mid_arr) = len(x_arr)-1
    z_mid_arr (np.array): Mid positions in z direction, len(z_mid_arr) = len(z_arr)-1
    """
    if len(x_arr) == 1:
        if len(z_arr) == 1:
            # if only 1x1 brick: no mid positions
            return None, None
        else:
            # for 1xn brick: no x-mid positions
            delta_z = z_arr[1] - z_arr[0]
            z_mid_arr = z_arr[:-1] + delta_z / 2.
            return None, z_mid_arr
    elif len(z_arr) == 1:
        # for nx1 brick: no z-mid positions
        delta_x = x_arr[1] - x_arr[0]
        x_mid_arr = x_arr[:-1] + delta_x / 2.
        return x_mid_arr, None
    else:
        delta_x = x_arr[1] - x_arr[0]
        delta_z = z_arr[1] - z_arr[0]
        #
        x_mid_arr = x_arr[:-1] + delta_x / 2.
        z_mid_arr = z_arr[:-1] + delta_z / 2.
    return x_mid_arr, z_mid_arr


def make_u_nodes(p, nx, nz, b, r_stud, h_stud, mesh_size):
    """Calculate necessary displacements of nodes of the brick cavities such 
    that there is no contact penetration between the stud and the cavity
    nodes.

    Args:
        p (Abaqus part): The part whose cavities should be widened
        nx (int): Number of studs in x direction
        nz (int): Number of studs in z direction
        b (float): Distance between two studs
        r_stud (float): Stud radius
        h_stud (float): Stud height

    Returns:
        ndarray: The node label, the radial displacement, and the negative angle as [[i_node,u_r,-u_ang],...]
    """
    # calculate the displacements of nodes that may contact studs
    # to establish initial contact
    def get_overlap_nodes(c_arr, x0, z0, r_stud, h_stud, mesh_size):
        """Calculate radial displacement to resolve contact penetration.

        Args:
            c_arr (ndarray): Array of node labels and node coordinates [[n_i,x_i,y_i,z_i],...]
            x0 (float): x center of stud
            z0 (float): z center of stud
            r_stud (float): Stud radius
            h_stud (float): Stud height

        Returns:
            list: The node label, the radial displacement, and the negative angle as [[i_node,u_r,-u_ang],...]
        """
        # for one stud: get Nodes that overlap
        # (including their penetration and direction)
        c_distance = c_arr - [0, x0, z0, 0]
        c_radius = (c_distance[:, 1] ** 2 + c_distance[:, 2] ** 2) ** 0.5
        c_angle = np.arctan2(c_distance[:, 2], c_distance[:, 1])
        # select overlapping nodes
        cyl_sel = (c_radius < r_stud) * (c_distance[:, -1] <= h_stud + mesh_size/2.)
        # obtain node number, penetration, angle (rad)
        n_arr = c_arr[cyl_sel][:, 0]
        r_arr = c_radius[cyl_sel] - r_stud
        ang_arr = c_angle[cyl_sel] * 180 / np.pi
        return list(np.column_stack((n_arr, -r_arr, ang_arr)))

    # output the nodes for the initial displacements
    nodes_cont_bot = p.sets['cont-bot-touch'].nodes
    
    nodes_cont_bot = np.array([np.array([i.label] + list(i.coordinates), dtype=float)
                               for i in nodes_cont_bot])
    
    move_nodes = []
    for j in range(nz):
        for i in range(nx):
            move_nodes += get_overlap_nodes(nodes_cont_bot, b * i, b * j,
                                            r_stud, h_stud, mesh_size)
    move_nodes = np.array(move_nodes)

    # create node sets
    for n_label_i in move_nodes[:, 0]:
        n_label = int(n_label_i)
        p.Set(name='x-n' + str(n_label), nodes=p.nodes[n_label - 1:n_label])
    
    # for debugging
    #np.savetxt('nodes_to_move_{}x{}.dat'.format(nz, nx), move_nodes)
    return move_nodes


# Find interacting Lego(r) parts ----------------------------------------------------

def make_brick_lists(brick, h, h_top):
    """Make arrays of the (ix,iz) positions of the studs of a brick and the y_range of the side wall of the brick

    Args:
        brick (dict): Dictionary that contains the brick `type`, `nx`, and `nz`
        h (float): Height of Lego brick without stud
        h_top (float): Height of Lego stud

    Returns:
        (ndarray, ndarray): Array of (ix,iz) Positions of studs and array of the y range of the brick without stud.
    """
    # create list of stud positions and the y range of the brick
    # (without moving the brick to position in model)
    b_type = brick['type']
    nx, nz = brick['nx'], brick['nz']

    # list of stud positions (in (8 mm))
    brick_list = ()
    for ix in range(1, nx+1):
        for iz in range(1, nz+1):
            brick_list += ((ix, iz),)
    
    # y range (without studs)
    if b_type == 'regular':
        y_range = (0, h)
    elif b_type == 'base-plate':
        y_range = (-h_top, 0)
    else:
        y_range = (0, h/3.)
    return np.array(brick_list), np.array(y_range)


def make_parts(bricks, parts, b, h, h_top, if_plot=0):
    """Add the stud positions and y ranges of all parts to the `parts` dictionary

    Args:
        bricks (dict): Definition of the bricks used in the model
        parts (dict): Definition of the positions of the bricks in the model (the `parts` sub-dictionary of the assembly)
        b (float): Distance between two studs
        h (float): Height of Lego brick without stud
        h_top (float): Height of Lego stud
        if_plot (int, optional): If all Lego parts should be plotted. Defaults to 0.

    Returns:
        dict: assembly: `parts` dictionary with `stud_list` and `y_range` in global coordinate system added
    """
    # load parts from 'parts' dictionary:
    # add y_range and stud_list to dictionary
    if if_plot:
        fig, (ax1,ax2) = plt.subplots(nrows=2, sharex=True)

    # make all parts
    for i_part in parts.keys():
        i_brick = parts[i_part]['brick_id']
        loc = parts[i_part]['loc']
        brick_list, y_range = make_brick_lists(bricks[i_brick], h, h_top)
        brick_list += (int(loc[0]/b), int(loc[2]/b))
        y_range += loc[1]

        # add information to ass dict
        parts[i_part]['stud_list'] = np.copy(brick_list)
        parts[i_part]['y_range'] = np.copy(y_range)
        
        # rough plot of the part positions (y over x and z over x)
        if if_plot:
            ax1.plot(list(brick_list[:,0])*2, [y_range[0]]*len(brick_list)+
                     [y_range[1]]*len(brick_list), 'x-', markersize=20)
            ax2.plot(brick_list[:,0], brick_list[:,1], 'x-', markersize=20)
    
    if if_plot:
        ax1.set_ylabel('y (9.6 mm)')
        ax2.set_ylabel('z (8 mm)')
        plt.xlabel('x (8 mm)')
        plt.show()
    return parts


def check_two_parts(part1, part2, h_stud, b):
    """Check relative positions between two bricks and returns the lower cavities of part2 to be widened and the  surfaces for contact between part1 and part2.

    Args:
        part1 (dict): Dictionary defining part1: Range of y coordinate and positions of studs
        part2 (dict): Dictionary defining part2: Range of y coordinate and positions of studs
        h_stud (float): Stud height
        b (float): Stud distance

    Returns:
        (str, list, list, int)
        y_rel (str): The relative position of part 2 to part 1. Can be 'top', 'bottom', 'side' or 'none'.
        widen_list (list): (xi,z1) positions of cavities that should be widened for part2
        side_list (list): [[str1,str2]] strings for contact surface names of part1 and part2, respectively
        if_stud_contact (int): One if at least one cavity of part2 is widened, otherwise 0
    """
    # if part1 and part2 are identical: return empty lists
    if part1['instance'][0] == part2['instance'][0]:
        return ('none', [], [])
    
    s_list1, s_list2 = part1['stud_list'], part2['stud_list']
    loc = part1['loc']
    y_range1, y_range2 = part1['y_range'], part2['y_range']

    # check for overlap in y-direction    
    if y_range1[1] + h_stud > y_range2[0] and y_range1[1] <= y_range2[0]+TOL:
        # part2 on top of part1
        y_rel = 'top'
    elif y_range2[1] + h_stud > y_range1[0] and y_range2[1] <= y_range1[0]+TOL:
        # part1 on top of part2
        y_rel = 'bot'
    elif ((y_range2[1] > y_range1[0] + TOL and y_range2[1] < y_range1[1] + TOL) or
          (y_range2[0] < y_range1[1] + TOL and y_range2[0] > y_range1[0] - TOL)):
        # possible side contact of part1 and part2
        y_rel = 'side'
    else:
        y_rel = 'none'
    
    # check if overlap in x,z plane (what cavity positions of part2 should be widened?)
    widen_list = []
    if_stud_contact = 0

    if y_rel == 'bot':
        for s1 in s_list1:
            for s2 in s_list2:
                # check bottom positions that should be widened in first part
                if list(s1) == list(s2):
                    widen_list += [list(s1 - (int(loc[0] / b), int(loc[2] / b)))]
                    if_stud_contact = 1

    # check for side contact in x,z plane
    side_list = []

    if y_rel == 'side':
        for s1 in s_list1:
            for s2 in s_list2:
                # only in positive direction: because checked from both sides
                # (side contact not relevant for overlapping bricks)
                if (s2 - s1)[0] == 1 and (s2-s1)[1] == 0:
                    side_list += [['x1', 'x0']]
                elif (s2 - s1)[1] == 1 and (s2-s1)[0] == 0:
                    side_list += [['z1', 'z0']]
    
    return y_rel, widen_list, side_list, if_stud_contact


def make_widen_cont(assembly, lego_geom, if_plot=0):
    """Check what faces of what parts are in contact (for implicit load step) and what lower cavity positions of what bricks are widened.

    Args:
        assembly (dict): The setup of the Lego set defining bricks used (sub-dictionary `bricks`), part positions (`parts`), boundary conditions (`bc`), and loads applied on sets (`loads_rp`)
        lego_geom (dict): Dictionary containing the general Lego dimensions.
        if_plot (int, optional): If all Lego parts should be plotted. Defaults to 0.

    Returns:
        (dict_widen, cont_list)
        dict_widen (dict): (xi,z1) positions of cavities that should be widened for all bricks
        cont_list (list): [i_brick1, i_brick2, [str1, str2]]: Contact should be defined between the bricks with index `i_brick1` and `i_brick2` with their surfaces with string `str1` and `str2`, respectively. 

    Examples:
        lego_geom: {'b, b_gap, b_wall': (8., 0.1, 1.6),
                    'h, h_stud, h_top': (9.6, 1.7, 1.5),
                    'inside small': {'r': 1.6},
                    'inside big': {'r': 3.3, 't': 0.9},
                    'delta_r': 0.05}
        assembly: {'name': 'wall1x2-down-dev',
                'bricks': {1: {'type': 'regular', 'nx': 2, 'nz': 1}},
                'parts': {1: {'brick_id': 1, 'loc': (0, 0, 0)},
                        2: {'brick_id': 1, 'loc': (16, 0, 0)},
                        3: {'brick_id': 1, 'loc': (32, 0, 0)},
                        4: {'brick_id': 1, 'loc': (48, 0, 0)},
                        5: {'brick_id': 1, 'loc': (8, -9.6, 0)},
                        6: {'brick_id': 1, 'loc': (24, -9.6, 0)},
                        7: {'brick_id': 1, 'loc': (40, -9.6, 0)}}}
        dict_widen: {1: [[2, 1]], 2: [[1, 1], [2, 1]], 3: [[1, 1], [2, 1]], 4: [[1, 1]],
                     5: [], 6: [], 7: []}
        cont_list: [[1, 2, ['x1', 'x0']], [1, 5, ['bot', 'top']], [2, 3, ['x1', 'x0']], [2, 5, ['bot', 'top']],
                    [2, 6, ['bot', 'top']], [3, 4, ['x1', 'x0']], [3, 6, ['bot', 'top']], [3, 7, ['bot', 'top']],
                    [4, 7, ['bot', 'top']], [5, 6, ['x1', 'x0']], [6, 7, ['x1', 'x0']]]
    """
    b, _, _ = lego_geom['b, b_gap, b_wall']
    h, h_stud, h_top = lego_geom['h, h_stud, h_top']

    # flatten the cavity lists that should be widened
    def flat_cavity(lists_0):
        out_list = []
        for i in lists_0:
            out_list += i
        return out_list

    # obtain all stud positions for widening (where studs will stick) and contact pairs
    parts = make_parts(assembly['bricks'], assembly['parts'], b, h, h_top, if_plot)
    
    dict_widen = {}
    cont_list = []

    # check two parts each and fill the dict_widen dictionary and the cont list
    for k_part in parts.keys():
        dic_temp = {j:check_two_parts(parts[k_part], i, h_stud, b) for j, i in parts.items() if j != k_part}
        #print(dic_temp)
        widen_list = [i[1] for i in dic_temp.values() if i[1] != []]
        dict_widen[k_part] = flat_cavity(widen_list)
        # (first part, second part, faces)
        cont_list += [[k_part, i_part, j[2][0]] for i_part, j in dic_temp.items() if j[2] != []]
        # add stud contact
        # (contact identifyers: 'x0', 'x1', 'z0', 'z1', 'bot', 'top')
        cont_list += [[k_part, i_part, ['bot','top']] for i_part, j in dic_temp.items() if j[-1] == 1]
    
    #with open('check-widen.json', 'w') as f:
    #    json.dump({'widen':dict_widen, 'cont_list':cont_list}, f)
    
    #raise ValueError('')
    
    return dict_widen, cont_list


# ABAQUS evaluation functions -------------------------------------------------------

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


def color_bricks(vp1,instances,rigid_color='ffffff'):
    """Assign colors to the displayed instances either from the defined colors in the `assembly` dictionary or as random colors.

    Args:
        vp1 (Abaqus viewport): The viewport, needed to set the display colors
        instances (dict): The instances of the model or the odb
    """

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
    
    # load the assembly dictionary from file
    with open('_dict-assembly.json', 'r') as f:
        assembly = json.load(f)['assembly']
    
    def get_color(part,cm_set,c_default='ffffff'):  # 7ed986
        """Get the hex code of color from col_str: can be a string with the color name or the hex code itself

        Args:
            part (dict): The part definition in the `assembly` dictionary
            cm_set (set): Set of color names in `cm_lego`
            c_default (str, optional): Default color in its hex code. Defaults to 'ffffff'.
        """
        # get the hex code of color from col_str
        if 'c' in part.keys():
            # use color string stated in part
            col_str = part['c']

            if col_str in cm_set:
                # if col_str in cm_set: use hex code from cm_lego
                return cm_lego[col_str]
            elif len(col_str) == 6:
                # if six digits: directly use as hex code
                return col_str
            else:
                # otherwise return default color
                return c_default
        else:
            # otherwise return default color
            return c_default

    # check if there are colors stated in the assembly dictionary
    if "c" in assembly['parts']['1'].keys():
        # if 'c' keys in the parts: use those colors for video
        cm_set = set(cm_lego.keys())
        print(assembly['parts'].items())
        # add color string for the rigid part(s) 'd5ef5b' for light green
        col_list = [get_color(assembly['parts'][str(i_part+1)],cm_set) for i_part in range(len(assembly['parts']))]+[rigid_color]*4

    else:
        # if no colors are defined in parts of assembly dictionary:
        #    use random colors
        try:
            # if file exists, use the file for color scheme
            with open('_color-file.dat', 'r') as f:
                col_list = f.read().split(', ')
            print('loaded')
        except:
            # if file does not exist: shuffle colors and write file
            # use the color list 5 times such that more than 11 bricks can be used
            col_list = list(cm_lego.values())
            random.shuffle(col_list)
            with open('_color-file.dat', 'w') as f:
                f.write(', '.join(col_list))

    # if brick is divided into two parts, it should still have the same color
    inst_names = instances.keys()
    i = 0
    cmap_dict = {inst_names[0]:(True, '#' + col_list[i], 'Default', '#' + col_list[i])}
    i_name_temp = inst_names[0]
    #
    for inst_name in inst_names[1:]:
        if inst_name[:-4] != i_name_temp[:-4]:
            i += 1
        cmap_dict[inst_name] = (True, '#' + col_list[i], 'Default', '#' + col_list[i])
        i_name_temp = inst_name
    
    # set the color mapping
    cmap = vp1.colorMappings['Part instance']
    vp1.setColor(colorMapping=cmap)
    vp1.disableMultipleColors()
    
    cmap.updateOverrides(overrides=cmap_dict)
    cmap = vp1.colorMappings['Part instance']
    vp1.setColor(colorMapping=cmap)
    vp1.disableMultipleColors()
    return


def print_assembly(model, job_name):
    """For debugging. Print the assembly of the model to a png image named {job_name}-ass.png
    """
    ass = model.rootAssembly

    # display the assembly
    vp1 = session.viewports['Viewport: 1']
    vp1.setValues(displayedObject=ass)
    vp1.assemblyDisplay.setValues(optimizationTasks=OFF, geometricRestrictions=OFF,
                                  stopConditions=OFF)
    
    # select the iso view and fit model into the viewport
    vp1.view.setValues(session.views['Iso'])
    vp1.view.fitView()

    # assign colors to the bricks
    color_bricks(vp1,model.rootAssembly.instances)

    # Print the image to a png file
    session.printOptions.setValues(vpDecorations=OFF, reduceColors=False)
    session.printToFile(fileName=job_name + '-ass', format=PNG, canvasObjects=(vp1,))
    return


def create_video(job_name):
    """Create a video of the deformed model in the load step. The deformed
    bricks are plotted in different colors and the mesh is not shown. Separate 
    images are written to the subfolder `img-video` and then converted to an 
    animated gif using `imagemagick`. The resolution of the images is set to
    960 x 960 pixels.
    """
    odb = session.openOdb(job_name + '.odb')
    vp1 = session.viewports['Viewport: 1']
    vp1.setValues(displayedObject=odb)

    # resore viewport to set the aspect ratio of the window for image
    vp1.makeCurrent()
    vp1.restore()
    vp1.setValues(width=160, height=160)
    vp1.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF,))
    vp1.odbDisplay.setFrame(step=0, frame=-1)

    # possibility to use a manual view
    session.View(name='User-1', nearPlane=119.04, farPlane=181.08, width=74.219, height=37.301,
                 projection=PERSPECTIVE, cameraPosition=(102.76, -116.14, 52.059), viewOffsetY=0,
                 cameraUpVector=(-0.34279, 0.52332, 0.78015), cameraTarget=(23.063, -2.4859, 3.6211),
                 viewOffsetX=0, autoFit=ON)
    
    # for the pumpkin
    session.View(name='User-1', nearPlane=221.05, farPlane=323.41, width=203.45, 
                height=99.922, projection=PERSPECTIVE, cameraPosition=(104.03, 92.217, 
                251.82), cameraUpVector=(-0.13447, 0.85535, -0.5003), cameraTarget=(31.195, 
                36.737, -5.3113), viewOffsetX=0, viewOffsetY=0, autoFit=ON)
    vp1.view.setValues(session.views['User-1'])

    # just use the iso view and fit model into the window
    #vp1.view.setValues(session.views['Iso'])
    #vp1.view.fitView()

    vp1.viewportAnnotationOptions.setValues(triad=OFF, state=OFF, legendBackgroundStyle=MATCH,
                                            annotations=OFF, compass=OFF, title=OFF, legend=OFF)
    #
    session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')

    session.printOptions.setValues(vpDecorations=OFF, vpBackground=ON, reduceColors=True)
    vp1.odbDisplay.commonOptions.setValues(visibleEdges=FREE)
    # create folder for the images
    make_dir('img-video', if_clear=1, if_change=0)
    vp1.odbDisplay.display.setValues(plotState=(DEFORMED,))
    #
    # assign colors to the bricks
    color_bricks(vp1,odb.rootAssembly.instances)
    
    # set the rezolution of the animated gif (160*6)
    session.pngOptions.setValues(imageSize=(160*6, 160*6))
    
    # print png images for all frames in the load step
    for i_frame in range(len(odb.steps['load'].frames)):
        try:        # for implicit load
            vp1.odbDisplay.setFrame(step=3, frame=i_frame)
        except:     # for explicit load
            vp1.odbDisplay.setFrame(step=0, frame=i_frame)
        # print image to files anim-000.png, anim-002.png, ...
        session.printToFile(fileName='img-video/anim-' + str(i_frame).zfill(3),
                            format=PNG, canvasObjects=(vp1,))
    
    # new version with ffmpeg (works only ok for reducedColors=True)
    try:
        os.remove('anim-'+job_name+'.gif')
    except:
        None
    os.system('ffmpeg -framerate 30 -i img-video/anim-%03d.png anim-'+job_name+'.gif')

    make_dir('img-video', if_clear=1, if_change=0)  
    vp1.maximize()
    
    return
    # merge the images into an animated gif: use a delay of 50 ms (-delay 5) 
    # & repeat (-loop 0)
    os.system('magick convert -dispose Background -delay 5 img-video/anim-*.png -loop 0 ' +
              job_name + '-anim.gif')
    make_dir('img-video', if_clear=1, if_change=0)
    vp1.maximize()
    return


def get_ho(job_name):
    """Load and export the history output of the model from the odb.

    Yields
    ------
    {job_name}-ho-ass.dat (ASCII file): All relevant energy outputs are written over the time into that file: The total strain energy ALLSE, the total kinetic energy ALLKE, the total work of external forces ALLWK, and the total frictional dissipation ALLFD.
    {job_name}-ho-node{i}.dat (ASCII files): For each reference point `i`, such a file is written into a .dat file that contains the time, displacements u, and reaction forces rf.
    """
    odb = session.openOdb(job_name + '.odb')
    step1 = odb.steps['load']

    # obtain the history regions of the reference points
    hr_node = [hr for hr in step1.historyRegions.values() if 'Node' in hr.name]

    # go through the reference point history regions and load the history output
    for i_hr, hr in enumerate(hr_node):
        ho_dic = {}
        ho_dic['time'] = np.array(hr.historyOutputs.values()[0].data)[:, 0]
        for ho_name, ho in hr.historyOutputs.items():
            ho_dic[ho_name] = np.array(ho.data)[:, 1]
        # write results to dat file
        with open(job_name + '-ho-node' + str(i_hr).zfill(2) + '.dat', 'w') as f:
            f.write(', '.join(ho_dic.keys()) + '\n')
            np.savetxt(f, np.column_stack(ho_dic.values()), delimiter=', ')

    # access energy in the model
    hr_ass = [hr for hr in step1.historyRegions.values() if 'ASSEMBLY' in hr.name][0]
    ho_ass = np.column_stack((np.array(hr_ass.historyOutputs['ALLSE'].data)[:,0],
                              np.array(hr_ass.historyOutputs['ALLSE'].data)[:,1],
                              np.array(hr_ass.historyOutputs['ALLKE'].data)[:,1],
                              np.array(hr_ass.historyOutputs['ALLWK'].data)[:,1],
                              np.array(hr_ass.historyOutputs['ALLFD'].data)[:,1]))

    # write energies to .dat file
    with open(job_name + '-ho-ass.dat', 'w') as f:
        f.write('time, ALLSE, ALLKE, ALLWK, ALLFD\n')
        np.savetxt(f, ho_ass, delimiter=', ')
    return


# ABAQUS helper functions -----------------------------------------------------------

def make_rigid_part(model, i_rigid, rigid_par, mesh_size):
    """Create rigid part for load (in `loads_rigid`)

    Args:
        model (Abaqus model)
        i_rigid (int): Number of load for the name of rigid part
        rigid_par (dict): Definition of the rigid part geometry, can bei either a sphere or a cylinder: {'shape':'sphere', 'loc':(0,0,0), 'radius':4., 'u':(0,0,0)}, {'shape':'cyl', 'loc':(0,0,0), 'radius':4., 'u':(0,0,-4), 'dir':(0,1,0), 'len':15.}
        mesh_size (float): Mesh size of the rigid part

    Returns:
        model instance: The instance of the created rigid part
    """
    # create sketch, part, and a reference point in the middle
    s_rigid = model.ConstrainedSketch(name='rigid', sheetSize=200.0)
    
    if rigid_par['shape'] == 'cyl':
        p_name = 'rigid-cyl-' + str(i_rigid)
        s_rigid.CircleByCenterPerimeter(center=(0, 0), point1=(rigid_par['radius'], 0))
        p_rigid = model.Part(name=p_name, dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
        p_rigid.BaseShellExtrude(sketch=s_rigid, depth=rigid_par['len'])
        
        rp = p_rigid.ReferencePoint(point=(0, 0, rigid_par['len'] / 2.))
    else:
        p_name = 'rigid-sphere-' + str(i_rigid)
        s_rigid.ConstructionLine(point1=(0, 0), point2=(0, 1))
        s_rigid.ArcByCenterEnds(center=(0, 0), point1=(0, rigid_par['radius']),
                                point2=(0, -rigid_par['radius']), direction=CLOCKWISE)
        p_rigid = model.Part(name=p_name, dimensionality=THREE_D,
                             type=DISCRETE_RIGID_SURFACE)
        p_rigid.BaseShellRevolve(sketch=s_rigid, angle=360.0, flipRevolveDirection=OFF)
        
        rp = p_rigid.ReferencePoint(point=(0, 0, 0))
    
    # seed and mesh rigid part
    p_rigid.seedPart(size=mesh_size)
    p_rigid.generateMesh()
    
    # create instance and move to given position
    ass = model.rootAssembly
    inst_rigid = ass.Instance(name=p_name, part=p_rigid, dependent=ON)
    
    if rigid_par['shape'] == 'cyl':
        ass.translate(instanceList=(p_name,), vector=(0, 0, -rigid_par['len'] / 2.))
        if rigid_par['dir'] == (0, 1, 0):
            ass.rotate(instanceList=(p_name,), axisPoint=(0, 0, 0),
                     axisDirection=(1, 0, 0), angle=90)
        elif rigid_par['dir'] == (1, 0, 0):
            ass.rotate(instanceList=(p_name,), axisPoint=(0, 0, 0),
                     axisDirection=(0, 1, 0), angle=90)
    
    ass.translate(instanceList=(inst_rigid.name,), vector=rigid_par['loc'])

    # define RP set and the contact surface
    p_rigid.Set(name='RP', referencePoints=(p_rigid.referencePoints[rp.id],))
    ass.Surface(side1Faces=inst_rigid.faces, name='rigid-contact')

    # apply mass in reference point
    if 'm' in rigid_par.keys():
        if rigid_par['shape'] == 'cyl':
            i11 = 0.5*rigid_par['m']*rigid_par['radius']**2
            i22, i33 = 1/12.*rigid_par['m']*rigid_par['len']**2
        else:
            i11, i22, i33 = [2/5.*rigid_par['m']*rigid_par['radius']**2]*3
        
        p_rigid.engineeringFeatures.PointMassInertia(name='Inertia-1', region=p_rigid.sets['RP'],
                                                     mass=rigid_par['m'], i11=i11, i22=i22, i33=i33, 
                                                     alpha=0.0, composite=0.0)

    return inst_rigid


def make_contact(model, mu, pos_dict={}, cont_pairs=[], is_expl=1, if_aug_lag=0):
    """Creates the contact property and interactions in the Abaqus model.

    Args:
        model (Abaqus model)
        mu (float): The friction coefficient (both static and dynamic)
        pos_dict (dict, optional): Dictionary that defines the Lego brick positions. Defaults to {}.
        cont_pairs (list, optional): [n1,n2,[str_1,str_2]],...] with n1 and n2 the instance numbers of the bricks and str_1 and str_2 the surface names for the contact. Only relevant for is_expl==0. Defaults to [].
        is_expl (int, optional): If the model is explicit (1) or implicit (0). Defaults to 1.
        if_aug_lag (int, optional): If 1, use Augumented Lagrange contact with default properties. Otherwise, Penalty Contact is used. Defaults to 0.
    
    Example:
        cont_pairs : [[1, 2, ['x1', 'x0']], [1, 5, ['bot', 'top']], [2, 3, ['x1', 'x0']], [2, 5, ['bot', 'top']], [2, 6, ['bot', 'top']],]]
    """
    def make_std_interaction(model,name,s1,s2,ip='cont-prop'):
        """Create the interaction, keywords have changed in Abaqus versions
        """
        try:
            # earlier Abaqus versions
            model.SurfaceToSurfaceContactStd(name=name, createStepName='contact', master=s1, slave=s2,
                                             sliding=FINITE, interactionProperty=ip)
        except:
            # Abaqus version 2022: master/slave is called main/secondary
            model.SurfaceToSurfaceContactStd(name=name, createStepName='contact', main=s1, secondary=s2,
                                             sliding=FINITE, interactionProperty=ip)
        return

    # make contact properties
    cp = model.ContactProperty('cont-prop')
    if not if_aug_lag:
        cp.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON,
                          constraintEnforcementMethod=DEFAULT)
    else:
        cp.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, contactStiffness=DEFAULT,
                          contactStiffnessScaleFactor=1.0, clearanceAtZeroContactPressure=0.0,
                          constraintEnforcementMethod=AUGMENTED_LAGRANGE)
    if mu == 0.:
        cp.TangentialBehavior(formulation=FRICTIONLESS)
    else:
        cp.TangentialBehavior(formulation=PENALTY, table=((mu,),),
                              maximumElasticSlip=FRACTION, fraction=0.005)
    
    # make interactions
    if is_expl:
        gen_cont = model.ContactExp(name='general-contact', createStepName='Initial')
        gen_cont.includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)
        gen_cont.contactPropertyAssignments.appendInStep(stepName='Initial',
                                                         assignments=((GLOBAL, SELF, 'cont-prop'),))
    else:
        # use cont_pairs. if pos_dict[i_inst2]['instance'][0] == pos_dict[i_inst2]['instance'][1]: 
        # side cont only once, otherwise twice!!
        for i_1, i_2, (str_1, str_2) in cont_pairs:
            inst_1_bot, inst_1_top = pos_dict[i_1]['instance']
            inst_2_bot, inst_2_top = pos_dict[i_2]['instance']
            #
            if str_1 == 'bot':
                make_std_interaction(model, 'cont-bt-'+str(i_1)+'-'+str(i_2), inst_2_top.surfaces['contact-top'],
                                     inst_1_bot.surfaces['contact-bot'])
            else:
                # side contact: x,z. Gets weird for 2-part bricks: 4 contact combinations
                make_std_interaction(model, 'cont-side-'+str(i_1)+'-'+str(i_2)+'bb', inst_1_bot.surfaces[str_1],
                                     inst_2_bot.surfaces[str_2])
                if inst_1_bot != inst_1_top:
                    make_std_interaction(model, 'cont-side-' + str(i_1) + '-' + str(i_2)+'tb',
                                         inst_1_top.surfaces[str_1], inst_2_bot.surfaces[str_2])
                    if inst_2_bot != inst_2_top:
                        make_std_interaction(model, 'cont-side-' + str(i_1) + '-' + str(i_2)+'tt',
                                             inst_1_top.surfaces[str_1], inst_2_top.surfaces[str_2])
                if inst_2_bot != inst_2_top:
                    make_std_interaction(model, 'cont-side-' + str(i_1) + '-' + str(i_2)+'bt',
                                         inst_1_bot.surfaces[str_1], inst_2_top.surfaces[str_2])
    return


def make_sections(model, mat_props):
    """
    Create linear elastic material with (Young's modulus E,Poisson's ratio nu, density) = mat_props. Create and assign solid section to all parts in the model.
    """
    # create the lin.elastic material for ABS
    mat = model.Material(name='ABS')
    mat.Elastic(table=(mat_props[:2],))
    mat.Density(table=((mat_props[-1],),))

    # define solid section and assign to all parts
    model.HomogeneousSolidSection(name='ABS', material='ABS',
                                  thickness=None)
    for p in model.parts.values():
        p.SectionAssignment(region=p.sets['all'], sectionName='ABS',
                            thicknessAssignment=FROM_SECTION)
    return


def load_impl_to_expl(model, job_name):
    """
    Load the deformed parts of the initial implicit model with result files `job_name` to the model `model` and create instances from the loaded parts. Returns the list `inst_names` that contains all instance names.
    """
    # open odb and obtain names of all instances (no assembly stuff like RP, then the model crashes)
    odb = session.openOdb(job_name + '.odb')
    a = model.rootAssembly
    inst_names = [i for i in odb.rootAssembly.instances.keys() if i != 'ASSEMBLY']
    
    # load all parts and create instances out of them
    for inst_name in inst_names:
        # files necessary: .odb, .res, .prt, .mdl, ...
        p_temp = model.PartFromOdb(name=inst_name, instance=inst_name, odb=odb,
                                   shape=DEFORMED, step=-1, frame=-1)
        a.Instance(name=p_temp.name, part=p_temp, dependent=ON)
    odb.close()
    return inst_names


def make_load_rps(model, loads_rp):
    """Create reference points and rigidly connect them to the corresponding set of a part.

    Args:
        model (Abaqus model)
        loads_rp (dict): Definition what part and what set should be used, e.g. {1:{'part_id':1, 'set_name':'STUD-12', 'ux':3.},...}
    """
    ass = model.rootAssembly
    for i, load_val in loads_rp.items():

        # create RP and rigid body coupling
        set_name = load_val['set_name']

        # take the right instance for the coupling to the RP
        # (can be either top or bottom part)
        inst_i_top = ass.instances['BRICK' + str(load_val['part_id']).zfill(2) + '-TOP']
        
        if set_name in inst_i_top.sets.keys():
            load_set = inst_i_top.sets[set_name]
        else:
            inst_i_bot = ass.instances['BRICK' + str(load_val['part_id']).zfill(2) + '-BOT']
            load_set = inst_i_bot.sets[set_name]
        
        try:  # in the implicit model: faces exist
            load_face = load_set.faces[0]
            rp_position = load_face.getCentroid()[0]
        except:  # in the explicit model: loaded mesh only has nodes and elements
            elem_bb = load_set.elements.getBoundingBox()
            rp_position = ((elem_bb['high'][0] + elem_bb['low'][0]) / 2.,
                           (elem_bb['high'][1] + elem_bb['low'][1]) / 2.,
                           (elem_bb['high'][2] + elem_bb['low'][2]) / 2.)

        # create RP and set out of it
        rp = ass.ReferencePoint(point=rp_position)
        ass.Set(name='RP-' + str(i), referencePoints=(ass.referencePoints[rp.id],))
        
        # create rigid body coupling
        model.RigidBody(name='RigBody-RP-' + str(i), pinRegion=load_set,
                        refPointRegion=(ass.referencePoints[rp.id],))
    return


def make_loads(model, assembly, explicit_par, mesh_size, t_step, is_expl=1, if_prestep=0):
    """Create the loads in the implicit or explicit model on the reference points that are connected to sets. Furthermore, load the reference points of rigid bodies when `if_prestep==0`. The load steps must exist already.

    Args:
        model (Abaqus model)
        assembly (dict): Arrangement of the bricks in the model, including boundary conditions and loads on reference points that are connected to sets.
        explicit_par (dict): Load parameters of the model; Uses `if_acc` (if the displacement should be applied using a constant acceleration instead of a constant velocity) and `loads_rigid` (the dictionary that defines the rigid bodies)
        mesh_size (float): The global mesh size of the model, needed for meshing the rigid parts
        t_step (float): Time to be used for the load step, needed for calculating the acceleration
        is_expl (int, optional):  Wether the load step is an explicit step.. Defaults to 1.
        if_prestep (int, optional): If the function is called for creating the loads in the clamping steps (`widen`, `contact`, `free`) or from the load step. Because for an explicit load, the reference points of the loads_rp should be fixed in the `free` step to avoid rigid body motion.. Defaults to 0.
    """
    ass = model.rootAssembly
    
    if is_expl:
        if_acc = explicit_par['is_acc']
        rigid_loads = explicit_par['loads_rigid']
    #
    # rigid body constraints at studs (RPs for `loads_rp` already exist and are 
    # in sets 'RP-{i}' with i being the index of the `loads_rp` dictionary)
    for i, load_val in assembly['loads_rp'].items():

        # create lists and fill with UNSET keyword (or FREED for implicit load step)
        load_list = [UNSET] * 6
        load_list0 = [UNSET] * 6
        load_list_impl = [FREED] * 6
        #
        for i_load, load_str in enumerate(('ux', 'uy', 'uz', 'rotx', 'roty', 'rotz')):
            if load_str in load_val.keys():
                load_list[i_load] = load_val[load_str]
                load_list_impl[i_load] = load_val[load_str]
                load_list0[i_load] = 0

        def multi_load(load, fac):
            """ Only multiply if load value is not UNSET
            """
            if load == UNSET:
                return UNSET
            else:
                return load * fac
        
        step_ho = 'load'

        # the RPs have already been created and are in assembly sets named 'RP-{i}'
        if is_expl and if_prestep == 0:
            if if_acc:
                acc_fac = 2 / t_step ** 2
                model.AccelerationBC(name='load-RP-' + str(i), createStepName='load', region=ass.sets['RP-' + str(i)],
                                     a1=multi_load(load_list[0],acc_fac), a2=multi_load(load_list[1],acc_fac),
                                     a3=multi_load(load_list[2],acc_fac), ar1=multi_load(load_list[3],acc_fac),
                                     ar2=multi_load(load_list[4],acc_fac), ar3=multi_load(load_list[5],acc_fac))
            else:
                model.VelocityBC(name='load-RP-' + str(i), createStepName='load', region=ass.sets['RP-' + str(i)],
                                 v1=multi_load(load_list[0], 1/t_step), v2=multi_load(load_list[1], 1/t_step),
                                 v3=multi_load(load_list[2], 1/t_step), vr1=multi_load(load_list[3], 1/t_step),
                                 vr2=multi_load(load_list[4], 1/t_step), vr3=multi_load(load_list[5], 1/t_step))
        else:
            # first create, than set values in step
            bc_temp = model.DisplacementBC(name='load-RP-' + str(i), createStepName='free', region=ass.sets['RP-' + str(i)],
                                           u1=load_list0[0], u2=load_list0[1], u3=load_list0[2], ur1=load_list0[3],
                                           ur2=load_list0[4], ur3=load_list0[5])
            if is_expl == 0:
                bc_temp.setValuesInStep(stepName='load', u1=load_list_impl[0], u2=load_list_impl[1], u3=load_list_impl[2],
                                        ur1=load_list_impl[3], ur2=load_list_impl[4], ur3=load_list_impl[5])
            else:
                step_ho = 'free'
        #
        ho_var = ('U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3', 'RF1', 'RF2', 'RF3', 'RM1', 'RM2', 'RM3')
        if is_expl and if_prestep == 0:
            model.HistoryOutputRequest(name='ho-RP-' + str(i), createStepName=step_ho, numIntervals=500,
                                       variables=ho_var,
                                       region=ass.sets['RP-' + str(i)])
        else:
            model.HistoryOutputRequest(name='ho-RP-' + str(i), createStepName=step_ho,
                                       variables=ho_var,
                                       region=ass.sets['RP-' + str(i)])
    
    # create loads on rigid bodies
    if if_prestep == 0 and is_expl:
        for i, load_val in rigid_loads.items():
            # create rigid part and instance
            inst_rigid = make_rigid_part(model, i, load_val, mesh_size)
            
            if is_expl:
                if 'v0' in load_val.keys():
                    # initial velocity
                    model.Velocity(name='v-ini-'+str(i), region=inst_rigid.sets['RP'],
                                   distributionType=MAGNITUDE, velocity1=load_val['v0'][0], velocity2=load_val['v0'][1],
                                   velocity3=load_val['v0'][2], omega=0)
                if 'u' in load_val.keys():
                    # displacement is applied via constant velocity or acceleration
                    if if_acc:
                        model.AccelerationBC(name='load-rigid-' + str(i), createStepName='load',
                                             region=inst_rigid.sets['RP'], a1=2 * load_val['u'][0] / t_step ** 2,
                                             a2=2 * load_val['u'][1] / t_step ** 2,
                                             a3=2 * load_val['u'][2] / t_step ** 2, ar1=0, ar2=0, ar3=0)
                    else:
                        model.VelocityBC(name='load-rigid-' + str(i), createStepName='load', vr1=0, vr2=0, vr3=0,
                                         region=inst_rigid.sets['RP'], v1=load_val['u'][0] / t_step,
                                         v2=load_val['u'][1] / t_step,
                                         v3=load_val['u'][2] / t_step)
            else:
                # direct displacement in implicit load step (however, does not work because 
                # of general contact problem)
                model.DisplacementBC(name='load-rigid-' + str(i), createStepName='load',
                                     region=inst_rigid.sets['RP'], u1=load_val['u'][0], u2=load_val['u'][1],
                                     u3=load_val['u'][2], ur1=0, ur2=0, ur3=0)
            
            # write history output of RPs into the odb file
            ho_vars = ('U1','U2','U3','V1','V2','V3','RF1','RF2','RF3')
            if is_expl:
                model.HistoryOutputRequest(name='ho-RP-' + str(i), createStepName='load',
                                           variables=ho_vars, numIntervals=500,
                                           region=inst_rigid.sets['RP'])
            else:
                model.HistoryOutputRequest(name='ho-RP-' + str(i), createStepName='load',
                                           variables=ho_vars,
                                           region=inst_rigid.sets['RP'])
    return


# ABAQUS main functions -------------------------------------------------------------

def make_abq_brick(model, lego_geom, nz, nx, b_type, mesh_size, is_tet=1):
    """Create the part(s) for one Lego brick. If tetrahedral elements are used (is_tet==1), only one part p_top is created, for hexahedral elements, the two parts p_top, p_bot are created.

    Args:
        model (Abaqus model)
        lego_geom (dict): Dictionary containing the general Lego dimensions.
        nz (int): Number of studs in x direction
        nx (int): Number of studs in z direction
        b_type (str): Brick type, can be 'regular', 'plate', 'tile', or 'base-plate'. The height of the 'brick' without the stud is 3 times higher than for the 'plate' and the 'tile'.
        mesh_size (float): The global mesh size to use for the part(s)
        is_tet (int, optional): What types of elements to use. is_tet==1: linear tetrahedral elements, is_tet == 0: linear hexahedral elements. Defaults to 1.

    Returns:
        p_bot (Abaqus part): The bottom part of the Lego brick, if necessary. Only != None for is_tet==0 and/or b_type == 'base_plate'.
        p_top (Abaqus part): The top part of the Lego brick.
        move_nodes (list): List of (node label, radial displacement, negative angle) elements
    """
    # read all parameters from the dictionaries
    b, b_gap, b_wall = lego_geom['b, b_gap, b_wall']
    h, h_stud, h_top = lego_geom['h, h_stud, h_top']
    r_in_small = lego_geom['inside small']['r']
    r_in_big, t_in_big = lego_geom['inside big']['r'], lego_geom['inside big']['t']
    # cavity wall properties
    t_rib_small, h_rib_small = lego_geom['inside small']['t_rib'], lego_geom['inside small']['h_rib']
    t_rib_big, h_rib_big = lego_geom['inside big']['t_rib'], lego_geom['inside big']['h_rib']
    
    delta_r = lego_geom['delta_r']
    #
    r_stud = b/2. - b_gap - b_wall + delta_r
    #
    brick_str = b_type + str(nz) + 'x' + str(nx)

    # Origin is at 1,1 stud
    x_arr = np.linspace(0, (nx - 1) * b, nx)
    z_arr = np.linspace(0, (nz - 1) * b, nz)

    # get outer coordinates of the brick (x-z plane)
    out_coord = ((-b / 2. + b_gap, -b / 2. + b_gap), (x_arr[-1] + b / 2. - b_gap, z_arr[-1] + b / 2. - b_gap))
    
    # get the mid positions between the studs
    x_mid_arr, z_mid_arr = make_middle_pos(x_arr, z_arr)
    
    # top part: find out about the height and if more than one part is needed
    if_full_part = 1
    # create only one part if tet or no delta_r
    if b_type == 'tile':
        # for tiles: always in one part
        h_total = h/3. + h_stud
    elif b_type == 'base-plate' or (is_tet == 0 and delta_r != 0.):
        # for (hex part where delta_r > 0) and baseplate: top part
        h_total = h_top + h_stud
        if_full_part = 0
    elif b_type == 'plate':
        # for plate & tet: one part with 1/3 of h
        h_total = h/3. + h_stud
    else:
        h_total = h + h_stud
    
    print(b_type, ' if_full_part:', if_full_part, ' h_total=', h_total)
    
    # create block for top part
    s_top = model.ConstrainedSketch(name='top-' + brick_str, sheetSize=200.0)
    s_top.rectangle(point1=out_coord[0], point2=out_coord[1])
    p_top = model.Part(name=brick_str + '-TOP', dimensionality=THREE_D,
                       type=DEFORMABLE_BODY)
    p_top.BaseSolidExtrude(depth=h_total, sketch=s_top)

    # cut from top to get studs
    t_cut = p_top.MakeSketchTransform(sketchPlane=p_top.faces.findAt(coordinates=(0, 0, h_total)),
                                      sketchUpEdge=p_top.edges.findAt(
                                          coordinates=(out_coord[1][0], 0, h_total)),
                                      sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                                      origin=(0, 0, 0))
    
    s_cut = model.ConstrainedSketch(name='top-cut-' + brick_str, sheetSize=200.0, transform=t_cut)
    s_cut.rectangle(point1=out_coord[0], point2=out_coord[1])
    
    # if tile, do not use circles for cut from top
    if b_type != 'tile':
        for x in x_arr:
            for y in z_arr:
                s_cut.CircleByCenterPerimeter(center=(x, y), point1=(x + r_stud, y))
    
    p_top.CutExtrude(sketchPlane=p_top.faces.findAt(coordinates=(0, 0, h_total)),
                     sketchUpEdge=p_top.edges.findAt(coordinates=(out_coord[1][0], 0, h_total)),
                     sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s_cut, depth=h_stud,
                     flipExtrudeDirection=OFF)
    
    p_top.Set(name='cont-top-touch', faces=p_top.faces.getByBoundingBox(zMin=h_total - h_stud - TOL, xMin=-r_stud - TOL,
                                                                        xMax=out_coord[1][0] - TOL))
    p_top.Surface(name='contact-top', side1Faces=p_top.faces.getByBoundingBox(zMin=h_total - h_stud - TOL))

    # partition part at bottom of studs
    if b_type != 'tile':
        p_top.PartitionCellByPlanePointNormal(point=p_top.vertices.findAt(coordinates=(r_stud, 0, h_total - h_stud)),
                                              normal=p_top.edges.findAt(
                                                  coordinates=(out_coord[1][0], out_coord[1][1], h_top / 2. + TOL)),
                                              cells=p_top.cells)

    # sets and surfaces for top part
    p_top.Set(name='all', cells=p_top.cells[:])
    p_top.Set(name='BOTTOM', faces=p_top.faces.getByBoundingBox(zMax=TOL))
    p_top.Set(name='TOP-FACES', faces=p_top.faces.getByBoundingBox(zMin=h_stud + h_top - TOL))

    # create side surfaces of top part
    p_top.Surface(name='x0', side1Faces=p_top.faces.getByBoundingBox(xMax=out_coord[0][0]+TOL))
    p_top.Surface(name='x1', side1Faces=p_top.faces.getByBoundingBox(xMin=out_coord[1][0] - TOL))
    p_top.Surface(name='z0', side1Faces=p_top.faces.getByBoundingBox(yMax=out_coord[0][1] + TOL))
    p_top.Surface(name='z1', side1Faces=p_top.faces.getByBoundingBox(yMin=out_coord[1][1] - TOL))

    # if there is also a bottom part, create tie surfaces
    if not if_full_part:
        p_top.Set(name='tie', faces=p_top.faces.getByBoundingBox(zMax=TOL))
        p_top.Surface(name='tie', side1Faces=p_top.faces.getByBoundingBox(zMax=TOL))
    #
    # create sets STUD-ij for loads_rp: i: ix and i: iz of the stud
    if b_type != 'tile':
        for i_z in range(1, nz + 1):
            for i_x in range(1, nx + 1):
                p_top.Set(name='STUD-' + str(i_x) + str(i_z),
                          faces=p_top.faces.getByBoundingBox(zMin=h_total - TOL - TOL,
                                                             yMin=-b / 2 + b * (i_z - 1), yMax=-b / 2 + b * i_z,
                                                             xMin=-b / 2 + b * (i_x - 1), xMax=-b / 2 + b * i_x))
    
    # create sketch for cut from bottom
    s_bot = model.ConstrainedSketch(name='bot-' + brick_str, sheetSize=200.0)
    s_bot.rectangle(point1=np.array(out_coord[0]) + b_wall, point2=np.array(out_coord[1]) - b_wall)

    # so that this sketch exists also for the 1x1 tiles
    model.ConstrainedSketch(name='bot-' + brick_str+'-wo-walls', objectToCopy=s_bot)
    
    # draw the inner walls, here!
    # using t_rib_small, h_rib_small, t_rib_big, h_rib_big
    # -------------------------------------------------------------
    
    def trim_curves(s,point_list):
        """Trim a curve at positions point_list
        """
        for point in point_list:
            s.autoTrimCurve(curve1=s.geometry.findAt(point), point1=point)
        return

    # the inner cylinders
    if type(z_mid_arr) != np.ndarray:
        if type(x_mid_arr) == np.ndarray:

            # positions for the inner walls
            y0 = out_coord[0][1]+b_wall
            y1 = (r_in_small**2-t_rib_small**2/4)**0.5

            # nx1 brick: just the small cylinders: directly in the sketch
            for i_x,x in enumerate(x_mid_arr):
                s_bot.CircleByCenterPerimeter(center=(x, 0), point1=(x + r_in_small, 0))
            
            # save sketch seperately: for cut from bottom
            model.ConstrainedSketch(name='bot-' + brick_str+'-wo-walls', objectToCopy=s_bot)
            
            # only put bottom walls when there is an even number of studs in this direction
            if b_type == 'regular' and nx%2 == 0:
                for i_x,x in enumerate(x_mid_arr):
                    # draw the inner walls if their thickness != 0
                    if (i_x-1)%2 == 0 and t_rib_small != 0.:
                        s_bot.Line(point1=(x-t_rib_small/2,y0), point2=(x-t_rib_small/2,-y1))
                        s_bot.Line(point1=(x+t_rib_small/2,y0), point2=(x+t_rib_small/2,-y1))
                        s_bot.Line(point1=(x-t_rib_small/2,-y0), point2=(x-t_rib_small/2,y1))
                        s_bot.Line(point1=(x+t_rib_small/2,-y0), point2=(x+t_rib_small/2,y1))
                        
                        # trim the curves of the outer line and the inner circle
                        trim_curves(s_bot,((x,y0),(x,-y0),(x,-r_in_small),(x,r_in_small)))
    else:
        # 1xn brick: just the small cylinders: directly in the sketch
        if type(x_mid_arr) != np.ndarray:
            
            for i_y,y in enumerate(z_mid_arr):
                # second point of circle up because of later trim curves
                s_bot.CircleByCenterPerimeter(center=(0, y), point1=(0, y + r_in_small))
            
            # save sketch seperately: for cut from bottom
            model.ConstrainedSketch(name='bot-' + brick_str+'-wo-walls', objectToCopy=s_bot)

            # positions for the inner walls
            x0 = out_coord[0][0]+b_wall
            x1 = (r_in_small**2-t_rib_small**2/4)**0.5
            
            # only put bottom walls when there is an even number of studs in this direction
            if b_type == 'regular' and nz%2 == 0:
                for i_y,y in enumerate(z_mid_arr):
                    # draw the inner walls if their thickness != 0
                    if (i_y-1)%2 == 0 and t_rib_small != 0.:
                        s_bot.Line(point1=(x0,y-t_rib_small/2), point2=(-x1,y-t_rib_small/2))
                        s_bot.Line(point1=(x0,y+t_rib_small/2), point2=(-x1,y+t_rib_small/2))
                        s_bot.Line(point1=(-x0,y-t_rib_small/2), point2=(x1,y-t_rib_small/2))
                        s_bot.Line(point1=(-x0,y+t_rib_small/2), point2=(x1,y+t_rib_small/2))
                        
                        # trim the curves of the outer line and the inner circle
                        trim_curves(s_bot,((x0,y),(-x0,y),(-r_in_small,y),(r_in_small,y)))
        else:
            # nx>1 and nz>1: bigger part: copy the sketch first,because 
            for i_x,x in enumerate(x_mid_arr):
                for i_y,y in enumerate(z_mid_arr):
                    # draw the inner tubes                    
                    s_bot.CircleByCenterPerimeter(center=(x, y), point1=(x + r_in_big/2**0.5, y + r_in_big/2**0.5))
                    s_bot.CircleByCenterPerimeter(center=(x, y), point1=(x + r_in_big - t_in_big, y))
            
            # save sketch seperately: for cut from bottom
            model.ConstrainedSketch(name='bot-' + brick_str+'-wo-walls', objectToCopy=s_bot)
            
            if b_type == 'regular':
                for i_x,x in enumerate(x_mid_arr):
                    for i_y,y in enumerate(z_mid_arr):
                        # outer coordinates
                        c0 = out_coord[0][0]+b_wall
                        c1 = (r_in_big**2-t_rib_big**2/4)**0.5
                        xmax = out_coord[1][0]-b_wall
                        ymax = out_coord[1][1]-b_wall        
                        # check if the wall goes to the outer wall or to the next tube

                        # only put bottom walls when there is an even number of studs in this direction
                        if nx%2 == 0:
                            # every second tube in x direction
                            if (i_x-1)%2 == 0 and t_rib_big != 0.:
                                # if bottom tube: wall to tube
                                if i_y == 0:
                                    s_bot.Line(point1=(x-t_rib_big/2,c0), point2=(x-t_rib_big/2,y-c1))
                                    s_bot.Line(point1=(x+t_rib_big/2,c0), point2=(x+t_rib_big/2,y-c1))
                                    trim_curves(s_bot,((x,c0),(x,y-r_in_big)))
                                # if top tube: tube to wall
                                if i_y+1 == len(z_mid_arr):
                                    s_bot.Line(point1=(x-t_rib_big/2,y+c1), point2=(x-t_rib_big/2,ymax))
                                    s_bot.Line(point1=(x+t_rib_big/2,y+c1), point2=(x+t_rib_big/2,ymax))
                                    
                                    trim_curves(s_bot,((x,y+r_in_big),(x,ymax)))
                                else:
                                    # if connection in the middle: tube-tube
                                    s_bot.Line(point1=(x-t_rib_big/2,y+c1), point2=(x-t_rib_big/2,z_mid_arr[i_y+1]-c1))
                                    s_bot.Line(point1=(x+t_rib_big/2,y+c1), point2=(x+t_rib_big/2,z_mid_arr[i_y+1]-c1))
                                    trim_curves(s_bot,((x,y+r_in_big),(x,z_mid_arr[i_y+1]-r_in_big)))
                        
                        # only put bottom walls when there is an even number of studs in this direction
                        if nz%2 == 0:
                            # every second tube in y direction
                            if (i_y-1)%2 == 0 and t_rib_big != 0.:
                                # if outer left tube: wall to tube
                                if i_x == 0:
                                    s_bot.Line(point1=(c0,y-t_rib_big/2), point2=(x-c1,y-t_rib_big/2))
                                    s_bot.Line(point1=(c0,y+t_rib_big/2), point2=(x-c1,y+t_rib_big/2))
                                    trim_curves(s_bot,((c0,y),(x-r_in_big,y)))
                                # if most right tube: tube to wall
                                if i_x+1 == len(x_mid_arr):
                                    s_bot.Line(point1=(x+c1,y-t_rib_big/2), point2=(xmax,y-t_rib_big/2))
                                    s_bot.Line(point1=(x+c1,y+t_rib_big/2), point2=(xmax,y+t_rib_big/2))
                                    trim_curves(s_bot,((xmax,y),(x+r_in_big,y)))
                                else:
                                    # if connection in the middle: tube-tube
                                    s_bot.Line(point1=(x+c1,y-t_rib_big/2), point2=(x_mid_arr[i_x+1]-c1,y-t_rib_big/2))
                                    s_bot.Line(point1=(x+c1,y+t_rib_big/2), point2=(x_mid_arr[i_x+1]-c1,y+t_rib_big/2))
                                    trim_curves(s_bot,((x+r_in_big,y),(x_mid_arr[i_x+1]-r_in_big,y)))
                            
    # for extruding bottom part: add the outer rectangle to sketch
    if not if_full_part:
        s_bot.rectangle(point1=out_coord[0], point2=out_coord[1])
        
    # use the sketch for cut from bottom if there is only one part
    if if_full_part:
        if b_type != 'base-plate':
            # make a sketch for the bottom face and load the cut from 
            # bottom sketch from above
            t = p_top.MakeSketchTransform(sketchPlane=p_top.faces.findAt(coordinates=(TOL, TOL, 0)),
                                      sketchUpEdge=p_top.edges.findAt(coordinates=(out_coord[1][0], -TOL, 0)),
                                      sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0,0,0))
            s_cut = model.ConstrainedSketch(name='bot-' + brick_str + '-cut', sheetSize=200.0, transform=t)
            
            # only leave the walls if 'brick'
            if b_type == 'regular':
                s_cut.retrieveSketch(sketch=s_bot)
            else:
                s_cut.retrieveSketch(sketch=model.sketches['bot-' + brick_str+'-wo-walls'])
            #
            g = s_cut.geometry
            hl = s_cut.ConstructionLine(point1=(0, 0), point2=(1, 0))
            s_cut.mirror(mirrorLine=hl, objectList=g.values())

            # cut from bottom
            p_top.CutExtrude(sketchPlane=p_top.faces.findAt(coordinates=(TOL, TOL, 0)),
                         sketchUpEdge=p_top.edges.findAt(coordinates=(out_coord[1][0], -TOL, 0)),
                         sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s_cut,
                         depth=h_total - h_stud - h_top, flipExtrudeDirection=OFF)

        # partition in z direction
        if b_type != 'base-plate':
            cut2 = p_top.DatumPlaneByPrincipalPlane(offset=h_total - h_stud - h_top, principalPlane=XYPLANE)
            p_top.PartitionCellByDatumPlane(cells=p_top.cells[:], datumPlane=p_top.datums[cut2.id])

        if b_type == 'regular':
            cut0 = p_top.DatumPlaneByPrincipalPlane(offset=h_stud + mesh_size, principalPlane=XYPLANE)
            p_top.PartitionCellByDatumPlane(cells=p_top.cells[:], datumPlane=p_top.datums[cut0.id])
    
    #print('part of type ' + b_type + ' created :-)')
    ass = model.rootAssembly

    # bottom part
    if not if_full_part and b_type != 'base-plate':
        p_bot = model.Part(name=brick_str + '-BOT', dimensionality=THREE_D,
                       type=DEFORMABLE_BODY)
        
        # just use the cut from bottom sketch from above
        if b_type == 'regular':
            # full height h
            p_bot.BaseSolidExtrude(depth=h - h_top, sketch=s_bot)
            p_bot.Set(name='tie', faces=p_bot.faces.getByBoundingBox(zMin=h - h_top - TOL))
            p_bot.Surface(name='tie', side1Faces=p_bot.faces.getByBoundingBox(zMin=h - h_top - TOL))
        else:
            # one one third of the height for 'plate'
            p_bot.BaseSolidExtrude(depth=h / 3. - h_top, sketch=s_bot)
            p_bot.Set(name='tie', faces=p_bot.faces.getByBoundingBox(zMin=h / 3. - h_top - TOL))
            p_bot.Surface(name='tie', side1Faces=p_bot.faces.getByBoundingBox(zMin=h / 3. - h_top - TOL))
    
    # only bricks can have walls in lowe cavities
    if b_type == 'regular':
        # cut from bottom, if needed
        if if_full_part:
            p_cut = p_top
        else:
            p_cut = p_bot

        if_cut = 0

        # cut here again if the height of inner walls h < h-h_top --------------------
        # check if 1xn or nx1 brick: then use h_rib_small
        if (nx == 1 and nz > 3) or (nz == 1 and nx > 3):
            h_rib = h_rib_small
            t_rib = t_rib_small
            if_cut = 1
        elif (nx > 1 and nz > 3) or (nz > 0 and nx > 3):
            h_rib = h_rib_big
            t_rib = t_rib_big
            if_cut = 1
        
        if if_cut and h_rib < h-h_top-TOL and t_rib != 0.:
            t = p_cut.MakeSketchTransform(sketchPlane=p_cut.faces.findAt(coordinates=(out_coord[0][0]+TOL, out_coord[0][1]+TOL, 0)),
                                    sketchUpEdge=p_cut.edges.findAt(coordinates=(out_coord[1][0], -TOL, 0)),
                                    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0,0,0))
            s_cut = model.ConstrainedSketch(name='bot-' + brick_str + '-cut2', sheetSize=200.0, transform=t)
            s_cut.retrieveSketch(sketch=model.sketches['bot-' + brick_str+'-wo-walls'])
            g = s_cut.geometry
            hl = s_cut.ConstructionLine(point1=(0, 0), point2=(1, 0))
            s_cut.mirror(mirrorLine=hl, objectList=g.values())

            # cut from bottom again
            p_cut.CutExtrude(sketchPlane=p_cut.faces.findAt(coordinates=(out_coord[0][0]+TOL, out_coord[0][1]+TOL, 0)),
                        sketchUpEdge=p_cut.edges.findAt(coordinates=(out_coord[1][0], -TOL, 0)),
                        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s_cut,
                        depth=h-h_top-h_rib, flipExtrudeDirection=OFF)
            
            # partition in this cut plane
            p0 = p_cut.DatumPlaneByPrincipalPlane(offset=h-h_top-h_rib, principalPlane=XYPLANE)
            p_cut.PartitionCellByDatumPlane(cells=p_cut.cells[:], datumPlane=p_cut.datums[p0.id])

    # mesh the top part
    p_top.seedPart(size=mesh_size)

    if is_tet:
        p_top.setMeshControls(regions=p_top.cells[:], elemShape=TET, technique=FREE)
    #
    p_top.generateMesh()
    if not is_tet:
        # there is usually some contact penetration for hex elements and gen. contact
        # in explicit model
        p_top.setElementType(regions=p_top.sets['all'],
                            elemTypes=(mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD),
                                        mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN,
                                        secondOrderAccuracy=OFF, hourglassControl=ENHANCED),
                                        mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)))
    else:
        p_top.setElementType(regions=p_top.sets['all'],elemTypes=(mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT,
                            secondOrderAccuracy=ON, distortionControl=ON,
                            lengthRatio=0.1),))
    
    if not if_full_part and b_type != 'base-plate':
        # Sets and Surfaces
        p_bot.Set(name='all', cells=p_bot.cells[:])
        p_bot.Set(name='BOTTOM', faces=p_bot.faces.getByBoundingBox(zMax=TOL))

        # partition the bottom part (only for brick)
        if b_type == 'regular':
            cut0 = p_bot.DatumPlaneByPrincipalPlane(offset=h_stud + mesh_size, principalPlane=XYPLANE)
            p_bot.PartitionCellByDatumPlane(cells=p_bot.cells[:], datumPlane=p_bot.datums[cut0.id])
        
        # mesh the bottom part
        p_bot.seedPart(size=mesh_size)
        p_bot.generateMesh()
        p_bot.setElementType(regions=p_bot.sets['all'], elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,
                                                                                 kinematicSplit=AVERAGE_STRAIN,
                                                                                 secondOrderAccuracy=OFF,
                                                                                 hourglassControl=ENHANCED),
                                                                   mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD),
                                                                   mesh.ElemType(elemCode=C3D4,
                                                                                 elemLibrary=STANDARD)))
        p_contact = p_bot

        # create side surfaces for bottom part
        p_bot.Surface(name='x0', side1Faces=p_bot.faces.getByBoundingBox(xMax=out_coord[0][0] + TOL))
        p_bot.Surface(name='x1', side1Faces=p_bot.faces.getByBoundingBox(xMin=out_coord[1][0] - TOL))
        p_bot.Surface(name='z0', side1Faces=p_bot.faces.getByBoundingBox(yMax=out_coord[0][1] + TOL))
        p_bot.Surface(name='z1', side1Faces=p_bot.faces.getByBoundingBox(yMin=out_coord[1][1] - TOL))
    else:
        # in only one part, the top one is used also for bottom contact
        p_bot = None
        p_contact = p_top
    
    # create contact surface for bottom cavities
    p_contact.Set(name='cont-bot-touch0', faces=p_contact.faces.getByBoundingBox(zMax=h_stud + mesh_size + TOL,
                                                                         xMin=out_coord[0][0] + b_wall - TOL,
                                                                         xMax=out_coord[1][0] - b_wall + TOL))
    p_contact.Set(name='cont-bot-neg-widen', faces=p_contact.faces.getByBoundingBox(zMax=TOL) +
                  p_contact.faces.getByBoundingBox(zMin=h/3.-h_top-TOL, zMax=h/3.-h_top+TOL) +
                  p_contact.faces.getByBoundingBox(zMin=h_stud + mesh_size - TOL, zMax=h_stud + mesh_size + TOL))

    p_contact.Set(name='cont-bot-neg-cont', faces=p_contact.faces.getByBoundingBox(zMax=TOL) +
                  p_contact.faces.getByBoundingBox(zMin=h_stud + mesh_size - TOL, zMax=h_stud + mesh_size + TOL))

    p_contact.SetByBoolean(name='cont-bot-touch', operation=DIFFERENCE,
                       sets=(p_contact.sets['cont-bot-touch0'], p_contact.sets['cont-bot-neg-widen'],))
    
    p_contact.SetByBoolean(name='cont-bot', operation=DIFFERENCE,
                       sets=(p_contact.sets['cont-bot-touch0'], p_contact.sets['cont-bot-neg-cont'],))
    
    del p_contact.sets['cont-bot-touch0']
    del p_contact.sets['cont-bot-neg-widen']
    del p_contact.sets['cont-bot-neg-cont']
    
    p_contact.Surface(name='contact-bot', side1Faces=p_contact.faces.getByBoundingBox(zMax=TOL)
                                                     + p_contact.sets['cont-bot'].faces)
    
    # calculate the node displacements using p,nx,ny,r_stud,h_stud
    if b_type != 'base-plate':
        if if_full_part:
            move_nodes = make_u_nodes(p_top, nx, nz, b, r_stud, h_stud, mesh_size)
        else:
            move_nodes = make_u_nodes(p_bot, nx, nz, b, r_stud, h_stud, mesh_size)
    else:
        # empty array for base-plate: No movements needed there
        move_nodes = np.array([])
    return (p_bot, p_top), move_nodes


def make_assembly(assembly, lego_geom, mesh_size, is_tet):
    """Create all needed parts, make instances of them in the assembly and move them to their desired position, create tie constraints between parts (only for hex mesh)

    Args:
        assembly (dict): The setup of the Lego set containing bricks, brick positions, boundary conditions, and loads on sets.
        lego_geom (dict): Basic Lego brick geometry parameters
        mesh_size (float): Global mesh size for all parts
        is_tet (int): If linear tetrahedral elements should be used (linear hexahedral elements, otherwise)

    Returns:
        (Abaqus model, dict): Dictionary based on the 'parts' dict of the assembly, but with added items 'instance':(i_top, i_bot)
    """
    # load parameters from dictionaries
    h, h_stud, h_top = lego_geom['h, h_stud, h_top']
    b, _, _ = lego_geom['b, b_gap, b_wall']
    
    brick_dict = assembly['bricks']
    pos_dict = assembly['parts']
    
    # initialize the model and change its name to model_name
    Mdb()
    mdb.models.changeKey(fromName='Model-1', toName=assembly['name'])
    model = mdb.models[assembly['name']]
    
    # create the brick parts, put the parts and their nodes to move
    # into dictionary `brick_dict`
    for i_brick, brick_i in brick_dict.items():
        #
        b_type, nx, nz = (brick_i['type'], brick_i['nx'], brick_i['nz'])
        dict_i = brick_dict[i_brick]
        #
        dict_i['part'], dict_i['move_nodes'] = make_abq_brick(model, lego_geom, nz, nx, b_type,
                                                              mesh_size, is_tet)
    # assign sections
    make_sections(model, lego_geom['E, nu, dens'])
    
    ass = model.rootAssembly
    
    # create all instances: one part can be used more than once
    # ('BRICK01', 'BRICK02', ...)
    for i_inst, inst_i in pos_dict.items():
        # number of brick used, and its position
        i_brick, loc = (inst_i['brick_id'], inst_i['loc'])
        
        # get brick parameters and part name from that
        brick_i = brick_dict[i_brick]
        b_type, nx, nz = (brick_i['type'], brick_i['nx'], brick_i['nz'])
        brick_str = b_type + str(nz) + 'x' + str(nx)

        # in `brick_dict`: (p_bot, p_top): if p_bot==None: only one part
        if brick_i['part'][0] is None:
            inst = ass.Instance(name='BRICK' + str(i_inst).zfill(2) + '-TOP',
                                part=model.parts[brick_str+'-TOP'], dependent=ON)
            
            # for use later in the interactions: top and bottom instance the same
            pos_dict[i_inst]['instance'] = (inst, inst)

            # so that brick lies in the x-z plane and its upper left stud lies in the origin
            ass.translate(instanceList=(inst.name,), vector=(0, -b*(nz-1), 0))
            ass.rotate(instanceList=(inst.name,), axisPoint=(0, 0, 0),
                       axisDirection=(1, 0, 0), angle=-90) 

            # translate instance to its position given by the user
            ass.translate(instanceList=(inst.name,), vector=tuple(loc))
            if b_type == 'base-plate':
                ass.translate(instanceList=(inst.name,), vector=(0, -h_top, 0))
        else:
            # here there are actually two parts per brick (except for base-plates)
            i_top = ass.Instance(name='BRICK' + str(i_inst).zfill(2) + '-TOP',
                               part=model.parts[brick_str + '-TOP'], dependent=ON)

            # so that brick lies in the x-z plane and its upper left stud lies in the origin
            ass.translate(instanceList=(i_top.name,), vector=(0, -b * (nz-1), 0))
            ass.rotate(instanceList=(i_top.name,), axisPoint=(0, 0, 0),
                       axisDirection=(1, 0, 0), angle=-90)
            
            # move in y direction such that bottom is at y == 0 (except for base-plate)
            ass.translate(instanceList=(i_top.name,), vector=tuple(loc))
            if b_type == 'regular':
                ass.translate(instanceList=(i_top.name,), vector=(0, h - h_top, 0))
            elif b_type != 'base-plate':
                ass.translate(instanceList=(i_top.name,), vector=(0, h / 3. - h_top, 0))
            
            # base palate can have hex. elements, but does not need bottom part
            if b_type != 'base-plate':
                i_bot = ass.Instance(name='BRICK' + str(i_inst).zfill(2) + '-BOT',
                                   part=model.parts[brick_str + '-BOT'], dependent=ON)
                
                # translate bottom instance to its position
                ass.translate(instanceList=(i_bot.name,), vector=(0, -b * (nz-1), 0))
                ass.rotate(instanceList=(i_bot.name,), axisPoint=(0, 0, 0),
                           axisDirection=(1, 0, 0), angle=-90)
                ass.translate(instanceList=(i_bot.name,), vector=tuple(loc))
                pos_dict[i_inst]['instance'] = (i_bot, i_top)

                # tie-constraints between top and bottom instance (for hex. elements)
                model.Tie(name='tie-brick-' + str(i_inst), master=i_bot.surfaces['tie'],
                          slave=i_top.surfaces['tie'], positionToleranceMethod=COMPUTED,
                          adjust=OFF, tieRotations=OFF, thickness=ON)
            else:
                # for base-plate: top and bottom instance is the same
                pos_dict[i_inst]['instance'] = (i_top, i_top)
    return model, pos_dict


def make_model_load(model_name0, assembly, explicit_par, t_step, if_prestep=0, nf_expl=80):
    """Create either the explicit model that loads the Lego set and runs & evaluates the model (for load_par['is_expl']==1) or creates an additional implicit load step (for load_par['is_expl']==0).

    Args:
        model_name0 (str): Name of the Abaqus model
        assembly (dict): The setup of the Lego set containing bricks, brick positions, boundary conditions, and loads on sets.
        explicit_par (dict): Load parameters of the model
        t_step (float): Time duration of the step.
        if_prestep (int, optional): If the function is called for creating the loads in the clamping steps (`widen`, `contact`, `free`) or from the load step. Because for an explicit load, the reference points of the loads_rp should be fixed in the `free` step to avoid rigid body motion. Defaults to 0.
        nf_expl (int, optional): Number of field output frames in the explicit model. Defaults to 80.
    """
    # load parameters from the input dictionaries
    mu = assembly['mu']
    mesh_size = assembly['mesh_size']

    if explicit_par == {}:
        is_expl = 0
        is_tet = 0
        load_str = ''
    else:
        is_expl = 1
        is_tet = 1
        ms_time = explicit_par['mass_scale_t']
        load_str = explicit_par['load_str']

    model_name = model_name0 + load_str

    # if prestep of explicit load: fix the reference points that will be later loaded in the 'free' step
    # (for implicit load step, this is done separately): do not run model yet
    if is_expl and if_prestep:
        make_loads(mdb.models[model_name0], assembly, explicit_par, mesh_size, 1, is_expl=1, if_prestep=1)
        return

    if is_expl:
        # use the time duration of the step in the model name
        model_name += '-'+str(round(t_step*1000,3)).zfill(5).replace('.', '_')+'ms'
        model = mdb.Model(name=model_name, modelType=STANDARD_EXPLICIT)
        a = model.rootAssembly

        # load the mesh of the deformed parts into the model & make general contact
        inst_names = load_impl_to_expl(model, model_name0)
        make_contact(model, mu, is_expl=1)

        # create the loads on reference points that are connected to sets
        make_load_rps(model, assembly['loads_rp'])

        # load initial state (deformations and stresses)
        model.InitialState(updateReferenceConfiguration=OFF,
                           fileName=model_name0, endStep=LAST_STEP, endIncrement=STEP_END,
                           name='from-impl', createStepName='Initial',
                           instances=tuple(a.instances[i] for i in inst_names))
        
        step = model.ExplicitDynamicsStep(name='load', previous='Initial',
                                          timePeriod=t_step, improvedDtMethod=ON)
        
        # request field output in `nf_expl` frames
        model.fieldOutputRequests['F-Output-1'].setValues(numIntervals=nf_expl)
        
        # set mass scaling to reduce copmutational cost
        if ms_time != 0:
            step.setValues(massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 0.0,
                                         ms_time, BELOW_MIN, 0, 0, 0.0, 0.0, 0, None),),
                           improvedDtMethod=ON)
        
        # tie the top and bottom bricks together
        for inst_bot_name in [i for i in inst_names if 'BOT' in i]:

            # connect top and bottom part using tie constraints, if necessary
            model.Tie(name='tie-part-' + str(inst_bot_name[:-4]), master=a.instances[inst_bot_name].sets['TIE'],
                  slave=a.instances[inst_bot_name.replace('BOT', 'TOP')].sets['TIE'],
                  positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=OFF, thickness=ON)
        
        # create the loads in the explicit model: RP of rigid bodies
        make_loads(model, assembly, explicit_par, mesh_size, t_step, is_expl=1)
        
        # apply the boundary conditions `bc`. Try which of the both parts contains the set
        # that should be coupled to a RP
        for i, fix_par in assembly['bc'].items():
            if is_tet:
                inst_i = a.instances['BRICK' + str(fix_par['part_id']).zfill(2)+'-TOP']
            else:
                try:
                    inst_i = a.instances['BRICK' + str(fix_par['part_id']).zfill(2) + '-BOT']
                except:
                    # when base-plate
                    inst_i = a.instances['BRICK' + str(fix_par['part_id']).zfill(2) + '-TOP']
            model.DisplacementBC(name='fix-part-' + str(i).zfill(2), createStepName='Initial',
                                 region=inst_i.sets[fix_par['set_name'].upper()], u1=0, u2=0, u3=0)
        
        # for debugging: write the load parameters and step time to a json file
        with open(model_name + '-par.json', 'w') as f:
            json.dump({'explicit_par': explicit_par, 't_step': t_step}, f)
        
        # run & evaluate the explcit model
        run_model(model, model_name)
        get_ho(model_name)
        create_video(model_name)
    else:
        # implicit model: create the load step after the free step, apply loads and run
        model = mdb.models.values()[0]
        model.StaticStep(name='load', previous='free', maxNumInc=1000,
                         initialInc=0.05, minInc=1e-5, maxInc=0.05, nlgeom=ON)
        make_loads(model, assembly, explicit_par, mesh_size, 1, is_expl=0, if_prestep=0)
    
        # for debugging: write the load parameters and step time to a json file
        with open(model_name + '-par.json', 'w') as f:
            json.dump({'explicit_par': explicit_par, 't_step': t_step}, f)
    return


def make_model_ini(model_name, assembly, lego_geom, explicit_par, mesh_size, mu, is_expl, is_tet):
    """Create and run the initial model with the clamping steps. For an implicit load step, the load step is added and the whole model is run and evaluated. For an explicit load step, only the clamping steps are run and the results of that are later loaded into the explicit model.

    Args:
        model_name (Abaqus Model)
        assembly (dict): The setup of the Lego set containing bricks, brick positions, boundary conditions, and loads on sets.
        lego_geom (dict): Basic Lego brick geometry parameters
        explicit_par (dict): Load parameters of the model
        mesh_size (float): Global mesh size for all parts
        mu (float): Coefficient of friction between the Lego bricks
        is_expl (int): If the load should be applied in an explicit step. For `is_expl==1`, the step `free` is the last step and a restart request is written for loading the results into the explicit model. For `is_expl==0`, the step `load` is added after the step `free` and loads are created.
        is_tet (int): If linear tetrahedral elements should be used (linear hexahedral elements, otherwise)
    """
    # create the assembly
    model, pos_dict = make_assembly(assembly, lego_geom, mesh_size, is_tet)

    # obtain contact faces between bricks and cavities that should be widened
    widen_cav, cont_pairs = make_widen_cont(assembly, lego_geom, if_plot=0)
    
    # Steps: widen cavities, establish contact with holding the bricks,
    #        then let them free
    model.StaticStep(name='widen', previous='Initial', maxNumInc=1000,
                     initialInc=1, minInc=1e-08, maxInc=1, nlgeom=ON)
    model.StaticStep(name='contact', previous='widen', maxNumInc=1000,
                     initialInc=1, minInc=1e-08, maxInc=1, nlgeom=ON)
    
    # free step
    free_step = model.StaticStep(name='free', previous='contact', maxNumInc=1000,
                                 initialInc=1, minInc=1e-08, maxInc=1, nlgeom=ON)
    if is_expl:
        # write restart request so that the explicit model can restart the analysis
        free_step.Restart(frequency=0, numberIntervals=1, overlay=OFF, timeMarks=OFF)

    # create the boundary conditions `bc`
    for i, fix_par in assembly['bc'].items():
        inst_i_bot, inst_i_top = pos_dict[fix_par['part_id']]['instance']
        try:
            model.DisplacementBC(name='fix-brick-' + str(i).zfill(2), createStepName='Initial',
                                 region=inst_i_bot.sets[fix_par['set_name']], u1=0, u2=0, u3=0)
        except:
            model.DisplacementBC(name='fix-brick-' + str(i).zfill(2), createStepName='Initial',
                                 region=inst_i_top.sets[fix_par['set_name']], u1=0, u2=0, u3=0)

    # widen the cavities so the studs fit in there
    for i_inst in pos_dict.keys():
        brick_i = assembly['bricks'][pos_dict[i_inst]['brick_id']]
        inst_i_bot, inst_i_top = pos_dict[i_inst]['instance']

        # only widen bricks where at least one stud is sticking (from `widen_cav`)
        if len(widen_cav[i_inst]) > 0:
            move_nodes = brick_i['move_nodes']
            for n_label_i, dr_i, ang_i in move_nodes:
                n_label = int(n_label_i)
                # get cartesian coordinates from angular and radial displacement
                ux, uz = dr_i * np.cos(ang_i * np.pi / 180.), -dr_i * np.sin(ang_i * np.pi / 180.)
                bc_temp = model.DisplacementBC(name='u0-inst' + str(i_inst) + 'n' + str(n_label),
                                               createStepName='widen', region=inst_i_bot.sets['x-n' + str(n_label)],
                                               u1=ux, u2=0, u3=uz)
                bc_temp.deactivate('contact')
        
        # fix the top faces while establishing contact (discontinue in `free` step)
        top_temp = model.DisplacementBC(name='fix-top-i' + str(i_inst), createStepName='Initial',
                                        region=inst_i_top.sets['TOP-FACES'], u1=0, u2=0, u3=0)
        top_temp.deactivate('free')

    # print the assembly for debugging
    print_assembly(model, model_name)

    # define contact & reference points for loads
    make_contact(model, mu, pos_dict, cont_pairs, is_expl=0)
    make_load_rps(model, assembly['loads_rp'])
    
    # implicit model: already add the load step
    make_model_load(model_name, assembly, explicit_par, 1, if_prestep=1)

    # run the initial model (expl.) or full model (impl.)
    run_model(model, model_name)

    # is_expl==0: implicit model already run and evaluated in function make_model_load
    if not is_expl:
        # evaluate the model, if implicit: because load has been already applied
        get_ho(model_name)
        create_video(model_name)
    return


def make_model(assembly0, explicit_par={}, lego_geom=lego_geom, is_new=1, n_frames_expl=80):
    """Creates, runs, and evaluates Lego model in Abaqus. Parameters should be given in the N-mm-s system.

    Args:
        assembly0 (dict): The setup of the Lego set defining bricks used (sub-dictionary `bricks`), brick positions (`pos`), boundary conditions (`bc`), and loads applied on sets (`loads_rp`)
        explicit_par (dict): Dictionary that defines the load parameters in the model
        lego_geom (dict): Dictionary containing the general Lego dimensions and elastic parameters.
        is_new (int, optional): If the implicit calculation (steps `widen`, `contact`, `free`) should be performed anew or already exists. Only relevant for an explicit load step. Defaults to 1.
        n_frames_expl (int, optional): Number of output frames in the explicit model. Defaults to 80.
    
    Examples:
        assembly0: {'name':'case0-widen', 
                       'bricks':{1:{'type':'base-plate', 'nx':1, 'nz':1}, 2:{'type':'brick', 'nx':1, 'nz':1}},
                       'parts':{1:{'brick_id':1, 'loc':(0,0,0)},
                                2:{'brick_id':2, 'loc':(0,0,0)}},
                       'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
                       'loads_rp':{1:{'part_id':2, 'set_name':'STUD-11', 'uy':2}}}
        explicit_par: {'mass_scale_t': 0, 't_step': 0.0005, 'is_acc': 0,
                          'load_str': '-plate', 'loads_rigid': {}}
        lego_geom: {'b, b_gap, b_wall': (8., 0.1, 1.6),
                              'h, h_stud, h_top': (9.6, 1.7, 1.5),
                              'inside small': {'r': 1.6, 't_rib':1.1, 'h_rib':8.1},
                              'inside big': {'r': 3.3, 't': 0.9, 't_rib':0.9, 'h_rib': 6.8},
                              'delta_r': 0.05,
                              'E, nu, dens': (2200., 0.35, 1e-9)}
    """
    remove_files(DIR0)
    model_name = assembly0['name']

    mesh_size = assembly0['mesh_size']
    mu = assembly0['mu']

    if explicit_par == {}:
        is_expl = 0
        is_tet = 0
    else:
        is_expl = 1
        is_tet = 1
        t_steps = explicit_par['t_step']

    # information will be written into assembly, assembly0 should stay unchanged
    assembly = copy.deepcopy(assembly0)

    solver_type = {1: 'expl', 0:'impl'}

    # create the name for the model directory
    run_dir = model_name + '-' + solver_type[is_expl] + '-mesh' + str(int(mesh_size * 100)).zfill(3)+'mm'

    # change into that directory
    make_dir(run_dir, if_clear=is_new, if_change=1)

    # write input into json file
    with open('_dict-assembly.json', 'w') as f:
        json.dump({'lego_geom': lego_geom, 'assembly': assembly}, f)

    if is_new:
        # run the initial model for establishing the clamping, implicit: already also run the load step
        make_model_ini(model_name, assembly, lego_geom, explicit_par, mesh_size, mu, is_expl, is_tet)

    # load the deformed and stressed Lego bricks from the implicit clamping model to the explicit model,
    # apply the loads, and run & evaluate the explicit model
    if is_expl:
        if type(t_steps) != list:
            # t_steps is a scalar: only call once
            make_model_load(model_name, assembly, explicit_par, t_steps, nf_expl=n_frames_expl)
        else:
            # t_steps is a list: repeatedly run explicit models with varying time step length
            for t_i in t_steps:
                make_model_load(model_name, assembly, explicit_par, t_i, nf_expl=n_frames_expl)
    #
    remove_files(DIR0+'/'+run_dir)
    os.chdir(DIR0)
    return
