# Pumpkin Lego model, halloween 2024
# -----------------------------------------------------------------------

from brickfem import make_model, make_dir, create_video

dy, dx, dz = 9.6, 8., 8.

assembly_pumpkin = {
                'name':'pumpkin',
                'bricks':{
                        1:{'type':'regular', 'nx':1, 'nz':2},
                        2:{'type':'regular', 'nx':2, 'nz':2},
                        3:{'type':'regular', 'nx':3, 'nz':2},
                        4:{'type':'regular', 'nx':4, 'nz':2},
                        5:{'type':'base-plate', 'nx':10, 'nz':2},
                        },
                'parts':{
                        # base plate
                        1:{'brick_id':5, 'loc':(-dx, 0, 0), 'c':'Grey'},
                        
                        2:{'brick_id':4, 'loc':(0, 0, 0), 'c':'Orange'},

                        3:{'brick_id':4, 'loc':(4*dx, 0, 0), 'c':'Orange'},
                        4:{'brick_id':3, 'loc':(-dx, dy, 0), 'c':'Orange'},
                        5:{'brick_id':4, 'loc':(2*dx, dy, 0), 'c':'Black'},
                        6:{'brick_id':3, 'loc':(6*dx, dy, 0), 'c':'Orange'},
                        7:{'brick_id':3, 'loc':(-2*dx, 2*dy, 0), 'c':'Orange'},

                        8:{'brick_id':1, 'loc':(dx, 2*dy, 0), 'c':'Black'},
                        9:{'brick_id':2, 'loc':(2*dx, 2*dy, 0), 'c':'Orange'},
                        10:{'brick_id':2, 'loc':(4*dx, 2*dy, 0), 'c':'Orange'},
                        11:{'brick_id':1, 'loc':(6*dx, 2*dy, 0), 'c':'Black'},
                        12:{'brick_id':3, 'loc':(7*dx, 2*dy, 0), 'c':'Orange'},

                        13:{'brick_id':1, 'loc':(-2*dx, 3*dy, 0), 'c':'Orange'},
                        14:{'brick_id':4, 'loc':(-dx, 3*dy, 0), 'c':'Orange'},
                        15:{'brick_id':2, 'loc':(3*dx, 3*dy, 0), 'c':'Black'},
                        16:{'brick_id':3, 'loc':(5*dx, 3*dy, 0), 'c':'Orange'},
                        17:{'brick_id':2, 'loc':(8*dx, 3*dy, 0), 'c':'Orange'},

                        18:{'brick_id':2, 'loc':(-2*dx, 4*dy, 0), 'c':'Orange'},
                        19:{'brick_id':2, 'loc':(0, 4*dy, 0), 'c':'Black'},
                        20:{'brick_id':4, 'loc':(2*dx, 4*dy, 0), 'c':'Orange'},
                        21:{'brick_id':2, 'loc':(6*dx, 4*dy, 0), 'c':'Black'},
                        22:{'brick_id':2, 'loc':(8*dx, 4*dy, 0), 'c':'Orange'},

                        23:{'brick_id':2, 'loc':(-dx, 5*dy, 0), 'c':'Orange'},
                        24:{'brick_id':2, 'loc':(dx, 5*dy, 0), 'c':'Orange'},
                        25:{'brick_id':2, 'loc':(3*dx, 5*dy, 0), 'c':'Orange'},
                        26:{'brick_id':2, 'loc':(5*dx, 5*dy, 0), 'c':'Orange'},
                        27:{'brick_id':2, 'loc':(7*dx, 5*dy, 0), 'c':'Orange'},

                        28:{'brick_id':2, 'loc':(0, 6*dy, 0), 'c':'Orange'},
                        29:{'brick_id':2, 'loc':(2*dx, 6*dy, 0), 'c':'Orange'},
                        30:{'brick_id':2, 'loc':(4*dx, 6*dy, 0), 'c':'Orange'},
                        31:{'brick_id':2, 'loc':(6*dx, 6*dy, 0), 'c':'Orange'},

                        32:{'brick_id':2, 'loc':(3*dx, 7*dy, 0), 'c':'Brown'},

                        },
                'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
                'loads_rp':{}, 'mesh_size':0.8, 'mu':0.2}

radius = 8

explicit_par_snowman = {
              'mass_scale_t': 14e-8, 't_step': 0.00005, 'is_acc': 0, # 0.001
              'load_str': '-XXmps-06mm',
              'loads_rigid': {1:{'shape':'sphere', 'radius':radius, 
                                 'loc':(3.5*dx, dy*3.5, -dz/2-radius),
                                 'u':(0,0,6) # 40
                                 }}
              }

explicit_par_snowman = {
              'mass_scale_t': 14e-8, 't_step': 0.001, 'is_acc': 0, # 0.001
              'load_str': '-40mps-40mm-higher',
              'loads_rigid': {1:{'shape':'sphere', 'radius':radius, 
                                 'loc':(3.5*dx, dy*4, -dz/2-radius),
                                 'u':(0,0,40) # 40
                                 }}
              }

make_model(assembly_pumpkin, explicit_par_snowman, is_new=0,n_frames_expl=100)

# only video evaluation, e.g. for changed view
#make_dir('pumpkin-expl-mesh080mm', if_clear=0, if_change=1)
#create_video('pumpkin-40mps-15mm-001_0ms')