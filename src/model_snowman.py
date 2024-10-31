"""Running test case 3 of the Lego model
"""

from brickfem import make_model

# Big snowman of Sebastian
# -----------------------------------------------------------------------

def_y = 9.6
def_x = 8
def_z = def_x

assembly_snowman = {
                'name':'snowman',
                'bricks':{
                        1:{'type':'regular', 'nx':4, 'nz':2},
                        2:{'type':'regular', 'nx':6, 'nz':2},
                        3:{'type':'regular', 'nx':2, 'nz':4},
                        4:{'type':'regular', 'nx':2, 'nz':4},
                        5:{'type':'regular', 'nx':2, 'nz':2},
                        6:{'type':'regular', 'nx':1, 'nz':2},
                        7:{'type':'base-plate', 'nx':8, 'nz':2},
                        8:{'type':'plate', 'nx':4, 'nz':2},
                        9:{'type':'plate', 'nx':8, 'nz':2},
                        10:{'type':'regular', 'nx':2, 'nz':2},
                        11:{'type':'regular', 'nx':6, 'nz':2}
                        },
                'parts':{
                        # base plate
                        1:{'brick_id':7, 'loc':(-def_x*2, 0, 0), 'c':'White'},
                        # brick 1
                        2:{'brick_id':1, 'loc':(0, 0, 0), 'c':'White'},
                        # brick 2
                        3:{'brick_id':10, 'loc':(-def_x, def_y, 0), 'c':'White'},
                        4:{'brick_id':10, 'loc':(def_x, def_y, 0), 'c':'Black'},
                        5:{'brick_id':10, 'loc':(3*def_x, def_y, 0), 'c':'White'},
                        6:{'brick_id':10, 'loc':(def_x, def_y*2, 0), 'c':'White'},
                        7:{'brick_id':10, 'loc':(def_x, def_y*3, 0), 'c':'Black'},
                        # bricks that should have twice the height 
                        # ------ brick 5
                        8:{'brick_id':5, 'loc':(-def_x*2, def_y*2, 0), 'c':'White'},
                        9:{'brick_id':5, 'loc':(-def_x*2, def_y*3, 0), 'c':'White'},
                        10:{'brick_id':5, 'loc':(def_x*4, def_y*2, 0), 'c':'White'},
                        11:{'brick_id':5, 'loc':(def_x*4, def_y*3, 0), 'c':'White'},
                        # --- brick 6
                        12:{'brick_id':6, 'loc':(0, def_y*2, 0), 'c':'White'},
                        13:{'brick_id':6, 'loc':(0, def_y*3, 0), 'c':'White'},
                        14:{'brick_id':6, 'loc':(def_x*3, def_y*2, 0), 'c':'White'},
                        15:{'brick_id':6, 'loc':(def_x*3, def_y*3, 0), 'c':'White'},
                        # ------ end of 5 and 6
                        # brick 1
                        16:{'brick_id':1, 'loc':(-def_x*2, def_y*4, 0), 'c':'White'},
                        17:{'brick_id':1, 'loc':(def_x*2, def_y*4, 0), 'c':'White'},
                        # brick 2
                        18:{'brick_id':10, 'loc':(-def_x, def_y*5, 0), 'c':'White'},
                        19:{'brick_id':10, 'loc':(def_x, def_y*5, 0), 'c':'Black'},
                        20:{'brick_id':10, 'loc':(def_x*3, def_y*5, 0), 'c':'White'},
                        # brick 9
                        21:{'brick_id':9, 'loc':(0, def_y*6, 0), 'c':'Red'},
                        # brick 8
                        22:{'brick_id':8, 'loc':(0, def_y*(6+(1./3.)), 0), 'c':'White'},
                        # brick 2
                        23:{'brick_id':10, 'loc':(-def_x, def_y*(6+(2./3.)), 0), 'c':'White'},
                        24:{'brick_id':10, 'loc':(def_x*3, def_y*(6+(2./3.)), 0), 'c':'White'},
                        # brick 4
                        25:{'brick_id':4, 'loc':(def_x, def_y*(6+(2./3.)), 0), 'c':'Orange'},
                        # brick 8
                        26:{'brick_id':8, 'loc':(-def_x*2, def_y*(8+(2./3.)), 0), 'c':'Brown'},
                        27:{'brick_id':8, 'loc':(def_x*2, def_y*(8+(2./3.)), 0), 'c':'Brown'},
                        28:{'brick_id':8, 'loc':(0, def_y*(8+(3./3.)), 0), 'c':'Brown'},
                        # brick 1
                        29:{'brick_id':1, 'loc':(0, def_y*(8+(4./3.)), 0), 'c':'Dark Nougat'},
                        # brick 8
                        30:{'brick_id':8, 'loc':(0, def_y*(9+(4./3.)), 0), 'c':'Brown'},
                        31:{'brick_id':11, 'loc':(-def_x,def_y*(7+(2./3.)), 0), 'c':'White'}
                        },
                'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
                'loads_rp':{}, 'mesh_size':0.8, 'mu':0.2}

radius = 8

explicit_par_snowman = {
              'mass_scale_t': 14e-8, 't_step': 0.0005, 'is_acc': 0,
              'load_str': '-20mps-300f',
              'loads_rigid': {1:{'shape':'sphere', 'radius':radius, 
                                 'loc':(-def_x-4-radius,def_y*(7+(2./3.)), 4),
                                 'v0':(40000,0,0), 'm':1200e-6
                                 }}
              }

make_model(assembly_snowman, explicit_par_snowman, is_new=0, n_frames_expl=100)
