"""Running test case 3 of the Lego model
"""

from brickfem import make_model

# Case 3: tower of 2x2 bricks
# -----------------------------------------------------------------------
assembly_case_3 = {'name':'case-3-tower6_2x2',
                  'bricks':{1:{'type':'base-plate', 'nx':2, 'nz':2},
                            2:{'type':'regular', 'nx':2, 'nz':2}},
                  'parts':{1:{'brick_id':1, 'loc':(0,0,0), 'c':'Yellow'},
                         2:{'brick_id':2, 'loc':(0,0,0), 'c':'Red'},
                         3:{'brick_id':2, 'loc':(0,9.6,0), 'c':'Yellow'},
                         4:{'brick_id':2, 'loc':(0,2*9.6,0), 'c':'Red'},
                         5:{'brick_id':2, 'loc':(0,3*9.6,0), 'c':'Yellow'},
                         6:{'brick_id':2, 'loc':(0,4*9.6,0), 'c':'Red'},
                         7:{'brick_id':2, 'loc':(0,5*9.6,0), 'c':'Yellow'}},
                  'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
                  'loads_rp':{}, 'mesh_size':0.75, 'mu':0.2}

explicit_par_3 = {'t_step': 0.0005, 'is_acc': 0, 'mass_scale_t': 0, 
                  'load_str': '', 
                  'loads_rigid': {1:{'shape':'sphere',  'loc':(-8.001,9.6*4.5,4),
                                     'radius':4., 'u':(20,0,0)}}}

make_model(assembly_case_3, explicit_par_3, is_new=1)

# same model, explicit load step with 0.2 ms instead of 0.5 ms
explicit_par_3['t_step'] = 0.0002
make_model(assembly_case_3, explicit_par_3, is_new=0)
