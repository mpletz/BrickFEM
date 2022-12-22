"""Running test case 1 of the Lego model
"""

from brickfem import make_model

# Case 1: pulling two 1x1 bricks apart
# -----------------------------------------------------------------------
assembly_case_1 = {'name':'case-1-pull-1x1',
                   'bricks':{1:{'type':'base-plate', 'nx':1, 'nz':1},
                            2:{'type':'regular', 'nx':1, 'nz':1}},
                   'parts':{1:{'brick_id':1, 'loc':(0,0,0), 'c':'Yellow'},
                            2:{'brick_id':2, 'loc':(0,0,0), 'c':'Red'}},
                   'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
                   'loads_rp':{1:{'part_id':2, 'set_name':'STUD-11', 'uy':1.2}},
                   'mesh_size':0.5, 'mu':0.2}

make_model(assembly_case_1)

# load explicitly in a separate computation
explicit_par_1 = {'t_step': 0.0005, 'is_acc': 0, 'mass_scale_t': 0,
                  'load_str': '', 'loads_rigid': {}}

assembly_case_1['loads_rp'][1]['uy'] = 2

make_model(assembly_case_1, explicit_par_1)