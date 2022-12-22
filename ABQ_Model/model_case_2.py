"""Running test case 2 of the Lego model
"""

from brickfem import make_model

# Case 2: 2x4 brick on base-plate
# -----------------------------------------------------------------------
assembly_case_2 = {'name':'case-2-bend2x4',
                   'bricks':{1:{'type':'base-plate', 'nx':2, 'nz':2},
                             2:{'type':'regular', 'nx':4, 'nz':2}},
                   'parts':{1:{'brick_id':1, 'loc':(0,0,0), 'c':'Yellow'},
                            2:{'brick_id':2, 'loc':(0,0,0), 'c':'Red'}},
                   'bc':{1:{'part_id':1, 'set_name':'BOTTOM'}},
                   'loads_rp':{1:{'part_id':2, 'set_name':'STUD-41', 'uy':-3}},
                   'mesh_size':0.75, 'mu':0.2}


explicit_par_2 = {'mass_scale_t': 0, 't_step': 0.0005, 'is_acc': 0,
                  'load_str': '', 'loads_rigid': {}}

make_model(assembly_case_2, explicit_par_2)