"""Running test case 3 of the Lego model
"""

from brickfem import make_model

# Case 3: tower of 2x2 bricks
# -----------------------------------------------------------------------
assembly_case_4 = {'name':'crack-bendH6',
                  'bricks':{1:{'type':'base-plate', 'nx':1, 'nz':1},
                            2:{'type':'plate', 'nx':2, 'nz':1},
                            3:{'type':'plate', 'nx':1, 'nz':1},
                            4:{'type':'tile', 'nx':2, 'nz':1}},
                  'parts':{
                         5:{'brick_id':1, 'loc':(-8*5,3.2,0)},
                         6:{'brick_id':1, 'loc':(8*4,3.2,0)},

                         9:{'brick_id':3, 'loc':(-8*5,3.2,0)},
                         10:{'brick_id':2, 'loc':(-8*4,3.2,0)},
                         11:{'brick_id':2, 'loc':(-8*2,3.2,0)},
                         12:{'brick_id':2, 'loc':(0,3.2,0)},
                         13:{'brick_id':2, 'loc':(8*2,3.2,0)},
                         14:{'brick_id':3, 'loc':(8*4,3.2,0)},
                         15:{'brick_id':2, 'loc':(-8*5,6.4,0)},
                         16:{'brick_id':2, 'loc':(-8*3,6.4,0)},
                         17:{'brick_id':2, 'loc':(-8*1,6.4,0)},
                         18:{'brick_id':2, 'loc':(8*1,6.4,0)},
                         19:{'brick_id':2, 'loc':(8*3,6.4,0)},

                         20:{'brick_id':3, 'loc':(-8*5,9.6,0)},
                         21:{'brick_id':2, 'loc':(-8*4,9.6,0)},
                         22:{'brick_id':2, 'loc':(-8*2,9.6,0)},
                         23:{'brick_id':4, 'loc':(0,9.6,0)},
                         24:{'brick_id':2, 'loc':(8*2,9.6,0)},
                         25:{'brick_id':3, 'loc':(8*4,9.6,0)},

                         26:{'brick_id':2, 'loc':(-8*5,6.4*2,0)},
                         27:{'brick_id':2, 'loc':(-8*3,6.4*2,0)},
                         28:{'brick_id':2, 'loc':(-8*1,6.4*2,0)},
                         29:{'brick_id':2, 'loc':(8*1,6.4*2,0)},
                         30:{'brick_id':2, 'loc':(8*3,6.4*2,0)}
                         
                         },
                  'bc':{},
                  'loads_rp':{1:{'part_id':5, 'set_name':'BOTTOM', 'uy':0,
                                 'uz':0, 'rotx':0, 'roty':0},
                              2:{'part_id':6, 'set_name':'BOTTOM','uy':0,
                                 'uz':0, 'rotx':0, 'roty':0},},
                  'mesh_size':0.75, 'mu':0.2}

"""
1:{'part_id':5, 'set_name':'BOTTOM'},
                        2:{'part_id':6, 'set_name':'BOTTOM'}
"""

# Schlagpendel v: 20 mm in 5 ms: 4 m/s :-)

explicit_par_4 = {'t_step': 0.005, 'is_acc': 0, 'mass_scale_t': 0, 
                  'load_str': '', 
                  'loads_rigid': {1:{'shape':'cyl',  'loc':(-4,17+6.7,0), 
                                     'dir':(0,0,1), 'len':10.,
                                     'radius':6., 'u':(0,-20,0)}}}

make_model(assembly_case_4, explicit_par_4, is_new=1)

# same model, explicit load step with 0.2 ms instead of 0.5 ms
#explicit_par_4['t_step'] = 0.0002
#make_model(assembly_case_4, explicit_par_4, is_new=0)
