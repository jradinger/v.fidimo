#!/usr/bin/env python
# -*- coding: utf-8 -*-
############################################################################
#
# MODULE:       FIDIMO Fish Dispersal Model for Vector River Networks for GRASS 7
#
# AUTHOR(S):    Johannes Radinger
#               
# VERSION:      V0.1 Beta
#
# DATE:         2013-04-11
#
#############################################################################
#%Module
#% description: Calculating fish dispersal in a river network from source populations with species specific dispersal parameters
#% keyword: Fish Dispersal Model
#%End
#%option
#% key: input
#% type: string
#% gisprompt: old,vector,vector
#% description: River network (Vector output from v.stream.order)
#% required: yes
#% guisection: Network parameters
#%end
#%option
#% key: strahler_col
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Column name indicating Strahler order of each stream segment
#% guisection: Network parameters
#%end
#%option
#% key: shreve_col
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Column name indicating Shreve order of each stream segment
#% guisection: Network parameters
#%end
#%option
#% key: network_col
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Column name indicating Network ID each stream segment belongs to
#% guisection: Network parameters
#%end
#%option
#% key: barriers
#% type: string
#% gisprompt:old,vector,vector
#% description: Barrier point file (vector map)
#% required: no
#% guisection: Barrier parameters
#%end
#%option
#% key: passability_col
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Column name indicating passability rate (0-1) of each barrier
#% guisection: Barrier parameters
#%end
#%option
#% key: threshold
#% type: string
#% key_desc: distance
#% description: Snapping distance threshold of barriers (in m)
#% required: no
#% guisection: Barrier populations
#%end
#%Option
#% key: l
#% type: integer
#% required: no
#% multiple: no
#% description: Fish Length [mm] (valid range=39-810)
#% guisection: Dispersal parameters
#%end
#%Option
#% key: ar
#% type: double
#% required: no
#% multiple: no
#% description: Aspect Ratio of Caudal Fin (valid range 0.51 - 2.29)
#% guisection: Dispersal parameters
#%end
#%Option
#% key: t
#% type: integer
#% required: no
#% multiple: no
#% description: Time interval for model step [days]
#% guisection: Dispersal parameters
#% options: 1-3285
#% answer: 30
#%end
#% key: p
#% type: double
#% required: no
#% multiple: no
#% description: Share of the stationary component (valid range 0 - 1)
#% answer:0.67 
#% guisection: Dispersal parameters
#%end
#%Option
#% key: seed_fishmove
#% type: integer
#% required: no
#% multiple: no
#% description: fixed seed for calculating dispersal parameters using 'fishmove'
#% guisection: Dispersal parameters
#%End
#%Option
#% key: statistical_interval
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Statistical Intervals (level=0.95)
#% guisection: Output
#% options:no,confidence_interval,prediction_interval
#% answer:no
#%end


###########################################
############# Import Libraries ############
###########################################

import os
import sys
from grass.script import core as grass
import math
import sqlite3

from igraph import *

from itertools import repeat

from timeit import default_timer as timer

from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import vector as db
from grass.pygrass.modules.shortcuts import general as g

# import required numpy/scipy modules
import numpy
from scipy import stats
from scipy import optimize

###########################################
############# Set of Functions ############
###########################################

def import_vector(  input_map, #input vector name
					output_map, #output vector name
					columns): # dictionary {'target_name':'input_name'} 
					
	''' Import of vector maps for stream network and barriers.
	Copies selected input map to the specified output map. Only
	specified columns are copied and renamed if necessary'''
	
	grass.run_command("g.copy", 
		overwrite=True,
		quiet=True,
		vector="%s,%s" %(input_map,output_map))
	old_columns = grass.read_command("db.columns",
		table=output_map).splitlines()
	
	if sum([x not in old_columns for x in columns.values()])>0:
		grass.fatal(_("At least one specified column does not exist in input map"))
	
	if sorted(columns.values()+["cat"])==sorted(old_columns):
		pass
	else:
		grass.run_command("v.db.dropcolumn",
			quiet=True,
			map=output_map,
			columns=[x for x in old_columns if x not in columns.values()+["cat"]])
			
	for i in columns:
		if columns[i]!=i:
			grass.run_command("v.db.renamecolumn",
				quiet=True,
				map=output_map,
				column=",".join([columns[i],i]))
	
	grass.run_command("v.db.addcolumn",
		quiet=True,
		map=output_map,
		columns="orig_cat INTEGER")
	grass.run_command("v.db.update",
			quiet=True,
			map=output_map,
			column="orig_cat",
			query_column="cat")
		


def fidimo_network( input,
					strahler_col,
					shreve_col,
					network_col,
					barriers=None,
					passability_col=None,
					threshold=25):
	
	'''GRASS Network based on the preprocessed vector input of streams
	and barriers (if present) is calculated (v.net). 
	Barriers (if present) are snapped to stream network (within threshold) and reaches 
	are split at barriers. Nodes are connected to streams network and 
	lengths of reaches are updated. Nodes and barriers are stored in attribute
	table and can be distinuished by the v_type column (1: barrier, 2: node)'''
	
	grass.message(_("Preparing river network for FIDIMO"))
	
	# Input stream network      
	streams_col_dict = {"strahler":strahler_col,"shreve":shreve_col,"network":network_col} #dictionary {'target_name':'input_name'}
	import_vector(input_map=input, columns=streams_col_dict, output_map="streams_tmp")
	
	# Get original length of reaches
	grass.run_command("v.db.addcolumn",
		quiet=True,
		map="streams_tmp",
		columns="orig_length DOUBLE")
	grass.run_command("v.to.db",
		quiet=True,
		map="streams_tmp",
		columns="orig_length",
		option="length",
		units="meters")
	
	# Connect barriers (if present) and nodes to network
	if barriers:
		grass.message(_("Connect barriers and nodes to network"))
		# Input barriers
		passability_col_dict = {"passability":passability_col} #dictionary {'target_name':'input_name'}
		import_vector(input_map=barriers, columns=passability_col_dict, output_map="barriers_tmp")
		
		# Get network ID for each barrier
		grass.run_command("v.db.addcolumn",
			quiet=True,
			map="barriers_tmp",
			columns="network INTEGER")
		grass.run_command("v.what.vect",
			map="barriers_tmp",
			column="network",
			query_map="streams_tmp",
			query_column="network",
			dmax=threshold)
		
		# Create network
		v.net(overwrite=True,
			quiet=True,
			flags="s",
			input="streams_tmp",
			points="barriers_tmp",
			output="fidimo_net1",
			operation="connect",
			threshold=threshold)
		v.net(overwrite=True,
			quiet=True,
			flags="c",
			input="fidimo_net1",
			output="fidimo_net2",
			operation="nodes")
	else:
		grass.message(_("Connect nodes to network"))
		v.net(overwrite=True,
			quiet=True,
			flags="c",
			input="streams_tmp",
			output="fidimo_net2",
			operation="nodes")
		
	# Get table with orig attributes for all (new split) network edges
	v.category(overwrite=True,
		quiet=True,
		input="fidimo_net2",
		layer=3,
		type="line",
		output="fidimo_net3",
		option="add")
	grass.run_command("v.db.addtable",
		quiet=True,
		map="fidimo_net3",
		table="edges",
		layer=3,
		column="orig_cat INTEGER, strahler INTEGER, shreve INTEGER, network INTEGER, orig_length DOUBLE, edge_length DOUBLE, from_orig_v INTEGER, to_orig_v INTEGER, from_orig_e INTEGER, to_orig_e INTEGER")
	for i in ["orig_cat","strahler","shreve","network","orig_length"]:
		grass.run_command("v.to.db",
			quiet=True,
			map="fidimo_net3",
			layer=3,
			option="query",
			columns=i,
			query_layer=1,
			query_colum=i)
	
	# Get lengths of (new split) network edges
	grass.run_command("v.to.db",
		quiet=True,
		map="fidimo_net3",
		layer=3,
		columns="edge_length",
		option="length",
		units="meters")
	
	#database-connection of current mapset
	mapset_db_settings = dict(x.split(": ") for x in grass.read_command("db.connect",
			flags="p").splitlines())
	if mapset_db_settings["driver"]!="sqlite":
		grass.fatal(_("Database driver of current mapset is not 'sqlite'."))
	mapset_database = sqlite3.connect(mapset_db_settings["database"])
	mapset_db = mapset_database.cursor()
	
	
	# Get for each edge the orig cat for the start (from_orig) and end point (to_orig)
	mapset_db.execute('''CREATE TEMP TABLE edges_tmp 
							(cat INTEGER, from_orig_v INTEGER, to_orig_v INTEGER)''')
	e = [(int(x.split()[0]),int(x.split()[1]),int(x.split()[2])) for x in grass.read_command("v.net",
		quiet=True,
		input="fidimo_net3",
		operation="report",
		arc_layer=3).splitlines()]
	mapset_db.executemany("INSERT INTO edges_tmp (cat, from_orig_v, to_orig_v) VALUES (?,?,?)", e)
	mapset_db.execute('''UPDATE edges SET 
								from_orig_v = (SELECT from_orig_v FROM edges_tmp WHERE cat=edges.cat),
								to_orig_v = (SELECT to_orig_v FROM edges_tmp WHERE cat=edges.cat)
							WHERE EXISTS (SELECT cat FROM edges_tmp WHERE cat=edges.cat)''')
	mapset_database.commit()
	mapset_database.close()
	
	# Add table for vertices in layer 2
	grass.run_command("v.db.addtable",
		quiet=True,
		map="fidimo_net3",
		table="vertices",
		layer=2,
		column="v_type INTEGER, orig_cat INTEGER")
	
	
	# Create and populate column to distinguish between barrier / node
	if barriers:
		grass.run_command("v.db.update",
			quiet=True,
			map="fidimo_net3",
			layer=2,
			column="v_type",
			value=1,
			where="cat IN (SELECT cat FROM barriers_tmp)")
		grass.run_command("v.db.update",
			quiet=True,
			map="fidimo_net3",
			layer=2,
			column="v_type",
			value=2,
			where="cat NOT IN (SELECT cat FROM barriers_tmp)")
	else:
		grass.run_command("v.db.update",
			quiet=True,
			map="fidimo_net3",
			layer=2,
			column="v_type",
			value=2)
		
	# Get barrier orig_cat
	grass.run_command("v.db.update",
		quiet=True,
		map="fidimo_net3",
		layer=2,
		column="orig_cat",
		query_column="cat")
		
	# Update network id for each barrier
	grass.run_command("v.db.join",
			quiet=True,
			map="fidimo_net3",
			layer=2,
			column="orig_cat",
			other_table="barriers_tmp",
			other_column="orig_cat",
			subset_columns="network")
	
	grass.message(_("Final networks prepared for FIDIMO"))
	
	#removing networks and left over maps
	#g.remove





def set_fidimo_db(fidimo_db_path):
	
	''' Create fidimo database (sqlite) and tables for 
	vertices and edges'''
	
	# If database exists it will be first removed
	try:
		os.remove(fidimo_db_path)
	except OSError:
		pass
	
	grass.message(_("Creating FIDIMO Database and copying edges and vertices"))
	
	fidimo_database = sqlite3.connect(fidimo_db_path)
	fidimo_db = fidimo_database.cursor()
	
	# Copy edges and vertices to FIDIMO DB ####
	grass.run_command("db.copy",
		overwrite=True,
		from_table="vertices",
		to_database=fidimo_db_path,
		to_table="vertices")
	grass.run_command("db.copy",
		overwrite=True,
		from_table="edges",
		to_database=fidimo_db_path,
		to_table="edges")
		
	# Create fidimo_distance table
	fidimo_db.execute('''CREATE TABLE fidimo_distance (fidimo_distance_id INTEGER PRIMARY KEY, source INTEGER, target INTEGER, distance DOUBLE, direction INTEGER, 
					from_orig_e INTEGER, to_orig_e INTEGER, from_orig_v INTEGER, to_orig_v INTEGER, source_edge_length DOUBLE, target_edge_length DOUBLE,
					upr_limit DOUBLE, lwr_limit DOUBLE, source_strahler INTEGER, source_shreve INTEGER, 
					target_shreve INTEGER, network INTEGER)''')
					
	# Create barriers_along table
	fidimo_db.execute('''CREATE TABLE barriers_along (source INTEGER, target INTEGER, barrier INTEGER)''')          
					
		
	fidimo_database.commit()
	fidimo_database.close()




def add_midpoints_fidimo_db(fidimo_db_path):
	
	'''Add midpoints (vertices) to each river reach to calculate
	distances between (midpoints of) river reaches.
	Each midpoint has its own id that is linked to the cat of the 
	river reach. Function is used in fidimo_distance'''
		
	grass.message(_("Add midpoints (vertices) for each river reach"))
	
	fidimo_database = sqlite3.connect(fidimo_db_path)
	fidimo_db = fidimo_database.cursor()
	
	# Get midpoint vertices (reach id of grass edges) and insert in fidimo_db 
	fidimo_db.execute('''SELECT cat FROM edges''')
	e_cat = [(x[0]) for x in fidimo_db.fetchall()]
	
	fidimo_db.execute('''SELECT cat FROM vertices''')
	v_cat = [(x[0]) for x in fidimo_db.fetchall()]
	
	# add midpoints with category values that start with the minimum value + 1 of the existing verices-cats
	# E.g. if the existing vertices have cats 1....5, then the midpoints cats (v_type=3) will start with 6...
	fidimo_db.executemany("INSERT INTO vertices (cat,v_type,orig_cat) VALUES (?,?,?)",
		zip(range(max(v_cat)+1,max(v_cat)+1+len(e_cat)),[3]*len(e_cat),e_cat))
	
	# Create index on vertices.cat
	fidimo_db.execute('''CREATE INDEX v_index_1 ON vertices (cat)''')
	fidimo_db.execute('''CREATE INDEX v_index_2 ON vertices (orig_cat)''')
	
	
	#### Alter edges in DB ####
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN part INTEGER''')
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN from_v INTEGER''')
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN to_v INTEGER''')
	
	fidimo_db.execute('''UPDATE edges SET part=1''')
	
	# Doubling edge entries to get upstream and downstream part of midpoint
	fidimo_db.execute('''INSERT INTO edges (cat,orig_cat,strahler,shreve,network,orig_length,edge_length,
								part,from_orig_v,to_orig_v) 
		SELECT cat,orig_cat,strahler,shreve,network,orig_length,edge_length,
								2 AS part,from_orig_v,to_orig_v FROM edges''')
											
	# Get cat of vertices for from and to columns while splitting each reach at midpoint into two parts
	fidimo_db.execute('''UPDATE edges SET from_v = (SELECT cat FROM vertices 
		WHERE vertices.cat=edges.from_orig_v AND vertices.v_type!=3) WHERE edges."part"=1;''')
	fidimo_db.execute('''UPDATE edges SET from_v = (SELECT cat FROM vertices 
	WHERE vertices.orig_cat=edges.cat AND vertices.v_type=3) WHERE edges."part"=2;''')
	
	fidimo_db.execute('''UPDATE edges SET to_v = (SELECT cat FROM vertices 
	WHERE vertices.orig_cat=edges.cat AND vertices.v_type=3) WHERE edges."part"=1;''')
	fidimo_db.execute('''UPDATE edges SET to_v = (SELECT cat FROM vertices 
	WHERE vertices.cat=edges.to_orig_v AND vertices.v_type!=3) WHERE edges."part"=2;''')
	
	#Add index to double-cat column
	fidimo_db.execute('''CREATE INDEX e_index_1 ON edges (part,cat)''')
	fidimo_db.execute('''CREATE INDEX e_index_2 ON edges (orig_cat)''')
	
	
	fidimo_database.commit()
	fidimo_database.close()


def fidimo_distance(fidimo_db_path):
	'''This function fetches all vertices/edges from the
	corresponding tables of the fidimo database and builds
	the directed graph/network based on the igraph library.
	
	From the directed graph distance between the midpoints
	of the river reaches are calculated and stored in the 
	fidimo_distance table along with other attributes of each
	from/to river reach.
	
	The output is the populated fidimo_distance table with all
	necessary attributes and values needed to calculate the
	dispersal kernel probabilities '''
	
	# connect to database
	fidimo_database = sqlite3.connect(fidimo_db_path)
	fidimo_db = fidimo_database.cursor()
	
	# Add midpoints to calculate distance between single river reaches
	add_midpoints_fidimo_db(fidimo_db_path)
	
	# Update main distance matrix
	grass.message(_("Update main distance matrix (fidimo_distance) between all connected river reaches..."))
	
	fidimo_db.execute('SELECT DISTINCT network FROM edges')
	network_id = [x[0] for x in fidimo_db.fetchall()]
	
	for i in network_id:
		grass.message(_("   ...processing network ID: %s" %str(i)))
		
		# Fetch all edges
		fidimo_db.execute('SELECT from_v, to_v, edge_length FROM edges WHERE network = ?', str(i))
		edges_attributes = fidimo_db.fetchall()
		edges = [(x[0],x[1]) for x in edges_attributes]
		
		# Get vertices and unique ids for graph (starting with 0)
		vertices = list(set([x[y] for x in edges_attributes for y in [0,1]]))
		g_vertices = range(len(vertices))
		vertices_dict = dict(zip(vertices,g_vertices))
		inv_vertices_dict = {v:k for k, v in vertices_dict.iteritems()}
		
		# Get new unique igraph ids for vertices in edges list
		g_edges = [(vertices_dict[x[0]],vertices_dict[x[1]]) for x in edges_attributes]
		
		# Fetch all vertices for specific network that are midpoints of river reaches 
		fidimo_db.execute('SELECT cat FROM vertices WHERE v_type=3 AND orig_cat IN (SELECT cat from edges WHERE network = ? )',str(i))
		midpoints = [x[0] for x in fidimo_db.fetchall()]
		g_midpoints = [vertices_dict[x] for x in midpoints]
		
		# Build Graph
		g = Graph(g_edges,directed=True)
		
		# Add length to edges (half length as each edge is split by the midpoint)
		g.es["half_length"] = [(x[2])/2.0 for x in edges_attributes]
		
		### Calculate shortest paths
		grass.message(_("Calculating shortest paths"))
		distance_mat_upstream = g.shortest_paths(source=g_midpoints,
						target=g_midpoints, weights="half_length",mode="IN")
		distance_mat_downstream = g.shortest_paths(source=g_midpoints,
						target=g_midpoints, weights="half_length",mode="OUT")
		distance_mat_all = g.shortest_paths(source=g_midpoints,
						target=g_midpoints, weights="half_length",mode="ALL")
				
		# Distinguish between paths up- and downstream [1: downstream, 2:upstram, 3:neither (down-up combination), 4: both directions (e.g. where source=target and dist=0)
		grass.message(_("Calculating directions between reaches"))
		direction_mat = numpy.select([  numpy.isinf(distance_mat_upstream)&~numpy.isinf(distance_mat_downstream),
										numpy.isinf(distance_mat_downstream)&~numpy.isinf(distance_mat_upstream),
										numpy.isinf(distance_mat_upstream)&numpy.isinf(distance_mat_downstream),
										~numpy.isinf(distance_mat_upstream)&~numpy.isinf(distance_mat_downstream)],
									[1,2,3,4],default=-99999)
		
		grass.message(_("Updating distance matrix for specific network"))
		paths_to_db = zip(      [x for item in midpoints for x in repeat(item,len(midpoints))],  #  list of all source/from midpoints of reaches
								midpoints*len(midpoints), #  list of all target/to midpoints of reaches
								[item for sublist in distance_mat_all for item in sublist], # list of distances
								[item for sublist in direction_mat for item in sublist], # list of directions
								[i]*(len(midpoints)**2)) # network
		
		fidimo_db.executemany("INSERT INTO fidimo_distance (source,target,distance,direction,network) VALUES (?,?,?,?,?)",paths_to_db)
		fidimo_database.commit()
		
		# Barriers along each path
		fidimo_db.execute('''SELECT cat FROM vertices WHERE v_type=1 and network=?''', str(i))
		barriers_list = [x[0] for x in fidimo_db.fetchall()]
		
		if barriers_list:
			grass.message(_("Processing barriers along networks paths"))
			
			g_barriers = [vertices_dict[x] for x in barriers_list]
			
			for j in range(len(midpoints)):
				#j=0
				vertices_paths = g.get_shortest_paths(v=g_midpoints[j],to=g_midpoints,
								weights="half_length",output="vpath",mode="ALL")
				paths_with_barriers = [[x for x in L if x in g_barriers] for L in vertices_paths]
				barriers_along = [(midpoints[j], inv_vertices_dict[b], inv_vertices_dict[e]) for a, b in zip(paths_with_barriers, g_midpoints) for e in a]
		
				fidimo_db.executemany('''INSERT INTO barriers_along (source,target,barrier) VALUES (?,?,?)''',barriers_along)
		
	grass.message(_("All networks processed"))
	
	##### Update values for main fidimo table (e.g. lengths, stream order,...)
	start = timer()
	
	#Create index on fidimo_distance
	#grass.message(_("Creating DB Index on fidimo_distance"))
	#fidimo_db.execute('''CREATE INDEX fidimo_distance_index_source ON fidimo_distance (source)''')
	#fidimo_db.execute('''CREATE INDEX fidimo_distance_index_target ON fidimo_distance (target)''')
	
	# Get cats of original vertices (from-to)
	grass.message(_("Updating original vertex categories to fidimo_distance"))
	fidimo_db.execute('''UPDATE fidimo_distance SET 
								from_orig_v = (SELECT orig_cat FROM vertices WHERE cat=fidimo_distance.source),
								to_orig_v = (SELECT orig_cat FROM vertices WHERE cat=fidimo_distance.target)''')
	fidimo_database.commit()
	
	#Create index on fidimo_distance
	#grass.message(_("Creating DB Index on fidimo_distance"))
	#fidimo_db.execute('''CREATE INDEX fidimo_distance_index_from_orig_v ON fidimo_distance (from_orig_v)''')
	#fidimo_db.execute('''CREATE INDEX fidimo_distance_index_to_orig_v ON fidimo_distance (to_orig_v)''')
	
	# Get cats of original edges (from-to)
	grass.message(_("Updating original river reach (edges) categories fidimo_distance"))
	fidimo_db.execute('''UPDATE fidimo_distance SET 
								from_orig_e = (SELECT orig_cat FROM edges WHERE cat=fidimo_distance.from_orig_v),
								to_orig_e = (SELECT orig_cat FROM edges WHERE cat=fidimo_distance.to_orig_v)''')
	
	# Get edge lengths, stream order and network id for source (and target) reach
	grass.message(_("Updating addtional attributes (e.g. stream order) in fidimo_distance"))
	fidimo_db.execute('''UPDATE fidimo_distance SET
								source_edge_length = (SELECT edge_length FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1),
								target_edge_length = (SELECT edge_length FROM edges WHERE cat=fidimo_distance.to_orig_v AND part=1),
								source_strahler = (SELECT strahler FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1),
								source_shreve = (SELECT shreve FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1),
								target_shreve = (SELECT shreve FROM edges WHERE cat=fidimo_distance.to_orig_v AND part=1),
								network = (SELECT network FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1)''')
	
	# Get upr and lwr limit for later integration based on the reach lengths and distance
	grass.message(_("Updating addtional attributes (upper and lower distance limit) in fidimo_distance"))
	fidimo_db.execute('''UPDATE fidimo_distance SET
								lwr_limit = distance-(target_edge_length/2),
								upr_limit = distance+(target_edge_length/2)''')
	
	fidimo_database.commit()
	fidimo_database.close()
	
	end = timer()
	
	grass.message(_("Time elapsed: %s" %str(end-start)))



def sigma_calc( l,
				ar,
				t,
				statistical_interval,
				fishmove_seed=None):
	'''This function calculates dispersal distances sigma_stat and sigma_mob
	for each stream order and for given input: l, ar, t'''
	
	
	# Regression model to calculate dispersal kernel parameter sigma_stat and sigma_mob
	if statistical_interval == "no":
		#log(sigma_stat) = -10.57 + 1.64*log(L) + 0.94*AR + 1.14*sqrt(SO) + 0.43*log(T)
		#log(sigma_mob) = -7.48 + 1.45*log(L) + 0.58*AR + 1.51*sqrt(SO) + 0.55*log(T)
		# Parameter coefficients calculated from fishmove(L=200,AR=1.25,rep=5000,seed=999)
		sigma_dict = {
			"stat" : {
				"fit" : dict(zip(range(1,10), [math.exp(-10.56720 + 1.643237*math.log(l) + 0.9647056*ar + 1.142661*math.sqrt(i) + 0.4273618*math.log(t)) for i in (1,2,3,4,5,6,7,8,9)]))},
			"mob" : {
				"fit" : dict(zip(range(1,10), [math.exp(-7.480260 + 1.445552*math.log(l) + 0.5820339*ar + 1.507528*math.sqrt(i) + 0.5527266*math.log(t)) for i in (1,2,3,4,5,6,7,8,9)]))}}
				
		return sigma_dict
		
	else:
		import rpy2.robjects as robjects # import required rpy2 module
		from rpy2.robjects.packages import importr
		fm = importr('fishmove')
		
		##### Calculating 'fishmove' depending on species or L & AR
		# Statistical interval
		if "prediction" in statistical_interval:
			interval = "prediction"
		elif "confidence" in statistical_interval:
			interval = "confidence"
		
		#Set fixed seed if specified
		if fishmove_seed:
			seed = ",seed="+str(fishmove_seed)
		else:
			seed = ""
		
		so_rvector = robjects.IntVector((1,2,3,4,5,6,7,8,9))
		 
		# Calculate 'fishmove' and store sigma values in pandas df
		fishmove = eval("fm.fishmove(L=l,AR=ar,SO=so_rvector,T=t,interval=interval,rep=200%s)"%(seed))
		fishmove = fishmove[1]
		sigma_dict = {
			"stat" : {
				"fit" : dict(zip(range(1,10), list(fishmove.rx("fit",'sigma_stat',1,1,so_rvector,1)))),
				"lwr" : dict(zip(range(1,10), list(fishmove.rx("lwr",'sigma_stat',1,1,so_rvector,1)))),
				"upr" : dict(zip(range(1,10), list(fishmove.rx("upr",'sigma_stat',1,1,so_rvector,1))))},
			"mob" : {"fit" : dict(zip(range(1,10), list(fishmove.rx("fit",'sigma_mob',1,1,so_rvector,1)))),
				"lwr" : dict(zip(range(1,10), list(fishmove.rx("lwr",'sigma_mob',1,1,so_rvector,1)))),
				"upr" : dict(zip(range(1,10), list(fishmove.rx("upr",'sigma_mob',1,1,so_rvector,1))))}}
		
		return sigma_dict
		

def fidimo_kernel_cdf_truncation(x,sigma_stat,sigma_mob,truncation,p):
	'''This function calculates the 
	maximum distance (cutting distance) 
	based on truncation criterion.'''
	return p * stats.norm.cdf(x, loc=0, scale=sigma_stat) + (1-p) * stats.norm.cdf(x, loc=0, scale=sigma_mob) - truncation


def fidimo_kernel_cdf(x,sigma_stat,sigma_mob,p):
	'''This function calculates the dispersal probability 
	based on leptokurtic dispersal kernels for river fish (Radinger and Wolter, 2014).
	Dispersal probability is calculated based on the cumulative density function
	between lower and upper distance.'''
	return (p * stats.norm.cdf(x, loc=0, scale=sigma_stat) + (1-p) * stats.norm.cdf(x, loc=0, scale=sigma_mob))


def fidimo_source_pop(	input,
						source_col,
						output,
						fidimo_db_path):
	'''This function appends source populations to distance matrix and creates output map'''
	
	# connect to database
	fidimo_database = sqlite3.connect(fidimo_db_path)
	fidimo_db = fidimo_database.cursor()
	
	# Check if source_col exists in input
	###
	
	# Import source populations
	source_col_dict = {"source_pop":source_col} #dictionary {'target_name':'input_name'}
	import_vector(input_map=input, columns=source_col_dict, output_map=output)
	
	# Copy source populations to FIDIMO DB ####
	grass.run_command("db.copy",
		overwrite=True,
		from_table=output,
		to_database=fidimo_db_path,
		to_table="fidimo_source_pop")
	
	# Creating index on fidimo_source_pop
	fidimo_db.execute('''CREATE INDEX fidimo_source_pop ON fidimo_source_pop (orig_cat)''')
	
	#check if cat of cats of fidimo_source_pop == from_orig of fidimo_distance
	fidimo_db.execute('SELECT DISTINCT from_orig_e FROM fidimo_distance')
	from_orig_e_fidimo_distance = [x[0] for x in fidimo_db.fetchall()]
	fidimo_db.execute('SELECT DISTINCT orig_cat FROM fidimo_source_pop')
	orig_cat_fidimo_source_pop = [x[0] for x in fidimo_db.fetchall()]
	if sorted(from_orig_e_fidimo_distance) != sorted(orig_cat_fidimo_source_pop):
		grass.fatal(_("Vector input map of source populations must must match vector input map that has been used for calculating fidimo distance matrix"))
	#else:
	#	print ("Vector input map of source populations matches vector input map that has been used for calculating fidimo distance matrix")
		
	# Join table fidimo_distance with fidimo_source_pop
	grass.message(_("Updating source population in fidimo_distance"))
	fidimo_db.execute('''ALTER TABLE fidimo_distance ADD COLUMN source_pop DOUBLE''')
	fidimo_db.execute('''UPDATE fidimo_distance SET
								source_pop = (SELECT source_pop FROM fidimo_source_pop WHERE orig_cat=fidimo_distance.from_orig_e)''')
	


def fidimo_probability(	input,
						output,
						fidimo_db_path,
						sigma_dict,
						p,
						statistical_interval):
	'''This function calculates dispersal kernel probabilities
	for each river reach and corresponding distances between reaches'''
	
	# connect to database
	fidimo_database = sqlite3.connect(fidimo_db_path)
	fidimo_db = fidimo_database.cursor()
	
	# Get number of runs for statistical intervals
	if statistical_interval == "no":
		nrun=["fit"]
		# Create a new sqlite table in fidimo_db to collect results
		fidimo_db.execute('''CREATE TABLE fidimo_prob (fidimo_distance_id INTEGER, fidimo_prob DOUBLE)''')
	else:
		nrun=["fit","lwr","upr"]
		# Create a new sqlite table in fidimo_db to collect results
		fidimo_db.execute('''CREATE TABLE fidimo_prob (fidimo_distance_id INTEGER, fidimo_prob DOUBLE, fidimo_prob_lwr DOUBLE, fidimo_prob_upr DOUBLE)''')
	
	# Commit changes
	fidimo_database.commit()
	
	# Get stream orders (ASC) of source populations where source_pop > 0
	fidimo_db.execute('''SELECT DISTINCT source_strahler FROM fidimo_distance WHERE source_pop > 0 AND source_pop IS NOT NULL AND source_pop != ""''')
	so_list = sorted([int(x[0]) for x in fidimo_db.fetchall()])
	return so_list
	
	if max(so_list) > 9:
		grass.fatal(_("Stream network has stream orders > 9. Please consider smaller stream network"))
	
	for i in so_list:
		#i=3
		
		# Calculate maximum distance (cutting distance) based on truncation criterion 
		if truncation != "inf":
			truncation = float(truncation)
			max_dist = float(optimize.zeros.newton(fidimo_kernel_cdf_truncation, 1.,  
				args=(sigma_dict["stat"]["fit"][i],sigma_dict["mob"]["fit"][i],truncation,p)))
		else
			fidimo_db.execute('''SELECT max(distance) FROM fidimo_distance WHERE source_strahler = ?''', str(i))
   			max_dist = [x[0] for x in fidimo_db.fetchall()][0]*2 # Get a distance that is definintly (2x) larger than max distance
   			
		direction_dict = {"upstream":1,"downstream":2,"source":3,"undef":4}
		
		for j in direction_dict: # loop over directions (up- vs. downstream)
			''' calculate fidimo probability in chunks based on source_strahler
			MAIN PART: leptokurtic probability density kernel based on fishmove '''
			
			grass.debug(_("Begin with core of fidimo, application of fishmove"))
			
			# Do not calculate any fidimo probability for direction that are combinations of down- and upstream (as in r.fidimo)
			if j=="undef":
				continue
			
			# Get all IDs for cases where distance, direction, source_pop and source_strahler match criteria
			fidimo_db.execute('''SELECT fidimo_distance_id FROM fidimo_distance 
				WHERE source_strahler = ?
					AND distance <= ? 
					AND direction = ?
					AND source_pop > 0''', (int(i),max_dist,int(direction_dict(j)))
			
			# Split all cases into chunks of max. chunk_size		
			row_ids = [x[0] for x in fidimo_db.fetchall()]
			chunk_size=500
			row_ids_chunks = [row_ids[x:x+chunk_size] for x in xrange(0, len(row_ids), chunk_size)]
			
			for k in row_ids_chunks
				fidimo_db.execute('SELECT * FROM fidimo_distance WHERE fidimo_distance_id IN (%s)' %','.join('?'*len(k)), tuple(k))
				fidimo_distance_colnames = dict(zip([x[0] for x in fidimo_db.description], range(0,len(fidimo_db.description))))
				fidimo_distance_array = scipy.array(fidimo_db.fetchall())
				
				# Create results array to collect fidimo_prob results
				fidimo_result_array = k
				
				for l in nrun:
					# CDF(upr) - CDF(lwr)
					if j=="source":
						fidimo_prob = 	(fidimo_kernel_cdf(	x = fidimo_distance_array[:,col_names["upr_limit"]],
															sigma_stat = sigma_dict["stat"][l][i],
															sigma_mob = sigma_dict["mob"][l][i],
															p = p) - 
										fidimo_kernel_cdf(	x = 0,
															sigma_stat = sigma_dict["stat"][l][i],
															sigma_mob = sigma_dict["mob"][l][i],
															p = p))*2.0
					else:
						fidimo_prob = 	fidimo_kernel_cdf(	x = fidimo_distance_array[:,col_names["upr_limit"]],
															sigma_stat = sigma_dict["stat"][l][i],
															sigma_mob = sigma_dict["mob"][l][i],
															p = p) - 
										fidimo_kernel_cdf(	x = fidimo_distance_array[:,col_names["lwr_limit"]],
															sigma_stat = sigma_dict["stat"][l][i],
															sigma_mob = sigma_dict["mob"][l][i],
															p = p)
					
					# if upstream direction than correct for network splits by a shreve stream order approach
					if j=="upstream":
						fidimo_prob = fidimo_prob * (fidimo_distance_array[:,col_names["target_shreve"]]]*1.0 / fidimo_distance_array[:,col_names["source_shreve"]]])
					
					# Stack fidimo_prob with fidimo_distance id into result array
					fidimo_result_array = numpy.concatenate(
											fidimo_result_array, # fidimo_distance_id
											fidimo_prob) # fidimo_prob for single nrun
				
				if statistical_interval == "no":
 					fidimo_db.executemany("INSERT INTO fidimo_prob VALUES (?,?)", fidimo_result_array)
 				else:
 					fidimo_db.executemany("INSERT INTO fidimo_prob VALUES (?,?,?,?)", fidimo_result_array)
 				
 				# Commit changes
				fidimo_database.commit()
				
	# Create key in final fidimo_prob
	
	# Join fidimo_prob with fidimo_distance
	
def barrier_correction()
	
		 

###########################################
########## Run FIDIMO from input ##########
###########################################

# Stream settings
barriers="fidimo_testbarriers"
passability_col="passability"
input="fidimo_testnetwork"
threshold=25
strahler_col="strahler"
shreve_col="shreve"
network_col="network"
passability_col="passability"

# Dispersal settings
l=200
ar=1.25
t=30
statistical_interval="no"

#fidimo database connection
fidimo_db_path='/home/radinger/Documents/FIDIMO_DB/fidimo_sqlite.db'

# Source population settings
source_col="sources_test"


# Create fidimo network from input data (river shape, barrier points)
fidimo_network(input=input,
				strahler_col=strahler_col,
				shreve_col=shreve_col,
				network_col=network_col,
				barriers=barriers,
				passability_col=passability_col,
				threshold=threshold)

# Set up fidimo_db and copying edges and vertices to fidimo_db
set_fidimo_db(fidimo_db_path=fidimo_db_path)
					
# Calculate distance between single river reaches
fidimo_distance(fidimo_db_path = fidimo_db_path)

# Calculate dispersal distances
my_sigma_dict = sigma_calc(	l=l,
							ar=ar,
							t=t,
							statistical_interval=statistical_interval,
							fishmove_seed=999)

# Append source populatoins
fidimo_source_pop(input=input,
				source_col=source_col,
				output="fidimo_test_output",
				fidimo_db_path=fidimo_db_path)

# Calcuate fidimo probability
fidimo_probability(input=input,
				output="fidimo_test_output",
				fidimo_db_path=fidimo_db_path,
				sigma_dict=my_sigma_dict,
				p=0.67,
				statistical_interval=statistical_interval)




