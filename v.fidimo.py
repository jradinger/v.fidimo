#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################
############# Import Libraries ############
###########################################

import os
import sys
from grass.script import core as grass
import math
import sqlite3
import numpy

from igraph import *

from itertools import repeat

from timeit import default_timer as timer

from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import vector as db
from grass.pygrass.modules.shortcuts import general as g

###########################################
############# Set of Functions ############
###########################################

def import_vector(	input_map, #input vector name
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
		


def fidimo_network(	input,
					strahler_col,
					shreve_col,
					network_col,
					fidimo_db_path,  
					barriers=None,
					passability_col=None,
					threshold=25):

	'''GRASS Network based on the preprocessed vector input of streams
	and barriers (if present) is calculated (v.net). 
	
	Barriers (if present) are snapped to stream network (within threshold) and reaches 
	are split at barriers. Nodes are connected to streams network and 
	lengths of reaches are updated. Nodes and barriers are stored in attribute
	table and can be distinuished by the v_type column (1: barrier, 2: node)
	'''
	grass.message(_("Preparing river network for FIDIMO"))
	
	# Input stream network		
	streams_columns = {"strahler":strahler_col,"shreve":shreve_col,"network":network_col} #dictionary {'target_name':'input_name'}
	import_vector(input_map=input, columns=streams_columns, output_map="streams_tmp")
		
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
	grass.message(_("Connect barriers (if present) and nodes to network"))
	if barriers:
		# Input barriers
		barriers_columns = {"passability":passability_col} #dictionary {'target_name':'input_name'}
		import_vector(input_map=barriers, columns=barriers_columns, output_map="barriers_tmp")

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
		column="orig_cat INTEGER, strahler INTEGER, shreve INTEGER, network INTEGER, orig_length DOUBLE, edge_length DOUBLE")
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
	
	# Copying edges and vertices to fidimo_db
	set_fidimo_db(fidimo_db_path)

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
		
	# Create fidimo_main table
	fidimo_db.execute('''CREATE TABLE fidimo_main (source INTEGER, target INTEGER, distance DOUBLE, direction INTEGER, 
					from_orig INTEGER, to_orig INTEGER, source_edge_length DOUBLE, target_edge_length DOUBLE,
					upr_limit DOUBLE, lwr_limit DOUBLE, source_strahler INTEGER, source_shreve INTEGER, 
					target_shreve INTEGER, network INTEGER)''')
					
	# Create barriers_along table
	fidimo_db.execute('''CREATE TABLE barriers_along (source INTEGER, target INTEGER, barrier INTEGER)''')			
					
		
	fidimo_database.commit()
	fidimo_database.close()




def add_midpoints_fidimo_db(fidimo_db_path):

	'''Add midpoints (vertices) to each river reach to calculate
	afterwards distances between (midpoints of) river reaches.
	Each midpoint has its own id that is linked to the cat of the 
	river reach'''
	
	grass.message(_("Add midpoints (vertices) for each river reach"))

	fidimo_database = sqlite3.connect(fidimo_db_path)
	fidimo_db = fidimo_database.cursor()

	# Get midpoint vertices(reach id of grass edges) and insert in fidimo_db 
	fidimo_db.execute('''SELECT cat FROM edges''')
	e_cat = [(x[0]) for x in fidimo_db.fetchall()]
	
	fidimo_db.execute('''SELECT cat FROM vertices''')
	v_cat = [(x[0]) for x in fidimo_db.fetchall()]
	
	fidimo_db.executemany("INSERT INTO vertices (cat,v_type,orig_cat) VALUES (?,?,?)",
		zip(range(max(v_cat)+1,max(v_cat)+1+len(e_cat)),[3]*len(e_cat),e_cat))

	# Create index on vertices.cat
	fidimo_db.execute('''CREATE INDEX v_index ON vertices (cat)''')

	#### Alter edges in DB ####
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN part INTEGER''')
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN from_orig INTEGER''')
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN to_orig INTEGER''')
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN from_v INTEGER''')
	fidimo_db.execute('''ALTER TABLE edges ADD COLUMN to_v INTEGER''')

	fidimo_db.execute('''CREATE TEMP TABLE edges_tmp 
							(cat INTEGER, from_orig INTEGER, to_orig INTEGER)''')
	e = [(int(x.split()[0]),int(x.split()[1]),int(x.split()[2])) for x in grass.read_command("v.net",
		quiet=True,
		input="fidimo_net3",
		operation="report",
		arc_layer=3).splitlines()]
	fidimo_db.executemany("INSERT INTO edges_tmp (cat, from_orig, to_orig) VALUES (?,?,?)", e)

	fidimo_db.execute('''UPDATE edges SET 
								from_orig = (SELECT from_orig FROM edges_tmp WHERE cat=edges.cat),
								to_orig = (SELECT to_orig FROM edges_tmp WHERE cat=edges.cat)
							WHERE EXISTS (SELECT cat FROM edges_tmp WHERE cat=edges.cat)''')


	fidimo_db.execute('''UPDATE edges SET part=1''')

	# Doubling edge entries to get upstream and downstream part of midpoint
	fidimo_db.execute('''INSERT INTO edges (cat,orig_cat,strahler,shreve,network,orig_length,edge_length,
								part,from_orig,to_orig) 
		SELECT cat,orig_cat,strahler,shreve,network,orig_length,edge_length,
								2 AS part,from_orig,to_orig FROM edges''')
											
	# Get cat of vertices for from and to columns while splitting each reach at midpoint into two parts
	fidimo_db.execute('''UPDATE edges SET from_v = (SELECT cat FROM vertices 
		WHERE vertices.cat=edges.from_orig AND vertices.v_type!=3) WHERE edges."part"=1;''')
	fidimo_db.execute('''UPDATE edges SET from_v = (SELECT cat FROM vertices 
	WHERE vertices.orig_cat=edges.cat AND vertices.v_type=3) WHERE edges."part"=2;''')

	fidimo_db.execute('''UPDATE edges SET to_v = (SELECT cat FROM vertices 
	WHERE vertices.orig_cat=edges.cat AND vertices.v_type=3) WHERE edges."part"=1;''')
	fidimo_db.execute('''UPDATE edges SET to_v = (SELECT cat FROM vertices 
	WHERE vertices.cat=edges.to_orig AND vertices.v_type!=3) WHERE edges."part"=2;''')
	
	#Add index to double-cat column
	fidimo_db.execute('''CREATE INDEX edges_index ON edges (part,cat)''')

	fidimo_database.commit()
	fidimo_database.close()


def fidimo_distance(fidimo_db_path):
	'''This function fetches all vertices/edges from the
	corresponding tables of the fidimo database and builds
	the directed graph/network based on the igraph library.
	
	From the directed graph distance between the midpoints
	of the river reaches are calculated and stored in the 
	fidimo_main table along with other attributes of each
	from/to river reach.
	
	The output is the populated fidimo_main table with all
	necessary attributes and values needed to calculate the
	dispersal kernel probabilities '''
	
	# connect to database
	fidimo_database = sqlite3.connect(fidimo_db_path)
	fidimo_db = fidimo_database.cursor()
	
	# Add midpoints to calculate distance between single river reaches
	add_midpoints_fidimo_db(fidimo_db_path)
	
	# Update main distance matrix
	grass.message(_("Update main distance matrix (fidimo_main) between all connected river reaches..."))

	fidimo_db.execute('SELECT DISTINCT network FROM edges')
	network_id = [x[0] for x in fidimo_db.fetchall()]

	for i in network_id:
		grass.message(_("	...processing network ID: %s" %str(i)))

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
		direction_mat = numpy.select([	numpy.isinf(distance_mat_upstream)&~numpy.isinf(distance_mat_downstream),
										numpy.isinf(distance_mat_downstream)&~numpy.isinf(distance_mat_upstream),
										numpy.isinf(distance_mat_upstream)&numpy.isinf(distance_mat_downstream),
										~numpy.isinf(distance_mat_upstream)&~numpy.isinf(distance_mat_downstream)],
									[1,2,3,4],default=-99999)
		
		grass.message(_("Updating distance matrix for specific network"))
		paths_to_db = zip(		[x for item in midpoints for x in repeat(item,len(midpoints))],  #  list of all source/from midpoints of reaches
								midpoints*len(midpoints), #  list of all target/to midpoints of reaches
								[item for sublist in distance_mat_all for item in sublist], # list of distances
								[item for sublist in direction_mat for item in sublist], # list of directions
								[i]*(len(midpoints)**2)) # network

	
		fidimo_db.executemany("INSERT INTO fidimo_main (source,target,distance,direction,network) VALUES (?,?,?,?,?)",paths_to_db)
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
	# Get cats of original river reaches (from-to)
	grass.message(_("Updating original river reach categories to fidimo_main"))
	fidimo_db.execute('''UPDATE fidimo_main SET 
								from_orig = (SELECT orig_cat FROM vertices WHERE cat=fidimo_main.source),
								to_orig = (SELECT orig_cat FROM vertices WHERE cat=fidimo_main.target)''')
	
	# Get edge lengths, stream order and network id for source (and target) reach
	grass.message(_("Updating addtional attributes (e.g. stream order) in fidimo_main"))
	fidimo_db.execute('''UPDATE fidimo_main SET
								source_edge_length = (SELECT edge_length FROM edges WHERE cat=fidimo_main.from_orig AND part=1),
								target_edge_length = (SELECT edge_length FROM edges WHERE cat=fidimo_main.to_orig AND part=1),
								source_strahler = (SELECT strahler FROM edges WHERE cat=fidimo_main.from_orig AND part=1),
								source_shreve = (SELECT shreve FROM edges WHERE cat=fidimo_main.from_orig AND part=1),
								target_shreve = (SELECT shreve FROM edges WHERE cat=fidimo_main.to_orig AND part=1),
								network = (SELECT network FROM edges WHERE cat=fidimo_main.from_orig AND part=1)''')

	# Get upr and lwr limit for integration based on the reach lengths and distance
	grass.message(_("Updating addtional attributes (e.g. distance) in fidimo_main"))
	fidimo_db.execute('''UPDATE fidimo_main SET
								lwr_limit = distance-target_edge_length/2,
								upr_limit = distance+target_edge_length/2''')
								
								
	fidimo_database.commit()
	fidimo_database.close()
	
	end = timer()

	grass.message(_("Time elapsed: %s" %str(end-start)))
	

###########################################
########## Run FIDIMO from input ##########
###########################################

barriers="barriers"
passability_col="passability"
input="streams_fidimo"
input="stream_network_order_test2@vstreamorder"
threshold=25
strahler_col="strahler"
shreve_col="shreve"
network_col="network"
passability_col="passability"

### Get and set DB settings
# get DB settings
old_db = grass.read_command("db.connect", flags="p")

# Set GRASS DB
env = grass.gisenv()
gisdbase = env['GISDBASE']
location = env['LOCATION_NAME']
mapset = env['MAPSET']

### Check if GRASS DB exists
# otherwise create connection


fidimo_network(input=input,
				strahler_col=strahler_col,
				shreve_col=shreve_col,
				network_col=network_col,
				fidimo_db_path = '/home/radinger/Documents/FIDIMO_DB/fidimo_sqlite.db',
				barriers=barriers,
				passability_col=passability_col,
				threshold=threshold)
					
fidimo_distance(fidimo_db_path = '/home/radinger/Documents/FIDIMO_DB/fidimo_sqlite.db')



################################################
#layout=g.layout("tree")
#plot(g, layout=layout)
	
#vert_cat = [x.split("|")[0] for x in vert_grass]
#vert_names = [x.split("|")[3] for x in vert_grass]
	

		

#e = [(int(x.split()[0]),int(x.split()[1]),int(x.split()[2])) for x in grass.read_command("v.net",
#	input="fidimo_net2",
#	operation="report").splitlines()]
#e_test = [(1,1,2),(1,5,3),(2,31,1),]

#python '/home/radinger/Documents/FIDIMO_DB/fidimo_test.py' 
