#!/usr/bin/env python
# -*- coding: utf-8 -*-
############################################################################
#
# MODULE:       FIDIMO Fish Dispersal Model for Vector River Networks for GRASS 7
#
# AUTHOR(S):    Johannes Radinger
#
# VERSION:      V0.0 Beta
#
# DATE:         2016-08-30
#
#############################################################################

#%module
#% description: calculates fish dispersal in a river network from source populations with species specific dispersal parameters
#% keyword: vector
#% keyword: database
#% keyword: network
#% keyword: fish dispersal model
#%end

#%option G_OPT_F_INPUT
#% key: fidimo_db_path
#% description: FIDIMO Database directory
#% required: yes
#%end

#%option G_OPT_V_MAP
#% key: input
#% description: Vector river network (e.g. output from v.stream.order)
#% required: no
#% guisection: Network parameters
#% guidependency: strahler_col,shreve_col,network_col,source_col
#%end
#%option G_OPT_DB_COLUMN
#% key: strahler_col
#% description: Column name indicating Strahler order of each stream segment
#% required: no
#% guisection: Network parameters
#%end
#%option G_OPT_DB_COLUMN
#% key: shreve_col
#% description: Column name indicating Shreve order of each stream segment
#% required: no
#% guisection: Network parameters
#%end
#%option G_OPT_DB_COLUMN
#% key: network_col
#% description: Column name indicating Network ID each stream segment belongs to
#% required: no
#% guisection: Network parameters
#%end
#%option
#% key: barriers
#% type: string
#% gisprompt:old,vector,vector
#% description: Barrier point file (vector map)
#% required: no
#% guisection: Network parameters
#% guidependency: passability_col
#%end
#%option G_OPT_DB_COLUMN
#% key: passability_col
#% description: Column name indicating passability rate (0-1) of each barrier
#% required: no
#% guisection: Network parameters
#%end
#%option
#% key: threshold
#% type: string
#% key_desc: distance
#% description: Snapping distance threshold of barriers (in m)
#% required: no
#% answer: 25
#% guisection: Network parameters
#%end

#%option G_OPT_F_INPUT
#% key: source_pop_csv
#% description: csv-File (comma-separeted) indicating source populations (columns: reach_ID,source_col,p_col)
#% required: no
#% guisection: Source populations
#%end
#%option

#% key: l
#% type: integer
#% required: no
#% multiple: no
#% description: Fish Length [mm] (valid range=39-810)
#% guisection: Dispersal parameters
#%end
#%option
#% key: ar
#% type: double
#% required: no
#% multiple: no
#% description: Aspect Ratio of Caudal Fin (valid range 0.51 - 2.29)
#% guisection: Dispersal parameters
#%end
#%option
#% key: t
#% type: integer
#% required: no
#% multiple: no
#% description: Time interval for model step [days]
#% guisection: Dispersal parameters
#% options: 1-3285
#% answer: 30
#%end
#%option
#% key: p
#% type: double
#% required: no
#% multiple: no
#% description: Share of the stationary component (valid range 0 - 1)
#% guisection: Dispersal parameters
#% answer: 0.67
#%end
#%option
#% key: seed_fishmove
#% type: integer
#% required: no
#% multiple: no
#% description: fixed seed for calculating dispersal parameters using 'fishmove'
#% guisection: Dispersal parameters
#%end

#%option G_OPT_V_OUTPUT
#% key: output
#% description: Output vector map to store results
#% required: no
#% guisection: Output
#%end
#%option
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

#%flag
#% key: u
#% description: Values for source populations are provided in relative units (i.e. fish/probability per m)
#%end
#%flag
#% key: r
#% description: Calculate realisation of dispersal instead of probabilities
#%end
#%flag
#% key: k
#% description: Keep all temporal maps
#%end
#%flag
#% key: n
#% description: Only generate network (incl. barriers) and set up FIDIMO database. Don't calculate dispersal.
#%end
#%flag
#% key: s
#% description: Only update source populations. Network must have been generated before. Dispersal is not calculated.
#%end
#%flag
#% key: f
#% description: Only calculate dispersal based on previosely generated network and previousely defined source populations.
#%end
#%flag
#% key: m
#% description: Print out metadata only and exit
#%end

#%rules
#% requires_all: barriers,passability_col,threshold
#%end


###########################################
############# Import Libraries ############
###########################################

import os
import sys
import atexit
from grass.script import core as grass
import math
import sqlite3
import csv
import gc

import igraph
from igraph import *

from itertools import repeat

import datetime
from timeit import default_timer as timer

from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import vector as db
from grass.pygrass.modules.shortcuts import general as g

# import required numpy/scipy modules
import numpy
import scipy
from scipy import stats
from scipy import optimize

FIDIMO_version = "0.0"


###########################################
####### Set quite and overwrite ###########
###########################################
quiet = True
if grass.verbosity() > 2:
    quiet = False
overwrite = grass.overwrite()

###########################################
############# Define Cleanup ##############
###########################################


tmp_map_vect = None


def cleanup():
    if tmp_map_vect and not flags['k']:
        grass.run_command("g.remove",
                          flags='f',
                          type='vector',
                          name=[f + str(os.getpid()) for f in tmp_map_vect],
                          quiet=quiet)


###########################################
############# Set of Functions ############
###########################################

def import_vector(input_map,  # input vector name
                  output_map,  # output vector name
                  columns):  # dictionary {'target_name':'input_name'}
    ''' Import of vector maps for stream network and barriers.
    Copies selected input map to the specified output map. Only
    specified columns are copied and renamed if necessary'''
    
    grass.run_command("g.copy",
                      overwrite=True,
                      quiet=quiet,
                      vector="%s,%s" % (input_map, output_map))
    old_columns = grass.read_command("db.columns",
                                     table=output_map).splitlines()
    
    if sum([x not in old_columns for x in columns.values()]) > 0:
        #grass.fatal(_("At least one specified column does not exist in input map"))
        raise ValueError(
            "At least one specified column does not exist in input map")
        
    if sorted(columns.values() + ["cat"]) == sorted(old_columns):
        pass
    else:
        grass.run_command("v.db.dropcolumn",
                          quiet=quiet,
                          map=output_map,
                          columns=[x for x in old_columns if x not in columns.values() + ["cat"]])
        
    for i in columns:
        if columns[i] != i:
            grass.run_command("v.db.renamecolumn",
                              quiet=quiet,
                              map=output_map,
                              column=",".join([columns[i], i]))
            
    grass.run_command("v.db.addcolumn",
                      quiet=quiet,
                      map=output_map,
                      columns="orig_cat INTEGER")
    grass.run_command("v.db.update",
                      quiet=quiet,
                      map=output_map,
                      column="orig_cat",
                      query_column="cat")


def create_fidimo_db(fidimo_db_path):
    ''' Create FIDIMO DB'''
    
    # HERE a check if db already exists and launch an error if so and if overwrite flag is not set
    ########
    
    # If database exists it will be first removed
    try:
        os.remove(fidimo_db_path)
    except OSError:
        pass
        
    #grass.message(_("Creating FIDIMO Database and copying edges and vertices"))
    print("Creating FIDIMO Database")
    
    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()
    
    # Create table for meta data
    fidimo_db.execute(
        '''CREATE TABLE meta (parameter VARCHAR(45), value VARCHAR(300))''')
    fidimo_db.executemany("INSERT INTO meta (parameter) VALUES (?)",
                          [(x,)for x in ["Fidimo DB", "Projection", "Network name", "Total reaches n",
                                         "Barriers n", "Source populations n", "Distance matrix created", "Source populations imported",
                                         "Fidimo probabilities calculated", "Last modified", "GRASS setup", "FIDIMO version",
                                         "Scipy version", "Numpy version", "igraph version"]])
    
    # Update metadata
    fidimo_db.execute(
        '''UPDATE meta SET value=? WHERE parameter="Fidimo DB"''', (fidimo_db_path,))
    fidimo_db.execute(
        '''UPDATE meta SET value=? WHERE parameter="FIDIMO version"''', (FIDIMO_version,))
    fidimo_db.execute(
        '''UPDATE meta SET value=? WHERE parameter="Scipy version"''', (scipy.__version__,))
    fidimo_db.execute(
        '''UPDATE meta SET value=? WHERE parameter="Numpy version"''', (numpy.version.version,))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="GRASS setup"''', (", ".join(
        grass.read_command("g.version", flags="g").splitlines()[0:3]),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Projection"''',
                      (grass.read_command("g.proj", flags="jf")[:300],))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    fidimo_database.commit()
    fidimo_database.close()


def fidimo_network(input,
                   output,
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
    table and can be distinuished by the v_type column (1: barrier, 2: node)'''
    
    #grass.message(_("Preparing river network for FIDIMO"))
    print("Preparing river network for FIDIMO")
    
    # Input stream network
    streams_col_dict = {"strahler": strahler_col, "shreve": shreve_col,
                        "network": network_col}  # dictionary {'target_name':'input_name'}
    import_vector(input_map=input, columns=streams_col_dict,
                  output_map="streams_tmp" + str(os.getpid()))
    
    # Get original length of reaches
    grass.run_command("v.db.addcolumn",
                      quiet=quiet,
                      map="streams_tmp" + str(os.getpid()),
                      columns="orig_length DOUBLE")
    grass.run_command("v.to.db",
                      quiet=quiet,
                      map="streams_tmp" + str(os.getpid()),
                      columns="orig_length",
                      option="length",
                      units="meters")
    
    # Connect barriers (if present) and nodes to network
    if barriers:
        #grass.message(_("Connect barriers and nodes to network"))
        print("Connect barriers and nodes to network")
        # Input barriers
        # dictionary {'target_name':'input_name'}
        passability_col_dict = {"passability": passability_col}
        import_vector(input_map=barriers, columns=passability_col_dict,
                      output_map="barriers_tmp" + str(os.getpid()))
        
        # Get network ID for each barrier
        grass.run_command("v.db.addcolumn",
                          quiet=quiet,
                          map="barriers_tmp" + str(os.getpid()),
                          columns="network INTEGER")
        grass.run_command("v.what.vect",
                          map="barriers_tmp" + str(os.getpid()),
                          column="network",
                          query_map="streams_tmp" + str(os.getpid()),
                          query_column="network",
                          dmax=threshold)
        
        # Create network
        v.net(overwrite=True,
              quiet=quiet,
              flags="s",
              input="streams_tmp" + str(os.getpid()),
              points="barriers_tmp" + str(os.getpid()),
              output="fidimo_net1_tmp" + str(os.getpid()),
              operation="connect",
              threshold=threshold)
        v.net(overwrite=True,
              quiet=quiet,
              flags="c",
              input="fidimo_net1_tmp" + str(os.getpid()),
              output="fidimo_net2_tmp" + str(os.getpid()),
              operation="nodes")
    else:
        #grass.message(_("Connect nodes to network"))
        print("Connect nodes to network")
        v.net(overwrite=True,
              quiet=quiet,
              flags="c",
              input="streams_tmp" + str(os.getpid()),
              output="fidimo_net2_tmp" + str(os.getpid()),
              operation="nodes")
        
    # Get table with orig attributes for all (new split) network edges
    v.category(overwrite=True,
               quiet=quiet,
               input="fidimo_net2_tmp" + str(os.getpid()),
               layer=3,
               type="line",
               output="fidimo_net3_tmp" + str(os.getpid()),
               option="add")
    grass.run_command("v.db.addtable",
                      quiet=quiet,
                      map="fidimo_net3_tmp" + str(os.getpid()),
                      table="edges",
                      layer=3,
                      column="orig_cat INTEGER, strahler INTEGER, shreve INTEGER, network INTEGER, orig_length DOUBLE, edge_length DOUBLE, from_orig_v INTEGER, to_orig_v INTEGER, from_orig_e INTEGER, to_orig_e INTEGER")
    for i in ["orig_cat", "strahler", "shreve", "network", "orig_length"]:
        grass.run_command("v.to.db",
                          quiet=quiet,
                          map="fidimo_net3_tmp" + str(os.getpid()),
                          layer=3,
                          option="query",
                          columns=i,
                          query_layer=1,
                          query_colum=i)
        
    # Get lengths of (new split) network edges
    grass.run_command("v.to.db",
                      quiet=quiet,
                      map="fidimo_net3_tmp" + str(os.getpid()),
                      layer=3,
                      columns="edge_length",
                      option="length",
                      units="meters")
    
    # database-connection of current mapset
    mapset_db_settings = dict(x.split(": ") for x in grass.read_command("db.connect",
                                                                        flags="p").splitlines())
    if mapset_db_settings["driver"] != "sqlite":
        #grass.fatal(_("Database driver of current mapset is not 'sqlite'."))
        raise ValueError("Database driver of current mapset is not 'sqlite'.")
    mapset_database = sqlite3.connect(mapset_db_settings["database"])
    mapset_db = mapset_database.cursor()
    
    # Get for each edge the orig cat for the start (from_orig) and end point
    # (to_orig)
    mapset_db.execute('''CREATE TEMP TABLE edges_tmp 
              (cat INTEGER, from_orig_v INTEGER, to_orig_v INTEGER)''')
    e = [(int(x.split()[0]), int(x.split()[1]), int(x.split()[2])) for x in grass.read_command("v.net",
                                                                                               quiet=quiet,
                                                                                               input="fidimo_net3_tmp" +
                                                                                               str(os.getpid(
                                                                                               )),
                                                                                               operation="report",
                                                                                               arc_layer=3).splitlines()]
    mapset_db.executemany(
        "INSERT INTO edges_tmp (cat, from_orig_v, to_orig_v) VALUES (?,?,?)", e)
    mapset_db.execute('''UPDATE edges SET 
                from_orig_v = (SELECT from_orig_v FROM edges_tmp WHERE cat=edges.cat),
                to_orig_v = (SELECT to_orig_v FROM edges_tmp WHERE cat=edges.cat)
              WHERE EXISTS (SELECT cat FROM edges_tmp WHERE cat=edges.cat)''')
    mapset_database.commit()
    mapset_database.close()
    
    # Add table for vertices in layer 2
    grass.run_command("v.db.addtable",
                      quiet=quiet,
                      map="fidimo_net3_tmp" + str(os.getpid()),
                      table="vertices",
                      layer=2,
                      column="v_type INTEGER, orig_cat INTEGER")
    
    # Create and populate column to distinguish between barrier / node
    if barriers:
        grass.run_command("v.db.update",
                          quiet=quiet,
                          map="fidimo_net3_tmp" + str(os.getpid()),
                          layer=2,
                          column="v_type",
                          value=1,
                          where="cat IN (SELECT cat FROM barriers_tmp%s)" % (str(os.getpid()),))
        grass.run_command("v.db.update",
                          quiet=quiet,
                          map="fidimo_net3_tmp" + str(os.getpid()),
                          layer=2,
                          column="v_type",
                          value=2,
                          where="cat NOT IN (SELECT cat FROM barriers_tmp%s)" % (str(os.getpid()),))
    else:
        grass.run_command("v.db.update",
                          quiet=quiet,
                          map="fidimo_net3_tmp" + str(os.getpid()),
                          layer=2,
                          column="v_type",
                          value=2)
        
    # Get barrier orig_cat
    grass.run_command("v.db.update",
                      quiet=quiet,
                      map="fidimo_net3_tmp" + str(os.getpid()),
                      layer=2,
                      column="orig_cat",
                      query_column="cat")
    
    # Update network id for each barrier
    if barriers:
        grass.run_command("v.db.join",
                          quiet=quiet,
                          map="fidimo_net3_tmp" + str(os.getpid()),
                          layer=2,
                          column="orig_cat",
                          other_table="barriers_tmp" + str(os.getpid()),
                          other_column="orig_cat",
                          subset_columns="network")
    else:
        grass.run_command("v.db.addcolumn",
                      quiet=quiet,
                      map="fidimo_net3_tmp" + str(os.getpid()),
                      layer=2,
                      columns='''network INT''')
    
    # Copy final network to output map and adapt associated tables
    grass.run_command("v.extract",
                      quiet=quiet,
                      overwrite=overwrite,
                      flags="t",
                      input="fidimo_net3_tmp" + str(os.getpid()),
                      layer=3,
                      output="output_tmp" + str(os.getpid()))
    grass.run_command("v.category",
                      quiet=quiet,
                      overwrite=True,
                      input="output_tmp" + str(os.getpid()),
                      layer="3,1",
                      output=output,
                      option="chlayer")
    
    # Check if output table exists and remove if it exists
    tables_list = grass.read_command("db.tables", flags="p").splitlines()
    if "fidimo_output" in tables_list:
        grass.run_command("db.droptable",
                          quiet=quiet,
                          flags="f",
                          table="fidimo_output")
        
    grass.run_command("db.copy",
                      quiet=quiet,
                      overwrite=True,
                      from_table="edges",
                      to_table="fidimo_output")
    
    grass.run_command("v.db.connect",
                      quiet=quiet,
                      overwrite=True,
                      map=output,
                      layer=1,
                      table="fidimo_output")
    
    output_columns = grass.read_command("db.columns",
                                        table="fidimo_output").splitlines()
    grass.run_command("v.db.dropcolumn",
                      quiet=quiet,
                      map=output,
                      layer=1,
                      columns=[x for x in output_columns if x not in ["cat", "orig_cat"]])
    grass.run_command("v.db.addcolumn",
                      quiet=quiet,
                      map=output,
                      layer=1,
                      columns='''reach_length DOUBLE, source_pop DOUBLE, fidimo_result DOUBLE, fidimo_result_lwr DOUBLE, 
    fidimo_result_upr DOUBLE, rel_fidimo_result DOUBLE, rel_fidimo_result_lwr DOUBLE, rel_fidimo_result_upr DOUBLE''')
    
    #grass.message(_("Final networks prepared for FIDIMO"))
    print("Final networks prepared for FIDIMO")
    
    # Update metadata
    print("Update Metadata")
    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()
    fidimo_db.execute(
        '''UPDATE meta SET value=? WHERE parameter="Network name"''', (input,))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    if barriers:
        n_barriers = int(grass.read_command("v.info", flags="t",
                                            map="fidimo_net1_tmp" + str(os.getpid())).splitlines()[1].split("=")[1])
        fidimo_db.execute(
            '''UPDATE meta SET value=? WHERE parameter="Barriers n"''', (n_barriers,))
    else:
        fidimo_db.execute(
            '''UPDATE meta SET value=0 WHERE parameter="Barriers n"''')
    
    fidimo_database.commit()
    fidimo_database.close()


def set_fidimo_db(fidimo_db_path):
    ''' Create and tables for 
    vertices, edges, fidimo distance and barriers in Fidimo DB'''
    
    #grass.message(_("Creating FIDIMO Database and copying edges and vertices"))
    print("Copying edges and vertices to FIDIMO database")
    
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
    fidimo_db.execute(
        '''CREATE TABLE barriers_along (source INTEGER, target INTEGER, barrier INTEGER)''')
    
    # Update metadata
    print("Update Metadata")
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    fidimo_database.commit()
    fidimo_database.close()


def add_midpoints_fidimo_db(fidimo_db_path):
    '''Add midpoints (vertices) to each river reach to calculate
    distances between (midpoints of) river reaches.
    Each midpoint has its own id that is linked to the cat of the 
    river reach. Function is used in fidimo_distance'''
    
    #grass.message(_("Add midpoints (vertices) for each river reach"))
    print("Add midpoints (vertices) for each river reach")
    
    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()
    
    # Get midpoint vertices (reach id of grass edges) and insert in fidimo_db
    fidimo_db.execute('''SELECT cat FROM edges''')
    e_cat = [(x[0]) for x in fidimo_db.fetchall()]
    
    fidimo_db.execute('''SELECT cat FROM vertices''')
    v_cat = [(x[0]) for x in fidimo_db.fetchall()]
    
    # add midpoints with category values that start with the minimum value + 1 of the existing verices-cats
    # E.g. if the existing vertices have cats 1....5, then the midpoints cats
    # (v_type=3) will start with 6...
    fidimo_db.executemany("INSERT INTO vertices (cat,v_type,orig_cat) VALUES (?,?,?)",
                          zip(range(max(v_cat) + 1, max(v_cat) + 1 + len(e_cat)), [3] * len(e_cat), e_cat))
    
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
    
    # Get cat of vertices for from and to columns while splitting each reach
    # at midpoint into two parts
    fidimo_db.execute('''UPDATE edges SET from_v = (SELECT cat FROM vertices 
    WHERE vertices.cat=edges.from_orig_v AND vertices.v_type!=3) WHERE edges."part"=1;''')
    fidimo_db.execute('''UPDATE edges SET from_v = (SELECT cat FROM vertices 
  WHERE vertices.orig_cat=edges.cat AND vertices.v_type=3) WHERE edges."part"=2;''')
    
    fidimo_db.execute('''UPDATE edges SET to_v = (SELECT cat FROM vertices 
  WHERE vertices.orig_cat=edges.cat AND vertices.v_type=3) WHERE edges."part"=1;''')
    fidimo_db.execute('''UPDATE edges SET to_v = (SELECT cat FROM vertices 
  WHERE vertices.cat=edges.to_orig_v AND vertices.v_type!=3) WHERE edges."part"=2;''')
    
    # Add index to double-cat column
    fidimo_db.execute('''CREATE INDEX e_index_1 ON edges (part,cat)''')
    fidimo_db.execute('''CREATE INDEX e_index_2 ON edges (orig_cat)''')
    fidimo_db.execute('''CREATE INDEX e_index_3 ON edges (cat)''')
    
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
    #grass.message(_("Update main distance matrix (fidimo_distance) between all connected river reaches..."))
    print(
        "Update main distance matrix (fidimo_distance) between all connected river reaches...")
    
    fidimo_db.execute('SELECT DISTINCT network FROM edges')
    network_id = [x[0] for x in fidimo_db.fetchall()]
    
    for i in network_id:
        #grass.message(_("   ...processing network ID: %s" %str(i)))
        print("   ...processing network ID: %s" % str(i))
        
        # Fetch all edges
        fidimo_db.execute(
            'SELECT from_v, to_v, edge_length FROM edges WHERE network = ?', str(i))
        edges_attributes = fidimo_db.fetchall()
        edges = [(x[0], x[1]) for x in edges_attributes]
        
        # Get vertices and unique ids for graph (starting with 0)
        vertices = list(set([x[y] for x in edges_attributes for y in [0, 1]]))
        g_vertices = range(len(vertices))
        vertices_dict = dict(zip(vertices, g_vertices))
        inv_vertices_dict = {v: k for k, v in vertices_dict.iteritems()}
        
        # Get new unique igraph ids for vertices in edges list
        g_edges = [(vertices_dict[x[0]], vertices_dict[x[1]])
                   for x in edges_attributes]
                   
        # Fetch all vertices for specific network that are midpoints of river
        # reaches
        fidimo_db.execute(
            'SELECT cat FROM vertices WHERE v_type=3 AND orig_cat IN (SELECT cat from edges WHERE network = ? )', str(i))
        midpoints = [x[0] for x in fidimo_db.fetchall()]
        g_midpoints = [vertices_dict[x] for x in midpoints]
        
        # Build Graph
        g = Graph(g_edges, directed=True)
        
        # Add length to edges (half length as each edge is split by the
        # midpoint)
        g.es["half_length"] = [(x[2]) / 2.0 for x in edges_attributes]
        
        # Create chunks for memory-saving processing        
        chunk_size_midpoints = int(100 * round(float(10E6/len(midpoints))/100))
        g_midpoints_chunks = [g_midpoints[x:x + chunk_size_midpoints]
                              for x in xrange(0, len(g_midpoints), chunk_size_midpoints)]
           
        # Calculate shortest pathss
        #grass.message(_("Calculating shortest paths"))
        print("Calculating shortest paths (chunk size: "+str(chunk_size_midpoints)+")...")
                                      
        for k in range(len(g_midpoints_chunks)):
            print("...chunk: "+str(k+1)+" of "+str(len(g_midpoints_chunks)))
            distance_mat_upstream = g.shortest_paths(source=g_midpoints_chunks[k],
                                                     target=g_midpoints, weights="half_length", mode="IN")
            distance_mat_downstream = g.shortest_paths(source=g_midpoints_chunks[k],
                                                       target=g_midpoints, weights="half_length", mode="OUT")
            distance_mat_all = g.shortest_paths(source=g_midpoints_chunks[k],
                                                target=g_midpoints, weights="half_length", mode="ALL")
            
            # Distinguish between paths up- and downstream [1: downstream, 2:upstream, 3:neither (down-up combination), 4: both directions (e.g. where source=target and dist=0)
            #grass.message(_("Calculating directions between reaches"))
            print("Calculating directions between reaches")
            
            direction_mat = numpy.select([numpy.isinf(distance_mat_upstream) & ~numpy.isinf(distance_mat_downstream),
                                          numpy.isinf(distance_mat_downstream) & ~numpy.isinf(
                                              distance_mat_upstream),
                                          numpy.isinf(distance_mat_upstream) & numpy.isinf(
                                              distance_mat_downstream),
                                          ~numpy.isinf(distance_mat_upstream) & ~numpy.isinf(distance_mat_downstream)],
                                         [1, 2, 3, 4], default=-99999)
            
            # Collect garbage to free memory
            del distance_mat_upstream
            del distance_mat_downstream
            gc.collect()
            
            #grass.message(_("Updating distance matrix for specific network"))
            print("Updating distance matrix for specific network and chunk")
            paths_to_db = zip([x for item in g_midpoints_chunks[k] for x in repeat(item, len(midpoints))],  # list of all source/from midpoints of reaches
                              # list of all target/to midpoints of reaches
                              midpoints * len(g_midpoints_chunks[k]),
                              # list of distances
                              [item for sublist in distance_mat_all for item in sublist],
                              # list of directions
                              [item for sublist in direction_mat for item in sublist],
                              [i] * (len(midpoints)*len(g_midpoints_chunks[k])))  # network
            print("Updating distance matrix for specific network and chunk 2")
            fidimo_db.executemany(
                "INSERT INTO fidimo_distance (source,target,distance,direction,network) VALUES (?,?,?,?,?)", paths_to_db)
            fidimo_database.commit()
                        
            # Collect garbage to free memory
            del paths_to_db
            del direction_mat
            gc.collect()
        
        # Barriers along each path (if exist)
        fidimo_db.execute(
            '''SELECT cat FROM vertices WHERE v_type=1 and network=?''', str(i))
        barriers_list = [x[0] for x in fidimo_db.fetchall()]
        
        if barriers_list:
            #grass.message(_("Processing barriers along networks paths"))
            print("Processing barriers along networks paths")
            
            g_barriers = [vertices_dict[x] for x in barriers_list]
            
            for j in xrange(len(midpoints)):
                # j=0
                vertices_paths = g.get_shortest_paths(v=g_midpoints[j], to=g_midpoints,
                                                      weights="half_length", output="vpath", mode="ALL")
                paths_with_barriers = [
                    [x for x in L if x in g_barriers] for L in vertices_paths]
                barriers_along = [(midpoints[j], inv_vertices_dict[b], inv_vertices_dict[
                                   e]) for a, b in zip(paths_with_barriers, g_midpoints) for e in a]
                
                fidimo_db.executemany(
                    '''INSERT INTO barriers_along (source,target,barrier) VALUES (?,?,?)''', barriers_along)
                
    #grass.message(_("All networks processed"))
    print("All networks processed")
    
    # Update values for main fidimo table (e.g. lengths, stream order,...)
    start = timer()
    
    # Get number of rows of fidimo_distance and chunk it into pieces of 10E5
    fidimo_db.execute('''SELECT max(rowid) FROM fidimo_distance''')
    max_fidimo_distance_rowid = [x[0] for x in fidimo_db.fetchall()][0]
    fidimo_distance_rowid_chunks = [[x+1,x + 10E5] for x in xrange(0, max_fidimo_distance_rowid, int(10E5))]
    
    print("Updating fidimo_distance database...")
    for k in range(len(fidimo_distance_rowid_chunks)):
        print("...chunk: "+str(k+1)+" of "+str(len(fidimo_distance_rowid_chunks)))
        
        # Get cats of original vertices (from-to), Get cats of original edges (from-to)
        #grass.message(_("Updating original vertex categories to fidimo_distance"))
        print("Updating original vertex categories to fidimo_distance")
        fidimo_db.execute('''UPDATE fidimo_distance SET 
                      from_orig_v = (SELECT orig_cat FROM vertices WHERE cat=fidimo_distance.source),
                      to_orig_v = (SELECT orig_cat FROM vertices WHERE cat=fidimo_distance.target)
                      WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
        
        print("Updating original river reach (edges) categories to fidimo_distance")
        fidimo_db.execute('''UPDATE fidimo_distance SET 
                      from_orig_e = (SELECT orig_cat FROM edges WHERE cat=fidimo_distance.from_orig_v),
                      to_orig_e = (SELECT orig_cat FROM edges WHERE cat=fidimo_distance.to_orig_v)
                      WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
      
        # Get edge lengths, stream order and network id for source (and target) reach
        #grass.message(_("Updating addtio nal attributes (e.g. stream order) in fidimo_distance"))
        print(
            "Updating addtional attributes (e.g. stream order) in fidimo_distance")
        fidimo_db.execute('''UPDATE fidimo_distance SET
                    source_edge_length = (SELECT edge_length FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1),
                    target_edge_length = (SELECT edge_length FROM edges WHERE cat=fidimo_distance.to_orig_v AND part=1),
                    source_strahler = (SELECT strahler FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1),
                    source_shreve = (SELECT shreve FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1),
                    target_shreve = (SELECT shreve FROM edges WHERE cat=fidimo_distance.to_orig_v AND part=1),
                    network = (SELECT network FROM edges WHERE cat=fidimo_distance.from_orig_v AND part=1)
                    WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
        
        # Get upr and lwr limit for later integration based on the reach lengths and distance
        #grass.message(_("Updating addtional attributes (upper and lower distance limit) in fidimo_distance"))
        print(
            "Updating addtional attributes (upper and lower distance limit) in fidimo_distance")
        fidimo_db.execute('''UPDATE fidimo_distance SET
                    lwr_limit = distance-(target_edge_length/2),
                    upr_limit = distance+(target_edge_length/2)
                    WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
        
    end = timer()
    
    #grass.message(_("Time elapsed: %s" %str(end-start)))
    print("Time elapsed: %s" % str(end - start))
    
    # Update metadata
    print("Update Metadata")
    fidimo_db.execute('SELECT COUNT(*) FROM vertices WHERE v_type==3')
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Total reaches n"''', (str(
        [x[0] for x in fidimo_db.fetchall()][0]),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Distance matrix created"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    fidimo_db.execute(
        '''UPDATE meta SET value=? WHERE parameter="igraph version"''', (igraph.__version__,))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    fidimo_database.commit()
    fidimo_database.close()


def sigma_calc(l,
               ar,
               t,
               statistical_interval,
               seed_fishmove=None):
    '''This function calculates dispersal distances sigma_stat and sigma_mob
    for each stream order and for given input: l, ar, t'''

    # Regression model to calculate dispersal kernel parameter sigma_stat and
    # sigma_mob
    if statistical_interval == "no":
        # log(sigma_stat) = -10.57 + 1.64*log(L) + 0.94*AR + 1.14*sqrt(SO) + 0.43*log(T)
        # log(sigma_mob) = -7.48 + 1.45*log(L) + 0.58*AR + 1.51*sqrt(SO) + 0.55*log(T)
        # Parameter coefficients calculated from
        # fishmove(L=200,AR=1.25,rep=5000,seed=999)
        sigma_dict = {
            "stat": {
                "fit": dict(zip(range(1, 10), [math.exp(-10.56720 + 1.643237 * math.log(l) + 0.9647056 * ar + 1.142661 * math.sqrt(i) + 0.4273618 * math.log(t)) for i in (1, 2, 3, 4, 5, 6, 7, 8, 9)]))},
            "mob": {
                "fit": dict(zip(range(1, 10), [math.exp(-7.480260 + 1.445552 * math.log(l) + 0.5820339 * ar + 1.507528 * math.sqrt(i) + 0.5527266 * math.log(t)) for i in (1, 2, 3, 4, 5, 6, 7, 8, 9)]))}}

        return sigma_dict

    else:
        import rpy2.robjects as robjects  # import required rpy2 module
        from rpy2.robjects.packages import importr
        fm = importr('fishmove')

        # Calculating 'fishmove' depending on species or L & AR
        # Statistical interval
        if "prediction" in statistical_interval:
            interval = "prediction"
        elif "confidence" in statistical_interval:
            interval = "confidence"

        # Set fixed seed if specified
        if seed_fishmove:
            seed = ",seed=" + str(seed_fishmove)
        else:
            seed = ""

        so_rvector = robjects.IntVector((1, 2, 3, 4, 5, 6, 7, 8, 9))

        # Calculate 'fishmove' and store sigma values in pandas df
        fishmove = eval(
            "fm.fishmove(L=l,AR=ar,SO=so_rvector,T=t,interval=interval,rep=200%s)" % (seed))
        fishmove = fishmove[1]
        sigma_dict = {
            "stat": {
                "fit": dict(zip(range(1, 10), list(fishmove.rx("fit", 'sigma_stat', 1, 1, so_rvector, 1)))),
                "lwr": dict(zip(range(1, 10), list(fishmove.rx("lwr", 'sigma_stat', 1, 1, so_rvector, 1)))),
                "upr": dict(zip(range(1, 10), list(fishmove.rx("upr", 'sigma_stat', 1, 1, so_rvector, 1))))},
            "mob": {"fit": dict(zip(range(1, 10), list(fishmove.rx("fit", 'sigma_mob', 1, 1, so_rvector, 1)))),
                    "lwr": dict(zip(range(1, 10), list(fishmove.rx("lwr", 'sigma_mob', 1, 1, so_rvector, 1)))),
                    "upr": dict(zip(range(1, 10), list(fishmove.rx("upr", 'sigma_mob', 1, 1, so_rvector, 1))))}}

        return sigma_dict


def fidimo_kernel_cdf_truncation(x, sigma_stat, sigma_mob, truncation, p):
    '''This function calculates the 
    maximum distance (cutting distance) 
    based on truncation criterion.'''
    return p * stats.norm.cdf(x, loc=0, scale=sigma_stat) + (1 - p) * stats.norm.cdf(x, loc=0, scale=sigma_mob) - truncation


def fidimo_kernel_cdf(x, sigma_stat, sigma_mob, p):
    '''This function calculates the dispersal probability 
    based on leptokurtic dispersal kernels for river fish (Radinger and Wolter, 2014).
    Dispersal probability is calculated based on the cumulative density function
    between lower and upper distance.'''
    return (p * stats.norm.cdf(x, loc=0, scale=sigma_stat) + (1 - p) * stats.norm.cdf(x, loc=0, scale=sigma_mob))


def fidimo_source_pop( source_pop_csv,
                       fidimo_db_path,
                       realisation):
    '''This function appends source populations to distance matrix and creates output map'''
    
    # connect to database
    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()
    
    # Delete fidimo_source_pop if exists
    fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_source_pop_tmp''')
    fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_source_pop''')
    
    fidimo_db.execute("CREATE TABLE fidimo_source_pop_tmp (cat, source_pop, p);")
    
    # Read CSV and check if three columns (reach cat, source_col, p)
    with open(source_pop_csv,'rb') as csv_file:
        source_pop_csv_read = csv.DictReader(csv_file,fieldnames=("cat", "source_pop", "p")) # comma is default delimiter,header is assumed 
        # Skip first header line
        next(source_pop_csv_read)
        # csv to list of tuples for input into db
        source_pop_csv_read_to_db = [(int(i['cat']), float(i['source_pop']), float(i['p'])) for i in source_pop_csv_read]
  
    fidimo_db.executemany("INSERT INTO fidimo_source_pop_tmp (cat, source_pop, p) VALUES (?, ?, ?);", source_pop_csv_read_to_db)
    fidimo_database.commit()
    
    fidimo_db.execute(
        '''CREATE INDEX fidimo_source_pop_tmp_index ON fidimo_source_pop_tmp (cat)''')
    
    # Check if source pop and p are within a valid range
    fidimo_db.execute('SELECT DISTINCT p FROM fidimo_source_pop_tmp')
    p_fidimo_source_pop_tmp = [x[0] for x in fidimo_db.fetchall()]
    if (min(p_fidimo_source_pop_tmp)<0) or (max(p_fidimo_source_pop_tmp)>1):
        raise ValueError(
            "Values for p must be decimal numbers between 0 and 1. NAs not allowed")
    fidimo_db.execute('SELECT DISTINCT source_pop FROM fidimo_source_pop_tmp')
    source_pop_fidimo_source_pop_tmp = [x[0] for x in fidimo_db.fetchall()]
    if min(source_pop_fidimo_source_pop_tmp)<0:
        raise ValueError(
            "Values for source populations must be positive integer or decimal numbers. NAs not allowed")
        
    # check if cats of fidimo_source_pop_tmp == orig_cat of edges
    fidimo_db.execute('SELECT DISTINCT orig_cat FROM edges')
    orig_cat_edges = [x[0] for x in fidimo_db.fetchall()]
    fidimo_db.execute('SELECT DISTINCT cat FROM fidimo_source_pop_tmp')
    cat_fidimo_source_pop_tmp = [x[0] for x in fidimo_db.fetchall()]
    if sorted(orig_cat_edges) != sorted(cat_fidimo_source_pop_tmp):
        #grass.fatal(_("Vector input map of source populations must must match vector input map that has been used for calculating fidimo distance matrix"))
        raise ValueError(
            "IDs (categories, cat) of source populations must must match cats of vector input map that has been used for calculating fidimo distance matrix")
    else:
        print(
            "IDs (categories, cat) of source populations match cats of vector input map that has been used for calculating fidimo distance matrix")
        
    # Create source_pop table in FIDIMO DB
    fidimo_db.execute('''CREATE TABLE fidimo_source_pop AS SELECT
              cat AS cat, orig_cat AS orig_cat, edge_length AS edge_length
              FROM edges WHERE edges.part=1;''')
    fidimo_db.execute(
        '''ALTER TABLE fidimo_source_pop ADD COLUMN source_pop DOUBLE''')
    
    # Update source populations from input
    fidimo_db.execute(
        'UPDATE fidimo_source_pop SET source_pop = (SELECT source_pop FROM fidimo_source_pop_tmp WHERE orig_cat=fidimo_source_pop_tmp.cat)')
    fidimo_database.commit()
    
    # Delete temporary input table fidimo_source_pop_tmp and corresponding
    # index if exists
    fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_source_pop_tmp''')
    fidimo_db.execute('''DROP INDEX IF EXISTS fidimo_source_pop_tmp_index''')
    
    # If source populations are provided in relative units (e.g. CPUE, density per meter, probability per m)
    # then caculate first absolute values
    # if flags['u']:
    # fidimo_db.execute('UPDATE fidimo_source_pop SET source_pop = source_pop*edge_length')
    
    # Correct for barrier splits
    # Get edges/reaches that are split by barriers (e.g. with more than 2
    # entries for orig_cat)
    fidimo_db.execute(
        'SELECT orig_Cat FROM fidimo_source_pop GROUP BY orig_cat HAVING COUNT(*) > 1;')
    split_reach_cats = [x[0] for x in fidimo_db.fetchall()]
    for i in split_reach_cats:
        if realisation == True:
            fidimo_db.execute(
                'UPDATE fidimo_source_pop SET source_pop = round((source_pop/(SELECT SUM(edge_length) FROM fidimo_source_pop WHERE orig_cat=?))*edge_length,0)  WHERE orig_cat=?', (i, i))
        else:
            fidimo_db.execute(
                'UPDATE fidimo_source_pop SET source_pop = (source_pop/(SELECT SUM(edge_length) FROM fidimo_source_pop WHERE orig_cat=?))*edge_length  WHERE orig_cat=?', (i, i))
            
    # Setting all reaches that are not defined as source populations to 0
    fidimo_db.execute(
        '''UPDATE fidimo_source_pop SET source_pop=0.0 WHERE source_pop IS NULL OR source_pop='';''')
    fidimo_database.commit()
    
    # First check if source_pop already exists in fidimo distance and create
    # if not exist
    fidimo_db.execute('SELECT * FROM fidimo_distance LIMIT 1')
    if "source_pop" not in [x[0] for x in fidimo_db.description]:
        fidimo_db.execute(
            '''ALTER TABLE fidimo_distance ADD COLUMN source_pop DOUBLE''')
    
    # Get number of rows of fidimo_distance and chunk it into pieces of 10E5
    fidimo_db.execute('''SELECT max(rowid) FROM fidimo_distance''')
    max_fidimo_distance_rowid = [x[0] for x in fidimo_db.fetchall()][0]
    fidimo_distance_rowid_chunks = [[x+1,x + 10E5] for x in xrange(0, max_fidimo_distance_rowid, int(10E5))]   
    
    # Join table fidimo_distance with fidimo_source_pop
    #grass.message(_("Updating source population in fidimo_distance"))
    print("Updating source population in fidimo_distance")
    for k in range(len(fidimo_distance_rowid_chunks)):
        print("...chunk: "+str(k+1)+" of "+str(len(fidimo_distance_rowid_chunks))) 
        fidimo_db.execute('''UPDATE fidimo_distance SET 
                      source_pop = (SELECT source_pop FROM fidimo_source_pop WHERE cat=fidimo_distance.from_orig_v)
                      WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
    
    # Update metadata
    print("Update Metadata")
    fidimo_db.execute(
        'SELECT COUNT(*) FROM fidimo_source_pop WHERE source_pop>0')
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Source populations n"''', (str(
        [x[0] for x in fidimo_db.fetchall()][0]),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Source populations imported"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    fidimo_database.commit()
    fidimo_database.close()


def fidimo_probability( fidimo_db_path,
                        sigma_dict,
                        p,
                        statistical_interval,
                        truncation="inf"):
    '''This function calculates dispersal kernel probabilities
    for each river reach and corresponding distances between reaches'''

    # connect to database
    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()

    # Get stream orders (ASC) of source populations where source_pop > 0
    fidimo_db.execute(
        '''SELECT DISTINCT source_strahler FROM fidimo_distance WHERE source_pop > 0''')
    so_list = sorted([x[0] for x in fidimo_db.fetchall()])
    if None in so_list:
        raise ValueError(
            "At least one stream reach with a source population has no value for Strahler stream order")
    elif max([int(x) for x in so_list]) > 9:
        #grass.fatal(_("Stream network has stream orders > 9. Please consider smaller stream network"))
        raise ValueError(
            "Stream network has stream orders > 9. Please consider smaller stream network")
    else:
        so_list = [int(x) for x in so_list]

    # First check if fidimo_prob already exists in fidimo distance and create
    # if not exist
    fidimo_db.execute('SELECT * FROM fidimo_distance LIMIT 1')
    if "basic_fidimo_prob" not in [x[0] for x in fidimo_db.description]:
        fidimo_db.execute(
            '''ALTER TABLE fidimo_distance ADD COLUMN basic_fidimo_prob DOUBLE''')
        fidimo_db.execute(
            '''ALTER TABLE fidimo_distance ADD COLUMN basic_fidimo_prob_lwr DOUBLE''')
        fidimo_db.execute(
            '''ALTER TABLE fidimo_distance ADD COLUMN basic_fidimo_prob_upr DOUBLE''')

    # Before any calculations set fidimo_prob in fidimo_distance to ''
    fidimo_db.execute(
        '''UPDATE fidimo_distance SET basic_fidimo_prob=NULL,basic_fidimo_prob_lwr=NULL,basic_fidimo_prob_upr=NULL;''')
    fidimo_database.commit()

    # Delete fidimo_prob if exists
    fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_prob''')

    # Get number of runs for statistical intervals
    if statistical_interval == "no":
        nrun = ["fit"]
    else:
        nrun = ["fit", "lwr", "upr"]

    # Create a new sqlite table in fidimo_db to collect results
    fidimo_db.execute('''CREATE TABLE fidimo_prob ( fidimo_distance_id INTEGER, 
                          basic_fidimo_prob DOUBLE, 
                          basic_fidimo_prob_lwr DOUBLE, 
                          basic_fidimo_prob_upr DOUBLE)''')

    # Commit changes
    fidimo_database.commit()

    #grass.debug(_("Begin with core of fidimo, application of fishmove"))
    print("Begin with core of fidimo, application of fishmove")

    for i in so_list:
        # Calculate maximum distance (cutting distance) based on truncation
        # criterion (based on p=0.67)
        if truncation != "inf":
            max_dist = float(optimize.zeros.newton(fidimo_kernel_cdf_truncation, 1.,
                                                   args=(sigma_dict["stat"]["fit"][i], sigma_dict["mob"]["fit"][i], float(truncation), 0.67)))
        else:
            max_dist = 1e8  # Max distance = 100000 km

        direction_dict = {"downstream": 1,
                          "upstream": 2, "undef": 3, "source": 4}

        for j in direction_dict:  # loop over directions (up- vs. downstream)
            ''' calculate fidimo probability in chunks based on source_strahler
            MAIN PART: leptokurtic probability density kernel based on fishmove '''

            # Do not calculate any fidimo probability for direction that are
            # combinations of down- and upstream (as in r.fidimo)
            if j == "undef":
                continue

            print("Calculation of Fidimo probability for direction: " +
                  j + " and stream order: " + str(i))

            # Get all IDs for cases where distance, direction, source_pop and
            # source_strahler match criteria
            fidimo_db.execute('''SELECT fidimo_distance_id FROM fidimo_distance 
        WHERE source_strahler = ?
          AND distance <= ? 
          AND direction = ?
          AND source_pop > 0''', (i, max_dist, direction_dict[j]))

            # Split all cases into chunks of max. chunk_size
            row_ids = [x[0] for x in fidimo_db.fetchall()]
            chunk_size = 500
            row_ids_chunks = [row_ids[x:x + chunk_size]
                              for x in xrange(0, len(row_ids), chunk_size)]

            for k in range(len(row_ids_chunks)):
                print("Processing chunk: " + str(k + 1) +
                      " of " + str(len(row_ids_chunks)))
                chunk = row_ids_chunks[k]
                fidimo_db.execute('SELECT * FROM fidimo_distance WHERE fidimo_distance_id IN (%s)' %
                                  ','.join('?' * len(chunk)), tuple(chunk))
                fidimo_distance_colnames = dict(
                    zip([x[0] for x in fidimo_db.description], range(0, len(fidimo_db.description))))
                fidimo_distance_array = scipy.array(fidimo_db.fetchall())

                # Create results array to collect fidimo_prob results
                fidimo_result_array = scipy.array(chunk)

                for l in nrun:
                    # CDF(upr) - CDF(lwr)
                    if j == "source":
                        basic_fidimo_prob = (fidimo_kernel_cdf(x=fidimo_distance_array[:, fidimo_distance_colnames["upr_limit"]].astype(float),
                                                               sigma_stat=sigma_dict[
                            "stat"][l][i],
                            sigma_mob=sigma_dict[
                            "mob"][l][i],
                            p=p) - fidimo_kernel_cdf( x=0,
                                                      sigma_stat=sigma_dict[
                                                          "stat"][l][i],
                                                      sigma_mob=sigma_dict[
                                                          "mob"][l][i],
                                                      p=p)) * 2.0
                    elif j == "upstream" or j == "downstream":
                        basic_fidimo_prob = fidimo_kernel_cdf(x=fidimo_distance_array[:, fidimo_distance_colnames["upr_limit"]].astype(float),
                                                              sigma_stat=sigma_dict[
                            "stat"][l][i],
                            sigma_mob=sigma_dict[
                            "mob"][l][i],
                            p=p) - fidimo_kernel_cdf( x=fidimo_distance_array[:, fidimo_distance_colnames["lwr_limit"]].astype(float),
                                                      sigma_stat=sigma_dict[
                                                          "stat"][l][i],
                                                      sigma_mob=sigma_dict[
                                                          "mob"][l][i],
                                                      p=p)

                    # if upstream direction than correct for network splits by
                    # a shreve stream order approach
                    if j == "upstream":
                        basic_fidimo_prob = basic_fidimo_prob * (fidimo_distance_array[:, fidimo_distance_colnames["target_shreve"]].astype(
                            float) / fidimo_distance_array[:, fidimo_distance_colnames["source_shreve"]].astype(float))

                    # Set all values that are 0 (i.e. to small to be calculated
                    # to a the minimum float number possible)
                    basic_fidimo_prob[
                        basic_fidimo_prob == 0] = sys.float_info.min

                    # Stack fidimo_prob with fidimo_distance id into result
                    # array
                    fidimo_result_array = numpy.column_stack(
                        (fidimo_result_array,  # fidimo_distance_id
                         basic_fidimo_prob))  # fidimo_prob for single nrun

                if statistical_interval == "no":
                    fidimo_db.executemany("INSERT INTO fidimo_prob (fidimo_distance_id,basic_fidimo_prob) VALUES (?,?)", map(
                        tuple, fidimo_result_array.tolist()))
                else:
                    fidimo_db.executemany("INSERT INTO fidimo_prob (fidimo_distance_id,basic_fidimo_prob,basic_fidimo_prob_lwr,basic_fidimo_prob_upr) VALUES (?,?,?,?)", map(
                        tuple, fidimo_result_array.tolist()))

                # Commit changes
                fidimo_database.commit()

    # Create key in final fidimo_prob
    fidimo_db.execute(
        '''CREATE INDEX fidimo_prob_index ON fidimo_prob (fidimo_distance_id)''')

    # Join fidimo_prob with fidimo_distance
    print(
        "Update fidimo_distance with basic (unweighted) fidimo probabilities")
    fidimo_db.execute('''UPDATE fidimo_distance SET 
                basic_fidimo_prob = (SELECT basic_fidimo_prob FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                basic_fidimo_prob_lwr = (SELECT basic_fidimo_prob_lwr FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                basic_fidimo_prob_upr = (SELECT basic_fidimo_prob_upr FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id)
                 WHERE EXISTS (SELECT fidimo_distance_id FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id)''')

    print("Update Metadata")
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Fidimo probabilities calculated"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))

    # Commit changes
    fidimo_database.commit()

    # Close database connection
    fidimo_database.close()


def fidimo_realisation( realisation,
                        fidimo_db_path):
    '''This function either weights the fidimo probabilities with the value of the value of
    the source population or calculates realised fish count (real integer values) from the
    probabilities using the multinomial distribution'''

    # connect to database
    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()

    # First check if fidimo_result already exists in fidimo distance and
    # create if not exist
    fidimo_db.execute('SELECT * FROM fidimo_distance LIMIT 1')
    if "fidimo_result" not in [x[0] for x in fidimo_db.description]:
        fidimo_db.execute(
            '''ALTER TABLE fidimo_distance ADD COLUMN fidimo_result DOUBLE''')
        fidimo_db.execute(
            '''ALTER TABLE fidimo_distance ADD COLUMN fidimo_result_lwr DOUBLE''')
        fidimo_db.execute(
            '''ALTER TABLE fidimo_distance ADD COLUMN fidimo_result_upr DOUBLE''')

    # Before any calculations set fidimo_result in fidimo_distance to ''
    fidimo_db.execute(
        '''UPDATE fidimo_distance SET fidimo_result=NULL,fidimo_result_lwr=NULL,fidimo_result_upr=NULL;''')
    fidimo_database.commit()

    if realisation == True:
        print(
            "Calculate realisation (i.e. fish counts that disperse from each source population)")

        # Get number of runs for statistical intervals
        fidimo_db.execute('''SELECT SUM(basic_fidimo_prob) AS basic_fidimo_prob,
                  SUM(basic_fidimo_prob_lwr) AS basic_fidimo_prob_lwr, 
                  SUM(basic_fidimo_prob_upr) AS basic_fidimo_prob_upr 
              FROM fidimo_distance''')
        nrun = zip([x[0] for x in fidimo_db.description], [
                   x != None for x in fidimo_db.fetchall()[0]])

        # Realisation using the mutlinomial distribution to obtain indiviual counts
        # Get ids of all source reaches that have source_pop>0
        fidimo_db.execute(
            'SELECT DISTINCT from_orig_v,source_pop FROM fidimo_distance WHERE source_pop>0')
        source_populations = [[x[0], x[1]] for x in fidimo_db.fetchall()]

        # Insert results of realisation in temporary table and then update
        # fidimo_distance from that table
        fidimo_db.execute('''CREATE TEMP TABLE realisation_result_tmp 
              (fidimo_distance_id INTEGER, fidimo_result DOUBLE, fidimo_result_lwr DOUBLE, fidimo_result_upr DOUBLE)''')

        for i in source_populations:
            fidimo_db.execute(
                '''SELECT * FROM fidimo_distance WHERE from_orig_v = ?''', (i[0],))
            fidimo_distance_colnames = dict(
                zip([x[0] for x in fidimo_db.description], range(0, len(fidimo_db.description))))
            fidimo_distance_array = scipy.array(fidimo_db.fetchall())
            # Array to collect results
            realised_fidimo_result = fidimo_distance_array[
                :, fidimo_distance_colnames["fidimo_distance_id"]].astype(int)
            for j in nrun:
                if j[1] == False:
                    continue
                basic_fidimo_prob = numpy.nan_to_num(
                    fidimo_distance_array[:, fidimo_distance_colnames[j[0]]].astype(float))

                # if options['seed2']:
                # numpy.random.seed(int(optionss['seed2']))
                realised_fidimo_result_i = numpy.random.multinomial(int(i[1]),
                                                                    (basic_fidimo_prob / numpy.sum(basic_fidimo_prob)))

                realised_fidimo_result = numpy.column_stack((realised_fidimo_result,
                                                             realised_fidimo_result_i.astype(int)))

            if sum([x[1] for x in nrun]) == 1:
                fidimo_db.executemany('''INSERT INTO realisation_result_tmp 
            (fidimo_distance_id, fidimo_result) 
            VALUES (?,?)''', map(tuple, realised_fidimo_result.tolist()))
            elif sum([x[1] for x in nrun]) == 3:
                fidimo_db.executemany('''INSERT INTO realisation_result_tmp 
            (fidimo_distance_id, fidimo_result, fidimo_result_lwr, fidimo_result_upr) 
            VALUES (?,?,?,?)''', map(tuple, realised_fidimo_result.tolist()))
            else:
                print(
                    "Realisation cannot be calculated due to some statistical intervals")

        # Join fidimo_prob with fidimo_distance
        print(
            "Update fidimo_distance with realised fidimo results (i.e. fish counts)")
        fidimo_db.execute('''UPDATE fidimo_distance SET 
                fidimo_result = (SELECT fidimo_result FROM realisation_result_tmp WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                fidimo_result_lwr = (SELECT fidimo_result_lwr FROM realisation_result_tmp WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                fidimo_result_upr = (SELECT fidimo_result_upr FROM realisation_result_tmp WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id)
                 WHERE EXISTS (SELECT fidimo_distance_id FROM realisation_result_tmp WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id)''')

    else:
        # Multiply by value of inital source population
        print(
            "Update fidimo_distance with source-population-weighted fidimo probability")
        fidimo_db.execute('''UPDATE fidimo_distance SET
                  fidimo_result = basic_fidimo_prob*source_pop,
                  fidimo_result_lwr = basic_fidimo_prob_lwr*source_pop,
                  fidimo_result_upr = basic_fidimo_prob_upr*source_pop''')

    # Update metadata
    print("Update Metadata")
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))

    # Commit changes
    fidimo_database.commit()

    # Close database connection
    fidimo_database.close()


def fidimo_summarize( output,
                      fidimo_db_path):
    '''This function summarizes the fidimo result for each target reach
    and writes results back to output vector map'''

    # connect to database
    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()

    # First, delete summary_fidimo_prob if exists
    fidimo_db.execute('''DROP TABLE IF EXISTS summary_fidimo_result''')

    # Summarize fidimo prob for each target reach
    fidimo_db.execute('''CREATE TABLE summary_fidimo_result AS SELECT
    to_orig_v, 
    target_edge_length AS reach_length,
    CAST(SUM(fidimo_result) AS DOUBLE) AS fidimo_result, 
    CAST(SUM(fidimo_result_lwr) AS DOUBLE) AS fidimo_result_lwr,
    CAST(SUM(fidimo_result_upr) AS DOUBLE) AS fidimo_result_upr 
    FROM fidimo_distance GROUP BY to_orig_v;''')

    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN source_pop DOUBLE''')
    fidimo_db.execute('''UPDATE summary_fidimo_result SET
                  source_pop = (SELECT source_pop FROM fidimo_source_pop WHERE cat=summary_fidimo_result.to_orig_v)''')
    # Commit changes
    fidimo_database.commit()

    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN rel_fidimo_result DOUBLE''')
    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN rel_fidimo_result_lwr DOUBLE''')
    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN rel_fidimo_result_upr DOUBLE''')

    fidimo_db.execute('''UPDATE summary_fidimo_result SET
                  rel_fidimo_result = fidimo_result/reach_length,
                  rel_fidimo_result_lwr = fidimo_result_lwr/reach_length,
                  rel_fidimo_result_upr = fidimo_result_upr/reach_length''')
    # Update metadata
    print("Update Metadata")
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))

    # Commit changes
    fidimo_database.commit()

    # Close database connection
    fidimo_database.close()

    # Attach fidimo results to attribute table of output map
    # database-connection of current mapset
    mapset_db_settings = dict(x.split(": ") for x in grass.read_command("db.connect",
                                                                        flags="p").splitlines())
    if mapset_db_settings["driver"] != "sqlite":
        #grass.fatal(_("Database driver of current mapset is not 'sqlite'."))
        raise ValueError("Database driver of current mapset is not 'sqlite'.")
    mapset_database = sqlite3.connect(mapset_db_settings["database"])
    mapset_db = mapset_database.cursor()

    mapset_db.execute("ATTACH DATABASE ? AS fidimo_db", (fidimo_db_path,))
    mapset_db.execute('''UPDATE fidimo_output SET 
                reach_length = (SELECT summary_fidimo_result.reach_length FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat),
                source_pop = (SELECT summary_fidimo_result.source_pop FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat),
                fidimo_result = (SELECT summary_fidimo_result.fidimo_result FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat), 
                fidimo_result_lwr = (SELECT summary_fidimo_result.fidimo_result_lwr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat), 
                fidimo_result_upr = (SELECT summary_fidimo_result.fidimo_result_upr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat),
                rel_fidimo_result = (SELECT summary_fidimo_result.rel_fidimo_result FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat), 
                rel_fidimo_result_lwr = (SELECT summary_fidimo_result.rel_fidimo_result_lwr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat), 
                rel_fidimo_result_upr = (SELECT summary_fidimo_result.rel_fidimo_result_upr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v=fidimo_output.cat);''')
    # Commit changes
    mapset_database.commit()
    # Close database connection
    mapset_database.close()


def fidimo_mapping(output):
    grass.run_command("v.colors",
                      flags="g",
                      map=output,
                      use="attr",
                      column="rel_fidimo_result",
                      color="bgyr")

    # grass.run_command("d.vect",
    # map=output)


def print_metadata(fidimo_db_path):
    ''' Print out metadata of FIDIMO database defined
    in fidimo_db_path'''

    fidimo_database = sqlite3.connect(fidimo_db_path)
    fidimo_db = fidimo_database.cursor()

    # Fetch meta data
    fidimo_db.execute('''SELECT * FROM meta''')
    meta_out = "\n".join([(str(x[0]) + ": " + str(x[1]))
                          for x in fidimo_db.fetchall()])

    fidimo_database.close()

    print(meta_out)


# def barrier_correction()


def main():
    ########## Run FIDIMO from input ##########

    input = options['input']
    output = options['output']
    strahler_col = options['strahler_col']
    shreve_col = options['shreve_col']
    network_col = options['network_col']
    barriers = options['barriers']
    fidimo_db_path = options['fidimo_db_path']
    passability_col = options['passability_col']
    threshold = options['threshold']

    source_pop_csv = options['source_pop_csv']
    l = options['l']
    ar = options['ar']
    t = options['t']
    statistical_interval = options['statistical_interval']
    seed_fishmove = options['seed_fishmove']

    ############ DEFINITION CLEANUP TEMPORARY FILES ##############
    # global variables for cleanup
    global tmp_map_vect
    tmp_map_vect = ['streams_tmp', 'barriers_tmp', 'fidimo_net1_tmp',
                    'fidimo_net2_tmp', 'fidimo_net3_tmp', 'output_tmp']

    ############ Start with FIDIMO modules ##############
    # Print out Metadata and exit main function
    if flags['m']:
        print_metadata(fidimo_db_path=options['fidimo_db_path'])
        return None

    if not (flags['s'] or flags['f']):
        # Set up fidimo_db
        create_fidimo_db(fidimo_db_path=fidimo_db_path)

        # Create fidimo network from input data (river shape, barrier points)
        fidimo_network(input=input,
                       output=output,
                       strahler_col=strahler_col,
                       shreve_col=shreve_col,
                       network_col=network_col,
                       barriers=barriers,
                       fidimo_db_path=fidimo_db_path,
                       passability_col=passability_col,
                       threshold=threshold)

        # Copying edges and vertices etc to fidimo_db
        set_fidimo_db(fidimo_db_path=fidimo_db_path)

        # Calculate distance between single river reaches
        fidimo_distance(fidimo_db_path=fidimo_db_path)

    if not (flags['n'] or flags['f']):
        # Append source populations
        fidimo_source_pop(source_pop_csv=source_pop_csv,
                          fidimo_db_path=fidimo_db_path,
                          realisation=flags['r'])

    if not (flags['n'] or flags['s']):
        # Calculate dispersal distances
        my_sigma_dict = sigma_calc( l=int(l),
                                    ar=float(ar),
                                    t=int(t),
                                    statistical_interval=statistical_interval,
                                    seed_fishmove=seed_fishmove)
        # Calcuate fidimo probability
        fidimo_probability(fidimo_db_path=fidimo_db_path,
                           sigma_dict=my_sigma_dict,
                           p=0.67,
                           statistical_interval=statistical_interval)

        # Calculate realisation or multiplication by value of initial source
        # population
        fidimo_realisation( realisation=flags['r'],
                            fidimo_db_path=fidimo_db_path)

        # Get sum of fidimo results for each target reach
        fidimo_summarize(output=output,
                         fidimo_db_path=fidimo_db_path)

        # Map fidimo result to display
        fidimo_mapping(output)

    return 0

if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    sys.exit(main())
