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
#% key: fidimo_dir
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
#% key: upstream_pass_col
#% description: Column name indicating upstream passability rate (0-1) of each barrier
#% required: no
#% guisection: Network parameters
#%end
#%option G_OPT_DB_COLUMN
#% key: downstream_pass_col
#% description: Column name indicating downstream passability rate (0-1) of each barrier
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
#%option
#% key: truncation
#% type: integer
#% key_desc: truncation
#% description: Truncation distance threshold (m). Connections between reaches with distance>truncation will be ignored.
#% required: no
#% answer: -1
#% guisection: Output
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
#% key: b
#% description: Consider barriers (if set) as fully passable.
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
#%flag
#% key: t
#% description: Flag for testing single modules
#%end

#%rules
#% requires_all: barriers,passability_col,threshold
#%end


###########################################
############# Import Libraries ############
###########################################

import os
import sys
import shutil
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

###########################################
############# Define Cleanup ##############
###########################################

############ DEFINITION CLEANUP TEMPORARY FILES ##############
# global variables for cleanup
global tmp_map_vect
tmp_map_vect = []

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
        grass.fatal(
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


def create_fidimo_db( fidimo_dir):
    ''' Create FIDIMO DB'''
     
    # If fidimo_dir exists it will be first removed
    if os.path.exists(fidimo_dir):
        if grass.overwrite()==True:
            shutil.rmtree(fidimo_dir)
            grass.warning(_("FIDIMO dir already exists and will be overwritten"))
        else:
            grass.fatal("FIDIMO directory already exists. Please use overwrite-flag to overwrite the existing FIDIMO directory")
    
    os.makedirs(fidimo_dir)
        
    #grass.message(_("Creating FIDIMO Database and copying edges and vertices"))
    grass.message(_("Creating FIDIMO Database"))
    
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
    fidimo_db = fidimo_database.cursor()
    
    # Create table for meta data
    fidimo_db.execute(
        '''CREATE TABLE meta (parameter VARCHAR(45), value VARCHAR(300))''')
    fidimo_db.executemany("INSERT INTO meta (parameter) VALUES (?)",
                          [(x,)for x in ["Fidimo directory", "Projection", "Network name", "Total reaches n",
                                         "Barriers n", "Source populations n", "Distance matrix created", "Source populations imported",
                                         "Fidimo probabilities calculated", "Last modified", "GRASS setup", "FIDIMO version",
                                         "Scipy version", "Numpy version", "igraph version"]])
    
    # Update metadata
    fidimo_db.execute(
        '''UPDATE meta SET value=? WHERE parameter="Fidimo directory"''', (fidimo_dir,))
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
                   strahler_col,
                   shreve_col,
                   network_col,
                   fidimo_dir,
                   barriers=None,
                   upstream_pass_col=None,
                   downstream_pass_col=None,
                   threshold=25):
    '''GRASS Network based on the preprocessed vector input of streams
    and barriers (if present) is calculated (v.net). 
    Barriers (if present) are snapped to stream network (within threshold) and reaches 
    are split at barriers. Nodes are connected to streams network and 
    lengths of reaches are updated. Nodes and barriers are stored in attribute
    table and can be distinuished by the v_type column (1: barrier, 2: node)'''
    
    grass.message(_("Preparing river network for FIDIMO"))
       
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
        grass.message(_("Connect barriers and nodes to network"))
        # Input barriers
        # dictionary {'target_name':'input_name'}
        passability_col_dict = {"upstream_passability": upstream_pass_col,
                                  "downstream_passability":downstream_pass_col}
        import_vector(input_map=barriers, columns=passability_col_dict,
                      output_map="barriers_tmp" + str(os.getpid()))
        
        # Get network ID for each barrier
        grass.run_command("v.db.addcolumn",
                          quiet=quiet,
                          map="barriers_tmp" + str(os.getpid()),
                          columns="network INTEGER")
        grass.run_command("v.what.vect",
                          quiet=quiet,
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
        
        # Update files to remove at cleanup
        tmp_map_vect.extend(['barriers_tmp', 'fidimo_net1_tmp'])
        
    else:
        grass.message(_("Connect nodes to network"))
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
        grass.fatal("Database driver of current mapset is not 'sqlite'.")
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
                      overwrite=True,
                      flags="t",
                      input="fidimo_net3_tmp" + str(os.getpid()),
                      layer=3,
                      output="output1_tmp" + str(os.getpid()))
    grass.run_command("v.category",
                      quiet=quiet,
                      overwrite=True,
                      input="output1_tmp" + str(os.getpid()),
                      layer="3,1",
                      output="output2_tmp" + str(os.getpid()),
                      option="chlayer")
    
    # Check if output table exists and remove if it exists
    tables_list = grass.read_command("db.tables", flags="p").splitlines()
    if "output2_tmp" + str(os.getpid()) in tables_list:
        grass.run_command("db.droptable",
                          quiet=quiet,
                          flags="f",
                          table="output2_tmp" + str(os.getpid()))
        
    grass.run_command("db.copy",
                      quiet=quiet,
                      overwrite=True,
                      from_table="edges",
                      to_table="output2_tmp" + str(os.getpid()))
    
    grass.run_command("v.db.connect",
                      quiet=quiet,
                      overwrite=True,
                      map="output2_tmp" + str(os.getpid()),
                      layer=1,
                      table="output2_tmp" + str(os.getpid()))
    
    output_columns = grass.read_command("db.columns",
                                        table="output2_tmp" + str(os.getpid())).splitlines()
    grass.run_command("v.db.dropcolumn",
                      quiet=quiet,
                      map="output2_tmp" + str(os.getpid()),
                      layer=1,
                      columns=[x for x in output_columns if x not in ["cat", "orig_cat"]])
    grass.run_command("v.db.addcolumn",
                      quiet=quiet,
                      map="output2_tmp" + str(os.getpid()),
                      layer=1,
                      columns='''length DOUBLE, source_pop DOUBLE, fidimo DOUBLE, fidimo_lwr DOUBLE, 
    fidimo_upr DOUBLE, rel DOUBLE, rel_lwr DOUBLE, rel_upr DOUBLE''')
    
    grass.run_command("v.out.ogr",
                      quiet=quiet,
                      overwrite=grass.overwrite(),
                      input="output2_tmp" + str(os.getpid()),
                      output=os.path.join(fidimo_dir,"fidimo_network"),
                      format="ESRI_Shapefile",
                      output_layer="fidimo_network")
    
    grass.message(_("Final networks prepared for FIDIMO"))
    
    # Update files to remove at cleanup
    tmp_map_vect.extend(['streams_tmp','fidimo_net2_tmp','fidimo_net3_tmp',
                    'output1_tmp','output2_tmp'])
    
    # Update metadata
    grass.verbose(_("Updating Metadata"))
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
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


def set_fidimo_db(fidimo_dir):
    ''' Create and tables for 
    vertices, edges, fidimo distance and barriers in Fidimo DB'''
    
    grass.message(_("Copying edges and vertices to FIDIMO database"))
    
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
    fidimo_db = fidimo_database.cursor()
    
    # Copy edges and vertices to FIDIMO DB ####
    grass.run_command("db.copy",
                      overwrite=True,
                      from_table="vertices",
                      to_database=os.path.join(fidimo_dir,"fidimo_database.db"),
                      to_table="vertices")
    grass.run_command("db.copy",
                      overwrite=True,
                      from_table="edges",
                      to_database=os.path.join(fidimo_dir,"fidimo_database.db"),
                      to_table="edges")
    
    # Copy barriers if exists
    grass.message(_("Copying barriers to FIDIMO database"))
    grass.run_command("db.copy",
                      overwrite=True,
                      from_table="barriers_tmp"+ str(os.getpid()),
                      to_database=os.path.join(fidimo_dir,"fidimo_database.db"),
                      to_table="barriers",
                      where="network IS NOT NULL")
    
    # Create fidimo_distance table
    fidimo_db.execute('''CREATE TABLE fidimo_distance (fidimo_distance_id INTEGER PRIMARY KEY, source INTEGER, target INTEGER, distance DOUBLE, direction INTEGER, 
          from_orig_e INTEGER, to_orig_e INTEGER, from_orig_v INTEGER, to_orig_v INTEGER, source_edge_length DOUBLE, target_edge_length DOUBLE,
          upr_limit DOUBLE, lwr_limit DOUBLE, source_strahler INTEGER, source_shreve INTEGER, 
          target_shreve INTEGER, network INTEGER)''')
    
    # Update metadata
    grass.verbose(_("Updating Metadata"))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    fidimo_database.commit()
    fidimo_database.close()


def add_midpoints_fidimo_db(fidimo_dir):
    '''Add midpoints (vertices) to each river reach to calculate
    distances between (midpoints of) river reaches.
    Each midpoint has its own id that is linked to the cat of the 
    river reach. Function is used in fidimo_distance'''
    
    grass.message(_("Add midpoints (vertices) for each river reach"))
    
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
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


def barrier_passability(  g_barriers,
                          g,
                          v,
                          to):
    '''This calculates the cumulative (i.e. multiplicative) passability
    of barriers from a single node (v) to a list of other nodes (to).
    It requires therefore a dictionary of barriers with passability rates
    and the graph of the network '''
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    
        upstream_barriers = [[x for x in L if x in g_barriers] for L in g.get_shortest_paths(v=v, to=to, output="vpath", mode="OUT")] 
        downstream_barriers = [[x for x in L if x in g_barriers] for L in g.get_shortest_paths(v=v, to=to, output="vpath", mode="IN")] 
        all_barriers = [[x for x in L if x in g_barriers] for L in g.get_shortest_paths(v=v, to=to, output="vpath", mode="ALL")] 
  
    unique_upstream_g_barriers = set([item for sublist in upstream_barriers for item in sublist])
    unique_downstream_g_barriers = set([item for sublist in downstream_barriers for item in sublist])
  
    barrier_pass = [None] * len(all_barriers)
    for i in xrange(len(all_barriers)):
        if not upstream_barriers[i] and not downstream_barriers[i] and not all_barriers[i]: #source=target (i.e no barriers inbetween, passability is unrestricted=1)
      barrier_pass[i] = 1.0 
    elif not upstream_barriers[i] and not downstream_barriers[i] and all_barriers[i]: # direction is first downstream then upstream
        barrier_pass[i] = round(numpy.prod(
          [g_barriers[x][1] if x in unique_downstream_g_barriers else g_barriers[x][0] for x in all_barriers[i]]
          ),8)
    elif upstream_barriers[i] and not downstream_barriers[i] and all_barriers[i]:
        barrier_pass[i] = round(numpy.prod([g_barriers[x][0] for x in upstream_barriers[i]]),8)
    elif not upstream_barriers[i] and downstream_barriers[i] and all_barriers[i]:
        barrier_pass[i] = round(numpy.prod([g_barriers[x][1] for x in downstream_barriers[i]]),8)
    else:
        barrier_pass[i] = -9999
        
    return barrier_pass



def fidimo_distance(fidimo_dir,
                    truncation=-1):
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
    
    fidimo_distance_timer1 = timer()
    
    # connect to database
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
    fidimo_db = fidimo_database.cursor()
    
    # Add midpoints to calculate distance between single river reaches
    add_midpoints_fidimo_db(fidimo_dir)
    
    # Set distance threshold
    if int(truncation) < 0:
        max_dist = 1e8 # Max distance = 100000 km
    else:
        max_dist = int(truncation)  
    
    # Update main distance matrix
    #grass.message(_("Update main distance matrix (fidimo_distance) between all connected river reaches..."))
    grass.message(_(
        "Updating main distance matrix (fidimo_distance) between all connected river reaches..."))
    
    fidimo_db.execute('SELECT DISTINCT network FROM edges')
    network_id = [x[0] for x in fidimo_db.fetchall()]
    
    for i in network_id:
        grass.message(_("   ...processing network ID: %s" % str(i)))
        
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
        midpoints_chunks = [midpoints[x:x + chunk_size_midpoints]
                              for x in xrange(0, len(midpoints), chunk_size_midpoints)]
           
        # Calculate shortest pathss
        grass.message(_("Calculating shortest paths (chunk size: "+str(chunk_size_midpoints)+")..."))
                                      
        for k in range(len(midpoints_chunks)):
            grass.message(_("...chunk: "+str(k+1)+" of "+str(len(midpoints_chunks))))
            
            g_midpoints_chunk_k = [vertices_dict[x] for x in midpoints_chunks[k]]
            
            distance_mat_upstream = g.shortest_paths(source=g_midpoints_chunk_k,
                                                     target=g_midpoints, weights="half_length", mode="IN")
            distance_mat_downstream = g.shortest_paths(source=g_midpoints_chunk_k,
                                                       target=g_midpoints, weights="half_length", mode="OUT")
            distance_mat_all = g.shortest_paths(source=g_midpoints_chunk_k,
                                                target=g_midpoints, weights="half_length", mode="ALL")
            
            # Distinguish between paths up- and downstream [1: downstream, 2:upstream, 3:neither (down-up combination), 4: both directions (e.g. where source=target and dist=0)
            grass.message(_("Calculating directions between reaches"))
            
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
            grass.message(_("Updating distance matrix for specific network and chunk - part 1"))
            paths_to_db = numpy.array([[x for item in midpoints_chunks[k] for x in repeat(item, len(midpoints))],  # list of all source/from midpoints of reaches
                              # list of all target/to midpoints of reaches
                              midpoints * len(midpoints_chunks[k]),
                              # list of distances
                              [item for sublist in distance_mat_all for item in sublist],
                              # list of directions
                              [item for sublist in direction_mat for item in sublist],
                              [i] * (len(midpoints)*len(midpoints_chunks[k]))])  # network
            
            # Only select those network connection that are below the threshold distance
            paths_to_db = zip(*paths_to_db[:,paths_to_db[2,:]<max_dist])
            
            grass.message(_("Updating distance matrix for specific network and chunk - part 2"))
            fidimo_db.executemany(
                "INSERT INTO fidimo_distance (source,target,distance,direction,network) VALUES (?,?,?,?,?)", paths_to_db)
            fidimo_database.commit()
                        
            # Collect garbage to free memory
            del paths_to_db
            del direction_mat
            gc.collect()
        
        # Barriers along each path (if exist)
        fidimo_db.execute(
            '''SELECT COUNT(*) FROM barriers WHERE network=?''', str(i))
        n_barriers_network_i = [x[0] for x in fidimo_db.fetchall()]
        
        if n_barriers_network_i>0 and not flags["b"]: #add b-flag here
            grass.message(_("Processing barriers along networks paths"))
            
            # Get all barriers and passabilities for that specific network
            fidimo_db.execute(
                '''SELECT cat,upstream_passability,downstream_passability FROM barriers WHERE network=?''', str(i))
            barriers_dict = {int(x[0]):[float(x[1]),float(x[0])] for x in fidimo_db.fetchall()}

            g_barriers = {vertices_dict[x]:barriers_dict[x] for x in barriers_dict} ####!!!
            
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
                
    grass.message(_("All networks processed"))
    
    fidimo_distance_timer2 = timer()
    
    # Get number of rows of fidimo_distance and chunk it into pieces of 10E5
    fidimo_db.execute('''SELECT max(rowid) FROM fidimo_distance''')
    max_fidimo_distance_rowid = [x[0] for x in fidimo_db.fetchall()][0]
    fidimo_distance_rowid_chunks = [[x+1,x + 10E5] for x in xrange(0, max_fidimo_distance_rowid, int(10E5))]
    
    grass.message(_("Updating fidimo_distance database..."))
    for k in range(len(fidimo_distance_rowid_chunks)):
        grass.message(_("...chunk: "+str(k+1)+" of "+str(len(fidimo_distance_rowid_chunks))))
        
        # Get cats of original vertices (from-to), Get cats of original edges (from-to)
        #grass.message(_("Updating original vertex categories to fidimo_distance"))
        grass.message(_("Updating original vertex categories to fidimo_distance"))
        fidimo_db.execute('''UPDATE fidimo_distance SET 
                      from_orig_v = (SELECT orig_cat FROM vertices WHERE cat=fidimo_distance.source),
                      to_orig_v = (SELECT orig_cat FROM vertices WHERE cat=fidimo_distance.target)
                      WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
        
        grass.message(_("Updating original river reach (edges) categories to fidimo_distance"))
        fidimo_db.execute('''UPDATE fidimo_distance SET 
                      from_orig_e = (SELECT orig_cat FROM edges WHERE cat=fidimo_distance.from_orig_v),
                      to_orig_e = (SELECT orig_cat FROM edges WHERE cat=fidimo_distance.to_orig_v)
                      WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
      
        # Get edge lengths, stream order and network id for source (and target) reach
        #grass.message(_("Updating addtio nal attributes (e.g. stream order) in fidimo_distance"))
        grass.message(_(
            "Updating addtional attributes (e.g. stream order) in fidimo_distance"))
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
        grass.message(_(
            "Updating addtional attributes (upper and lower distance limit) in fidimo_distance"))
        fidimo_db.execute('''UPDATE fidimo_distance SET
                    lwr_limit = distance-(target_edge_length/2),
                    upr_limit = distance+(target_edge_length/2)
                    WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
        
    fidimo_distance_timer3 = timer()
    
    grass.verbose(_("Time elapsed fidimo_distance: %s and %s" %(str(fidimo_distance_timer2 - fidimo_distance_timer1),
      str(fidimo_distance_timer3 - fidimo_distance_timer1))))
    
    # Update metadata
    grass.verbose(_("Updating Metadata"))
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
    
    grass.verbose(_("Calculating dispersal distances (i.e. sigma dictionary)"))
    
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
                       fidimo_dir,
                       realisation):
    '''This function appends source populations to distance matrix'''
    
    # connect to database
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
    fidimo_db = fidimo_database.cursor()
    
    # Check if FIDIMO DB already contains source populations and check for overwrite flag
    fidimo_db.execute('''SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='fidimo_source_pop';''')
    if sum([x[0] for x in fidimo_db.fetchall()])==1:
        if grass.overwrite()==True:
            # Delete source populations if exists
            fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_source_pop_tmp''')
            fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_source_pop''')
            grass.warning(_("FIDIMO database already contains source populations which will be overwritten"))
        else:
            grass.fatal("FIDIMO database already contains source populations. Please use overwrite-flag to overwrite the existing source populations")
       
    # Read CSV and check if three columns (reach cat, source_col, p)
    with open(source_pop_csv,'rb') as csv_file:
        source_pop_csv_read = csv.DictReader(csv_file,fieldnames=("cat", "source_pop", "p")) # comma is default delimiter,header is assumed 
        # Skip first header line
        next(source_pop_csv_read)
        # csv to list of tuples for input into db
        source_pop_csv_read_to_db = [(int(i['cat']), float(i['source_pop']), float(i['p'])) for i in source_pop_csv_read]
  
    # Check if source pop and p are within a valid range
    n_source_pop = sum([i[1]>0 for i in source_pop_csv_read_to_db])
    if n_source_pop==0:
        grass.fatal(
            "Source population csv contains no source populations > 0")
        
    p_fidimo_source_pop_tmp = [i[2] for i in source_pop_csv_read_to_db]
    if (min(p_fidimo_source_pop_tmp)<0) or (max(p_fidimo_source_pop_tmp)>1):
        grass.fatal(
            "Values for p must be decimal numbers between 0 and 1. NAs not allowed")
        
    source_pop_fidimo_source_pop_tmp = [i[1] for i in source_pop_csv_read_to_db]
    if min(source_pop_fidimo_source_pop_tmp)<0:
        grass.fatal(
            "Values for source populations must be positive integer or decimal numbers. NAs not allowed")
    
    # check if cats of fidimo_source_pop_tmp == orig_cat of edges
    fidimo_db.execute('SELECT DISTINCT orig_cat FROM edges')
    orig_cat_edges = [x[0] for x in fidimo_db.fetchall()]
    cat_fidimo_source_pop_tmp = [i[0] for i in source_pop_csv_read_to_db]
    if sorted(orig_cat_edges) != sorted(cat_fidimo_source_pop_tmp):
        #grass.fatal(_("Vector input map of source populations must must match vector input map that has been used for calculating fidimo distance matrix"))
        grass.fatal(
            "IDs (categories, cat) of source populations must must match cats of vector input map that has been used for calculating fidimo distance matrix")
    else:
        grass.verbose(_(
            "IDs (categories, cat) of source populations match cats of vector input map that has been used for calculating fidimo distance matrix"))
 
   
    fidimo_db.execute("CREATE TABLE fidimo_source_pop_tmp (cat, source_pop, p);")
 
    # Insert source populations into a tmp table          
    fidimo_db.executemany("INSERT INTO fidimo_source_pop_tmp (cat, source_pop, p) VALUES (?, ?, ?);", source_pop_csv_read_to_db)
    fidimo_database.commit()
    
    fidimo_db.execute(
        '''CREATE INDEX fidimo_source_pop_tmp_index ON fidimo_source_pop_tmp (cat)''')
      
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
        'SELECT orig_cat FROM fidimo_source_pop GROUP BY orig_cat HAVING COUNT(*) > 1;')
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
    
    # Add index to source_pop
    fidimo_db.execute('''DROP INDEX IF EXISTS fidimo_source_pop_index''')
    fidimo_db.execute(
        '''CREATE INDEX fidimo_source_pop_index ON fidimo_source_pop (cat)''')
    
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
    grass.message(_("Updating %s source populations in fidimo_distance" %(str(n_source_pop))))
    for k in range(len(fidimo_distance_rowid_chunks)):
        grass.message(_("...chunk: "+str(k+1)+" of "+str(len(fidimo_distance_rowid_chunks)))) 
        fidimo_db.execute('''UPDATE fidimo_distance SET 
                      source_pop = (SELECT source_pop FROM fidimo_source_pop WHERE cat=fidimo_distance.from_orig_v)
                      WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        fidimo_database.commit()
    
    # Update metadata
    grass.verbose(_("Updating Metadata"))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Source populations n"''',(str(n_source_pop),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Source populations imported"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    fidimo_database.commit()
    fidimo_database.close()


def fidimo_probability( fidimo_dir,
                        sigma_dict,
                        p,
                        statistical_interval):
    '''This function calculates dispersal kernel probabilities
    for each river reach and corresponding distances between reaches'''
    
    grass.message(_("Calculation of FIDIMO probablities based on species-specific dispersal parameters"))
    
    # connect to database
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
    fidimo_db = fidimo_database.cursor()
    
    # Get stream orders (ASC) of source populations where source_pop > 0
    fidimo_db.execute(
        '''SELECT DISTINCT source_strahler FROM fidimo_distance WHERE source_pop > 0;''')
    so_list = sorted([x[0] for x in fidimo_db.fetchall()])
    
    grass.debug(_("Test 0 fidimo_probability"),1)
    
    if None in so_list:
        grass.fatal(
            "At least one stream reach with a source population has no value for Strahler stream order")
    elif max([int(x) for x in so_list]) > 9:
        grass.fatal(
            "Stream network has stream orders > 9. Please consider smaller stream network")
    else:
        so_list = [int(x) for x in so_list]
        
    grass.debug(_("Test 1 fidimo_probability"),1)
        
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
    
            
    else:
        grass.debug(_("Test 2 fidimo_probability"),1)
        
        # Before any calculations set fidimo_prob in fidimo_distance to ''
        fidimo_db.execute('''UPDATE fidimo_distance SET basic_fidimo_prob=NULL,basic_fidimo_prob_lwr=NULL,basic_fidimo_prob_upr=NULL;''')
    
    fidimo_database.commit()
    
    # Delete fidimo_prob if exists
    fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_prob''')
    
    grass.debug(_("Test 3 fidimo_probability"),1)
    
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
    
    grass.message(_("Begin with core of fidimo, application of fishmove"))
    
    for i in so_list:
        #i=2
        direction_dict = {"downstream": 1,
                          "upstream": 2, "undef": 3, "source": 4}
        
        for j in direction_dict:  # loop over directions (up- vs. downstream)
            ''' calculate fidimo probability in chunks based on source_strahler
            MAIN PART: leptokurtic probability density kernel based on fishmove '''
            
            # Do not calculate any fidimo probability for direction that are
            # combinations of down- and upstream (as in r.fidimo)
            if j == "undef":
                continue
                
            grass.message(_("Calculation of Fidimo probability for direction: " +
                  j + " and stream order: " + str(i)+"..."))
            
            # Create temp table to collect rows for a specific stream order
            fidimo_db.execute('''DROP TABLE IF EXISTS fidimo_prob_calculation_tmp''')
            
            fidimo_db.execute('''CREATE TABLE fidimo_prob_calculation_tmp AS SELECT
              fidimo_distance_id AS fidimo_distance_id,
              upr_limit AS upr_limit,
              lwr_limit AS lwr_limit,
              source_shreve AS source_shreve,
              target_shreve AS target_shreve
              FROM fidimo_distance
              WHERE source_strahler = ?
              AND direction = ?
              AND source_pop > 0''', (i, direction_dict[j]))
                      
            # Split all cases into chunks of max. chunk_size
            fidimo_db.execute('''SELECT max(rowid) FROM fidimo_prob_calculation_tmp''')
            max_fidimo_prob_calculation_tmp_rowid = [x[0] for x in fidimo_db.fetchall()][0]
            fidimo_prob_calculation_rowid_chunks = [[x+1,x + 500] for x in xrange(0, max_fidimo_prob_calculation_tmp_rowid, int(500))]
            
            for k in range(len(fidimo_prob_calculation_rowid_chunks)):
                #k=1
                grass.message(_("...chunk: " + str(k + 1) +
                      " of " + str(len(fidimo_prob_calculation_rowid_chunks))))
                               
                #chunk = fidimo_prob_calculation_rowid_chunks[k]
                
                fidimo_db.execute('''SELECT * FROM fidimo_prob_calculation_tmp 
                  WHERE rowid BETWEEN %s and %s;'''%(fidimo_prob_calculation_rowid_chunks[k][0],fidimo_prob_calculation_rowid_chunks[k][1]))
                fidimo_prob_calculation_colnames = dict(
                    zip([x[0] for x in fidimo_db.description], range(0, len(fidimo_db.description))))
                fidimo_prob_calculation_array = scipy.array(fidimo_db.fetchall())
                #fidimo_db.execute('SELECT * FROM fidimo_distance WHERE fidimo_distance_id IN (%s)' %
                 #                 ','.join('?' * len(chunk)), tuple(chunk))
                #fidimo_distance_colnames = dict(
                #    zip([x[0] for x in fidimo_db.description], range(0, len(fidimo_db.description))))
                #fidimo_distance_array = scipy.array(fidimo_db.fetchall())
                
                # Create results array to collect fidimo_prob results
                fidimo_result_array = fidimo_prob_calculation_array[:, fidimo_prob_calculation_colnames["fidimo_distance_id"]].astype(int)
                
                for l in nrun:
                    #l="fit"
                    # CDF(upr) - CDF(lwr)
                    if j == "source":
                        basic_fidimo_prob = (fidimo_kernel_cdf(x=fidimo_prob_calculation_array[:, fidimo_prob_calculation_colnames["upr_limit"]].astype(float),
                                                              sigma_stat=sigma_dict["stat"][l][i],
                                                              sigma_mob=sigma_dict["mob"][l][i],
                                                              p=p) - 
                                            fidimo_kernel_cdf( x=0,
                                                      sigma_stat=sigma_dict["stat"][l][i],
                                                      sigma_mob=sigma_dict["mob"][l][i],
                                                      p=p)) * 2.0
                    elif j == "upstream" or j == "downstream":
                        basic_fidimo_prob = fidimo_kernel_cdf(x=fidimo_prob_calculation_array[:, fidimo_prob_calculation_colnames["upr_limit"]].astype(float),
                                                      sigma_stat=sigma_dict["stat"][l][i],
                                                      sigma_mob=sigma_dict["mob"][l][i],
                                                      p=p) - fidimo_kernel_cdf( x=fidimo_prob_calculation_array[:, fidimo_prob_calculation_colnames["lwr_limit"]].astype(float),
                                                      sigma_stat=sigma_dict["stat"][l][i],
                                                      sigma_mob=sigma_dict["mob"][l][i],
                                                      p=p)
                                                      
                    # if upstream direction than correct for network splits by
                    # a shreve stream order approach
                    if j == "upstream":
                        basic_fidimo_prob = basic_fidimo_prob * (fidimo_prob_calculation_array[:, fidimo_prob_calculation_colnames["target_shreve"]].astype(
                            float) / fidimo_prob_calculation_array[:, fidimo_prob_calculation_colnames["source_shreve"]].astype(float))
                        
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
    grass.message(_(
        "Updating fidimo_distance with basic (unweighted) fidimo probabilities"))
    fidimo_db.execute('''UPDATE fidimo_distance SET 
                basic_fidimo_prob = (SELECT basic_fidimo_prob FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                basic_fidimo_prob_lwr = (SELECT basic_fidimo_prob_lwr FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                basic_fidimo_prob_upr = (SELECT basic_fidimo_prob_upr FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id)
                 WHERE EXISTS (SELECT fidimo_distance_id FROM fidimo_prob WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id)''')
    
    grass.message(_("Updating Metadata"))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Fidimo probabilities calculated"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    # Commit changes
    fidimo_database.commit()
    
    # Close database connection
    fidimo_database.close()


def fidimo_realisation( realisation,
                        statistical_interval,
                        fidimo_dir):
    '''This function either weights the fidimo probabilities with the value of the value of
    the source population or calculates realised fish count (real integer values) from the
    probabilities using the multinomial distribution'''
    
    # connect to database
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
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
        
    else:
        # Before any calculations set fidimo_result in fidimo_distance to ''
        grass.debug(_("Test 1 fidimo_realisation"),1)
        fidimo_db.execute(
        '''UPDATE fidimo_distance SET fidimo_result=NULL,fidimo_result_lwr=NULL,fidimo_result_upr=NULL;''')
    
    fidimo_database.commit()
    grass.debug(_("Test 2 fidimo_realisation"),1)
    
    if realisation == True:
        grass.message(_(
            "Calculate realisation (i.e. fish counts that disperse from each source population)"))
        
        # Get number of runs for statistical intervals
        if statistical_interval == "no":
            nrun = zip(["basic_fidimo_prob", "basic_fidimo_prob_lwr", "basic_fidimo_prob_upr"],[True,False,False])
        else:
            nrun = zip(["basic_fidimo_prob", "basic_fidimo_prob_lwr", "basic_fidimo_prob_upr"],[True,True,True])
        
        # Realisation using the mutlinomial distribution to obtain indiviual counts
        # Get ids of all source reaches that have source_pop>0
        fidimo_db.execute(
            'SELECT DISTINCT cat,source_pop FROM fidimo_source_pop WHERE source_pop>0')
        source_populations = [[x[0], x[1]] for x in fidimo_db.fetchall()]
        
        # Insert results of realisation in temporary table and then update
        # fidimo_distance from that table
        fidimo_db.execute('''CREATE TEMP TABLE realisation_result_tmp 
              (fidimo_distance_id INTEGER, fidimo_result DOUBLE, fidimo_result_lwr DOUBLE, fidimo_result_upr DOUBLE)''')
        
        # Create index on fidimo_distance to improve query based on from_orig_v
        fidimo_db.execute('''CREATE INDEX fidimo_distance_index_from_orig_v ON fidimo_distance (from_orig_v)''')
        
        for i in source_populations:
            fidimo_db.execute(
                '''SELECT fidimo_distance_id, 
                basic_fidimo_prob, 
                basic_fidimo_prob_lwr, 
                basic_fidimo_prob_upr
                FROM fidimo_distance WHERE from_orig_v = ?''', (i[0],))
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
                grass.error(_(
                    "Realisation cannot be calculated due to some statistical intervals"))
        
        # Drop index from fidimo distance
        fidimo_db.execute('''DROP INDEX IF EXISTS fidimo_distance_index_from_orig_v''')
        
        # Create index on realisation_result_tmp
        fidimo_db.execute('''CREATE INDEX realisation_result_tmp_index ON realisation_result_tmp (fidimo_distance_id)''')

        # Get number of rows of fidimo_distance and chunk it into pieces of 10E5
        fidimo_db.execute('''SELECT max(rowid) FROM fidimo_distance''')
        max_fidimo_distance_rowid = [x[0] for x in fidimo_db.fetchall()][0]
        fidimo_distance_rowid_chunks = [[x+1,x + 10E5] for x in xrange(0, max_fidimo_distance_rowid, int(10E5))]
        
        # Join fidimo_prob with fidimo_distance
        grass.message(_("Updating fidimo_distance with realised fidimo results (i.e. fish counts)..."))
        for k in range(len(fidimo_distance_rowid_chunks)):
            grass.message(_("...chunk: "+str(k+1)+" of "+str(len(fidimo_distance_rowid_chunks))))        
            fidimo_db.execute('''UPDATE fidimo_distance SET 
                    fidimo_result = (SELECT fidimo_result FROM realisation_result_tmp WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                    fidimo_result_lwr = (SELECT fidimo_result_lwr FROM realisation_result_tmp WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id),
                    fidimo_result_upr = (SELECT fidimo_result_upr FROM realisation_result_tmp WHERE fidimo_distance_id=fidimo_distance.fidimo_distance_id)
                    WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        
    else:
        # Get number of rows of fidimo_distance and chunk it into pieces of 10E5
        fidimo_db.execute('''SELECT max(rowid) FROM fidimo_distance''')
        max_fidimo_distance_rowid = [x[0] for x in fidimo_db.fetchall()][0]
        fidimo_distance_rowid_chunks = [[x+1,x + 10E5] for x in xrange(0, max_fidimo_distance_rowid, int(10E5))]
      
        # Join fidimo_prob with fidimo_distance
        grass.message(_("Updating fidimo_distance with source-population-weighted fidimo probability..."))
        for k in range(len(fidimo_distance_rowid_chunks)):
            grass.message(_("...chunk: "+str(k+1)+" of "+str(len(fidimo_distance_rowid_chunks))))  # Multiply by value of inital source population
            fidimo_db.execute('''UPDATE fidimo_distance SET
                      fidimo_result = basic_fidimo_prob*source_pop,
                      fidimo_result_lwr = basic_fidimo_prob_lwr*source_pop,
                      fidimo_result_upr = basic_fidimo_prob_upr*source_pop
                      WHERE rowid BETWEEN %s and %s;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))
        
    # Update metadata
    grass.verbose(_("Updating Metadata"))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))
    
    # Commit changes
    fidimo_database.commit()
    
    # Close database connection
    fidimo_database.close()


def fidimo_summarize( output,
                      fidimo_dir):
    '''This function summarizes the fidimo result for each target reach
    and writes results back to output vector map'''
    
    grass.message(_("FIDIMO summarize"))
    
    time_fidimo_summarize1 = timer()

    # connect to database
    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
    fidimo_db = fidimo_database.cursor()

    # First, delete summary_fidimo_prob if exists
    fidimo_db.execute('''DROP TABLE IF EXISTS summary_fidimo_result''')
    
    # Create index on fidimo_distance (to_orig_v)
    grass.debug(_("Test 1 fidimo_summarize"),1)
    fidimo_db.execute('''DROP INDEX IF EXISTS fidimo_distance_index_to_orig_v''')
    fidimo_db.execute('''CREATE INDEX fidimo_distance_index_to_orig_v ON fidimo_distance (to_orig_v)''')

    # Split summarize into chunks
    fidimo_db.execute('''SELECT max(rowid) FROM fidimo_distance''')
    max_fidimo_distance_rowid = [x[0] for x in fidimo_db.fetchall()][0]
    fidimo_distance_rowid_chunks = [[x+1,x + 10E5] for x in xrange(0, max_fidimo_distance_rowid, int(10E5))]
 
    # Create temp table to collect subresults
    fidimo_db.execute('''DROP TABLE IF EXISTS summary_fidimo_result_tmp''')
    fidimo_db.execute('''DROP INDEX IF EXISTS summary_fidimo_result_tmp_index_to_orig_v''')
    fidimo_db.execute('''CREATE TEMP TABLE summary_fidimo_result_tmp 
              (to_orig_v INTEGER, reach_length DOUBLE, fidimo_result DOUBLE, fidimo_result_lwr DOUBLE, fidimo_result_upr DOUBLE)''')
 
    grass.message(_("Summing up fidimo result per target reach..."))
    for k in range(len(fidimo_distance_rowid_chunks)):
        grass.message(_("...chunk: "+str(k+1)+" of "+str(len(fidimo_distance_rowid_chunks))))

        # Summarize fidimo prob for each target reach
        fidimo_db.execute('''INSERT INTO summary_fidimo_result_tmp
          (to_orig_v,reach_length,fidimo_result,fidimo_result_lwr,fidimo_result_upr)
          SELECT
          to_orig_v, 
          target_edge_length AS reach_length,
          CAST(SUM(fidimo_result) AS DOUBLE) AS fidimo_result, 
          CAST(SUM(fidimo_result_lwr) AS DOUBLE) AS fidimo_result_lwr,
          CAST(SUM(fidimo_result_upr) AS DOUBLE) AS fidimo_result_upr 
          FROM fidimo_distance 
          WHERE direction!=3 AND
          rowid BETWEEN %s and %s
          GROUP BY to_orig_v;'''%(fidimo_distance_rowid_chunks[k][0],fidimo_distance_rowid_chunks[k][1]))

        # Commit changes
        fidimo_database.commit()
     
    # Summing up subresults
    fidimo_db.execute('''CREATE INDEX summary_fidimo_result_tmp_index_to_orig_v ON summary_fidimo_result_tmp (to_orig_v)''')
    fidimo_db.execute('''CREATE TABLE summary_fidimo_result AS SELECT
    to_orig_v AS to_orig_v, 
    reach_length AS reach_length,
    CAST(SUM(fidimo_result) AS DOUBLE) AS fidimo_result, 
    CAST(SUM(fidimo_result_lwr) AS DOUBLE) AS fidimo_result_lwr,
    CAST(SUM(fidimo_result_upr) AS DOUBLE) AS fidimo_result_upr 
    FROM summary_fidimo_result_tmp GROUP BY to_orig_v;''')

    # Drop index on fidimo_distance (to_orig_v)
    fidimo_db.execute('''DROP INDEX IF EXISTS summary_fidimo_result_tmp_index_to_orig_v''')

    grass.debug(_("Test 2 fidimo_summarize"),1)

    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN source_pop DOUBLE''')
    fidimo_db.execute('''UPDATE summary_fidimo_result SET
                  source_pop = (SELECT source_pop FROM fidimo_source_pop WHERE cat=summary_fidimo_result.to_orig_v)''')
    # Commit changes
    fidimo_database.commit()
    
    grass.debug(_("Test 3 fidimo_summarize"),1)

    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN rel_fidimo_result DOUBLE''')
    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN rel_fidimo_result_lwr DOUBLE''')
    fidimo_db.execute(
        '''ALTER TABLE summary_fidimo_result ADD COLUMN rel_fidimo_result_upr DOUBLE''')
    
    grass.debug(_("Test 4 fidimo_summarize"),1)

    fidimo_db.execute('''UPDATE summary_fidimo_result SET
                  rel_fidimo_result = fidimo_result/reach_length,
                  rel_fidimo_result_lwr = fidimo_result_lwr/reach_length,
                  rel_fidimo_result_upr = fidimo_result_upr/reach_length''')
    # Update metadata
    grass.verbose(_("Updating Metadata"))
    fidimo_db.execute('''UPDATE meta SET value=? WHERE parameter="Last modified"''',
                      (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),))

    # Commit changes
    fidimo_database.commit()
    # Close database connection
    fidimo_database.close()
    
    time_fidimo_summarize2 = timer()
    
    # Import network map from fidimo_dir
    grass.run_command("v.in.ogr",
                      quiet=quiet,
                      overwrite=grass.overwrite(),
                      input=os.path.join(fidimo_dir,"fidimo_network/fidimo_network.shp"),
                      output=output,
                      key="cat")

    # Attach fidimo results to attribute table of output map
    # database-connection of current mapset
    mapset_db_settings = dict(x.split(": ") for x in grass.read_command("db.connect",
                                                                        flags="p").splitlines())
    if mapset_db_settings["driver"] != "sqlite":
        grass.fatal("Database driver of current mapset is not 'sqlite'.")
    mapset_database = sqlite3.connect(mapset_db_settings["database"])
    mapset_db = mapset_database.cursor()
    
    grass.message(_("Updating attribute table of output map"))
    mapset_db.execute("ATTACH DATABASE ? AS fidimo_db", (os.path.join(fidimo_dir,"fidimo_database.db"),))
    mapset_db.execute('''UPDATE {output} SET 
                length = (SELECT summary_fidimo_result.reach_length FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat),
                source_pop = (SELECT summary_fidimo_result.source_pop FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat),
                fidimo = (SELECT summary_fidimo_result.fidimo_result FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat), 
                fidimo_lwr = (SELECT summary_fidimo_result.fidimo_result_lwr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat), 
                fidimo_upr = (SELECT summary_fidimo_result.fidimo_result_upr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat),
                rel = (SELECT summary_fidimo_result.rel_fidimo_result FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat), 
                rel_lwr = (SELECT summary_fidimo_result.rel_fidimo_result_lwr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat), 
                rel_upr = (SELECT summary_fidimo_result.rel_fidimo_result_upr FROM fidimo_db.summary_fidimo_result summary_fidimo_result 
                  WHERE summary_fidimo_result.to_orig_v={output}.cat);'''.format(output=output))
    # Commit changes
    mapset_database.commit()
    # Close database connection
    mapset_database.close()
    
    time_fidimo_summarize3 = timer()
   
    #grass.message(_("Time elapsed: %s" %str(end-start)))
    grass.verbose(_("Time for fidimo_summarize : %s and %s" % (str(time_fidimo_summarize2 - time_fidimo_summarize1),
      str(time_fidimo_summarize3 - time_fidimo_summarize1))))



def fidimo_mapping(output):
    grass.run_command("v.colors",
                      flags="g",
                      map=output,
                      use="attr",
                      column="rel_fidimo_result",
                      color="bgyr")

    # grass.run_command("d.vect",
    # map=output)


def print_metadata(fidimo_dir):
    ''' Print out metadata of FIDIMO database defined
    in fidimo_dir'''

    fidimo_database = sqlite3.connect(os.path.join(fidimo_dir,"fidimo_database.db"))
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
    fidimo_dir = options['fidimo_dir']
    upstream_pass_col = options['upstream_pass_col']
    downstream_pass_col = options['downstream_pass_col']
    threshold = options['threshold']
    truncation = options['truncation']

    source_pop_csv = options['source_pop_csv']
    l = options['l']
    ar = options['ar']
    t = options['t']
    statistical_interval = options['statistical_interval']
    seed_fishmove = options['seed_fishmove']
    
    ############### Testing single modules ##################
    ######### exit main after module is finished ############
    if flags['t']:
        # Testing single modules and exit main after module is finished
        # Set up fidimo_db
        create_fidimo_db(fidimo_dir=fidimo_dir)

        # Create fidimo network from input data (river shape, barrier points)
        fidimo_network(input=input,
                       strahler_col=strahler_col,
                       shreve_col=shreve_col,
                       network_col=network_col,
                       barriers=barriers,
                       fidimo_dir=fidimo_dir,
                       passability_col=passability_col,
                       threshold=threshold)
        return None



    ############ Start with FIDIMO modules ##############
    # Print out Metadata and exit main function
    if flags['m']:
        print_metadata(fidimo_dir=options['fidimo_dir'])
        return None

    if not (flags['s'] or flags['f']):
        # Set up fidimo_db
        create_fidimo_db(fidimo_dir=fidimo_dir)

        # Create fidimo network from input data (river shape, barrier points)
        fidimo_network(input=input,
                       strahler_col=strahler_col,
                       shreve_col=shreve_col,
                       network_col=network_col,
                       barriers=barriers,
                       fidimo_dir=fidimo_dir,
                       passability_col=passability_col,
                       threshold=threshold)

        # Copying edges and vertices etc to fidimo_db
        set_fidimo_db(fidimo_dir=fidimo_dir)

        # Calculate distance between single river reaches
        fidimo_distance(fidimo_dir=fidimo_dir,
                        truncation=truncation)

    if not (flags['n'] or flags['f']):
        # Append source populations
        fidimo_source_pop(source_pop_csv=source_pop_csv,
                          fidimo_dir=fidimo_dir,
                          realisation=flags['r'])

    if not (flags['n'] or flags['s']):
        # Calculate dispersal distances
        my_sigma_dict = sigma_calc( l=int(l),
                                    ar=float(ar),
                                    t=int(t),
                                    statistical_interval=statistical_interval,
                                    seed_fishmove=seed_fishmove)
        # Calcuate fidimo probability
        fidimo_probability(fidimo_dir=fidimo_dir,
                           sigma_dict=my_sigma_dict,
                           p=0.67,
                           statistical_interval=statistical_interval)

        # Calculate realisation or multiplication by value of initial source
        # population
        fidimo_realisation( realisation=flags['r'],
                            statistical_interval=statistical_interval,
                            fidimo_dir=fidimo_dir)

        # Get sum of fidimo results for each target reach
        fidimo_summarize(output=output,
                         fidimo_dir=fidimo_dir)

        # Map fidimo result to display
        #fidimo_mapping(output) # does not work yet

    return 0

if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    sys.exit(main())
