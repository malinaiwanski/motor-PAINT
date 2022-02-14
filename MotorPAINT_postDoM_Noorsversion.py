# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 16:16:22 2022
Last updated; 25 Jan 2022

@author: Noor
"""

import numpy as np
import pandas as pd
import math
from scipy import interpolate as interp

import sys as sys
import os
import tkinter as tk
from tkinter import filedialog

import matplotlib.pylab as plt

from datetime import datetime

#To write a log file
import logging as logging



class Motor_paint: 
    """Perform motor-PAINT analysis using a 
    results table with track IDs in format produced by DoM plugin"""
    
    def __init__(self):
        self.logfileopen='no'
        self._select_file()
        
        self.log(getstarted=True)

        
        self._set_default_parameters()
        self.report_parameters(write_to_log=False)
        
    def log(self, writetolog='', getstarted=False, done= False):
        """Creating and Configuring Logger. Function allows statements to be written to a logfile. 
        
        :param writelog: string to log to the logfile.log
        :type writelog: str
        :param getstarted: Initiates the creation of a log-file, call this part of the function only once at the beginning
        :type getstarted: bool
        :param done: closes the handlers so that the log file can be closed, call this after analysis has finished
        :type done:bool
            
            
        """
       
        if getstarted & done: 
            print('Error- Either getstarted OR done should be True, not both!')
        if getstarted:
            if self.logfileopen == 'yes':
                print('logfile is already open')
            else: 
                for handler in logging.root.handlers[:]:
                    logging.root.removeHandler(handler)
                logging.basicConfig(filename = "logfile.log", 
                            filemode = "w", 
                            level = logging.INFO)
                self.logger = logging.getLogger()
                self.logger.info(f"{datetime.now()}")
                self.logger.info(f'This is a logfile for motor-PAINT analysis of : {self.basename}, Data directory: {self.datadir}')
                self.logfileopen='yes'
        if done: 
            #This empties the file?
            print('This option does not work yet')
            #logging.shutdown()
        
        if getstarted == False & done ==False: 
            self.logger.info(writetolog)

      

        
    def _select_file(self):
        """
        Allows user to select a file in pop-up window, saves the directory. 

        """
        #create the root window
        root = tk.Tk()
        root.withdraw()
        
        filetypes = (("CSV Files","*.csv"),)
        
        #Detected particles that were DomLinked
        raw = filedialog.askopenfilename(\
                                         title = 'Choose Linked Particle Table OR Drift Corrected Linked Particle Table format:.csv',\
                                           filetypes=filetypes  ) 
                        
        
        self.datadir=raw
        self.basename = str(os.path.basename(raw))
        self.savename = self.basename[:-4]
        print (self.savename)
        print(raw)
        
        #Set the current working directory
        os.chdir(os.path.dirname(raw))
        
        root.destroy()
        
    def _set_default_parameters(self):
        """ ..."""
        self.mintracklength = 4 #5 # default @ 4 because you loose lots of BGD without afecting the tracks you see too much minimum length of the tracks you want to include minimum is 3 because you need at least two segments to see if your track is straight
        self.maxtracklength = 2000 #maximum length of tracks you want to include
        self.min_displacement = 150 #nm 
        self.minspeed = 100 #nm/s included in the vectormap
        self.maxspeed = 1500 #nm/s included in the vectormap #was at 15000
        self.exptime = 0.060 # sec was at 0.1
        self.max_angle = 'None' #75# maximum angle that can be observed in a track default 75

    def report_parameters(self, print_to_screen=True, write_to_log=True):
        if print_to_screen:
            print('The current parameter values are listed. You can redefine each parameter')   
            print(f"mintracklength : {self.mintracklength} frames \nmaxtracklength : {self.maxtracklength} frames")
            print(f"min_displacement : {self.min_displacement} nm \nminspeed : {self.minspeed} nm/s ") 
            print(f"maxspeed : {self.maxspeed} nm/s \nexptime : {self.exptime} sec ")
            print(f"max_angle : {self.max_angle} degrees"  )
        if write_to_log: 
            #write this to log file
            self.log('Parameter values used for filtering the Results table are listed.')   
            self.log(f"mintracklength : {self.mintracklength} frames \nmaxtracklength : {self.maxtracklength} frames \nmin_displacement : {self.min_displacement} nm \nminspeed : {self.minspeed} nm/s \nmaxspeed : {self.maxspeed} nm/s \nexptime : {self.exptime} sec \nmax_angle : {self.max_angle} degrees")
            #self.log(f"min_displacement : {self.min_displacement} nm \nminspeed : {self.minspeed} nm/s ") 
            #self.log(f"maxspeed : {self.maxspeed} nm/s \nexptime : {self.exptime} sec ")
            #self.log(f"max_angle : {self.max_angle} degrees"  )
        
    def filter_data(self):
        """
        Calculate track displacements and other track parameters and filters based on given parameter values
        
            
        
        """
        self.report_parameters(print_to_screen=False)
        
        #Read in the data set
        dfs=pd.read_csv(self.datadir)

        #Put original columnnames in a list for later when we want to restore the table to the orginial format
        self.og_cols=list(dfs.columns)
        
        number_of_tracks=len(dfs['Track_ID'].unique())
        self.log(f'The number of tracks before filtering: {number_of_tracks}')
        
        #filter out the longest and shortest tracks
        dfs=dfs[ (self.maxtracklength > dfs.Track_Length) & (dfs.Track_Length >= self.mintracklength)]

        #Make calculations per track by setting track ID as index
        dfs=dfs.set_index(['Track_ID'])
        
        number_of_tracks=len(dfs.index.unique())
        self.log(f'The number of tracks after filtering for track length: {number_of_tracks}')

        #Caluclate the x and y displacement over the whole track

        dfs['dX_track_(nm)'] = dfs.groupby('Track_ID')['X_(nm)'].nth(-1)-dfs.groupby('Track_ID')['X_(nm)'].nth(0) #dx= x_end-x_start
        dfs['dY_track_(nm)'] = dfs.groupby('Track_ID')['Y_(nm)'].nth(-1)-dfs.groupby('Track_ID')['Y_(nm)'].nth(0) #dy
        dfs['Nettrack_Distance_(nm)'] = np.sqrt(dfs['dX_track_(nm)']**2 + dfs['dY_track_(nm)']**2)

        #Calculate the (eucledean) distance between each consecuative point within a track

        dfs['dx']=dfs['X_(nm)']-dfs.groupby(['Track_ID'])['X_(nm)'].shift()
        dfs['dy']=dfs['Y_(nm)']-dfs.groupby(['Track_ID'])['Y_(nm)'].shift()
        dfs['Step_Distance_(nm)'] = np.sqrt(dfs['dx']**2 + dfs['dy']**2)

        #Calculate speed and netto displacement per track for filtering
        dfs['Speed_(nm/s)'] =(dfs.groupby(['Track_ID'])['Step_Distance_(nm)'].sum()/(dfs.groupby(['Track_ID']).size()-1)/self.exptime)
        self.dfs=dfs[ (dfs['Speed_(nm/s)'] > self.minspeed) & (self.maxspeed > dfs['Speed_(nm/s)']) &\
       (dfs['Nettrack_Distance_(nm)'] > self.min_displacement)]
            
        number_of_tracks=len(self.dfs.index.unique())
        self.log(f'The number of tracks after filtering for traveled distance per track and speed: {number_of_tracks}')
        
    def split_tracks(self):
        #remove track_ID as index
        dfs=self.dfs.reset_index()

    #Devide into frout tables, based on the overall direction of the track
    #The Orgin (0,0) is in the upper left corner of the image, so that Left Bottom (LB) equals dx<0 &dy >0
 
        Results_RT = dfs[ (dfs['dX_track_(nm)'] > 0) & (dfs['dY_track_(nm)'] < 0)]
        Results_RB = dfs[ (dfs['dX_track_(nm)'] > 0) & (dfs['dY_track_(nm)'] > 0)]
        Results_LT = dfs[ (dfs['dX_track_(nm)'] < 0) & (dfs['dY_track_(nm)'] < 0)]
        Results_LB = dfs[ (dfs['dX_track_(nm)'] < 0) & (dfs['dY_track_(nm)'] > 0)]


        save_name=['Results_RT','Results_RB','Results_LT','Results_LB']

        for i,table in enumerate([Results_RT,Results_RB,Results_LT,Results_LB]):
            #Delete the columns
            for table_col in list(table.columns): 
                if table_col not in self.og_cols: 
                    del table[f'{table_col}']
            
            #Set order of cols as in DoM Results table
            table=table[self.og_cols]
    
            #Save as CSV in current working directory
            print(fr'{save_name[i]}.csv') 
            self.log(fr'{save_name[i]}.csv')
            table.to_csv(path_or_buf=fr'{save_name[i]}.csv', index=False)
            
        self.log(f' Saved split tabled to {self.datadir}')