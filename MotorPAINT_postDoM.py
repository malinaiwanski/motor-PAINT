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
        self.mintracklength = 6 #5 # default @ 4 because you loose lots of BGD without afecting the tracks you see too much minimum length of the tracks you want to include minimum is 3 because you need at least two segments to see if your track is straight
        self.maxtracklength = 2000 #maximum length of tracks you want to include
        self.min_displacement = 200 #nm
        self.minspeed = 100 #nm/s included in the vectormap
        self.maxspeed = 1500 #nm/s included in the vectormap #was at 15000
        self.exptime = 0.060 # sec was at 0.1
        
        self.max_angle = 50 # maximum angle that can be observed in a track default 75
        self.stepsize= 3 #in frames for displacement vectors used to calculate angle along tracks
        self.upperbound = 0.9
        self.lowerbound =0.1 #Confinement ratio Net displacement over total dispalcement


    def report_parameters(self, print_to_screen=True, write_to_log=True):
        if print_to_screen:
            print('The current parameter values are listed. You can redefine each parameter')
            print(f"mintracklength : {self.mintracklength} frames \nmaxtracklength : {self.maxtracklength} frames")
            print(f"min_displacement : {self.min_displacement} nm \nminspeed : {self.minspeed} nm/s ")
            print(f"maxspeed : {self.maxspeed} nm/s \nexptime : {self.exptime} sec ")
            print(f"max_angle : {self.max_angle} degrees \nstepsize displacement vectors:{self.stepsize} Note these are only applied when you call the function to filter for angle speretly"  )
            print(f"Confinement ratio min-max:{self.lowerbound,self.upperbound}")
        if write_to_log:
            #write this to log file
            self.log('Parameter values used for filtering the Results table are listed.')
            self.log(f"mintracklength : {self.mintracklength} frames \nmaxtracklength : {self.maxtracklength} frames \nmin_displacement : {self.min_displacement} nm \nminspeed : {self.minspeed} nm/s \nmaxspeed : {self.maxspeed} nm/s \nexptime : {self.exptime} sec \nConfinement ratio min-max:{self.lowerbound,self.upperbound}")
            #self.log(f"min_displacement : {self.min_displacement} nm \nminspeed : {self.minspeed} nm/s ")
            #self.log(f"maxspeed : {self.maxspeed} nm/s \nexptime : {self.exptime} sec ")
            #self.log(f"max_angle : {self.max_angle} degrees"  )


    def filter_data(self):
        """
        Calculate track displacements and other track parameters and filters based on given parameter values



        """
        self.report_parameters(print_to_screen=False)

        #Read in the data set
        dfs=pd.read_csv(self.datadir, error_bad_lines=False, warn_bad_lines=True)
        dfs.sort_values(['Track_ID','Particle_ID'], inplace=True)

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
        dfs['Displacement_(nm)'] = np.sqrt(dfs['dX_track_(nm)']**2 + dfs['dY_track_(nm)']**2)

        #Calculate the (eucledean) distance between each consecuative point within a track

        dfs['dx']=dfs['X_(nm)']-dfs.groupby(['Track_ID'])['X_(nm)'].shift()
        dfs['dy']=dfs['Y_(nm)']-dfs.groupby(['Track_ID'])['Y_(nm)'].shift()
        dfs['Step_Distance_(nm)'] = np.sqrt(dfs['dx']**2 + dfs['dy']**2)

        #Calculate speed and netto displacement per track for filtering
        dfs['Speed_(nm/s)'] =(dfs.groupby(['Track_ID'])['Step_Distance_(nm)'].sum()/(dfs.groupby(['Track_ID']).size()-1)/self.exptime)
        dfs=dfs[ (dfs['Speed_(nm/s)'] > self.minspeed) & (self.maxspeed > dfs['Speed_(nm/s)']) &\
       (dfs['Displacement_(nm)'] > self.min_displacement)]

        #number_of_tracks=len(dfs['Track_ID'].unique())
        #self.log(f'The number of tracks after filtering for traveled distance per track and speed: {number_of_tracks}')


        dfs['Net/total']= (dfs['Displacement_(nm)']/dfs.groupby(['Track_ID'])['Step_Distance_(nm)'].sum())
        self.dfs=dfs[(dfs['Net/total'] > self.lowerbound)&(dfs['Net/total'] < self.upperbound)]

        #self.log(f'The number of tracks after filtering for Confinement: {number_of_tracks}')
        
    def filter_for_angle(self, minduration = 3, actually_filter=False): 
        print('Filtering for angle....')
        
        self.log(f'\nmax_angle : {self.max_angle} degrees \nstepsize displacement vectors:{self.stepsize}')
        
        self.dfs['dx2']=self.dfs['X_(nm)']-self.dfs.groupby(['Track_ID'])['X_(nm)'].shift(self.stepsize)
    
        self.dfs['dy2']=self.dfs['Y_(nm)']-self.dfs.groupby(['Track_ID'])['Y_(nm)'].shift(self.stepsize)
    
        vectors=self.dfs[['dx2', 'dy2']].to_numpy()
        magnitude=np.sqrt((vectors[:,0]**2+vectors[:,1]**2))
        #normalize vectors
        vn=vectors/magnitude[:,None]
        #make a shifted copy

        vn_shifted=np.delete(vn,0, 0)
        newrow=[np.nan,np.nan]
        vn_shifted=np.vstack([vn_shifted,newrow])


        #Calcualte the angle using the dot product
        angle = np.arccos(np.sum(vn*vn_shifted, axis=1))
        #Convert to degrees
        degree=np.rad2deg(angle)
        #Add to the data frame
        self.dfs['Angle_degree_shifted']=degree.tolist()
    
        #Shift angle values up (without changing the index) such that the angle value between the first and second vector is in the first row
        #etc the last rows of track ID will have Nan value
        self.dfs['Angle_degree'] = self.dfs['Angle_degree_shifted'].shift(-self.stepsize)


        self.dfs['Angle_degree'].fillna(value = self.dfs['Angle_degree_shifted'].shift(-1), inplace=True)

        self.dfs.drop(['Angle_degree_shifted','dy2', 'dx2' ], axis = 1, inplace = True)

        g = self.dfs.groupby('Track_ID').cumcount(ascending=False)

        #Pad the edge with the angles as if you are walking backward on the track and fill the other two Nans with the mean as you cannot caluculate 
        #an angle with only one vector, but if you do not fill in an angel the point will be removed later in the code
        #This code as it is works only for a stepsize of 3


        if self.stepsize-3 !=0: 
            print('Code works currently only for a stepsize of 3!')
    
        pd.options.mode.chained_assignment = None  # default='warn'
        print('Warning turned of : pd.options.mode.chained_assignment = None')

        self.dfs['Angle_degree'][g ==  1] = self.dfs['Angle_degree'][g ==  2]
        self.dfs['Angle_degree'][g ==  2] = self.dfs['Angle_degree'][g ==  3]
        self.dfs['Angle_degree'][g ==  3] = self.dfs['Angle_degree'][g ==  0] =  (self.dfs['Angle_degree'][g ==  5]+self.dfs['Angle_degree'][g ==  4])/2
    

        #Creates a column which states whether the angle is below max_angle (True) or above (False)
        self.dfs['under_max_angle'] = self.dfs['Angle_degree']<self.max_angle
        #Creates another column which gives consequative rows with the same boolean the same number
        self.dfs['streak_id'] = (self.dfs['under_max_angle'] != self.dfs.groupby('Track_ID')['under_max_angle'].shift(1)).cumsum()
        
        
        if actually_filter:
            #here tracks are actually removed based on length and angle
            
            #Only keep streteched where consequative boolean value is the same for longer than 2 frames
            self.dfs=self.dfs.groupby(['streak_id']).filter(lambda x: len(x)>minduration )
    
            #Remove all the long stretches where the angle is above max angle
            self.dfs=self.dfs[self.dfs.under_max_angle]
            
            print(f'Filtered for angle, all stretches of tracks below {self.max_angle} were removed and strecthec shorted than {minduration} frames as well, the index is reset such that each stretch becomes a track on its own')
            self.dfs = self.dfs.reset_index().set_index('streak_id')
            self.dfs.rename(columns = {'Track_ID':'Old_Track_ID'}, inplace = True)
            self.dfs.index.names = ['Track_ID']
            
        #Caluclate the x and y displacement over the whole stretch which is a part of a whole track

        self.dfs['dX_stretch_(nm)'] = self.dfs.groupby('Track_ID')['X_(nm)'].nth(-1)-self.dfs.groupby('Track_ID')['X_(nm)'].nth(0) #dx= x_end-x_start
        self.dfs['dY_stretch_(nm)'] = self.dfs.groupby('Track_ID')['Y_(nm)'].nth(-1)-self.dfs.groupby('Track_ID')['Y_(nm)'].nth(0) #dy
        self.dfs['Displacement_stretch_(nm)'] = np.sqrt(self.dfs['dX_stretch_(nm)']**2 + self.dfs['dY_stretch_(nm)']**2)

        #Calculate the (eucledean) distance between each consecuative point within a track

        self.dfs['dx_s']=self.dfs['X_(nm)']-self.dfs.groupby(['Track_ID'])['X_(nm)'].shift()
        self.dfs['dy_s']=self.dfs['Y_(nm)']-self.dfs.groupby(['Track_ID'])['Y_(nm)'].shift()
        self.dfs['Step_Distance_stretch(nm)'] = np.sqrt(self.dfs['dx_s']**2 + self.dfs['dy_s']**2)

        #Calculate speed and netto displacement per track for filtering
        self.dfs['Speed_stretch(nm/s)'] =(self.dfs.groupby(['Track_ID'])['Step_Distance_stretch(nm)'].sum()/(self.dfs.groupby(['Track_ID']).size()-1)/self.exptime)
        #self.dfs=self.dfs[ (self.dfs['Speed_(nm/s)'] > self.minspeed) & (self.maxspeed > self.dfs['Speed_stretch(nm/s)']) &\
       #(self.dfs['Displacement_(nm)'] > self.min_displacement)]

        number_of_tracks=len(self.dfs.index.unique())
        print(f'The number of stretches after filtering for angle and splitting tracks: {number_of_tracks}')
        self.log(f'The number of stretches after filtering for angle and splitting tracks: {number_of_tracks}')


        self.dfs['Net/total_stretch']= (self.dfs['Displacement_stretch_(nm)']/self.dfs.groupby(['Track_ID'])['Step_Distance_stretch(nm)'].sum())
        #self.self.dfs=self.dfs[(self.dfs['Net/total_stretch'] > self.lowerbound)&(self.dfs['Net/total_stretch'] < self.upperbound)]
        
        
        return True

    def divide_two_directions(self, soma_coords=None):

        """ Given the x and y coordinate of the middle of the soma an extra column is added
        which states wether the track is towards (plus in = True ) the soma coordinates or away from (False) """
        self.soma_coords=soma_coords
        self.dfs['Xn']=self.dfs['X_(px)']-self.soma_coords[0][0]
        self.dfs['Yn']=self.dfs['Y_(px)']-self.soma_coords[0][1]
        #Displacement is towards the soma (plus end in) if the end of the track is closer to the soma coordinate than
        #the beginning i.e. the  difference with the soma coordinate is smaller at the end of the track
        self.dfs['X_towards_soma'] = abs(self.dfs.groupby('Track_ID')['Xn'].nth(-1))<abs(self.dfs.groupby('Track_ID')['Xn'].nth(0)) #x_end<x_start
        self.dfs['Y_towards_soma'] = abs(self.dfs.groupby('Track_ID')['Yn'].nth(-1))<abs(self.dfs.groupby('Track_ID')['Yn'].nth(0))

        #If the displacement along the y positioin/displacement along x direction is smaller than 1: use the change in x coordinates to determine
        #wether the track is plus in, otherwise the displacement in the y direction is used

        # use track displacement as conditon :
        #dfs['dY_track_(nm)'] = dfs.groupby('Track_ID')['Y_(nm)'].nth(-1)-dfs.groupby('Track_ID')['Y_(nm)'].nth(0)

        self.dfs['Plus_IN']=pd.DataFrame.where(cond=(abs(self.dfs['dY_track_(nm)']/self.dfs['dX_track_(nm)'])<1),\
                                                 self=self.dfs['X_towards_soma'], other=self.dfs['Y_towards_soma'])
        #clean up dataframe
        self.dfs=self.dfs.drop(['Xn', 'Yn','X_towards_soma','Y_towards_soma'], axis = 1)

    def split_tracks(self, two_directions=True, save_tables = True):
            #remove track_ID as index
        dfs=self.dfs.reset_index()

    #Devide into frout tables, based on the overall direction of the track
    #The Orgin (0,0) is in the upper left corner of the image, so that Left Bottom (LB) equals dx<0 &dy >0
        Results_filtered = dfs
        Results_RT = dfs[ (dfs['dX_track_(nm)'] > 0) & (dfs['dY_track_(nm)'] < 0)]
        Results_RB = dfs[ (dfs['dX_track_(nm)'] > 0) & (dfs['dY_track_(nm)'] > 0)]
        Results_LT = dfs[ (dfs['dX_track_(nm)'] < 0) & (dfs['dY_track_(nm)'] < 0)]
        Results_LB = dfs[ (dfs['dX_track_(nm)'] < 0) & (dfs['dY_track_(nm)'] > 0)]

        save_name=['Results_filtered' , 'Results_RT','Results_RB','Results_LT','Results_LB','Results_PLUSIN','Results_PLUSOUT']
        table_list=[Results_RT,Results_RB,Results_LT,Results_LB]

        if two_directions:
            Results_PLUSIN = dfs[ (dfs['Plus_IN'] ==True) ]
            number_of_localisations_IN=len(Results_PLUSIN.index.unique())
            print('number of localisations toward soma', round(number_of_localisations_IN,2))
            self.log(f'number of localisations toward soma {round(number_of_localisations_IN,2)}')    
            Results_PLUSOUT = dfs[ (dfs['Plus_IN'] ==False) ]
            number_of_localisations_OUT=len(Results_PLUSOUT.index.unique())
            print('number of localisations away from soma', round(number_of_localisations_OUT,2))
            self.log(f'number of localisations away from soma {round(number_of_localisations_OUT,2)}')
            
            plusin_overtotal = round(number_of_localisations_IN/(number_of_localisations_OUT+number_of_localisations_IN), 2)

            print(f'Ratio in over total {plusin_overtotal}')
            self.log(writetolog=f'Ratio in over total {plusin_overtotal}')

            table_list=[Results_filtered, Results_RT,Results_RB,Results_LT,Results_LB, Results_PLUSIN,Results_PLUSOUT]

        else:
            save_name=save_name[0:5]


        if save_tables: 

            for i,table in enumerate(table_list):
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
