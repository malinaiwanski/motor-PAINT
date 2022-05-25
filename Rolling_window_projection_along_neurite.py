# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:45:59 2022

@author: Hanna
"""

import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter



class projection_3:
    
    def __init__(self):
        self.df              = None
        self.change          = 'X' #Needed for median calculation if the change in X direction is bigger than Y direction the X coordnates get sorted and median for Y is calcluated
        self.window_width    = 3000 #length of line to project onto (in nm)
        self.n               = 0 #counts the number of windows
        self.overlap_window  = 0 #set to 1 for half a window overlap
        self.counter_offset  = 0
        self.pixelsize       = 65 #nm
        self.stored_window   = {} #note it is restarted for every neurite!
        
        #Coordinates of a point in the soma of the cell
        self.soma_coords = None
        #This class should be initiated per neurite
        #Number of neurite 
        self.i = None
        
        
        
        
        
    def define_neurite_axis(self):
        
        if self.change =='X':
            plt.plot(self.df['X_(px)'],self.df['Y_(px)'], 'o', ms=0.1 )
            self.df['Y_medianfilt']=median_filter(self.df['Y_(px)'], 1800, mode='nearest')
            plt.plot(self.df['X_(px)'],self.df['Y_medianfilt'],  'ko', ms=0.8)
            plt.show()  

            if abs(self.df['X_(px)'].iloc[-1]-self.soma_coords[0][0]) < abs(self.df['X_(px)'].iloc[0]-self.soma_coords[0][0]): 
            # If last point in the dataframe lies closer to the the soma, reverse the datat frame such that distance calculation along the median starts from the soma
                self.df=self.df[::-1]

            ###### Calculate the distance along the median line ###########
            self.df['Step_distance_(px)'] = np.sqrt((self.df['X_(px)']-self.df['X_(px)'].shift())**2 + (self.df['Y_medianfilt']-self.df['Y_medianfilt'].shift())**2)
        
        elif self.change == 'Y': 
            plt.plot(self.df['X_(px)'],self.df['Y_(px)'], 'o', ms=0.1 )
            self.df['X_medianfilt']=median_filter(self.df['X_(px)'], 1800, mode='nearest')
            plt.plot(self.df['X_medianfilt'], self.df['Y_(px)'], 'ko', ms=0.8)
            plt.show()   

            if abs(self.df['Y_(px)'].iloc[-1]-self.soma_coords[0][1]) < abs(self.df['Y_(px)'].iloc[0]-self.soma_coords[0][1]): 
            # If last point in the dataframe lies closer to the the soma, reverse the datat frame such that distance calculation along the median starts from the soma
                self.df=self.df[::-1]

            ###### Calculate the distance along the median line ###########
            self.df['Step_distance_(px)'] = np.sqrt((self.df['Y_(px)']-self.df['Y_(px)'].shift())**2 + (self.df['X_medianfilt']-self.df['X_medianfilt'].shift())**2)
        
        
        #To make sure the distance value along the line is not exactly the same for consequative rows 
        #(would cause problem later in the code)

        self.df['Step_distance_(px)']=self.df['Step_distance_(px)']+0.0000001
        #remove NaN for the first row and replace by zero 
        self.df['Step_distance_(px)']   = self.df['Step_distance_(px)'].fillna(0)
        self.df['Distance_(px)']        = self.df['Step_distance_(px)'].cumsum()
        self.df['Distance_median_(nm)'] = self.df['Distance_(px)']*self.pixelsize
        
        
    def rollBy(self,basis,window,overlap_window,*args,**kwargs):
        #note that basis must be sorted in order for this to work properly
        self.window         = window
        self.overlap_window = overlap_window
        windows_min         = basis.min()
    
        windows_max         = basis.max()
        window_starts       = np.arange(windows_min, windows_max, window)
        
        #For overlapping windows
        if self.overlap_window != 0: 
            if self.overlap_window == 1: 
                window_starts1 = np.arange(windows_min+(window/(self.overlap_window+1)), windows_max, window)
                
            if self.overlap_window == 2: 
                window_starts1 = np.arange(windows_min+(window/(self.overlap_window+1)), windows_max, window)
                window_starts2 = np.arange(windows_min+((self.overlap_window*window)/(self.overlap_window+1)), windows_max, window)
                window_starts2 = pd.Series(window_starts2, index = window_starts2)
        
            window_starts1 = pd.Series(window_starts1, index = window_starts1)
            
        window_starts = pd.Series(window_starts, index = window_starts)
        
        
        #indexed_what = pd.Series(what.values,index=basis.values)
        indexed_what = pd.DataFrame(self.df.values,index=basis.values, columns=self.df.columns)
        
        
        #make columns for coordinates of the line to project onto
        for f in range(0, self.overlap_window+1):
            
                self.df[f'{f}P_x']   = np.nan
                self.df[f'{f}P_y']   = np.nan
                self.df[f'{f}Q_x']   = np.nan
                self.df[f'{f}Q_y']   = np.nan
        

        def applyToWindow(val):
            # using slice_indexer rather that what.loc [val:val+window] allows
            # window limits that are not specifically in the index

            indexer = indexed_what.index.slice_indexer(val,val+window,1)

            chunk = indexed_what.iloc[indexer]

            return self.func(data=chunk, *args,**kwargs)
        
        self.overlap_window_counter = 0 
        window_starts.apply(applyToWindow)
                                           
                
        if self.overlap_window == 1:
            self.overlap_window_counter = 1 #Counter for round of windowing
            self.n=1 #1, 3, 5, ... counter for individual plots
            window_starts1.apply(applyToWindow)
            
        if self.overlap_window == 2: 
            self.overlap_window_counter = 1
            self.n = 1
            window_starts1.apply(applyToWindow)
            self.overlap_window_counter = 2
            self.n = 2
            window_starts2.apply(applyToWindow)
                
             
        
    def func(self, data):
        
               
        # A and B form a line along the median

        if self.change == 'X':
            A = (data.loc[data.index[0],'X_(px)'],data.loc[data.index[0],'Y_medianfilt'])
            B = (data.loc[data.index[-1],'X_(px)'],data.loc[data.index[-1],'Y_medianfilt'])
        elif self.change == 'Y':
            A = (data.loc[data.index[0],'X_medianfilt'], data.loc[data.index[0],'Y_(px)'])
            B = (data.loc[data.index[-1],'X_medianfilt'], data.loc[data.index[-1],'Y_(px)'])
        
        

        delta = np.array(A) - np.array(B) # should be windowsize? 

        distance = np.linalg.norm(delta)

        width = self.window_width/self.pixelsize #windowheight = 2000/pixelsize
        rect_x = 0.5*width*delta[0]/distance

        rect_y = 0.5*width*delta[1]/distance
        
        #calculate corner points of the window
        r1 = (A[0]-rect_y, A[1]+rect_x)
        r2 = (A[0] + rect_y, A[1]- rect_x)
        r3 = (B[0]-rect_y, B[1]+rect_x)
        r4 = (B[0] + rect_y, B[1]- rect_x)
        

        #Line to project is the midline boardered by P and Q
        #B is always further along than A, so r1 and r2 form the start boundry and r3 and r4 the end boundry of the window
        P = ((r1[0]+r3[0])/2, (r1[1]+r3[1])/2)
        Q = ((r2[0]+r4[0])/2, (r2[1]+r4[1])/2)
        

        #select all the points within the window
        slope_r12, intercept_r12 = np.polyfit([r1[0], r2[0]],[r1[1], r2[1]],1)
        slope_r13, intercept_r13 = np.polyfit([r1[0], r3[0]],[r1[1], r3[1]],1)
        slope_r34, intercept_r34 = np.polyfit([r3[0], r4[0]],[r3[1], r4[1]],1)
        slope_r24, intercept_r24 = np.polyfit([r2[0], r4[0]],[r2[1], r4[1]],1)
         
            
        if A[1]<B[1]:
            condition  = (self.df['Y_(px)']>(slope_r12*self.df['X_(px)']+intercept_r12)) & \
            (self.df['Y_(px)']<(slope_r34*self.df['X_(px)']+intercept_r34)) 
        else:
            condition  = (self.df['Y_(px)']<(slope_r12*self.df['X_(px)']+intercept_r12)) & \
            (self.df['Y_(px)']>(slope_r34*self.df['X_(px)']+intercept_r34))
            
            
        if A[0]>B[0]:
            condition = condition & (self.df['Y_(px)']>(slope_r24*self.df['X_(px)']+intercept_r24)) &\
            (self.df['Y_(px)']<(slope_r13*self.df['X_(px)']+intercept_r13))
        else: 
            condition = condition & (self.df['Y_(px)']<(slope_r24*self.df['X_(px)']+intercept_r24)) &\
            (self.df['Y_(px)']>(slope_r13*self.df['X_(px)']+intercept_r13))
        
        #fill the data frame columns with the line coordinates 
        #on which these points should be projected
         
        f=self.overlap_window_counter
        
        #Note as points can occur in sequential windows this way of storing info overwrites previous
        self.df.loc[condition, f'{f}P_x'] = P[0]
        self.df.loc[condition, f'{f}P_y'] = P[1]
        self.df.loc[condition, f'{f}Q_x'] = Q[0]
        self.df.loc[condition, f'{f}Q_y'] = Q[1]
        
    
        #self.df.loc[condition, 'delta']=delta
        
      
        self.df[f'projection{f}'] = ((self.df['X_(px)']-self.df[f'{f}P_x'])*(self.df[f'{f}Q_x']-self.df[f'{f}P_x'])+\
        (self.df['Y_(px)']-self.df[f'{f}P_y'])*(self.df[f'{f}Q_y']-self.df[f'{f}P_y']))/width
        
        self.df[f'projection_(nm){f}']=self.pixelsize*self.df[f'projection{f}']
        
        
        A_ = data.loc[data.index[0],'Distance_median_(nm)']
        B_ = data.loc[data.index[-1],'Distance_median_(nm)']
        
        if self.n>12:
            
            if B_-A_ > self.window-(self.window/100):
         
                self.stored_window[f'Neurite{self.i}_window_{self.n}'] =\
                    [self.df[['Distance_median_(nm)', f'{f}P_x', f'projection_(nm){f}', 'Plus_IN']][(self.df[f'{f}P_x']==P[0])], [A_, B_]] 


                self.plot_window(condition, r1, r2, r3, r4, A, B, P, Q)
                print(self.n)
                self.n+=1 + self.overlap_window 
            
        else: 
            if self.n<10:
                self.stored_window[f'Neurite{self.i}_window_0{self.n}'] =\
                [self.df[[ f'{f}P_x', f'projection_(nm){f}', 'Plus_IN']][(self.df[f'{f}P_x']==P[0])], [A_, B_]] 


            else: 
                self.stored_window[f'Neurite{self.i}_window_{self.n}'] =\
                [self.df[['Distance_median_(nm)', f'{f}P_x', f'projection_(nm){f}', 'Plus_IN']][(self.df[f'{f}P_x']==P[0])], [A_, B_]] 


            self.plot_window(condition, r1, r2, r3, r4, A, B, P, Q)
        
            print(self.n)
            self.n+=1 + self.overlap_window 
            

        
    def plot_window(self, condition,  r1, r2, r3, r4, A, B, P, Q):
        
        
        fig1, ax = plt.subplots()    
        
        #Plot bounding window
        points_x=[r1[0],r2[0],r4[0], r3[0], r1[0], A[0], B[0]]
        points_y=[r1[1],r2[1],r4[1],r3[1], r1[1], A[1], B[1]]
        ax.plot(points_x,points_y, color='black')
        #ax.scatter(points_x,points_y, c=['#87CEFA','#778899','#B0C4DE','#FFFFE0','#87CEFA', '#9932CC', '#8B0000'])
        ax.plot([P[0], Q[0]],[P[1], Q[1]], color='firebrick' , label='projectionline' )
        #ax.plot([A[0],B[0]],[A[1],B[1]], label='Line along median points')

        #plot points meeting the condition, should be in the window, with a cross
        ax.plot(self.df['X_(px)'][condition], self.df['Y_(px)'][condition], 'o', color='slateblue', markersize=5)

        #plot all the other points as well 

        ax.plot(self.df['X_(px)'],self.df['Y_(px)'], 'o',color='slategrey', ms=0.1 )

        #Plot the midline
        if self.change == 'X':
                ax.plot(self.df['X_(px)'],self.df['Y_medianfilt'],  'ko', ms=0.8)

        elif self.change == 'Y':
                ax.plot(self.df['X_medianfilt'], self.df['Y_(px)'],  'ko', ms=0.8)


       # annotations=['r1', 'r2', 'r4','r3','r1' ,'A', 'B']
        #for n, label in enumerate(annotations):
            #  plt.annotate(label, (points_x[n], points_y[n]))
        ax.legend()
        
         
        if self.n<10: 
            fig1.savefig(os.getcwd() + fr"\\Neurite{self.i}_window_0{self.n}.jpeg")
        else:
            fig1.savefig(os.getcwd() + fr"\\Neurite{self.i}_window_{self.n}.jpeg")
        
        plt.show() 
        
        
        plt.close(fig1) 
            
              
            

        return True
    
#Reasoning behind projection
# Dot product of vectors u and v =|u||v|cos(q), projection of u onto v is |u|cos(q)=(u.v)/|v|=(ux*vx+uy*vy)/|v|
#if v is a line spanned by points A=(ax,ay) and B=(bx,by) and u is a vector from P to a random point (x,y) then 
#v=(Q-P)=[(qx-px),(qy-py)], |v|=sqrt((qx-px)**2+(qy-py)**2) and u =(point-P)=[(x-px),(y-py)], this gives
#(ux*vx+uy*vy)/|v|= ((x-px)*(qx-px)+(y-py)*(qy-py))/|v| 
# |v| is Length AB          
     