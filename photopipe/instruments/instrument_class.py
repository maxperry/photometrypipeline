import numpy as np
import sys
from abc import ABCMeta, abstractmethod

class instrument(object):

    __metaclass__ = ABCMeta
    
    def __init__(self, name, camnum):
        self.name = name
        self.camname  = ['C'+str(i) for i in np.arange(camnum)]
        
        self.flatname = 'flat'
        self.biasname = 'bias'
        self.darkname = 'dark'
        self.objname  = 'img'
        
        self.ftype_post = {self.objname: 'o', self.flatname: 'f', 
            self.biasname: 'b', self.darkname: 'd'}
                
    ########
    @abstractmethod
    def has_cam_bias(self, idx):
        # Input: index of camera
        # Output: boolean saying if indexed camera needs bias calibration
        pass

    @abstractmethod
    def has_cam_dark(self, idx):
        # Input: index of camera
        # Output: boolean saying if indexed camera needs dark calibration
        pass

    @abstractmethod
    def change_header_keywords(self, h, cam):
        pass

    @abstractmethod
    def slice(self, cam):
        # Input: camera name
        # Output: slicing of specified camera
        pass

    @abstractmethod
    def is_cam_split(self, idx):
        # Input: index of camera
        # Output: boolean saying if indexed camera is a split filter camera
        pass

    @abstractmethod
    def get_cam_sat(self, h, idx):
        # Input: header, index of camera
        # Output: saturation of indexed camera
        pass
        
    @abstractmethod
    def get_cam_gain(self, h, idx):
        # Input: header, index of camera
        # Output: gain of indexed camera
        pass

    @abstractmethod
    def get_exptime(self, h):
        # Input: header
        # Output: exposure time of frame
        pass

    @abstractmethod
    def get_filter(self, h, cam):
        # Input: header, camera name (ex. 'C0', 'C3a' for split filter camera)
        # Output: filter name for file
        pass
        
    @abstractmethod
    def possible_filters(self):
        # Input: nothing
        # Output: instruments' list of possible filters
        pass
        
    @abstractmethod
    def get_centered_filter(self,h,idx):
        # Input: header, index of camera
        # Output: filter name at center of indexed camera
        pass

    @abstractmethod
    def original_file_format(self):
        # Input: nothing
        # Output: file name format (ex. '????????T??????C??.fits')
        pass

    @abstractmethod
    def change_file_names(self, files):
        # Input: files to change name
        # Output: change file names in same directory to match 
        #         '????????T??????C{}{}.fits'.format(cam_i, ftype_post) format
        pass

    ########