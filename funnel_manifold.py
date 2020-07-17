# -*- coding: utf-8 -*-

"""
Created on Fri Jul 17 15:03:41 2020
"""

from utilities import sign
import math

class funnel():
 
    def __init__(self):
        self.suggestedFrame = {
            "yMin": -1.5,
            "yMax": 1.5,
            "xMin": -1,
            "xMax": +1,
            }

    def x_of_y(self,y):
        x= math.exp(y-1)
        return(x)

    def keep_in_funnel(self,r):
        assert(len(r)==2)
        if abs(r[0])>self.x_of_y(r[1]):
            r[0] =sign(r[0])*self.x_of_y(r[1])

    #def showAll():
        
        