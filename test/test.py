##Class test : Vector class

import numpy as np

class Vector(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __add__(self,othervector):
        added_x = self.x + othervector.x
        added_y = self.y + othervector.y
        return Vector(added_x,added_y)
    def __sub__(self,othervector):
        sub_x = self.x - othervector.x
        sub_y = self.y - othervector.y
        return Vector(sub_x, sub_y)
    def mag(self):
        return (self.x**2+self.y**2)**0.5
    def deltax(self, othervector):
        return othervector.x - self.x
    def deltay(self,othervector):
        return othervector.y - self.y
    #def dotprod(self,othervector):
    #    return mag(self)*mag(othervector)*cos

bob = Vector(1,2)
alice = Vector(4,5)
charles = bob+alice

print charles.deltax(bob)
