#!/usr/bin/env python
# FILENAME: foo.py

from ctypes import *

foo = CDLL('./foo.so')

myprint = foo.myprint
myprint.argtypes = [POINTER(c_char)] 
myprint.restype = c_char_p 
res = myprint('hello')
print res

add = foo.add
add.argtypes = [c_float, c_float]
add.restype = c_float
print add(1.3, 4.2)
"""
hifi = CDLL('./libhifi_fu.so')
a = hifi.hifi_m
a("a")
"""
