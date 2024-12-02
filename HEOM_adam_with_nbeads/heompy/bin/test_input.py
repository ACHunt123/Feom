#!/usr/bin/env python
from heom.input import InputObj

file_input = '/home/ap837/code/heompy/input/input'

inp = InputObj(file_input)

for key,value in inp.__dict__.items():
    print(key,value)
