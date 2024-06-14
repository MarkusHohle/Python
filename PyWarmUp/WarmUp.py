# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 01:57:32 2024

@author: MMH_user
"""
#### a few things interviewers would ask you on the 1st technical interview!##

#########some Python warm-ups before interview#################################
###############################################################################

class functions:
    
    def f_default(a, b = 2):
        
        print(a+b)
        
    def f_optional(a, b = 2, *c):
        
        if c:#tuple
            for c in c:
                print(a**c)
            
        else:
            print('no c available')
            
    def f_keyword(a, b = 2, *c, **d):
        
        if 'hi' in d:#dict
            print(d['hi'])

###############################################################################
###############################################################################
class CountInstances:
    
    ct = 0
    
    def __init__(self):
        CountInstances.ct += 1
        
    @classmethod
    def reset(cls):
        cls.ct = 0

###############################################################################
###############################################################################
class Super():
    
    class C1():
    
        def f1(self,a,b):
            res = a + b
            
            #print(res)
            
            return res
            
    class C2(C1):
        
        def f2(self,a,b,c,d):
            res = super().f1(a,b)
            print(c+d+res)

###############################################################################
###############################################################################
class Encoder():
    
    def __init__(self):
        
        NT   = ['A', 'C', 'G', 'T']
        Code = [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]] 
        
        Dict = {key: value for key, value in zip(NT,Code)}
        
        self.Enc = lambda Sequence: [Dict[s] for s in Sequence]
        
        return(self.Enc)
        
    def Encode(self, Sequence):
        
        print(self.Enc(Sequence))

###############################################################################
###############################################################################
class Map(Encoder): #combining inheritance, super & maping
    
    def __init__(self):
        self.Enc = super().__init__()
    
    def fun(self, s):
        NT   = ['A', 'C', 'G', 'T']
        Code = [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]] 
        Dict = {key: value for key, value in zip(NT,Code)}
        
        return Dict[s]

    def runboth(self,S):
        
        print(list(map(self.fun, S)))#mapping actual function
        print(list(map(self.Enc, S)))#mapping anonymus function
###############################################################################
###############################################################################
class Recursion():
    
    def MinusOne(self,n):
        
        n -= 1
        
        if n > 0:
            print('n = ' + str(n))
            self.MinusOne(n)
            
        else:
            
            print('n = ' + str(n) + ' done!')
            
###############################################################################
###############################################################################
import time
import multiprocessing as mp
from datetime import datetime
from multiprocessing import cpu_count

class Multi():
    
    def __init__(self):
        self.n = cpu_count()
    
    def DoSomething(self, sleeptime = 3):
        
        time.sleep(sleeptime)# sleep for sleeptime seconds
        print('Good Morning!')
        
    def Serial(self):
        t1 = datetime.now()
        
        for i in range(self.n): 
            self.DoSomething()
        
        t2 = datetime.now()
        
        print('total time sequential = ' + str(t2-t1))

    def Parallel(self):
        t1 = datetime.now()
        
        Processes = [mp.Process(target = self.DoSomething, args = ()) for i in range(self.n)]
        for p in Processes:
            p.start()
        for p in Processes:
            p.join()
        
        t2 = datetime.now()
        
        print('total time parallel = ' + str(t2-t1))

            
###############################################################################
###############################################################################
class DecoratorAsWrapper():
    
    #this function just measures time
    
    def my_timer(my_function):
        def get_args(*args,**kwargs):
            t1 = time.monotonic()
            results = my_function(*args,**kwargs)
            t2 = time.monotonic()
            dt = t2 - t1
            print("Total runtime: " + str(dt))
            return results
        return get_args
    
    def __init__(self):
        self.n = cpu_count()
    
    def DoSomething(self, sleeptime = 3):
        
        time.sleep(sleeptime)# sleep for sleeptime seconds
        print('Good Morning!')
    
    @my_timer
    def Serial(self):
        #t1 = datetime.now()
        
        for i in range(self.n): 
            self.DoSomething()
        
        #t2 = datetime.now()
        
        #print('total time sequential = ' + str(t2-t1))

    @my_timer
    def Parallel(self):
        #t1 = datetime.now()
        
        Processes = [mp.Process(target = self.DoSomething, args = ()) for i in range(self.n)]
        for p in Processes:
            p.start()
        for p in Processes:
            p.join()
        
        #t2 = datetime.now()
        
        #print('total time parallel = ' + str(t2-t1))

###############################################################################
###############################################################################