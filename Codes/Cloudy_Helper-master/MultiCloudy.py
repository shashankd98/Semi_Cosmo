#!/usr/bin/env python3

import sys, os, math, threading, subprocess, queue

########## Initialization ##########

# Directory of CLOUDY
CLOUDY = '~/CLOUDY/source/cloudy.exe'
# How many theads you want to use? Don't set this number larger than your core number.
Max_Threads = 4
# How many points you want to calculate?
Model_Num = 6
# What is the name of your model?
Model_Name = 'test'

def InputFile(sequence, model, inputscript):
    # Write CLOUDY input script here
    # Here gives an example show a calculation with 6 points, let log10(temperature)+log10(hden) = 10
    # How to calculate temperature, sequence begin from 0
    temperature = 4+sequence
    # How to calcuate hydrogen density
    hden = 10 - temperature
    # Add and only add one \n at end of each line
    inputscript.write('coronal equilibrium {}\n'.format(temperature))
    inputscript.write('hden {}\n'.format(hden))
    inputscript.write('stop zone 1\n')
    inputscript.write('set dr 0\n')
    # Do NOT change format of name of the saving file, change extension is OK
    inputscript.write('save cooling "grid{0:09d}_{1}.col" last no hash\n'.format(sequence, model))
    inputscript.write('save overview "grid{0:09d}_{1}.ovr" last no hash\n'.format(sequence, model))
    return 0

# Extensions of files you save.
Extensions = ['col','ovr']
    
########## End of initialization ##########

class CalThread(threading.Thread):
    def __init__(self, sequence, model):
        threading.Thread.__init__(self)
        self.sequence = sequence
        self.model = model

    def run(self):
        inputscript = open('grid{0:09d}_{1}.in'.format(self.sequence, self.model), 'w')
        InputFile(self.sequence, self.model, inputscript)
        inputscript.close()
        subprocess.call('{0} -r grid{1:09d}_{2}'.format(CLOUDY, self.sequence, self.model), shell=True)

Error_Queue = queue.Queue()

def GatherFile(suffix, model):
    Label = 0 # 0 is false, 1 is true
    Output = open('./{0}.{1}'.format(model, suffix), 'w')
    blankline = 0
    for i in range(Model_Num):
        if(not os.path.exists('./grid{0:09d}_{1}.{2}'.format(i, model, suffix))):
            Error_Queue.put(i)
            continue
            
        SingleFile = open('./grid{0:09d}_{1}.{2}'.format(i, model, suffix), 'r')
        Lines = SingleFile.readlines()
        SingleFile.close()
        if(len(Lines) < 2):
            Error_Queue.put(i)
            blankline += 1
            if(len(Lines) == 1 and Label == 0):
                Output.write(Lines[0])
                Label = 1
        else:
            if(Label == 0):
                Output.write(Lines[0])
                Label = 1
            while( blankline !=0 ):
                Output.write('\n')
                blankline -= 1
            for j in range(1, len(Lines)):
                Output.write(Lines[j])
        os.remove('grid{0:09d}_{1}.{2}'.format(i, model, suffix))
        

class GatherThread(threading.Thread):
    def __init__(self, suffix, model):
        threading.Thread.__init__(self)
        self.suffix = suffix
        self.model = model

    def run(self):
        GatherFile(self.suffix, self.model)

def main():

    # establish and launch calculation threads
    threads_c = []
    for i in range(Model_Num):
        thread_i = CalThread(i, Model_Name)
        threads_c += [thread_i]
        thread_i.start()
        while( len(threads_c) == Max_Threads ):
            for j in range(Max_Threads):
                if( not threads_c[j].is_alive() ):
                    del threads_c[j]
                    break
        if( i + 1 == Model_Num ):
            for k in threads_c:
                k.join()

    # establish and launch gather threads
    threads_g = []
    for suffix in Extensions:
        thread_gather = GatherThread(suffix, Model_Name)
        threads_g = [thread_gather]
        thread_gather.start()
    for g in threads_g:
        g.join()
    # gather output
    stdoPath = './output/'
    if( not os.path.exists(stdoPath) ):
        os.makedirs(stdoPath)
    subprocess.call('mv *.out {}'.format(stdoPath), shell=True)

    error = []
    while(not Error_Queue.empty()):
        error += [Error_Queue.get()]
    error = list(set(error))

    if(len(error) != 0):
        print('Something is wrong, the input scripts which have problem is kept.')
    for i in range(Model_Num):
        if( i not in error):
            os.remove('./grid{0:09d}_{1}.in'.format(i, Model_Name))

    
    return 0

if __name__ == "__main__":
    main()
