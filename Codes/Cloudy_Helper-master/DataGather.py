#!/usr/bin/env python3

######################################################################### 
# DataGather.py can merge temporary files generate by 'grid' command.   # 
# This is used when Cloudy fails to merge those files together.         # 
# Synopsis:  ./DataGather.py (original_data_path) file1 file2 ...       #
# original_data_path is a optional argument, it must end with '/'       #
# its default value is './original_data/',                              # 
# file1, file2 and ... are the files you want to save.                  # 
# We don't suggest mergint standard output of Cloudy together,          # 
# since that file may very big.                                         # 
#########################################################################

import sys, threading, os, queue

Missing_Queue = queue.PriorityQueue()
Blank_Queue = queue.PriorityQueue()
Info_Queue = queue.Queue()

# 1 means labels locates at single file, 2 means labels are at grid files, 0 means labels are missing
def CheckColumnLabel( SubFilePath, FileName ):
    if( os.path.exists('{0}{1}'.format(SubFilePath, FileName))):
        File = open('{0}{1}'.format(SubFilePath, FileName), 'r')
        Line = File.readline()
        File.close()
        if( len(Line) != 0 ):
            return 1
    if( os.path.exists('{0}grid000000000_{1}'.format( SubFilePath, FileName ))):
        return 2
    else:
        return 0

 # LabelStatus comes from CheckColumnLabel, return number of files gathered and gathered data (a list of string)  
def GatherFile( SubFilePath, FileName, LabelStatus ):
    Data = []
    if( LabelStatus == 1 ):
        LabelFile = open('{0}{1}'.format(SubFilePath, FileName), 'r')
        LabelLine = LabelFile.readlines()
        LabelFile.close()
        for i in LabelLine:
            Data.append(i)

    FileExist = 1
    FileNum = 0
    FileSequence = 0
    CurrentExistFile = 0
    PreviousExistFile = -1
    Check = 0
                        
    while( FileExist == 1 ) :
        if( os.path.exists('{0}grid{1:09d}_{2}'.format(SubFilePath, FileSequence, FileName)) ):
            ## Check if files are missing
            CurrentExistFile = FileSequence
            while( (CurrentExistFile - PreviousExistFile) != 1 ):
                Data.append('\n')
                PreviousExistFile = 1 + PreviousExistFile
                Missing_Queue.put(['{}'.format(FileName),'grid{0:09d}_{1}'.format(PreviousExistFile, FileName)])
            PreviousExistFile = FileSequence
                
            SubFile = open('{0}grid{1:09d}_{2}'.format(SubFilePath, FileSequence, FileName), 'r')
            Lines = SubFile.readlines()
            SubFile.close()
            if (len(Lines) == 0):
                Data.append('\n')
                Blank_Queue.put(['{}'.format(FileName), 'grid{0:09d}_{1}'.format(FileSequence, FileName)])
            else:
                for LineNum in range(0, len(Lines)):
                    Data.append(Lines[LineNum])
            FileNum += 1
            Check = 0
        else:
            Check += 1
            if( Check == 1000 ):
                FileExist = 0
        FileSequence += 1

    return [FileNum, Data]

        
def ExportData( Data, FileName ):
    File = open('./{}'.format(FileName), 'w')
    for i in Data:
        File.write('{}'.format(i))
    File.close()        

class Process(threading.Thread):
    def __init__(self, SubFilePath, FileName):
        threading.Thread.__init__(self)
        self.SubFilePath = SubFilePath
        self.FileName = FileName

    def run(self):
        print('Thread \'{}\' start!'.format(self.FileName))

        if( self.FileName.find('.out') == -1):
            LabelStatus = CheckColumnLabel(self.SubFilePath, self.FileName)
        else:
            LabelStatus = 2
        GatheredInfo = GatherFile(self.SubFilePath, self.FileName, LabelStatus)

        if( GatheredInfo[0] == 0 ):
            print('ERROR: I cannot find any {} file, please check path!'.format(self.FileName))
            return 0
            
        if( LabelStatus == 0 ):
             Info_Queue.put('{0} column labels may be missed!\nCheck gathered file to see whether column labels exist.'.format(self.FileName))
             
        ExportData( GatheredInfo[1], self.FileName )
        print('Thread \'{0}\' finished! {1} files are merged together'.format(self.FileName, GatheredInfo[0]))
     
        return 0

        
def main():
    # Check whether arguments are provided
    if( len( sys.argv )<2 ):
        sys.exit('ERROR: I don\'t get all arguments I need!\nPlease provide file names.')

    # Default Data Path
    Path = './original_data/'

    # Get Argumetns
    Arguments = sys.argv[1:]
    Files = []

    for ArgvNum in range( len(Arguments) ):
        if( Arguments[ArgvNum].endswith('/') ):
            Path = Arguments[ArgvNum]
        else:
            Files.append( Arguments[ArgvNum] )

    if( not os.path.isdir(Path) ):
        sys.exit('ERROR: I cannot find where are the data files!\nPlease input path to folder which contains the data files and make sure it end with \'/\'.')
    
    threads = []
    for i in Files:
        thread_i = Process(Path,i)
        threads += [thread_i]
        thread_i.start()

    for t in threads:
        t.join()
        
    # Export error info
    if(not Missing_Queue.empty()):
        print('\nFollowing files are missing:')
    while(not Missing_Queue.empty()):
        MissFileTempList = Missing_Queue.get()
        print('{}'.format(MissFileTempList[1]))

    if(not Blank_Queue.empty()):
        print('\nFollowing files are blank:')
    while(not Blank_Queue.empty()):
        BlankFileTempList = Blank_Queue.get()
        print('{}'.format(BlankFileTempList[1]))

    if(not Info_Queue.empty()):
        print()
    while(not Info_Queue.empty()):
        print(Info_Queue.get())

if __name__ == '__main__':
    main()
