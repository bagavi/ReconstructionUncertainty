from Bio import SeqIO
import sys, numpy, math, datetime
import CommonFunctions
from CommonFunctions import WriteArrayinFile
from operator import itemgetter

class Genome:
    Name = ""
    DNAList= []
    DNA_current = []
    ReadLength_Considered = 100
    RepeatsHashDictionary = dict()
    Breakpoint = 100
    Filename = ""
    SideLengths = 1
    Gap = 1
    def __init__(self, Filename, Gap):
        self.Gap = Gap
        self.Filename = Filename
        Handle = open(Filename)
        for seq_record in SeqIO.parse(Handle, "fasta"):
            self.DNAList.append(str(seq_record.seq))
        print( "Length of the DNA is", len(self.DNAList[0]))
        
    def __IncreamentDictElement(self, key ):
        self.RepeatsHashDictionary[ key  ] = self.RepeatsHashDictionary.get(key, 0) + 1
    
    def __getDictValue(self, key ):
        return( self.RepeatsHashDictionary.get( tuple(key), 0) )
    
    def __getLength_L_read(self, position = 0):
        Extra = max( position -  len(self.DNA_current) + self.ReadLength_Considered, 0 )
        if Extra > 0 :
            return ( self.DNA_current[position:] + self.DNA_current[:Extra] )
        else:
            return ( self.DNA_current[position: position + self.ReadLength_Considered])
    
    def getUncertainty(self, MaxLength):
        Summary = []
        for length in range(1,MaxLength):
            Summary += [ self.RepeatsofLengthL(length) ]
        CommonFunctions.ReWriteArrayinFile(Summary, "Summary_"+self.Filename[:-6]+".csv")
#     
    def RepeatsofLengthL(self, length = 100):
        self.ReadLength_Considered = length
        self.RepeatsHashDictionary = dict()
        self.DNA_current = self.DNAList[0]
        Reads = []
        print("Read Length", length)
        for position in range(self.SideLengths, len(self.DNA_current) - self.ReadLength_Considered - self.SideLengths - 1872915):
            if position % 1000000 == 0:
                print("In position", position, "DNA left", len(self.DNA_current) - position )           
            Reads += [ [ self.DNA_current[position - self.SideLengths:                             position] ,
                         self.DNA_current[position                   :position + self.ReadLength_Considered] ,
                         self.DNA_current[position + self.ReadLength_Considered : position + self.ReadLength_Considered + self.SideLengths],
                         position
                         ] ]
        print("Sorting Read Array")
        Reads.sort(key=itemgetter(1), reverse=False)
        CountInfo = []
        ReadInfo = []
        print("Counting Repeats")
        CurrentString = Reads[0][1]
        Repeat = 0
        Titles = [ "Gene", "Gene Length", "Read Length", "Read", "Neighbhours", "Repeat"]
        WriteInfo = [ [Titles ] ]
        LeftNeighbhors = []
        RightNeighbors = []
        counter = 0
        TempReads = []
        for read in Reads[1:]:
            counter +=1
#             if counter%50000 == 0:
#                 print("Counter", counter, "Reads", len(Reads) )
            if read[1] == CurrentString:
                Repeat += 1
                LeftNeighbhors += read[0]
                RightNeighbors += read[2]
            else:
                if Repeat > 0:
                    if len( set(LeftNeighbhors) ) != 1 or len( set(RightNeighbors) ) != 1:
#                         print("YAY", "LeftNeighbhours", set(LeftNeighbhors), "Read", read, "Right Neighbors", set(RightNeighbors) )
                        CountInfo += [ Repeat *math.log(len( set(RightNeighbors) ), 2)]
                        ReadInfo += [ CurrentString ]
                        WriteInfo += [ [ self.Filename, len(self.DNA_current), self.ReadLength_Considered, read[1], -1, Repeat ] ]
#                 elif len( set(LeftNeighbhors) ) == 1 and len( set(RightNeighbors) ) == 1 and Repeat > 0:
#                     print( "Saved!" ,Repeat)
                
                LeftNeighbhors = []
                RightNeighbors = []
                Repeat = 0
                CurrentString = read[1]
     #   WriteArrayinFile(TempReads, "TempRead.csv")                
        CountInfo.sort(key=None, reverse=True)
        CommonFunctions.ReWriteArrayinFile(WriteInfo,'RepeatData.csv')
        Summary = [ self.Filename, len(self.DNA_current), self.ReadLength_Considered, len(CountInfo), sum(CountInfo), datetime.datetime.now() ]
        if False:
            print("Length of the DNA is", len(self.DNA_current))
            print("Number of reads repeating of length", self.ReadLength_Considered," is", len(CountInfo))
            print("Time", datetime.datetime.now())
        print("Uncertainity", sum(CountInfo))
        
        return(Summary)
if len( sys.argv )> 1 :
    MaxReadLength = int( float( sys.argv[1] ) )
else:
    MaxReadLength = int( 10)
    
if len( sys.argv )> 2 :
    Gap = int( float( sys.argv[2] ) )
else:
    Gap = int(5)

Staphylococcus = Genome("StaphylococcusAureus.fasta", Gap)
Staphylococcus.getUncertainty(MaxReadLength)

# Rhodobacter = Genome("RhodobacterSphaeroides.fasta")
# 
#     def RepeatsofgivenLength(self, length = 100):
#         self.ReadLength_Considered = length
#         self.RepeatsHashDictionary = dict()
#         self.DNA_current = self.DNAList[0]
#         unit = 1
#         for position in range( len(self.DNA_current) ):
#             if position%5000 == 0:
#                 print("In position", position, "Unit", unit, "Breakpoint", self.Breakpoint, "DNA left", len(self.DNA_current) - position )
#                 unit +=1
#                 if unit == self.Breakpoint:
#                     break
#             Read = tuple( self.__getLength_L_read(position))
#             ReadString = ''.join(Read)
#             self.__IncreamentDictElement(ReadString)
#         
#         self.AnalyzeRepeatStructure()
# Rhodobacter.RepeatsofgivenLength(100)