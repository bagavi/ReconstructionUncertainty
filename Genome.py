from Bio import SeqIO
import sys, numpy, math, datetime
import CommonFunctions
from CommonFunctions import WriteArrayinFile
from operator import itemgetter
from collections import Counter

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
    
    def getUncertainty(self):
        Summary = []
        length = 88
        while True:
            length += self.Gap
            Answer = self.RepeatsofLengthL(length)
            print("Answer", Answer)
            Summary += [ Answer ]
            if Answer[-3] < 1:
                break
        CommonFunctions.WriteArrayinFile(Summary, "Precise2_Summary_"+self.Filename[:-6]+".csv")
#   
    

    def factlog(self,n):
        return( n*(math.log(n,2) - 1.44269) + 0.5*math.log(2*math.pi*n)) #math.log(math.e,2) = 1.4426950408889634
#         for i in range(1,n):
#             sum += math.log(i,2)
#         return sum
    
    def RepeatsofLengthL(self, length = 100):
        self.ReadLength_Considered = length
        self.RepeatsHashDictionary = dict()
        self.DNA_current = self.DNAList[0]
        print(len(self.DNA_current))
      #  a = input("FFF")
        Reads = []
        print("Read Length", length)
        for position in range(self.SideLengths, len(self.DNA_current) - self.ReadLength_Considered - self.SideLengths ):
            if position % 1000000 == 0:
                print("In position", position, "DNA left", len(self.DNA_current) - position )           
            Reads +=[   [ 
                         self.DNA_current[position - self.SideLengths:                             position] ,
                         self.DNA_current[position                   :position + self.ReadLength_Considered] ,
                         self.DNA_current[position + self.ReadLength_Considered : position + self.ReadLength_Considered + self.SideLengths],
                         position
                        ]
                    ]
        print("Sorting Read Array")
        Reads = sorted(Reads,key=itemgetter(1) )
        CountInfo = []
        ReadInfo = []
        print("Counting Repeats")
        CurrentString = Reads[0][1]
        Repeat = 0
        Titles = [ "Gene", "Gene Length", "Read Length", "Read", "Neighbhours", "Repeat"]
        Is_less_than_critical_length = False
        Reason = "No reason"
        Position_of_repeat_less_than_2 = []
        LeftNeighbhors = []
        RightNeighbors = []
        RepeatPositions = []
        Repeat = 0

        for read in Reads[1:]:
            if read[1] == CurrentString:
                Repeat += 1
                LeftNeighbhors += read[0]
                RightNeighbors += read[2]
                RepeatPositions += [ read[3] ]
            else:
                if Repeat > 1 and len( set(RightNeighbors) ) != 1:
                    #Checking for l_critical
                    if Repeat > 2:
                        Is_less_than_critical_length = True
                        Reason = "Triple Repeats"
                    else:
                        Position_of_repeat_less_than_2 += RepeatPositions
                        Reason = "Interleaved Repeats"
                    Count_stats = Counter(RightNeighbors).values()
                    try:
                        Uncertainty = self.factlog(sum(Count_stats))
                        for i in Count_stats:
                            Uncertainty -= self.factlog(i)
                    except:
                        Uncertainty = 2*(sum(Count_stats))
                    CountInfo += [ Uncertainty ]
                LeftNeighbhors = []
                RightNeighbors = []
                RepeatPositions = []
                Repeat = 0
            
                CurrentString = read[1]
        
        if Position_of_repeat_less_than_2 != sorted(Position_of_repeat_less_than_2):
#            print( Position_of_repeat_less_than_2, sorted(Position_of_repeat_less_than_2))
            Is_less_than_critical_length = True
        Summary = [ self.Filename, len(self.DNA_current), self.ReadLength_Considered, len(CountInfo), sum(CountInfo), Is_less_than_critical_length, str(datetime.datetime.now()) ]
        if False:
            print("Length of the DNA is", len(self.DNA_current))
            print("Number of reads repeating of length", self.ReadLength_Considered," is", len(CountInfo))
            print("Time", datetime.datetime.now())
        print("Uncertainity", sum(CountInfo), "Is less than critical length", Is_less_than_critical_length, "reason", Reason, "Repeat position", Position_of_repeat_less_than_2)
        return(Summary)

filename = "Buchnera_aphidicola.fasta"
Gene = Genome(filename, 3)
Gene.getUncertainty()


#RhodobacterSphaeroides = Genome("RhodobacterSphaeroides.fasta", 3)
#RhodobacterSphaeroides.getUncertainty()
