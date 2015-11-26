from Bio import SeqIO

class Genome:
    Name = ""
    DNAList= []
    DNA_current = []
    ReadLength_Considered = 100
    RepeatsHashDictionary = dict()
    
    def __init__(self, filename=""):
        Handle = open(filename)
        for seq_record in SeqIO.parse(Handle, "fasta"):
            self.DNAList.append(seq_record.seq)
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
        
    def RepeatsofgivenLength(self, length = 100):
        self.ReadLength_Considered = length
        self.RepeatsHashDictionary = dict()
        self.DNA_current = self.DNAList[0]
        for position in range( len(self.DNA_current) ):
            if position%5000 == 0:
                print(" In position", position)
            Read = tuple( self.__getLength_L_read(position))
            ReadString = ''.join(Read)
            self.__IncreamentDictElement(ReadString)
        
        self.AnalyzeRepeatStructure()
        
    def AnalyzeRepeatStructure(self):
        RepeatValues = list( self.RepeatsHashDictionary.values())
        RepeatValues.sort()
Staphylococcus = Genome("StaphylococcusAureus.fasta")
Staphylococcus.RepeatsofgivenLength(100)

# Rhodobacter = Genome("RhodobacterSphaeroides.fasta")
# Rhodobacter.RepeatsofgivenLength(100)