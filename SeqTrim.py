#!/usr/bin/env python

import sys

class Convert(object):
    """
        The class Convert provides functions for converting between different 
        Phred quality score encodings and allows retrieving ASCII characters
        corresponding to integer Phred scores.
    """
    def __init__(self):
        # Map to associate raw Phred score to ASCII+33 (Sanger encoding)
        self.Phred33 = {
                0:'!',
                1:'"',
                2:'#',
                3:'$',
                4:'%',
                5:'&',
                6:"'",
                7:'(',
                8:')',
                9:'*',
                10:'+',
                11:',',
                12:'-',
                13:'.',
                14:'/',
                15:'0',
                16:'1',
                17:'2',
                18:'3',
                19:'4',
                20:'5',
                21:'6',
                22:'7',
                23:'8',
                24:'9',
                25:':',
                26:';',
                27:'<',
                28:'=',
                29:'>',
                30:'?',
                31:'@',
                32:'A',
                33:'B',
                34:'C',
                35:'D',
                36:'E',
                37:'F',
                38:'G',
                39:'H',
                40:'I',
                41:'J',
                42:'K',
                43:'L',
                44:'M',
                45:'N',
                46:'O',
                47:'P',
                48:'Q',
                49:'R',
                50:'S',
                51:'T',
                52:'U',
                53:'V',
                54:'W',
                55:'X',
                56:'Y',
                57:'Z',
                58:'[',
                59:'\\',
                60:']',
                61:'^',
                62:'_',
                63:'`',
                64:'a',
                65:'b',
                66:'c',
                67:'d',
                68:'e',
                69:'f',
                70:'g',
                71:'h',
                72:'i',
                73:'j',
                74:'k',
                75:'l',
                76:'m',
                77:'n',
                78:'o',
                79:'p',
                80:'q',
                81:'r',
                82:'s',
                83:'t',
                84:'u',
                85:'v',
                86:'w',
                87:'x',
                88:'y',
                89:'z',
                90:'{',
                91:'|',
                92:'}',
                93:'~'
        }

        # Map to associate raw ASCII+64 to phred score
        self.Phred64 = {
                '@':0,
                'A':1,
                'B':2,
                'C':3,
                'D':4,
                'E':5,
                'F':6,
                'G':7,
                'H':8,
                'I':9,
                'J':10,
                'K':11,
                'L':12,
                'M':13,
                'N':14,
                'O':15,
                'P':16,
                'Q':17,
                'R':18,
                'S':19,
                'T':20,
                'U':21,
                'V':22,
                'W':23,
                'X':24,
                'Y':25,
                'Z':26,
                '[':27,
                '\\':28,
                ']':29,
                '^':30,
                '_':31,
                '`':32,
                'a':33,
                'b':34,
                'c':35,
                'd':36,
                'e':37,
                'f':38,
                'g':39,
                'h':40,
                'i':41,
                'j':42,
                'k':43,
                'l':44,
                'm':45,
                'n':46,
                'o':47,
                'p':48,
                'q':49,
                'r':50,
                's':51,
                't':52,
                'u':53,
                'v':54,
                'w':55,
                'x':56,
                'y':57,
                'z':58,
                '{':59,
                '|':60,
                '}':61,
                '~':62
        }

    def Convert(self, QualLine):
        """ Converts Illumina >1.3<1.8 phred encoding to Sanger format. """
        NewQualLine = []
        for Value in QualLine:
            try:
                Key = self.Phred64[Value]
                NewQualLine.append(self.Phred33[Key])
            except KeyError:
                sys.exit('Invalid Format. Are you sure you are trying to convert between the right formats?')
        return ''.join(NewQualLine)
    
    def PhredToASCII(self, PhredScore):
        """ Takes phred score as argument and returns corresponding ASCII character (Sanger encoding). """
        try:
            return self.Phred33[PhredScore]
        except KeyError:
            sys.exit('Cannot retrieve ASCII character. Valid integers range from 0 to 93!')

class Trim(Convert):
    """ 
        Methods for trimming fastq sequences and corresponding quality scores from 5 and 3 
        prime ends. To instantiate, provide 3 arguments, Nucleotide sequence and correspondig
        quality string in Sanger Phred encoding, and minimum Phred score to accept. 
    """
    def __init__(self, Sequence, Quality, QScore):
        """ Instantiation of the class requires 3 arguments:
        Sequence, Quality line, and minimum acceptable Phred score. """
        try:
            assert type(Sequence) == str
            self.Sequence = Sequence
        except AssertionError:
            sys.exit('Sequence not a string!')
        try:
            assert type(Quality) == str
            self.Quality = Quality
        except AssertionError:
            sys.exit('Quality line not a string!')
        try:
            assert 0 <= QScore <= 93 and type(QScore) == int
            super(Trim, self).__init__()
            self.QScore = super(Trim, self).PhredToASCII(QScore)
#            THIS Works as an alternative to the 2 lines above:
#            self.Convert = Convert()
#            self.QScore = self.Convert.PhredToASCII(QScore)
        except AssertionError:
            self.QScore = QScore
            sys.exit('Quality score not an integer or out of range (0->93)!')

    def FivePrime(self, Crawl=0):
        """ Trim from the 5 prime end of a sequence. Trims until first nucleotide passing threshold is encountered.
        By setting Crawl to be a positive integer the search can be extended (eg, Crawl=5 will make the search continue
        for 5 nucleotides downstream of the first one that passed phred scores threshold). """
        try:
            assert Crawl >= 0 and Crawl <= len(self.Quality)
            if Crawl == 0:
                for i in range(len(self.Quality)):
                    if self.Quality[i] >= self.QScore:
                        Start = i + 1
                        break
                self.Sequence = self.Sequence[i:]
                self.Quality = self.Quality[i:]
            else:
                # Initialize Start and Trim, so that they are available
                Start = 0
                Trim = Start
                # Find the first base below Phred threshold
                for i in range(len(self.Quality)):
                    if self.Quality[i] >= self.QScore:
                        Start = i + 1
                        break
                # Move beyond first base that failed test and do a more
                # aggressive trim
                Count = 0
                Trim = Start - 1
                for j in self.Quality[Start:(Crawl + Start)]:
                    Count += 1
                    if j < self.QScore:
                        Trim = (Start + Count)
                    elif Count == Crawl:
                        break
                # Trim the sequence
                self.Sequence = self.Sequence[Trim:]
                self.Quality = self.Quality[Trim:]
        except AssertionError:
            if Crawl > len(self.Quality):
                sys.exit('Crawl variable passed to function must be <= length of sequence.')
            else:
                sys.exit('Crawl variable passed to function must be >= 0 and <= length of sequence.')

    def ThreePrime(self, Crawl=0):
        """ Trim from the 3 prime end of a sequence. Trims until first nucleotide passing threshold is encountered.
        By setting Crawl to be a positive integer the search can be extended (eg, Crawl=5 will make the search continue
        for 5 nucleotides upstream of the first one that passed phred scores threshold). """
        try:
            assert Crawl >= 0 and Crawl <= len(self.Quality)
            if Crawl == 0:
                for i in range(len(self.Quality)-1,-1,-1):
                    if self.Quality[i] >= self.QScore:
                        Start = i + 1
                        break
                self.Sequence = self.Sequence[0:(i + 1)]
                self.Quality = self.Quality[0:(i + 1)]
            else:
                # Initialize Start and Trim, so that they are available
                Start = (len(self.Quality))
                Trim = Start
                # Find the first base below Phred threshold
                for i in range(len(self.Quality)-1,-1,-1):
                    if self.Quality[i] >= self.QScore:
                        Start = (len(self.Quality) - i) + 1
                        break
                # Move beyond first base that failed test and do a more
                # aggressive trim
                Count = 0
                Trim = Start - 2
                RevSeq = self.Sequence[::-1]
                RevQual = self.Quality[::-1]
                for j in RevQual[Start:(Start + Crawl)]:
                    Count += 1
                    if j < self.QScore:
                        Trim = (Start + Count)
                    elif Count == Crawl:
                        break
                # Trim the sequence
                self.Sequence = RevSeq[Trim:][::-1]
                self.Quality = RevQual[Trim:][::-1]
        except AssertionError:
            if Crawl > len(self.Quality):
                sys.exit('Crawl variable passed to function must be <= length of sequence.')
            else:
                sys.exit('Crawl variable passed to function must be >= 0 and <= length of sequence.')
    
    def GlobalTrim(self, PercentBases=1):
        """ Counts the number of nucleotides below QScore threshold. PercentBases specifies the
        acceptable fraction of nucleotides below the Phred score threshold. Sequence is clipped
        5-3 prime from first unacceptable nucleotide if too many nuclotides fail quality check. """
        try:
            assert PercentBases >= 0 and type(PercentBases) == int
            Count = 0
            Trim = None
            # Check how many bases fail QScore test
            for i in range(len(self.Quality)):
                if self.Quality[i] < self.QScore:
                    Count += 1
                    if Trim == None:
                        Trim = i
            if Count >= (Count/PercentBases*100):
                self.Sequence = self.Sequence[:Trim]
                self.Quality = self.Quality[:Trim]
        except AssertionError:
            sys.exit('Provide number of acceptable nucleotides below Phred threshold as integer argument.')

    def MinLength(self, Length=None):
        """ Sets Sequence and Quality line to None if they are below threshold length. This may break
        your downstream pipe if you are not careful how you handle returns from this function. Safeguard 
        by using a statement like "if None in trim.Retrieve(): pass" when using MinLength function. """
        try:
            assert Length >= 0 and type(Length) == int
            if len(self.Sequence) < Length:
                self.Sequence = None
                self.Quality = None
        except AssertionError:
            if type(Length) != int:
                sys.exit('Provide min length as integer argument.')
            else:
                sys.exit('Min length needs to be > 0.')
        
    def Retrieve(self):
        """ Retrieve the sequence and quality line. Returns a list. """
        try:
            if len(self.Sequence) == 0:
                self.Sequence = None
                self.Quality = None
            return self.Sequence, self.Quality
        except TypeError:
            return self.Sequence, self.Quality
