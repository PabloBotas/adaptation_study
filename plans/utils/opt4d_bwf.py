#!/usr/bin/python3

class Bwf:
    """
    Class for bwf files.
    
    The __init__ constructor fills all the fields from the file string.
    """
    
    ## Init method / constructor
    def __init__(self, file):
        # print('Calling __init__')
        self.file = str(file)
        self.totalWeight = float()
        self.number_spots = int()
        self.w = list()
        self.x = list()
        self.y = list()
        self.read()

    def read(self):
        with open(self.file, 'r') as f:
            for i, line in enumerate(f):
                try:
                    line = line.strip('\n')
                    line = line.lstrip()
                    if not line.startswith('#'):
                        _,w,_,_,_,_,_,x,y  = [float(j) for j in line.split()]
                        self.w.append(w)
                        self.x.append(x)
                        self.y.append(y)
                        self.totalWeight += w
                        self.number_spots += 1
                except ValueError:
                    print("ValueError on line {0}:".format(i+1))
                    print("\t{0}".format(line))
                    raise

    ## String transformation overloading
    def __str__(self):
        # print('Calling __str__')
        out = "N. of spots: {0}".format(self.number_spots) + "\n"\
              "totalWeight: {0}".format(self.totalWeight) + "\n"
        return out

