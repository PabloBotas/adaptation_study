#!/usr/bin/python3

# import tramp_spot

class Tramp:
    """
    Class for tramp files.
    
    The __init__ constructor fills all the fields from the file string.
    """
    
    ## Init method / constructor
    def __init__(self, file):
        # print('Calling __init__')
        self.file = str(file)
        self.patient_id             = str()
        self.patient_first_name     = str()
        self.patient_middle_initial = str()
        self.patient_last_name      = str()
        self.astroid_id             = str()
        self.course_name            = str()
        self.beam_name              = str()
        self.gantry                 = float()
        self.couch_rotation         = float()
        self.gigaproton_total       = float()
        self.number_spots           = int()
        self.nlayers                = int()
        self.phases                 = list()
        self.spotspacing            = float()
        self.spots                  = list()
        self.read()
    
    def read(self):
        with open(self.file, 'r') as f:
            for i, line in enumerate(f):
                try:
                    line = line.strip('\n')
                    line = line.lstrip()
                    if line.startswith('#'):
                        if "patient_i"  in line: self.patient_id             = line.split(' ')[-1]
                        if "patient_f"  in line: self.patient_first_name     = line.split(' ')[-1]
                        if "patient_m"  in line: self.patient_middle_initial = line.split(' ')[-1]
                        if "patient_l"  in line: self.patient_last_name      = line.split(' ')[-1]
                        if "astroid"    in line: self.astroid_id             = line.split(' ')[-1]
                        if "course"     in line: self.course_name            = line.split(' ')[-1]
                        if "beam"       in line: self.beam_name              = line.split(' ')[-1]
                        if "gantry"     in line: self.gantry                 = float(line.split(' ')[-1])
                        if "couch"      in line: self.couch_rotation         = float(line.split(' ')[-1])
                        # if "gigaproton" in line: self.gigaproton_total       = float(line.split(' ')[-1])
                        # if "rows_total" in line: self.number_spots           = int(line.split(' ')[-1])
                    else:
                        time   = float(0)
                        phase  = float(0)
                        energy = float(line.split()[0]) # MeV
                        x      = float(line.split()[1]) # mm
                        y      = float(line.split()[2]) # mm
                        ngp    = float(line.split()[3]) # gigaprotons (10^9 protons)
                        self.gigaproton_total += ngp
                        self.number_spots     += 1
                        self.spots.append(Spot(energy, x, y, ngp, time, phase))
                except ValueError:
                    print("ValueError on line {0}:".format(i+1))
                    print("\t{0}".format(line))
                    raise
    
    def getHeader(self, nspots=None, ngps=None):
        if not nspots or not ngps:
            nspots = self.number_spots
            ngps = self.gigaproton_total
        header = "# patient_id {0}".format(self.patient_id) + "\n" + \
                 "# patient_first_name {0}".format(self.patient_first_name) + "\n" + \
                 "# patient_middle_initial {0}".format(self.patient_middle_initial) + "\n" + \
                 "# patient_last_name {0}".format(self.patient_last_name) + "\n" + \
                 "# astroid_id {0}".format(self.astroid_id) + "\n" + \
                 "# course_name {0}".format(self.course_name) + "\n" + \
                 "# beam_name {0}".format(self.beam_name) + "\n" + \
                 "# gantry {0}".format(self.gantry) + "\n" + \
                 "# couch_rotation {0}".format(self.couch_rotation) + "\n" + \
                 "# gigaproton_total {0}".format(ngps) + "\n" + \
                 "# rows_total {0}".format(nspots) + "\n" + \
                 "# E(MeV) X(mm) Y(mm) N(Gp)\n"
        return header
    
    
    def getPhases(self):
        temp = list()
        for i in range(0, self.number_spots):
            if self.spots[i].phase not in temp:
                temp.append(self.spots[i].phase)
        self.phases = temp
        return self.phases
    
    def setNumberLayers(self, num=None):
        if num:
            try:
                self.nlayers = int(num)
            except ValueError:
                print('Please, enter a number')
                raise
        else:
            self.nlayers = 1
            for i in range(1, self.number_spots):
                if self.spots[i].energy != self.spots[i-1].energy:
                    self.nlayers += 1
    
    def getNumberLayers(self):
        if not self.nlayers:
            self.setNumberLayers()
        return self.nlayers
    
    
    ## String transformation overloading
    def __str__(self):
        # print('Calling __str__')
        out = "N. of spots: {0}".format(self.number_spots) + "\n"\
              "Gigaprotons: {0}".format(self.gigaproton_total) + "\n"\
              "First and last spot:" + "\n"\
              "\t{0}".format(self.spots[0]) + "\n"\
              "\t{0}".format(self.spots[-1])
        return out


class Spot:
    """
    Class for Spot.
    
    A spot is defined inside a tramp file as a line specifying:
        - Energy (MeV)
        - X position (mm)
        - Y position (mm)
        - Number of gigaprotons (10^9 protons)
    This class additionally adds support for time information, initializing it to 0.
    Provides a custom conversion to string for better printing.
    """
    
    ## Init method / constructor
    def __init__(self, energy, x, y, ngp, time=0, phase=0):
        # print('Calling __init__')
        self.energy = float(energy) # MeV
        self.x      = float(x) # mm
        self.y      = float(y) # mm
        self.ngp    = float(ngp) # gigaprotons (10^9 protons)
        self.time   = float(time)
        self.phase  = int(phase)
    
    ## Initialization from line
    @classmethod
    def fromline(cls, line):
        """
        Method to initialize a line
        """
        time   = float(0)
        phase  = float(0)
        line   = line.strip('\n')
        energy = float(line.split('\t')[0]) # MeV
        x      = float(line.split('\t')[1]) # mm
        y      = float(line.split('\t')[2]) # mm
        ngp    = float(line.split('\t')[3]) # gigaprotons (10^9 protons)
        return cls(energy, x, y, ngp, time, phase)
    
    ## String transformation overloading
    def __str__(self):
        # print('Calling __str__')
        return "energy = {0} x = {1} y = {2} ngp = {3} time = {4} phase = {5}".format(self.energy, self.x, self.y, self.ngp, self.time, self.phase)

