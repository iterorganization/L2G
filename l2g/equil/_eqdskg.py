import os
import re

import logging

log = logging.getLogger(__name__)

class NanValuesInEqdskFile(Exception):
    pass

class FileDoesNotExist(Exception):
    pass

class EQDSKIO(object):
    """EQDSK interface implementation

    The following quantities are read from a ``EQDSK G``.

    Attributes:
        PSIRZ (arr): Array of plasma psi values in units Web/rad.
        CURRENT (float): Plasma current in Ampere
        RDIM (float): Horizontal dimension in meter of computational box
        ZDIM (float): Vertical dimension in meter of computational box
        NW (int): Number of horizontal R grid points
        NH (int): Number of vertical Z grid points
        LIMITR (int): Number of limiter points
        RLIM (arr): R of surrounding limiter contour in meter
        ZLIM (arr): Z of surrounding limiter contour in meter
        NBBBS (int): Number of boundary points
        RBBBS (arr): R of boundary points in meter
        ZBBBS (arr): Z of boundary points in meter
        RMAXIS (float): R of magnetic axis in meter
        ZMAXIS (float): Z of magnetic axis in meter
        FPOL (arr): Poloidal current function in m-T, F=RB_T on flux grid
        PRES (arr): Plasma pressure in nt/m^2 on uniform flux grid
        FFPRIM (arr): FF'(Psi) in (mT)^2/(Weber/rad) on uniform flux grid
        PPRIME (arr): P'(Psi) in (nt/m^2) / (Weber/rad) on uniform flux grid
        QPSI (arr): q values on uniform flux grid from axis to boundary
        SIMAG (arr): poloidal flux at magnetic axis in Weber /rad
        SIBRY (arr): poloidal flux at the plasma boundary in Weber/rad
        RCENTR (float): R in meter of vacuum toroidal magnetic field BCENTR
        BCENTR (float): Vacuum toroidal magnetic field in Tesla at RCENTR
        RLEFT (float): Minimum R in meter of rectangular computational box
        ZMID (float): Z of center of computational box in meter

    They are all accessible with the following example

    .. code-block:: python

       import l2g.equil
       x = l2g.equil.EQDSKIO()
       eqdskFile = '/path/to/eqdsk/format/file'
       # If SMITER is located on your machine
       x.openAndRead(eqdskFile)

       # If SMITER is used remotely you have to manually pass the contents
       # of the EQDSK file. In this case:
       with open(eqdskFile) as f:
           data = f.read()
       ok = x.read(data)
       # Then set the name of the eqdsk file
       if ok:
           # Reading Successful
           eqdskName = eqdskFile.rpslit('/', 1)[-1] # Get the name of the file
           x.setName(eqdskFile)
       else:
           print('Trouble reading EQDSK file %s' % eqdskFile)

       # If the file was read without a problem, you can access, i.e, the
       # geometry of the limiter
       x.RLIM
       x.ZLIM

    """

    # Patterns for regexing

    # Search for Not A Number expressions.
    nanPatt = re.compile(r'(?i)\bnan\b')

    # Integer pattern for reading values separated with space
    integerPattern = r'(?:^|(?<=\s))-?[0-9]+(?:$|(?=\s))'
    intPatt = re.compile(f'.*{integerPattern}.*')

    # Pattern for reading scientific numbers
    scientificPattern = r'[+-]?[0-9]\.[0-9]{6,9}[eE][+-][0-9]{2}'

    genPatt = re.compile(f'(?:{integerPattern}|{scientificPattern})')


    def __init__(self, file=''):
        self.name = ''
        self.resetValues()

        if file:
            ok = self.openAndRead(file)
            if not ok:
                log.error(f'Failed to load file: {file}')

    def resetValues(self) -> None:
        self.INPUT_EQDSK = ''
        self.HEADER = ''

        self.PSIRZ = [] # Plasma psi WEBER

        self.CURRENT: float = 0

        self.RDIM: float = 0 # Size dimension for R and Z
        self.ZDIM: float = 0

        self.NW: int = 0  # Dimension of grid
        self.NH: int = 0

        self.LIMITR: int = 0  # Number of limiter points
        self.RLIM: list = []  # rlim
        self.ZLIM: list = []  # zlim

        self.NBBBS: int = 0 # Number of Boundary points
        self.RBBBS: list = []  # rbbbs
        self.ZBBBS: list = []  # zbbs

        self.RMAXIS: float = 0 # Position of toroidal magnetic
        self.ZMAXIS: float = 0 # field
        self.FPOL: list = []
        self.PRES: list = []
        self.FFPRIM: list = []
        self.PPRIME: list = []
        self.QPSI: list = []

        self.SIMAG: float = 0  # Flux value at magnetix axis
        self.SIBRY: float = 0  # Flux value on boundary

        self.RCENTR: float = 0  # Reference value: position R, mag field B
        self.BCENTR: float = 0  # B value in RCENTR

        self.RLEFT: float = 0  # Minimum of rectangular computational box.
        self.ZMID: float = 0  # Vertical dimension of computational box
        """
        Or in another case:

        The grid goes from :

        RLEFT to RLEFT + RDIM
        -ZDIM / 2  to  ZDIM / 2
        """

        self.eqdskString = ''
        self.successfullRead = False

    def getName(self) -> str:
        """Returns the name.

        Returns:
            name (str): String which has whites-paces replaced with `_`.
        """
        return self.name.replace(' ', '_')

    def setName(self, name: str) -> None:
        """Set a name for eqdsk.. It should be the same as the name displayed
        to the **Salome Object** in **Object browser**.

        Arguments:
            name (str): Name for the case.
        """
        self.name = name

    def openAndRead(self, file_path: str) -> bool:
        """Opens and reads the contents of file at file_path.

        Arguments:
            file_path (str): Path to EQDSK-G format file.
        """
        if not os.access(file_path, os.F_OK | os.R_OK):
            raise FileDoesNotExist(f'File {file_path} does not exist or lacking permission to read it!')

        data = ''
        with open(file_path, 'r') as f:
            data = f.read()
        # Setting name
        self.setName(os.path.basename(file_path))
        return self.read(data)

    def read(self, eqdskString: str) -> bool:
        """Passing string instead of a filename, to avoid any troubles if this
        is used on a server.

        This function parses EQDSK G format.
        """
        self.resetValues()
        self.eqdskString = eqdskString
        # Reset the log

        LINES = eqdskString.splitlines()

        HEADER = LINES[0]
        self.HEADER = HEADER
        # regexHeader = '^.+[0-9] +?([0-9]+) +?([0-9]+)'
        regexHeader = r'\s[0-9]\s+([0-9]+)\s+([0-9]+)'

        result = re.findall(regexHeader, HEADER)
        # self.log("Reading 1st line of EQDSK.")
        if result:
            nw, nh = result[0]
            self.NW = int(nw)
            self.NH = int(nh)
        else:
            return False
        # self.log('Done reading 1st line.')
        # self.log("Reading 2nd line of EQDSK")
        i = 1
        result = self.readOneLine(LINES[i])

        if not result:
            # self.log("[INFO] 2nd line must be a comment")
            while 1:
                # self.log("[INFO] Filtering the comments.")
                i += 1

                result = self.readOneLine(LINES[i])
                if result:
                    break

                # Precaution, if there are 100 wrong lines, just end the
                # reading
                if i > 100:
                    return False
        try:
            self.RDIM = result[0]
            self.ZDIM = result[1]
            self.RCENTR = result[2]
            self.RLEFT = result[3]
            self.ZMID = result[4]

            i += 1
            result = self.readOneLine(LINES[i])
            self.RMAXIS = result[0]
            self.ZMAXIS = result[1]
            self.SIMAG = result[2]
            self.SIBRY = result[3]
            self.BCENTR = result[4]

            i += 1
            result = self.readOneLine(LINES[i])
            self.CURRENT = result[0]

            msg = "The 2nd value on line 4 of the EQDSK is not equal to " \
                  "the flux at magnetic axis read from the 3rd line"
            assert self.SIMAG == result[1], msg
            DUMMY = result[2]
            msg = "The 4th value on line 4 of the EQDSK is not equal to the"\
                  " magnetic axis R position, read from the 3rd line!"
            assert self.RMAXIS == result[3], msg
            DUMMY = result[4]
            # self.log("Done reading 4th line.")
            i += 1
            # self.log("Reading 5th line of EQDSK")
            result = self.readOneLine(LINES[i])
            msg = "The 1st value on line 5 is not equal to the magnetic " \
                  " axis Z position, read from the 3rd line!"
            assert self.ZMAXIS == result[0], msg
            DUMMY = result[1]
            msg = " The 3rd value on line 5 is not equal to the flux at " \
                  " boundary, read from the 3rd line"
            assert self.SIBRY == result[2], msg
            DUMMY = result[3]
            DUMMY = result[4]
            # self.log("Done reading 5th line.")

            line_index = i + 1
            P1, P2 = self.readValues(LINES[line_index:])

            _I = 0
            self.FPOL = P1[_I: _I + self.NW]
            _I += self.NW

            self.PRES = P1[_I: _I + self.NW]
            _I += self.NW

            self.FFPRIM = P1[_I: _I + self.NW]
            _I += self.NW

            self.PPRIME = P1[_I: _I + self.NW]
            _I += self.NW

            self.PSIRZ = []
            for i in range(self.NH):
                self.PSIRZ.append(P1[_I: _I + self.NW])
                _I += self.NW

            self.QPSI = P1[_I: _I + self.NW]
            _I += self.NW

            _I = 0
            self.NBBBS = int(P2[_I])
            _I += 1

            self.LIMITR = int(P2[_I])
            _I += 1

            self.RBBBS = []
            self.ZBBBS = []
            for i in range(self.NBBBS):
                self.RBBBS.append(P2[_I])
                self.ZBBBS.append(P2[_I + 1])
                _I += 2

            self.RLIM = []
            self.ZLIM = []
            for i in range(self.LIMITR):
                self.RLIM.append(P2[_I])
                self.ZLIM.append(P2[_I + 1])
                _I += 2
        except NanValuesInEqdskFile as e:
            log.error("NaN values in EQDSK file!")
            self.successfullRead = False
            return False
        except Exception as e:
            log.error("Failed to read EQDSK")
            self.successfullRead = False
            return False
        self.successfullRead = True
        return True

    def readOneLine(self, line: str) -> list:
        """Reads the float values from one line and returns a list of 5
        elements. If the line does not have 5 elements, it is considered a
        comment line and returns an empty list.

        Returns:
            out (Union[int, List[Float]]): Returns -1 if a line contains an
                alpha at the beginning, else it returns a list of values.
        """

        if line.lstrip()[0].isalpha():
            return []

        pattern = r'-?[0-9]{1}\.[0-9]{9}[eE][-+]?[0-9]{2}'
        result = re.findall(pattern, line)

        if len(result) == 5:
            # 5 float values on line.
            out = [float(el) for el in result]
            return out
        else:
            # Manually check the values
            out = []
            for el in line.split():
                try:
                    out.append(float(el))
                except Exception as e:
                    pass
            if len(out) == 5:
                # 5 values. All ok
                return out
            else:
                # Must be a comment line
                return []

    def readValues(self, LINES) -> tuple[list, list]:
        """Reads values from EQDSK lines.
        """
        # Part one consists of values before the section of limiter and
        # boundary geometry. Basically meaning reading values until finding
        # two integer values.
        part_1 = []
        # Part two are all values regarding plasma and boundary geometry.
        # HOPEFULLY, the lines start with two INTEGERS and not somewhere in
        # between.
        part_2 = []

        SWITCH = 1

        for line in LINES:
            if SWITCH and self.intPatt.match(line):
                SWITCH = 0

            nan_check = self.nanPatt.findall(line)
            if len(nan_check) > 0:
                raise NanValuesInEqdskFile

            result = self.genPatt.findall(line)
            if SWITCH:
                part_1 += [float(el) for el in result]
            else:
                part_2 += [float(el) for el in result]

        return part_1, part_2

    def readOk(self) -> bool:
        return self.successfullRead

    def generateText(self) -> str:
        """Generates an EQDSK format text.
        """
        # return self.eqdskString

        # The problem with EQDSK G format is that is no longer a consistent
        # standard. Therefore the utility eqdsk.py is mainly used for reading
        # the value, maybe validating the format and using the read value
        # to plot data. Therefore when we want to save the contents of the
        # eqdskString, the eqdsk is no longer generated by the EQDSK G format
        # but the contents of the loaded eqdsk file is just printed.

        #---------------------------------------------------------------------
        if not self.successfullRead:
            return self.eqdskString

        self.generatedText = ''

        self.generatedText += self.HEADER + '\n'

        # First line
        DUMMY = 0.0
        self._PerLineCount = 1


        self.createRealLine(self.RDIM, self.ZDIM, self.RCENTR,
                            self.RLEFT, self.ZMID)

        self.createRealLine(self.RMAXIS,
                            self.ZMAXIS,
                            self.SIMAG, self.SIBRY,
                            self.BCENTR)

        self.createRealLine(self.CURRENT, self.SIMAG, DUMMY,
                            self.RMAXIS, DUMMY)

        self.createRealLine(self.ZMAXIS, DUMMY,
                            self.SIBRY, DUMMY, DUMMY)

        self.createRealLine(*self.FPOL)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.PRES)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.FFPRIM)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.PPRIME)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        for psi_width in self.PSIRZ:
            self.createRealLine(*psi_width)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'
        self._PerLineCount = 1

        self.createRealLine(*self.QPSI)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'

        self.createIntLine(self.NBBBS, self.LIMITR)
        self.generatedText += '\n'
        self._PerLineCount = 1

        alternating_array = [None] * 2 * self.NBBBS
        alternating_array[::2] = self.RBBBS
        alternating_array[1::2] = self.ZBBBS
        self.createRealLine(*alternating_array)
        self.generatedText += '\n'

        self._PerLineCount = 1
        alternating_array = [None] * 2 * self.LIMITR
        alternating_array[::2] = self.RLIM
        alternating_array[1::2] = self.ZLIM
        self.createRealLine(*alternating_array)

        if not self.generatedText.endswith('\n'):
            self.generatedText += '\n'

        return self.generatedText

    def createRealLine(self, *args) -> None:
        Format = '%.9E'
        for el in args:
            stringReal = Format % el
            self.generatedText += ' ' * (16 - len(stringReal)) + stringReal
            if self._PerLineCount % 5 == 0:
                self.generatedText += '\n'
            self._PerLineCount += 1

    def createIntLine(self, *args):
        for el in args:
            self.generatedText += '   ' + str(el)

    def mmMoveInR(self, mm):
        """ Radially moves the values of eqdsk for *mm* millimeters.

        All R relevant quantities:
            RLIM
            RBBBS
            RLEFT
            RMAXIS

        Note that these quantities are not necessarily in millimeters!
        """

        if not self.successfullRead:
            return

        orderLim = '%.1E' % max(self.RLIM)
        orderLim = int(orderLim[-2:])

        scale = 1.0
        if orderLim == 0:
            # Meters
            scale = 1e-3
        elif orderLim == 2:
            # Centimeters
            scale = 1e-2
        elif orderLim == 3:
            # Millimeters
            scale = 1
        else:
            # print('Strange units as RLIM!')
            return

        for i in range(len(self.RLIM)):
            self.RLIM[i] += scale * mm


        orderBoundary = '%.1E' % max(self.RBBBS)
        orderBoundary = int(orderBoundary[-2:])

        scale = 1.0

        if orderBoundary == 0:
            scale = 1e-3
        elif orderBoundary == 2:
            scale = 1e-2
        elif orderBoundary == 3:
            scale = 1
        else:
            # print('Strange Units at RBBBS!')
            return

        for i in range(len(self.RBBBS)):
            self.RBBBS[i] += scale * mm

        orderLeftR = '%.1E' % self.RLEFT
        orderLeftR = int(orderLeftR[-2:])

        scale = 1.0

        if orderLeftR == 0:
            scale = 1e-3
        elif orderLeftR == 2:
            scale = 1e-2
        elif orderLeftR == 3:
            scale = 1
        else:
            # print('Strange Units at RLEFT!')
            return

        self.RLEFT += scale * mm

        orderMagAxisR = '%.1E' % self.RLEFT
        orderMagAxisR = int(orderMagAxisR[-2:])

        scale = 1.0

        if orderMagAxisR == 0:
            scale = 1e-3
        elif orderMagAxisR == 2:
            scale = 1e-2
        elif orderMagAxisR == 3:
            scale = 1
        else:
            # print('Strange Units at RMAXIS!')
            return

        self.RMAXIS += scale * mm

    def fluxOffset(self, offset):
        """Offsetting the flux quantities.
        SIMAG += offset
        SIBRY = offset
        PSIRZ += offset
        """
        if not self.successfullRead:
            return

        self.SIMAG += offset
        self.SIBRY = offset

        X = len(self.PSIRZ)
        Y = len(self.PSIRZ[0])

        for i in range(X):
            for j in range(Y):
                self.PSIRZ[i][j] += offset

    def convertTo_m(self, a):
        """Determine the size unit. If the units are in mm, divide the
        values by 1000.

        Arguments:
            a (array): Array of floats holding the R or Z values.
        """
        import numpy as np
        order = int(np.log10(max([abs(el) for el in a])))
        if order == 2:
            # Hopefully centimeters
            return [el / 100.0 for el in a]
        elif order == 3:
            # Hopefully mm
            return [el / 1000.0 for el in a]
        else:
            return a

    def getRBBBS(self):
        """Gets the R points for LCFS.

        """
        return self.convertTo_m(self.RBBBS)

    def getZBBBS(self):
        """Gets the Z points for LCFS.

        Determine the size unit. If the units are in meters, multiply the
        values by 1000.
        """
        return self.convertTo_m(self.ZBBBS)

    def getRLIM(self):
        """Gets the R points for wall silhouette.
        """
        return self.convertTo_m(self.RLIM)

    def getZLIM(self):
        """Gets the Z points for wall silhouette.
        """
        return self.convertTo_m(self.ZLIM)
