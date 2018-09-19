import numpy, struct, os

"""
    Not written by Greg Blumberg (OU/CIMMS)

    I got this file from one of Dave Turner's colleagues (Aronne Merrelli)
"""

class panel_file:
    """
    Object class defining LBLRTM fortran formatted binary file.

    Since there is no documentation for these files, the assumed format 
    was cobbled together mainly from previous MATLAB codes (lbl_read.m 
    from various SSEC people), and then some reverse engineering and/or 
    guessing for other parts. The header is still not correctly read in, 
    but I think otherwise this is correct.

    Usage:
    panel_data = panel_file('FileName')

    Will load just the header; this does the various format checking which 
    is necessary to read the file, but does not attempt to load the arrays.
    This might be preferable when in doubt of the actual format of a file.

    To also load data:
    panel_data = panel_file('FileName', do_load_data = True)

    printing panel data to console gives a summary of the file contents:

    In [12]: panel_data  = lblrtm_utils.panel_file('TAPE10')
    In [13]: print panel_data
    -------> print(panel_data)
    filename    : TAPE10
    panel_format: single
    n_panels    : 370
    n_pts       : 887497
    vmin        : 500.000000
    vmax        : 700.000000
    dv          : 0.000225353

    The header can be printed in a similar fashion:

    In [14]: print panel_data.hdr
    -------> print(panel_data.hdr)
    user_id     :   blah blah blah
    secant      : 0.999611
    p_avg       : 10.622692
    t_avg       : 234.389465
    molecule id, col_dens:
           H2O    1.07833e+17
           CO2    8.65037e+18
    etc...

    The different header fields are object attributes:
    In [16]: print panel_data.hdr.secant
    -------> print(panel_data.hdr.secant)
    0.999610960484

    If the data is loaded, then object attributes v, data1 are defined 
    for a single panel file, and in addition, data2 for a double panel 
    file. These attributes will be 1-D numpy ndarrays:

    In [18]: panel_data.load_data()
    In [19]: type(panel_data.v), panel_data.v.shape
    Out[19]: (<type 'numpy.ndarray'>, (887497,))


    """

    def __init__(self, filename, do_load_data = False, 
                 valid_avg_P = [1e-4,1e4], valid_avg_T = [100,400]):
        
        self.filename = filename
        self._valid_avg_P = valid_avg_P
        self._valid_avg_T = valid_avg_T

        [junk_type, float_type, int_type] = self._check_field_type()

        self._data_types = {'junk':junk_type, 'int':int_type, 
                           'float':float_type}

        self._f = open(self.filename, 'rb')
        try:
            self._fread_hdr()
        except IOError as theError:
            self._f.close()
            raise theError

        self._hdr_size = self._f.tell()
        self._f.close()

        try:
            self._check_panel_format()
        except IOError as theError:
            self._f.close()
            raise theError

        self._check_length()

        if do_load_data:
            self.load_data()

    def load_data(self):

        self._f = open(self.filename, 'rb')
        self._f.seek(self._hdr_size)

        # allocate numpy arrays for efficiency
        panel_dtype = numpy.dtype( [ ('v1', 'float64'), 
                                     ('v2', 'float64'), 
                                     ('dv', 'float64'), 
                                     ('n_pts', 'int') ] )
        self.panel_hdrs = numpy.zeros(self.n_panels, panel_dtype)
        self.data1 = numpy.zeros(self.n_pts, 'float32')
        if self.panel_format == 'double':
            self.data2 = numpy.zeros(self.n_pts, 'float32')

        ct = 0
        for n in range(self.n_panels):
            panel_hdr = self._fread_panel_hdr()
            self.panel_hdrs['v1'][n] = panel_hdr[0]
            self.panel_hdrs['v2'][n] = panel_hdr[1]
            self.panel_hdrs['dv'][n] = panel_hdr[2]
            self.panel_hdrs['n_pts'][n] = panel_hdr[3]
            panel_data = self._fread_panel(self.panel_format, 
                                           panel_hdr[3], 'read')
            self.data1[ct:ct+panel_hdr[3]] = panel_data[0]
            if self.panel_format == 'double':
                self.data2[ct:ct+panel_hdr[3]] = panel_data[1]
            ct = ct + panel_hdr[3]

        self.vmin = self.panel_hdrs['v1'][0]
        self.vmax = self.panel_hdrs['v2'][-1]
        self.dv = self.panel_hdrs['dv'][0]
        self.v = numpy.linspace(self.vmin, self.vmax, self.n_pts)


    def _check_field_type(self):
        """
        Attempts to check the panel file's type, either 32 bit v 64 bit, 
        and the size of the junk fields. Needed since the LBLRTM panel 
        file does not self-describe in any way. This is done by a brute 
        force method where the first several fields are attempted to be 
        read in a number of ways, and the first way that returns believable 
        values is taken to be the correct file type.
        This is a private function, and should not typically be accessed 
        by users.
        """

        # the 4 tested possibilities are:
        #
        # file_type 1: Junk fields are 4 bytes, default floats are 32 bit.
        # file_type 2: Junk fields are 4 bytes, default floats are 64 bit.
        # file_type 3: Junk fields are 8 bytes, default floats are 32 bit.
        # file_type 4: Junk fields are 8 bytes, default floats are 64 bit.
        #
        # It is possible that file types 3 & 4 are solely due to gfortran 
        # compiler problems on 64 bit linux. From my limited experience 
        # attempting to compile LBLRTM on a few machines:
        # Using the Intel ifort compiler produces files of type 1 or 2 only. 
        # Using gfortran on 32 bit Mac OS X produces file type 1 
        # (have been unable to compile it with default double precision).
        # Using f77 on 32 bit linux produces file type 1.
        # Using gfortran on 64 bit linux produces file type 3 and 4.

        self._f = open(self.filename, 'rb')
        # reads max number of bytes we will possibly need.
        test_bytes = self._f.read(112)
        self._f.close()

        # now check each case - note that I hardcoded the length of the 
        # unpacked data in each case; this was simpler to implement, 
        # since I only have 4 cases.
        # the tests are done on test_data[3] and test_data[4] which are 
        # supposed to be the average pressure and temperature in the 
        # profile data.

        # Case 1
        junk_type = 'i'
        float_type = 'f'
        int_type = 'i'
        if self._valid_header_bytes(test_bytes, junk_type, float_type):
            return [junk_type, float_type, int_type]
        # Case 2
        junk_type = 'i'
        float_type = 'd'
        int_type = 'q'
        if self._valid_header_bytes(test_bytes, junk_type, float_type):
            return [junk_type, float_type, int_type]
        # Case 3
        junk_type = 'q'
        float_type = 'f'
        int_type = 'i'
        if self._valid_header_bytes(test_bytes, junk_type, float_type):
            return [junk_type, float_type, int_type]
        # Case 4
        junk_type = 'q'
        float_type = 'd'
        int_type = 'q'
        if self._valid_header_bytes(test_bytes, junk_type, float_type):
            return [junk_type, float_type, int_type]

        # if we make it here, no tests worked - so raise an error
        raise IOError("Failed to determine field size in LBLRTM Panel file")


    def _check_panel_format(self):
        # Brute force approach to attempt to read double panel data to see 
        # we get something sensible. If nonsense, then assume it is a single 
        # panel, and try that - if neither works, throw an error.

        self._f = open(self.filename, 'rb')
        self._f.seek(self._hdr_size)

        panel_hdr = self._fread_panel_hdr()
        if not self._valid_panel_hdr(panel_hdr):
            raise IOError('Panel header not successfully read')
        # this is a little clumsy - tell size of header by how much file 
        # pointer moves after panel header read
        self._panel_hdr_size = self._f.tell() - self._hdr_size
        

        # now left with 2 ugly cases
        # First, there is only one panel in the entire file. 
        # In this case, check the size of the file versus the size of the 
        # header plus a single or double format panel.

        file_size = os.stat(self.filename).st_size
        single_panel_size = self._fread_panel('single',panel_hdr[3],'no')
        double_panel_size = self._fread_panel('double',panel_hdr[3],'no')
        if file_size == (self._hdr_size + single_panel_size + 
                         self._panel_hdr_size):
            self.panel_format = 'single'
            self._f.close()
            return
        if file_size == (self._hdr_size + double_panel_size + 
                         self._panel_hdr_size):
            self.panel_format = 'double'
            self._f.close()
            return

        # files seem to have a "null" panel header at the end, sometimes?
        # so check for that.
        if file_size == (self._hdr_size + single_panel_size + 
                         2*self._panel_hdr_size):
            self.panel_format = 'single'
            self._f.close()
            return
        if file_size == (self._hdr_size + double_panel_size + 
                         2*self._panel_hdr_size):
            self.panel_format = 'double'
            self._f.close()
            return

        # Second case - there are multiple panels. In this case, try to read 
        # the second panel's header to see if we can get sensible values.
        self._f.seek(self._hdr_size)
        panel_hdr = self._fread_panel_hdr()
        self._f.seek(single_panel_size, 1)
        panel_hdr = self._fread_panel_hdr()
        single_panel_valid = self._valid_panel_hdr(panel_hdr)

        # this can cause a read-over-EOF error; when this occurred, the 
        # above method did succeed (I think this produced an error only 
        # when the file was between 1 and 2 full single panels in length)
        if not single_panel_valid:
            self._f.seek(self._hdr_size)
            panel_hdr = self._fread_panel_hdr()
            self._f.seek(double_panel_size, 1)
            panel_hdr = self._fread_panel_hdr()
            double_panel_valid = self._valid_panel_hdr(panel_hdr)
        else:
            double_panel_valid = False

        if single_panel_valid:
            self.panel_format = 'single'
            self._f.close()
            return
        elif double_panel_valid:
            self.panel_format = 'double'
            self._f.close()
            return
        else:
            self._f.close()
            raise IOError("Could not determine panel format")


    def _check_length(self):
        # compute number of panels, by sequentially reading the file.
        # Unfortunately, I do not think there is anyway to check this 
        # except for reading through the file. Previous versions attempted 
        # to take a shortcut by assuming the panel structure was always 
        # n L-sized panels (L is usually 2400, but this is not assumed), 
        # with one panel at the end with the remainder;
        # however, since that time I discovered that at some points it 
        # will write the remainder-sized panel at the *front* of the file, 
        # which of course screws everything up.

        # Now, the method is to just repeat reading panels until finding 
        # a panel with -99 points (sometimes seen as an "EOF-like" marker; 
        # reading to EOF; or reading a panel with 0 points, or reaching EOF.

        self._f = open(self.filename, 'rb')
        self._f.seek(self._hdr_size)

        # read panels until "end condition" is met: either nonpositive 
        # number of points in panel, or reached EOF. This code 
        # should even catch a case of no panels (if the initial n_pts 
        # count is already nonpositive)

        file_size = os.stat(self.filename).st_size
        n_panels = 0
        total_n_pts = 0
        
        at_EOF = self._f.tell() == file_size
        at_last_panel = False

        while not at_last_panel and not at_EOF:

            panel_hdr = self._fread_panel_hdr()
            n_pts_in_panel = panel_hdr[3]
            if n_pts_in_panel > 0:
                total_n_pts += n_pts_in_panel
                n_panels += 1
                # seek argument means this should just do a seek, and 
                # not read data, so should be hopefully fast.
                panel_size = self._fread_panel(self.panel_format, 
                                               n_pts_in_panel, 'seek')
            else:
                at_last_panel = True

            # check if at EOF
            at_EOF = self._f.tell() == file_size


        self.n_panels = n_panels
        self.n_pts = total_n_pts

        self._f.close()

    def _valid_header_bytes(self, test_bytes, 
                            test_junk_type, test_float_type):
        unpack_fmt = '=' + test_junk_type + \
            '80sd' + 2*test_float_type
        test_data = struct.unpack(unpack_fmt, 
                                  test_bytes[0:struct.calcsize(unpack_fmt)])
        return \
            (test_data[3] > self._valid_avg_P[0]) & \
            (test_data[3] < self._valid_avg_P[1]) & \
            (test_data[4] > self._valid_avg_T[0]) & \
            (test_data[4] < self._valid_avg_T[1])
        
    def _valid_panel_hdr(self, panel_hdr):
        # these checks are fully hardcoded - the wavenumber values (first and 
        # second values) must be between 1 and 50000, the wavenum increment 
        # (dv, third value) must be between 1e-9 and 10.0, and the number of 
        # points (np), must be between 0 and 10000. The last value seems 
        # to be 2400 in practice, at most ... so just guessing that 10000
        # makes sense as an upper limit.
        # also, check that the wnum's are in the right order (should be wnum 
        # min followed my wnum max, so wnum2 > wnum1.)
        return \
            (panel_hdr[0] > 1) & (panel_hdr[0] < 50000) & \
            (panel_hdr[1] > 1) & (panel_hdr[1] < 50000) & \
            (panel_hdr[2] > 1e-9) & (panel_hdr[2] <= 10.0) & \
            (panel_hdr[3] > 0) & (panel_hdr[3] < 10000) & \
            (panel_hdr[1] > panel_hdr[0])

    def _fread_hdr(self):
        self.hdr = panel_file_hdr(self._f, self._data_types)

    def _fread_panel_hdr(self):
        # unpack string data into binary, and discard the first 
        # and last values since those are junk values.
        unpack_fmt = '=' + self._data_types['junk'] + \
            'dd' + self._data_types['float'] + self._data_types['int'] + \
            self._data_types['junk']
        hdr_size = struct.calcsize(unpack_fmt)
        data = struct.unpack(unpack_fmt, self._f.read(hdr_size))
        return data[1:-1]

    def _fread_panel(self, panel_type, panel_len, read_file):
        if panel_type == 'single':
            panel_fmt = \
                self._data_types['junk'] + \
                panel_len*self._data_types['float'] + \
                self._data_types['junk']
        else:
            panel_fmt = \
                self._data_types['junk'] + \
                panel_len*self._data_types['float'] + \
                self._data_types['junk'] + \
                self._data_types['junk'] + \
                panel_len*self._data_types['float'] + \
                self._data_types['junk']
        panel_size = struct.calcsize(panel_fmt)
        if read_file == 'no':
            return panel_size
        elif read_file == 'seek':
            self._f.seek(panel_size, 1)
            return panel_size
        else:
            raw_bytes = self._f.read(panel_size)
            data = struct.unpack(panel_fmt, raw_bytes)
            # get rid of junk bytes
            if panel_type == 'single':
                data = [data[1:panel_len+1]]
            else:
                data = [data[1:panel_len+1], data[panel_len+3:2*panel_len+3]]
            return data        

    def __repr__(self):
        return \
            '{0:12s}: {1:s}\n'.format('filename', self.filename) + \
            '{0:12s}: {1:s}\n'.format('panel_format', self.panel_format) + \
            '{0:12s}: {1:d}\n'.format('n_panels', self.n_panels) + \
            '{0:12s}: {1:d}\n'.format('n_pts', self.n_pts) + \
            '{0:12s}: {1:f}\n'.format('vmin', self.hdr.v1) + \
            '{0:12s}: {1:f}\n'.format('vmax', self.hdr.v2) + \
            '{0:12s}: {1:g}\n'.format('dv', self.hdr.dv)


class panel_file_hdr:

    def __init__(self, _f, data_types):

        unpack_fmt = '=' + data_types['junk'] + '80sd' + \
            2*data_types['float']
        raw_bytes = _f.read(struct.calcsize(unpack_fmt))
        data = struct.unpack(unpack_fmt, raw_bytes)
        self.user_id = data[1]
        self.secant = data[2]
        self.p_avg = data[3]
        self.t_avg = data[4]

        unpack_fmt = 64*'8s'
        raw_bytes = _f.read(struct.calcsize(unpack_fmt))
        self.molecule_id = struct.unpack(unpack_fmt, raw_bytes)

        unpack_fmt = 64*data_types['float']
        raw_bytes = _f.read(struct.calcsize(unpack_fmt))
        self.mol_col_dens = struct.unpack(unpack_fmt, raw_bytes)

        unpack_fmt = 2*data_types['float'] + \
            'dd' + 2*data_types['float'] + \
            11*data_types['int'] + 2*data_types['float'] + \
            6*data_types['int'] + data_types['float']
        raw_bytes = _f.read(struct.calcsize(unpack_fmt))
        data = struct.unpack(unpack_fmt, raw_bytes)
        self.broad_dens = data[0]
        self.dv = data[1]
        self.v1 = data[2]
        self.v2 = data[3]
        self.t_bound = data[4]
        self.emis_bound = data[5]
        self.lblrtm_flag = dict()
        self.lblrtm_flag['hirac'] = data[6]
        self.lblrtm_flag['lblf4'] = data[7]
        self.lblrtm_flag['xscnt'] = data[8]
        self.lblrtm_flag['aersl'] = data[9]
        self.lblrtm_flag['emit'] = data[10]
        self.lblrtm_flag['scan'] = data[11]
        self.lblrtm_flag['plot'] = data[12]
        self.lblrtm_flag['path'] = data[13]
        self.lblrtm_flag['jrad'] = data[14]
        self.lblrtm_flag['test'] = data[15]
        self.lblrtm_flag['merge'] = data[16]
        self.lblrtm_flag['scnid'] = data[17]
        self.lblrtm_flag['hwhm'] = data[18]
        self.lblrtm_flag['idabs'] = data[19]
        self.lblrtm_flag['atm'] = data[20]
        self.lblrtm_flag['layr1'] = data[21]
        self.lblrtm_flag['nlayr'] = data[22]
        self.n_mol = data[23]
        self.layer = data[24]
        self.yi1 = data[25]

        unpack_fmt = '=8s8s6s8s4s46s' + data_types['junk']
        raw_bytes = _f.read(struct.calcsize(unpack_fmt))
        data = struct.unpack(unpack_fmt, raw_bytes)
        # drop last junk word
        self.yid = data[0:-1]

    def __repr__(self):

        molecule_id_col_str = ''
        for n in range(len(self.mol_col_dens)):
            molecule_id_col_str = molecule_id_col_str + \
                '    {0:9s} {1:g}\n'.format(self.molecule_id[n], 
                                           self.mol_col_dens[n])

        lblrtm_flag_str = 'lblrtm_flag:\n'
        for item in self.lblrtm_flag.items():
            lblrtm_flag_str = lblrtm_flag_str + \
                '    {0:12s}: {1:g}\n'.format(item[0], item[1])

        yid_fmt_str = '{0:12s}:\n{1:>12s}\n{2:>12s}\n{3:>10s}' + \
            '\n{4:>12s}\n{5:>8s}\n{5:>50s}\n'

        return \
            '{0:12s}: {1:s}\n'.format('user_id', self.user_id) + \
            '{0:12s}: {1:f}\n'.format('secant', self.secant) + \
            '{0:12s}: {1:f}\n'.format('p_avg', self.p_avg) + \
            '{0:12s}: {1:f}\n'.format('t_avg', self.t_avg) + \
            'molecule id, col_dens:\n' + \
            molecule_id_col_str + \
            '{0:12s}: {1:g}\n'.format('broad_dens', self.broad_dens) + \
            '{0:12s}: {1:g}\n'.format('dv', self.dv) + \
            '{0:12s}: {1:f}\n'.format('v1', self.v1) + \
            '{0:12s}: {1:f}\n'.format('v2', self.v2) + \
            '{0:12s}: {1:f}\n'.format('t_bound', self.t_bound) + \
            '{0:12s}: {1:f}\n'.format('emis_bound', self.emis_bound) + \
            lblrtm_flag_str + \
            '{0:12s}: {1:g}\n'.format('n_mol', self.n_mol) + \
            '{0:12s}: {1:g}\n'.format('yi1', self.yi1) + \
            yid_fmt_str.format('yid', self.yid[0], self.yid[1], self.yid[2], 
                               self.yid[3], self.yid[4], self.yid[5])

