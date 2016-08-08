"""
Core functions and classes for BCG search.

redshift_bias_corr function:
    When looking for photoz members of a cluster, this returns a corrected cluster redshift for photoz

contourapl function:
    Create a contour array that can be used by aplpy show_lines() method to overlay ds9 *_fine.con fine contour files

RetrieveW1Cutouts function:
    Retrieve cutouts from CADC W1 database based on template

Catalog class:
    General usage, comma-separated catalog object

Cluster class:
    A class to handle cluster parameters (X-ray, mass, etc) and functions such a sources catalog and others.

Galaxy class:
    Do not exist yet
    TO INCLUDE:
        - Get things from EzGal
            - k (k+e) correction to a given redshift
            - Stellar mass based on some IMF and SFH
        - And more!
            - Interactive GalFit profile fit?
"""


def redshift_bias_corr(red, field):
    """
    When looking for photoz members of a cluster, this returns a corrected cluster redshift for photoz

    red: input redshift of the cluster (confirmed from specz usually)
    field: north or south field, photoz bias is diffferent in both

    returns: a "corrected" redshift for photoz search
    """

    import numpy as np

    data = np.genfromtxt('/Users/seb-air/Dropbox/Astro/PhD/general_catalogs/redshift_bias_over.txt', comments='#',
                         dtype='str', delimiter='\t')

    bins = []
    over = []
    if field in ['north']:
        for i in range(len(data[:, 0])):
            bins.append(float(data[:, 0][i]))
            over.append(float(data[:, 1][i]))
    elif field in ['south']:
        for i in range(len(data[:, 0])):
            bins.append(float(data[:, 2][i]))
            over.append(float(data[:, 3][i]))

    for i in range(len(bins[:-1])):
        if bins[i] < red <= bins[i + 1]:
            int_over = np.interp(red, [bins[i], bins[i+1]], [over[i], over[i+1]])

    return red + int_over


def contourapl(file):
    """
    Create a contour array that can be used by aplpy show_lines() method to overlay ds9 *_fine.con fine contour files.

    file: contour file

    Returns: A numpy array composed of ra-dec arrays for each vertices of the contour. This way they can be used by AplPy
    """
    import numpy as np                                  # Import numpy and fileinput
    import fileinput

    ra = []                                             # Create empty vectors for ra, dec and contour
    dec = []
    con_arr = []
    same = 0                                            # Initialize a swtich. If the contour level changes,
                                                        #   it starts another ra-dec array

    for line in fileinput.input(file):
        line_str = line.split()                         # Loop through the lines, splitting them
        if not line.strip():                            # If the line is empty, change in contour level
            same = 0                                    # Not the same contour = 0
            ra.pop()                                    # Remove the last ra-dec couple because it's the end
            dec.pop()                                       # of a line segment
            con_arr.append(np.array([ra, dec]))         # Append the final level array
            ra = []                                     # Empty ra-dec
            dec = []
        else:
            if same == 1:                               # If the line is on the same contour level
                ra.append(float(line_str[0]))           # Enter ra-dec couple twice. Once for the end of the line
                ra.append(float(line_str[0]))               # the other for the start of the next line
                dec.append(float(line_str[1]))
                dec.append(float(line_str[1]))
            elif same == 0:                             # Only used for the first line
                ra.append(float(line_str[0]))           # The first point doesn't need to be there twice
                dec.append(float(line_str[1]))
                same = 1                                # Switch to same level so the next vertices are
                                                            # written twice
    return con_arr                                      # Return the array of level arrays


def RetrieveW1Cutouts(output, ra, dec, radius, field, filter):
    """
    Retrieve cutouts from CADC W1 database based on template :

                                                                                                              v- cutout RA/DEC
    http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout?uri=ad:CFHTSG/D2.I.fits&cutout=Circle+ICRS+150.570478+2.172356+0.01
                                                                                ^- field + filter                                ^- cutout radius in degrees

    output: Name of output file
    ra: Cutout central RA
    dec: Cutout central DEC
    radius: Cutout radius in degrees
    field: W1 subfield, e.g. W1+3+1
    filter: W1 filter [U,G,R,I,Z]
    """
    import urllib
    import os

    url_base = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout?uri=ad:CFHTSG/' # Build the url parts
    url_field = field.replace('+', '%2B') + '.'
    url_filter = str(filter).upper() + '.fits&'
    ra = '%s' % ra                                                          # Make Ra/Dec/Radius into strings
    dec = '%s' % dec
    radius = '%s' % radius
    cutout_url = 'cutout=Circle+ICRS+' + ra + '+' + dec + '+' + radius
    complete_url = url_base + url_field + url_filter + cutout_url           # Combine url

    if os.path.isfile(output):                                              # If the file exists, append a _2 an continue
        output_base = output[:-5]
        output = output_base + '_2.fits'

    print complete_url                                                      # Print URL
    urllib.urlretrieve(complete_url, output)                                # Retrieve data
    urllib.urlcleanup()                                                     # Clean URL (possibly not necessary)


class catalog:
    """
    General usage, comma-separated catalog object

    Expected format:
        #,header key, value
        #,header key, value
        ...
        col1, col2, col3,...
        data1, data2, data3,...

    WriteCatalog: Saves a catalog with header and column names according to format
    LoadCatalog: Loads a catalog in the expected format into a dictionary
    AppendCatalog: Append a line or column to a catalog object dictionary
    LookupCatalog: Lookup a specific key-value combo and returns the index
    ExtractIndex: Extract a sub-dictionary of a single index source
    """

    def __init__(self):
        """
        Create catalog object
        """
        self.data = []                      # dictionary of data with column keys as names
        self.columns_keys = []              # names of columns in order for dictionary (so the output is in the right order)
        self.header = []                    # dictionary of header parameter
        self.header_keys = []               # names of header parameters in order (so the output is in the right order)

    def WriteCatalog(self, filename, write_columns=False, write_header=False):
        """
        Write catalog

        filename: Name of the output file
        write_columns: Write name of columns in dictionary in order
        write_header: write the header info
        """
        import csv
        import numpy as np

        with open(filename, 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')                              # Comma separated file
            if write_header:                                                            # Write header
                for i in self.header.keys():                                            # Loop through keys
                    csvwriter.writerow(['#'] + ['%s' % i] +['%s' % (self.header[i])])   # Write header
            if write_columns:                                                           # Write columns
                csvwriter.writerow(self.columns_keys)                                   # Columns in order

            if np.shape(self.data[self.columns_keys[0]]) is ():
                row = []  # Create row from dictionary
                for j in self.columns_keys:  # following prescribed order
                    row.append(self.data[j])
                csvwriter.writerow(row)
            else:
                for i in range(len(self.data[self.columns_keys[0]])):                        # Write data
                    row = []                                                                # Create row from dictionary
                    for j in self.columns_keys:                                             #   following prescribed order
                        row.append(self.data[j][i])
                    csvwriter.writerow(row)
        csvfile.close()

    def LoadCatalog(self, filename):
        """
        Load catalog of expected format

        filename: Name of file to be loaded
        """
        import numpy as np
        import csv

        head_key = []                                           # Value of header params
        head_val = []

        with open(filename, 'rb') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',')

            previous_line_header = 0                            # Used to figure out if after header
            previous_line_data = 0
            num_row = 0
            for row in csvreader:
                if row[0] in ['#']:                             # Header line
                    head_key.append(row[1])                     # Append header key in order
                    head_val.append(row[2])                     #   as well as value
                    previous_line_header = 1
                elif previous_line_header == 1:                 # First line after header is column
                    row_temp = []
                    for i in row:
                        row_temp.append(i.strip())
                    self.columns_keys = row_temp                     # Store row in order
                    previous_line_header = 0
                elif previous_line_data == 0:                   # First data line
                    row_temp = []
                    for i in row:
                        row_temp.append(i.strip())
                    dat = row_temp                                   # Start array
                    previous_line_data = 1
                    num_row += 1
                elif previous_line_data == 1:                   # Append to the array
                    row_temp = []
                    for i in row:
                        row_temp.append(i.strip())
                    dat = np.vstack((dat, row_temp))
                    num_row += 1
        csvfile.close()

        if len(head_key) > 0:                                   # If there is a header, make it into a dict
            self.header_keys = head_key
            self.header = dict(zip(head_key, head_val))
        if len(dat) > 0:                                        # If there is data, arrange it into a dict from
            self.data = dict()                                  #   columns and data
            for i in range(len(self.columns_keys)):
                if num_row > 1:
                    self.data[self.columns_keys[i]] = dat[:, i]
                else:
                    self.data[self.columns_keys[i]] = dat[i]

    def AppendCatalog(self, append_line=False, append_column=False, column_name=False):
        """
        Append a line or column to an existing catalog

        append_line: A vector containing the line to be appended, must be same lenght as the number of columns
        append_column: A vector containing the column to be appended, must be the same lenght as the number of lines
                            and must be associated with a column_name
        column_name: A name for the appended column
        """
        import numpy as np

        if append_line:                                                 # If append line, check the lenght is ok
            if len(append_line) == len(self.columns_keys):
                dat = self.data
                for i, j in zip(append_line, self.columns_keys):        # Append to all columns
                    dat[j] = np.append(dat[j], i)
                self.data = dat
            else:                                                       # Print error if lenght is not ok
                print 'Line is not the right lenght. Line is %i and should be %i' % (len(append_line), len(self.columns_keys))
                raise StandardError
        elif append_column:                                             # If append column, check the lenght
            if not column_name:                                         # Raise error if column_name not specified
                print 'Specify column name'
                raise StandardError
            if len(append_column) == len(self.data[self.columns_keys[0]]):
                dat = self.data
                dat[column_name] = append_column                        # Add new column and name to data dict
                self.data = dat
            else:                                                       # Raise error if lenght is not ok
                print 'Column is not the right lenght. Column is %i and should be %i' % (len(append_column), len(self.data[self.columns_keys[0]]))
                raise StandardError
            if column_name:                                             # Add the column name to the column params to keep
                col = self.columns_keys                                 #   in order when saving
                col = np.hstack((col, column_name))
                self.columns_keys = col

    def LookupCatalog(self, key, value):
        """
        Lookup a specific key-value combo and returns the index

        key: Dict key to look into
        value: value of key looked for

        returns: index of the combo
        """
        if key in self.columns_keys:                                    # If key is in the structure
            if value in self.data[key]:                                 # If value is in the key vector
                vector = list(self.data[key])                           # Transform vector in list
                return vector.index(value)                              # Return index from list
            else:
                print 'Value is not in catalog[key]'                    # If value not in catalog[key], raise error
                raise StandardError
        else:
            print 'Key not in catalog'                                  # If key is not in catalog, raise error
            raise StandardError

    def ExtractIndex(self, index):
        """
        Extract a sub-dictionary of a single index source

        index: index of line in dictionary

        returns: a dictionary containing the values for each original keys at [index]
        """
        ret_dic = dict()                                                # Create empty dictionary
        for i in self.columns_keys:                                     # Loop through keys
            ret_dic[i] = self.data[i][index]                            # Add data and keys to empty dic
        return ret_dic                                                  # Return dictionary for specific index


class Cluster:
    """
    A class to handle cluster parameters (X-ray, mass, etc) and functions such a sources catalog and others.

    __init__: Load cluster related infos into the object
    BuildSourceCatalog: Build and save a catalog of all sources within some radius of the cluster centre
    SourceCatalog: Return the sources as a catalog object
    BuildMembersCatalog: Build and save a catalog of member sources within some radius of the cluster centre
    MembersCatalog: Return the members as a catalog object
    FindBCG: Find brightest cluster galaxy in a given filter within a certain radius from the x-ray centroid and saves the information
    BCGInfo: Return BCG info as a catalog object
    MakeImage: Saves an image of the field with sources and contours overlaid
    """

    def __init__(self, tag):
        """
        Initialize various cluster-related parameters
        """
        from astropy.cosmology import WMAP9 as cosmo
        import os

        cluster_info_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c1_clusters.csv'    # Cluster basic parameters file path
        cluster_info = catalog()                                    # Create catalog for global cluster params
        cluster_info.LoadCatalog(filename=cluster_info_path)        # Load global cluster params
        ind = cluster_info.LookupCatalog(key='tag', value=tag)      # Find index of cluster in global params file
        cluster_basics = cluster_info.ExtractIndex(index=ind)       # Load params values relevant to present cluster

        self.tag = tag                                                      # Cluster tag
        self.type = cluster_basics['type']                                  # c1, c2 class
        self.ra = float(cluster_basics['ra'])                               # cluster ra
        self.dec = float(cluster_basics['dec'])                             # cluster dec
        self.status = cluster_basics['status']                              # redshift status: undef, tentative, confirmed, provisional, photometric
        self.quality = float(cluster_basics['quality'])                     # cluster quality
        self.redshift = float(cluster_basics['redshift'])                   # cluster redshift
        self.scale = cosmo.kpc_proper_per_arcmin(self.redshift).value / 60. # kpc per arcsec
        self.cutout = cluster_basics['cutout']                              # available cutout of cluster

        if tag[0] in ['n', 'N']:                                            # defining north of south field
            self.field = 'north'
        elif tag[0] in ['s', 'S']:
            self.field = 'south'

        if os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/contours/' + self.tag + '_fine.con'):         # looking for x-ray contours
            self.contours = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/contours/' + self.tag + '_fine.con'         # defining contours path
        else:
            self.contours = False

        if os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/sources/' + self.tag + '_sources.txt'):         # looking for sources catalog
            self.sources = True
        else:
            self.sources = False

        if os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/members/' + self.tag + '_members.txt'):         # looking for members catalog
            self.members = True
        else:
            self.members = False

        if os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/bcgs/' + self.tag + '_bcg.txt'):         # looking for bcg info
            self.bcg = True
        else:
            self.bcg = False

    def BuildSourcesCatalog(self, rad=1000., overwrite=False):
        """
        Build and save a catalog of all sources within some radius of the cluster centre. In the North,
        this is based on the Terapix T007 CFHTLS-wide photoz catalog.

        rad: radius in kpc at the cluster redshift
        overwrite: overwrite pre-existing catalog
        """
        import fileinput                                                     # File in/out management package
        import numpy as np
        import datetime
        import os

        data_path = '/Volumes/Storage_bay/CFHTLS/photozCFHTLS-W1_270912.out'       # Path to photoz catalog
        if not os.path.isfile(data_path):
            data_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/photozCFHTLS-W1_270912.out' # Alternative path
        output_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/sources/' + self.tag + '_sources.txt' # Output path of sources

        cl_ra = self.ra                         # Cluster RA
        cl_dec = self.dec                       # Cluster DEC
        kpc_deg = self.scale * 3600.            # kpc per degree
        rad_deg = rad / kpc_deg                 # Source radius in degrees

        if os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/sources/' + self.tag + '_sources.txt') and not overwrite:
            print 'Source file already exists'                  # If source file exists and not overwrite, raise error
            raise StandardError
        elif os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/sources/' + self.tag + '_sources.txt') and overwrite:
            print 'Source file already exists'                  # If source file exists and overwrite, print message, continue
            print 'Overwriting...'

        print 'Building sources catalog within %4.0f kpc of %s from %s' % (rad, self.tag, data_path)    # Building catalog message
        print '...'

        cat = catalog()                 # Create catalog object for data
        cat.header = {'run': datetime.date.today(), 'radius': rad}  # Catalog header information
        cat.columns_keys = ['Id', 'ra', 'dec', 'flag', 'StarGal', 'r2', 'photoz', 'zPDF', 'zPDF_l68', 'zPDF_u68',
                            'chi2_zPDF', 'mod', 'ebv', 'NbFilt', 'zMin', 'zl68', 'zu68', 'chi2_best', 'zp_2', 'chi2_2',
                            'mods', 'chis', 'zq', 'chiq', 'modg', 'U', 'G', 'R', 'I', 'Y', 'Z', 'eU', 'eG', 'eR', 'eI',
                            'eY', 'eZ', 'MU', 'MG', 'MR', 'MI', 'MY', 'MZ'] # Keys in order

        for line in fileinput.FileInput(data_path):                             # Line by line scan of data file
            src_ra = float(line.split(' ')[2])                                  # Extract RA and DEC of source
            src_dec = float(line.split(' ')[3])

            if -rad_deg < src_ra - cl_ra < rad_deg and -rad_deg < src_dec - cl_dec < rad_deg:
                dist = np.sqrt((src_ra - cl_ra) ** 2. + (src_dec - cl_dec) ** 2.)   # If source within defined rad
                if dist < rad_deg:
                    row = filter(None, line.split(' '))                             # Split line, removing whitespaces
                    if not cat.data:                                                # If no data yet, create dict
                        cat.data = dict(zip(cat.columns_keys, row[:-1]))
                    else:
                        cat.AppendCatalog(append_line=row[:-1])                     # Else, append
        fileinput.close()                                                           # Close data file
        cat.WriteCatalog(output_path, write_header=True, write_columns=True)        # Write catalog to output path
        print 'Done'

    def SourcesCatalog(self):
        """
        Return the sources as a catalog object
        """
        cat_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/sources/' + self.tag + '_sources.txt'    # Source catalog path
        cat = catalog()                             # Create catalog object
        cat.LoadCatalog(filename=cat_path)          # Load catalog data
        return cat                                  # Return object

    def BuildMembersCatalog(self, rad=1000., reject_stars=False, class_max=0.5, photoz_sigma=1., photoz_bias_corr=True, overwrite=False):
        """
        Build and save a catalog of member sources within some radius of the cluster centre. In the North,
        this is based on the Terapix T007 CFHTLS-wide photoz catalog.

        rad : radius in kpc at the cluster redshift
        reject_stars:
        class_max:
        photoz_sigma:
        overwrite:
        """
        import fileinput  # File in/out management package
        import numpy as np
        import datetime
        import os

        data_path = '/Volumes/Storage_bay/CFHTLS/photozCFHTLS-W1_270912.out'
        if not os.path.isfile(data_path):
            data_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/photozCFHTLS-W1_270912.out'
        output_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/members/' + self.tag + '_members.txt'

        cl_ra = self.ra                 # Cluster RA
        cl_dec = self.dec               # Cluster DEC
        kpc_deg = self.scale * 3600.    # kpc per degree
        rad_deg = rad / kpc_deg         # Radius input in degrees

        if self.field in ['north']:
            del_red = photoz_sigma * 0.04 * (1. + self.redshift) # Correct for sigma larger at higher redshift. CFHTLS T0007. sig_z/(1+z) = 0.04 for i < 22.5
        elif self.field in ['south']:
            del_red = photoz_sigma * 0.06 * (1. + self.redshift) # Benitez 2000, RMS sig_z/(1+z) = 0.06

        if os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/members/' + self.tag + '_members.txt') and not overwrite:
            print 'Member file already exists'          # If not overwrite and file exists, raise error
            raise StandardError
        elif os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/members/' + self.tag + '_members.txt') and overwrite:
            print 'Member file already exists'          # If overwrite and file exists, print message
            print 'Overwriting...'

        print 'Building members catalog within %4.0f kpc of %s from %s' % (rad, self.tag, data_path)
        print '...'                                     # Print information on catalog

        cat = catalog()                                 # Create catalog object
        cat.header = {'run': datetime.date.today(), 'radius': rad, 'reject_star' : reject_stars, 'class_max' : class_max,
                      'photoz_sigma': photoz_sigma, 'photoz_bias_corr' : photoz_bias_corr}  # Create header
        cat.columns_keys = ['Id', 'ra', 'dec', 'flag', 'StarGal', 'r2', 'photoz', 'zPDF', 'zPDF_l68', 'zPDF_u68',
                            'chi2_zPDF', 'mod', 'ebv', 'NbFilt', 'zMin', 'zl68', 'zu68', 'chi2_best', 'zp_2', 'chi2_2',
                            'mods', 'chis', 'zq', 'chiq', 'modg', 'U', 'G', 'R', 'I', 'Y', 'Z', 'eU', 'eG', 'eR', 'eI',
                            'eY', 'eZ', 'MU', 'MG', 'MR', 'MI', 'MY', 'MZ'] # Keys in order

        if not photoz_bias_corr:                        # Define redshift range if not photoz bias correction applied
            red_min_photoz = self.redshift - del_red
            red_max_photoz = self.redshift + del_red
        elif photoz_bias_corr:                          # If photoz correction applied, get effective photoz
            red_min_photoz = redshift_bias_corr(self.redshift, self.field) - del_red
            red_max_photoz = redshift_bias_corr(self.redshift, self.field) + del_red

        for line in fileinput.FileInput(data_path):     # Loop through lines in data file
            src_ra = float(line.split(' ')[2])          # Load source RA, DEC, photoz and star class
            src_dec = float(line.split(' ')[3])
            src_photoz = float(line.split(' ')[7])
            src_class = float(line.split(' ')[5])

            if not reject_stars:                        # If not star_class cutoff
                if -rad_deg < src_ra - cl_ra < rad_deg and -rad_deg < src_dec - cl_dec < rad_deg \
                        and red_min_photoz <= src_photoz <= red_max_photoz:
                    dist = np.sqrt((src_ra - cl_ra) ** 2. + (src_dec - cl_dec) ** 2.)
                    if dist < rad_deg:                  # If source if within photoz range and radius
                        row = filter(None, line.split(' ')) # Split line and remove whitespaces
                        if not cat.data:                # If catalog empty, add keys
                            cat.data = dict(zip(cat.columns_keys, row[:-1]))
                        else:                           # If catlaog exists, append
                            cat.AppendCatalog(append_line=row[:-1])
            elif reject_stars:                          # If reject stars cutoff
                if -rad_deg < src_ra - cl_ra < rad_deg and -rad_deg < src_dec - cl_dec < rad_deg \
                        and red_min_photoz <= src_photoz <= red_max_photoz and src_class < class_max:
                    dist = np.sqrt((src_ra - cl_ra) ** 2. + (src_dec - cl_dec) ** 2.)
                    if dist < rad_deg:                  # If source is in photoz, class_star and radius
                        row = filter(None, line.split(' ')) # Split line and remove whitespaces
                        if not cat.data:                # If catalog empty, add keys
                            cat.data = dict(zip(cat.columns_keys, row[:-1]))
                        else:                           # If catalog exists, append
                            cat.AppendCatalog(append_line=row[:-1])
        fileinput.close()                               # Close data file
        cat.WriteCatalog(output_path, write_header=True, write_columns=True)    # Write catalog to ouput path
        print 'Done'

    def MembersCatalog(self):
        """
        Return the members as a catalog object
        """
        cat_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/members/' + self.tag + '_members.txt' # Path of members catalog
        cat = catalog()                         # Create catalog object
        cat.LoadCatalog(filename=cat_path)      # Load members
        return cat                              # Return catalog object

    def FindBCG(self, max_rad=350., filter='Z', overwrite=False, xbrightest=False):
        """
        Find brightest cluster galaxy in a given filter within a certain radius from the x-ray centroid and
        saves the information

        max_rad:
        filter:
        overwrite:
        xbrightest: Find additional brightest galaxies up to "x"th brightest
        """
        import numpy as np
        import datetime
        import os

        if not xbrightest:
            output_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/bcgs/' + self.tag + '_bcg.txt'    # Output path of BCG data
        elif xbrightest:
            output_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/xbcgs/' + self.tag + '_xbcg.txt'  # Output path of BCG data

        kpc_deg = self.scale * 3600.    # kpc per degree
        rad_deg = max_rad / kpc_deg     # Lookup radius in degrees

        members = self.MembersCatalog().data    # Put members data in variable

        members_dist_ind = []                   # Index of members within BCG selection radius

        for i in range(len(members['ra'])):     # Loop through members
            src_ra = float(members['ra'][i])    # Members RA and DEC
            src_dec = float(members['dec'][i])

            if np.sqrt((src_ra - self.ra)**2. + (src_dec - self.dec)**2.) < rad_deg:
                members_dist_ind.append(i)      # If member i is within raidus, append its index

        mag_vec_str = self.MembersCatalog().data[filter][members_dist_ind]  # Extract magnitude of members within radius from specified filter
        mag_vec = []
        for i in range(len(mag_vec_str)):               # Make BCG mag vector in float format
            mag_vec.append(float(mag_vec_str[i]))
        for i in range(len(mag_vec)):                   # Change any non real magnitude like -99 and such to 1e8
            if mag_vec[i] < -50.:
                mag_vec[i] = 1.e8

        if not xbrightest:      # If only looking for the BCG
            bcg_index = members_dist_ind[np.argmin(mag_vec)]        # Find BCg index, aka brightest galaxy in filter

            if os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/bcgs/' + self.tag + '_bcg.txt') and not overwrite:
                print 'BCG file already exists'     # If not overwrite and files exists, raise error
                raise StandardError
            elif os.path.isfile('/Users/seb-air/Google Drive/Astro/PhD/BCG_find/bcgs/' + self.tag + '_bcg.txt') and overwrite:
                print 'BCG file already exists'     # If overwrite and file exists, print message
                print 'Overwriting...'

            print 'Finding BCG'                     # Print message
            print '...'

            cat = catalog()                         # Create catalog object for BCG data
            cat.header = {'run': datetime.date.today(), 'max_radius': max_rad, 'filter' : filter}   # Define header variables
            cat.columns_keys = ['Id', 'ra', 'dec', 'flag', 'StarGal', 'r2', 'photoz', 'zPDF', 'zPDF_l68', 'zPDF_u68',
                                'chi2_zPDF', 'mod', 'ebv', 'NbFilt', 'zMin', 'zl68', 'zu68', 'chi2_best', 'zp_2', 'chi2_2',
                                 'mods', 'chis', 'zq', 'chiq', 'modg', 'U', 'G', 'R', 'I', 'Y', 'Z', 'eU', 'eG', 'eR', 'eI',
                                    'eY', 'eZ', 'MU', 'MG', 'MR', 'MI', 'MY', 'MZ'] # Keys in order

            bcg_vec = []                                                    # Vector for BCG data
            for i in cat.columns_keys:                                      # Loop through members keys
                bcg_vec.append(self.MembersCatalog().data[i][bcg_index])    # Append BCG data

            cat.data = dict(zip(cat.columns_keys, bcg_vec))                 # Create dictionary with BCG data

            cat.WriteCatalog(output_path, write_header=True, write_columns=True)    # Write BCG catalog
            print 'Done'

        elif xbrightest:            # If looking for many BCGs
            bcg_index = []
            while len(bcg_index) < xbrightest:                          # Find the index of the next BCG up to x
                bcg_index.append(members_dist_ind[np.argmin(mag_vec)])
                mag_vec[np.argmin(mag_vec)] = 1e8

            if os.path.isfile(
                                    '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/xbcgs/' + self.tag + '_xbcg.txt') and not overwrite:
                print 'xBCG file already exists'  # If not overwrite and files exists, raise error
                raise StandardError
            elif os.path.isfile(
                                    '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/xbcgs/' + self.tag + '_xbcg.txt') and overwrite:
                print 'xBCG file already exists'  # If overwrite and file exists, print message
                print 'Overwriting...'

            print 'Finding %i BCGs' % xbrightest  # Print message
            print '...'

            cat = catalog()  # Create catalog object for BCG data
            cat.header = {'run': datetime.date.today(), 'max_radius': max_rad, 'filter': filter}  # Define header variables
            cat.columns_keys = ['Id', 'ra', 'dec', 'flag', 'StarGal', 'r2', 'photoz', 'zPDF', 'zPDF_l68', 'zPDF_u68',
                                'chi2_zPDF', 'mod', 'ebv', 'NbFilt', 'zMin', 'zl68', 'zu68', 'chi2_best', 'zp_2', 'chi2_2',
                                'mods', 'chis', 'zq', 'chiq', 'modg', 'U', 'G', 'R', 'I', 'Y', 'Z', 'eU', 'eG', 'eR', 'eI',
                                'eY', 'eZ', 'MU', 'MG', 'MR', 'MI', 'MY', 'MZ']  # Keys in order

            for ind in bcg_index:               # Create xbcg catalog
                bcg_vec = []  # Vector for BCG data
                for i in cat.columns_keys:  # Loop through members keys
                    bcg_vec.append(self.MembersCatalog().data[i][ind])  # Append BCG data
                if not cat.data:  # If catalog empty, add keys
                    cat.data = dict(zip(cat.columns_keys, bcg_vec))
                else:  # If catlaog exists, append
                    cat.AppendCatalog(append_line=bcg_vec)

            cat.WriteCatalog(output_path, write_header=True, write_columns=True)  # Write xBCG catalog
            print 'Done'


    def BCGInfo(self):
        """
        Return BCG info as a catalog object
        """
        cat_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/bcgs/' + self.tag + '_bcg.txt'   # BCG path
        cat = catalog()                     # Create catalog object for BCG data
        cat.LoadCatalog(filename=cat_path)  # Load BCg data
        return cat                          # Return BCG catalog

    def xBCGsInfo(self):
        """
        Return xBCGs info as a catalog object
        """
        cat_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/xbcgs/' + self.tag + '_xbcg.txt'   # BCG path
        cat = catalog()                     # Create catalog object for BCG data
        cat.LoadCatalog(filename=cat_path)  # Load BCg data
        return cat                          # Return BCG catalog

    def MakeImage(self, out_path=False, filter='r', show_contours=False, show_sources=False, show_members=False, show_BCG=False,
                  show_xBCGs=False, contours_color='blue', sources_color='red', members_color='green', BCG_color='yellow',
                  show_scalebar=False, scalebar_lenght=500.):
        """
        Saves an image of the field with sources and contours overlaid

        out_path: Output path and name of the image
        filter: Filter to use (determine where the code looks up for fits)
        show_contours: Show X-ray contours is they exist
        show_sources: Show sources from catalog
        show_members: Show members from catalog
        show_BCG: Show BCG position
        show_xBCGs: Show x BCGs with numbers
        contours_color: Contours color
        sources_color: Color of sources
        members_color: Color of members
        BCG_color: Color of BCG
        show_scalebar: Show scalebar or not
        scalebar_lenght: Scalebar lenght in kpc
        """
        import aplpy

        if not out_path:            # Default path if output path not specified
            out_path = '/Users/seb-air/Desktop/' + self.tag + '_image.pdf'

        cutout_path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/cutouts/' + filter + '/' + self.tag + '_cutout.fits'  # General path for cutouts

        fig = aplpy.FITSFigure(cutout_path)     # Create aplpy figure
        fig.show_grayscale(invert=True)         # Show cutout in greyscale
        if show_contours and self.contours:     # If contours exists, add contours
            fig.show_lines(contourapl(self.contours), color=contours_color, linewidth=0.5)
        if show_sources:                        # If sources, add circles for each
            fig.show_circles(self.SourcesCatalog().data['ra'].astype(float), self.SourcesCatalog().data['dec'].astype(float), 3./3600., color=sources_color, linewidth=0.5)
        if show_members:                        # If members, add circles foe each
            fig.show_circles(self.MembersCatalog().data['ra'].astype(float), self.MembersCatalog().data['dec'].astype(float), 4./3600., color=members_color, linewidth=1)
        if show_BCG:                            # If BCG data, add cicle
            fig.show_circles([float(self.BCGInfo().data['ra'])], [float(self.BCGInfo().data['dec'])], 5./3600., color=BCG_color, linewidth=1)
        if show_scalebar:                       # If scalebar, add scalebar of defined lenght in kpc
            scalebar_degrees = scalebar_lenght / self.scale / 3600.
            fig.add_scalebar(scalebar_degrees)
            fig.scalebar.set(linestyle='solid', color='black', linewidth=3) # Deifne scalebar style
            fig.scalebar.set_label(" %3.0f kpc" % scalebar_lenght)          # Scalebar legend
            fig.scalebar.set_font(size='large', family='serif')

        if show_xBCGs:
            fig.show_circles(self.xBCGsInfo().data['ra'].astype(float), self.xBCGsInfo().data['dec'].astype(float), 5./3600., color=BCG_color, linewidth=1)
            for i in range(len(self.xBCGsInfo().data['ra'])):
                ra_xbcg = self.xBCGsInfo().data['ra'].astype(float)[i]
                dec_xbcg = self.xBCGsInfo().data['dec'].astype(float)[i]
                fig.add_label(ra_xbcg + 10./3600., dec_xbcg, i+1, )

        fig.save(out_path)                                                  # Save figure
        fig.close()                                                         # Close figure








############################################
################
# Old functions
###############
############################################

# def LoadClusterCatalog():
#     """
#     """
#     import numpy as np
#
#     path = '/Users/seb-air/Google Drive/Astro/PhD/BCG_find/catalog/c2_clusters.csv'
#     print 'Loading cluster catalog from : ', path
#     print '...'
#     print 'Done'
#     return catalog().LoadCatalog(filename=path)


# def GetClusterInfo(tag):
#     """
#     """
#     cluster_cat = LoadClusterCatalog()
#
#     tag_vector = list(cluster_cat['tag'])
#
#     if tag in tag_vector:
#         return cluster_cat[tag_vector.index(tag)]
#     else:
#         print 'Cluster tag not in the catalog'
#         raise LookupError