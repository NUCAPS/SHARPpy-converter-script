import xarray as xr
import numpy as np
import os, sys, glob, datetime
from csv import writer

# Define constants
Rd = 287.
g = 9.81
p_std = 1013.25
navog = 6.02214199e+23  # NIST Avogadro's number (molecules/mole)
mw_d = 28.9644  # gm/mole dry air
mw_wv = 18.0151  # gm/mole water
g_std = 980.664  # acceleration of gravity, cm/s^2
cdair = 1000.0 * p_std * navog / (mw_d * g_std)
wv_eps = mw_wv / mw_d
FILL_VAL = np.nan

# Define the home directory where the .csv and text files will be stored.
HOME = os.path.expanduser("~")

#########################################
#### Cleanup Old Text Files and CSVs ####
#########################################

def create_text_file_path():
    # Create the satellite directory path to store the text files.
    if not os.path.exists(os.path.join(HOME, '.sharppy', 'datasources', satName)):
        os.makedirs(os.path.join(HOME, '.sharppy', 'datasources', satName))

def remove_old_txt_csv():
    # Remove text files if they already exist
    textFiles = glob.glob(os.path.join(HOME, '.sharppy', 'datasources', satName, f'*{satName}*.txt'))
    for f in textFiles:
        os.remove(f)

    # Remove CSV if already exists
    pathCSV = os.path.join(HOME, '.sharppy', 'datasources', f'{satName}_case_study.csv')
    if os.path.isfile(pathCSV):
        os.remove(pathCSV)

def write_csv_header():
    # Store CSV header in a temporary list
    csvHeader = ['icao','iata','synop','name','state','country','lat','lon','elev','priority','srcid']
    write_obj = open(os.path.join(HOME, '.sharppy', 'datasources', f'{satName}_case_study.csv'), "w", newline='')
    csv_writer = writer(write_obj)
    csv_writer.writerow(csvHeader)
    write_obj.close()


################################
#### Science Code Functions ####
################################

# Find the surface within a sounding
def findSurface(pres, surfpres):
    diff = np.abs((pres - surfpres))
    mindiff = np.min(diff)
    clev = np.where(diff == mindiff)[0][0]
    surflev = clev
    if surfpres < pres[clev]:
        surflev = clev
    if surfpres > pres[clev] and mindiff >= 5.0:
        surflev = clev
    if surfpres > pres[clev] and mindiff > 5.0:
        surflev = clev + 1
    return surflev

# Find BOTLEV & BLMULT
def get_botlev_blmult(plev, psurf, nobs):
    botlev = np.zeros((nobs), dtype=float)
    blmult = np.zeros((nobs), dtype=float)

    for i in range(nobs):
        surflev = findSurface(plev, psurf[i])
        num = psurf[i] - plev[surflev - 1]
        denom = plev[surflev] - plev[surflev - 1]
        blmult[i] = num / denom
        botlev[i] = surflev
    return blmult, botlev

# Calculate surface water vapor column density using BLMULT
def calc_wvcd_sfc(botlev, blmult, nobs, p_layer, wvcd):
    nlev_100 = np.shape(p_layer)[0]
    wvcd_sfc = np.zeros((nobs), dtype=float)

    for i in range(nobs):
        if botlev[i] < nlev_100 - 1:
            sfc = int(botlev[i])
        if botlev[i] == nlev_100 - 1:
            sfc = nlev_100 - 1

        wvcd_sfc[i] = wvcd[i, int(sfc)] * blmult[i]
    return wvcd_sfc

# Calculate surface temperature using BLMULT
def calc_Tsfc(botlev, blmult, nobs, plev, temperature):
    nlev_100 = np.shape(plev)[0]
    tsfc = np.zeros((nobs), dtype=float)

    for i in range(nobs):
        if botlev[i] < nlev_100 - 1:
            sfc = int(botlev[i])
        if botlev[i] == nlev_100 - 1:
            sfc = nlev_100 - 1

        t_diff = temperature[i, int(sfc)] - temperature[i, int(sfc) - 1]
        tsfc[i] = temperature[i, int(sfc) - 1] + blmult[i] * t_diff
    return tsfc

# Insert pressure into final lists.
def insert_surface_pressure(botlev, nobs, plev, psurf):
    nlev_100 = np.shape(plev)[0]
    blmult_P_ALL = []

    for i in range(nobs):
        sfc = botlev[i]
        # Add surface values
        blmult_P_footprint = np.append(plev[0:int(sfc)], psurf[i])

        # Append each footprint to the larger array
        blmult_P_ALL.append(blmult_P_footprint)

    # Convert to arrays
    blmult_P_ALL = np.asarray(blmult_P_ALL, dtype="object")

    return blmult_P_ALL

# Insert surface temperature into final lists.
def insert_surface_temperature(botlev, nobs, plev, tsfc, temperature):
    nlev_100 = np.shape(plev)[0]
    blmult_T_ALL = []

    for i in range(nobs):
        sfc = botlev[i]
        if sfc + 1 == nlev_100 - 1:
            temperature[i, sfc + 1] = FILL_VAL
        if sfc + 1 > nlev_100 - 1:
            temperature[i, sfc + 1:nlev_100 - 1] = FILL_VAL

        lev_T = np.zeros(nlev_100, dtype=float)

        for j in range(int(sfc)):
            lev_T[j] = temperature[i, j]

        # Add surface values
        blmult_T_footprint = np.append(lev_T[0:int(sfc)], tsfc[i])

        # Append each footprint to the larger array
        blmult_T_ALL.append(blmult_T_footprint)

    # Convert to arrays
    blmult_T_ALL = np.asarray(blmult_T_ALL, dtype="object")

    return blmult_T_ALL

# Insert surface values into final lists.
def insert_surface_water_vapor(botlev, nobs, plev, wvcd_sfc, wvcd):
    nlev_100 = np.shape(plev)[0]
    blmult_wvcd_ALL = []

    for i in range(nobs):
        sfc = botlev[i]
        if sfc + 1 == nlev_100 - 1:
            wvcd[i, sfc + 1] = FILL_VAL
        if sfc + 1 > nlev_100 - 1:
            wvcd[i, sfc + 1:nlev_100 - 1] = FILL_VAL

        lev_wvcd = np.zeros(nlev_100, dtype=float)
        lev_wvcd[0] = wvcd[i, 0]

        for j in range(1, int(sfc)):
            lev_wvcd[j] = (wvcd[i,j-1] + wvcd[i,j]) / 2

        # Add surface values
        blmult_wvcd_footprint = np.append(lev_wvcd[0:int(sfc)], wvcd_sfc[i])

        # Append each footprint to the larger array
        blmult_wvcd_ALL.append(blmult_wvcd_footprint)

    # Convert to arrays
    blmult_wvcd_ALL = np.asarray(blmult_wvcd_ALL, dtype="object")

    return blmult_wvcd_ALL

# Convert water vapor column density to mixing ratio
def convert_cd2mr(nobs, blmult_wvcd_ALL, blmult_P_ALL, psurf, botlev):
    wvmr = []

    for i in range(nobs):
        nlev_NEW = np.shape(blmult_P_ALL[i])[0]
        wvmr_footprint = np.zeros((nlev_NEW), dtype=float)
        wvcd_footprint = blmult_wvcd_ALL[i]

        deltap = np.zeros((nlev_NEW), dtype=float)
        pres = blmult_P_ALL[i]
        deltap[0] = pres[0]
        deltap[1:nlev_NEW] = pres[1:nlev_NEW] - pres[0:nlev_NEW-1]

        for j in range(nlev_NEW):
            wvmr_footprint[j] = wvcd_footprint[j] / ((cdair * deltap[j] / p_std) / wv_eps)

        # Append footprint array into larger array
        wvmr.append(wvmr_footprint)

    # Convert larger array to numpy array
    wvmr = np.asarray(wvmr, dtype="object")
    return wvmr

# Calculate virtual temperature
def calc_virtual_temperature(nobs, blmult_P_ALL, wvmr, blmult_T_ALL):
    tv = []

    for i in range(nobs):
        nlev_NEW = np.shape(blmult_P_ALL[i])[0]
        tv_footprint = np.zeros((nlev_NEW), dtype=float)
        T_footprint = blmult_T_ALL[i]
        wvmr_footprint = wvmr[i]

        for j in range(nlev_NEW):
            tv_footprint[j] = T_footprint[j] * (1 + 0.608 * wvmr_footprint[j])

        # Append footprint array into larger array
        tv.append(tv_footprint)

    # Convert larger array to numpy array
    tv = np.asarray(tv, dtype="object")
    return tv

# Calculate mean sea level pressure
def calc_mslp(nobs, topography, tsfc, psurf):
    mslp = np.zeros((nobs), dtype=float)

    height_const = np.multiply(0.0065, topography)
    fraction = np.divide(height_const, np.add(tsfc, height_const))
    mslp = np.multiply(psurf, np.power((1 - fraction), -5.257))
    return mslp

# Calculate geopotential Height
def calc_geopotential_height(nobs, blmult_P_ALL, mslp, tv):
    z = []

    for i in range(nobs):
        nlev_NEW = np.shape(blmult_P_ALL[i])[0]
        z_footprint = np.zeros((nlev_NEW), dtype=float)
        mtv = np.zeros((nlev_NEW), dtype=float)
        tv_footprint = tv[i]
        tvsfc = tv_footprint[nlev_NEW - 1]
        P_footprint = blmult_P_ALL[i]

        for j in range(nlev_NEW):
            mtv[j] = (tvsfc + tv_footprint[j]) / 2
            z_footprint[j] = ((Rd * mtv[j]) / g) * np.log(mslp[i] / P_footprint[j])

        # Append footprint array into larger array
        z.append(z_footprint)

    # Convert larger array to numpy array
    z = np.asarray(z, dtype="object")
    return z

# Convert mixing ratio to dew point
def calc_dewpoint(nobs, blmult_P_ALL, blmult_T_ALL, wvmr):
    dew_point = []

    for i in range(nobs):
        nlev_NEW = np.shape(blmult_P_ALL[i])[0]
        RH_footprint = np.zeros((nlev_NEW), dtype=float)
        dew_point_footprint = np.zeros((nlev_NEW), dtype=float)
        T_footprint = blmult_T_ALL[i]
        P_footprint = blmult_P_ALL[i]
        wvmr_footprint = wvmr[i]

        for j in range(nlev_NEW):
            # Substitute mixing ratio (MR) for specific humidity (SHx) since they are approximately equal.
            RH_footprint[j] = find_RH(P_footprint[j], T_footprint[j], wvmr_footprint[j])
            dew_point_footprint[j] = find_dewpoint(T_footprint[j], RH_footprint[j])

        # Append footprint array into larger array
        dew_point.append(dew_point_footprint)

    # Convert larger array to numpy array
    dew_point = np.asarray(dew_point, dtype="object")
    return dew_point

# Find relative humidity for single profile.
def find_RH(P, TEMP_K, SHx):
    "Calculate Relative Humidity(0 to 100) from Pressure(hPa), Temp(K), and Specific Humidity."
    a = 22.05565
    b = 0.0091379024
    c = 6106.396
    epsilonx1k = 622.0

    shxDenom = SHx * 0.378
    shxDenom += epsilonx1k

    tDenom = -b*TEMP_K
    tDenom += a
    tDenom -= c/TEMP_K

    RH = P * SHx
    RH /= shxDenom
    RH /= np.exp(tDenom)

    RH = RH * 1000
    return RH

# Find dewpoint for single profile.
def find_dewpoint(TEMP_K, RH):
   "Calculate dewpoint(C) from temperature(K) and relative humidity(0 to 100)."
   b=0.0091379024*TEMP_K
   b += 6106.396/TEMP_K
   b -= np.log(RH/100)
   val = b*b
   val -= 223.1986
   val = np.sqrt(val)
   DpT = b-val
   DpT /= 0.0182758048
   DpT = DpT - 273.15
   return DpT

# Find cloud top fraction and pressure for each layer
def find_ctf_ctp(nobs, cloud_top_fraction, cloud_top_pressure):
    ctf_low = []
    ctf_high = []
    ctp_low = []
    ctp_high = []

    for i in range(nobs):
        ctf_low_footprint = np.zeros((2), dtype=float)
        ctf_high_footprint = np.zeros((2), dtype=float)
        ctp_low_footprint = np.zeros((2), dtype=float)
        ctp_high_footprint = np.zeros((2), dtype=float)

        # Populate the footprint arrays
        ctf_low_footprint = cloud_top_fraction[i, 1]
        ctf_high_footprint = cloud_top_fraction[i, 0]
        ctp_low_footprint = cloud_top_pressure[i, 1]
        ctp_high_footprint = cloud_top_pressure[i, 0]

        # Multiply each element in ctf_low_footprint and ctf_high_footprint by 100 to get %
        ctf_low_footprint = ctf_low_footprint * 100
        ctf_high_footprint = ctf_high_footprint * 100

        # Round to nearest integer
        ctf_low_footprint = int(np.round(ctf_low_footprint, decimals=0))
        ctf_high_footprint = int(np.round(ctf_high_footprint, decimals=0))
        ctp_low_footprint = int(np.round(ctp_low_footprint, decimals=0))
        ctp_high_footprint = int(np.round(ctp_high_footprint, decimals=0))

        # Append footprint array into larger list
        ctf_low.append(ctf_low_footprint)
        ctf_high.append(ctf_high_footprint)
        ctp_low.append(ctp_low_footprint)
        ctp_high.append(ctp_high_footprint)

    # Convert larger lists to numpy arrays
    ctf_low = np.asarray(ctf_low, dtype="object")
    ctf_high = np.asarray(ctf_high, dtype="object")
    ctp_low = np.asarray(ctp_low, dtype="object")
    ctp_high = np.asarray(ctp_high, dtype="object")

    return ctf_low, ctf_high, ctp_low, ctp_high


# Once all netCDFs have been downloaded, process them.
def Process(FILES):
    for FILE in FILES:
        print(f'Now processing file: {FILE}')
        # Extract the satellite identifier from netCDF filename.
        satName = FILE.split('_')[2]

        # Create nc object
        nc = xr.open_dataset(FILE, decode_times=False)
        temperature = np.array(nc.Temperature)
        wvcd = np.array(nc.H2O)
        p_layer = np.array(nc.Effective_Pressure[0, :])
        plev = np.array(nc.Pressure[0, :])
        psurf = np.array(nc.Surface_Pressure)
        nobs = len(nc.Latitude)
        topography = np.array(nc.Topography)
        time = np.array(nc.Time)

        # Find the total cloud top fraction and cloud top pressure for each footprint.
        # This will get written to the .csv file.
        cloud_top_fraction = np.array(nc.Cloud_Top_Fraction)
        cloud_top_pressure = np.array(nc.Cloud_Top_Pressure)

        # Remove NaNs or -9999
        cloud_top_fraction = cloud_top_fraction[~np.isnan(cloud_top_fraction)]
        cloud_top_fraction = cloud_top_fraction[cloud_top_fraction != -9999]
        cloud_top_pressure = cloud_top_pressure[~np.isnan(cloud_top_pressure)]
        cloud_top_pressure = cloud_top_pressure[cloud_top_pressure != -9999]

        # Reshape the arrays
        cloud_top_fraction = np.reshape(cloud_top_fraction, (-1, 2))
        cloud_top_pressure = np.reshape(cloud_top_pressure, (-1, 2))

        # Make new arrays of the cloud top pressure and fraction at each layer.
        ctf_low, ctf_high, ctp_low, ctp_high = find_ctf_ctp(nobs, cloud_top_fraction, cloud_top_pressure)

        blmult, botlev = get_botlev_blmult(plev, psurf, nobs)

        # Find wvcd and temperature at surface using BLMULT
        wvcd_sfc = calc_wvcd_sfc(botlev, blmult, nobs, p_layer, wvcd)
        tsfc = calc_Tsfc(botlev, blmult, nobs, plev, temperature)

        # Insert surface values into temperature, water vapor and pressure arrays.
        blmult_P_ALL = insert_surface_pressure(botlev, nobs, plev, psurf)
        blmult_T_ALL = insert_surface_temperature(botlev, nobs, plev, tsfc, temperature)
        blmult_wvcd_ALL = insert_surface_water_vapor(botlev, nobs, plev, wvcd_sfc, wvcd)

        # Derive the other variables
        wvmr = convert_cd2mr(nobs, blmult_wvcd_ALL, blmult_P_ALL, psurf, botlev)
        tv = calc_virtual_temperature(nobs, blmult_P_ALL, wvmr, blmult_T_ALL)
        mslp = calc_mslp(nobs, topography, tsfc, psurf)
        z = calc_geopotential_height(nobs, blmult_P_ALL, mslp, tv)
        dew_point = calc_dewpoint(nobs, blmult_P_ALL, blmult_T_ALL, wvmr)
        #######################################

        # Convert temperature to Celsius before writing to text files.
        for x in range(len(blmult_T_ALL)):
            blmult_T_ALL[x] = blmult_T_ALL[x] - 273.15

        ########################
        # Create the CSV
        srcid = []

        if satName == 'm01' or satName == 'm02' or satName == 'm03':
            scenes = nc.Number_of_Dice.values
        elif satName == 'j01' or satName == 'npp':
            scenes = nc.Number_of_CrIS_FORs.values
        elif satName == 'aq0':
            scenes = nc.Number_of_AIRS_FORs.values

        for FOR in scenes:
            if satName == 'm01' or satName == 'm02' or satName == 'm03':
                tmp = nc.sel(Number_of_Dice=FOR, drop=True)
            elif satName == 'j01' or satName == 'npp':
                tmp = nc.sel(Number_of_CrIS_FORs=FOR, drop=True)
            elif satName == 'aq0':
                tmp = nc.sel(Number_of_AIRS_FORs=FOR, drop=True)

            temps = blmult_T_ALL[FOR]
            dewPoint = dew_point[FOR]
            press = blmult_P_ALL[FOR]
            Z = z[FOR]
            lat = np.round(tmp.Latitude.values, 2)
            lon = np.round(tmp.Longitude.values, 2)
            sfcHgt = int(tmp.Topography.values)
            qc_flag = int(tmp.Quality_Flag.values)

            # Find the date/time info of each footprint
            epoch_time = time[FOR] / 1000 # convert milliseconds to seconds.
            footprint_time = str(datetime.datetime.utcfromtimestamp(epoch_time))
            year = footprint_time[2:4]
            month = footprint_time[5:7]
            day = footprint_time[8:10]
            hour = footprint_time[11:13]
            minute = footprint_time[14:16]
            seconds = footprint_time[17:19]
            date = f'{year}{month}{day}'
            timestamp = f'{hour}{minute}{seconds}'

            # Convert QC flag to a color based on which sensors pass/fail.
            if qc_flag == 0: # successful retrieval
                qc_flag = 'green'
            elif qc_flag == 1 or qc_flag == 17: # ir+mw failed, mw-only passed
                qc_flag = 'yellow'
            elif qc_flag == 9 or qc_flag == 25: # ir+mw failed and mw-only failed
                qc_flag = 'red'

            # Write the srcid and elev to the CSV
            ID = FOR + 1
            ID = '{:003d}'.format(ID)
            srcid.append(f'{date}_{timestamp}_{ID}_{satName}')

            # Append the name, lat, lon, elev and srcid to the .csv file.
            csvEntry = ['','','',f'{str(srcid[FOR])}','','',f'{str(lat)}',f'{str(lon)}',f'{str(sfcHgt)}',f'{str(qc_flag)}',f'{str(srcid[FOR])}']
            write_obj = open(os.path.join(HOME, '.sharppy', 'datasources', f'{satName}_case_study.csv'), "a+", newline='')
            csv_writer = writer(write_obj)
            csv_writer.writerow(csvEntry)
            write_obj.close()

            ########################
            # Create the Text Files
            headers = []
            headers.append("%TITLE%")
            headers.append(f"STC {str(date)}/{str(timestamp)} {str(lat)},{str(lon)}")
            headers.append("")
            vars = ["LEVEL", "HGHT", "TEMP", "DWPT", "WDIR", "WSPD"]

            line=''
            for var in vars:
                line = f'{line}       {var}'

            headers.append(line)
            headers.append("-------------------------------------------------------------------")
            headers.append('%RAW%')

            text_file = open(os.path.join(HOME, '.sharppy', 'datasources', satName, f'{date}_{timestamp}_{str(FOR+1).zfill(3)}_{satName}.txt'), "w")

            for header in headers:
                text_file.write(f'{header}\n')

            for i in np.arange(len(temps)-1, -1, -1):
                HGHT = Z[i]
                TEMP = temps[i]
                LEVEL = press[i]
                DWPT = dewPoint[i]
                WDIR = -9999.00
                WSPD = -9999.00

                if (DWPT > TEMP):
                    DWPT = TEMP

                LEVEL = '{:.2f}'.format(round(LEVEL, 2))
                HGHT = '{:.2f}'.format(round(HGHT, 2))
                TEMP = '{:.2f}'.format(round(TEMP, 2))
                DWPT = '{:.2f}'.format(round(DWPT, 2))
                WDIR = '{:.2f}'.format(round(WDIR, 2))
                WSPD = '{:.2f}'.format(round(WSPD, 2))
                text_file.write(f"{LEVEL.rjust(12, ' ')},{HGHT.rjust(10, ' ')},{TEMP.rjust(10, ' ')},{DWPT.rjust(10, ' ')},{WDIR.rjust(10, ' ')},{WSPD.rjust(10, ' ')}\n")
            text_file.write('%END%')
            text_file.write('\n')
            text_file.write(f'ctf_low {ctf_low[FOR]},ctf_high {ctf_high[FOR]},ctp_low {ctp_low[FOR]},ctp_high {ctp_high[FOR]}')
            text_file.close()


######################
######## MAIN ########
######################
if __name__ == '__main__':
    # j01 = NOAA-20
    # npp = Suomi-NPP
    # m02 = MetOp-A
    # m01 = MetOp-B
    # m03 = MetOp-C
    # aq0 = Aqua

    # Assign command line arguments to a list.
    satNames = sys.argv[1:]

    for satName in satNames:
        # Check that user-supplied arguments are valid.
        if satName == 'j01' or satName == 'aq0' or satName == 'npp' or satName == 'm01' or satName == 'm02' or satName == 'm03':
            create_text_file_path()
            remove_old_txt_csv()
            write_csv_header()

            # Process the text files
            FILES = glob.glob(f'NUCAPS-EDR*{satName}*.nc')
            Process(FILES)
        else:
            print(f'{satName} is not a valid satellite identifier.  Valid arguments are "j01", "aq0", "npp", "m01", "m02", "m03".')
            print('Continuing with next argument.')
            continue

    data_path = os.path.join(HOME, '.sharppy', 'datasources')
    print(f'Script has completed.  Your data has been stored under {data_path}')
