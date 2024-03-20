# Downloading and Running Offline NUCAPS Data in SHARPpy
This tutorial explains how to run NUCAPS in SHARPpy in offline mode. The steps here outline how to generate NUCAPS text files from netCDF files that are stored locally.  Anaconda and SHARPpy are assumed to be already installed before trying the steps below.

## Downloading Data for Case Studies

1. Order the NUCAPS Environmental Data Records (EDR) data from [NOAA CLASS](class.noaa.gov), under the [JPSS Sounder Products (JPSS_SND)](https://www.avl.class.noaa.gov/saa/products/search?sub_id=0&datatype_family=JPSS_SND&submit.x=28&submit.y=11) dropdown menu. If you are unfamiliar with CLASS, refer to the [CLASS Tutorial document](https://weather.msfc.nasa.gov/nucaps/resources_training.html) on how to download data. More recently, NUCAPS data are distributed on the [cloud via the NOAA NODD Program](https://registry.opendata.aws/noaa-jpss/).

## Process NUCAPS Data Locally for Case Studies

### Creating SHARPpy-formatted text files from NUCAPS netCDF files

2. Download and save the [sharppy_offline_netcdf_converter.py](https://github.com/NUCAPS/SHARPpy-converter-script) script to an easily accessible directory. This script will convert the netCDF files to a format that SHARPpy can read.  Install the xarray and netcdf4 Python libraries needed by the script.

```bash
conda install xarray
conda install h5netcdf
```

3. After downloading the EDR files, move them to the same directory as *sharppy_offline_netcdf_converter.py*.

4. Run the *sharppy_offline_netcdf_converter.py* script from the command line to process all the netCDFs in the current directory. This creates the sounding text and location csv files. These files will be saved in */home/{user}/.sharppy/datasources*.  Add the appropriate satellite identifier(s) as arguments after the script name and press ENTER.

```bash
# "j01" for NOAA20
# "j02" for NOAA21
# "npp" for Suomi-NPP
# "m02" for Metop-A
# "m01" for Metop-B
# "m03" for Metop-C
# "aq0" for Aqua

# For example, run the script on NOAA-20 EDRs only...
python sharppy_offline_netcdf_converter.py j01
```

### Updating SHARPpy to point to the case study files

5. Change your directory to */home/{user}/SHARPpy/datasources* and open *case_study.xml*.  In *case_study.xml*, uncomment the datasource tag(s) you'll be using and change the URL to point where the text files reside.

For NOAA-20, the code will look like:

```xml
<datasource name="NUCAPS Case Study NOAA-20" ensemble="false" observed="true">
    <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/j01/{srcid}.txt" format="nucaps" >
        <time first="0" range="0" delta="0" offset="0" delay="1" cycle="1200" archive="12" start="-" end="-"/>
        <points csv="j01_case_study.csv" />
    </outlet>
</datasource>
```

For NOAA-21, the code will look like:

```xml
<datasource name="NUCAPS Case Study NOAA-21" ensemble="false" observed="true">
    <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/j02/{srcid}.txt" format="nucaps" >
        <time first="0" range="0" delta="0" offset="0" delay="1" cycle="1200" archive="12" start="-" end="-"/>
        <points csv="j02_case_study.csv" />
    </outlet>
</datasource>
```

For Suomi-NPP data, the code will look like:

```xml
<datasource name="NUCAPS Case Study Suomi-NPP" ensemble="false" observed="true">
    <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/npp/{srcid}.txt" format="nucaps" >
        <time first="0" range="0" delta="0" offset="0" delay="1" cycle="1200" archive="12" start="-" end="-"/>
        <points csv="npp_case_study.csv" />
    </outlet>
</datasource>
```

For Aqua data, the code will look like:

```xml
<datasource name="NUCAPS Case Study Aqua" ensemble="false" observed="true">
    <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/aq0/{srcid}.txt" format="nucaps" >
        <time first="0" range="0" delta="0" offset="0" delay="1" cycle="1200" archive="12" start="-" end="-"/>
        <points csv="aq0_case_study.csv" />
    </outlet>
</datasource>
```

For MetOp-A data, the code will look like:

```xml
<datasource name="NUCAPS Case Study MetOp-A" ensemble="false" observed="true">
    <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/m02/{srcid}.txt" format="nucaps" >
        <time first="0" range="0" delta="0" offset="0" delay="1" cycle="1200" archive="12" start="-" end="-"/>
        <points csv="m01_case_study.csv" />
    </outlet>
</datasource>
```

For MetOp-B data, the code will look like:

```xml
<datasource name="NUCAPS Case Study MetOp-B" ensemble="false" observed="true">
    <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/m01/{srcid}.txt" format="nucaps" >
        <time first="0" range="0" delta="0" offset="0" delay="1" cycle="1200" archive="12" start="-" end="-"/>
        <points csv="m02_case_study.csv" />
    </outlet>
</datasource>
```

For MetOp-C data, the code will look like:

```xml
<datasource name="NUCAPS Case Study MetOp-C" ensemble="false" observed="true">
    <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/m03/{srcid}.txt" format="nucaps" >
        <time first="0" range="0" delta="0" offset="0" delay="1" cycle="1200" archive="12" start="-" end="-"/>
        <points csv="m03_case_study.csv" />
    </outlet>
</datasource>
```

6. For the XML changes to take effect, you need to reinstall SHARPpy.  Go to the folder that contains your SHARPpy install. For example, in the terminal, type:

```bash
cd /home/<user>/SHARPpy
python setup.py install
```

7. Lastly, launch the SHARPpy GUI by typing *sharppy* into the terminal:

```bash
sharppy
```
