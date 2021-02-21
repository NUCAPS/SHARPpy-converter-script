# Downloading and Running Offline NUCAPS Data in SHARPpy
This tutorial explains how to run NUCAPS in SHARPpy in offline mode. The steps here outline how to generate NUCAPS text files from netCDF files that are stored locally.
Anaconda and SHARPpy are assumed to be already installed before trying the steps below.

## Downloading data for Case Studies

1. Order the NUCAPS Environmental Data Records (EDR) data from [NOAA CLASS](class.noaa.gov), under the [JPSS Sounder Products (JPSS_SND)](https://www.avl.class.noaa.gov/saa/products/search?sub_id=0&datatype_family=JPSS_SND&submit.x=28&submit.y=11) dropdown menu. If you are unfamiliar with CLASS, refer to the [CLASS Tutorial document](https://weather.msfc.nasa.gov/nucaps/resources_training.html) on how to download data.

## Process NUCAPS Data Locally for Case Studies

### Creating SHARPpy-formatted text files from NUCAPS netCDF files

2. Download and save the [sharppy_offline_netcdf_converter.py](https://github.com/NUCAPS/SHARPpy) script and save to an easily accessible directory. This script will convert the netCDF files to a format that SHARPpy can read.

3. After downloading the EDR files, move them to the same directory as *sharppy_offline_netcdf_converter.py*.

4. Run the *sharppy_offline_netcdf_converter.py* script to process all the netCDFs in the current directory. This creates the sounding text and location csv files. These files will be saved in */home/{user}/.sharppy/datasources*.

### Updating SHARPpy to point to the case study files

5. Change your directory to */home/{user}/SHARPpy/datasources* and create a new xml file (i.e. *case_study.xml*) that will point to your local text file paths.
   * Do not modify standard.xml, only use it as a template.
   * Change *standard.xml* to *standard.xml_IGNORE* so that SHARPpy does not read from this file. The original *standard.xml* points to real-time data.

The code may look like:

```bash
cd /home/<user>/SHARPpy/datasources
cp standard.xml case_study.xml
mv standard.xml standard.xml_IGNORE
cd /home/<user>/.sharppy/datasources
mv standard.xml standard.xml_IGNORE
```

6. In *case_study.xml*, add the below datasource xml tag and change the url parameter to point to the text file locations. The "datasource name=" parameter should be either: "NUCAPS CONUS NOAA-20" or "NUCAPS CONUS Suomi-NPP" based on the data you are using.

For NOAA-20, the code will look like:

```xml
 <datasource name="NUCAPS CONUS NOAA-20" ensemble="false" observed="true">
     <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/j01/{srcid}.txt" format="spc" >
         <time first="0" range="48" delta="1" offset="6" delay="4" cycle="12" archive="24" start="-" end="-"/>
         <points csv="j01_case_study.csv" />
     </outlet>
 </datasource>
 ```

 For Suomi-NPP data, the code will look like:

 ```xml
 <datasource name="NUCAPS CONUS Suomi-NPP" ensemble="false" observed="true">
     <outlet name="STC" url="file:///home/<user>/.sharppy/datasources/npp/{srcid}.txt" format="spc" >
         <time first="0" range="48" delta="1" offset="6" delay="4" cycle="12" archive="24" start="-" end="-"/>
         <points csv="npp_case_study.csv" />
     </outlet>
 </datasource>
 ```

7. For the XML changes to take effect, you need to reinstall SHARPpy. Go to the folder that contains your SHARPpy install. For example, in the terminal, type:

```bash
cd /home/<user>/SHARPpy/
python setup.py install
```


8. Lastly, launch the SHARPpy GUI by typing *sharppy* into the terminal:

```bash
sharppy
```
