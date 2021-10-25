# -*- coding: utf-8 -*-

import os
import arcpy
from arcpy import env
import math
import datetime


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
        self.label = "Pukaskwa and Fire Hazard"
        self.alias = "pukaskwafirehazard"

        # List of tool classes associated with this toolbox
        self.tools = [ImportTool, CreateNormal, CreateClimateScenario, CreateInterpolations, Mapping, Layout]


class ImportTool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "1 Import"
        self.description = """Import Tool creates a Pukaskwa Fire Hazard file geodatabase and a Cartography feature dataset, modifies the weather observation 
                            text file and import it into the geodatabase, and import all the feature classes from the source data folder."""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""       
        #location of text file with weather station data
        Weather_txt = arcpy.Parameter(
            displayName="Input Weather Stations Data",
            name="Weather_txt",
            datatype="DETextfile",
            parameterType="Required",
            direction="Input")

        #source data folder with shapefiles and feature classes for cartography
        Source_data = arcpy.Parameter(
            displayName="Source Data Folder",
            name="Source_data",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        #new file geodatabase location
        Geodatabase_folder = arcpy.Parameter(
            displayName="New Geodatabase Location",
            name="Geodatabase Location",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        #file geodatabase output parameter
        Geodatabase = arcpy.Parameter(
            displayName="Pukaskwa Fire Hazard Geodatabase",
            name="Geodatabase",
            datatype="DEWorkspace",
            parameterType="Derived",
            direction="Output")

        parameters = [Weather_txt, Source_data, Geodatabase_folder, Geodatabase]
        return parameters


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""
        #assigning parameters to variables
        Weather_txt = parameters[0].valueAsText
        Source_data = parameters[1].valueAsText
        Geodatabase_folder = parameters[2].valueAsText

        #environment settingS
        env.overwriteOutput = True
        env.workspace = Geodatabase_folder

        #take the weather observations text file, create a new Month field and output the modified data in a new text file named Raw_Weather.txt
        with open(Weather_txt, "r") as infile:
            with open("Raw_Weather.txt", "w") as outfile:
                for row in infile:
                    striprow = row.rstrip()
                    splitrow = striprow.split(',')
                    if splitrow[0] == "WX_STN_NUMBER":
                        outfile.write(striprow + ',MONTH\n')
                        continue
                    datesplit = splitrow[15].split('/')

                    newrow = (f"{striprow},{datesplit[0]}\n")
                    outfile.write(newrow)
        arcpy.AddMessage("Raw Weather text file with Month field created.")

        #create new file geodatabase named Pukaskwa_FireHazard
        GDB_name = "Pukaskwa_FireHazard"
        FireHazard_GDB = arcpy.management.CreateFileGDB(Geodatabase_folder, GDB_name, "CURRENT")
        arcpy.AddMessage("File Geodatabase created.")

        #create feature dataset Cartography in the Pukaskwa_FireHazard file geodatabase
        Dataset_name = "Cartography"
        Spatial_reference = arcpy.SpatialReference(26916)   #NAD83 / UTM Zone 16N
        Cartography_dataset = arcpy.management.CreateFeatureDataset(FireHazard_GDB, Dataset_name, Spatial_reference)
        arcpy.AddMessage("Feature Dataset created.")
        
        #import Raw Weather text file to Pukaskwa_FireHazard geodatabase
        arcpy.conversion.TableToGeodatabase('Raw_Weather.txt', FireHazard_GDB)
        arcpy.AddMessage("Weather stations data imported.")

        #import all the shapefiles and feature classes in the Pukaskwa_FireHazard geodatabase that are contained in the source data folder
        FC_list = []
        walk = arcpy.da.Walk(Source_data, datatype='FeatureClass')
        for dirpath, dirnames, filenames in walk:
            for filename in filenames:
                FC_list.append(os.path.join(dirpath, filename))
        arcpy.conversion.FeatureClassToGeodatabase(FC_list, Cartography_dataset)
        arcpy.AddMessage("Feature classes from the source data folder imported.")

        #setting derived output parameter
        Output_workspace = os.path.join(Geodatabase_folder, "Pukaskwa_FireHazard.gdb")
        parameters[3].value = Output_workspace

        return None




class CreateNormal(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "2 Create Normal Weather"
        self.description = """Create Normal Weather Tool creates a Weather Stations feature dataset, a new Weather Stations feature class with new geometry and a 
                        new Station ID field and High Hazard fields for each month from April to October. The tool also removes rows in the Raw Weather Table where 
                        the Fire Weather Index, Temperature and Build-Up Index fields are null, it then split the Raw Weather Table into individual tables for each 
                        weather station (Normal Tables). Finally the tool calculates the number of normal days with high hazard for each month at each stations and 
                        add the values in the Weather Stations feature class."""
        self.canRunInBackground = False


    def getParameterInfo(self):
        """Define parameter definitions"""
        #workspace geodatabase
        Wkspace_geodatabase = arcpy.Parameter(
            displayName="Pukaskwa FireHazard Geodatabase",
            name="Wkspace_geodatabase",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
    
        #file geodatabase output parameter
        Geodatabase = arcpy.Parameter(
            displayName="Pukaskwa Fire Hazard Geodatabase",
            name="Geodatabase",
            datatype="DEWorkspace",
            parameterType="Derived",
            direction="Output")

        parameters = [Wkspace_geodatabase, Geodatabase]
        return parameters


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""
        #assigning parameter to variable
        Wkspace_geodatabase = parameters[0].valueAsText

        #environment settingS
        env.overwriteOutput = True
        env.workspace = Wkspace_geodatabase

        #create feature dataset WeatherStations
        Dataset_name = 'WeatherStations'
        Spatial_reference = arcpy.SpatialReference(26916)   #NAD83 / UTM Zone 16N
        Weatherstations_dataset = arcpy.management.CreateFeatureDataset(Wkspace_geodatabase, Dataset_name, Spatial_reference)
        arcpy.AddMessage("Feature Dataset created.")

        #create weather stations feature class in WeatherStations feature dataset
        Featureclass_name = 'WeatherStations'
        Geometry_type = "Point"
        Weatherstations_fc = arcpy.management.CreateFeatureclass(Weatherstations_dataset, Featureclass_name, Geometry_type)
        arcpy.AddMessage("Weather Stations feature class created.")

        #create list of station IDs containing list of their easting and northing values using the most recent entries in the Raw Weather table
        DUB = [667800, 5354800]
        LUR = [585800, 5366600]
        PUK = [552000, 5383000]
        WAW = [666000, 5316000]
        WHR = [627200, 5383500]
        Weatherstations_list = [DUB, LUR, PUK, WAW, WHR]

        #create point geometry objects in the WeatherStations feature class for each weather station
        for weatherstation in Weatherstations_list:
            point = arcpy.Point(weatherstation[0], weatherstation[1])
            with arcpy.da.InsertCursor(Weatherstations_fc, "SHAPE@") as cursor:
                cursor.insertRow([point])
        arcpy.AddMessage("Weather Stations geometry created.")

        #add Station_ID field in the WeatherStations feature class
        Field_type = "Text"
        ID_field = "Station_ID"
        arcpy.management.AddField(Weatherstations_fc, ID_field, Field_type)

        #insert Station_ID value in Station_ID field
        Stationnames_list = ["DUB", "LUR", "PUK", "WAW", "WHR"]
        index = 0
        with arcpy.da.UpdateCursor(Weatherstations_fc, ID_field) as cursor:
            for row in cursor:
                row[0] = Stationnames_list[index]
                cursor.updateRow(row)
                index += 1
        arcpy.AddMessage("Station ID field in the Weather Stations feature class created and populated.")

        #removing rows from Raw Weather table where the Fire Weather Index is null (4315 rows)
        Table = "Raw_Weather"
        Field_name = "FWI"
        with arcpy.da.UpdateCursor(Table, Field_name) as cursor:
            for row in cursor:
                if row[0] == None:
                    cursor.deleteRow()

        #removing rows from Raw Weather table where the Temp is null (69 rows)
        Field_name = "TEMP"
        with arcpy.da.UpdateCursor(Table, Field_name) as cursor:
            for row in cursor:
                if row[0] == None:
                    cursor.deleteRow()
        arcpy.AddMessage("Rows with null values deleted.")

        #changing values to 0 from Raw Weather table where the DMC and BUI is null (4 rows)
        Field_names = ['DMC', 'BUI']
        with arcpy.da.UpdateCursor(Table, Field_names) as cursor:
            for row in cursor:
                if row[0] == None:
                    dmc = 0
                    bui = 0
                    Row = [dmc, bui]
                    cursor.updateRow(Row)

        #create tables for each weather station based on the Raw_Weather table
        Table_names = ["DUB_Normal", "LUR_Normal", "PUK_Normal", "WAW_Normal", "WHR_Normal"]
        Stationnames_list = ["DUB", "LUR", "PUK", "WAW", "WHR"]
        Table = "Raw_Weather"
        Field_name = "WSTNID"
        Delim_field = arcpy.AddFieldDelimiters(Table, Field_name)
        index = 0
        for name in Table_names:
            SQL_where = f"{Delim_field} = '{Stationnames_list[index]}'"
            arcpy.conversion.TableToTable(Table, Wkspace_geodatabase, name, SQL_where)
            index += 1
        arcpy.AddMessage("Tables for each weather station created and populated.")

        #add fields for each month for the normal days with FWI values >= 10 to the Weather Station feature class
        New_fields = ["HH_Norm_Apr", "HH_Norm_May", "HH_Norm_Jun", "HH_Norm_Jul","HH_Norm_Aug","HH_Norm_Sep","HH_Norm_Oct"]
        for name in New_fields:
            Field_type = "DOUBLE" 
            arcpy.management.AddField(Weatherstations_fc, name, Field_type)
        arcpy.AddMessage("Normal High Hazard field added in the Weather Stations feature class.")

        #calculate how many rows in each weather station tables have a FWI values >= 10 and how many rows there are for each month
        Table_names = ["DUB_Normal", "LUR_Normal", "PUK_Normal", "WAW_Normal", "WHR_Normal"]
        for table in Table_names:
            Field_names = ["FWI", "MONTH"]
            High_hazard = [0,0,0,0,0,0,0,0,0,0,0,0]
            Days_total = [0,0,0,0,0,0,0,0,0,0,0,0]
            with arcpy.da.SearchCursor(table, Field_names) as cursor:
                for row in cursor:
                    if row[0] >= 10:
                        High_hazard[row[1]-1] += 1
                    Days_total[row[1]-1] += 1
            
            #calculate the average number of day per month with FWI >= 10
            Days_permonth = [31,28,31,30,31,30,31,31,30,31,30,31]
            HHrate_list = []
            FullHHrate_list = []
            Month_list = range(12)
            for month in Month_list:
                #to avoid dividing by 0
                if Days_total[month] == 0:
                    HH_rate = 0
                    HHrate_list.append(HH_rate)
                    continue
                HH_rate = High_hazard[month]/Days_total[month] * Days_permonth[month]
                HHrate_list.append(HH_rate)
            FullHHrate_list.append(HHrate_list)

            #adding value from the high hazard rate list to the Weather Stations feature class attribute table
            FC = Weatherstations_fc
            Update_fields = ["HH_Norm_Apr", "HH_Norm_May", "HH_Norm_Jun", "HH_Norm_Jul","HH_Norm_Aug","HH_Norm_Sep","HH_Norm_Oct"]
            Field_name = 'Station_ID'
            Delim_field = arcpy.AddFieldDelimiters(FC, Field_name)
            SQL_where = f"{Delim_field} = '{table[0:3]}'"
            Row = (HHrate_list[3], HHrate_list[4], HHrate_list[5], HHrate_list[6], HHrate_list[7], HHrate_list[8], HHrate_list[9])
            with arcpy.da.UpdateCursor(FC, Update_fields, SQL_where) as cursor:
                for row in cursor:
                    cursor.updateRow(Row)
        arcpy.AddMessage(f"Fields for the number of normal days of high hazard per month for each weather stations updated.")

        #setting derived output parameter
        parameters[1].value = Wkspace_geodatabase

        return None




class CreateClimateScenario(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "3 Create Climate Scenario"
        self.description = """Create Climate Scenario Tool creates new tables for each weather station and add new fields for new weather values. Using the increase 
                        of temperature specified, the tool recalculates the temperature, relative humidity, fire weather codes and fire weather indices. The tool 
                        then add new high hazard fields in the weather station feature class and calculates a new values for the number of days of high hazard for 
                        each month."""
        self.canRunInBackground = False


    def getParameterInfo(self):
        """Define parameter definitions"""
        #geodatabase workspace
        Wkspace_geodatabase = arcpy.Parameter(
            displayName="Pukaskwa FireHazard Geodatabase",
            name="Wkspace_geodatabase",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
    
        #year of the forecasted climate scenario
        Year = arcpy.Parameter(
            displayName="Year of forecast",
            name="Year",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        #increase of temperature for Dubreuilville (DUB)
        NewTemp_DUB = arcpy.Parameter(
            displayName="Increase of average temperature for Dubreuilville (DUB)",
            name="NewTemp_DUB",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        #increase of temperature for Lurch Creek (LUR)
        NewTemp_LUR = arcpy.Parameter(
            displayName="Increase of average temperature for Lurch Creek (LUR)",
            name="NewTemp_LUR",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        #increase of temperature for Pukaskwa (PUK)
        NewTemp_PUK = arcpy.Parameter(
            displayName="Increase of average temperature for Pukaskwa (PUK)",
            name="NewTemp_PUK",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        #increase of temperature for Wawa (WAW)
        NewTemp_WAW = arcpy.Parameter(
            displayName="Increase of average temperature for Wawa (WAW)",
            name="NewTemp_WAW",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        #increase of temperature for White River (WHR)
        NewTemp_WHR = arcpy.Parameter(
            displayName="Increase of average temperature for White River (WHR)",
            name="NewTemp_WHR",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        #file geodatabase output parameter
        Geodatabase = arcpy.Parameter(
            displayName="Pukaskwa Fire Hazard Geodatabase",
            name="Geodatabase",
            datatype="DEWorkspace",
            parameterType="Derived",
            direction="Output")

        parameters = [Wkspace_geodatabase, Year, NewTemp_DUB, NewTemp_LUR, NewTemp_PUK, NewTemp_WAW, NewTemp_WHR, Geodatabase]
        return parameters


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""

        #assigning parameter to variable
        Wkspace_geodatabase = parameters[0].valueAsText
        Year = parameters[1].valueAsText
        NewTemp_DUB = float(parameters[2].valueAsText)
        NewTemp_LUR = float(parameters[3].valueAsText)
        NewTemp_PUK = float(parameters[4].valueAsText)
        NewTemp_WAW = float(parameters[5].valueAsText)
        NewTemp_WHR = float(parameters[6].valueAsText)

        #environment settingS
        env.overwriteOutput = True
        env.workspace = Wkspace_geodatabase
        
        #create table of weather data for each weather station for the climate scenario
        Weatherstations_list = ['DUB', 'LUR', 'PUK', 'WAW', 'WHR']
        Newtemp_list = [NewTemp_DUB, NewTemp_LUR, NewTemp_PUK, NewTemp_WAW, NewTemp_WHR]
        index = 0
        Tablename_list = []
        for weatherstation in Weatherstations_list:
            Input_table = Weatherstations_list[index] + '_Normal'
            Table_name = f"{Weatherstations_list[index]}_{Year}"
            index += 1
            Tablename_list.append(Table_name)
            arcpy.conversion.TableToTable(Input_table, Wkspace_geodatabase, Table_name)
        arcpy.AddMessage('New tables for each weather station created.')

        #add weather fields in each weather station table
        for table in Tablename_list:
            Input_table = table
            Fieldname_list = [f"TEMP_{Year}", f"DEW_{Year}", f"RH_{Year}", f"WIND_{Year}", f"RAIN_{Year}", f"FFMC_{Year}", f"DMC_{Year}", f"DC_{Year}", f"ISI_{Year}", f'BUI_{Year}', f'FWI_{Year}']
            Field_type = "FLOAT"
            for fieldname in Fieldname_list:
                arcpy.management.AddField(Input_table, fieldname, Field_type)
        arcpy.AddMessage('New weather and fire hazard fields added to the new climate scenario tables.')

        #add new temperature values in the new TEMP field
        Newtemp_list = [NewTemp_DUB, NewTemp_LUR, NewTemp_PUK, NewTemp_WAW, NewTemp_WHR]
        index = 0
        for table in Tablename_list:
            Input_table = table
            Field_names = ['TEMP', f"TEMP_{Year}"]
            with arcpy.da.UpdateCursor(Input_table, Field_names) as cursor:
                    for row in cursor:
                        New_value = row[0] + Newtemp_list[index]
                        Row = [row[0], New_value]
                        cursor.updateRow(Row)
            index += 1

        #add dew points values in the new DEW field
        for table in Tablename_list:
            Input_table = table
            Field_names = ["TEMP", "REL_HUM", f"DEW_{Year}"]
            with arcpy.da.UpdateCursor(Input_table, Field_names) as cursor:
                    for row in cursor:
                        New_value = 243.04*(math.log(row[1]/100)+((17.625*row[0])/(243.04+row[0])))/(17.625-math.log(row[1]/100)-((17.625*row[0])/(243.04+row[0])))
                        Row = [row[0], row[1], New_value]
                        cursor.updateRow(Row)

        #add relative humidity values in the new RH field
        for table in Tablename_list:
            Input_table = table
            Field_names = [f"TEMP_{Year}", f"DEW_{Year}", f"RH_{Year}"]
            with arcpy.da.UpdateCursor(Input_table, Field_names) as cursor:
                    for row in cursor:
                        if row[1] == 0:
                            row[1] = 0.01
                        New_value = 100*(math.exp((17.625*row[1])/(243.04+row[1]))/math.exp((17.625*row[0])/(243.04+row[0])))
                        Row = [row[0], row[1], New_value]
                        cursor.updateRow(Row)

        #transfer values for wind speed and rain to new WIND and RAIN fields
        for table in Tablename_list:
            Input_table = table
            Field_names = ["WIND_SPEED", "RAIN", f"WIND_{Year}", f"RAIN_{Year}"]
            with arcpy.da.UpdateCursor(Input_table, Field_names) as cursor:
                    for row in cursor:
                        New_value_wind = row[0]
                        New_value_rain = row[1]
                        Row = [row[0], row[1], New_value_wind, New_value_rain]
                        cursor.updateRow(Row)
        arcpy.AddMessage("Weather fields for each climate scenario table calculated and populated.")

        ### Forest Fire Weather Index System Python Code from
        ### Wang, Y.; Anderson, K.R.; Suddaby, R.M. 2005
        ### Updated source code for calculation fire danger indes in the Canadian Forest Fire Weather Index System
        ### Natural Resources Canada, Canadian Forest Service. Northern Forest Centre
        ### Edmonton, AB. Inf. Rep. NOR-X-424

        ### ORIGINAL VARIABLES IN 1988 FORTRAN CODE ###

        # T = Noon temperature (°C)
        # H = Noon relative humidity (%)
        # W = Noon wind speed (km/h)
        # Ro = Rain fall in open over previous 24 hours, at noon (mm)
        # Rf = Effective rain fall for calculating FFMC
        # Re = Effective rain fall for calculating DMC
        # Rd = Effective rain fall for calculating DC

        # Mo = Fine Fuel Moisture Content from previous day
        # Mr = Fine Fuel Moisture Content after rain
        # M = Fine Fuel Moisture Content after drying
        # Ed = Fine Fuel equlibrium moisture content (EMC) for drying
        # Ew = Fine Fuel EMC for wetting
        # Ko = Intermediate step in calculation of Kd
        # Kd = Log drying rate, FFMC drying rate, FFMC in (M)/day
        # K1 = Intermediate step in calculation of Kw
        # Kw = Natural log wetting rate, ln (M)/day
        # Fo = Previous day's FFMC
        # F = Fine Fuel Moisture Code

        # DMo = Duff Moisture Content from previous day
        # DMr = Duff Moisture content after rain
        # DM = Duff moisture content after drying
        # K = Log drying rate in DMC, ln (M)/day
        # Le = Effective day length in DMC, hours
        # B = Slope variable in DMC rain effect
        # Po = Previous day's DMC
        # Pr = DMC after rain
        # P = DMC

        # Q = Moisture equivalent of DC, units of 0.254 mm
        # Qo = Moisture equivalent of previous day's DC
        # Qr = Moisture equivalent after rain
        # V = Potential evapotranspiration, units of 0.254 mm water/day
        # Lf = Day-length adjustement in DC
        # Do = Previous day's DC
        # Dr = DC after rain
        # D = DC

        # Fw = Wind function
        # Ff = Fine fuel moisture function
        # Fd = Duff moisture function
        # R = Initial spread index
        # U = Buildup index
        # B = Intermediate FWI
        # S = Final FWI

        """ Define Class FWI Class first"""
        class FWICLASS:
            def __init__(self, temp, rhum, wind, prcp):
                self.h = rhum
                self.t = temp
                self.w = wind
                self.p = prcp

            def FFMCcalc(self, ffmc0):
                mo = (147.2*(101.0 - ffmc0))/(59.5 + ffmc0)                                                                             #*Eq. 1*#
                if (self.p > 0.5):
                    rf = self.p - 0.5                                                                                                   #*Eq. 2*#
                    if (mo > 150.0):
                        mo = (mo+42.5*rf*math.exp(-100.0/(251.0-mo))*(1.0-math.exp(-6.93/rf)))+(0.0015*(mo-150.0)**2)*math.sqrt(rf)     #*Eq. 3b*#
                    elif mo <= 150.0:
                        mo = mo+42.5*rf*math.exp(-100.0/(251.0-mo))*(1.0-math.exp(-6.93/rf))                                            #*Eq. 3a*#
                    if (mo > 250.0):
                        mo = 250.0
                ed = 0.942*(self.h**0.679) + (11.0*math.exp((self.h-100.0)/10.0))+0.18*(21.1-self.t)*(1.0-1.0/math.exp(0.1150*self.h))  #*Eq. 4*#
                if (mo < ed):
                    ew = 0.618*(self.h**0.753)+(10.0*math.exp((self.h-100.0)/10.0))+0.18*(21.1-self.t)*(1.0-1.0/math.exp(0.115*self.h)) #*Eq. 5*#
                    if (mo <= ew):
                        kl = 0.424*(1.0-((100.0-self.h)/100.0)**1.7)+(0.0694*math.sqrt(self.w))*(1.0-((100.0-self.h)/100.0)**8)         #*Eq. 7a*#
                        kw = kl*(0.581*math.exp(0.0365*self.t))                                                                         #*Eq. 7b*#
                        m = ew-(ew - mo)/10.0**kw                                                                                       #*Eq. 9*#                                   
                    elif mo > ew:
                        m = mo
                elif (mo == ed):
                    m = mo
                elif mo > ed:
                    kl = 0.424*(1.0-(self.h/100.0)**1.7)+(0.0694*math.sqrt(self.w))*(1.0-(self.h/100.0)**8)                             #*Eq. 6a*#
                    kw = kl*(0.581*math.exp(0.0365*self.t))                                                                             #*Eq. 6b*#
                    m = ed+(mo-ed)/10.0**kw                                                                                             #*Eq. 8*#
                ffmc = (59.5*(250.0-m))/(147.2+m)                                                                                       #*Eq. 10*#
                if (ffmc > 101.0):
                    ffmc - 101.0
                if (ffmc <= 0.0):
                    ffmc = 0.0
                return ffmc

            def DMCcalc(self, dmc0, mth):
                el = [6.5,7.5,9.0,12.8,13.9,13.9,12.4,10.9,9.4,8.0,7.0,6.0]
                t = self.t
                if (t < -1.1):
                    t = -1.1
                rk = 1.894*(t+1.1)*(100.0-self.h)*(el[mth-1]*0.0001)                                                                    #*Eqs. 16 & 17*#
                if self.p > 1.5:
                    ra = self.p
                    rw = 0.92*ra-1.27                                                                                                   #*Eq. 11*#
                    wmi = 20.0 + 280.0/math.exp(0.023*dmc0)                                                                             #*Eq. 12*#
                    if dmc0 <= 33.0:
                        b = 100.0/(0.5+0.3*dmc0)                                                                                        #*Eq. 13a*#
                    elif dmc0 > 33.0:
                        if dmc0 <= 65.0:
                            b = 14.0-1.3*math.log(dmc0)                                                                                 #*Eq. 13b*#
                        elif dmc0 > 65.0:
                            b = 6.2 * math.log(dmc0)-17.2                                                                                   #*Eq. 13c*#
                    wmr = wmi+(1000*rw)/(48.77+b*rw)                                                                                     #*Eq. 14*#
                    pr = 43.43*(5.6348-math.log(wmr-20.0))                                                                              #*Eq. 15*#
                elif self.p <= 1.5:
                    pr = dmc0
                if (pr<0.0):
                    pr = 0.0
                dmc = pr + rk
                if (dmc <= 1.0):
                    dmc = 1.0
                return dmc

            def DCcalc(self, dc0, mth):
                global dc
                fl = [-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5.0,2.4,0.4,-1.6,-1.6]
                t = self.t
                if (t < -2.8):
                    t = -2.8
                pe = (0.36*(t+2.8)+fl[mth-1])/2                                                                                         #*Eq. 22*#        
                if pe <= 0.0:
                    pe = 0.0
                if (self.p > 2.8):
                    ra = self.p
                    rw = 0.83*ra-1.27                                                                                                   #*Eq. 18*#
                    smi = 800.0*math.exp(-dc0/400.0)                                                                                    #*Eq. 19*#
                    dr = dc0 - 400.0*math.log(1.0+((3.937*rw)/smi))                                                                     #*Eqs. 20 & 21*#
                    if (dr > 0.0):
                        dc = dr + pe
                elif self.p <= 2.8:
                    dc = dc0 + pe
                return dc

            def ISIcalc(self, ffmc):
                mo = 147.2*(101.0-ffmc)/(59.5+ffmc)                                                                                     #*Eq. 1*#
                ff = 19.115*math.exp(mo*-0.1386)*(1.0+(mo**5.31)/49300000.0)                                                            #*Eq. 25*#
                isi = ff*math.exp(0.05039*self.w)                                                                                       #*Eq. 26*#
                return isi

            def BUIcalc(self, dmc, dc):
                if dmc <= 0.4*dc:
                    bui = (0.8*dc*dmc)/(dmc+0.4*dc)                                                                                     #*Eq. 27a*#
                else:
                    bui = dmc-(1.0-0.8*dc/(dmc+0.4*dc))*(0.92+(0.0114*dmc)**1.7)                                                        #*Eq. 27b*#        
                if bui < 0.0:
                    bui = 0.0
                return bui

            def FWIcalc(self, isi, bui):
                if bui <= 80.0:
                    bb = 0.1*isi*(0.626*bui**0.809+2.0)                                                                                 #*Eq. 28a*#
                else:
                    bb = 0.1*isi*(1000.0/(25.0+108.64/math.exp(0.023*bui)))                                                             #*Eq. 28b*#        
                if (bb <= 1.0):
                    fwi = bb
                else:
                    fwi = math.exp(2.72*(0.434*math.log(bb))**0.647)                                                                    #*Eq. 30a*#        
                return fwi

        ##calculate new fire hazard codes and indices
        for table in Tablename_list:
            Input_table = table
            Field_names = ['DFS_WOBS_ARCHIVE_WX_YEAR', 'DFS_WOBS_ARCHIVE_WX_DATE', 'MONTH', 'FFMC', 'DMC', 'DC', 'ISI', 'BUI', 'FWI', f"TEMP_{Year}", f"RH_{Year}", f"WIND_{Year}", f"RAIN_{Year}", f"FFMC_{Year}", f"DMC_{Year}", f"DC_{Year}", f"ISI_{Year}", f'BUI_{Year}', f'FWI_{Year}']
            
            #setting time variables
            year = 0
            date = datetime.date.today()
            One_day = datetime.timedelta(days=+1)

            with arcpy.da.UpdateCursor(Input_table, Field_names) as cursor:
                for row in cursor:
                    #if year change, starting indices remain the same because it is a new season and row before is unrelated to this row
                    if year != row[0]:
                        year = row[0]
                        ffmc0 = row[3]
                        dmc0 = row[4]
                        dc0 = row[5]
                        isi = row[6]
                        bui = row[7]
                        fwi = row[8]
               
                        Row = [row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], ffmc0, dmc0, dc0, isi, bui, fwi]
                        cursor.updateRow(Row)

                    
                    #if row before is not from the day before, starting indices remain the same as the normal (otherwise it creates unrealistic values and 
                    #cummulatively high overestimation of hazard)
                    if date + One_day != row[1]:
                        date = row[1]
                        ffmc0 = row[3]
                        dmc0 = row[4]
                        dc0 = row[5]
                        isi = row[6]
                        bui = row[7]
                        fwi = row[8]
               
                        Row = [row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], ffmc0, dmc0, dc0, isi, bui, fwi]
                        cursor.updateRow(Row)

                    #calculating codes and indices with NRCan FWICLASS
                    else:
                        year = row[0]
                        date = row[1]
                        mth = row[2]
                        temp = row[9]
                        rhum = row[10]
                        wind = row[11]
                        prcp = row[12]

                        if rhum > 100.0:
                            rhum = 100.0

                        fwisystem = FWICLASS(temp, rhum, wind, prcp)
                        ffmc = fwisystem.FFMCcalc(ffmc0)
                        dmc = fwisystem.DMCcalc(dmc0, mth)
                        dc = fwisystem.DCcalc(dc0, mth)
                        isi = fwisystem.ISIcalc(ffmc)
                        bui = fwisystem.BUIcalc(dmc, dc)
                        fwi = fwisystem.FWIcalc(isi, bui)
                    
                        ffmc0 = ffmc
                        dmc0 = dmc
                        dc0 = dc

                        Row = [year, row[1], mth, row[3], row[4], row[5], row[6], row[7], row[8], temp, rhum, wind, prcp, round(ffmc, 1), round(dmc, 1), round(dc, 1), round(isi, 1), round(bui, 1), round(fwi, 1)]
                        cursor.updateRow(Row)
        arcpy.AddMessage(f"Fire hazards fields for each climate scenario tables calculated and populated.")

        #add fields for each month for the normal days with FWI values >= 10 to the Weather Station feature class
        Weatherstations_fc = r'WeatherStations\WeatherStations'
        New_fields = [f"HH_{Year}_Apr", f"HH_{Year}_May", f"HH_{Year}_Jun", f"HH_{Year}_Jul",f"HH_{Year}_Aug",f"HH_{Year}_Sep",f"HH_{Year}_Oct"]
        for name in New_fields:
            Field_type = "DOUBLE" 
            arcpy.management.AddField(Weatherstations_fc, name, Field_type)
        arcpy.AddMessage(f"High Hazard fields for year {Year} added in the Weather Stations feature class.")

        #calculate how many rows in each weather station tables have a FWI values >= 10 and how many rows there are for each month
        Table_names = [f"DUB_{Year}", f"LUR_{Year}", f"PUK_{Year}", f"WAW_{Year}", f"WHR_{Year}"]
        for table in Table_names:
            Field_names = [f"FWI_{Year}", "MONTH"]
            High_hazard = [0,0,0,0,0,0,0,0,0,0,0,0]
            Days_total = [0,0,0,0,0,0,0,0,0,0,0,0]
            with arcpy.da.SearchCursor(table, Field_names) as cursor:
                for row in cursor:
                    if row[0] >= 10:
                        High_hazard[row[1]-1] += 1
                    Days_total[row[1]-1] += 1
            
            #calculate the average number of day per month with FWI >= 10
            Days_permonth = [31,28,31,30,31,30,31,31,30,31,30,31]
            HHrate_list = []
            FullHHrate_list = []
            Month_list = range(12)
            for month in Month_list:
                #to avoid dividing by 0
                if Days_total[month] == 0:
                    HH_rate = 0
                    HHrate_list.append(HH_rate)
                    continue
                HH_rate = High_hazard[month]/Days_total[month] * Days_permonth[month]
                HHrate_list.append(HH_rate)
            FullHHrate_list.append(HHrate_list)

            #adding value from the high hazard rate list to the Weather Stations feature class attribute table
            FC = Weatherstations_fc
            Update_fields = [f"HH_{Year}_Apr", f"HH_{Year}_May", f"HH_{Year}_Jun", f"HH_{Year}_Jul",f"HH_{Year}_Aug",f"HH_{Year}_Sep",f"HH_{Year}_Oct"]
            Field_name = 'Station_ID'
            Delim_field = arcpy.AddFieldDelimiters(FC, Field_name)
            SQL_where = f"{Delim_field} = '{table[0:3]}'"
            Row = (HHrate_list[3], HHrate_list[4], HHrate_list[5], HHrate_list[6], HHrate_list[7], HHrate_list[8], HHrate_list[9])
            with arcpy.da.UpdateCursor(FC, Update_fields, SQL_where) as cursor:
                for row in cursor:
                    cursor.updateRow(Row)
        arcpy.AddMessage(f"Fields for the number of days of high hazard per month for the year {Year} for each weather stations updated.")

        #setting derived output parameter
        parameters[7].value = Wkspace_geodatabase

        return None




class CreateInterpolations(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "4 Create Interpolations"
        self.description = """Create Interpolation Tool creates weighted inversed distance interpolations for each field in the Weather Station feature class. The values represent the number of high hazard days."""
        self.canRunInBackground = False


    def getParameterInfo(self):
        """Define parameter definitions"""
        #geodatabase workspace
        Wkspace_geodatabase = arcpy.Parameter(
            displayName="Pukaskwa FireHazard Geodatabase",
            name="Wkspace_geodatabase",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
    
        #file geodatabase output parameter
        Geodatabase = arcpy.Parameter(
            displayName="Pukaskwa Fire Hazard Geodatabase",
            name="Geodatabase",
            datatype="DEWorkspace",
            parameterType="Derived",
            direction="Output")

        parameters = [Wkspace_geodatabase, Geodatabase]
        return parameters


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""      
        #assigning parameter to variable
        Wkspace_geodatabase = parameters[0].valueAsText

        #environment settings
        env.overwriteOutput = True
        env.workspace = Wkspace_geodatabase
        env.extent = arcpy.Extent(531700, 5274400, 687600, 5435300)
        
        #making a list of high hazard fields in the Weather Stations feature class
        fc = r'WeatherStations\WeatherStations'
        fields = arcpy.ListFields(fc, 'HH*')
        HHfields_list = []
        for field in fields:
            HHfields_list.append(field.name)

        #create interpolations for each field in the Weather Station feature class using weighted inversed distance interpolation
        for field in HHfields_list:
            Inpoint_fc = fc
            z_field = field
            outIDW = arcpy.sa.Idw(fc, z_field)
            outIDW.save(field)
            arcpy.AddMessage(f"Interpolation surface for {field} created.")

        #setting derived output parameter
        parameters[1].value = Wkspace_geodatabase
        
        return None




class Mapping(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "5 Mapping"
        self.description = """Mapping Tool remove all the layers already present on the maps, it adds each interpolation layer on their respective climate scenario 
                            maps and symbolizes them using consistent class breaks and colour ramp. The tool adds layer files to create a basemap on each map, 
                            then it adds and symbolizes the Weather Station feature class."""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        #gworkspace geodatabase
        Wkspace_geodatabase = arcpy.Parameter(
            displayName="Pukaskwa FireHazard Geodatabase",
            name="Wkspace_geodatabase",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        #Folder containing layer files
        Layer_folder = arcpy.Parameter(
            displayName="Layer Files Folder",
            name="Layer_folder",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
    
        #file geodatabase output parameter
        Geodatabase = arcpy.Parameter(
            displayName="Pukaskwa Fire Hazard Geodatabase",
            name="Geodatabase",
            datatype="DEWorkspace",
            parameterType="Derived",
            direction="Output")

        parameters = [Wkspace_geodatabase, Layer_folder, Geodatabase]
        return parameters


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""      

        #assigning parameters to variables
        Wkspace_geodatabase = parameters[0].valueAsText
        Layer_folder = parameters[1].valueAsText
        Geodatabase = parameters[2].valueAsText

        #environment settingS
        env.overwriteOutput = True
        env.workspace = Wkspace_geodatabase

        #create a set for map names, the final set contains an element 'Norm' for the normal and additional element for each climate scenario
        FC_location = "WeatherStations/WeatherStations"
        fc = os.path.join(Wkspace_geodatabase, FC_location)
        Map_set = set()
        Fields_list = arcpy.ListFields(fc, 'HH*')
        for field in Fields_list:
            Field_name = field.name.split('_')
            Map_set.add(Field_name[1])

        #transfer the map name elements from the Map_set set to a list
        Mapname_list = []
        for element in Map_set:
            Mapname_list.append(element)

        #create project reference
        aprx = arcpy.mp.ArcGISProject("Current")
        
        #create layer file references
        BoundaryPUK = os.path.join(Layer_folder, "BoundaryPUK.lyrx")
        Road = os.path.join(Layer_folder, "Road.lyrx")
        Water = os.path.join(Layer_folder, "Water.lyrx")
        WaterPUK = os.path.join(Layer_folder, "WaterPUK.lyrx")

        #list of all layers to be used for the basemap
        Layer_list = [ Road, Water, WaterPUK, BoundaryPUK]

        #for each element in the map set, extract a map and create a map object
        count = 0
        for map in Map_set:
            try:
                m = aprx.listMaps()[count]
                m.name = f"HH_{Mapname_list[count]}"
            #empty maps need to be created in the GUI for the Mapping Tool to work successfully
            except:
                arcpy.AddWarning("Pukaskwa Fire Hazard Model needs 6 empty maps.")
                arcpy.AddError("CREATE MAPS IN THE PROJECT: The project needs to contain as many empty maps as there are 'Create Climate Scenario' tools plus one in the model.")
                sys.exit()

            #remove all layers that are already in the maps 
            #(for newly created maps: the topographic basemap or if the tool as already been run and other layers are still present)
            layers = m.listLayers()
            for layer in layers:
                m.removeLayer(layer)

            #walking through the workspace geodatabase to create a list of all the interpolations 
            Raster_list = []
            walk = arcpy.da.Walk(Wkspace_geodatabase, datatype='RasterDataset')
            for dirpath, dirnames, filenames in walk:
                for filename in filenames:
                    #creating list of interpolation surfaces for the normal or each climate scenario
                    if map in filename:
                        Complete_filename = f"{filename}/Band_1"
                        path = os.path.join(Wkspace_geodatabase, Complete_filename)
                        Raster_list.append(path)
            
            #adding interpolations on their respective maps
            for path in Raster_list:
                m.addDataFromPath(path)

            #finding the minimum value in the interpolation rasters to applied the proper symbolization
            #(if no values are in the lower class break, symbolization is not working properly)
            Interpolation_rasters = m.listLayers('HH*')
            for raster in Interpolation_rasters:
                vegras = arcpy.sa.Raster(raster.dataSource)
                minimum = 30
                for i,j in vegras:
                    if vegras[i,j] < minimum:
                        minimum = vegras[i,j]

                #for symbology, if minimum raster value is lower than 5, apply 6 class breaks
                if minimum < 5:
                    sym = raster.symbology
                    sym.updateColorizer('RasterClassifyColorizer')

                    sym.colorizer.breakCount = 6

                    #list of labels for the legend and colours for the colour ramp
                    Labels_list = ("0 - 5 days", "6 - 10 days", "11 - 15 days", "16 - 20 days", "21 - 25 days", "26 - 30 days")
                    Colour_list = [[0, 114, 6, 100], [101, 171, 20, 100], [196, 227, 29, 100], [249, 223, 26, 100], [225, 154, 11, 100], [252, 59, 9, 100]]
                    Break_value = 5
                    Label_index = 0

                    #applying class break value, labels and colours
                    for brk in sym.colorizer.classBreaks:
                        brk.upperBound = Break_value
                        brk.label = f"{Labels_list[Label_index]}"
                        brk.color = {'RGB' : Colour_list[Label_index]}
                        Break_value += 5
                        Label_index += 1
                    
                    raster.symbology = sym
    
                #for symbology, if minimum raster value is lower than 10, apply 5 class breaks
                elif minimum < 10:
                    sym = raster.symbology
                    sym.updateColorizer('RasterClassifyColorizer')

                    sym.colorizer.breakCount = 5

                    #list of labels for the legend and colours for the colour ramp
                    Labels_list = ("6 - 10 days", "11 - 15 days", "16 - 20 days", "21 - 25 days", "26 - 30 days")
                    Colour_list = [[101, 171, 20, 100], [196, 227, 29, 100], [249, 223, 26, 100], [225, 154, 11, 100], [252, 59, 9, 100]]
                    Break_value = 10
                    Label_index = 0

                    #applying class break value, labels and colours
                    for brk in sym.colorizer.classBreaks:
                        brk.upperBound = Break_value
                        brk.label = f"{Labels_list[Label_index]}"
                        brk.color = {'RGB' : Colour_list[Label_index]}
                        Break_value += 5
                        Label_index += 1
                    
                    raster.symbology = sym

                #for symbology, if minimum raster value is lower than 15, apply 4 class breaks
                elif minimum < 15:
                    sym = raster.symbology
                    sym.updateColorizer('RasterClassifyColorizer')

                    sym.colorizer.breakCount = 4

                    #list of labels for the legend and colours for the colour ramp
                    Labels_list = ("11 - 15 days", "16 - 20 days", "21 - 25 days", "26 - 30 days")
                    Colour_list = [[196, 227, 29, 100], [249, 223, 26, 100], [225, 154, 11, 100], [252, 59, 9, 100]]
                    Break_value = 15
                    Label_index = 0

                    #applying class break value, labels and colours
                    for brk in sym.colorizer.classBreaks:
                        brk.upperBound = Break_value
                        brk.label = f"{Labels_list[Label_index]}"
                        brk.color = {'RGB' : Colour_list[Label_index]}
                        Break_value += 5
                        Label_index += 1
                    
                    raster.symbology = sym

                #for symbology, if minimum raster value is lower than 20, apply 3 class breaks
                elif minimum < 20:
                    sym = raster.symbology
                    sym.updateColorizer('RasterClassifyColorizer')

                    sym.colorizer.breakCount = 3

                    #list of labels for the legend and colours for the colour ramp
                    Labels_list = ("16 - 20 days", "21 - 25 days", "26 - 30 days")
                    Colour_list = [[249, 223, 26, 100], [225, 154, 11, 100], [252, 59, 9, 100]]
                    Break_value = 20
                    Label_index = 0

                    #applying class break value, labels and colours
                    for brk in sym.colorizer.classBreaks:
                        brk.upperBound = Break_value
                        brk.label = f"{Labels_list[Label_index]}"
                        brk.color = {'RGB' : Colour_list[Label_index]}
                        Break_value += 5
                        Label_index += 1
                    
                    raster.symbology = sym

                #for symbology, if minimum raster value is lower than 25, apply 2 class breaks
                elif minimum < 25:
                    sym = raster.symbology
                    sym.updateColorizer('RasterClassifyColorizer')

                    sym.colorizer.breakCount = 2

                    #list of labels for the legend and colours for the colour ramp
                    Labels_list = ("21 - 25 days", "26 - 30 days")
                    Colour_list = [[225, 154, 11, 100], [252, 59, 9, 100]]
                    Break_value = 25
                    Label_index = 0

                    #applying class break value, labels and colours
                    for brk in sym.colorizer.classBreaks:
                        brk.upperBound = Break_value
                        brk.label = f"{Labels_list[Label_index]}"
                        brk.color = {'RGB' : Colour_list[Label_index]}
                        Break_value += 5
                        Label_index += 1
                    
                    raster.symbology = sym

            #add the layer files to the map
            for layer in Layer_list:
                lyrx = arcpy.mp.LayerFile(layer)
                m.addLayer(lyrx)
            
            #add the weather station feature class to the map and symbolize
            FC_location = "WeatherStations/WeatherStations"
            fc = os.path.join(Wkspace_geodatabase, FC_location)
            m.addDataFromPath(fc)
            Weather_stations = m.listLayers('Weather*')[0]
            sym = Weather_stations.symbology
            orange = {"RGB": [230, 152, 0, 100]}
            sym.renderer.symbol.color = orange
            sym.renderer.symbol.size = 8
            Weather_stations.symbology = sym
            count += 1

            arcpy.AddMessage(f"{m.name} Map created.")           


        #importing layout template to the project for next tool 
        #(the layer template is situated in the same folder as the layer files)
        Layout_template = os.path.join(Layer_folder, "FWILayout.pagx")
        aprx.importDocument(Layout_template)
        arcpy.AddMessage("Layout template file added to the project.")    

        #save mapping changes and delete project reference
        aprx.save()
        del aprx

        #setting derived output parameter
        parameters[2].value = Wkspace_geodatabase
    
        return None




class Layout(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "6 Layout"
        self.description = """Layout Tool is assigning maps with specific interpolation layer visiblility turned on for each map frames in a layout template. 
                            It creates a layout for each month with changing subtitles and then export it in the specified format in a specified folder."""
        self.canRunInBackground = False


    def getParameterInfo(self):
        """Define parameter definitions"""
        #workspace geodatabase
        Wkspace_geodatabase = arcpy.Parameter(
            displayName="Pukaskwa FireHazard Geodatabase",
            name="Wkspace_geodatabase",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")
    
        #climate scenario name string
        Scenario_name = arcpy.Parameter(
            displayName="Climate Scenario Name",
            name="Scenario_name",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        #output file format (PDF, JPEG OR PNG)
        Output_format = arcpy.Parameter(
            displayName="Output Format",
            name="Output_format",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        Output_format.filter.type = "ValueList"
        Output_format.filter.list = ['JPEG', 'PDF', 'PNG']

        #Export maps folder 
        Output_folder = arcpy.Parameter(
            displayName="Exported Map Folder",
            name="Output_folder",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        #file geodatabase output parameter
        Geodatabase = arcpy.Parameter(
            displayName="Pukaskwa Fire Hazard Geodatabase",
            name="Geodatabase",
            datatype="DEWorkspace",
            parameterType="Derived",
            direction="Output")

        parameters = [Wkspace_geodatabase, Scenario_name, Output_format, Output_folder, Geodatabase]
        return parameters


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return


    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""

        #assigning parameter to variable
        Wkspace_geodatabase = parameters[0].valueAsText
        Scenario_name = parameters[1].valueAsText
        Output_format = parameters[2].valueAsText
        Output_folder = parameters[3].valueAsText
        Geodatabase = parameters[4].valueAsText

        #create project reference
        aprx = arcpy.mp.ArcGISProject("Current")

        #creating layout reference
        try:
            layout = aprx.listLayouts("FWI*")[0]
        except:
            arcpy.AddWarning("The Layout Tool need to have the FWILayout file available in the project.")
            arcpy.AddError("OPEN THE LAYOUT FILE IN THE PROJECT.")
            sys.exit()

        #create a set for map names, the final set contains an element 'Norm' for the normal and additional element for each climate scenario
        FC_location = "WeatherStations/WeatherStations"
        fc = os.path.join(Wkspace_geodatabase, FC_location)
        Map_set = set()
        Fields_list = arcpy.ListFields(fc, 'HH*')
        for field in Fields_list:
            Field_name = field.name.split('_')
            Map_set.add(Field_name[1])

        #transfer the map name elements from the Map_set set to a list
        Mapname_list = []
        for element in Map_set:
            Mapname_list.append(element)
        count = 0

        #listsof month abbreviation and month complete name
        Months_abbreviation = ['Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct']
        Months_list = ["April", "May", "June", "July", "August", "September", "October"]

        #creating list of map and turning off the visibility of all the interpolations to ensure that only one will have its visibility turn on at a time
        for map in Mapname_list:
            Mapframesname_list = []
            m = aprx.listMaps("*" + map + "*")[0]
            Layers = m.listLayers("HH*")
            for layer in Layers:
                layer.visible = False

        #producing map for each month
        for month in Months_abbreviation:
            #creating a list of interpolation layers in each map, 
            #turning the interpolation visibility on for the current month map and all the other interpolation off
            for map in Mapname_list:
                m = aprx.listMaps("*" + map + "*")[0]
                Layers = m.listLayers("HH*")
                for layer in Layers:
                    if month in layer.name:
                        layer.visible = True
                    else:
                        layer.visible =  False

            #extracting the month name from the month list
            for month2 in Months_list:
                if month in month2:
                    Month_name = month2
            #extracting the subtitle text element and changing the text to write the layout month name and climate scenario
            Subtitle = layout.listElements("TEXT_ELEMENT", 'Subtitle')[0]
            Subtitle.text = f"{Month_name}, {Scenario_name}, Pukaskwa National Park"

            #placing the maps in their proper maps frames in the layout
            for map2 in Mapname_list:
                Map_frame = layout.listElements("MAPFRAME_ELEMENT", "*" + map2 + "*")[0]
                map3 = aprx.listMaps("*" + Map_frame.name[5:9] + "*")[0]
                Map_frame.map = map3

            #if output format is JPEG, export to JPEG format
            if Output_format == "JPEG":
                layout.exportToJPEG(os.path.join(Output_folder, month))
            
            #if output format is PDF, export to PDF format
            elif Output_format == "PDF":
                layout.exportToPDF(os.path.join(Output_folder, month))

            #if output format is PNG, export to PNG format
            elif Output_format == "PNG":
                layout.exportToPNG(os.path.join(Output_folder, month))
 
            arcpy.AddMessage(f"{Output_format} layout for the month of {Month_name} exported.")


        #save mapping changes and delete project reference
        aprx.save()
        del aprx

        #setting derived output parameter
        parameters[3].value = Wkspace_geodatabase

        return None