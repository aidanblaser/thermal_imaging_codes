#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 11:31:44 2021

@author: aidanblaser
"""

import numpy as np
import subprocess as sp
import cv2 as cv
from PIL import Image

class ThermalImage:
    """Thermal Image class."""

    def __init__(self, image_path, camera_manufacturer, color_map="jet", thermal_np=None):
        """Base Class for Thermal Images.
        Args:
            image_path (str): Path of image to be loaded
            camera_manufacturer (str, optional): Which type of thermal camera was the image captured from.
                Supported values: ["flir","dji"].
            color_map (str, optional): The default colour map to be used when loading the image. Defaults to "jet".
            thermal_np (np.ndarray): Initialize directly with temp array.
        """
        self.image_path = image_path
        # Convert the string false colour map to an opencv object
        self.cmap = getattr(cv, f"COLORMAP_{color_map.upper()}")

        # Load the temperature matrix, sensor matrix and image metadata. First two are loaded as np arrays
        if camera_manufacturer.lower() == "flir":
            self.thermal_np, self.raw_sensor_np, self.meta = self.extract_temperatures_flir()
        elif camera_manufacturer.lower() == "dji":
            self.thermal_np, self.raw_sensor_np, self.meta = self.extract_temperatures_dji()
        elif camera_manufacturer.lower() == "pass":
            pass
        else:
            #logger.error(f"Cannot handle data from camera manufacturer {camera_manufacturer}")
            pass
        if thermal_np is not None:
            self.thermal_np = thermal_np

        self.global_min_temp = np.min(self.thermal_np)
        self.global_max_temp = np.max(self.thermal_np)

    def extract_temperatures_dji(self):
            """Extracts the DJI-encoded thermal image as 2D floating-point numpy array with temperatures in degC."""
            # read image metadata for the dji camera images
            exif_binary = "exiftool"
            meta = sp.Popen(
                (f'"{exif_binary}"  "{self.image_path}"'),
                shell=True,
                stdout=sp.PIPE,
            ).communicate()[0]
            meta = meta.decode("utf8")
    
            meta = {
                "Emissivity": 0.98,
                "ObjectDistance": 5,
                "AtmosphericTemperature": 25,
                "ReflectedApparentTemperature": 25,
                "IRWindowTemperature": 25,
                "IRWindowTransmission": 0.95,
                "RelativeHumidity": 50,
                "PlanckR1": 14906.216,#21106.77,
                "PlanckB": 1396.5,#1501,
                "PlanckF": 1,#1,
                "PlanckO": -7261,#-7340,
                "PlanckR2":0.010956882,#0.012545258,
            }
    
            ## calculating the raw sensor values
    
            im = Image.open(self.image_path)
            # concatenate APP3 chunks of data
            a = im.applist[3][1]
            for i in range(4, 14):
                a += im.applist[i][1]
            # create image from bytes
            img = Image.frombytes("I;16L", (640, 512), a)
    
            raw_sensor_np = np.array(img)
    
            ## extracting the temperatures from thermal images
            thermal_np = sensor_vals_to_temp(raw_sensor_np, **meta)
    
            return thermal_np, raw_sensor_np, meta
    
    
def sensor_vals_to_temp(
            raw,
            Emissivity='hello',
            ObjectDistance=5,
            AtmosphericTemperature=25,
            ReflectedApparentTemperature=20,
            IRWindowTemperature=20,
            IRWindowTransmission=1,
            RelativeHumidity=50,
            PlanckR1=14364.633,#21106.77,
            PlanckB= 1385.4,#1501,
            PlanckF=1,#1,
            PlanckO=-5753,#-7340,
            PlanckR2=0.010603162,#0.012545258,
            **kwargs,
            ):
            """Convert raw values from the thermographic sensor sensor to temperatures in °C. Tested for Flir cams."""
            # this calculation has been ported to python from https://github.com/gtatters/Thermimage/blob/master/R/raw2temp.R
            # a detailed explanation of what is going on here can be found there
            
            # constants
            ATA1 = 0.006569
            ATA2 = 0.01262
            ATB1 = -0.002276
            ATB2 = -0.00667
            ATX = 1.9
            
            # transmission through window (calibrated)
            emiss_wind = 1 - IRWindowTransmission
            refl_wind = 0
            # transmission through the air
            h2o = (RelativeHumidity / 100) * np.exp(
                1.5587
                + 0.06939 * (AtmosphericTemperature)
                - 0.00027816 * (AtmosphericTemperature) ** 2
                + 0.00000068455 * (AtmosphericTemperature) ** 3
            )
            tau1 = ATX * np.exp(-np.sqrt(ObjectDistance / 2) * (ATA1 + ATB1 * np.sqrt(h2o))) + (1 - ATX) * np.exp(
                -np.sqrt(ObjectDistance / 2) * (ATA2 + ATB2 * np.sqrt(h2o))
            )
            tau2 = ATX * np.exp(-np.sqrt(ObjectDistance / 2) * (ATA1 + ATB1 * np.sqrt(h2o))) + (1 - ATX) * np.exp(
                -np.sqrt(ObjectDistance / 2) * (ATA2 + ATB2 * np.sqrt(h2o))
            )
            # radiance from the environment
            raw_refl1 = PlanckR1 / (PlanckR2 * (np.exp(PlanckB / (ReflectedApparentTemperature + 273.15)) - PlanckF)) - PlanckO
            raw_refl1_attn = (1 - Emissivity) / Emissivity * raw_refl1  # Reflected component
            
            raw_atm1 = (
                PlanckR1 / (PlanckR2 * (np.exp(PlanckB / (AtmosphericTemperature + 273.15)) - PlanckF)) - PlanckO
            )  # Emission from atmosphere 1
            raw_atm1_attn = (1 - tau1) / Emissivity / tau1 * raw_atm1  # attenuation for atmospheric 1 emission
            
            raw_wind = (
                PlanckR1 / (PlanckR2 * (np.exp(PlanckB / (IRWindowTemperature + 273.15)) - PlanckF)) - PlanckO
            )  # Emission from window due to its own temp
            raw_wind_attn = (
                emiss_wind / Emissivity / tau1 / IRWindowTransmission * raw_wind
            )  # Componen due to window emissivity
            
            raw_refl2 = (
                PlanckR1 / (PlanckR2 * (np.exp(PlanckB / (ReflectedApparentTemperature + 273.15)) - PlanckF)) - PlanckO
            )  # Reflection from window due to external objects
            raw_refl2_attn = (
                refl_wind / Emissivity / tau1 / IRWindowTransmission * raw_refl2
            )  # component due to window reflectivity
            
            raw_atm2 = (
                PlanckR1 / (PlanckR2 * (np.exp(PlanckB / (AtmosphericTemperature + 273.15)) - PlanckF)) - PlanckO
            )  # Emission from atmosphere 2
            raw_atm2_attn = (
                (1 - tau2) / Emissivity / tau1 / IRWindowTransmission / tau2 * raw_atm2
            )  # attenuation for atmospheric 2 emission
            
            raw_obj = (
                raw / Emissivity / tau1 / IRWindowTransmission / tau2
                - raw_atm1_attn
                - raw_atm2_attn
                - raw_wind_attn
                - raw_refl1_attn
                - raw_refl2_attn
            )
            val_to_log = PlanckR1 / (PlanckR2 * (raw_obj + PlanckO)) + PlanckF
            if any(val_to_log.ravel() < 0):
                raise Exception("Image seems to be corrupted")
            # temperature from radiance
            return PlanckB / np.log(val_to_log) - 273.15

