## MEDSLIK-II NECCTON

MEDSLIK-II NECCTON is a customised version of the open-source MEDSLIK-II oil spill model documented [`MEDSLIK_II`](https://www.medslik-ii.org/). The code is based on Version 1.01 of the historical code, which works for the Mediterranean Sea.
MEDSLIK-II simulates oil spill trajectory and fate due to three main processes known collectively as weathering: viscous–gravity spreading, evaporation, and natural dispersion. The formation of water- in-oil emulsion is also considered. 
If oil arrives on the coast, the model simulates the adsorption of oil droplets into the coastal environment while considering the probability that oil may be washed back into the water.


### Installation

The currently supported architectures is Linux. The details of installation can be found in the Quick-Start Guide and Manual published at [`MEDSLIK_II`](https://www.medslik-ii.org/).

### Required data

`MEDSLIK-II NECCTON` requires the following data:

* Sea current and SST data: [MEDSEA_MULTIYEAR_PHY_006_004](https://data.marine.copernicus.eu/product/MEDSEA_MULTIYEAR_PHY_006_004/)
* Wind data: [ECMWF ERA5 Wind](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels)
* Bathymetry: [medf_.bath](https://www.medslik-ii.org/) 
* Coastline: [medf_.map](https://www.medslik-ii.org/)

### License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
For more details see [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html) 
