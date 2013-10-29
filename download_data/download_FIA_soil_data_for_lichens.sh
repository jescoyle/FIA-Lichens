#! /bin/bash



# This script downloads FIA Soil table for states where lichen plots are located





cd C:/BulkDownload/FIA_SOIL



# Download Plot data



for s in AK AL AZ CA CO CT GA ID IL IN MA MD ME MI MN MT NC NH NJ NV OR PA RI SC UT VA VT WA WI WV WY

do 

address=http://apps.fs.fed.us/fiadb-downloads/${s}_SOILS_LAB.ZIP
curl $address > ${s}_SOILS_LAB.ZIP

address=http://apps.fs.fed.us/fiadb-downloads/${s}_SOILS_SAMPLE_LOC.ZIP
curl $address > ${s}_SAMPLE_LOC.ZIP

address=http://apps.fs.fed.us/fiadb-downloads/${s}_SOILS_VISIT.ZIP
curl $address > ${s}_SOILS_VISIT.ZIP

done





