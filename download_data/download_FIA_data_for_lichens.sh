#! /bin/bash



# This script downloads FIA Plot, Subplot, and Condition data for eastern states into the DBDGS shared irectory.







cd /gpfs/nfs/share/ftp/priv/priv/Groups/DBDGS/FIA/RawData/




# Download Plot data

cd Plots

for s in AK AZ CA CO ID MT NV OR UT WA

do 

address=http://apps.fs.fed.us/fiadb-downloads/${s}_PLOT.CSV
curl $address > ${s}_PLOT.CSV

done





# Download Subplot data

cd ../Subplots


for s in AK AZ CA CO ID MT NV OR UT WA

do 

address=http://apps.fs.fed.us/fiadb-downloads/${s}_SUBPLOT.CSV
curl $address > ${s}_SUBPLOT.CSV

done



# Download Condition Data

cd ../Conditions



for s in AK AZ CA CO ID MT NV OR UT WA

do 

address=http://apps.fs.fed.us/fiadb-downloads/${s}_COND.CSV
curl $address > ${s}_COND.CSV

done


# Download Seedling Data

cd ../Seedlings


for s in AK AZ CA CO ID MT NV OR UT WA

do 

address=http://apps.fs.fed.us/fiadb-downloads/${s}_SEEDLING.CSV
curl $address > ${s}_SEEDLING.CSV

done


# Download Tree Data

cd ../Trees


for s in AK AZ CA CO ID MT NV OR UT WA

do 

address=http://apps.fs.fed.us/fiadb-downloads/${s}_TREE.ZIP
curl $address > ${s}_TREE.ZIP

done

for s in AK AZ CA CO ID MT NV OR UT WA

do 

unzip ${s}_TREE.ZIP

done

rm *.ZIP