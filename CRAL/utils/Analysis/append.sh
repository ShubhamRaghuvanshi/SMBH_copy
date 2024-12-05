#!/bin/bash

SimDir=/ptmp/mpa/sraghu/sink_nodrag/star_sn_sink

nlogfile=1
SimLogFile=$SimDir/sinkdata_back/sink_0000$nlogfile.csv
Sinkfile_append=$SimDir/sinkdata_new/sink_0000$nlogfile.csv
while [[ -f "$SimLogFile" && $nlogfile -lt 10 ]]
do 
  echo "SinkLogFile : $SimLogFile"
  cat "$SimLogFile" "$Sinkfile_append" > "$SimDir/sinkdata_all/sinkappend_0000$nlogfile.csv"
  mv "$SimDir/sinkdata_all/sinkappend_0000$nlogfile.csv" "$SimDir/sinkdata_all/sink_0000$nlogfile.csv"
  nlogfile=$((nlogfile + 1))
  SimLogFile=$SimDir/sinkdata_back/sink_0000$nlogfile.csv
  Sinkfile_append=$SimDir/sinkdata_new/sink_0000$nlogfile.csv
done 

nlogfile=10
SimLogFile=$SimDir/sinkdata_back/sink_000$nlogfile.csv
Sinkfile_append=$SimDir/sinkdata_new/sink_000$nlogfile.csv
while [[ -f "$SimLogFile" && $nlogfile -lt 100 ]]
do 
  echo "SinkLogFile : $SimLogFile"
  cat "$SimLogFile" "$Sinkfile_append" > "$SimDir/sinkdata_all/sinkappend_000$nlogfile.csv"
  mv "$SimDir/sinkdata_all/sinkappend_000$nlogfile.csv" "$SimDir/sinkdata_all/sink_000$nlogfile.csv"
  nlogfile=$((nlogfile + 1))
  SimLogFile=$SimDir/sinkdata_back/sink_000$nlogfile.csv
  Sinkfile_append=$SimDir/sinkdata_new/sink_000$nlogfile.csv
done 

nlogfile=100
SimLogFile=$SimDir/sinkdata_back/sink_00$nlogfile.csv
Sinkfile_append=$SimDir/sinkdata_new/sink_00$nlogfile.csv
while [[ -f "$SimLogFile" && $nlogfile -lt 1000 ]]
do 
  echo "SinkLogFile : $SimLogFile"
  cat "$SimLogFile" "$Sinkfile_append" > "$SimDir/sinkdata_all/sinkappend_00$nlogfile.csv"
  mv "$SimDir/sinkdata_all/sinkappend_00$nlogfile.csv" "$SimDir/sinkdata_all/sink_00$nlogfile.csv"
  nlogfile=$((nlogfile + 1))
  SimLogFile=$SimDir/sinkdata_back/sink_00$nlogfile.csv
  Sinkfile_append=$SimDir/sinkdata_new/sink_00$nlogfile.csv
done 
