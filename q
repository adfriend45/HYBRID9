cd SOURCE
rm HYBRID9.exe
make HYBRID9.exe
cd ..
cd EXECUTE
rm HYBRID9.exe
cd ..
cp SOURCE/HYBRID9.exe EXECUTE
cd EXECUTE
nice -19 ./HYBRID9.exe
