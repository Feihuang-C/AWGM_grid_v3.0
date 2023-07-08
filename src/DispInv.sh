#!/bin/bash
# 

ls bk.*.deg.*.dsp > dsplist
for dsp in `cat dsplist`
do
#azm=`echo $dsp | awk -F '.' '{print $4}'`
rfpath='../'                     # path of RF
rfile=stacksac2.0AJWA			 # RF file
intMod=../mod.d2				 # initial model
invScriptF=../jointinversion.sh	 # file1 for setting inversion parameters
jobi=$dsp.job 					 # file2 for setting inversion parameters

cp   $intMod ./mod.d         
cp   $rfpath/$rfile  ./
cp   $invScriptF ./
echo $rfile>rfnt

echo "0.00499999989 0.00499999989 0. 0.00499999989 0." > $jobi
echo "1 2 2 2 2 2 2 0 1 0" >> $jobi
echo "mod.d" >> $jobi
echo $dsp >> $jobi
echo "rfnt" >> $jobi

mv $jobi jobs.d
./jointinversion.sh > log
cp modl.out $dsp.modo3

ps2pdf figjnt1.eps $dsp.fit1.pdf

rm tmp*
rm *.eps
rm *PLT
rm  rfnt mod.d modl.out 
done
rm log
