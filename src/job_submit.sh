#!/bin/bash


if [ -f ./bklist.txt ];then
rm ./bklist.txt
fi

if [ -f input.lst ];then
rm input.lst
fi

echo 'making bklist.txt'
bklist='./bklist.txt'
ls -d 0* >> $bklist

echo 'making input.lst'
for diri in `cat $bklist`
do
cp ./DispInv.sh $diri/DispInv.sh
cd $diri
./DispInv.sh
rm ./DispInv.sh
cd ..
done
