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
echo "cd $diri; ./DispInv.sh" >> input.lst
done

mpifort run_shell_mpi.f90 -o run_shell_mpi

echo 'execute:  mpirun -np $num run_shell_mpi'







