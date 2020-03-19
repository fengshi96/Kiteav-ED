#!/bin/bash -x
make
Kz_max=1.0

for Kz in $(seq 1.0 0.1 $Kz_max)
do
	mkdir -p Kzz_$Kz
	cd Kzz_$Kz
	cp ../input_template.inp input.inp
	sed -i 's/kzkzkz/'$Kz'/g' input.inp

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#PBS -N ED24-$Kz
#PBS -l nodes=1:ppn=10
#PBS -l walltime=64:00:00 
#PBS -l vmem=150gb			
#PBS -l mem=80gb
#PBS -m a				
#PBS -M feng.934@osu.edu
#PBS -A PAS1528	
#PBS -j oe			
#PBS -q batch
hostname
#PBS -r n
module load intel
module load mkl
cd \$PBS_O_WORKDIR
time
../main input.inp &> runForinput.cout
time
EOF
)"
	echo "$rawjob" &> job.pbs
	echo -n "Kz = $Kz, "
	qsub job.pbs
	
	cd ..
done
