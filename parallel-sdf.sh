#!/bin/sh

SDFCMD=sdf
TIMEOUT=1
SDFTMPDIR=.sdftmp
SDFJOBS=sdf.jobs
SDFJOBSTMP=sdf.jobs.tmp
SDFPIDS=sdf.pids
SDFPIDSTMP=sdf.pids.tmp

#sets up an environment for an sdf process
run_sdf()
{
	CURDIR=`pwd`

	input_geometry=$1
        output_volume=$2
        subvoldim=$3
        flip_normals=$4
        shift 4 #we must shift because solaris sh does not support ${10} or greater
        minx=$1 
        miny=$2
        minz=$3
        maxx=$4
        maxy=$5
        maxz=$6

	rm -rf ${minx}_${miny}_${minz}_${maxx}_${maxy}_${maxz} 
	mkdir ${minx}_${miny}_${minz}_${maxx}_${maxy}_${maxz}
	cd ${minx}_${miny}_${minz}_${maxx}_${maxy}_${maxz} 
	$SDFCMD $input_geometry $output_volume $subvoldim $flip_normals $minx $miny $minz $maxx $maxy $maxz > sdf.log &

	#return to the directory we were in prior to the invocation of SDF
	cd $CURDIR
}

add_to_job_queue()
{
	input_geometry=$1
	output_volume=$2
	subvoldim=$3
	flip_normals=$4
	shift 4 #we must shift because solaris sh does not support ${10} or greater
	minx=$1
        miny=$2
        minz=$3
        maxx=$4
        maxy=$5
        maxz=$6

	echo $input_geometry $output_volume $subvoldim $flip_normals $minx $miny $minz $maxx $maxy $maxz >> $SDFJOBS 
}

if [ $# -ne 12 ]; then
	echo Usage: $0 \<input geometry\> \<output volume\> \<dim of each subvolume \(volume size will be dim^3\)\> \<flipNormals \(0 or 1\)\> \<minx\> \<miny\> \<minz\> \<maxx\> \<maxy\> \<maxz\> \<number of divisions along each axis\> \<max number of processors\>
	exit 1
fi

input_geometry=$1
output_volume=$2
subvoldim=$3
flip_normals=$4
shift 4 #we must shift because solaris sh does not support ${10} or greater
orig_xmin=$1
orig_ymin=$2
orig_zmin=$3
orig_xmax=$4
orig_ymax=$5
orig_zmax=$6
num_divs=$7
maxprocesses=$8

echo Num Divisions: $num_divs
echo Max Processes: $maxprocesses


xspan=`echo "scale=8; ( $orig_xmax - $orig_xmin ) / $num_divs" | bc`
yspan=`echo "scale=8; ( $orig_ymax - $orig_ymin ) / $num_divs" | bc`
zspan=`echo "scale=8; ( $orig_zmax - $orig_zmin ) / $num_divs" | bc`

# make a temp directory if it doesn't exist
if [ ! -d $SDFTMPDIR ]; then
	mkdir $SDFTMPDIR 
fi

cd $SDFTMPDIR 
rm -f $SDFJOBS
rm -f $SDFPIDS

#build the job queue
idx_i=0
while [ $idx_i -ne $num_divs ]; do
	idx_j=0
	while [ $idx_j -ne $num_divs ]; do
		idx_k=0
		while [ $idx_k -ne $num_divs ]; do
		
			minx=`echo "scale=8; $orig_xmin + $idx_i * $xspan" | bc`
			miny=`echo "scale=8; $orig_ymin + $idx_j * $yspan" | bc`
			minz=`echo "scale=8; $orig_zmin + $idx_k * $zspan" | bc`

			maxx=`echo "scale=8; $orig_xmin + ( $idx_i + 1 ) * $xspan" | bc`
			maxy=`echo "scale=8; $orig_ymin + ( $idx_j + 1 ) * $yspan" | bc`
			maxz=`echo "scale=8; $orig_zmin + ( $idx_k + 1 ) * $zspan" | bc`

			echo ${minx}_${miny}_${minz}_${maxx}_${maxy}_${maxz}
			add_to_job_queue $input_geometry $output_volume $subvoldim $flip_normals $minx $miny $minz $maxx $maxy $maxz

			idx_k=`expr $idx_k + 1`
		done
		idx_j=`expr $idx_j + 1`
	done
	idx_i=`expr $idx_i + 1`
done

#now go through the queue and invoke processes according to the max number of processes requested
processes=0
curjob=`head -n 1 $SDFJOBS`
while [ -n "$curjob" ]; do
	if [ $processes -lt $maxprocesses ]; then
		# run the job
		echo `date` : Running job: $curjob
		run_sdf $curjob
		echo $! >> $SDFPIDS
		processes=`expr $processes + 1`

		# now remove the job from the queue
		numlines=`wc -l $SDFJOBS | awk '{ printf("%s",$1); }'`
		tail -n `expr $numlines - 1` $SDFJOBS > $SDFJOBSTMP 
		mv $SDFJOBSTMP $SDFJOBS
	else
		sleep $TIMEOUT	
		rm -f $SDFPIDSTMP
		touch $SDFPIDSTMP

		#For every pid in the file, check to see if that process is still running.
		#If not, remove it from the list
		curpid=`head -n 1 $SDFPIDS`
		curpididx=0
		while [ -n "$curpid" ]; do
			
			if [ -n "`ps -p $curpid --no-headers`" ]; then 
				echo $curpid >> $SDFPIDSTMP
			fi
			
			#get the next pid... perhaps there is a nicer way to extract a single line from a file?
			curpididx=`expr $curpididx + 1`	
			numlines=`wc -l $SDFPIDS | awk '{ printf("%s",$1); }'`
			numleft=`expr $numlines - $curpididx`	
			curpid=`tail -n $numleft $SDFPIDS | head -n 1`
		done
		mv $SDFPIDSTMP $SDFPIDS

		numlines=`wc -l $SDFPIDS | awk '{ printf("%s",$1); }'`
		processes=$numlines
	fi

	curjob=`head -n 1 $SDFJOBS`
done

# move all output rawiv files to the root sdftmp directory
#find . -name "*.rawiv" -exec echo \{\} \; | 
#  awk '{ cmd = sprintf("echo -n %s$ ",$1); system(cmd); cmd = sprintf("echo %s | perl -pi -e s/\\\\//_/g;",$1); system(cmd); }' | 
#  awk -F $ '{ cmd = sprintf("mv %s %s",$1,$2); system(cmd); }'

# TODO: combine rawiv files into a single volume
