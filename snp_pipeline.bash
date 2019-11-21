#!/bin/bash
# Place your gatech userame in the below export
export NAME="pparekh32"
get_input () {
	# Function for doing your getopts
	#
	# Input: Getopts array
	while getopts "a:b:r:f:o:evizh*" option
        do
                case $option in
                a)reads1=$OPTARG;;
                b)reads2=$OPTARG;;
                r)ref=$OPTARG;;
		#have allocated reallign variable
                e)realign=1;;
		#have allocated for millsFile
		f)millsFile=$OPTARG;;
		#have allocated output variable
                o)output=$OPTARG;;
                #have allocated a value for gunzip
	       	z)gunzip=1;;
		#have allocated a function to verbose
		v)v=1;;
                #have allocated a function to index
		i)index=1;;
		#have allocated a function to help
		h)h=1;;
		*)echo "Wrong Option"
			exit 1;;
        	esac
		if [ $h == 1 ]; then
			echo "This is the help file for the snp_pipeline.bash program
			      The commands that you can pass for getopts are:
			      Commands Necessary for the Program to run:
					-a [options] :To give the location of reads1 file.
					-b [options] :To give the location of reads2 file.
					-r [options] :To give the location of the reference file.
					-o [options] :To define the output file name.
			      Optional Commands:
					-e :To realign the BAM sorted file.
					-f [options]: Necessary after calling -e to give the location of the millsFile.
					-z :To have your output file unzipped.
					-v :To switch on the verbose mode.
					-i :To index your BAM file.
					-h :To access the help page."
			exit 0
		fi
	done
}
check_files () {
	# Function for checking for presence of input files, reference genome,
	# and the output VCF file
	#
	# Input: File locations (string)
	# Output: True, if checks pass; False, if checks fail (bool)
	#taking three flags to collectively return the value 1 as well as print the missing files
	# taking the thrid file to check the output file and stop the program if the user doesn't want to rewrite the files
  	if [ "$v" -eq 1 ]; then
		echo "Checking if Reads 1, Reads 2, Reference files exist"
		echo "Checking if the output file with the same name exists"
	fi
	if [ -e "$reads1" ]; then
                flag=1
        else
		echo "Missing Reads 1 file"
		exit 1
       	fi
        if [ -e "$reads2" ]; then
                flag1=1
        else
                echo "Missing Reads 2 file"
		exit 1
	fi
	if [ -e "$ref" ]; then
		flag2=1
	else
                echo "Missing Reference file"
         	exit 1
	fi
	if [ -e "$millsFile" ]; then
		flag3=1
	else
		echo "Missing Mills File"
		exit 1
	fi	
	if [ -e "$output".vcf.gz ]; then
		echo "Output file already exists"
		echo "Enter 0 to exit or enter any number for the output file to be overwritten"
		read -r ch
		if [ "$ch" == 0 ]; then
			exit 1
		else 
			flag4=1
		fi
	fi
	if [ "$flag" -eq 1 ] && [ "$flag1" -eq 1 ] && [ "$flag2" -eq 1 ] && [ "$flag3" -eq 1 ] && [ "$flag4" -eq 1 ];then
		return 0
	fi
}
prepare_temp () {
	# Preparing your temporary directory
	#
	#
        if [ $v -eq 1 ]; then
		echo "Checking if indexed file exists, if it doesn't creates a new one"
	fi
	if [ -e "$(pwd)"/data/chr17.fa.ann ]; then
		echo "Indexed File exists"
	else
		bwa index "$ref"
	fi
}
mapping () {
	#Function for the mapping step of the SNP-calling pipeline
		#
	# Input: File locations (string), Verbose flag (bool)
	# Output: File locations (string)
	if [ $v -eq 1 ]; then
		echo "Checks if the sorted BAM file exists, and if not creates one"
	fi
	if [ -e "$(pwd)"/tmp/out1_sorted.bam ]; then
		echo "Sorted BAM file already exits"
	else
		bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' "$ref" "$reads1" "$reads2" >"$(pwd)"/tmp/out.sam
		samtools fixmate -O bam "$(pwd)"/tmp/out.sam "$(pwd)"/tmp/out1.bam
		#without realigning we need the output file to pass to the call variants
		samtools sort -O bam -o "$(pwd)"/tmp/out1_sorted.bam -T /tmp/lane_temp "$(pwd)"/tmp/out1.bam
	fi
}
improvement () {
	# Function for improving the number of miscalls
	#
	# Input: File locations (string)
	# Output: File locations (string0)
	if [ $v -eq 1 ]; then
		echo "Checking whether .fai and .dict files exists and if it doesn't exists, creates them"
	fi
	if [ -e "$(pwd)"/data/chr17.fai ] && [ -e "$(pwd)"/data/chr17.dict ]; then
		echo ".fai and .dict file exists"
	else
		samtools faidx "$ref" -o "$(pwd)"/data/chr17.fai
		samtools dict "$ref" -o "$(pwd)"/data/chr17.dict
	fi
	if [ "$index" -eq 1 ];then
		if [ $v -eq 1 ]; then
			echo "Indexing sorted BAM file"
		fi
		samtools index "$(pwd)"/tmp/out1_sorted.bam
	fi
	if [ $v -eq 1 ]; then
		echo "Checking if -r is called"
	fi
	if [ "$realign" -eq 1 ]; then
		if [ $v -eq 1 ]; then
			echo "Checking if realigned file exists and if not creates one"
		fi
		if [ -e "$(pwd)"/tmp/out_realigned.bam ]; then
			echo "Realigned file exists"
		else
			java -Xmx2g -jar "$(pwd)"/lib/GenomeAnalysisTK.jar -T RealignerTargetCreator -R "$ref" -I "$(pwd)"/tmp/out1_sorted.bam -o "$(pwd)"/tmp/out1.intervals --known "$millsFile" --log_to_file "$(pwd)"/output/pparekh32.log
			java -Xmx4g -jar "$(pwd)"/lib/GenomeAnalysisTK.jar -I "$(pwd)"/tmp/out1_sorted.bam  -R "$ref" -T IndelRealigner -targetIntervals "$(pwd)"/tmp/out1.intervals -o "$(pwd)"/tmp/out_realigned.bam -known "$millsFile" --log_to_file "$(pwd)"/output/pparekh32_1.log
		fi
	else
		echo "Not Realigning"
	fi
}
call_variants () {
	# Function to call variants
	#
	# Input: File locations (string)
	# Ouput: None
	#to get the file name of the output vcf file
	if [ $v -eq 1 ]; then
		echo "Checks if the realign option was called and takes the input file accordingly"
	fi
	if [ "$realign" -eq 1 ]; then
		if [ "$gunzip" -eq 1 ]; then
	                if [ $v -eq 1 ]; then
                        	 echo "Gunzipped output"
                 	fi
			bcftools mpileup -Ou -f "$ref" "$(pwd)"/tmp/out_realigned.bam | bcftools call -vmO z -o "$output".vcf.gz
		else
			bcftools mpileup -Ou -f "$ref" "$(pwd)"/tmp/out_realigned.bam | bcftools call -vmO v -o "$output".vcf
	   	fi
	else
		if [ "$gunzip" -eq 1 ]; then
			 if [ $v -eq 1 ]; then   
   				 echo "Gunzipped output"
                         fi
			bcftools mpileup -Ou -f "$ref" "$(pwd)"/tmp/out1_sorted.bam  | bcftools call -vmO z -o "$output".vcf.gz
		else
			 bcftools mpileup -Ou -f "$ref" "$(pwd)"/tmp/out1_sorted.bam  | bcftools call -vmO v -o "$output".vcf
		fi
	fi
}
main() {
	# Function that defines the order in which functions will be called
	# You will see this construct and convention in a lot of structured code.
	# Add flow control as you see appropriate
	get_input "$@"
	check_files # Add arguments here
	prepare_temp
	mapping # Add arguments here
	improvement # Add arguments here
	call_variants # Add arguments here
}

# Calling the main function
main "$@"
