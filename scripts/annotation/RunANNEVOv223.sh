#!/bin/bash

##############################################
## This is a script to run several QC tools:
##############################################

annevo_exe="python /data/tfg_users/mollai/software/ANNEVO/annotation.py"


######################
## STARTING MESSAGE ##
######################

echo "";
echo "===================================";
echo "Wrapper of the ANNEVO tool starts";
echo "===================================";
echo "";

##########################
## DEFINE HELP FUNCTION ##
##########################

Help()
{
        echo "Usage: RunANNEVO.sh -g <genome.fasta> -o <output.gff>";
        echo "       [-c <cuda_gpus>] [-t <threads>] [-B]";
        echo "";
        echo " Arguments:";
        echo "  -h <print_this_help>"
        echo "  -g <genome.fasta> [MANDATORY]";
        echo "  -o <output.gff> [MANDATORY]";
        echo "  -c <cuda_gpus> [Comma separated (e.g., 0,1). OPTIONAL, default: 0]";
        echo "  -t <threads> [OPTIONAL, default:8]";
        echo "  -m <annevo_model_file> [OPTIONAL, default:/data/tfg_users/lopezp/software/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt]";
        echo "  -l <linage> [{Mammalia,Insecta,Aves,Actinopteri,Magnoliopsida,Fungi}]"
        echo "  -B Enable boundary aware model (--num_classes 15 --boundary-aware --min_intron_length 20)"
        echo "";
}

###############################
## DEFINE ARGUMENTS FUNCTION ##
###############################

genome="";
outname="";
cuda="";
threads="";
models="";
lineage="";
baware=0;

while getopts ":hg:o:c:t:m:l:B" option; do
   case $option in
      h) ## Print help and exit
         Help
         exit;;
      g) ## Get the genome option
         genome=$OPTARG;;
      o) ## Get the outname option
         outname=$OPTARG;;
      c) ## Get the cuda
         cuda=$OPTARG;;
      t) ## Get the number of threads
         threads=$OPTARG;;
      m) ## Get the ANNEVO models
         models=$OPTARG;;
      l) ## Get lineage
         lineage=$OPTARG;;
      B) ## activate boundary aware
         baware=1;;
     \?) # Invalid option
         echo "";
         echo "Error: Invalid option. Please check help running the command:"
         echo "";
         echo "RunANNEVO -h";
         echo "";
         exit;;
   esac
done



###########################
## DEFINE DEFAULT VALUES ##
###########################

cuda="${cuda:-0}";
threads="${threads:-8}";
models="${models:-/data/tfg_users/mollai/software/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt}"
lineage="${lineage:-Magnoliopsida}"

## Print help if no argument is used

if [[ -z $genome ]]
then
        echo "";
        echo "No arguments used... printing usage";
        echo "";
        Help;
        exit;
else
        if ! [ -f $genome ]
        then
           echo "";
           echo "ERROR: <genome.fasta> $genome file does not exist";
           echo "";
           exit;

        else
           if [[ -z $outname ]]
           then
              echo "";
              echo "ERROR: <outname.gff> is a mandatory argument.";
              echo "";
              exit;
           fi
        fi
fi




## RUN ANNEVO

echo "1- RUN ANNEVO"
date; 

command="$annevo_exe --genome $genome --model_path $models --output "temp_"$outname --threads $threads -l $lineage"
if [[ $baware == 1 ]]
then 
    command=$command" --num_classes 15 --boundary-aware --min_intron_length 20";
fi

echo "";
echo "COMMAND TO RUN: CUDA_VISIBLE_DEVICES=$cuda $command";
echo ""; 

#CUDA_VISIBLE_DEVICES=$cuda $annevo_exe --genome $genome --model_path $models --output "temp_"$outname --threads $threads --num_classes 15 --boundary-aware --min_intron_length 20
CUDA_VISIBLE_DEVICES=$cuda $command;

echo "2- RUN AGAT TO CONVERT AND SORT THE OUTPUT"
date;

agat_convert_sp_gxf2gxf.pl -g "temp_"$outname -o "temp_formatted_"$outname
sed -r 's/transcript/mRNA/' "temp_formatted_"$outname > $outname

echo "DONE"
date;

