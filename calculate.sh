#!/bin/bash
set -o nounset
set -o pipefail

cd $1

ls | sort

module load samtools/1.9 > /dev/null

find -name *.bam -exec samtools view -H {} \; | grep '^@RG' | sort

find -name *.bam -exec samtools flagstat {} \; | sort