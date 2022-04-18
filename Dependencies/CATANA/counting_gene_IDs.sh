#!/usr/bin/env sh

# Author: Cheng-Kai Shiau <shiauck at gmail dot com>

if [ -n "$1" ]; then
   for i in `ls $1/*`;
      do echo $i;
         awk '{if($3=="gene"){print $9}}' $i | cut -d '=' -f 2 | cut -d '.' -f 1 | sort -u | wc -l
   done
else
   echo "\nUsage: $0 <folder name>\n";
fi

