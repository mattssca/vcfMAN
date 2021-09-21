#!/bin/bash

if ls in/SVs/*.gz &>/dev/null; 
then 
  gunzip in/SVs/*.gz; 
else 
  echo "vcf successfully extrcted from *.gz"
fi