#!/bin/bash

if ls in/small_variants/*.gz &>/dev/null; 
then 
  gunzip in/small_variants/*.gz; 
else 
  echo "vcf successfully extrcted from *.gz"
fi