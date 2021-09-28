#!/bin/bash

#combine png for report format
convert -quiet out/SVs/report.png  out/small_variants/report.png -append -quiet out/report.png