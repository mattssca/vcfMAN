#!/bin/bash
bedtools merge -i out/SVs/tables/bed/sv_del.bed > out/SVs/tables/bed/sv_del_merged.bed
bedtools merge -i out/SVs/tables/bed/sv_dup.bed > out/SVs/tables/bed/sv_dup_merged.bed
bedtools merge -i out/SVs/tables/bed/sv_ins.bed > out/SVs/tables/bed/sv_ins_merged.bed

rm out/SVs/tables/bed/sv_del.bed
rm out/SVs/tables/bed/sv_dup.bed
rm out/SVs/tables/bed/sv_ins.bed