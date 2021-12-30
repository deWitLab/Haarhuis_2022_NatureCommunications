macs2 bdgbroadcall -i ../md.bg -g 65 -l 500 -G 2000 -C 1.5 -c 2.4 -o md_K9tri_broad_g65_l500_G2000.bed
awk '{if ($3-$2 > 100000) print $1"\t"$2"\t"$3}' md_K9tri_broad_g65_l500_G2000.bed > md_K9tri_broad_g65_l500_G2000_min100K.bed
bedtools merge -i md_K9tri_broad_g65_l500_G2000_min100K.bed -d 100000 > md_k9tri_broad_g65_l500_G2000_min100K_merged10k.bed

macs2 bdgbroadcall -i ../wt.bg -g 65 -l 500 -G 2000 -C 1.5 -c 2.4 -o wt_K9tri_broad_g65_l500_G2000.bed
awk '{if ($3-$2 > 100000) print $1"\t"$2"\t"$3}' wt_K9tri_broad_g65_l500_G2000.bed > wt_K9tri_broad_g65_l500_G2000_min100K.bed
bedtools merge -i wt_K9tri_broad_g65_l500_G2000_min100K.bed -d 100000 > wt_k9tri_broad_g65_l500_G2000_min100K_merged10k.bed

macs2 bdgbroadcall -i ../wa.bg -g 65 -l 500 -G 2000 -C 1.5 -c 2.4 -o wa_K9tri_broad_g65_l500_G2000.bed
awk '{if ($3-$2 > 100000) print $1"\t"$2"\t"$3}' wa_K9tri_broad_g65_l500_G2000.bed > wa_K9tri_broad_g65_l500_G2000_min100K.bed
bedtools merge -i wa_K9tri_broad_g65_l500_G2000_min100K.bed -d 100000 > wa_k9tri_broad_g65_l500_G2000_min100K_merged10k.bed

macs2 bdgbroadcall -i ../dk.bg -g 65 -l 500 -G 2000 -C 1.5 -c 2.4 -o dk_K9tri_broad_g65_l500_G2000.bed
awk '{if ($3-$2 > 100000) print $1"\t"$2"\t"$3}' dk_K9tri_broad_g65_l500_G2000.bed > dk_K9tri_broad_g65_l500_G2000_min100K.bed
bedtools merge -i dk_K9tri_broad_g65_l500_G2000_min100K.bed -d 100000 > dk_k9tri_broad_g65_l500_G2000_min100K_merged10k.bed


macs2 bdgbroadcall -i ../ccnc.bg -g 65 -l 500 -G 2000 -C 1.5 -c 2.4 -o ccnc_K9tri_broad_g65_l500_G2000.bed
awk '{if ($3-$2 > 100000) print $1"\t"$2"\t"$3}' ccnc_K9tri_broad_g65_l500_G2000.bed > ccnc_K9tri_broad_g65_l500_G2000_min100K.bed
bedtools merge -i ccnc_K9tri_broad_g65_l500_G2000_min100K.bed -d 100000 > ccnc_k9tri_broad_g65_l500_G2000_min100K_merged10k.bed
