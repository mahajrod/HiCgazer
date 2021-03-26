bwa index draft.fa
generate_site_positions.py MboI draft draft.fa
./juicer/scripts/juicer.sh –g draft –s MboI –z draft.fa –y draft_MboI.txt –p
assembly