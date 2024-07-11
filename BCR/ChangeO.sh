AssignGenes.py igblast -s npc_merge_bcr.fasta -b /data/home/liuyang/software/ncbi-igblast-1.17.0 \
--organism human --loci ig --format blast

MakeDb.py igblast -i npc_merge_bcr_igblast.fmt7 -s npc_merge_bcr.fasta \
-r /data/home/liuyang/liuyang/RAW_data/1-SingleCell/3-NPC/analysis/NPC-BCR/Change-O/database/AIRR_Example/IMGT_Human_*.fasta \
--10x bcr_merge_10x.csv --extended

ParseDb.py select -d npc_merge_bcr_igblast_db-pass.tsv -f productive -u T
#ParseDb.py split -d npc_merge_bcr_igblast_db-pass.tsv -f productive

ParseDb.py select -d npc_merge_bcr_igblast_db-pass.tsv -f v_call j_call c_call -u "IGH" \
    --logic all --regex --outname heavy
	
cut_off = 0.12

DefineClones.py -d npc_merge_bcr_igblast_db-pass_parse-select.tsv --act set --model ham \
--norm len --dist 0.12

CreateGermlines.py -d npc_merge_bcr_igblast_db-pass_parse-select_clone-pass.tsv -g dmask --cloned \
    -r IMGT_Human_IGHV.fasta IMGT_Human_IGHD.fasta IMGT_Human_IGHJ.fasta
	
BuildTrees.py -d npc_merge_bcr_igblast_db-pass_parse-select_clone-pass_germ-pass.tsv --outname ex --log ex.log --collapse \
    --sample 3000 --igphyml --clean all --nproc 16
	
