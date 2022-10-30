gcc -o methylate_cg methylate_cg_region.c mcg_region/kmer2binary_cg.c mcg_region/viterbi_algorithm_cg.c -Lhtslib/ -Ihtslib/ -lhts -lm -lz
gcc -o methylate_ma methylate_ma_region.c ma_region/kmer2binary.c ma_region/viterbi_algorithm.c -Lhtslib/ -Ihtslib/ -lhts -lm -lz
gcc -o re_align.ma re_align_ma.c training/ma/kmer2binary.c training/ma/viterbi_algorithm.c -Lhtslib/ -Ihtslib/ -lhts -lm -lz
gcc -o re_align.mcg re_align_mcg.c training/mcg/kmer2binary.c training/mcg/viterbi_algorithm.c -Lhtslib/ -Ihtslib/ -lhts -lm -lz
gcc -o re_align.ma_mcg re_align_ma_mcg.c training/ma_mcg/kmer2binary.c training/ma_mcg/viterbi_algorithm.c -Lhtslib/ -Ihtslib/ -lhts -lm -lz
