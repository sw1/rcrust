#!/usr/bin/env Rscript

SCRIPTS <- '/usr/local/bin'
format_tree_and_trait_table <- file.path(SCRIPTS,'format_tree_and_trait_table.py')
ancestral_state_reconstruction <- file.path(SCRIPTS,'ancestral_state_reconstruction.py')
predict_traits <- file.path(SCRIPTS,'predict_traits.py')

cat('Formatting traits for 16s\n')
system2(format_tree_and_trait_table,
        args=c('-t','~/rcrust/testing_workflow/rtree.tree',
               '-i','~/rcrust/picrust_starting_files/IMG_16S_counts.tab',
               '-m','~/rcrust/gg_13_5_img.txt',
               '-o','~/rcrust/testing_workflow/out/16S/'))

cat('Formatting traits for ko\n')
system2(format_tree_and_trait_table,
        args=c('-t','~/rcrust/testing_workflow/rtree.tree',
               '-i','~/rcrust/picrust_starting_files/IMG_ko_counts.tab',
               '-m','~/rcrust/gg_13_5_img.txt',
               '-o','~/rcrust/testing_workflow/out/KEGG/'))

cat('Asr for 16s\n')
system2(ancestral_state_reconstruction,
        args=c('-t','~/rcrust/testing_workflow/out/16S/pruned_tree.newick',
               '-i','~/rcrust/testing_workflow/out/16S/trait_table.tab',
               '-o','~/rcrust/testing_workflow/out/asr/16S_asr_counts.tab'))

cat('Asr for ko\n')
system2(ancestral_state_reconstruction,
        args=c('-t','~/rcrust/testing_workflow/out/KEGG/pruned_tree.newick',
               '-i','~/rcrust/testing_workflow/out/KEGG/trait_table.tab',
               '-o','~/rcrust/testing_workflow/out/asr/KEGG_asr_counts.tab'))

cat('Predicting traits for 16s\n')
system2(predict_traits,
        args=c('-t','~/rcrust/testing_workflow/out/16S/reference_tree.newick',
               '-i','~/rcrust/testing_workflow/out/16S/trait_table.tab',
               '-r','~/rcrust/testing_workflow/out/asr/16S_asr_counts.tab',
               '-o','~/rcrust/testing_workflow/out/new_refs/16S_13_5_precalculated.tab'))

cat('Predicting traits for ko\n')
system2(predict_traits,
        args=c('-t','~/rcrust/testing_workflow/out/KEGG/reference_tree.newick',
               '-i','~/rcrust/testing_workflow/out/KEGG/trait_table.tab',
               '-r','~/rcrust/testing_workflow/out/asr/KEGG_asr_counts.tab',
               '-o','~/rcrust/testing_workflow/out/new_refs/ko_13_5_precalculated.tab'))

cat('Compressing output\n')
system('gzip ~/rcrust/testing_workflow/out/new_refs/16S_13_5_precalculated.tab')
system('gzip ~/rcrust/testing_workflow/out/new_refs/ko_13_5_precalculated.tab')
