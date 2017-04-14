SCRIPTS <- '/usr/local/bin'
format_tree_and_trait_table <- file.path(SCRIPTS,'format_tree_and_trait_table.py')
ancestral_state_reconstruction <- file.path(SCRIPTS,'ancestral_state_reconstruction.py')
predict_traits <- file.path(SCRIPTS,'predict_traits.py')

system2(format_tree_and_trait_table,
        args=c('-t','~/rcrust/tree.tree',
               '-i','~/rcrust/picrust_starting_files/IMG_16S_counts.tab',
               '-m','~/rcrust/gg_13_5_img.txt',
               '-o','~/rcrust/out/16S/'))

system2(format_tree_and_trait_table,
        args=c('-t','~/rcrust/tree.tree',
               '-i','~/rcrust/picrust_starting_files/IMG_ko_counts.tab',
               '-m','~/rcrust/gg_13_5_img.txt',
               '-o','~/rcrust/out/KEGG/'))

system2(ancestral_state_reconstruction,
        args=c('-t','~/rcrust/out/16S/pruned_tree.newick',
               '-i','~/rcrust/out/16S/trait_table.tab',
               '-o','~/rcrust/asr/16S_asr_counts.tab'))

system2(ancestral_state_reconstruction,
        args=c('-t','~/rcrust/out/KEGG/pruned_tree.newick',
               '-i','~/rcrust/out/KEGG/trait_table.tab',
               '-o','~/rcrust/asr/KEGG_asr_counts.tab'))

system2(predict_traits,
        args=c('-t','~/rcrust/out/16S/reference_tree.newick',
               '-i','~/rcrust/out/16S/trait_table.tab',
               '-r','~/rcrust/asr/16S_asr_counts.tab',
               '-o','~/rcrust/new_refs/16S_13_5_precalculated.tab'))

system2(predict_traits,
        args=c('-t','~/rcrust/out/KEGG/reference_tree.newick',
               '-i','~/rcrust/out/KEGG/trait_table.tab',
               '-r','~/rcrust/asr/KEGG_asr_counts.tab',
               '-o','~/rcrust/new_refs/ko_13_5_precalculated.tab'))

system('gzip ~/rcrust/new_refs/16S_13_5_precalculated.tab')
system('gzip ~/rcrust/new_refs/ko_13_5_precalculated.tab')
