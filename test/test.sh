../../samtools_merge/src/samtools_merge test_2.bam test_1.bam test_3.bam badmerge.bam
../src/untangle badmerge.bam test_header.sam test_2_reconstruct.bam test_1.reads:test_1.sam:test_1_reconstruct.bam test_3.reads:test_3.sam:test_3_reconstruct.bam
../../bridgebuilder/brunel/src/brunel test_header.sam test_1_reconstruct.bam:../../bridgebuilder/brunel/test/trans.txt test_2_reconstruct.bam test_3_reconstruct.bam out_reconstruct.bam
