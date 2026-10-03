# Reproducibility scope and limitations

This package reproduces the read-mapping validation, per-base coverage calculations, coverage statistics, primer-to-final-genome remapping, TP13 ABI quality trimming/alignment, and all coverage figures generated during the analysis.

The mapping workflow is a custom validation script, not BWA-MEM2, Bowtie2, SAMtools, or BCFtools. It uses exact 19-nt seed lookup followed by full-read Hamming mismatch scoring on a circularized reference. It does not model insertions/deletions and should not be presented as a general-purpose short-read aligner. It was used here because the reads are 150-bp sequences generated from predefined overlapping mitochondrial PCR amplicons and the goal was to check read support against the completed sequence.

CAP3 v10.2011 is retained in the manuscript as the primary assembly program because that is the assembly method recorded in the deposited GenBank record. The original intermediate CAP3 assembly files and any historical draft containing the reported run of Ns were not uploaded to this chat. Consequently, this package does not claim to reproduce the historical assembly process byte-for-byte or to recover the exact historical N-block coordinates. The 14,937–15,261 interval is explicitly an empirical lowest-coverage 325-bp window calculated from the mapped Illumina reads.

The robust reproducible Sanger interval in the notebook is the Biopython `abi-trim` high-quality core (651 bp), which aligns from final-genome position 14,968 through the circular origin to position 123. Exploratory alignment of untrimmed chromatogram sequence is quality/scoring dependent and is therefore not used as the primary manuscript coordinate.
