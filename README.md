# FISH-globe

## HCR probe design

### Pilot work flow

1. Get GTF file from Ensembl into tidy valr tibble.
2. Retain protein-coding/non-coding/rRNA gene_biotypes.
3. Remove any outlier genes that have both +/- strand attributes.
4. Split into list for each gene element.
5. Identify isoform-overlapping intervals using bedtools multiinter. (later used to toggle whether the user is happy to include probes spanning exon-intron junctions and/or isoform-specific junctions). Multiinter output also contains information which region overlaps with which isoforms. 
6. 
