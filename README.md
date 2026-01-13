<h2>Leonardo da Vinci Y Haplogroup Analysis Code</h2><br/>

This repository contains the bash scripts/tool commands (Tools_Scripts_bashCommandLine Directory) as well as the python jupyter notebooks (Haplogroup_Analysis_python_notebooks Directory) used in assigning Y haplogroup lineage information to short read illumina data taken from artifacts. </br>

<h3>Tools_Scripts_bashCommandLine</h3></br>
- part1_mappingSamples.sh: Paired short-read data was mapped using BWA MEM to the T2T CHM13 human genome, the Y chromosome was replaced with the GRCh38 Y chromosome. This was done so we could use that coordinate space.</br>
- part2_callVariants.sh: bcftools mpileup was used against the Poznik Y callable regions. SNVs 5 bases from INDELs were removed from the analysis. INDEL records are filtered in python notebooks.</br>

<h3>Haplogroup_Analysis_python_notebooks</h3>
- variants_part1_read_vcfs.ipynb: Code used to compare the variant calls to the ISOGG database of Y haplogroup markers. Produces CSV files containing ancestral and derived Y haplogroup SNV counts. </br>
- variants_part2_visualize_data.ipynb: Python code to visualize/cluster artifact Y haplogroup derived SNV profiles (this code doesn't affect any results its simply to aid in visualization of the data in some informative way). 
