#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process TCOFFEE_DEFAULT{
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee default  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_default.aln
    """
} 
process TCOFFEE_QUICKALN{
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee quickaln  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode quickaln $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_quickaln.aln 
    """
} 
process TCOFFEE_MCOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee mcoffee  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode mcoffee $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_mcoffee.aln
    """
} 
process TCOFFEE_ACCURATE {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee accurate  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode accurate -blast LOCAL -pdb_db ${params.database_path} $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_accurate.aln
    """
} 
process TCOFFEE_FMCOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee fmcoffee  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode fmcoffee $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_fmcoffee.aln
    """
}  
process TCOFFEE_PSICOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee psicoffee  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    def cache_filter = params.blastOutdir != null ? "-cache ${params.blastOutdir}" : ''
    """
    t_coffee -seq $seqs -mode psicoffee -blast_server LOCAL -protein_db ${params.database_path} $cache_filter $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_psicoffee.aln
    """
} 
process TCOFFEE_EXPRESSO {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee expresso  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode expresso -pdb_type d -blast LOCAL -pdb_db ${params.database_path} $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_expresso.aln
    """
}           
process TCOFFEE_PROCOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee procoffee  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode procoffee $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_procoffee.aln
    """
}    
process TCOFFEE_3DCOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee 3dcoffee  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -method sap_pair $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_3dcoffee.aln
    """
} 
process TCOFFEE_TRMSD {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee trmsd  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:     // TODO -> check to output
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : '' 
    """
    t_coffee -seq $seqs -method mustang_pair $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_trmsd.aln

    t_coffee -other_pg trmsd -aln ${id}.tcoffee_trmsd.aln
    """
} 
process TCOFFEE_RCOFFEE {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee rcoffee  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode rcoffee $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_rcoffee.aln
    """
}        
process TCOFFEE_RCOFFEE_CONSAN {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee rcoffee_consan  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -mode rcoffee_consan $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_rcoffee_consan.aln
    """
} 
process TCOFFEE_3DALIGN {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee 3d_align  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -method sap_pair TMalign_pair -output fasta_aln $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_3d_align.aln
    """
} 
process TCOFFEE_3DMALIGN {
    container 'edgano/tcoffee:protocols'
    tag "Tcoffee 3d_align  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), file(seqs), file(template), file(library)
    val(fakeId)             //to ensure the process before when precompute Blast

    output:
    tuple val (id), path ("*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    def template_filter = template.name != 'input.2' ? "-template_file $template" : ''
    def libs_filter = library.name != 'input.3' ? "-lib $library" : ''
    """
    t_coffee -seq $seqs -method sap_pair TMalign_pair mustang_pair -output fasta_aln $template_filter $libs_filter ${params.params4tcoffee} -outfile ${id}.tcoffee_3dM_align.aln
    """
}  

