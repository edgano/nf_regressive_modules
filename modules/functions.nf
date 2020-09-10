#!/bin/bash nextflow

nxf_min_ver = "20.07.1"

def set_templates_path () {
    if( !nextflow.version.matches(">=${nxf_min_ver}") ) {
        println "It is advisable to run this workflow with Nextflow version $nxf_min_ver or greater -- You are running version $nextflow.version"
        path_templates = workflow.commandLine.contains('--pipeline') ? "$baseDir/modules/regressive_alignment/modules/templates" : "$baseDir/modules/templates"
    }
    else {
        path_templates = "${moduleDir}/templates"
    }
    return path_templates
}

