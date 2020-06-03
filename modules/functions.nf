#!/bin/bash nextflow

def set_templates_path () {
    if( !nextflow.version.matches('20.04.1-edge+') ) {
        println "It is advisable to run this workflow with Nextflow version 20.04.1-edge or greater -- You are running version $nextflow.version"
        path_templates = projectDir.name == "nf-benchmark" ? "$projectDir/modules/regressive_alignment/modules/templates" : "$projectDir/modules/templates"
    }
    else {
        path_templates = "${moduleDir}/templates"
    }
    return path_templates
}

