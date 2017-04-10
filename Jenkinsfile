#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.2') _

def builders = [:]

timestamps {

builders['linux'] = {
    node('envinf1') {
        dir('ogs') { checkoutWithTags() }
        load 'ogs/scripts/jenkins/linux.groovy'
    }
}

builders['mingw'] = {
    node('docker') {
        dir('ogs') { checkoutWithTags() }
        load 'ogs/scripts/jenkins/mingw.groovy'
    }
}

parallel builders

} // end timestamps
