#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.2') _

def builders = [:]

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

node {
    step([$class: 'GitHubCommitStatusSetter', reposSource:
        [$class: 'ManuallyEnteredRepositorySource',
        url: 'https://github.com/ufz/ogs5.git']])
}

properties([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator',
    artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
