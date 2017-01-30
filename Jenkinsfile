#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.2') _

node('master') {
    checkout scm

    parallel linux: { load 'scripts/jenkins/linux.groovy' },
    mingw: { load 'scripts/jenkins/mingw.groovy' }

    step([$class: 'GitHubCommitStatusSetter', reposSource:
        [$class: 'ManuallyEnteredRepositorySource',
        url: 'https://github.com/ufz/ogs5.git']])
}

properties([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator',
    artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
