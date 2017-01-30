node('master') {
    checkout scm

    parallel linux: { load 'scripts/jenkins/linux.groovy' },
    mingw: { load 'scripts/jenkins/mingw.groovy' }

    step([$class: 'GitHubCommitStatusSetter'])
}

properties([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator',
    artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
